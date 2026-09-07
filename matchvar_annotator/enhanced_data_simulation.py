#!/usr/bin/env python3
"""增强的数据模拟模块：在 GeneTranscript 上叠加突变谱感知的 SNV 生成。"""

from __future__ import annotations

import logging
from typing import Any, Dict, List, Optional

from .mutation_spectrum_generator import (
    ContextAwareVariantGenerator,
    MutationProbability,
    create_enhanced_variant_generator,
)
from .variant_simulation import ExonExtractor, GeneTranscript

logger = logging.getLogger(__name__)


class EnhancedGeneTranscript(GeneTranscript):
    """继承传统仿真，增加基于突变谱的变异生成。"""

    def __init__(self, *args, **kwargs):
        self.genome_version = kwargs.pop("genome_version", None) or "hg19"
        super().__init__(*args, **kwargs)
        self.enhanced_generator: Optional[ContextAwareVariantGenerator] = None
        self._init_enhanced_generator()

    def _init_enhanced_generator(self) -> None:
        try:
            self.enhanced_generator = create_enhanced_variant_generator(
                genome_version=self.genome_version,
                genome_fasta=getattr(self, "fasta_file", None),
                genome_handle=getattr(self, "genome", None),
            )
            logger.info("增强的变异生成器初始化成功（genome=%s）", self.genome_version)
        except Exception as exc:
            logger.error("初始化增强变异生成器失败: %s", exc)
            self.enhanced_generator = None

    @classmethod
    def from_gtf(
        cls,
        gene_name: str,
        transcript_id: str,
        gtf_file: str,
        fasta_file: str,
        genome_version: str = "hg19",
    ):
        extractor = ExonExtractor(gtf_file, fasta_file)
        exons, chromosome, strand = extractor.extract_exons(gene_name, transcript_id)
        return cls(
            gene_name=gene_name,
            transcript_id=transcript_id,
            exons=exons,
            chromosome=chromosome,
            strand=strand,
            genome=extractor.genome,
            fasta_file=fasta_file,
            genome_version=genome_version,
        )

    def generate_spectrum_aware_variants(
        self,
        max_variants: int = 1000,
        min_probability: float = 0.001,
        use_constraint_filter: bool = True,
        sample_size: Optional[int] = None,
    ) -> Dict[str, Any]:
        if not self.enhanced_generator:
            logger.error("增强变异生成器未初始化")
            return {"variants": [], "total": 0, "error": "Generator not initialized"}

        try:
            variants = self.enhanced_generator.generate_context_aware_variants(
                gene_transcript=self,
                max_variants=max_variants,
                min_probability=min_probability,
                use_constraint_filter=use_constraint_filter,
                sort_by_probability=False,
            )
            if sample_size and sample_size < len(variants):
                variants = self.enhanced_generator.sample_variants_by_probability(
                    variants, sample_size
                )
            formatted_variants = self._format_spectrum_variants(variants)
            logger.info("生成了 %d 个基于突变谱的变异", len(formatted_variants))
            return {
                "variants": formatted_variants,
                "total": len(formatted_variants),
                "spectrum_aware": True,
                "constraint_filtered": use_constraint_filter,
            }
        except Exception as exc:
            logger.error("生成基于突变谱的变异失败: %s", exc)
            return {"variants": [], "total": 0, "error": str(exc)}

    def _format_spectrum_variants(self, variants: List[MutationProbability]) -> List[Dict]:
        formatted: List[Dict] = []
        for i, variant in enumerate(variants):
            try:
                if self.strand == "-":
                    hgvs_ref = self._reverse_complement(variant.ref_base)
                    hgvs_alt = self._reverse_complement(variant.alt_base)
                else:
                    hgvs_ref = variant.ref_base
                    hgvs_alt = variant.alt_base

                hgvs = f"c.{variant.context.cds_position}{hgvs_ref}>{hgvs_alt}"
                is_synonymous = self.is_synonymous(
                    variant.context.cds_position - 1,
                    hgvs_alt,
                )
                codon_start = ((variant.context.cds_position - 1) // 3) * 3
                codon_pos = (variant.context.cds_position - 1) % 3
                original_codon = self.cds_sequence[codon_start:codon_start + 3]
                mutant_codon = original_codon[:codon_pos] + hgvs_alt + original_codon[codon_pos + 1:]
                is_stop = self.CODON_TABLE.get(mutant_codon) == "*"
                score = self._calculate_snv_score(
                    cds_pos=variant.context.cds_position - 1,
                    is_syn=is_synonymous,
                    is_stop=is_stop,
                    codon_pos=codon_pos + 1,
                )
                formatted.append({
                    "type": "SNV",
                    "hgvs": hgvs,
                    "cds_position": variant.context.cds_position,
                    "original": variant.ref_base,
                    "mutant": variant.alt_base,
                    "synonymous": is_synonymous,
                    "is_stop_codon": is_stop,
                    "genomic_pos": variant.context.genomic_position - 1,
                    "ref_seq": variant.ref_base,
                    "alt_seq": variant.alt_base,
                    "hgvs_ref": hgvs_ref,
                    "hgvs_alt": hgvs_alt,
                    "mutation_probability": variant.probability,
                    "trinucleotide_context": variant.context.trinucleotide_context,
                    "is_constrained_gene": variant.is_constrained,
                    "gnomad_constraint_score": variant.gnomad_constraint,
                    "Total_Score": round(score, 2),
                })
            except Exception as exc:
                logger.warning("格式化变异 %s 失败: %s", i, exc)
                continue
        return formatted

    def generate_hybrid_variants(
        self,
        spectrum_ratio: float = 0.7,
        max_variants: Optional[int] = 1000,
        min_probability: float = 0.001,
        use_constraint_filter: bool = True,
    ) -> Dict[str, Any]:
        if max_variants is None:
            spectrum_count = None
            traditional_count = None
        else:
            spectrum_count = int(max_variants * spectrum_ratio)
            traditional_count = max_variants - spectrum_count

        results: Dict[str, Any] = {
            "spectrum_variants": [],
            "traditional_variants": [],
            "total": 0,
            "hybrid": True,
        }

        if spectrum_count is None or spectrum_count > 0:
            spectrum_result = self.generate_spectrum_aware_variants(
                max_variants=spectrum_count,
                min_probability=min_probability,
                use_constraint_filter=use_constraint_filter,
            )
            results["spectrum_variants"] = spectrum_result.get("variants", [])

        if traditional_count is None or traditional_count > 0:
            snv_result = self.generate_snv_variants(
                synonymous=True,
                include_stop_codon=True,
            )
            splice_result = self.generate_splice_variants(
                max_intron_offset=20,
                min_intron_offset=1,
                include_classic_sites=True,
            )
            if traditional_count is None:
                results["traditional_variants"] = snv_result + splice_result
            else:
                snv_limit = int(traditional_count * 0.7)
                splice_limit = traditional_count - snv_limit
                results["traditional_variants"] = (
                    (snv_result[:snv_limit] if snv_limit > 0 else [])
                    + (splice_result[:splice_limit] if splice_limit > 0 else [])
                )

        all_variants = results["spectrum_variants"] + results["traditional_variants"]
        results["total"] = len(all_variants)
        results["all_variants"] = all_variants
        logger.info(
            "生成了混合变异: %d 个基于突变谱, %d 个传统枚举",
            len(results["spectrum_variants"]),
            len(results["traditional_variants"]),
        )
        return results


def create_enhanced_gene_transcript(
    gene_name: str,
    transcript_id: str,
    exons: List[Dict],
    chromosome: str = "chr1",
    strand: str = "+",
    genome_version: str = "hg19",
    utr5: Optional[Dict] = None,
    utr3: Optional[Dict] = None,
    genome=None,
    fasta_file: Optional[str] = None,
) -> EnhancedGeneTranscript:
    return EnhancedGeneTranscript(
        gene_name=gene_name,
        transcript_id=transcript_id,
        exons=exons,
        chromosome=chromosome,
        strand=strand,
        utr5=utr5,
        utr3=utr3,
        genome=genome,
        fasta_file=fasta_file,
        genome_version=genome_version or "hg19",
    )


def create_enhanced_gene_transcript_from_gtf(
    gene_name: str,
    transcript_id: Optional[str] = None,
    gtf_file: str = None,
    fasta_file: str = None,
    genome_version: str = "hg19",
) -> Optional[EnhancedGeneTranscript]:
    if not gtf_file or not fasta_file or not transcript_id:
        logger.error("from_gtf 需要 gtf_file、fasta_file 和 transcript_id")
        return None
    try:
        return EnhancedGeneTranscript.from_gtf(
            gene_name=gene_name,
            transcript_id=transcript_id,
            gtf_file=gtf_file,
            fasta_file=fasta_file,
            genome_version=genome_version,
        )
    except Exception as exc:
        logger.error("从 GTF 创建增强基因转录本失败: %s", exc)
        return None
