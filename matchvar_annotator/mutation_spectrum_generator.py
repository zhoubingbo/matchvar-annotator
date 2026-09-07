#!/usr/bin/env python3
"""基于人类突变谱模型的变异生成器。

整合三核苷酸上下文突变率 + gnomAD 约束分数。参考序列优先使用
GeneTranscript.genome / FASTA，无参考时回退到 CDS 序列推断。
"""

from __future__ import annotations

import json
import logging
import os
from dataclasses import dataclass
from typing import Dict, List, Optional

import numpy as np

logger = logging.getLogger(__name__)


@dataclass
class MutationContext:
    """突变上下文信息"""
    trinucleotide_context: str
    center_base: str
    flanking_bases: str
    genomic_position: int
    cds_position: int
    chromosome: str


@dataclass
class MutationProbability:
    """突变概率信息"""
    ref_base: str
    alt_base: str
    probability: float
    context: MutationContext
    gnomad_constraint: Optional[float] = None
    is_constrained: bool = False


_RC = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def _reverse_complement(seq: str) -> str:
    return seq.translate(_RC)[::-1]


def _fetch_fasta_slice(fasta, chromosome: str, start0: int, end0: int) -> Optional[str]:
    names = (
        chromosome,
        f"chr{chromosome}" if not str(chromosome).startswith("chr") else chromosome,
        str(chromosome).replace("chr", "", 1) if str(chromosome).startswith("chr") else chromosome,
    )
    for name in names:
        try:
            rec = fasta[name][start0:end0]
            seq = str(getattr(rec, "seq", rec)).upper()
            if len(seq) == (end0 - start0):
                return seq
        except Exception:
            continue
    return None


class MuRaLModel:
    """基于序列上下文的三核苷酸突变率模型。"""

    def __init__(self, model_file: Optional[str] = None):
        self.model_file = model_file
        self.trinucleotide_rates = self._load_mural_model()

    def _load_mural_model(self) -> Dict[str, Dict[str, float]]:
        if self.model_file and os.path.exists(self.model_file):
            try:
                return self._load_from_file()
            except Exception as exc:
                logger.warning("加载 MuRaL 模型文件失败，改用默认速率: %s", exc)
        return self._get_default_mural_rates()

    def _load_from_file(self) -> Dict[str, Dict[str, float]]:
        with open(self.model_file, "r", encoding="utf-8") as fh:
            data = json.load(fh)
        if not isinstance(data, dict):
            raise ValueError("MuRaL model file must be a JSON object")
        return data

    def _get_default_mural_rates(self) -> Dict[str, Dict[str, float]]:
        rates: Dict[str, Dict[str, float]] = {}
        nucleotides = ["A", "T", "C", "G"]
        for n1 in nucleotides:
            for n2 in nucleotides:
                for n3 in nucleotides:
                    context = n1 + n2 + n3
                    rates[context] = {}
                    for alt in nucleotides:
                        if alt != n2:
                            mutation = f"{n2}>{alt}"
                            rates[context][mutation] = self._get_context_specific_rate(context, mutation)
        return rates

    def _get_context_specific_rate(self, context: str, mutation: str) -> float:
        base_rates = {
            "C>T": 0.3,
            "T>C": 0.1,
            "C>A": 0.1,
            "T>A": 0.05,
            "C>G": 0.05,
            "T>G": 0.05,
            "A>T": 0.1,
            "A>C": 0.1,
            "A>G": 0.1,
            "G>T": 0.1,
            "G>C": 0.1,
            "G>A": 0.1,
        }
        base_rate = base_rates.get(mutation, 0.05)
        if context[1] == "C" and context[2] == "G" and mutation == "C>T":
            base_rate *= 10
        if context in {"TCA", "TCT", "TCC", "TCG"} and mutation == "C>T":
            base_rate *= 2
        return base_rate

    def get_mutation_probability(self, context: MutationContext, alt_base: str) -> float:
        if context.center_base == alt_base:
            return 0.0
        mutation = f"{context.center_base}>{alt_base}"
        trinucleotide = context.trinucleotide_context
        if trinucleotide in self.trinucleotide_rates:
            return self.trinucleotide_rates[trinucleotide].get(mutation, 0.0)
        return 0.0


class GnomADConstraintFilter:
    """gnomAD 约束分数过滤器；找不到文件时使用内置默认分。"""

    def __init__(self, constraint_file: Optional[str] = None):
        self.constraint_file = constraint_file
        self.constraint_scores = self._load_constraint_scores()

    def _candidate_files(self) -> List[str]:
        names = [
            os.path.join("gnomad", "hg19", "constraint_scores.json"),
            os.path.join("gnomad", "constraint_scores.json"),
        ]
        roots: List[str] = []
        env = os.environ.get("MATCHVAR_RESOURCES_DIR")
        if env:
            roots.append(env)
        pkg_root = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
        roots.extend([
            os.path.join(pkg_root, "resources"),
            os.path.join(os.getcwd(), "resources"),
        ])
        out: List[str] = []
        if self.constraint_file:
            out.append(self.constraint_file)
        for root in roots:
            for name in names:
                out.append(os.path.join(root, name))
        return out

    def _load_constraint_scores(self) -> Dict[str, float]:
        for path in self._candidate_files():
            if not path or not os.path.exists(path):
                continue
            try:
                with open(path, "r", encoding="utf-8") as fh:
                    scores = json.load(fh)
                if isinstance(scores, dict) and scores:
                    logger.info("加载 gnomAD 约束分数: %s (%d genes)", path, len(scores))
                    return {str(k): float(v) for k, v in scores.items()}
            except Exception as exc:
                logger.warning("加载 gnomAD 约束分数失败 %s: %s", path, exc)
        logger.info("使用默认 gnomAD 约束分数")
        return self._get_default_constraint_scores()

    def _get_default_constraint_scores(self) -> Dict[str, float]:
        return {
            "NIPBL": 0.032,
            "TP53": 0.469,
            "RB1": 0.10,
            "APC": 0.161,
            "MSH2": 0.334,
            "STK11": 0.245,
            "BRCA2": 0.635,
            "MLH1": 0.575,
            "MSH6": 0.498,
            "PTEN": 0.507,
            "CDH1": 0.430,
            "FANCM": 0.593,
            "BRCA1": 0.915,
            "VHL": 0.927,
            "PMS2": 1.266,
            "ATM": 0.710,
            "CHEK2": 1.530,
            "PALB2": 1.006,
            "BARD1": 1.358,
            "BRIP1": 0.786,
            "RAD51C": 1.491,
            "RAD51D": 1.218,
            "NBN": 1.010,
            "MRE11A": 0.839,
            "RAD50": 0.873,
            "FANCC": 1.043,
            "FANCG": 1.091,
            "FANCA": 1.366,
        }

    def is_constrained_gene(self, gene_name: str, threshold: float = 0.35) -> bool:
        if gene_name in self.constraint_scores:
            return self.constraint_scores[gene_name] < threshold
        return False

    def get_constraint_score(self, gene_name: str) -> float:
        return self.constraint_scores.get(gene_name, 0.5)


class ContextAwareVariantGenerator:
    """基于序列上下文和突变谱的变异生成器。"""

    def __init__(
        self,
        genome_version: str = "hg19",
        genome_fasta: Optional[str] = None,
        genome_handle=None,
        mural_model: Optional[MuRaLModel] = None,
        gnomad_filter: Optional[GnomADConstraintFilter] = None,
    ):
        self.genome_version = genome_version or "hg19"
        self.genome_fasta = genome_fasta
        self._fasta = genome_handle
        self.mural_model = mural_model or MuRaLModel()
        self.gnomad_filter = gnomad_filter or GnomADConstraintFilter()

        if self._fasta is None and genome_fasta and os.path.exists(genome_fasta):
            try:
                from pyfaidx import Fasta
                self._fasta = Fasta(genome_fasta)
            except Exception as exc:
                logger.warning("打开 FASTA 失败，将使用 CDS 回退: %s", exc)

    def extract_trinucleotide_context(
        self,
        chromosome: str,
        position: int,
        gene_transcript=None,
        cds_pos: Optional[int] = None,
    ) -> MutationContext:
        """提取三核苷酸上下文。position 为 1-based 基因组坐标。"""
        start0 = position - 2
        end0 = position + 1
        trinucleotide = None

        fasta_sources = []
        if gene_transcript is not None:
            genome = getattr(gene_transcript, "genome", None)
            if genome is not None:
                fasta_sources.append(genome)
        if self._fasta is not None:
            fasta_sources.append(self._fasta)

        for fasta in fasta_sources:
            trinucleotide = _fetch_fasta_slice(fasta, chromosome, start0, end0)
            if trinucleotide:
                break

        if (not trinucleotide or len(trinucleotide) != 3) and gene_transcript is not None and cds_pos is not None:
            trinucleotide = self._trinucleotide_from_cds(gene_transcript, cds_pos)

        if not trinucleotide or len(trinucleotide) != 3:
            raise ValueError(f"无法提取完整的三核苷酸上下文: {trinucleotide}")

        return MutationContext(
            trinucleotide_context=trinucleotide,
            center_base=trinucleotide[1],
            flanking_bases=trinucleotide[0] + trinucleotide[2],
            genomic_position=position,
            cds_position=0,
            chromosome=chromosome,
        )

    @staticmethod
    def _trinucleotide_from_cds(gene_transcript, cds_pos: int) -> str:
        seq = str(getattr(gene_transcript, "cds_sequence", "") or "").upper()
        if not seq or cds_pos < 0 or cds_pos >= len(seq):
            return ""
        left = seq[cds_pos - 1] if cds_pos > 0 else "N"
        center = seq[cds_pos]
        right = seq[cds_pos + 1] if cds_pos + 1 < len(seq) else "N"
        tri = left + center + right
        strand = getattr(gene_transcript, "strand", "+")
        if strand == "-":
            return _reverse_complement(tri)
        return tri

    def generate_context_aware_variants(
        self,
        gene_transcript,
        max_variants: int = 1000,
        min_probability: float = 0.001,
        use_constraint_filter: bool = True,
        sort_by_probability: bool = True,
    ) -> List[MutationProbability]:
        variants: List[MutationProbability] = []
        constraint_score = 0.5
        if use_constraint_filter:
            constraint_score = self.gnomad_filter.get_constraint_score(gene_transcript.gene_name)
            if constraint_score < 0.35:
                logger.info("基因 %s 是约束基因 (LOEUF=%.3f)", gene_transcript.gene_name, constraint_score)
            else:
                logger.info("基因 %s 不是约束基因", gene_transcript.gene_name)

        for cds_pos in range(gene_transcript.cds_length):
            try:
                genomic_pos = gene_transcript.get_genomic_coord(cds_pos)
                context = self.extract_trinucleotide_context(
                    gene_transcript.chromosome,
                    genomic_pos + 1,
                    gene_transcript=gene_transcript,
                    cds_pos=cds_pos,
                )
                context.cds_position = cds_pos + 1

                for alt_base in ("A", "T", "C", "G"):
                    if alt_base == context.center_base:
                        continue
                    probability = self.mural_model.get_mutation_probability(context, alt_base)
                    if probability < min_probability:
                        continue
                    is_constrained = False
                    if use_constraint_filter:
                        is_constrained = self.gnomad_filter.is_constrained_gene(gene_transcript.gene_name)
                    variants.append(
                        MutationProbability(
                            ref_base=context.center_base,
                            alt_base=alt_base,
                            probability=probability,
                            context=context,
                            gnomad_constraint=constraint_score,
                            is_constrained=is_constrained,
                        )
                    )
            except Exception as exc:
                logger.warning("处理 CDS 位置 %s 时出错: %s", cds_pos, exc)
                continue

        if sort_by_probability:
            variants.sort(key=lambda x: x.probability, reverse=True)
        else:
            variants.sort(key=lambda x: x.context.cds_position)

        if max_variants:
            variants = variants[:max_variants]

        logger.info("生成了 %d 个基于上下文的变异", len(variants))
        return variants

    def sample_variants_by_probability(
        self,
        variants: List[MutationProbability],
        sample_size: int,
    ) -> List[MutationProbability]:
        if not variants:
            return []
        probabilities = [v.probability for v in variants]
        total_prob = sum(probabilities)
        if total_prob == 0:
            return []
        normalized = [p / total_prob for p in probabilities]
        sampled_indices = np.random.choice(
            len(variants),
            size=min(sample_size, len(variants)),
            replace=False,
            p=normalized,
        )
        return [variants[i] for i in sampled_indices]


def create_enhanced_variant_generator(
    genome_version: str = "hg19",
    genome_fasta: Optional[str] = None,
    genome_handle=None,
) -> ContextAwareVariantGenerator:
    return ContextAwareVariantGenerator(
        genome_version=genome_version,
        genome_fasta=genome_fasta,
        genome_handle=genome_handle,
        mural_model=MuRaLModel(),
        gnomad_filter=GnomADConstraintFilter(),
    )
