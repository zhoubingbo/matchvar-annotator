#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
MATCHVAR Pipeline Module

End-to-end pipeline for variant simulation, annotation, and evaluation:
1. Simulate variants from GTF + gene/transcript → VCF with biological scores
2. Annotate simulated VCF with external prediction tools
3. Calculate auROC scores comparing predictions vs simulated truth
4. Generate publication-quality visualizations
"""

import os
import sys
import logging
import tempfile
import argparse
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Any
from datetime import datetime
import json

import pandas as pd
import numpy as np

from .variant_simulation import GeneTranscript, ExonExtractor
from .table_matchvar import TableAnnotator

logger = logging.getLogger(__name__)


def simulate_variants(gtf_file: str,
                       fasta_file: str,
                       gene_name: str,
                       transcript_id: str,
                       output_vcf: str,
                       variant_types: Optional[List[str]] = None,
                       max_indel_length: int = 5,
                       min_indel_length: int = 1,
                       synonymous: bool = True,
                       include_stop_codon: bool = True,
                       max_splice_offset: int = 20,
                       min_splice_offset: int = 1,
                       include_classic_splice_sites: bool = True,
                       max_variants: Optional[int] = None) -> Tuple[GeneTranscript, Dict]:
    """
    Extract exon data from GTF/FASTA, generate variants, and export to VCF.

    Args:
        gtf_file: Path to GTF annotation file
        fasta_file: Path to reference genome FASTA
        gene_name: Target gene symbol
        transcript_id: Target transcript ID
        output_vcf: Path for output VCF file
        variant_types: Types of variants to simulate (default: ['SNV', 'insertion', 'deletion'])
        max_indel_length: Maximum indel length (default: 15)
        min_indel_length: Minimum indel length (default: 1)
        synonymous: Whether to generate synonymous SNVs (default: True)
        include_stop_codon: Whether to include stop-codon variants (default: True)
        max_splice_offset: Maximum intron offset for splice sites (default: 20)
        min_splice_offset: Minimum intron offset for splice sites (default: 1)
        include_classic_splice_sites: Include classic splice sites (±1/±2) (default: True)
        max_variants: Maximum total variants to generate (default: no limit)

    Returns:
        GeneTranscript instance ready for further operations.
    """
    # Step 1: Extract exon/CDS structure from GTF and reference genome
    extractor = ExonExtractor(gtf_file, fasta_file)
    exons, chromosome, strand = extractor.extract_exons(gene_name, transcript_id)

    # Step 2: Build GeneTranscript object (pass genome and fasta_file so
    # reference-base fallbacks use the user-supplied FASTA, not a hardcoded path)
    transcript = GeneTranscript(
        gene_name=gene_name,
        transcript_id=transcript_id,
        exons=exons,
        chromosome=chromosome,
        strand=strand,
        genome=extractor.genome,
        fasta_file=fasta_file,
    )

    # Step 3: Generate variants
    variants = transcript.generate_all_variants(
        variant_types=variant_types,
        max_indel_length=max_indel_length,
        min_indel_length=min_indel_length,
        synonymous=synonymous,
        include_stop_codon=include_stop_codon,
        max_splice_offset=max_splice_offset,
        min_splice_offset=min_splice_offset,
        include_classic_splice_sites=include_classic_splice_sites,
        max_variants=max_variants,
    )

    # Step 4: Export to VCF
    transcript.export_to_vcf(variants, output_vcf)

    return transcript, variants


class MatchingPipeline:
    """
    End-to-end pipeline for variant simulation and annotation evaluation

    This class orchestrates the complete workflow:
    GTF + Gene/Transcript → Variant Simulation → Table Annotation → auROC Evaluation → Visualization
    """

    def __init__(self,
                  gtf_file: str,
                  fasta_file: str,
                  gene_name: str,
                  transcript_id: str,
                  database_dir: str,
                  output_dir: str,
                  protocols: Optional[List[str]] = None,
                  operations: Optional[List[str]] = None,
                  variant_types: Optional[List[str]] = None,
                  buildver: str = 'hg19',
                  threads: int = 4,
                  max_indel_length: int = 6,
                  min_indel_length: int = 1,
                  synonymous: bool = True,
                  include_stop_codon: bool = True,
                  max_splice_offset: int = 20,
                  min_splice_offset: int = 1,
                  include_classic_splice_sites: bool = True,
                  max_variants: Optional[int] = None):
        """
        Initialize the matching pipeline

        Args:
            gtf_file: Path to GTF annotation file
            fasta_file: Path to reference genome FASTA
            gene_name: Target gene symbol
            transcript_id: Target transcript ID
            database_dir: Directory containing annotation databases
            output_dir: Directory for all output files
            protocols: List of annotation protocols (default: ['refGene'])
            operations: List of operations (default: ['g'])
            variant_types: Types of variants to simulate (default: ['SNV', 'insertion', 'deletion'])
            buildver: Genome version (hg19/hg38)
            threads: Number of threads for annotation
            max_indel_length: Maximum indel length (default: 15)
            min_indel_length: Minimum indel length (default: 1)
            synonymous: Whether to generate synonymous SNVs (default: True)
            include_stop_codon: Whether to include stop-codon variants (default: True)
            max_splice_offset: Maximum intron offset for splice sites (default: 20)
            min_splice_offset: Minimum intron offset for splice sites (default: 1)
            include_classic_splice_sites: Include classic splice sites (±1/±2) (default: True)
            max_variants: Maximum total variants to generate (default: no limit)
        """
        self.gtf_file = self._validate_file(gtf_file, "GTF")
        self.fasta_file = self._validate_file(fasta_file, "FASTA")
        self.gene_name = gene_name
        self.transcript_id = transcript_id
        self.database_dir = self._validate_dir(database_dir, "database")
        self.output_dir = self._ensure_dir(output_dir)

        self.protocols = protocols or ['refGene']
        self.operations = operations or ['g']
        self.variant_types = variant_types or ['SNV', 'insertion', 'deletion']
        self.buildver = buildver
        self.threads = threads
        self.max_indel_length = max_indel_length
        self.min_indel_length = min_indel_length
        self.synonymous = synonymous
        self.include_stop_codon = include_stop_codon
        self.max_splice_offset = max_splice_offset
        self.min_splice_offset = min_splice_offset
        self.include_classic_splice_sites = include_classic_splice_sites
        self.max_variants = max_variants

        self.transcript: Optional[GeneTranscript] = None
        self.variants: Optional[Dict] = None
        self.simulated_vcf: Optional[str] = None
        self.annotated_tsv: Optional[str] = None

        self._setup_logging()

        logger.info("Initialized MatchingPipeline")
        logger.info(f"  Gene: {gene_name} ({transcript_id})")
        logger.info(f"  Output directory: {output_dir}")

    def _validate_file(self, path: str, file_type: str) -> str:
        if not os.path.exists(path):
            raise FileNotFoundError(f"{file_type} file not found: {path}")
        return os.path.abspath(path)

    def _validate_dir(self, path: str, dir_type: str) -> str:
        if not os.path.isdir(path):
            raise NotADirectoryError(f"{dir_type} directory not found: {path}")
        return os.path.abspath(path)

    def _ensure_dir(self, path: str) -> str:
        Path(path).mkdir(parents=True, exist_ok=True)
        return os.path.abspath(path)

    def _setup_logging(self):
        log_file = os.path.join(self.output_dir, f"pipeline_{datetime.now().strftime('%Y%m%d_%H%M%S')}.log")
        logging.basicConfig(
            level=logging.INFO,
            format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
            handlers=[
                logging.StreamHandler(sys.stdout),
                logging.FileHandler(log_file, encoding='utf-8')
            ]
        )

    def run(self) -> Dict[str, Any]:
        """
        Execute the complete pipeline

        Returns:
            Dictionary containing pipeline results and metrics
        """
        logger.info("=" * 60)
        logger.info("Starting MATCHVAR Pipeline")
        logger.info("=" * 60)

        # Step 1: Variant Simulation
        logger.info("\n[STEP 1/4] Variant Simulation")
        self._run_simulation()

        # Step 2: Annotation
        logger.info("\n[STEP 2/4] Table Annotation")
        self._run_annotation()

        # Step 3: auROC Calculation
        logger.info("\n[STEP 3/4] auROC Evaluation")
        scores = self._calculate_auroc_scores()

        # Step 4: Visualization
        logger.info("\n[STEP 4/4] Visualization")
        figures = self._generate_visualizations(scores)

        results = {
            'gene_name': self.gene_name,
            'transcript_id': self.transcript_id,
            'simulated_vcf': self.simulated_vcf,
            'annotated_tsv': self.annotated_tsv,
            'total_variants': len(self.variants) if self.variants else 0,
            'auroc_scores': scores,
            'figures': figures,
            'output_dir': self.output_dir
        }

        self._save_summary(results)
        logger.info("\n" + "=" * 60)
        logger.info("Pipeline completed successfully!")
        logger.info("=" * 60)

        return results

    def _run_simulation(self):
        """Run variant simulation"""
        self.simulated_vcf = os.path.join(self.output_dir, f"{self.gene_name}_simulated.vcf")

        # Run simulation and get both transcript and variants
        self.transcript, variants_dict = simulate_variants(
            gtf_file=self.gtf_file,
            fasta_file=self.fasta_file,
            gene_name=self.gene_name,
            transcript_id=self.transcript_id,
            output_vcf=self.simulated_vcf,
            variant_types=self.variant_types,
            max_indel_length=self.max_indel_length,
            min_indel_length=self.min_indel_length,
            synonymous=self.synonymous,
            include_stop_codon=self.include_stop_codon,
            max_splice_offset=self.max_splice_offset,
            min_splice_offset=self.min_splice_offset,
            include_classic_splice_sites=self.include_classic_splice_sites,
            max_variants=self.max_variants,
        )

        # Use the variants directly from simulation (no need to regenerate)
        self.variants = variants_dict
        # Extract variant lists (excluding the 'total' key if present)
        variant_lists = {k: v for k, v in variants_dict.items() if k != 'total'}
        total = sum(len(v) for v in variant_lists.values()) if isinstance(variant_lists, dict) else 0
        logger.info(f"Generated {total} variants")

    def _run_annotation(self):
        """Run table annotation"""
        # The TableAnnotator will create a file with genome version suffix
        self.annotated_tsv = os.path.join(self.output_dir, f"{self.gene_name}_annotated.tsv")
        actual_output = os.path.join(self.output_dir, f"{self.gene_name}_annotated.{self.buildver}_multianno.tsv")

        annotator = TableAnnotator(
            queryfile=self.simulated_vcf,
            dbloc=self.database_dir,
            outfile=os.path.join(self.output_dir, f"{self.gene_name}_annotated"),
            buildver=self.buildver,
            protocol=','.join(self.protocols),
            operation=','.join(self.operations),
            thread=self.threads,
            vcfinput=True,
            otherinfo=True,
            remove=True
        )

        annotator.run_annotation()

        # Check if the expected output file exists
        if not os.path.exists(actual_output):
            raise FileNotFoundError(f"Annotation output not found: {actual_output}")
        
        # Update the annotated_tsv to point to the actual output file
        self.annotated_tsv = actual_output
        logger.info(f"Annotation completed: {self.annotated_tsv}")

    def _calculate_auroc_scores(self) -> Dict[str, Dict[str, Any]]:
        """
        Calculate auROC / auPRC scores for every score column found in the
        annotated TSV.

        This method is intentionally lenient:
        - does NOT require Total_Score to be present
        - does NOT silently drop columns with < 10 valid scores (only warns)
        - still returns an empty dict when appropriate
        """
        if not self.annotated_tsv or not os.path.exists(self.annotated_tsv):
            raise FileNotFoundError(f"Annotated file not found: {self.annotated_tsv}")

        df = pd.read_csv(self.annotated_tsv, sep='\t', low_memory=False)

        # ── ground-truth labels ────────────────────────────────────────────────
        y_true = self._get_true_labels(df, vcf_path=self.simulated_vcf)
        if y_true is None:
            return {}
        if len(y_true) == 0:
            logger.warning("Label array is empty – skipping auROC calculation")
            return {}

        n_pos = int(y_true.sum())
        logger.info(
            f"Ground truth: {n_pos} / {len(y_true)} = "
            f"{n_pos / max(len(y_true), 1):.1%} positive")

        if n_pos == 0:
            logger.warning(
                "All labels are negative (0 positives). "
                "auROC will be NaN – consider checking your label extraction.")
        elif n_pos == len(y_true):
            logger.warning("All labels are positive – degenerate case.")

        # ── score columns ─────────────────────────────────────────────────────
        from .visualization import extract_score_columns

        score_columns = extract_score_columns(df)
        logger.info(f"Score columns found: {score_columns}")

        if not score_columns:
            logger.warning("No numeric score columns found in annotated TSV")
            return {}

        # ── per-column metric computation ─────────────────────────────────────
        from sklearn.metrics import average_precision_score, roc_auc_score, roc_curve

        scores: Dict[str, Dict[str, Any]] = {}

        for col in score_columns:
            y_scores = pd.to_numeric(df[col], errors='coerce').values
            valid_mask = ~np.isnan(y_scores)
            n_valid = int(valid_mask.sum())

            if n_valid < 2:
                logger.warning(f"  '{col}': only {n_valid} valid score(s) – skipped")
                continue

            if n_valid < 10:
                logger.warning(
                    f"  '{col}': {n_valid} valid score(s) < 10 – "
                    f"computing but results may be unreliable")

            try:
                try:
                    auroc = roc_auc_score(y_true[valid_mask], y_scores[valid_mask])
                except ValueError:
                    auroc = float('nan')

                try:
                    auprc = average_precision_score(
                        y_true[valid_mask], y_scores[valid_mask])
                except ValueError:
                    auprc = float('nan')

                fpr_arr, tpr_arr, _ = roc_curve(
                    y_true[valid_mask], y_scores[valid_mask])

                scores[col] = {
                    'auroc': round(float(auroc), 4)
                    if not np.isnan(float(auroc)) else float('nan'),
                    'auprc': round(float(auprc), 4)
                    if not np.isnan(float(auprc)) else float('nan'),
                    'roc_curve': {
                        'fpr': fpr_arr.tolist(),
                        'tpr': tpr_arr.tolist(),
                    },
                    'n_variants': n_valid,
                }
                logger.info(
                    f"  {col}: AUROC={scores[col]['auroc']:.4f}  "
                    f"AUPRC={scores[col]['auprc']:.4f}  (n={n_valid})")

            except Exception as exc:
                logger.warning(f"  Failed to calculate metrics for '{col}': {exc}")

        return scores

    def _get_true_labels(self, df: pd.DataFrame,
                         vcf_path: Optional[str] = None) -> Optional[np.ndarray]:
        """
        Extract binary ground-truth labels from annotated dataframe.

        Tries, in priority order:
        1. Func.refGene / ExonicFunc.refGene (from gene annotation operations)
        2. Otherinfo VCF-INFO columns (--vcfinput --otherinfo path)
        3. VCF file directly (most reliable for --operations f runs)
        4. TYPE column (already-parsed VCF label column)
        """
        import re as _re
        _TYPE_RE  = _re.compile(r'(?:^|;)TYPE=([^;]+)')
        _FRM_RE   = _re.compile(r'(?:^|;)FRAMESHIFT=([^;]+)')
        _PATHO_T  = {'SPLICING', 'SPLICE_SITE', 'FRAMESHIFT'}

        def _lbl_from_info_str(s: str) -> int:
            vt = _TYPE_RE.search(s)
            fr = _FRM_RE.search(s)
            vtt = vt.group(1).upper().strip()  if vt else ''
            frt = fr.group(1).lower().strip()  if fr else 'false'
            return int(vtt in _PATHO_T or frt == 'true')

        try:
            # 1. Func.refGene
            if 'Func.refGene' in df.columns:
                func = df['Func.refGene'].fillna('').str.lower()
                y = func.str.contains('splicing|stopgain|stoploss|frameshift|nonsyn',
                                       na=False).astype(int).values
                logger.info(f"Extracted {y.sum()} positive labels from Func.refGene")
                return y

            # 2. ExonicFunc.refGene
            if 'ExonicFunc.refGene' in df.columns:
                exo = df['ExonicFunc.refGene'].fillna('').str.lower()
                y = exo.str.contains('splicing|stopgain|stoploss|frameshift|nonsynonymous',
                                     na=False).astype(int).values
                logger.info(f"Extracted {y.sum()} positive labels from ExonicFunc.refGene")
                return y

            # 3. TYPE column
            if 'TYPE' in df.columns:
                var_type = df['TYPE'].fillna('').str.upper()
                y = var_type.isin(['SPLICE_SITE', 'FRAMESHIFT', 'NONSENSE', 'STOPLOSS']).astype(int).values
                logger.info(f"Extracted {y.sum()} positive labels from TYPE column")
                return y

            # 4. Otherinfo columns — VCF INFO content (TYPE=…;FRAMESHIFT=…)
            # The pipeline passes -includeinfo -withfreq; the VCF INFO field is typically
            # in the last Otherinfo column.  We parse every Otherinfo column and pick
            # the one that actually contains pathogenic labels.
            for col in df.columns:
                if not col.startswith('Otherinfo'):
                    continue
                raw    = df[col].fillna('').astype(str)
                labels = [_lbl_from_info_str(v) for v in raw]
                if sum(labels) == 0:
                    continue   # not an INFO column — skip
                y = np.array(labels, dtype=np.int8)
                logger.info(
                    f"Extracted {y.sum()} positive labels from VCF INFO in "
                    f"column '{col}' ({y.sum()}/{len(y)} "
                    f"= {y.sum()/max(len(y),1):.1%})")
                return y

            # 5. VCF file directly
            if vcf_path and os.path.exists(vcf_path):
                y = _labels_from_vcf_info(vcf_path)
                if y is not None and len(y) > 0:
                    logger.info(
                        f"Extracted {y.sum()} positive labels from VCF INFO "
                        f"({y.sum()}/{len(y)} = "
                        f"{y.sum()/max(len(y),1):.1%})")
                    return y

            logger.error(
                "No label source found. Tried: Func.refGene / ExonicFunc.refGene / "
                "TYPE / Otherinfo (VCF INFO) / VCF file.")
            return None

        except Exception as exc:
            logger.error(f"Error extracting true labels: {exc}")
            return None


    def _generate_visualizations(self, scores: Dict[str, Dict]) -> Dict[str, str]:
        """
        Generate ROC curves, auROC comparison bar plots, and (when scores are
        empty) diagnostic / exploratory figures so the pipeline always produces
        visual output rather than silently returning nothing.

        Args:
            scores: Dictionary of metric scores per tool

        Returns:
            Dictionary mapping figure keys to file paths
        """
        from .visualization import (create_summary_figure,
                                    create_diagnostic_figures)

        fig_dir = os.path.join(self.output_dir, 'figures')
        Path(fig_dir).mkdir(parents=True, exist_ok=True)

        if scores:
            figures = create_summary_figure(
                tool_metrics=scores,
                gene_name=self.gene_name,
                output_dir=fig_dir,
                prefix=f"{self.gene_name}_performance",
            )
            logger.info(f"Generated {len(figures)} performance figures")
            return figures

        # ── No metrics – produce diagnostic figures instead of returning empty ──
        logger.info("No auROC scores available; producing diagnostic figures")
        tsv_path  = self.annotated_tsv
        vcf_path  = self.simulated_vcf  # set during _run_simulation

        diag = create_diagnostic_figures(
            annotated_tsv=tsv_path,
            gene_name=self.gene_name,
            output_dir=fig_dir,
            vcf_path=vcf_path,
        )
        logger.info(f"Generated {len(diag['figures'])} diagnostic figures")
        return {k: v for k, v in diag.get('figures', {}).items()}

    def _save_summary(self, results: Dict[str, Any]):
        """Save pipeline summary as JSON"""
        summary_file = os.path.join(self.output_dir, f"{self.gene_name}_pipeline_summary.json")

        summary = {
            'gene_name': results['gene_name'],
            'transcript_id': results['transcript_id'],
            'timestamp': datetime.now().isoformat(),
            'total_variants': results['total_variants'],
            'simulated_vcf': results['simulated_vcf'],
            'annotated_tsv': results['annotated_tsv'],
            'auroc_scores': results['auroc_scores'],
            'figures': {k: v for k, v in results['figures'].items() if k.endswith('_pdf')}
        }

        with open(summary_file, 'w', encoding='utf-8') as f:
            json.dump(summary, f, indent=2)

        logger.info(f"Summary saved: {summary_file}")

        # Also create a TSV statistics table
        stats_file = os.path.join(self.output_dir, f"{self.gene_name}_auroc_statistics.tsv")
        if results['auroc_scores']:
            stats_df = pd.DataFrame([
                {
                    'Tool': tool,
                    'AUROC': metrics['auroc'],
                    'AUPRC': metrics['auprc'],
                    'N_Variants': metrics['n_variants']
                }
                for tool, metrics in results['auroc_scores'].items()
            ])
            stats_df.to_csv(stats_file, sep='\t', index=False)
            logger.info(f"Statistics table saved: {stats_file}")


def _labels_from_vcf_info(vcf_path: str) -> Optional[np.ndarray]:
    """
    Parse VCF INFO and return per-position binary labels.
    Each VCF record contributes one label: 1 if TYPE is SPLICING/SPLICE_SITE/
    FRAMESHIFT or if FRAMESHIFT=true, 0 otherwise.
    When a genomic position has multiple records the label is the OR of all labels.
    """
    import re as _re
    _TYPE_RE = _re.compile(r'(?:^|;)TYPE=([^;]+)')
    _FRM_RE  = _re.compile(r'(?:^|;)FRAMESHIFT=([^;]+)')
    _PATHO   = {'SPLICING', 'SPLICE_SITE', 'FRAMESHIFT'}

    if not os.path.exists(vcf_path):
        logger.error(f"VCF file not found: {vcf_path}")
        return None

    pos_labels: Dict[str, int] = {}
    pos_order: List[str] = []

    try:
        with open(vcf_path, 'r') as fh:
            for line in fh:
                if line.startswith('#'):
                    continue
                parts = line.strip().split('\t')
                if len(parts) < 8:
                    continue
                info  = parts[7]
                vt_m  = _TYPE_RE.search(info)
                fr_m  = _FRM_RE.search(info)
                vt    = vt_m.group(1).upper().strip() if vt_m else ''
                fr    = fr_m.group(1).lower().strip() if fr_m else 'false'
                lbl   = int(vt in _PATHO or fr == 'true')
                if parts[1] not in pos_labels:
                    pos_labels[parts[1]] = lbl
                    pos_order.append(parts[1])
                elif lbl == 1:
                    pos_labels[parts[1]] = 1

        if not pos_order:
            logger.warning("No VCF records found while extracting labels")
            return None

        aligned = [pos_labels[p] for p in pos_order]
        logger.info(
            f"VCF direct labels: {sum(aligned)}/{len(aligned)} positive "
            f"({sum(aligned)/max(len(aligned),1):.1%}) across {len(pos_order)} positions")
        return np.array(aligned, dtype=np.int8)

    except Exception as exc:
        logger.error(f"Error reading VCF file for labels: {exc}")
        return None

    def _save_summary(self, results: Dict[str, Any]):
        """Save pipeline summary as JSON"""
        summary_file = os.path.join(self.output_dir, f"{self.gene_name}_pipeline_summary.json")

        summary = {
            'gene_name': results['gene_name'],
            'transcript_id': results['transcript_id'],
            'timestamp': datetime.now().isoformat(),
            'total_variants': results['total_variants'],
            'simulated_vcf': results['simulated_vcf'],
            'annotated_tsv': results['annotated_tsv'],
            'auroc_scores': results['auroc_scores'],
            'figures': {k: v for k, v in results['figures'].items() if k.endswith('_pdf')}
        }

        with open(summary_file, 'w', encoding='utf-8') as f:
            json.dump(summary, f, indent=2)

        logger.info(f"Summary saved: {summary_file}")

        # Also create a TSV statistics table
        stats_file = os.path.join(self.output_dir, f"{self.gene_name}_auroc_statistics.tsv")
        if results['auroc_scores']:
            stats_df = pd.DataFrame([
                {
                    'Tool': tool,
                    'AUROC': metrics['auroc'],
                    'AUPRC': metrics['auprc'],
                    'N_Variants': metrics['n_variants']
                }
                for tool, metrics in results['auroc_scores'].items()
            ])
            stats_df.to_csv(stats_file, sep='\t', index=False)
            logger.info(f"Statistics table saved: {stats_file}")


def run_pipeline(gtf_file: str, fasta_file: str, gene_name: str, transcript_id: str,
                 database_dir: str, output_dir: str, **kwargs) -> Dict[str, Any]:
    """
    Run pipeline directly from Python code

    Args:
        gtf_file: Path to GTF file
        fasta_file: Path to reference genome FASTA
        gene_name: Gene symbol
        transcript_id: Transcript ID
        database_dir: Database directory
        output_dir: Output directory
        **kwargs: Additional arguments for MatchingPipeline

    Returns:
        Dictionary with pipeline results
    """
    pipeline = MatchingPipeline(
        gtf_file=gtf_file,
        fasta_file=fasta_file,
        gene_name=gene_name,
        transcript_id=transcript_id,
        database_dir=database_dir,
        output_dir=output_dir,
        **kwargs
    )
    return pipeline.run()


def run_pipeline_from_args(args):
    """Run pipeline from parsed command-line arguments"""
    return run_pipeline(
        gtf_file=args.gtf,
        fasta_file=args.fasta,
        gene_name=args.gene,
        transcript_id=args.transcript,
        database_dir=args.database,
        output_dir=args.output_dir,
        protocols=args.protocol.split(',') if args.protocol else None,
        operations=args.operation.split(',') if args.operation else None,
        variant_types=args.variant_types.split(',') if args.variant_types else None,
        buildver=args.buildver,
        threads=args.threads
    )