#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
MATCHVAR Pipeline Module

End-to-end pipeline for variant simulation, annotation, and evaluation:
1. Simulate variants from GTF + gene/transcript -> VCF with biological scores
2. Annotate simulated VCF with external prediction tools
3. Calculate auROC scores comparing predictions vs simulated truth
4. Generate publication-quality visualizations
"""

import os
import sys
import logging
import argparse
from pathlib import Path
from types import SimpleNamespace
from typing import Dict, List, Optional, Tuple, Any
from datetime import datetime
import json

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from .variant_simulation import GeneTranscript, ExonExtractor
from .table_matchvar import TableAnnotator

logger = logging.getLogger(__name__)

# VCF INFO thresholds used by _lbl_from_info_str
_PATHO_VT = frozenset({'SPLICING', 'SPLICE_SITE', 'FRAMESHIFT'})


# ══════════════════════════════════════════════════════════════════════════════
# Standalone utilities
# ══════════════════════════════════════════════════════════════════════════════

def simulate_variants(
    gtf_file: str,
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
    max_variants: Optional[int] = None,
) -> Tuple[GeneTranscript, Dict]:
    """
    Extract exon data from GTF/FASTA, generate variants, and export to VCF.

    Args:
        gtf_file: Path to GTF annotation file
        fasta_file: Path to reference genome FASTA
        gene_name: Target gene symbol
        transcript_id: Target transcript ID
        output_vcf: Path for output VCF file
        variant_types: Variant types to simulate (default: ['SNV', 'insertion', 'deletion'])
        max_indel_length: Maximum indel length in bp (default: 5)
        min_indel_length: Minimum indel length in bp (default: 1)
        synonymous: Whether to generate synonymous SNVs (default: True)
        include_stop_codon: Whether to include stop-codon variants (default: True)
        max_splice_offset: Maximum intron offset for splice site variants (default: 20)
        min_splice_offset: Minimum intron offset for splice site variants (default: 1)
        include_classic_splice_sites: Include classic splice sites (+/-1, +/-2) (default: True)
        max_variants: Maximum total variants to generate (default: no limit)

    Returns:
        Tuple of (GeneTranscript instance, variants dict).
    """
    extractor = ExonExtractor(gtf_file, fasta_file)
    exons, chromosome, strand = extractor.extract_exons(gene_name, transcript_id)

    transcript = GeneTranscript(
        gene_name=gene_name,
        transcript_id=transcript_id,
        exons=exons,
        chromosome=chromosome,
        strand=strand,
        genome=extractor.genome,
        fasta_file=fasta_file,
    )

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

    transcript.export_to_vcf(variants, output_vcf)
    return transcript, variants


def _labels_from_vcf_info(vcf_path: str) -> Optional[np.ndarray]:
    """
    Parse VCF INFO column and return per-position binary labels.

    Each record is labelled 1 if TYPE is SPLICING/SPLICE_SITE/FRAMESHIFT or
    if FRAMESHIFT=true; 0 otherwise.  When multiple records share a position
    the label is the logical-OR of all record labels.
    """
    import re as _re
    _TYPE_RE = _re.compile(r'(?:^|;)TYPE=([^;]+)')
    _FRM_RE  = _re.compile(r'(?:^|;)FRAMESHIFT=([^;]+)')

    if not os.path.exists(vcf_path):
        logger.error(f" VCF file not found: {vcf_path}")
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
                lbl   = int(vt in _PATHO_VT or fr == 'true')
                if parts[1] not in pos_labels:
                    pos_labels[parts[1]] = lbl
                    pos_order.append(parts[1])
                elif lbl == 1:
                    pos_labels[parts[1]] = 1

        if not pos_order:
            logger.warning(" No VCF records found while extracting labels")
            return None

        aligned = [pos_labels[p] for p in pos_order]
        logger.info(
            f" VCF direct labels: {sum(aligned)}/{len(aligned)} positive "
            f"({sum(aligned) / max(len(aligned), 1):.1%}) across {len(pos_order)} positions")
        return np.array(aligned, dtype=np.int8)

    except Exception as exc:
        logger.error(f" Error reading VCF file for labels: {exc}")
        return None


def _run_annotation_for_vcf(
    simulated_vcf: str,
    output_dir: str,
    gene_name: str,
    database_dir: str,
    buildver: str,
    protocols: List[str],
    operations: List[str],
    threads: int,
) -> str:
    """
    Standalone annotation helper: annotate a simulated VCF and return the
    path to the resulting multi-anno TSV.

    Used by both MatchingPipeline._run_annotation() and the merge-output
    path in run_pipeline() so that no dummy MatchingPipeline object is needed.
    """
    actual_output = os.path.join(
        output_dir, f"{gene_name}_annotated.{buildver}_multianno.tsv")

    annotator = TableAnnotator(
        queryfile=simulated_vcf,
        dbloc=database_dir,
        outfile=os.path.join(output_dir, f"{gene_name}_annotated"),
        buildver=buildver,
        protocol=','.join(protocols),
        operation=','.join(operations),
        thread=threads,
        vcfinput=True,
        otherinfo=True,
        remove=True,
    )

    annotator.run_annotation()

    if not os.path.exists(actual_output):
        raise FileNotFoundError(f"Annotation output not found: {actual_output}")

    logger.info("=" * 70)
    logger.info(f" Annotation completed: {actual_output}")
    logger.info("=" * 70)
    return actual_output


# ══════════════════════════════════════════════════════════════════════════════
# Label extraction (module-level, shared by pipeline and visualization)
# ══════════════════════════════════════════════════════════════════════════════

def _get_true_labels(
    df: pd.DataFrame,
    vcf_path: Optional[str] = None,
) -> Optional[np.ndarray]:
    """
    Extract binary ground-truth labels from an annotated DataFrame.

    Priority order:
    1. Func.refGene / ExonicFunc.refGene
    2. Otherinfo VCF-INFO columns (VCF INFO embedded in TSV text)
    3. VCF file directly (most reliable for --operations f runs)
    4. TYPE column

    Returns:
        Integer ndarray of 0/1 labels, or None when no suitable source is found.
    """
    import re as _re
    _TYPE_RE = _re.compile(r'(?:^|;)TYPE=([^;]+)')
    _FRM_RE  = _re.compile(r'(?:^|;)FRAMESHIFT=([^;]+)')

    def _lbl_from_info_str(s: str) -> int:
        vt_m = _TYPE_RE.search(s)
        fr_m = _FRM_RE.search(s)
        vt   = vt_m.group(1).upper().strip() if vt_m else ''
        fr   = fr_m.group(1).lower().strip() if fr_m else 'false'
        return int(vt in _PATHO_VT or fr == 'true')

    # 1. Func.refGene
    if 'Func.refGene' in df.columns:
        func = df['Func.refGene'].fillna('').str.lower()
        y    = func.str.contains(
            'splicing|stopgain|stoploss|frameshift|nonsyn', na=False
        ).astype(np.int8).values
        logger.info(f" Extracted {y.sum()} positive labels from Func.refGene")
        return y

    # 2. ExonicFunc.refGene
    if 'ExonicFunc.refGene' in df.columns:
        exo  = df['ExonicFunc.refGene'].fillna('').str.lower()
        y    = exo.str.contains(
            'splicing|stopgain|stoploss|frameshift|nonsynonymous', na=False
        ).astype(np.int8).values
        logger.info(f" Extracted {y.sum()} positive labels from ExonicFunc.refGene")
        return y

    # 3. TYPE column
    if 'TYPE' in df.columns:
        var_type  = df['TYPE'].fillna('').str.upper()
        y         = var_type.isin(
            ['SPLICE_SITE', 'FRAMESHIFT', 'NONSENSE', 'STOPLOSS']
        ).astype(np.int8).values
        logger.info(f" Extracted {y.sum()} positive labels from TYPE column")
        return y

    # 4. Otherinfo columns – parse VCF INFO from embedded text
    for col in df.columns:
        if not col.startswith('Otherinfo'):
            continue
        raw    = df[col].fillna('').astype(str)
        labels = [_lbl_from_info_str(v) for v in raw]
        if sum(labels) == 0:
            continue   # not an INFO column – try next Otherinfo
        y = np.array(labels, dtype=np.int8)
        logger.info(
            f" Extracted {y.sum()} positive labels from VCF INFO in "
            f"column '{col}' ({y.sum()}/{len(y)} = "
            f"{y.sum() / max(len(y), 1):.1%})")
        return y

    # 5. VCF file directly
    if vcf_path and os.path.exists(vcf_path):
        y = _labels_from_vcf_info(vcf_path)
        if y is not None and len(y) > 0:
            logger.info(
                f" Extracted {y.sum()} positive labels from VCF INFO "
                f"({y.sum()}/{len(y)} = {y.sum() / max(len(y), 1):.1%})")
            return y

    logger.error(
        " No label source found. Tried: Func.refGene / ExonicFunc.refGene / "
        "TYPE / Otherinfo (VCF INFO) / VCF file.")
    return None


# ══════════════════════════════════════════════════════════════════════════════
# Module-level helpers – counting and summary
# ══════════════════════════════════════════════════════════════════════════════

def _count_total_variants(variants: Optional[Dict]) -> int:
    """Return total number of simulated variants excluding any 'total' key."""
    if not variants:
        return 0
    return sum(len(v) for k, v in variants.items()
               if k != 'total' and isinstance(v, list))


def _save_summary(results: Dict[str, Any]) -> None:
    """Save pipeline summary as JSON + AUROC statistics TSV."""
    output_dir   = results['output_dir']
    gene_name    = results['gene_name']
    summary_file = os.path.join(output_dir, f"{gene_name}_pipeline_summary.json")

    summary = {
        'gene_name':      results['gene_name'],
        'transcript_id':  results['transcript_id'],
        'timestamp':      datetime.now().isoformat(),
        'total_variants': results['total_variants'],
        'simulated_vcf':  results['simulated_vcf'],
        'annotated_tsv':  results['annotated_tsv'],
        'auroc_scores':   results['auroc_scores'],
        'figures': {k: v for k, v in results['figures'].items() if k.endswith('_pdf')},
    }
    with open(summary_file, 'w', encoding='utf-8') as f:
        json.dump(summary, f, indent=2)
    logger.info(f" Summary saved: {summary_file}")

    stats_file = os.path.join(output_dir, f"{gene_name}_auroc_statistics.tsv")
    if results['auroc_scores']:
        stats_df = pd.DataFrame([
            {
                'Tool':        tool,
                'AUROC':       metrics['auroc'],
                'AUPRC':       metrics['auprc'],
                'N_Variants':  metrics['n_variants'],
            }
            for tool, metrics in results['auroc_scores'].items()
        ])
        stats_df.to_csv(stats_file, sep='\t', index=False)
        logger.info(f" Statistics table saved: {stats_file}")


def _save_summary_for_result(
    results: Dict[str, Any],
    output_dir: str,
    gene_name: str,
) -> None:
    """Save a pipeline results dictionary as JSON + AUROC statistics TSV (merged-output path)."""
    summary_file = os.path.join(output_dir, f"{gene_name}_pipeline_summary.json")
    summary = {
        'gene_name':      results.get('gene_name', gene_name),
        'transcript_id':  results.get('transcript_id', ''),
        'timestamp':      datetime.now().isoformat(),
        'total_variants': results.get('total_variants', 0),
        'simulated_vcf':  results.get('simulated_vcf', ''),
        'annotated_tsv':  results.get('annotated_tsv', ''),
        'auroc_scores':   results.get('auroc_scores', {}),
        'figures': {k: v for k, v in results.get('figures', {}).items() if k.endswith('_pdf')},
    }
    with open(summary_file, 'w', encoding='utf-8') as f:
        json.dump(summary, f, indent=2)
    logger.info(f" Summary saved: {summary_file}")

    stats_file = os.path.join(output_dir, f"{gene_name}_auroc_statistics.tsv")
    if results.get('auroc_scores'):
        stats_df = pd.DataFrame([
            {
                'Tool':        tool,
                'AUROC':       metrics['auroc'],
                'AUPRC':       metrics['auprc'],
                'N_Variants':  metrics['n_variants'],
            }
            for tool, metrics in results['auroc_scores'].items()
        ])
        stats_df.to_csv(stats_file, sep='\t', index=False)
        logger.info(f" Statistics table saved: {stats_file}")


def _generate_visualizations_for_scores(
    scores: Dict[str, Any],
    output_dir: str,
    gene_name: str,
) -> Dict[str, str]:
    """Generate summary figures from pre-computed auROC scores (merged-output path)."""
    from .visualization import create_summary_figure, create_diagnostic_figures
    fig_dir = os.path.join(output_dir, 'figures')
    Path(fig_dir).mkdir(parents=True, exist_ok=True)
    if scores:
        return create_summary_figure(
            tool_metrics=scores,
            gene_name=gene_name,
            output_dir=fig_dir,
            prefix=f"{gene_name}_performance",
        )
    return {}


def _calculate_auroc_scores_core(
    annotated_tsv: Optional[str],
    vcf_path: Optional[str],
) -> Dict[str, Dict[str, Any]]:
    """
    Standalone auROC / auPRC computation — the heart of both
    MatchingPipeline._calculate_auroc_scores() and
    _calculate_auroc_scores_from_tsv().

    Args:
        annotated_tsv: Path to the annotated TSV file.
        vcf_path:      Path to the companion VCF (for label fallback).

    Returns:
        ``{col: {auroc, auprc, roc_curve, n_variants}}`` per score column.
    """
    if not annotated_tsv or not os.path.exists(annotated_tsv):
        raise FileNotFoundError(f"Annotated file not found: {annotated_tsv}")

    df = pd.read_csv(annotated_tsv, sep='\t', low_memory=False)

    # ground-truth labels
    y_true = _get_true_labels(df, vcf_path=vcf_path)
    if y_true is None:
        return {}
    if len(y_true) == 0:
        logger.warning(" Label array is empty – skipping auROC calculation")
        return {}

    n_pos = int(y_true.sum())
    logger.info(
        f" Ground truth: {n_pos} / {len(y_true)} = "
        f"{n_pos / max(len(y_true), 1):.1%} positive")

    if n_pos == 0:
        logger.warning(
            " All labels are negative (0 positives). "
            " auROC will be NaN – consider checking your label extraction.")
    elif n_pos == len(y_true):
        logger.warning(" All labels are positive – degenerate case.")

    # score columns
    from .visualization import extract_score_columns
    score_columns = extract_score_columns(df)
    logger.info(f" Score columns found: {score_columns}")

    if not score_columns:
        logger.warning(" No numeric score columns found in annotated TSV")
        return {}

    # per-column metric computation
    from sklearn.metrics import average_precision_score, roc_auc_score, roc_curve
    scores: Dict[str, Dict[str, Any]] = {}

    for col in score_columns:
        y_scores   = np.asarray(pd.to_numeric(df[col], errors='coerce'))
        valid_mask = ~np.isnan(y_scores)
        n_valid    = int(valid_mask.sum())

        if n_valid < 2:
            logger.warning(f" '{col}': only {n_valid} valid score(s) – skipped")
            continue

        if n_valid < 10:
            logger.warning(
                f" '{col}': {n_valid} valid score(s) < 10 – "
                f"computing but results may be unreliable")

        try:
            try:
                auroc = roc_auc_score(y_true[valid_mask], y_scores[valid_mask])
            except ValueError:
                auroc = float('nan')

            try:
                auprc = average_precision_score(y_true[valid_mask], y_scores[valid_mask])
            except ValueError:
                auprc = float('nan')

            fpr_arr, tpr_arr, _ = roc_curve(y_true[valid_mask], y_scores[valid_mask])

            scores[col] = {
                'auroc':      round(float(auroc), 4)
                              if not np.isnan(float(auroc)) else float('nan'),
                'auprc':      round(float(auprc), 4)
                              if not np.isnan(float(auprc)) else float('nan'),
                'roc_curve': {
                    'fpr': fpr_arr.tolist(),
                    'tpr': tpr_arr.tolist(),
                },
                'n_variants': n_valid,
            }
            logger.info(
                f" {col}: AUROC={scores[col]['auroc']:.4f}  "
                f"AUPRC={scores[col]['auprc']:.4f}  (n={n_valid})")

        except Exception as exc:
            logger.warning(f" Failed to calculate metrics for '{col}': {exc}")

    return scores


def _calculate_auroc_scores_from_tsv(
    annotated_tsv: str,
    vcf_path: str,
) -> Dict[str, Dict[str, Any]]:
    """Standalone auROC computation from file paths (merged-output path)."""
    return _calculate_auroc_scores_core(
        annotated_tsv=annotated_tsv,
        vcf_path=vcf_path,
    )


# ══════════════════════════════════════════════════════════════════════════════
# MatchingPipeline
# ══════════════════════════════════════════════════════════════════════════════

class MatchingPipeline:
    """
    End-to-end pipeline for variant simulation and annotation evaluation.

    Workflow:
        GTF + Gene/Transcript -> Variant Simulation -> Table Annotation ->
        auROC Evaluation -> Visualization
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
                 max_variants: Optional[int] = None,
                 simulation_strategy: str = 'traditional',
                 min_probability: float = 0.001,
                 constraint_filter: bool = False,
                 spectrum_ratio: float = 0.7,
                 score_correction_method: str = 'none',
                 score_correction_clusters: int = 3):
        """
        Initialize the matching pipeline.

        Args:
            gtf_file:   Path to GTF annotation file.
            fasta_file: Path to reference genome FASTA.
            gene_name:  Target gene symbol.
            transcript_id: Target transcript ID.
            database_dir: Directory containing annotation databases.
            output_dir:  Directory for all output files.
            protocols:  Annotation protocols (default: ['refGene']).
            operations: Operations list (default: ['g']).
            variant_types: Variant types to simulate (default: ['SNV','insertion','deletion']).
            buildver:   Genome build version (default: hg19).
            threads:    Threads for annotation (default: 4).
            max_indel_length: Maximum indel length in bp (default: 6).
            min_indel_length: Minimum indel length in bp (default: 1).
            synonymous: Include synonymous SNVs (default: True).
            include_stop_codon: Include stop-codon variants (default: True).
            max_splice_offset: Max intron offset for splice-site variants (default: 20).
            min_splice_offset: Min intron offset for splice-site variants (default: 1).
            include_classic_splice_sites: Include +/-1, +/-2 splice sites (default: True).
            max_variants: Maximum total simulated variants, None=unlimited.
            simulation_strategy: Variant simulation strategy ('traditional', 'spectrum', 'hybrid').
            min_probability: Min mutation probability for spectrum simulation.
            constraint_filter: Use gnomAD constraint filtering.
            spectrum_ratio: Ratio of spectrum vs traditional variants for hybrid mode.
            score_correction_method: ClinVar-based score correction method.
            score_correction_clusters: Number of clusters for K-means correction.
        """
        self.gtf_file          = self._validate_file(gtf_file, "GTF")
        self.fasta_file        = self._validate_file(fasta_file, "FASTA")
        self.gene_name         = gene_name
        self.transcript_id     = transcript_id
        self.database_dir      = self._validate_dir(database_dir, "database")
        self.output_dir        = self._ensure_dir(output_dir)

        self.protocols                   = protocols or ['refGene']
        self.operations                  = operations or ['g']
        self.variant_types               = variant_types or ['SNV', 'insertion', 'deletion']
        self.buildver                    = buildver
        self.threads                     = threads
        self.max_indel_length            = max_indel_length
        self.min_indel_length            = min_indel_length
        self.synonymous                  = synonymous
        self.include_stop_codon          = include_stop_codon
        self.max_splice_offset           = max_splice_offset
        self.min_splice_offset           = min_splice_offset
        self.include_classic_splice_sites = include_classic_splice_sites
        self.max_variants                = max_variants
        self.simulation_strategy         = simulation_strategy
        self.min_probability             = min_probability
        self.constraint_filter           = constraint_filter
        self.spectrum_ratio              = spectrum_ratio
        self.score_correction_method     = score_correction_method
        self.score_correction_clusters   = score_correction_clusters

        self.transcript:    Optional[GeneTranscript] = None
        self.variants:      Optional[Dict]           = None
        self.simulated_vcf: Optional[str]            = None
        self.annotated_tsv: Optional[str]            = None

        self._setup_logging()

        logger.info(" Initialized MatchingPipeline")
        logger.info(f" Gene: {gene_name} ({transcript_id})")
        logger.info(f" Output directory: {output_dir}")

    # ── Validation helpers ────────────────────────────────────────────────

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

    def _setup_logging(self) -> None:
        log_file = os.path.join(
            self.output_dir,
            f"pipeline_{datetime.now().strftime('%Y%m%d_%H%M%S')}.log")
        logging.basicConfig(
            level=logging.INFO,
            format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
            handlers=[
                logging.StreamHandler(sys.stdout),
                logging.FileHandler(log_file, encoding='utf-8'),
            ]
        )

    # ── Public API ────────────────────────────────────────────────────────

    def run(self) -> Dict[str, Any]:
        """Execute the complete pipeline. Returns a dict of results and metrics."""
        logger.info("=" * 70)
        logger.info(" Starting MATCHVAR Pipeline")
        logger.info("=" * 70)

        logger.info("==================== [STEP 1/4] Variant Simulation ====================")
        self._run_simulation()

        logger.info("==================== [STEP 2/4] Table Annotation ====================")
        self._run_annotation()

        logger.info("==================== [STEP 3/4] auROC Evaluation ====================")
        scores = self._calculate_auroc_scores()

        logger.info("==================== [STEP 4/4] Visualization ====================")
        figures = self._generate_visualizations(scores)

        total_variants = _count_total_variants(self.variants)
        results = {
            'gene_name':      self.gene_name,
            'transcript_id':  self.transcript_id,
            'simulated_vcf':  self.simulated_vcf,
            'annotated_tsv':  self.annotated_tsv,
            'total_variants': total_variants,
            'auroc_scores':   scores,
            'figures':        figures,
            'output_dir':     self.output_dir,
        }
        _save_summary(results)

        logger.info("=" * 70)
        logger.info(" Pipeline completed successfully!")
        logger.info("=" * 70)
        return results

    # ── Pipeline steps ────────────────────────────────────────────────────

    def _run_simulation(self) -> None:
        """Run variant simulation and write the simulated VCF."""
        self.simulated_vcf = os.path.join(self.output_dir, f"{self.gene_name}_simulated.vcf")

        # Use enhanced simulation for spectrum/hybrid strategies
        if self.simulation_strategy == 'spectrum':
            self._run_simulation_spectrum()
        elif self.simulation_strategy == 'hybrid':
            self._run_simulation_hybrid()
        else:
            self._run_simulation_traditional()

        variant_lists = {k: v for k, v in self.variants.items() if k != 'total'}
        total = sum(len(v) for v in variant_lists.values()) if isinstance(variant_lists, dict) else 0
        logger.info(f" Generated {total} variants")

    def _run_simulation_traditional(self) -> None:
        """Traditional variant simulation."""
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
        self.variants = variants_dict

    def _run_simulation_spectrum(self) -> None:
        """Spectrum-aware variant simulation using utils."""
        import sys
        sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
        from variant_simulation import ExonExtractor
        from utils.enhanced_data_simulation import EnhancedGeneTranscript
        
        try:
            extractor = ExonExtractor(self.gtf_file, self.fasta_file)
            exons, chromosome, strand = extractor.extract_exons(self.gene_name, self.transcript_id)
            
            # Convert exons format to match GeneTranscript expectations
            formatted_exons = []
            for exon in exons:
                formatted_exons.append({
                    'genomic_start': exon['genomic_start'],
                    'genomic_end': exon['genomic_end'],
                    'cds_sequence': exon['cds_sequence'],
                })
            
            enhanced_tx = EnhancedGeneTranscript(
                gene_name=self.gene_name,
                transcript_id=self.transcript_id,
                exons=formatted_exons,
                chromosome=chromosome,
                strand=strand,
                genome=extractor.genome,
                fasta_file=self.fasta_file,
            )
            
            variants_result = enhanced_tx.generate_spectrum_aware_variants(
                max_variants=self.max_variants,
                min_probability=self.min_probability,
                use_constraint_filter=self.constraint_filter,
            )
            variants_list = variants_result.get('variants', [])
            
            # If no variants generated, fallback to traditional
            if not variants_list:
                raise ValueError("No spectrum variants generated")
                
            # Convert list format to dict format expected by export_to_vcf
            variants_dict = {'snvs': [v for v in variants_list if v.get('type') == 'SNV'],
                            'insertions': [v for v in variants_list if v.get('type') == 'insertion'],
                            'deletions': [v for v in variants_list if v.get('type') == 'deletion'],
                            'splice_sites': [v for v in variants_list if v.get('type') == 'splice_site']}
            self.variants = variants_dict
            enhanced_tx.export_to_vcf(variants_dict, self.simulated_vcf)
            self.transcript = enhanced_tx
            
        except Exception as e:
            logger.warning(f"Spectrum simulation failed, falling back to traditional: {e}")
            self._run_simulation_traditional()

    def _run_simulation_hybrid(self) -> None:
        """Hybrid variant simulation using utils."""
        import sys
        sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
        from variant_simulation import ExonExtractor
        from utils.enhanced_data_simulation import EnhancedGeneTranscript
        
        try:
            extractor = ExonExtractor(self.gtf_file, self.fasta_file)
            exons, chromosome, strand = extractor.extract_exons(self.gene_name, self.transcript_id)
            
            # Convert exons format to match GeneTranscript expectations
            formatted_exons = []
            for exon in exons:
                formatted_exons.append({
                    'genomic_start': exon['genomic_start'],
                    'genomic_end': exon['genomic_end'],
                    'cds_sequence': exon['cds_sequence'],
                })
            
            enhanced_tx = EnhancedGeneTranscript(
                gene_name=self.gene_name,
                transcript_id=self.transcript_id,
                exons=formatted_exons,
                chromosome=chromosome,
                strand=strand,
                genome=extractor.genome,
                fasta_file=self.fasta_file,
            )
            
            variants_result = enhanced_tx.generate_hybrid_variants(
                spectrum_ratio=self.spectrum_ratio,
                max_variants=self.max_variants,
                min_probability=self.min_probability,
                use_constraint_filter=self.constraint_filter,
            )
            variants_list = variants_result.get('all_variants', [])
            
            if not variants_list:
                raise ValueError("No hybrid variants generated")
                
            variants_dict = {'snvs': [v for v in variants_list if v.get('type') == 'SNV'],
                            'insertions': [v for v in variants_list if v.get('type') == 'insertion'],
                            'deletions': [v for v in variants_list if v.get('type') == 'deletion'],
                            'splice_sites': [v for v in variants_list if v.get('type') == 'splice_site']}
            self.variants = variants_dict
            enhanced_tx.export_to_vcf(variants_dict, self.simulated_vcf)
            self.transcript = enhanced_tx
            
        except Exception as e:
            logger.warning(f"Hybrid simulation failed, falling back to traditional: {e}")
            self._run_simulation_traditional()

    def _run_annotation(self) -> None:
        """Run table annotation on the simulated VCF."""
        assert self.simulated_vcf, "Call _run_simulation before _run_annotation"
        self.annotated_tsv = _run_annotation_for_vcf(
            simulated_vcf=self.simulated_vcf,
            output_dir=self.output_dir,
            gene_name=self.gene_name,
            database_dir=self.database_dir,
            buildver=self.buildver,
            protocols=self.protocols,
            operations=self.operations,
            threads=self.threads,
        )
        logger.info(f" Annotation completed: {self.annotated_tsv}")

    def _calculate_auroc_scores(self) -> Dict[str, Dict[str, Any]]:
        """Calculate auROC / auPRC for every score column in the annotated TSV."""
        return _calculate_auroc_scores_core(
            annotated_tsv=self.annotated_tsv,
            vcf_path=self.simulated_vcf,
        )

    def _generate_visualizations(
        self,
        scores: Dict[str, Dict],
    ) -> Dict[str, str]:
        """
        Generate ROC / auROC / PR figures when scores are available,
        or diagnostic / exploratory figures otherwise.
        """
        from .visualization import (
            create_summary_figure,
            create_diagnostic_figures,
        )

        fig_dir = os.path.join(self.output_dir, 'figures')
        Path(fig_dir).mkdir(parents=True, exist_ok=True)

        if scores:
            figures = create_summary_figure(
                tool_metrics=scores,
                gene_name=self.gene_name,
                output_dir=fig_dir,
                prefix=f"{self.gene_name}_performance",
            )
            logger.info(f" Generated {len(figures)} performance figures")
            return figures

        logger.info(" No auROC scores available; producing diagnostic figures")
        assert self.annotated_tsv, "Run _run_annotation before _generate_visualizations"
        assert self.simulated_vcf, "Run _run_simulation before _generate_visualizations"
        diag = create_diagnostic_figures(
            annotated_tsv=self.annotated_tsv,
            gene_name=self.gene_name,
            output_dir=fig_dir,
            vcf_path=self.simulated_vcf,
        )
        logger.info(f" Generated {len(diag['figures'])} diagnostic figures")
        return {k: v for k, v in diag.get('figures', {}).items()}

    # ── Label extraction (delegates to module function) ───────────────────

    def _get_true_labels(
        self,
        df: pd.DataFrame,
        vcf_path: Optional[str] = None,
    ) -> Optional[np.ndarray]:
        """Extract binary ground-truth labels from an annotated DataFrame."""
        return _get_true_labels(df, vcf_path)


# ══════════════════════════════════════════════════════════════════════════════
# Multi-gene / merge-output entry point
# ══════════════════════════════════════════════════════════════════════════════

def run_pipeline(gtf_file: str, fasta_file: str, gene_name: str, transcript_id: str,
                  database_dir: str, output_dir: str, merge_output: bool = False,
                  simulation_strategy: str = 'traditional',
                  min_probability: float = 0.001,
                  constraint_filter: bool = False,
                  spectrum_ratio: float = 0.7,
                  score_correction_method: str = 'none',
                  score_correction_clusters: int = 3,
                  **kwargs) -> Dict[str, Any]:
    """
    Run pipeline directly from Python code.

    Args:
        gtf_file:      Path to GTF file.
        fasta_file:    Path to reference genome FASTA.
        gene_name:     Gene symbol(s), comma-separated for multiple genes.
        transcript_id: Transcript ID(s), comma-separated matching genes.
        database_dir:  Database directory.
        output_dir:    Output directory.
        merge_output:  If True, merge all genes into one VCF/annotation output;
                       False = separate files per gene (default).
        simulation_strategy: Variant simulation strategy.
        min_probability: Min mutation probability for spectrum simulation.
        constraint_filter: Use gnomAD constraint filtering.
        spectrum_ratio: Ratio of spectrum vs traditional variants.
        score_correction_method: ClinVar-based score correction method.
        score_correction_clusters: Number of clusters for K-means correction.
        **kwargs:      Additional keyword arguments forwarded to MatchingPipeline.

    Returns:
        Dictionary with pipeline results.
    """
    gene_list = [g.strip() for g in gene_name.split(",") if g.strip()]
    tx_list   = [t.strip() for t in transcript_id.split(",") if t.strip()]

    if len(gene_list) != len(tx_list):
        raise ValueError(
            f"Gene counts ({len(gene_list)}) != Transcript counts ({len(tx_list)})")

    # ── Single gene ──────────────────────────────────────────────────────────
    if len(gene_list) == 1:
        pipeline = MatchingPipeline(
            gtf_file=gtf_file,
            fasta_file=fasta_file,
            gene_name=gene_list[0],
            transcript_id=tx_list[0],
            database_dir=database_dir,
            output_dir=output_dir,
            simulation_strategy=simulation_strategy,
            min_probability=min_probability,
            constraint_filter=constraint_filter,
            spectrum_ratio=spectrum_ratio,
            score_correction_method=score_correction_method,
            score_correction_clusters=score_correction_clusters,
            **kwargs,
        )
        return pipeline.run()

    # ── Multiple genes ───────────────────────────────────────────────────────
    from collections import OrderedDict
    all_results: Dict[str, Any] = {}
    merged_vcf = os.path.join(output_dir, "multi_genes_simulated.vcf")
    vcf_lines: List[str] = []

    for gene, tx in zip(gene_list, tx_list):
        gene_out = output_dir if merge_output else os.path.join(output_dir, f"gene_{gene}")
        Path(gene_out).mkdir(parents=True, exist_ok=True)

        pipeline = MatchingPipeline(
            gtf_file=gtf_file,
            fasta_file=fasta_file,
            gene_name=gene,
            transcript_id=tx,
            database_dir=database_dir,
            output_dir=gene_out,
            simulation_strategy=simulation_strategy,
            min_probability=min_probability,
            constraint_filter=constraint_filter,
            spectrum_ratio=spectrum_ratio,
            score_correction_method=score_correction_method,
            score_correction_clusters=score_correction_clusters,
            **kwargs,
        )

        logger.info("=" * 70)
        logger.info(f" Running for gene: {gene} | {tx}")
        logger.info("=" * 70)
        res = pipeline.run()
        all_results[gene] = res

        if merge_output and os.path.exists(res["simulated_vcf"]):
            with open(res["simulated_vcf"], 'r') as f:
                for line in f:
                    line = line.rstrip("\n")
                    if line.startswith("##"):
                        if not vcf_lines:
                            vcf_lines.append(line)
                    elif line.startswith("#"):
                        if not any(l.startswith("#CHROM") for l in vcf_lines):
                            vcf_lines.append(line)
                    else:
                        vcf_lines.append(line)

    # ── Merge output ─────────────────────────────────────────────────────────
    merged_main_result: Optional[Dict[str, Any]] = None
    if merge_output:
        with open(merged_vcf, 'w') as f:
            f.write("\n".join(vcf_lines) + "\n")

        logger.info("========== Merged annotation ==========")
        merged_annotated_tsv = _run_annotation_for_vcf(
            simulated_vcf=merged_vcf,
            output_dir=output_dir,
            gene_name="merged_genes",
            database_dir=database_dir,
            buildver=kwargs.get('buildver', 'hg19'),
            protocols=kwargs.get('protocols', ['refGene']),
            operations=kwargs.get('operations', ['g']),
            threads=kwargs.get('threads', 4),
        )
        scores  = _calculate_auroc_scores_from_tsv(merged_annotated_tsv, merged_vcf)
        figs    = _generate_visualizations_for_scores(scores, output_dir, "merged_genes")
        total_variants = sum(r.get('total_variants', 0) for r in all_results.values())

        merged_main_result = {
            "gene_name":        "merged_genes",
            "transcript_id":    "merged_transcripts",
            "simulated_vcf":    merged_vcf,
            "annotated_tsv":    merged_annotated_tsv,
            "total_variants":   total_variants,
            "auroc_scores":     scores,
            "figures":          figs,
            "all_gene_results": all_results,
        }
        _save_summary_for_result(merged_main_result, output_dir, "merged_genes")
        return merged_main_result

    return {"mode": "separate_output", "all_gene_results": all_results}


def run_pipeline_from_args(args: argparse.Namespace) -> Dict[str, Any]:
    """Run pipeline from parsed command-line arguments (moved to pipeline.py for centralisation)."""
    return run_pipeline(
        gtf_file=args.gtf,
        fasta_file=args.fasta,
        gene_name=args.gene,
        transcript_id=args.transcript,
        database_dir=args.database,
        output_dir=args.output_dir,
        merge_output=args.merge_output,
        variant_types=[v.strip() for v in args.variant_types.split(',')],
        protocols=[p.strip() for p in args.protocols.split(',')],
        operations=[o.strip() for o in args.operations.split(',')],
        buildver=args.buildver,
        threads=args.threads,
        max_indel_length=args.max_indel_length,
        min_indel_length=args.min_indel_length,
        synonymous=args.synonymous,
        include_stop_codon=args.include_stop_codon,
        max_splice_offset=args.max_splice_offset,
        min_splice_offset=args.min_splice_offset,
        include_classic_splice_sites=args.include_classic_splice_sites,
        max_variants=args.max_variants,
        simulation_strategy=getattr(args, 'simulation_strategy', 'traditional'),
        min_probability=getattr(args, 'min_probability', 0.001),
        constraint_filter=getattr(args, 'constraint_filter', False),
        spectrum_ratio=getattr(args, 'spectrum_ratio', 0.7),
        score_correction_method=getattr(args, 'score_correction_method', 'none'),
        score_correction_clusters=getattr(args, 'score_correction_clusters', 3),
    )


def run_clinvar_roc_evaluation(
    annotated_tsv: str,
    clinvar_csv: str,
    gene_name: str,
    output_dir: str,
    vcf_path: Optional[str] = None,
) -> Dict[str, Any]:
    """
    Run ClinVar-based ROC evaluation comparing Total_Score vs AlphaMissense.

    This is the main evaluation function for the matchvar-pipeline CLI command.
    It loads ClinVar data, extracts benign/pathogenic labels, matches against
    annotated variants, and produces ROC curves with comparison metrics.

    Args:
        annotated_tsv: Path to annotated multianno TSV file
        clinvar_csv: Path to ClinVar CSV with GeneSymbol, ClinicalSignificance, etc.
        gene_name: Target gene name for filtering ClinVar data
        output_dir: Directory for output figures and results
        vcf_path: Optional path to companion VCF for extracting Total_Score

    Returns:
        Dictionary with ROC analysis results and figure paths
    """
    from .evaluation import ClinVarProcessor, ROCAnalyzer, PublicationPlotGenerator

    Path(output_dir).mkdir(parents=True, exist_ok=True)
    fig_dir = Path(output_dir) / "figures"
    fig_dir.mkdir(parents=True, exist_ok=True)

    # Load ClinVar data and extract binary labels
    clinvar = ClinVarProcessor(clinvar_csv)
    clinvar.load_data()
    gene_data = clinvar.filter_by_gene(gene_name)

    # Read annotated TSV
    df = pd.read_csv(annotated_tsv, sep='\t', low_memory=False)

    # Merge Total_Score from VCF if needed and available
    if vcf_path and 'Total_Score' not in df.columns:
        from .visualization import merge_vcf_scores_into_tsv
        df = merge_vcf_scores_into_tsv(annotated_tsv, vcf_path)

    # Find matches between annotated variants and ClinVar
    matches = clinvar.find_matches(df, gene_data)

    if not matches:
        logger.warning(f"No ClinVar matches found for gene {gene_name}")
        return {"error": "No matches found", "gene_name": gene_name}

    logger.info(f"Found {len(matches)} ClinVar matches for {gene_name}")

    # Build comparison result for ROCAnalyzer
    comparison_result = {
        "matches": matches,
        "total_matches": len(matches),
        "total_simulated": len(df),
    }

    # Run ROC analysis
    roc_analyzer = ROCAnalyzer(output_dir=output_dir)
    roc_result = roc_analyzer.perform_roc_analysis(
        comparison_result, gene_name, save_results=True
    )

    # Generate plots
    plot_generator = PublicationPlotGenerator(output_dir=output_dir)
    analysis_result = {"basic_statistics": {"total_variants": len(df)}}
    plots = plot_generator.create_comprehensive_analysis_plot(
        analysis_result, comparison_result, roc_result, gene_name
    )

    return {
        "gene_name": gene_name,
        "total_matches": len(matches),
        "roc_result": roc_result,
        "figures": plots,
        "comparison_result": comparison_result,
    }
