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
import matplotlib.pyplot as plt

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
        logger.info(f"Annotation completed: {self.annotated_tsv}")

    def _load_and_preprocess_tsv(self, tsv_path: str) -> pd.DataFrame:
        """
        Load TSV and preprocess for calibration/evaluation.

        - Filters ClinicalSignificance to 6 allowed values
        - Extracts TOTAL_SCORE from Otherinfo11 where TYPE=SNV
        - Creates binary Label column
        """
        import re

        df = pd.read_csv(tsv_path, sep='\t', low_memory=False)

        valid_clinical = ['Benign/Likely benign', 'Likely benign', 'Benign',
                          'Pathogenic/Likely pathogenic', 'Pathogenic', 'Likely pathogenic']

        if 'ClinicalSignificance' in df.columns:
            df = df[df['ClinicalSignificance'].isin(valid_clinical)].copy()

        type_pattern = re.compile(r'TYPE=SNV')
        total_score_pattern = re.compile(r'TOTAL_SCORE=([0-9.+-eE]+)')

        def extract_total_score(row):
            for col in row.index:
                if col.startswith('Otherinfo') and pd.notna(row[col]):
                    text = str(row[col])
                    if type_pattern.search(text):
                        match = total_score_pattern.search(text)
                        if match:
                            return float(match.group(1))
            return np.nan

        df['TOTAL_SCORE'] = df.apply(extract_total_score, axis=1)

        def get_label(clinical_sig):
            if pd.isna(clinical_sig):
                return np.nan
            sig_lower = str(clinical_sig).lower()
            if 'benign' in sig_lower:
                return 0
            elif 'pathogenic' in sig_lower:
                return 1
            return np.nan

        df['Label'] = df['ClinicalSignificance'].apply(get_label)

        return df

    def _run_calibration(self, clinvar_validation_set_path: str):
        """
        Calibrate TOTAL_SCORE using gene-specific or global ClinVar data.

        Trains a LogisticRegression model on ClinVar variants and applies
        it to the simulated variants to obtain calibrated harmfulness scores.
        """
        df_sim = self._load_and_preprocess_tsv(self.annotated_tsv)

        if 'TOTAL_SCORE' not in df_sim.columns or df_sim['TOTAL_SCORE'].isna().all():
            raise ValueError("No TOTAL_SCORE values found in simulated variants")
        
        df_clinvar = self._load_and_preprocess_tsv(clinvar_validation_set_path)
        
        df_gene = df_clinvar[
            df_clinvar['Gene'].str.upper() == self.gene_name.upper()
        ].copy() if 'Gene' in df_clinvar.columns else pd.DataFrame()
        
        use_global = False
        if len(df_gene) < 10 or df_gene['Label'].nunique() < 2:
            logger.warning(
                f"Insufficient gene-specific ClinVar data for {self.gene_name} "
                f"(n={len(df_gene)}, unique labels={df_gene['Label'].nunique()}). "
                "Falling back to global ClinVar model."
            )
            use_global = True
            train_df = df_clinvar[df_clinvar['Label'].notna()].copy()
        else:
            train_df = df_gene[df_gene['Label'].notna()].copy()
        
        if len(train_df) < 10 or train_df['Label'].nunique() < 2:
            raise ValueError(
                "Insufficient training data for calibration (need ≥10 samples with both classes)"
            )
        
        X_train = train_df['TOTAL_SCORE'].values.reshape(-1, 1)
        y_train = train_df['Label'].values
        
        model = LogisticRegression(random_state=42, max_iter=1000)
        model.fit(X_train, y_train)
        logger.info(
            f"Trained {'gene-specific' if not use_global else 'global'} "
            f"LogisticRegression on {len(train_df)} ClinVar variants"
        )
        
        valid_mask = df_sim['TOTAL_SCORE'].notna()
        df_sim.loc[valid_mask, 'calibrated_score'] = model.predict_proba(
            df_sim.loc[valid_mask, 'TOTAL_SCORE'].values.reshape(-1, 1)
        )[:, 1]
        
        self.calibrated_scores = df_sim['calibrated_score'].values
        self.df_calibrated = df_sim
        
        calibrated_output = os.path.join(
            self.output_dir, f"{self.gene_name}_calibrated_scores.tsv"
        )
        df_sim.to_csv(calibrated_output, sep='\t', index=False)
        logger.info(f"Calibrated scores saved to: {calibrated_output}")

    def _run_evaluation(self):
        """
        Compute auROC for am_pathogenicity and calibrated scores,
        and generate publication-quality ROC comparison figure.
        """
        from sklearn.metrics import roc_auc_score, roc_curve
        
        if not hasattr(self, 'df_calibrated') or self.df_calibrated is None:
            raise RuntimeError(
                "Run _run_calibration first before evaluation"
            )
        
        df = self.df_calibrated
        valid_mask = df['Label'].notna()
        
        if valid_mask.sum() == 0:
            raise ValueError("No valid labels for evaluation")
        
        y_true = df.loc[valid_mask, 'Label'].values
        
        scores_to_eval = {}
        
        if 'am_pathogenicity' in df.columns:
            am_scores = pd.to_numeric(df.loc[valid_mask, 'am_pathogenicity'], errors='coerce').values
            if not np.isnan(am_scores).all():
                scores_to_eval['am_pathogenicity'] = am_scores
        
        if 'calibrated_score' in df.columns:
            cal_scores = pd.to_numeric(df.loc[valid_mask, 'calibrated_score'], errors='coerce').values
            if not np.isnan(cal_scores).all():
                scores_to_eval['calibrated_score'] = cal_scores
        
        auroc_results = {}
        roc_curves = {}
        
        for name, scores in scores_to_eval.items():
            valid_idx = ~np.isnan(scores)
            if valid_idx.sum() >= 2:
                try:
                    auroc = roc_auc_score(y_true[valid_idx], scores[valid_idx])
                    fpr, tpr, _ = roc_curve(y_true[valid_idx], scores[valid_idx])
                    auroc_results[name] = auroc
                    roc_curves[name] = {'fpr': fpr, 'tpr': tpr}
                    logger.info(f"{name} auROC: {auroc:.4f}")
                except Exception as e:
                    logger.warning(f"Failed to compute auROC for {name}: {e}")
        
        self.auroc_results = auroc_results
        
        fig, ax = plt.subplots(figsize=(8, 8), dpi=300)
        
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        
        ax.plot([0, 1], [0, 1], color='gray', lw=1, ls=':', alpha=0.7)
        
        colors = {'calibrated_score': '#D62728', 'am_pathogenicity': '#1F77B4'}
        linestyles = {'calibrated_score': '-', 'am_pathogenicity': '--'}
        linewidths = {'calibrated_score': 2.5, 'am_pathogenicity': 2.0}
        
        for name, color in [('calibrated_score', '#D62728'), ('am_pathogenicity', '#1F77B4')]:
            if name in roc_curves:
                fpr = roc_curves[name]['fpr']
                tpr = roc_curves[name]['tpr']
                auroc = auroc_results[name]
                label_name = 'Calibrated Score' if name == 'calibrated_score' else 'am_pathogenicity'
                ax.plot(fpr, tpr, color=color, lw=linewidths.get(name, 2.0),
                        ls=linestyles.get(name, '-'),
                        label=f"{label_name}  (AUC={auroc:.3f})")
        
        ax.set_xlim(0.0, 1.0)
        ax.set_ylim(0.0, 1.05)
        ax.set_xlabel('False Positive Rate', fontweight='bold')
        ax.set_ylabel('True Positive Rate', fontweight='bold')
        ax.set_title(f'ROC Curve Comparison for {self.gene_name}', fontweight='bold', fontsize=14)
        ax.legend(loc='lower right', fontsize=10)
        ax.grid(True, alpha=0.3)
        
        pdf_path = os.path.join(self.output_dir, f"{self.gene_name}_ROC_comparison.pdf")
        fig.savefig(pdf_path, dpi=300, bbox_inches='tight', format='pdf')
        plt.close(fig)
        
        logger.info(f"ROC comparison figure saved: {pdf_path}")

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

# create labels from vcf 
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
    actual_output = os.path.join(output_dir, f"{gene_name}_annotated.{buildver}_multianno.tsv")

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

    logger.info(f"Annotation completed: {actual_output}")
    return actual_output


# old run pipeline
# def run_pipeline(gtf_file: str, fasta_file: str, gene_name: str, transcript_id: str,
#                  database_dir: str, output_dir: str, **kwargs) -> Dict[str, Any]:
#     """
#     Run pipeline directly from Python code

#     Args:
#         gtf_file: Path to GTF file
#         fasta_file: Path to reference genome FASTA
#         gene_name: Gene symbol
#         transcript_id: Transcript ID
#         database_dir: Database directory
#         output_dir: Output directory
#         **kwargs: Additional arguments for MatchingPipeline

#     Returns:
#         Dictionary with pipeline results
#     """
#     pipeline = MatchingPipeline(
#         gtf_file=gtf_file,
#         fasta_file=fasta_file,
#         gene_name=gene_name,
#         transcript_id=transcript_id,
#         database_dir=database_dir,
#         output_dir=output_dir,
#         **kwargs
#     )
#     return pipeline.run()

# new run pipeline
def run_pipeline(gtf_file: str, fasta_file: str, gene_name: str, transcript_id: str,
                 database_dir: str, output_dir: str, merge_output: bool = False, **kwargs) -> Dict[str, Any]:
    """
    Run pipeline directly from Python code

    Args:
        gtf_file: Path to GTF file
        fasta_file: Path to reference genome FASTA
        gene_name: Gene symbol(s), comma-separated for multiple
        transcript_id: Transcript ID(s), comma-separated matching genes
        database_dir: Database directory
        output_dir: Output directory
        merge_output: If True, merge all genes into one VCF/annotation; False = separate files
        **kwargs: Additional arguments for MatchingPipeline

    Returns:
        Dictionary with pipeline results
    """
    # Analyzing multiple genes/transcripts
    gene_list = [g.strip() for g in gene_name.split(",") if g.strip()]
    tx_list = [t.strip() for t in transcript_id.split(",") if t.strip()]

    if len(gene_list) != len(tx_list):
        raise ValueError(f"Gene counts({len(gene_list)}) != Transcript counts({len(tx_list)})")

    # ===================== Single gene =====================
    if len(gene_list) == 1:
        pipeline = MatchingPipeline(
            gtf_file=gtf_file,
            fasta_file=fasta_file,
            gene_name=gene_list[0],
            transcript_id=tx_list[0],
            database_dir=database_dir,
            output_dir=output_dir,
            **kwargs
        )
        return pipeline.run()

    # ===================== Multi genes =====================
    from collections import OrderedDict
    all_results = {}
    merged_vcf = os.path.join(output_dir, "multi_genes_simulated.vcf")
    vcf_lines = []

    # Run one by one
    for gene, tx in zip(gene_list, tx_list):
        if not merge_output:
            gene_out = os.path.join(output_dir, f"gene_{gene}")
            Path(gene_out).mkdir(parents=True, exist_ok=True)
        else:
            gene_out = output_dir

        pipeline = MatchingPipeline(
            gtf_file=gtf_file,
            fasta_file=fasta_file,
            gene_name=gene,
            transcript_id=tx,
            database_dir=database_dir,
            output_dir=gene_out,
            **kwargs
        )

        logger.info(f"\n=== Running for gene: {gene} | {tx} ===")
        res = pipeline.run()
        all_results[gene] = res

        # merge VCF
        if merge_output and os.path.exists(res["simulated_vcf"]):
            with open(res["simulated_vcf"], 'r') as f:
                for line in f:
                    line = line.rstrip("\n")
                    if line.startswith("##"):
                        # Keep ## meta-headers from the first file only
                        if not vcf_lines:
                            vcf_lines.append(line)
                    elif line.startswith("#"):
                        # Keep #CHROM header line from the first file only
                        if not any(l.startswith("#CHROM") for l in vcf_lines):
                            vcf_lines.append(line)
                    else:
                        # Always keep variant data lines
                        vcf_lines.append(line)

    # Write merged VCF and all annotation
    merged_main_result = None
    if merge_output:
        with open(merged_vcf, 'w') as f:
            f.write("\n".join(vcf_lines) + "\n")

        # Annotate the merged VCF using the standalone helper (no dummy pipeline needed)
        logger.info("\n=== Merged annotation ===")
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
        scores = _calculate_auroc_scores_from_tsv(merged_annotated_tsv, merged_vcf)
        figs = _generate_visualizations_for_scores(scores, output_dir, "merged_genes")
        merged_main_result = {
            "gene_name": "merged_genes",
            "transcript_id": "merged_transcripts",
            "simulated_vcf": merged_vcf,
            "annotated_tsv": merged_annotated_tsv,
            "auroc_scores": scores,
            "figures": figs,
            "all_gene_results": all_results
        }
        _save_summary_for_result(merged_main_result, output_dir, "merged_genes")
        return merged_main_result

    return {"mode": "separate_output", "all_gene_results": all_results}


def _calculate_auroc_scores_from_tsv(annotated_tsv: str, vcf_path: str) -> Dict:
    """Standalone auROC calculation from an annotated TSV path."""
    import tempfile
    tmp = type('TmpPipeline', (), {})()
    tmp.annotated_tsv = annotated_tsv
    tmp.simulated_vcf = vcf_path
    # Reuse the instance method via a temporary object
    return MatchingPipeline._calculate_auroc_scores(tmp)


def _generate_visualizations_for_scores(
    scores: Dict, output_dir: str, gene_name: str
) -> Dict[str, str]:
    """Standalone visualization generation from scores dict."""
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


def _save_summary_for_result(results: Dict[str, Any], output_dir: str, gene_name: str):
    """Standalone summary saver for a results dict."""
    summary_file = os.path.join(output_dir, f"{gene_name}_pipeline_summary.json")
    summary = {
        'gene_name': results.get('gene_name', gene_name),
        'transcript_id': results.get('transcript_id', ''),
        'timestamp': datetime.now().isoformat(),
        'total_variants': results.get('total_variants', 0),
        'simulated_vcf': results.get('simulated_vcf', ''),
        'annotated_tsv': results.get('annotated_tsv', ''),
        'auroc_scores': results.get('auroc_scores', {}),
        'figures': {k: v for k, v in results.get('figures', {}).items() if k.endswith('_pdf')}
    }
    with open(summary_file, 'w', encoding='utf-8') as f:
        json.dump(summary, f, indent=2)
    logger.info(f"Summary saved: {summary_file}")

    stats_file = os.path.join(output_dir, f"{gene_name}_auroc_statistics.tsv")
    if results.get('auroc_scores'):
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


def run_pipeline_from_args(args):
    """Run pipeline from parsed command-line arguments"""
    return run_pipeline(
        gtf_file=args.gtf,
        fasta_file=args.fasta,
        gene_name=args.gene,
        transcript_id=args.transcript,
        database_dir=args.database,
        output_dir=args.output_dir,
        merge_output=args.merge_output,
        protocols=args.protocols.split(',') if args.protocols else None,
        operations=args.operations.split(',') if args.operations else None,
        variant_types=args.variant_types.split(',') if args.variant_types else None,
        buildver=args.buildver,
        threads=args.threads
    )