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

from .vsimulator import simulate_variants, GeneTranscript
from .table_matchvar import TableAnnotator

logger = logging.getLogger(__name__)


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
                 threads: int = 4):
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

        self.transcript = simulate_variants(
            gtf_file=self.gtf_file,
            fasta_file=self.fasta_file,
            gene_name=self.gene_name,
            transcript_id=self.transcript_id,
            output_vcf=self.simulated_vcf,
            variant_types=self.variant_types
        )

        # Get variants by loading from VCF or regenerate
        variants_result = self.transcript.generate_all_variants(variant_types=self.variant_types)
        self.variants = variants_result
        # Extract variant lists (excluding the 'total' key)
        variant_lists = {k: v for k, v in variants_result.items() if k != 'total'}
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

    def _calculate_auroc_scores(self) -> Dict[str, Dict[str, float]]:
        """
        Calculate auROC scores for each annotation tool

        Returns:
            Dictionary: {tool_name: {'auroc': float, 'auprc': float, ...}}
        """
        if not self.annotated_tsv or not os.path.exists(self.annotated_tsv):
            raise FileNotFoundError(f"Annotated file not found: {self.annotated_tsv}")

        df = pd.read_csv(self.annotated_tsv, sep='\t', low_memory=False)

        if 'Total_Score' not in df.columns:
            logger.warning("Total_Score column not found, skipping auROC calculation")
            return {}

        y_true = self._get_true_labels(df)
        if y_true is None or len(y_true) == 0:
            logger.warning("No ground truth labels found")
            return {}

        scores = {}
        score_columns = [col for col in df.columns if 'Score' in col or 'score' in col.lower()]

        for col in score_columns:
            if col == 'Total_Score':
                continue

            y_scores = pd.to_numeric(df[col], errors='coerce').values
            valid_mask = ~np.isnan(y_scores)

            if valid_mask.sum() < 10:
                logger.debug(f"Skipping {col}: insufficient valid scores")
                continue

            from sklearn.metrics import roc_auc_score, average_precision_score, roc_curve

            try:
                auroc = roc_auc_score(y_true[valid_mask], y_scores[valid_mask])
                auprc = average_precision_score(y_true[valid_mask], y_scores[valid_mask])

                fpr, tpr, _ = roc_curve(y_true[valid_mask], y_scores[valid_mask])

                scores[col] = {
                    'auroc': round(float(auroc), 4),
                    'auprc': round(float(auprc), 4),
                    'roc_curve': {
                        'fpr': fpr.tolist(),
                        'tpr': tpr.tolist()
                    },
                    'n_variants': int(valid_mask.sum())
                }

                logger.info(f"  {col}: AUROC={auroc:.4f}, AUPRC={auprc:.4f}")

            except Exception as e:
                logger.warning(f"Failed to calculate metrics for {col}: {e}")

        return scores

    def _get_true_labels(self, df: pd.DataFrame) -> Optional[np.ndarray]:
        """
        Extract ground truth binary labels from annotated dataframe
        Labels derived from variant type: pathogenic variants (loss-of-function, splicing) = 1
        """
        try:
            if 'Func.refGene' in df.columns:
                func = df['Func.refGene'].fillna('').str.lower()
                y_true = func.str.contains('splicing|stopgain|stoploss|frameshift|nonsyn', na=False).astype(int).values
                logger.info(f"Extracted {y_true.sum()} positive labels from functional annotation")
                return y_true

            elif 'ExonicFunc.refGene' in df.columns:
                exonic_func = df['ExonicFunc.refGene'].fillna('').str.lower()
                y_true = exonic_func.str.contains('splicing|stopgain|stoploss|frameshift|nonsynonymous', na=False).astype(int).values
                logger.info(f"Extracted {y_true.sum()} positive labels from exonic function")
                return y_true

            elif 'TYPE' in df.columns:
                var_type = df['TYPE'].fillna('').str.upper()
                y_true = var_type.isin(['SPLICE_SITE', 'FRAMESHIFT', 'NONSENSE', 'STOPLOSS']).astype(int).values
                logger.info(f"Extracted {y_true.sum()} positive labels from variant type")
                return y_true

            else:
                logger.error("No suitable column for ground truth extraction")
                return None

        except Exception as e:
            logger.error(f"Error extracting true labels: {e}")
            return None

    def _generate_visualizations(self, scores: Dict[str, Dict]) -> Dict[str, str]:
        """
        Generate ROC curves and auROC comparison bar plots

        Args:
            scores: Dictionary of metric scores per tool

        Returns:
            Dictionary mapping figure names to file paths
        """
        if not scores:
            logger.warning("No scores available for visualization")
            return {}

        figures = {}

        from matplotlib import pyplot as plt
        import seaborn as sns
        sns.set_style("whitegrid")
        plt.rcParams['figure.dpi'] = 300

        fig_dir = os.path.join(self.output_dir, 'figures')
        Path(fig_dir).mkdir(parents=True, exist_ok=True)

        # Figure 1: ROC curves for all tools
        fig_roc, ax_roc = plt.subplots(figsize=(8, 8))

        for tool_name, metrics in scores.items():
            if 'roc_curve' in metrics:
                fpr = np.array(metrics['roc_curve']['fpr'])
                tpr = np.array(metrics['roc_curve']['tpr'])
                ax_roc.plot(fpr, tpr, lw=2,
                           label=f"{tool_name} (AUROC={metrics['auroc']:.3f})")

        ax_roc.plot([0, 1], [0, 1], 'k--', lw=1, alpha=0.7)
        ax_roc.set_xlabel('False Positive Rate', fontsize=12)
        ax_roc.set_ylabel('True Positive Rate', fontsize=12)
        ax_roc.set_title(f'ROC Curves - {self.gene_name}', fontsize=14, fontweight='bold')
        ax_roc.legend(loc='lower right', fontsize=10)
        ax_roc.grid(True, alpha=0.3)

        roc_path = os.path.join(fig_dir, f"{self.gene_name}_roc_curves.png")
        fig_roc.savefig(roc_path, dpi=300, bbox_inches='tight')
        figures['roc_curves'] = roc_path
        logger.info(f"ROC curves saved: {roc_path}")

        # Figure 2: auROC comparison bar plot
        fig_bar, ax_bar = plt.subplots(figsize=(10, 6))

        tools = list(scores.keys())
        aurocs = [scores[t]['auroc'] for t in tools]

        bars = ax_bar.barh(range(len(tools)), aurocs, color='steelblue', alpha=0.8)

        for i, (bar, val) in enumerate(zip(bars, aurocs)):
            ax_bar.text(val + 0.01, bar.get_y() + bar.get_height()/2,
                       f"{val:.3f}", va='center', fontsize=10)

        ax_bar.set_yticks(range(len(tools)))
        ax_bar.set_yticklabels(tools, fontsize=11)
        ax_bar.set_xlabel('auROC Score', fontsize=12)
        ax_bar.set_title(f'Annotation Tool Performance - {self.gene_name}',
                        fontsize=14, fontweight='bold')
        ax_bar.set_xlim(0, 1.05)
        ax_bar.grid(True, axis='x', alpha=0.3)

        bar_path = os.path.join(fig_dir, f"{self.gene_name}_auroc_comparison.png")
        fig_bar.savefig(bar_path, dpi=300, bbox_inches='tight')
        figures['auroc_comparison'] = bar_path
        logger.info(f"auROC comparison saved: {bar_path}")

        # Also save PDF versions for publication
        pdf_roc = roc_path.replace('.png', '.pdf')
        pdf_bar = bar_path.replace('.png', '.pdf')
        fig_roc.savefig(pdf_roc, format='pdf', bbox_inches='tight')
        fig_bar.savefig(pdf_bar, format='pdf', bbox_inches='tight')
        figures['roc_curves_pdf'] = pdf_roc
        figures['auroc_comparison_pdf'] = pdf_bar

        plt.close('all')

        return figures

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