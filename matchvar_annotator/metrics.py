#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Metrics Module for auROC Calculation and Comparison

Provides comprehensive evaluation metrics for variant annotation tools:
- Area Under ROC Curve (auROC)
- Area Under Precision-Recall Curve (auPRC)
- ROC curve data
- Statistical comparisons between tools
"""

import numpy as np
import pandas as pd
from typing import Dict, List, Optional, Tuple, Any
from sklearn.metrics import (
    roc_auc_score,
    average_precision_score,
    roc_curve,
    precision_recall_curve,
    confusion_matrix
)
import logging

logger = logging.getLogger(__name__)


class VariantMetricCalculator:
    """
    Calculate and compare performance metrics for variant annotation tools

    This calculator evaluates how well each annotation tool's scores correlate
    with ground truth pathogenicity labels derived from simulated variants.
    """

    def __init__(self, ground_truth_source: str = 'functional'):
        """
        Initialize metric calculator

        Args:
            ground_truth_source: Source of ground truth labels
                - 'functional': Use gene function annotation (splicing, stopgain, etc.)
                - 'variant_type': Use variant type classification
                - 'custom': Provide custom labels via get_labels method
        """
        self.ground_truth_source = ground_truth_source
        self.results: Dict[str, Dict] = {}

    def evaluate_tool(self,
                     tool_scores: np.ndarray,
                     tool_name: str,
                     y_true: np.ndarray,
                     variant_ids: Optional[List[str]] = None) -> Dict[str, Any]:
        """
        Calculate comprehensive metrics for a single annotation tool

        Args:
            tool_scores: Continuous scores from annotation tool (higher = more pathogenic)
            tool_name: Name of the annotation tool
            y_true: Binary ground truth labels (1 = pathogenic, 0 = benign)
            variant_ids: Optional variant identifiers

        Returns:
            Dictionary containing all calculated metrics
        """
        if len(tool_scores) != len(y_true):
            raise ValueError(f"Score length mismatch: {len(tool_scores)} vs {len(y_true)}")

        valid_mask = ~np.isnan(tool_scores)
        if valid_mask.sum() < 10:
            logger.warning(f"Insufficient valid scores for {tool_name}")
            return {}

        y_true_valid = y_true[valid_mask]
        scores_valid = tool_scores[valid_mask]

        if len(np.unique(y_true_valid)) < 2:
            logger.warning(f"{tool_name}: only one class present, skipping ROC AUC")
            return {}

        try:
            auroc = roc_auc_score(y_true_valid, scores_valid)
            auprc = average_precision_score(y_true_valid, scores_valid)

            fpr, tpr, thresholds = roc_curve(y_true_valid, scores_valid)

            precision, recall, pr_thresholds = precision_recall_curve(y_true_valid, scores_valid)

            # Youden's J statistic for optimal threshold
            j_scores = tpr - fpr
            optimal_idx = np.argmax(j_scores)
            optimal_threshold = thresholds[optimal_idx] if optimal_idx < len(thresholds) else 0.5

            metrics = {
                'tool_name': tool_name,
                'auroc': float(auroc),
                'auprc': float(auprc),
                'n_variants': int(valid_mask.sum()),
                'n_positive': int(y_true_valid.sum()),
                'optimal_threshold': float(optimal_threshold),
                'roc_curve': {
                    'fpr': fpr.tolist(),
                    'tpr': tpr.tolist(),
                    'thresholds': thresholds.tolist()
                },
                'pr_curve': {
                    'precision': precision.tolist(),
                    'recall': recall.tolist(),
                    'thresholds': pr_thresholds.tolist()
                }
            }

            self.results[tool_name] = metrics
            logger.info(f"{tool_name}: AUROC={auroc:.4f}, AUPRC={auprc:.4f}")

            return metrics

        except Exception as e:
            logger.error(f"Failed to calculate metrics for {tool_name}: {e}")
            return {}

    def extract_labels_from_annotation(self,
                                      annotation_df: pd.DataFrame,
                                      label_column: Optional[str] = None) -> np.ndarray:
        """
        Extract ground truth labels from annotated dataframe

        Args:
            annotation_df: Annotated variants dataframe
            label_column: Specific column to use (auto-detected if None)

        Returns:
            Binary array: 1 for pathogenic, 0 for benign
        """
        if label_column and label_column in annotation_df.columns:
            labels = annotation_df[label_column].fillna(0).astype(int).values
            logger.info(f"Using labels from column: {label_column}")
            return labels

        candidate_columns = [
            'Func.refGene', 'ExonicFunc.refGene', 'TYPE',
            'Func.ensGene', 'ExonicFunc.ensGene',
            'Func.knownGene', 'ExonicFunc.knownGene'
        ]

        for col in candidate_columns:
            if col in annotation_df.columns:
                logger.info(f"Auto-detected label column: {col}")
                values = annotation_df[col].fillna('').astype(str).str.lower()

                pathogenic_patterns = [
                    'splicing', 'stopgain', 'stoploss', 'frameshift',
                    'nonsynonymous', 'missense', 'nonframeshift'
                ]

                labels = np.array([
                    1 if any(pat in val for pat in pathogenic_patterns) else 0
                    for val in values
                ])

                logger.info(f"Extracted {labels.sum()} pathogenic / {len(labels)-labels.sum()} benign labels")
                return labels

        logger.error("No suitable label column found in annotation dataframe")
        return np.array([])

    def compare_tools(self) -> pd.DataFrame:
        """
        Compare all evaluated tools in a summary table

        Returns:
            DataFrame with auROC, auPRC, and other statistics
        """
        if not self.results:
            return pd.DataFrame()

        rows = []
        for tool, metrics in self.results.items():
            rows.append({
                'Tool': tool,
                'AUROC': metrics['auroc'],
                'AUPRC': metrics['auprc'],
                'N_Variants': metrics['n_variants'],
                'N_Pathogenic': metrics['n_positive'],
                'N_Benign': metrics['n_variants'] - metrics['n_positive']
            })

        df = pd.DataFrame(rows)
        df = df.sort_values('AUROC', ascending=False)

        return df

    def calculate_statistical_tests(self) -> pd.DataFrame:
        """
        Perform DeLong test for AUROC comparison

        Returns:
            DataFrame with pairwise p-values
        """
        # This would require implementing DeLong test or using external library
        # For now, return empty dataframe
        logger.info("Statistical comparison requires additional implementation")
        return pd.DataFrame()

    def save_results(self, output_prefix: str):
        """
        Save all metrics to files

        Args:
            output_prefix: Prefix for output files
        """
        # Save summary table
        summary_df = self.compare_tools()
        if not summary_df.empty:
            summary_file = f"{output_prefix}_auroc_summary.tsv"
            summary_df.to_csv(summary_file, sep='\t', index=False)
            logger.info(f"Saved metrics summary: {summary_file}")

        # Save full results as JSON
        import json
        results_file = f"{output_prefix}_detailed_metrics.json"
        with open(results_file, 'w') as f:
            json.dump(self.results, f, indent=2)
        logger.info(f"Saved detailed metrics: {results_file}")


def calculate_auroc_from_dataframe(df: pd.DataFrame,
                                   score_columns: List[str],
                                   label_column: str = 'label') -> Dict[str, float]:
    """
    Convenience function to calculate auROC for multiple score columns

    Args:
        df: DataFrame containing scores and labels
        score_columns: List of column names with prediction scores
        label_column: Name of column with binary true labels

    Returns:
        Dictionary mapping tool names to auROC scores
    """
    if label_column not in df.columns:
        raise ValueError(f"Label column '{label_column}' not found")

    y_true = df[label_column].values
    calculator = VariantMetricCalculator()
    scores = {}

    for col in score_columns:
        if col in df.columns:
            tool_scores = pd.to_numeric(df[col], errors='coerce').values
            metrics = calculator.evaluate_tool(tool_scores, col, y_true)
            if metrics:
                scores[col] = metrics['auroc']

    return scores


def bootstrap_confidence_interval(y_true: np.ndarray,
                                 y_scores: np.ndarray,
                                 n_bootstrap: int = 1000,
                                 alpha: float = 0.05) -> Tuple[float, float, float]:
    """
    Calculate bootstrap confidence interval for AUROC

    Args:
        y_true: True binary labels
        y_scores: Predicted scores
        n_bootstrap: Number of bootstrap samples
        alpha: Significance level (default 0.05 for 95% CI)

    Returns:
        (auroc, lower_bound, upper_bound)
    """
    from sklearn.utils import resample

    n_samples = len(y_true)
    auroc_base = roc_auc_score(y_true, y_scores)

    bootstrap_aurocs = []
    for _ in range(n_bootstrap):
        indices = resample(range(n_samples), n_samples=n_samples)
        y_true_bs = y_true[indices]
        y_scores_bs = y_scores[indices]

        if len(np.unique(y_true_bs)) >= 2:
            auroc_bs = roc_auc_score(y_true_bs, y_scores_bs)
            bootstrap_aurocs.append(auroc_bs)

    if bootstrap_aurocs:
        lower = np.percentile(bootstrap_aurocs, 100 * alpha / 2)
        upper = np.percentile(bootstrap_aurocs, 100 * (1 - alpha / 2))
        return auroc_base, lower, upper
    else:
        return auroc_base, np.nan, np.nan