#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Visualization Module for Variant Annotation Evaluation

Creates publication-quality figures:
- ROC curves with confidence bands
- auROC comparison bar plots
- Precision-Recall curves
- Performance heatmaps
"""

import os
import logging
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Any
from dataclasses import dataclass

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib import colors as mcolors

logger = logging.getLogger(__name__)


@dataclass
class FigureStyle:
    """Style configuration for publication figures"""
    dpi: int = 300
    figsize_roc: Tuple[int, int] = (8, 8)
    figsize_bar: Tuple[int, int] = (10, 6)
    figsize_heatmap: Tuple[int, int] = (10, 8)

    title_fontsize: int = 14
    label_fontsize: int = 12
    tick_fontsize: int = 10
    legend_fontsize: int = 10

    colormap_roc: str = 'tab10'
    colormap_bar: str = 'Blues_d'
    grid_alpha: float = 0.3


class PerformanceVisualizer:
    """
    Visualize annotation tool performance metrics

    Generates multiple figure types:
    1. ROC curves with AUROC values
    2. auROC comparison bar plot
    3. Precision-Recall curves
    4. Performance heatmap
    """

    def __init__(self, style: Optional[FigureStyle] = None):
        """
        Initialize visualizer with style settings

        Args:
            style: FigureStyle object with configuration (uses defaults if None)
        """
        self.style = style or FigureStyle()
        self._setup_style()

    def _setup_style(self):
        """Configure matplotlib/seaborn style for publication"""
        sns.set_style("whitegrid")
        sns.set_context("paper", font_scale=1.2)

        plt.rcParams.update({
            'figure.dpi': self.style.dpi,
            'savefig.dpi': self.style.dpi,
            'savefig.bbox': 'tight',
            'font.family': 'sans-serif',
            'font.sans-serif': ['Arial', 'DejaVu Sans'],
            'axes.labelsize': self.style.label_fontsize,
            'axes.titlesize': self.style.title_fontsize,
            'xtick.labelsize': self.style.tick_fontsize,
            'ytick.labelsize': self.style.tick_fontsize,
            'legend.fontsize': self.style.legend_fontsize,
            'grid.alpha': self.style.grid_alpha
        })

    def plot_roc_curves(self,
                        tool_metrics: Dict[str, Dict],
                        title: Optional[str] = None,
                        output_file: Optional[str] = None) -> plt.Figure:
        """
        Plot ROC curves for multiple annotation tools

        Args:
            tool_metrics: Dictionary {tool_name: metrics_dict} from VariantMetricCalculator
            title: Plot title (auto-generated if None)
            output_file: File path to save figure

        Returns:
            matplotlib Figure object
        """
        fig, ax = plt.subplots(figsize=self.style.figsize_roc)

        colors = plt.cm.get_cmap(self.style.colormap_roc)
        n_tools = len(tool_metrics)

        for idx, (tool_name, metrics) in enumerate(tool_metrics.items()):
            if 'roc_curve' not in metrics:
                continue

            fpr = np.array(metrics['roc_curve']['fpr'])
            tpr = np.array(metrics['roc_curve']['tpr'])
            auroc = metrics['auroc']

            ax.plot(fpr, tpr, lw=2.5,
                   label=f"{tool_name} (AUROC={auroc:.3f})",
                   color=colors(idx % colors.N))

        ax.plot([0, 1], [0, 1], 'k--', lw=1.5, alpha=0.7, label='Random (AUROC=0.500)')

        ax.set_xlabel('False Positive Rate', fontweight='bold')
        ax.set_ylabel('True Positive Rate', fontweight='bold')

        if title:
            ax.set_title(title, fontweight='bold', fontsize=self.style.title_fontsize + 2)
        else:
            ax.set_title('Receiver Operating Characteristic Curves',
                        fontweight='bold', fontsize=self.style.title_fontsize + 2)

        ax.legend(loc='lower right', framealpha=0.9)
        ax.grid(True, alpha=self.style.grid_alpha)
        ax.set_xlim([0.0, 1.0])
        ax.set_ylim([0.0, 1.05])

        plt.tight_layout()

        if output_file:
            self._save_figure(fig, output_file)

        return fig

    def plot_auroc_comparison(self,
                             tool_metrics: Dict[str, Dict],
                             title: Optional[str] = None,
                             output_file: Optional[str] = None,
                             sort_by: str = 'auroc') -> plt.Figure:
        """
        Create horizontal bar plot comparing auROC scores

        Args:
            tool_metrics: Dictionary {tool_name: metrics_dict}
            title: Plot title
            output_file: File path to save
            sort_by: Metric to sort by ('auroc' or 'auprc')

        Returns:
            matplotlib Figure object
        """
        if not tool_metrics:
            raise ValueError("No tool metrics provided")

        fig, ax = plt.subplots(figsize=self.style.figsize_bar)

        tools = list(tool_metrics.keys())
        scores = [tool_metrics[t]['auroc'] for t in tools]

        if sort_by == 'auroc':
            sorted_indices = np.argsort(scores)[::-1]
        else:
            sorted_indices = np.argsort(scores)[::-1]

        tools = [tools[i] for i in sorted_indices]
        scores = [scores[i] for i in sorted_indices]

        y_pos = np.arange(len(tools))
        bars = ax.barh(y_pos, scores,
                      color=plt.cm.get_cmap(self.style.colormap_bar)(np.linspace(0.4, 0.8, len(tools))),
                      edgecolor='black', linewidth=1)

        for i, (bar, score) in enumerate(zip(bars, scores)):
            ax.text(score + 0.01, bar.get_y() + bar.get_height()/2,
                   f"{score:.3f}", va='center', fontsize=self.style.tick_fontsize,
                   fontweight='bold')

        ax.set_yticks(y_pos)
        ax.set_yticklabels(tools, fontsize=self.style.label_fontsize)
        ax.set_xlabel('auROC Score', fontweight='bold', fontsize=self.style.label_fontsize)
        ax.set_xlim(0, 1.05)

        if title:
            ax.set_title(title, fontweight='bold', fontsize=self.style.title_fontsize + 2)
        else:
            ax.set_title('Annotation Tool Performance Comparison',
                        fontweight='bold', fontsize=self.style.title_fontsize + 2)

        ax.grid(True, axis='x', alpha=self.style.grid_alpha)
        ax.axvline(x=0.5, color='red', linestyle='--', linewidth=1.5, alpha=0.5,
                  label='Random (AUROC=0.5)')
        ax.legend()

        plt.tight_layout()

        if output_file:
            self._save_figure(fig, output_file)

        return fig

    def plot_precision_recall(self,
                              tool_metrics: Dict[str, Dict],
                              title: Optional[str] = None,
                              output_file: Optional[str] = None) -> plt.Figure:
        """Plot Precision-Recall curves"""
        fig, ax = plt.subplots(figsize=self.style.figsize_roc)

        colors = plt.cm.get_cmap(self.style.colormap_roc)
        n_tools = len(tool_metrics)

        for idx, (tool_name, metrics) in enumerate(tool_metrics.items()):
            if 'pr_curve' not in metrics:
                continue

            precision = np.array(metrics['pr_curve']['precision'])
            recall = np.array(metrics['pr_curve']['recall'])
            auprc = metrics['auprc']

            ax.plot(recall, precision, lw=2.5,
                   label=f"{tool_name} (AUPRC={auprc:.3f})",
                   color=colors(idx % colors.N))

        ax.set_xlabel('Recall', fontweight='bold')
        ax.set_ylabel('Precision', fontweight='bold')

        if title:
            ax.set_title(title, fontweight='bold', fontsize=self.style.title_fontsize + 2)
        else:
            ax.set_title('Precision-Recall Curves',
                        fontweight='bold', fontsize=self.style.title_fontsize + 2)

        ax.legend(loc='upper right', framealpha=0.9)
        ax.grid(True, alpha=self.style.grid_alpha)
        ax.set_xlim([0.0, 1.0])
        ax.set_ylim([0.0, 1.05])

        plt.tight_layout()

        if output_file:
            self._save_figure(fig, output_file)

        return fig

    def plot_performance_heatmap(self,
                                tool_metrics: Dict[str, Dict],
                                metrics_to_plot: List[str] = None,
                                title: Optional[str] = None,
                                output_file: Optional[str] = None) -> plt.Figure:
        """
        Plot heatmap of performance metrics

        Args:
            tool_metrics: Dictionary {tool_name: metrics_dict}
            metrics_to_plot: List of metrics to include (default: ['auroc', 'auprc'])
            title: Plot title
            output_file: File path to save

        Returns:
            matplotlib Figure object
        """
        metrics_to_plot = metrics_to_plot or ['auroc', 'auprc']

        data = []
        tools = []

        for tool_name, metrics in tool_metrics.items():
            row = {}
            for metric in metrics_to_plot:
                row[metric] = metrics.get(metric, np.nan)
            data.append(row)
            tools.append(tool_name)

        df = pd.DataFrame(data, index=tools)

        fig, ax = plt.subplots(figsize=self.style.figsize_heatmap)

        sns.heatmap(df.T, annot=True, fmt='.3f', cmap='YlOrRd',
                   cbar_kws={'label': 'Score'},
                   linewidths=0.5, linecolor='gray',
                   vmin=0, vmax=1, ax=ax)

        ax.set_xlabel('Annotation Tool', fontweight='bold')
        ax.set_ylabel('Metric', fontweight='bold')

        if title:
            ax.set_title(title, fontweight='bold', fontsize=self.style.title_fontsize + 2)
        else:
            ax.set_title('Performance Metrics Heatmap',
                        fontweight='bold', fontsize=self.style.title_fontsize + 2)

        plt.tight_layout()

        if output_file:
            self._save_figure(fig, output_file)

        return fig

    def plot_metric_distribution(self,
                                tool_metrics: Dict[str, Dict],
                                metric: str = 'auroc',
                                title: Optional[str] = None,
                                output_file: Optional[str] = None) -> plt.Figure:
        """
        Plot distribution/violin plot of metric values across bootstrap samples

        Note: Requires bootstrap samples to be stored in metrics
        """
        fig, ax = plt.subplots(figsize=(8, 6))

        # Extract bootstrap distributions if available
        has_bootstrap = any('bootstrap' in m for m in tool_metrics.values())

        if has_bootstrap:
            data = []
            tools = []
            for tool_name, metrics in tool_metrics.items():
                if 'bootstrap' in metrics:
                    tools.append(tool_name)
                    data.append(metrics['bootstrap'][metric])

            if data:
                ax.violinplot(data, vert=False)
                ax.set_yticks(range(1, len(tools) + 1))
                ax.set_yticklabels(tools)

        ax.set_xlabel(metric.upper())
        ax.set_ylabel('Annotation Tool')

        if title:
            ax.set_title(title)
        else:
            ax.set_title(f'{metric.upper()} Distribution')

        plt.tight_layout()

        if output_file:
            self._save_figure(fig, output_file)

        return fig

    def _save_figure(self, fig, output_file: str):
        """Save figure in multiple formats"""
        output_path = Path(output_file)

        # Ensure output directory exists
        output_path.parent.mkdir(parents=True, exist_ok=True)

        # Save as PNG
        fig.savefig(output_path, dpi=self.style.dpi, bbox_inches='tight')
        logger.info(f"Saved figure: {output_file}")

        # Also save as PDF for publication
        pdf_path = output_path.with_suffix('.pdf')
        fig.savefig(pdf_path, format='pdf', bbox_inches='tight')
        logger.info(f"Saved PDF: {pdf_path}")


def create_summary_figure(tool_metrics: Dict[str, Dict],
                         gene_name: str,
                         output_dir: str,
                         prefix: Optional[str] = None) -> Dict[str, str]:
    """
    Create a complete set of summary figures

    Args:
        tool_metrics: Dictionary of tool performance metrics
        gene_name: Name of target gene (for titles)
        output_dir: Directory to save figures
        prefix: Optional filename prefix

    Returns:
        Dictionary mapping figure types to file paths
    """
    prefix = prefix or f"{gene_name}_performance"
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    visualizer = PerformanceVisualizer()
    figures = {}

    # ROC curves
    roc_file = output_dir / f"{prefix}_roc.png"
    visualizer.plot_roc_curves(tool_metrics,
                              title=f'ROC Curves - {gene_name}',
                              output_file=str(roc_file))
    figures['roc'] = str(roc_file)

    # AUROC comparison
    bar_file = output_dir / f"{prefix}_auroc_bar.png"
    visualizer.plot_auroc_comparison(tool_metrics,
                                    title=f'Annotation Tool Performance - {gene_name}',
                                    output_file=str(bar_file))
    figures['bar'] = str(bar_file)

    # Precision-Recall curves
    pr_file = output_dir / f"{prefix}_pr.png"
    visualizer.plot_precision_recall(tool_metrics,
                                    title=f'Precision-Recall Curves - {gene_name}',
                                    output_file=str(pr_file))
    figures['pr'] = str(pr_file)

    # Heatmap
    heatmap_file = output_dir / f"{prefix}_heatmap.png"
    visualizer.plot_performance_heatmap(tool_metrics,
                                       title=f'Performance Summary - {gene_name}',
                                       output_file=str(heatmap_file))
    figures['heatmap'] = str(heatmap_file)

    logger.info(f"Generated {len(figures)} summary figures")
    return figures