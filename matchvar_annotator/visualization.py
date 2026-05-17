#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Visualization Module for Variant Annotation Evaluation

Creates publication-quality figures:
- ROC / PR curves with AUROC and AUPRC values
- auROC comparison bar plots
- Performance heatmaps
- Score-vs-truth scatter / box plots
- Variant-type distribution charts

Robustness features:
  * Falls back to VCF INFO TOTAL_SCORE when the TSV has no score column
  * Does not hard-skip on all-zero labels; still plots and warns
  * Does not drop scores just because they are less-than-10 non-null
  * Produces diagnostic / no-data placeholder figures rather than failing
"""

import json
import os
import re
import logging
import warnings
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import matplotlib
matplotlib.use('Agg')
import matplotlib.figure as mfigure
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib import colors as mcolors

logger = logging.getLogger(__name__)

# ── Constants ──────────────────────────────────────────────────────────────────
_FUNC_RE    = re.compile(r'splicing|stopgain|stoploss|frameshift|nonsyn',
                         re.IGNORECASE)
_EXFUNC_RE  = re.compile(r'splicing|stopgain|stoploss|frameshift|nonsynonymous',
                         re.IGNORECASE)
_TOTAL_SCORE_COL = 'Total_Score'
_AM_PATH_COL     = 'am_pathogenicity'


# ══════════════════════════════════════════════════════════════════════════════
# Helper – data extraction
# ══════════════════════════════════════════════════════════════════════════════

def extract_score_columns(df: pd.DataFrame) -> List[str]:
    """Return an ordered, priority-ranked list of candidate score columns."""
    # Priority 1 — explicit Total_Score written during variant simulation
    if _TOTAL_SCORE_COL in df.columns:
        series = pd.to_numeric(df[_TOTAL_SCORE_COL], errors='coerce')
        if series.count() > 0:
            return [_TOTAL_SCORE_COL]

    # Priority 2 — AlphaMissense pathogenicity score
    if _AM_PATH_COL in df.columns:
        series = pd.to_numeric(df[_AM_PATH_COL], errors='coerce')
        if series.count() > 0:
            return [_AM_PATH_COL]

    # Priority 3 — any column whose name signals a prediction score
    return [
        c for c in df.columns
        if re.search(r'(?i)score|pathogenicity|prediction|rank', c)
        and c not in (_TOTAL_SCORE_COL, _AM_PATH_COL)
    ]


def extract_true_labels(df: pd.DataFrame) -> Optional[np.ndarray]:
    """
    Extract binary ground-truth labels from an annotated TSV DataFrame.

    Checks, in order of preference:
      Func.refGene  →  ExonicFunc.refGene  →  TYPE

    Returns
    -------
    Optional[np.ndarray]
        Integer array of 0/1 labels, or ``None`` when no suitable column exists.
    """
    try:
        if 'Func.refGene' in df.columns:
            func    = df['Func.refGene'].fillna('').astype(str).str.lower()
            labels  = func.str.contains(_FUNC_RE.pattern,
                                        regex=True, na=False).astype(np.int8)
            n_pos   = int(labels.sum())
            logger.info(
                f"Extracted {n_pos} positive labels from Func.refGene "
                f"({n_pos}/{len(labels)} = {n_pos / max(len(labels), 1):.1%})")
            return labels.values

        if 'ExonicFunc.refGene' in df.columns:
            exonic  = df['ExonicFunc.refGene'].fillna('').astype(str).str.lower()
            labels  = exonic.str.contains(_EXFUNC_RE.pattern,
                                          regex=True, na=False).astype(np.int8)
            n_pos   = int(labels.sum())
            logger.info(
                f"Extracted {n_pos} positive labels from ExonicFunc.refGene "
                f"({n_pos}/{len(labels)} = {n_pos / max(len(labels), 1):.1%})")
            return labels.values

        if 'TYPE' in df.columns:
            var_type = df['TYPE'].fillna('').astype(str).str.upper()
            labels   = var_type.isin(
                ['SPLICE_SITE', 'FRAMESHIFT', 'NONSENSE', 'STOPLOSS']
            ).astype(np.int8)
            n_pos    = int(labels.sum())
            logger.info(
                f"Extracted {n_pos} positive labels from TYPE column "
                f"({n_pos}/{len(labels)} = {n_pos / max(len(labels), 1):.1%})")
            return labels.values

        logger.error(
            "No suitable label column found. "
            "Expected one of: Func.refGene, ExonicFunc.refGene, TYPE")
        return None

    except Exception as exc:
        logger.error(f"Error extracting true labels: {exc}")
        return None


def extract_scores_from_vcf(vcf_path: str) -> Optional[pd.Series]:
    """
    Parse the INFO column of a VCF and return a Series of ``TOTAL_SCORE``
    values keyed by 1-based genomic POS.

    Returns ``None`` when the file is missing or no scores are found.
    """
    if not vcf_path or not os.path.exists(vcf_path):
        return None

    scores: Dict[str, float] = {}
    try:
        with open(vcf_path, 'r') as fh:
            for line in fh:
                if line.startswith('#'):
                    continue
                parts = line.strip().split('\t')
                if len(parts) < 8:
                    continue
                pos  = parts[1]
                info = parts[7]
                m    = re.search(r'TOTAL_SCORE=([0-9.]+)', info)
                if m:
                    try:
                        scores[pos] = float(m.group(1))
                    except ValueError:
                        pass

        if scores:
            logger.info(
                f"Extracted TOTAL_SCORE from {len(scores)} VCF record(s)")
            return pd.Series(scores, name=_TOTAL_SCORE_COL)

    except Exception as exc:
        logger.warning(f"Could not extract TOTAL_SCORE from VCF: {exc}")

    return None


def merge_vcf_scores_into_tsv(tsv_path: str,
                               vcf_path: str) -> pd.DataFrame:
    """
    Inject a ``Total_Score`` column into the TSV by matching VCF TOTAL_SCORE
    values on the genomic ``Start`` position (both 1-based).

    If the TSV already has a ``Total_Score`` column, the file is returned
    unchanged.
    """
    df = pd.read_csv(tsv_path, sep='\t', low_memory=False)

    if _TOTAL_SCORE_COL in df.columns:
        return df

    if 'Start' not in df.columns:
        logger.warning("TSV lacks 'Start' column; cannot merge VCF scores")
        return df

    vcf_scores = extract_scores_from_vcf(vcf_path)
    if vcf_scores is None or vcf_scores.empty:
        return df

    df[_TOTAL_SCORE_COL] = df['Start'].astype(str).map(vcf_scores)
    n = int(df[_TOTAL_SCORE_COL].notna().sum())
    logger.info(f"Injected {n} Total_Score value(s) from VCF into TSV")
    return df


# ══════════════════════════════════════════════════════════════════════════════
# Metric computation (standalone, shared by pipeline.py and visualization.py)
# ══════════════════════════════════════════════════════════════════════════════

def calculate_metrics(
        annotated_tsv: str,
        score_columns: Optional[List[str]] = None,
        label_col: Optional[str] = None,
) -> Tuple[Dict[str, Dict[str, Any]], Optional[np.ndarray], Dict[str, Any]]:
    """
    Read an annotated TSV and compute auROC / auPRC for every score column.

    Parameters
    ----------
    annotated_tsv : str
        Path to the ``_multianno`` TSV.
    score_columns : list[str], optional
        Score columns to evaluate.  When ``None``, columns are auto-detected.
    label_col : str, optional
        Name of the column to use for ground-truth labels.
        When ``None``, auto-detected from ``Func.refGene`` /
        ``ExonicFunc.refGene`` / ``TYPE``.

    Returns
    -------
    metrics : dict
        ``{col: {auroc, auprc, roc_curve, pr_curve, n_variants}}``
    y_true : ndarray or None
        Binary label array, or ``None`` when no label column was found.
    diagnostics : dict
        Intermediate values useful for logging / debugging.
    """
    df = pd.read_csv(annotated_tsv, sep='\t', low_memory=False)
    diagnostics: Dict[str, Any] = {
        'tsv_columns':    list(df.columns),
        'n_rows':         len(df),
        'label_source':   None,
        'n_positive':     0,
        'n_negative':     0,
        'score_columns_found': [],
    }

    # ── labels ───────────────────────────────────────────────────────────────
    y_true = extract_true_labels(df)
    if y_true is not None:
        diagnostics['label_source'] = 'auto-detected'
    else:
        logger.warning("No ground-truth labels found – returning empty metrics")
        return {}, None, diagnostics

    n_pos = int(y_true.sum())
    n_neg = int(len(y_true) - n_pos)
    diagnostics['n_positive']  = n_pos
    diagnostics['n_negative']  = n_neg

    if n_pos == 0:
        logger.warning(
            "All labels are 0 (no positives).  auROC will still be computed "
            "but results may not be meaningful.")
    elif n_pos == len(y_true):
        logger.warning("All labels are 1 – degenerate case.")

    # ── score columns ────────────────────────────────────────────────────────
    if score_columns is None:
        score_columns = extract_score_columns(df)
        diagnostics['score_columns_found'] = score_columns

    if not score_columns:
        logger.warning("No score columns detected in the annotated TSV")
        return {}, y_true, diagnostics

    # ── metrics ──────────────────────────────────────────────────────────────
    from sklearn.metrics import (average_precision_score,
                                 precision_recall_curve,
                                 roc_auc_score, roc_curve)

    metrics: Dict[str, Dict[str, Any]] = {}

    for col in score_columns:
        raw  = pd.to_numeric(df[col], errors='coerce').values
        mask = ~np.isnan(raw)
        n_val = int(mask.sum())

        if n_val < 2:
            logger.warning(f"  '{col}': only {n_val} valid score(s) – skipped")
            continue
        if n_val < 10:
            logger.warning(
                f"  '{col}': {n_val} valid score(s) < 10 – "
                f"computing but results may be unreliable.")

        try:
            # AUROC / AUPRC  ─  sklearn raises ValueError when all labels are 0 or 1
            try:
                auroc_val = float(roc_auc_score(y_true[mask], raw[mask]))
            except ValueError:
                auroc_val = float('nan')

            try:
                auprc_val = float(average_precision_score(y_true[mask], raw[mask]))
            except ValueError:
                auprc_val = float('nan')

            # ROC curve
            try:
                fpr_arr, tpr_arr, _ = roc_curve(y_true[mask], raw[mask])
                roc_data = {'fpr': fpr_arr.tolist(),
                            'tpr': tpr_arr.tolist()}
            except Exception:
                roc_data = {'fpr': [0.0, 1.0], 'tpr': [0.0, 1.0]}

            # PR curve
            try:
                p_arr, r_arr, _ = precision_recall_curve(y_true[mask], raw[mask])
                pr_data = {'precision': p_arr.tolist(),
                           'recall':    r_arr.tolist()}
            except Exception:
                pr_data = {'precision': [], 'recall': []}

            metrics[col] = {
                'auroc':       round(auroc_val, 4)
                if not np.isnan(auroc_val) else float('nan'),
                'auprc':       round(auprc_val, 4)
                if not np.isnan(auprc_val) else float('nan'),
                'roc_curve':   roc_data,
                'pr_curve':    pr_data,
                'n_variants':  n_val,
            }
            logger.info(
                f"  '{col}': auROC={metrics[col]['auroc']:.4f}  "
                f"auPRC={metrics[col]['auprc']:.4f}  n={n_val}")

        except Exception as exc:
            logger.warning(f"  Failed to compute metrics for '{col}': {exc}")

    return metrics, y_true, diagnostics


def make_summary_diagnostics(
        metrics:      Dict[str, Any],
        y_true:       Optional[np.ndarray],
        diagnostics:  Dict[str, Any],
        output_dir:   str,
        gene_name:    str,
) -> Dict[str, Any]:
    """
    Write a JSON diagnostic summary of what happened during metric calculation.
    """
    Path(output_dir).mkdir(parents=True, exist_ok=True)
    summary: Dict[str, Any] = {
        'gene':              gene_name,
        'n_rows':            diagnostics.get('n_rows', 0),
        'n_positive':        diagnostics.get('n_positive',  0),
        'n_negative':        diagnostics.get('n_negative',  0),
        'label_source':      diagnostics.get('label_source', 'unknown'),
        'score_columns':     diagnostics.get('score_columns_found', []),
        'auroc_scores': {
            k: (round(v['auroc'], 4)
                if not isinstance(v.get('auroc'), float)
                   or not np.isnan(v['auroc'])
                else None)
            for k, v in metrics.items()
        },
        'auprc_scores': {
            k: (round(v['auprc'], 4)
                if not isinstance(v.get('auprc'), float)
                   or not np.isnan(v['auprc'])
                else None)
            for k, v in metrics.items()
        },
        'n_valid_per_score': {
            k: v['n_variants'] for k, v in metrics.items()
        },
    }
    json_path = os.path.join(output_dir, f"{gene_name}_diagnostics.json")
    with open(json_path, 'w', encoding='utf-8') as fh:
        json.dump(summary, fh, indent=2, default=str)
    logger.info(f"Diagnostic JSON saved: {json_path}")
    return summary


# ══════════════════════════════════════════════════════════════════════════════
# Visualizer class
# ══════════════════════════════════════════════════════════════════════════════

class FigureStyle:
    """Style configuration for publication figures."""

    dpi               : int  = 300
    figsize_roc       : Tuple[int, int] = (8,  8)
    figsize_bar       : Tuple[int, int] = (10, 6)
    figsize_heatmap   : Tuple[int, int] = (10, 8)
    figsize_hist      : Tuple[int, int] = (10, 6)

    title_fontsize    : int  = 14
    label_fontsize    : int  = 12
    tick_fontsize     : int  = 10
    legend_fontsize   : int  = 10

    colormap_roc      : str  = 'tab10'
    colormap_bar      : str  = 'Blues_d'
    grid_alpha        : float = 0.3


class PerformanceVisualizer:
    """
    Visualize annotation tool performance with publication-quality figures.

    Supported plot types
    --------------------
    * ROC curves                        :meth:`plot_roc_curves`
    * auROC horizontal bar comparison   :meth:`plot_auroc_comparison`
    * Precision-Recall curves           :meth:`plot_precision_recall`
    * Performance heatmap               :meth:`plot_performance_heatmap`
    * Score-vs-ground-truth scatter/box :meth:`plot_score_vs_truth`
    * Variant-type distribution         :meth:`plot_variant_distribution`
    """

    def __init__(self, style: Optional[FigureStyle] = None) -> None:
        self.style = style or FigureStyle()
        self._setup_style()

    # ── Internal helpers ──────────────────────────────────────────────────────

    def _setup_style(self) -> None:
        sns.set_style("whitegrid")
        sns.set_context("paper", font_scale=1.2)
        plt.rcParams.update({
            'figure.dpi':        self.style.dpi,
            'savefig.dpi':       self.style.dpi,
            'savefig.bbox':      'tight',
            'font.family':       'sans-serif',
            'font.sans-serif':   ['Arial', 'DejaVu Sans'],
            'axes.labelsize':    self.style.label_fontsize,
            'axes.titlesize':    self.style.title_fontsize,
            'xtick.labelsize':   self.style.tick_fontsize,
            'ytick.labelsize':   self.style.tick_fontsize,
            'legend.fontsize':   self.style.legend_fontsize,
            'grid.alpha':        self.style.grid_alpha,
        })

    def _save_figure(self, fig: mfigure.Figure, output_file: str) -> None:
        output_path = Path(output_file)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(output_path,        dpi=self.style.dpi, bbox_inches='tight')
        pdf = output_path.with_suffix('.pdf')
        fig.savefig(str(pdf),  format='pdf', bbox_inches='tight')
        logger.info(f"Saved figure: {output_file}")
        logger.info(f"Saved PDF:    {pdf}")

    def _empty_axes_placeholder(self,
                                fig_size: Tuple[int, int],
                                message: str) -> mfigure.Figure:
        """Return a one-subplot figure with a centered placeholder message."""
        fig, ax = plt.subplots(figsize=fig_size)
        ax.text(0.5, 0.5, message, ha='center', va='center',
                fontsize=14, color='gray', transform=ax.transAxes)
        ax.set_xticks([]); ax.set_yticks([])
        return fig

    # ── 1. ROC curves ─────────────────────────────────────────────────────────

    def plot_roc_curves(self,
                        tool_metrics:   Dict[str, Dict],
                        title:          Optional[str]        = None,
                        output_file:    Optional[str]        = None,
                        show_random:    bool                 = True,
                        ) -> mfigure.Figure:
        """Plot ROC curves for multiple annotation tools."""
        fig, ax = plt.subplots(figsize=self.style.figsize_roc)
        colors  = plt.cm.get_cmap(self.style.colormap_roc)
        n_tools = len(tool_metrics)

        if n_tools == 0:
            fig = self._empty_axes_placeholder(
                fig_size=(8, 8),
                message='No score data available\nfor ROC curves')
            fig.suptitle(
                title or 'ROC Curves',
                fontweight='bold', fontsize=self.style.title_fontsize + 2)
            if output_file:
                self._save_figure(fig, output_file)
            return fig

        for idx, (tool_name, metrics) in enumerate(tool_metrics.items()):
            if 'roc_curve' not in metrics:
                continue
            fpr   = np.array(metrics['roc_curve']['fpr'])
            tpr   = np.array(metrics['roc_curve']['tpr'])
            auroc = metrics['auroc']
            ax.plot(fpr, tpr, lw=2.5,
                    label=f"{tool_name}  (AUROC={auroc:.3f})",
                    color=colors(idx % colors.N))

        if show_random:
            ax.plot([0.0, 1.0], [0.0, 1.0], 'k--', lw=1.5, alpha=0.7,
                    label='Random  (AUROC=0.500)')

        ax.set_xlabel('False Positive Rate', fontweight='bold')
        ax.set_ylabel('True Positive Rate',  fontweight='bold')
        ax.set_xlim(0.0, 1.0)
        ax.set_ylim(0.0, 1.05)

        if title:
            ax.set_title(title, fontweight='bold',
                         fontsize=self.style.title_fontsize + 2)
        else:
            ax.set_title('Receiver Operating Characteristic Curves',
                         fontweight='bold',
                         fontsize=self.style.title_fontsize + 2)

        ax.legend(loc='lower right', framealpha=0.9)
        ax.grid(True, alpha=self.style.grid_alpha)
        plt.tight_layout()
        if output_file:
            self._save_figure(fig, output_file)
        return fig

    # ── 2. auROC comparison bar ───────────────────────────────────────────────

    def plot_auroc_comparison(self,
                              tool_metrics:  Dict[str, Dict],
                              title:         Optional[str]  = None,
                              output_file:   Optional[str]  = None,
                              sort_by:       str            = 'auroc',
                              ) -> mfigure.Figure:
        """Horizontal bar plot comparing auROC scores across annotation tools."""
        if not tool_metrics:
            fig = self._empty_axes_placeholder(
                fig_size=(10, 6), message='No metrics available')
            if output_file:
                self._save_figure(fig, output_file)
            return fig

        fig, ax   = plt.subplots(figsize=self.style.figsize_bar)
        tools     = list(tool_metrics.keys())
        au_scores = [tool_metrics[t].get('auroc', 0.0) for t in tools]

        order = np.argsort(au_scores)[::-1]
        tools = [tools[i]       for i in order]
        scores = [au_scores[i] for i in order]

        y_pos    = np.arange(len(tools))
        bar_cmap = plt.cm.get_cmap(self.style.colormap_bar)
        bar_colors = [
            bar_cmap(0.4 + 0.4 * i / max(len(tools) - 1, 1))
            for i in range(len(tools))
        ]

        bars = ax.barh(y_pos, scores, color=bar_colors,
                       edgecolor='black', linewidth=0.8)

        for bar, score in zip(bars, scores):
            ax.text(score + 0.005,
                    bar.get_y() + bar.get_height() / 2,
                    f"{score:.3f}", va='center',
                    fontsize=self.style.tick_fontsize,
                    fontweight='bold')

        ax.set_yticks(y_pos)
        ax.set_yticklabels(tools, fontsize=self.style.label_fontsize)
        ax.set_xlabel('auROC Score', fontweight='bold',
                      fontsize=self.style.label_fontsize)
        ax.set_xlim(0.0, 1.05)
        ax.axvline(x=0.5, color='red', linestyle='--',
                   linewidth=1.5, alpha=0.5, label='Random  (0.50)')
        ax.legend(fontsize=self.style.tick_fontsize)
        ax.grid(True, axis='x', alpha=self.style.grid_alpha)

        if title:
            ax.set_title(title, fontweight='bold',
                         fontsize=self.style.title_fontsize + 2)
        else:
            ax.set_title('Annotation Tool Performance Comparison',
                         fontweight='bold',
                         fontsize=self.style.title_fontsize + 2)

        plt.tight_layout()
        if output_file:
            self._save_figure(fig, output_file)
        return fig

    # ── 3. Precision-Recall curves ───────────────────────────────────────────

    def plot_precision_recall(self,
                              tool_metrics:  Dict[str, Dict],
                              title:         Optional[str] = None,
                              output_file:   Optional[str] = None,
                              ) -> mfigure.Figure:
        """Plot Precision-Recall curves for each annotation tool."""
        fig,   ax   = plt.subplots(figsize=self.style.figsize_roc)
        colors      = plt.cm.get_cmap(self.style.colormap_roc)
        n_tools     = len(tool_metrics)
        plotted     = False

        for idx, (tool_name, metrics) in enumerate(tool_metrics.items()):
            if 'pr_curve' not in metrics:
                continue
            precision = np.array(metrics['pr_curve']['precision'])
            recall    = np.array(metrics['pr_curve']['recall'])
            auprc     = metrics['auprc']
            ax.plot(recall, precision, lw=2.5,
                    label=f"{tool_name}  (AUPRC={auprc:.3f})",
                    color=colors(idx % colors.N))
            plotted = True

        if not plotted:
            fig = self._empty_axes_placeholder(
                fig_size=(8, 8), message='No PR-curve data available')
            fig.suptitle(
                title or 'Precision-Recall Curves',
                fontweight='bold', fontsize=self.style.title_fontsize + 2)
            if output_file:
                self._save_figure(fig, output_file)
            return fig

        ax.set_xlabel('Recall',    fontweight='bold')
        ax.set_ylabel('Precision', fontweight='bold')
        ax.set_xlim(0.0, 1.0)
        ax.set_ylim(0.0, 1.05)

        if title:
            ax.set_title(title, fontweight='bold',
                         fontsize=self.style.title_fontsize + 2)
        else:
            ax.set_title('Precision-Recall Curves',
                         fontweight='bold',
                         fontsize=self.style.title_fontsize + 2)

        ax.legend(loc='upper right', framealpha=0.9)
        ax.grid(True, alpha=self.style.grid_alpha)
        plt.tight_layout()
        if output_file:
            self._save_figure(fig, output_file)
        return fig

    # ── 4. Performance heatmap ────────────────────────────────────────────────

    def plot_performance_heatmap(self,
                                 tool_metrics:   Dict[str, Dict],
                                 metrics_to_plot: Optional[List[str]] = None,
                                 title:          Optional[str]  = None,
                                 output_file:    Optional[str]  = None,
                                 ) -> mfigure.Figure:
        """Plot a heatmap of performance metrics across all annotation tools."""
        metrics_to_plot = metrics_to_plot or ['auroc', 'auprc']

        if not tool_metrics:
            fig = self._empty_axes_placeholder(
                fig_size=(10, 8), message='No metrics data available')
            if output_file:
                self._save_figure(fig, output_file)
            return fig

        data, tools = [], []
        for tool_name, metrics in tool_metrics.items():
            data.append({m: metrics.get(m, float('nan'))
                         for m in metrics_to_plot})
            tools.append(tool_name)

        df_mat = pd.DataFrame(data, index=tools)

        fig, ax = plt.subplots(figsize=self.style.figsize_heatmap)
        sns.heatmap(df_mat.T, annot=True, fmt='.3f', cmap='YlOrRd',
                    cbar_kws={'label': 'Score'},
                    linewidths=0.5, linecolor='gray',
                    vmin=0, vmax=1, ax=ax)

        ax.set_xlabel('Annotation Tool', fontweight='bold')
        ax.set_ylabel('Metric', fontweight='bold')
        fig.suptitle(
            title or 'Performance Metrics Heatmap',
            fontweight='bold', fontsize=self.style.title_fontsize + 2, y=1.01)
        plt.tight_layout()
        if output_file:
            self._save_figure(fig, output_file)
        return fig

    # ── 5. Score-vs-ground-truth scatter/box ─────────────────────────────────

    def plot_score_vs_truth(self,
                            annotated_tsv:  str,
                            gene_name:      str,
                            output_dir:     str,
                            score_col:      str  = _TOTAL_SCORE_COL,
                            ) -> Optional[str]:
        """
        Scatter / box-plot of a prediction score grouped by binary truth label.

        The annotated TSV is read directly.  If ``score_col`` is absent from the
        TSV, a VCF companion file is located automatically and the score is
        injected from the VCF INFO field.

        Returns
        -------
        path to the saved figure, or ``None`` if no usable score data exists.
        """
        df = pd.read_csv(annotated_tsv, sep='\t', low_memory=False)

        # Inject TOTAL_SCORE from VCF if needed
        if score_col not in df.columns:
            # Guess VCF path from annotated TSV path
            guess = annotated_tsv.replace('_annotated.hg19_multianno.tsv',
                                         '_simulated.vcf')\
                                 .replace('_annotated.hg38_multianno.tsv',
                                         '_simulated.vcf')
            df = merge_vcf_scores_into_tsv(annotated_tsv, guess)

        if score_col not in df.columns:
            logger.warning(f"Score column '{score_col}' not found in TSV or VCF")
            return None

        numeric = pd.to_numeric(df[score_col], errors='coerce')
        if numeric.notna().sum() == 0:
            logger.warning(
                f"Column '{score_col}' has no non-NaN numeric values")
            return None

        y_true = extract_true_labels(df)
        if y_true is None:
            logger.warning("Cannot plot score-vs-truth without ground-truth labels")
            return None

        plot_df = pd.DataFrame({
            'score':          numeric.values,
            'ground_truth':   y_true.astype(str),
        }).dropna(subset=['score'])

        if plot_df.empty:
            logger.warning("No non-null score values for score-vs-truth plot")
            return None

        n_cats   = plot_df['ground_truth'].nunique()
        fig, axes = plt.subplots(1, 2, figsize=(14, 6))

        # (a) distribution histogram
        ax = axes[0]
        for lbl, grp in plot_df.groupby('ground_truth'):
            ax.hist(grp['score'], bins=min(20, len(grp)),
                    alpha=0.6, label=f'Label = {lbl}', density=True)
        ax.set_xlabel(score_col,          fontweight='bold')
        ax.set_ylabel('Density',          fontweight='bold')
        ax.set_title(
            f'{score_col}  Distribution by Ground Truth\n[{gene_name}]',
            fontweight='bold', fontsize=self.style.title_fontsize)
        ax.legend(fontsize=self.style.tick_fontsize)
        ax.grid(True, alpha=self.style.grid_alpha)

        # (b) box-plots
        ax = axes[1]
        order = sorted(plot_df['ground_truth'].unique().tolist())
        sns.boxplot(data=plot_df, x='ground_truth', y='score',
                    ax=ax, order=order)
        ax.set_xlabel('Ground Truth Label', fontweight='bold')
        ax.set_ylabel(score_col,           fontweight='bold')
        ax.set_title(
            f'{score_col}  vs Ground Truth  [{gene_name}]',
            fontweight='bold', fontsize=self.style.title_fontsize)
        ax.grid(True, axis='y', alpha=self.style.grid_alpha)

        plt.tight_layout()
        out_path = os.path.join(output_dir,
                                f"{gene_name}_{score_col}_vs_truth.png")
        self._save_figure(fig, out_path)
        return out_path

    # ── 6. Variant-type distribution ──────────────────────────────────────────

    def plot_variant_distribution(self,
                                  annotated_tsv:   str,
                                  gene_name:       str,
                                  output_dir:      str,
                                  type_col:        Optional[str] = None,
                                  ) -> Optional[str]:
        """
        Bar chart + (optionally) pie chart of variant-type counts.

        The type column is auto-detected from the TSV when ``type_col`` is
        ``None``.  Falls back to ``Func.refGene`` or ``ExonicFunc.refGene``.
        """
        df = pd.read_csv(annotated_tsv, sep='\t', low_memory=False)

        if type_col is None:
            for candidate in ['VarType', 'VarType.refGene',
                              'VarType.AlphaMissense']:
                if candidate in df.columns and df[candidate].notna().sum() > 0:
                    type_col = candidate
                    break

        if type_col is None:
            if 'Func.refGene' in df.columns and df['Func.refGene'].notna().sum() > 0:
                type_col = 'Func.refGene'
            elif 'ExonicFunc.refGene' in df.columns \
                    and df['ExonicFunc.refGene'].notna().sum() > 0:
                type_col = 'ExonicFunc.refGene'

        if type_col is None:
            logger.warning("No variant-type column found; skipping distribution plot")
            return None

        counts = df[type_col].fillna('NA').value_counts()
        if counts.empty:
            return None

        n_cats  = len(counts)
        show_pie = n_cats <= 8

        n_cols   = 2 if show_pie else 1
        fig, axes = plt.subplots(1, n_cols,
                                 figsize=(14 if show_pie else 10, 6))
        axes = np.atleast_1d(axes)

        cmap      = plt.cm.get_cmap('Set3')
        bar_colors = [cmap(0.15 + 0.70 * i / max(n_cats - 1, 1))
                      for i in range(n_cats)]

        # Bar chart
        ax = axes[0]  # type: ignore[index]
        counts.plot(kind='bar', ax=ax, color=bar_colors,
                    edgecolor='black', linewidth=0.6)
        ax.set_xlabel(type_col,          fontweight='bold')
        ax.set_ylabel('Count',           fontweight='bold')
        ax.set_title(
            f'Variant Type Distribution  [{gene_name}]\n'
            f'(source: {type_col})',
            fontweight='bold', fontsize=self.style.title_fontsize)
        ax.set_xticklabels(
            ax.get_xticklabels(), rotation=30, ha='right',
            fontsize=self.style.tick_fontsize)
        ax.grid(True, axis='y', alpha=self.style.grid_alpha)
        for container in ax.containers:
            ax.bar_label(container,
                         fontsize=self.style.tick_fontsize - 1, padding=3)

        # Pie chart (only when there are few categories)
        if show_pie:
            ax = axes[1]  # type: ignore[index]
            wedges, texts, autotexts = ax.pie(
                counts.values,
                labels      = counts.index.tolist(),
                autopct     = '%1.1f%%',
                startangle  = 90,
                colors      = bar_colors,
                textprops   = {'fontsize': self.style.tick_fontsize},
            )
            for at in autotexts:
                at.set_fontsize(self.style.tick_fontsize)
            ax.set_title('Proportion by Variant Type',
                         fontweight='bold', fontsize=self.style.title_fontsize)
            ax.axis('equal')

        plt.tight_layout()
        out_path = os.path.join(output_dir,
                                f"{gene_name}_variant_type_distribution.png")
        self._save_figure(fig, out_path)
        return out_path


# ══════════════════════════════════════════════════════════════════════════════
# High-level helper – called by pipeline.py
# ══════════════════════════════════════════════════════════════════════════════

def create_summary_figure(tool_metrics:    Dict[str, Dict],
                          gene_name:       str,
                          output_dir:      str,
                          prefix:          Optional[str]    = None,
                          ) -> Dict[str, str]:
    """
    Create a complete set of publication-quality summary figures.

    Parameters
    ----------
    tool_metrics : dict
        ``{tool_name: {'auroc': …, 'auprc': …, 'roc_curve': …, 'pr_curve': …}}``
        produced by :func:`calculate_metrics`.
    gene_name : str
        Gene symbol (used in figure titles).
    output_dir : str
        Directory to save all produced figures.
    prefix : str, optional
        Base filename without extension
        (default: ``"{gene_name}_performance"``).

    Returns
    -------
    Dict[str, str]
        Map of ``{figure_key → absolute_path}`` for every produced figure.
    """
    prefix   = prefix or f"{gene_name}_performance"
    out_dir  = Path(output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    visualizer = PerformanceVisualizer()
    figures: Dict[str, str] = {}

    # ROC curves
    roc_path = str(out_dir / f"{prefix}_roc.png")
    visualizer.plot_roc_curves(
        tool_metrics,
        title       = f'ROC Curves — {gene_name}',
        output_file = roc_path,
    )
    figures['roc'] = roc_path

    # auROC bar
    bar_path = str(out_dir / f"{prefix}_auroc_bar.png")
    visualizer.plot_auroc_comparison(
        tool_metrics,
        title       = f'Annotation Tool Performance — {gene_name}',
        output_file = bar_path,
    )
    figures['auroc_bar'] = bar_path

    # Precision-Recall curves
    pr_path = str(out_dir / f"{prefix}_pr.png")
    visualizer.plot_precision_recall(
        tool_metrics,
        title       = f'Precision-Recall Curves — {gene_name}',
        output_file = pr_path,
    )
    figures['pr'] = pr_path

    # Heatmap
    hm_path = str(out_dir / f"{prefix}_heatmap.png")
    visualizer.plot_performance_heatmap(
        tool_metrics,
        title       = f'Performance Summary — {gene_name}',
        output_file = hm_path,
    )
    figures['heatmap'] = hm_path

    logger.info(f"Generated {len(figures)} summary figures → {output_dir}")
    return figures


def create_diagnostic_figures(annotated_tsv:   str,
                              gene_name:       str,
                              output_dir:      str,
                              vcf_path:        Optional[str]      = None,
                              score_col:       str                = _TOTAL_SCORE_COL,
                              ) -> Dict[str, Any]:
    """
    Produce exploratory figures directly from an annotated TSV.

    Useful when ground-truth labels are degenerate or when a full evaluation is
    not yet possible.  Does **not** require any score columns to exist;
    diagnostic plots that can be drawn will be drawn, the rest are skipped with
    a warning.

    Parameters
    ----------
    annotated_tsv : str
        Path to the ``_multianno`` TSV.
    gene_name : str
        Gene symbol for figure titles.
    output_dir : str
        Directory to save produced figures.
    vcf_path : str, optional
        Path to the companion VCF.  When provided, ``TOTAL_SCORE`` is
        extracted from the VCF INFO field and injected into the TSV.
    score_col : str
        Name of the score column to plot.

    Returns
    -------
    dict  with keys  ``figures``  (path map) and  ``df``  (augmented DataFrame)
    """
    out_dir = Path(output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(annotated_tsv, sep='\t', low_memory=False)

    # Inject Total_Score from VCF if the column is absent
    if score_col not in df.columns and vcf_path:
        df = merge_vcf_scores_into_tsv(annotated_tsv, vcf_path)

    visualizer = PerformanceVisualizer()

    result: Dict[str, Any] = {'figures': {}, 'df': df}

    # Score-vs-truth
    plot_path = visualizer.plot_score_vs_truth(
        annotated_tsv, gene_name, output_dir, score_col=score_col)
    if plot_path:
        result['figures']['score_vs_truth'] = plot_path

    # Variant-type distribution
    dist_path = visualizer.plot_variant_distribution(
        annotated_tsv, gene_name, output_dir)
    if dist_path:
        result['figures']['variant_type_distribution'] = dist_path

    return result
