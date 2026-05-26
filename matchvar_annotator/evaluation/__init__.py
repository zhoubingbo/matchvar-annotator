#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Evaluation module for variant annotation comparison against ClinVar.

This module provides:
- ClinVar data loading and label extraction (benign/pathogenic)
- ROC/AUPRC calculation comparing Total_Score vs AlphaMissense
- Visualization of ROC curve comparison
"""

from .clinvar_processor import ClinVarProcessor
from .roc_analyzer import ROCAnalyzer
from .plot_generator import PublicationPlotGenerator

__all__ = [
    'ClinVarProcessor',
    'ROCAnalyzer',
    'PublicationPlotGenerator',
]