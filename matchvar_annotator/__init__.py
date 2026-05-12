#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
MATCHVAR Annotator Package

This is a Python package for functional annotation and analysis of genomic variants.
Provides complete MATCHVAR annotation functionality, including variant annotation, format conversion, coding change analysis, etc.

Main Features:
- Variant Annotation
- Format Conversion
- Coding Change Analysis
- Table Annotation
"""

__version__ = "1.2.0"
__author__ = "Bingbo Zhou"
__email__ = "zhoubingbo@hotmail.com"
__description__ = "MATCHVAR Annotator - Functional annotation and analysis of genomic variants"

# Import main classes
from .matchvar_annotator import MatchvarRunner
from .table_matchvar import TableAnnotator
from .convert2matchvar import Convert2Matchvar
from .coding_change import CodingChange
from .database_manager import DatabaseManager

# Import new modules
from .vsimulator import GeneTranscript, simulate_variants, ExonExtractor
from .pipeline import MatchingPipeline, run_pipeline, run_pipeline_from_args
from .metrics import VariantMetricCalculator, calculate_auroc_from_dataframe
from .visualization import PerformanceVisualizer, create_summary_figure

# Import database indexing functionality
try:
    from .build_tabix_indexes import (
        build_index_for_file,
        verify_tabix_index,
        diagnose_index_issues,
        discover_txt_files,
        discover_gz_files
    )
    TABIX_AVAILABLE = True
except ImportError:
    TABIX_AVAILABLE = False
    # Provide placeholder functions
    def build_index_for_file(*args, **kwargs):
        raise ImportError("pysam not available. Please install: pip install pysam")
    def verify_tabix_index(*args, **kwargs):
        raise ImportError("pysam not available. Please install: pip install pysam")
    def diagnose_index_issues(*args, **kwargs):
        raise ImportError("pysam not available. Please install: pip install pysam")
    def discover_txt_files(*args, **kwargs):
        raise ImportError("pysam not available. Please install: pip install pysam")
    def discover_gz_files(*args, **kwargs):
        raise ImportError("pysam not available. Please install: pip install pysam")

# Define public API
__all__ = [
    # Core classes
    'MatchvarRunner',
    'TableAnnotator',
    'Convert2Matchvar',
    'CodingChange',
    'DatabaseManager',

    # New modules
    'GeneTranscript',
    'simulate_variants',
    'ExonExtractor',
    'MatchingPipeline',
    'run_pipeline',
    'VariantMetricCalculator',
    'calculate_auroc_from_dataframe',
    'PerformanceVisualizer',
    'create_summary_figure',

    # Indexing functions
    'build_index_for_file',
    'verify_tabix_index',
    'diagnose_index_issues',
    'discover_txt_files',
    'discover_gz_files',

    # Version info
    '__version__',
    '__author__',
    '__email__',
    '__description__',

    # Flags
    'TABIX_AVAILABLE',
]

# Package-level configuration
import logging

# Set default logging level
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)

# Get package-level logger
logger = logging.getLogger(__name__)

def get_version():
    """Get package version"""
    return __version__

def get_author():
    """Get author information"""
    return f"{__author__} <{__email__}>"

def get_description():
    """Get package description"""
    return __description__
