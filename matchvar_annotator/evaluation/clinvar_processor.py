#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ClinVar data processor for extracting benign/pathogenic labels.
"""

import os
import logging
import pandas as pd
import numpy as np
from typing import Dict, List, Optional, Any, Tuple
from datetime import datetime

logger = logging.getLogger(__name__)


class ClinVarProcessor:
    """
    Process ClinVar data to extract binary labels for ROC evaluation.
    
    Labels:
        - benign: label = 0
        - pathogenic/likely_pathogenic: label = 1
        - uncertain/significance omitted (ignored)
    """
    
    BENIGN_TERMS = frozenset([
        'Benign', 'Likely benign', 'Benign/Likely benign'
    ])
    PATHOGENIC_TERMS = frozenset([
        'Pathogenic', 'Likely pathogenic', 'Pathogenic/Likely pathogenic'
    ])
    
    def __init__(self, clinvar_csv_path: str):
        """
        Initialize with path to ClinVar CSV file.
        
        Args:
            clinvar_csv_path: Path to ClinVar data CSV (must have GeneSymbol, 
                            ClinicalSignificance, Chr, Start, ReferenceAlleleVCF, 
                            AlternateAlleleVCF columns)
        """
        self.clinvar_csv_path = clinvar_csv_path
        self._data: Optional[pd.DataFrame] = None
        
    def load_data(self) -> pd.DataFrame:
        """Load ClinVar CSV data."""
        if not os.path.exists(self.clinvar_csv_path):
            raise FileNotFoundError(f"ClinVar file not found: {self.clinvar_csv_path}")
        
        self._data = pd.read_csv(self.clinvar_csv_path)
        logger.info(f"Loaded {len(self._data)} ClinVar records")
        return self._data
    
    def filter_by_gene(self, gene_name: str) -> pd.DataFrame:
        """Filter records by gene symbol."""
        if self._data is None:
            self.load_data()
        
        gene_data = self._data[self._data['GeneSymbol'] == gene_name].copy()
        logger.info(f"Found {len(gene_data)} records for gene {gene_name}")
        return gene_data
    
    def extract_binary_labels(self, df: pd.DataFrame) -> Tuple[pd.DataFrame, np.ndarray]:
        """
        Extract binary labels from ClinVar DataFrame.
        
        Args:
            df: DataFrame with ClinicalSignificance column
            
        Returns:
            Tuple of (filtered_df, labels) where labels are 0 (benign) or 1 (pathogenic)
        """
        import numpy as np
        
        if df.empty:
            return df, np.array([], dtype=np.int8)
        
        # Classify each record
        sig = df['ClinicalSignificance'].fillna('')
        
        is_benign = sig.isin(self.BENIGN_TERMS)
        is_pathogenic = sig.isin(self.PATHOGENIC_TERMS)
        
        # Filter to only clear benign/pathogenic
        mask = is_benign | is_pathogenic
        filtered = df[mask].copy()
        
        labels = np.where(is_benign[mask], 0, 1).astype(np.int8)
        
        logger.info(f"Extracted {len(labels)} clear (benign/pathogenic) labels: "
                   f"{sum(labels==0)} benign, {sum(labels==1)} pathogenic")
        
        return filtered, labels
    
    def get_position_key(self, row: pd.Series) -> str:
        """Create unique position key for matching (chr stripped, Ref/Alt aliases)."""
        chr_val = str(row.get('Chr', row.get('Chromosome', row.get('chromosome', ''))))
        chr_val = chr_val.replace('chr', '').replace('Chr', '')
        if chr_val.upper() == 'MT':
            chr_val = 'M'
        pos = row.get('Start', row.get('Position', row.get('genomic_pos', '')))
        ref = str(row.get('Ref', row.get('ref', row.get('ReferenceAlleleVCF', row.get('original', '')))))
        alt = str(row.get('Alt', row.get('alt', row.get('AlternateAlleleVCF', row.get('mutant', '')))))
        return f"{chr_val}:{pos}:{ref}>{alt}"
    
    def find_matches(
        self, 
        annotated_df: pd.DataFrame,
        clinvar_df: pd.DataFrame
    ) -> List[Dict[str, Any]]:
        """
        Find matching positions between annotated variants and ClinVar.
        """
        matches = []
        
        clinvar_pos_map: Dict[str, pd.Series] = {}
        for _, row in clinvar_df.iterrows():
            key = self.get_position_key(row)
            clinvar_pos_map[key] = row
        
        for _, var in annotated_df.iterrows():
            key = self.get_position_key(var)
            
            if key in clinvar_pos_map:
                clinvar_row = clinvar_pos_map[key]
                sig = clinvar_row.get('ClinicalSignificance', '')
                
                if sig in self.BENIGN_TERMS:
                    label = 0
                elif sig in self.PATHOGENIC_TERMS:
                    label = 1
                else:
                    continue
                
                matches.append({
                    'position': var.get('Start', var.get('genomic_pos', '')),
                    'ref': var.get('Ref', var.get('ref', '')),
                    'alt': var.get('Alt', var.get('alt', '')),
                    'label': label,
                    'clinical_significance': sig,
                    'total_score': var.get('Total_Score', var.get('total_score', 0.0)),
                    'am_pathogenicity': var.get('am_pathogenicity', None),
                })
        
        return matches