#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
MATCHVAR Database Manager
Provides database compression, indexing, validation and other functions
"""

import os
import logging
from typing import List, Dict, Optional, Tuple
from pathlib import Path

logger = logging.getLogger(__name__)

class DatabaseManager:
    """Database manager responsible for database compression, indexing and validation"""
    
    def __init__(self, humandb_dir: str, genome_version: str = "hg19"):
        """
        Initialize database manager
        
        Args:
            humandb_dir: Path to humandb directory
            genome_version: Genome version (hg19, hg38)
        """
        self.humandb_dir = os.path.abspath(humandb_dir)
        self.genome_version = genome_version
        
        # Check if directory exists
        if not os.path.isdir(self.humandb_dir):
            raise FileNotFoundError(f"humandb directory does not exist: {self.humandb_dir}")
        
        logger.info(f"Database manager initialized: {self.humandb_dir}, genome version: {self.genome_version}")
    
    def list_databases(self) -> Dict[str, List[str]]:
        """
        List all database files
        
        Returns:
            Dictionary containing database files in different states
        """
        databases = {
            'txt_files': [],
            'compressed_files': [],
            'indexed_files': [],
            'large_files': []
        }
        
        # 查找txt文件
        pattern = os.path.join(self.humandb_dir, f"{self.genome_version}_*.txt")
        txt_files = [f for f in os.listdir(self.humandb_dir) if f.startswith(f"{self.genome_version}_") and f.endswith('.txt')]
        databases['txt_files'] = [os.path.join(self.humandb_dir, f) for f in txt_files]
        
        # 查找压缩文件
        gz_files = [f for f in os.listdir(self.humandb_dir) if f.startswith(f"{self.genome_version}_") and f.endswith('.txt.gz')]
        databases['compressed_files'] = [os.path.join(self.humandb_dir, f) for f in gz_files]
        
        # 查找已索引文件
        tbi_files = [f for f in os.listdir(self.humandb_dir) if f.startswith(f"{self.genome_version}_") and f.endswith('.txt.gz.tbi')]
        databases['indexed_files'] = [os.path.join(self.humandb_dir, f) for f in tbi_files]
        
        # 查找大文件 (>5GB)
        for txt_file in databases['txt_files']:
            try:
                size_gb = os.path.getsize(txt_file) / (1024**3)
                if size_gb >= 5.0:
                    databases['large_files'].append(txt_file)
            except OSError:
                continue
        
        return databases
    
    def get_database_status(self) -> Dict[str, Dict]:
        """
        Get database status information
        
        Returns:
            Status information for each database file
        """
        status = {}
        databases = self.list_databases()
        
        # Check txt file status
        for txt_file in databases['txt_files']:
            filename = os.path.basename(txt_file)
            gz_file = txt_file + '.gz'
            tbi_file = gz_file + '.tbi'
            
            status[filename] = {
                'txt_exists': os.path.exists(txt_file),
                'compressed': os.path.exists(gz_file),
                'indexed': os.path.exists(tbi_file),
                'txt_size_gb': 0,
                'gz_size_gb': 0,
                'compression_ratio': 0
            }
            
            try:
                if os.path.exists(txt_file):
                    status[filename]['txt_size_gb'] = os.path.getsize(txt_file) / (1024**3)
                if os.path.exists(gz_file):
                    status[filename]['gz_size_gb'] = os.path.getsize(gz_file) / (1024**3)
                    if status[filename]['txt_size_gb'] > 0:
                        status[filename]['compression_ratio'] = (
                            status[filename]['gz_size_gb'] / status[filename]['txt_size_gb']
                        )
            except OSError:
                continue
        
        return status
    
    def build_indexes(self, 
                     min_size_gb: float = 5.0,
                     threads: int = 4,
                     force_rebuild: bool = False,
                     pre_sort: bool = False,
                     verify: bool = True) -> Dict[str, bool]:
        """
        Build database indexes
        
        Args:
            min_size_gb: Minimum file size threshold (GB)
            threads: Number of threads
            force_rebuild: Force rebuild
            pre_sort: Pre-sort
            verify: Verify indexes
            
        Returns:
            Build results for each file
        """
        try:
            from .build_tabix_indexes import main as build_main
            import sys
            import tempfile
            
            # Save original parameters
            original_argv = sys.argv.copy()
            
            # Build command line parameters
            cmd_args = [
                'build_tabix_indexes.py',
                '--humandb', self.humandb_dir,
                '--buildver', self.genome_version,
                '--min-size-gb', str(min_size_gb),
                '--threads', str(threads)
            ]
            
            if force_rebuild:
                cmd_args.append('--force-rebuild')
            if pre_sort:
                cmd_args.append('--pre-sort')
            if verify:
                cmd_args.append('--verify')
            
            # Set command line parameters
            sys.argv = cmd_args
            
            # Execute build
            build_main()
            
            # Restore original parameters
            sys.argv = original_argv
            
            # Return results
            status = self.get_database_status()
            results = {}
            for filename, info in status.items():
                results[filename] = info['indexed']
            
            return results
            
        except ImportError:
            logger.error("pysam not installed, cannot build indexes. Please install: pip install pysam")
            return {}
        except Exception as e:
            logger.error(f"Error building indexes: {e}")
            return {}
    
    def verify_indexes(self) -> Dict[str, bool]:
        """
        Verify all index files
        
        Returns:
            Verification results for each index file
        """
        try:
            from .build_tabix_indexes import verify_tabix_index
            
            results = {}
            databases = self.list_databases()
            
            for gz_file in databases['compressed_files']:
                filename = os.path.basename(gz_file)
                try:
                    results[filename] = verify_tabix_index(gz_file)
                    logger.info(f"Index verification {filename}: {'Success' if results[filename] else 'Failed'}")
                except Exception as e:
                    logger.error(f"Error verifying index {filename}: {e}")
                    results[filename] = False
            
            return results
            
        except ImportError:
            logger.error("pysam not installed, cannot verify indexes. Please install: pip install pysam")
            return {}
    
    def diagnose_indexes(self) -> Dict[str, Dict]:
        """
        Diagnose index issues
        
        Returns:
            Diagnosis results for each index file
        """
        try:
            from .build_tabix_indexes import diagnose_index_issues
            
            results = {}
            databases = self.list_databases()
            
            for gz_file in databases['compressed_files']:
                filename = os.path.basename(gz_file)
                try:
                    results[filename] = diagnose_index_issues(gz_file)
                    logger.info(f"Index diagnosis {filename}: Completed")
                except Exception as e:
                    logger.error(f"Error diagnosing index {filename}: {e}")
                    results[filename] = {'error': str(e)}
            
            return results
            
        except ImportError:
            logger.error("pysam not installed, cannot diagnose indexes. Please install: pip install pysam")
            return {}
    
    def get_compression_stats(self) -> Dict[str, Dict]:
        """
        Get compression statistics
        
        Returns:
            Compression statistics
        """
        stats = {
            'total_files': 0,
            'compressed_files': 0,
            'indexed_files': 0,
            'total_original_size_gb': 0,
            'total_compressed_size_gb': 0,
            'total_compression_ratio': 0,
            'space_saved_gb': 0
        }
        
        status = self.get_database_status()
        
        for filename, info in status.items():
            stats['total_files'] += 1
            
            if info['compressed']:
                stats['compressed_files'] += 1
                stats['total_original_size_gb'] += info['txt_size_gb']
                stats['total_compressed_size_gb'] += info['gz_size_gb']
            
            if info['indexed']:
                stats['indexed_files'] += 1
        
        if stats['total_original_size_gb'] > 0:
            stats['total_compression_ratio'] = (
                stats['total_compressed_size_gb'] / stats['total_original_size_gb']
            )
            stats['space_saved_gb'] = (
                stats['total_original_size_gb'] - stats['total_compressed_size_gb']
            )
        
        return stats
    
    def print_status_report(self):
        """Print status report"""
        print(f"\n=== MATCHVAR Database Status Report ===")
        print(f"Database directory: {self.humandb_dir}")
        print(f"Genome version: {self.genome_version}")
        
        # Basic statistics
        databases = self.list_databases()
        print(f"\nFile statistics:")
        print(f"  TXT files: {len(databases['txt_files'])}")
        print(f"  Compressed files: {len(databases['compressed_files'])}")
        print(f"  Indexed files: {len(databases['indexed_files'])}")
        print(f"  Large files(≥5GB): {len(databases['large_files'])}")
        
        # Compression statistics
        stats = self.get_compression_stats()
        print(f"\nCompression statistics:")
        print(f"  Original total size: {stats['total_original_size_gb']:.2f} GB")
        print(f"  Compressed total size: {stats['total_compressed_size_gb']:.2f} GB")
        print(f"  Compression ratio: {stats['total_compression_ratio']:.2%}")
        print(f"  Space saved: {stats['space_saved_gb']:.2f} GB")
        
        # Detailed status
        status = self.get_database_status()
        print(f"\nDetailed status:")
        for filename, info in status.items():
            status_str = []
            if info['compressed']:
                status_str.append("Compressed")  
            if info['indexed']:
                status_str.append("Indexed")
            
            print(f"  {filename}: {', '.join(status_str) if status_str else 'Not processed'}")
            if info['txt_size_gb'] > 0:
                print(f"    Size: {info['txt_size_gb']:.2f} GB")
                if info['compressed']:
                    print(f"    Compressed size: {info['gz_size_gb']:.2f} GB (Compression ratio: {info['compression_ratio']:.2%})")
        
        print("=" * 50)
