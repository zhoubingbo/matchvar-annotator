#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
MATCHVAR Database Management Command Line Tool
"""

import argparse
import sys
import os
import logging
from typing import Optional

from .database_manager import DatabaseManager

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

def main():
    """Main function"""
    parser = argparse.ArgumentParser(
        description='MATCHVAR Database Management Tool',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Example usage:
  # View database status
  matchvar-db status --humandb /path/to/humandb --buildver hg19
  
  # Build indexes
  matchvar-db index --humandb /path/to/humandb --buildver hg19 --threads 8
  
  # Verify indexes
  matchvar-db verify --humandb /path/to/humandb --buildver hg19
  
  # Diagnose index issues
  matchvar-db diagnose --humandb /path/to/humandb --buildver hg19
        """
    )
    
    # Global parameters
    parser.add_argument(
        '--humandb', '-d',
        default=os.path.join('resources', 'humandb'),
        help='Path to humandb directory'
    )
    
    parser.add_argument(
        '--buildver', '-b',
        default='hg19',
        help='Genome version (hg19, hg38)'
    )
    
    parser.add_argument(
        '--verbose', '-v',
        action='store_true',
        help='Verbose output'
    )
    
    # Subcommands
    subparsers = parser.add_subparsers(dest='command', help='Available commands')
    
    # status command
    status_parser = subparsers.add_parser('status', help='View database status')
    
    # index command
    index_parser = subparsers.add_parser('index', help='Build database indexes')
    index_parser.add_argument(
        '--min-size-gb',
        type=float,
        default=5.0,
        help='Minimum file size threshold (GB)'
    )
    index_parser.add_argument(
        '--threads', '-t',
        type=int,
        default=4,
        help='Number of threads'
    )
    index_parser.add_argument(
        '--force-rebuild',
        action='store_true',
        help='Force rebuild indexes'
    )
    index_parser.add_argument(
        '--pre-sort',
        action='store_true',
        help='Pre-sort'
    )
    index_parser.add_argument(
        '--verify',
        action='store_true',
        help='Verify indexes after building'
    )
    
    # verify command
    verify_parser = subparsers.add_parser('verify', help='Verify indexes')
    
    # diagnose command
    diagnose_parser = subparsers.add_parser('diagnose', help='Diagnose index issues')
    
    # stats command
    stats_parser = subparsers.add_parser('stats', help='Show compression statistics')
    
    args = parser.parse_args()
    
    # Set logging level
    if args.verbose:
        logging.getLogger().setLevel(logging.DEBUG)
    
    # Check command
    if not args.command:
        parser.print_help()
        sys.exit(1)
    
    try:
        # Create database manager
        db_manager = DatabaseManager(args.humandb, args.buildver)
        
        if args.command == 'status':
            db_manager.print_status_report()
            
        elif args.command == 'index':
            logger.info("Starting database index building...")
            results = db_manager.build_indexes(
                min_size_gb=args.min_size_gb,
                threads=args.threads,
                force_rebuild=args.force_rebuild,
                pre_sort=args.pre_sort,
                verify=args.verify
            )
            
            success_count = sum(1 for success in results.values() if success)
            total_count = len(results)
            
            print(f"\nIndex building completed:")
            print(f"  Success: {success_count}/{total_count}")
            print(f"  Failed: {total_count - success_count}/{total_count}")
            
            if args.verify:
                logger.info("Verifying indexes...")
                verify_results = db_manager.verify_indexes()
                verified_count = sum(1 for success in verify_results.values() if success)
                print(f"  Verification passed: {verified_count}/{len(verify_results)}")
            
        elif args.command == 'verify':
            logger.info("Verifying database indexes...")
            results = db_manager.verify_indexes()
            
            success_count = sum(1 for success in results.values() if success)
            total_count = len(results)
            
            print(f"\nIndex verification results:")
            print(f"  Passed: {success_count}/{total_count}")
            print(f"  Failed: {total_count - success_count}/{total_count}")
            
            for filename, success in results.items():
                status = "✅" if success else "❌"
                print(f"  {status} {filename}")
            
        elif args.command == 'diagnose':
            logger.info("Diagnosing database indexes...")
            results = db_manager.diagnose_indexes()
            
            for filename, diagnosis in results.items():
                print(f"\n{filename}:")
                print(f"  File exists: {diagnosis.get('file_exists', False)}")
                print(f"  Index exists: {diagnosis.get('tbi_exists', False)}")
                print(f"  File size: {diagnosis.get('file_size', 0):,} bytes")
                print(f"  Index size: {diagnosis.get('tbi_size', 0):,} bytes")
                print(f"  Can open: {diagnosis.get('can_open', False)}")
                print(f"  Contigs: {diagnosis.get('contigs', [])}")
                
                if 'error' in diagnosis:
                    print(f"  Error: {diagnosis['error']}")
                
                for query in diagnosis.get('sample_queries', []):
                    if query['success']:
                        print(f"  Query {query['query']}: {query['results']} results")
                    else:
                        print(f"  Query {query['query']}: Failed - {query['error']}")
            
        elif args.command == 'stats':
            stats = db_manager.get_compression_stats()
            
            print(f"\n=== Compression Statistics ===")
            print(f"Total files: {stats['total_files']}")
            print(f"Compressed files: {stats['compressed_files']}")
            print(f"Indexed files: {stats['indexed_files']}")
            print(f"Total original size: {stats['total_original_size_gb']:.2f} GB")
            print(f"Total compressed size: {stats['total_compressed_size_gb']:.2f} GB")
            print(f"Compression ratio: {stats['total_compression_ratio']:.2%}")
            print(f"Space saved: {stats['space_saved_gb']:.2f} GB")
            
    except FileNotFoundError as e:
        logger.error(f"File or directory does not exist: {e}")
        sys.exit(1)
    except Exception as e:
        logger.error(f"Error executing command: {e}")
        if args.verbose:
            import traceback
            traceback.print_exc()
        sys.exit(1)

if __name__ == '__main__':
    main()
