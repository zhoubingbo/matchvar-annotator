#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
MATCHVAR Annotator Command Line Interface
"""

import argparse
import sys
import os
import logging
from typing import List, Optional

from .matchvar_annotator import MatchvarRunner

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

def main():
    """Main function"""
    parser = argparse.ArgumentParser(
        description='MATCHVAR Annotator - Functional annotation tool for genomic variants',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Example usage:
  # Basic annotation
  matchvar-annotator input.vcf --protocol refGene,exac03 --operation g,f
  
  # Specify output file
  matchvar-annotator input.vcf --output result --protocol refGene,cytoBand,exac03
  
  # Use custom resource directory
  matchvar-annotator input.vcf --resources-dir /path/to/resources --genome-version hg38
  
  # Multi-threaded processing
  matchvar-annotator input.vcf --threads 8 --protocol refGene,ensGene,knownGene
        """
    )
    
    # Required arguments
    parser.add_argument(
        'input_file',
        help='Input variant file (VCF, BED, or MATCHVAR format)'
    )
    
    # Optional arguments
    parser.add_argument(
        '--output', '-o',
        help='Output file prefix (default: matchvar_result)',
        default='matchvar_result'
    )
    
    parser.add_argument(
        '--protocol', '-p',
        help='Comma-separated list of annotation protocols (e.g., refGene,exac03,avsift)',
        default='refGene,cytoBand,exac03,avsift,dbnsfp42a'
    )
    
    parser.add_argument(
        '--operation', '-op',
        help='Comma-separated list of operation types (g=gene annotation, r=region annotation, f=filter)',
        default='g,r,f,f,f'
    )
    
    parser.add_argument(
        '--resources-dir', '-r',
        help='Resource directory path (containing humandb database)',
        default=None
    )
    
    parser.add_argument(
        '--genome-version', '-g',
        help='Genome version (hg19, hg38)',
        default='hg19',
        choices=['hg19', 'hg38']
    )
    
    parser.add_argument(
        '--threads', '-t',
        type=int,
        help='Number of threads',
        default=4
    )
    
    parser.add_argument(
        '--timeout',
        type=int,
        help='Timeout in seconds',
        default=600
    )
    
    parser.add_argument(
        '--use-mane-transcript',
        action='store_true',
        help='Use MANE transcript filtering'
    )
    
    parser.add_argument(
        '--verbose', '-v',
        action='store_true',
        help='Verbose output'
    )
    
    parser.add_argument(
        '--version',
        action='version',
        version='%(prog)s 1.0.0'
    )
    
    args = parser.parse_args()
    
    # Set logging level
    if args.verbose:
        logging.getLogger().setLevel(logging.DEBUG)
    
    try:
        # Check if input file exists
        if not os.path.exists(args.input_file):
            logger.error(f"Input file does not exist: {args.input_file}")
            sys.exit(1)
        
        # Parse protocols and operations
        protocols = [p.strip() for p in args.protocol.split(',')]
        operations = [o.strip() for o in args.operation.split(',')]
        
        # Ensure protocol and operation counts match
        if len(protocols) != len(operations):
            logger.error("Number of protocols and operations must match")
            sys.exit(1)
        
        # Create MATCHVAR runner
        runner = MatchvarRunner(
            resources_dir=args.resources_dir,
            genome_version=args.genome_version,
            thread_count=args.threads
        )
        
        # Prepare additional arguments
        additional_args = {}
        if args.use_mane_transcript:
            additional_args['use_mane_transcript'] = True
        
        logger.info(f"Starting annotation of file: {args.input_file}")
        logger.info(f"Using protocols: {protocols}")
        logger.info(f"Using operations: {operations}")
        logger.info(f"Genome version: {args.genome_version}")
        logger.info(f"Number of threads: {args.threads}")
        
        # Run annotation
        result_df = runner.run_matchvar(
            input_file=os.path.basename(args.input_file),
            protocols=protocols,
            buildver=args.genome_version,
            output_prefix=args.output,
            additional_args=additional_args,
            timeout=args.timeout
        )
        
        if result_df is not None:
            logger.info(f"Annotation completed! Results contain {len(result_df)} rows of data")
            logger.info(f"Output file: {args.output}.{args.genome_version}_multianno.tsv")
        else:
            logger.error("Annotation failed!")
            sys.exit(1)
            
    except KeyboardInterrupt:
        logger.info("User interrupted operation")
        sys.exit(1)
    except Exception as e:
        logger.error(f"An error occurred: {e}")
        if args.verbose:
            import traceback
            traceback.print_exc()
        sys.exit(1)

if __name__ == '__main__':
    main()
