#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
MATCHVAR Simulate-Annotate Pipeline CLI

One-command pipeline for end-to-end variant simulation, annotation,
and performance evaluation.

Usage:
    matchvar-pipeline \\
        --gtf annotation.gtf \\
        --fasta genome.fa \\
        --gene BRCA1 \\
        --transcript NM_007294.4 \\
        --database /path/to/humandb \\
        --output_dir ./results

Python API:
    from matchvar_annotator import run_pipeline
    results = run_pipeline(
        gtf_file='...',
        fasta_file='...',
        gene_name='BRCA1',
        transcript_id='NM_007294.4',
        database_dir='...',
        output_dir='...'
    )
"""

import sys
import os
import argparse
import logging

# Add parent directory to path for package imports
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from matchvar_annotator import __version__
from matchvar_annotator.pipeline import MatchingPipeline, run_pipeline

logger = logging.getLogger(__name__)


# ══════════════════════════════════════════════════════════════════════════════
# Module-level helpers
# ══════════════════════════════════════════════════════════════════════════════

def _parse_bool(v: str) -> bool:
    """Parse a command-line string as a boolean value."""
    v = v.strip().lower()
    if v in ('true', 't', 'yes', 'y', '1'):
        return True
    if v in ('false', 'f', 'no', 'n', '0'):
        return False
    raise argparse.ArgumentTypeError(f"Invalid boolean value: {v!r}")


# ══════════════════════════════════════════════════════════════════════════════
# CLI entry point
# ══════════════════════════════════════════════════════════════════════════════

def main():
    """Main CLI entry point."""
    examples = """
Examples:
  # Basic usage with required arguments
  matchvar-pipeline \\
    --gtf /data/annotation.gtf \\
    --fasta /data/genome.fa \\
    --gene BRCA1 \\
    --transcript NM_007294.4 \\
    --database /data/humandb \\
    --output_dir ./brca1_results

  # Customize variant types and protocols
  matchvar-pipeline \\
    --gtf annotation.gtf \\
    --fasta genome.fa \\
    --gene TP53 \\
    --transcript NM_000546.6 \\
    --database humandb \\
    --output_dir tp53_results \\
    --variant-types SNV,insertion,deletion,splice_site \\
    --protocols refGene,exac03,avsift \\
    --operations g,f,f \\
    --threads 8 \\
    --buildver hg38

  # Using gene symbol only (auto-detects transcript)
  matchvar-pipeline \\
    --gtf gencode.v44.annotation.gtf \\
    --fasta GRCh38.primary_assembly.genome.fa \\
    --gene CFTR \\
    --transcript auto \\
    --database humandb_hg38 \\
    --output_dir cftr_results \\
    --buildver hg38
"""

    parser = argparse.ArgumentParser(
        prog='matchvar-pipeline',
        description='MATCHVAR End-to-End Pipeline: Variant Simulation -> Annotation -> Evaluation',
        epilog=examples,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    # ── Required arguments ───────────────────────────────────────────────────
    required = parser.add_argument_group('Required Arguments')
    required.add_argument('--gtf', required=True,
                         help='Path to GTF annotation file (can be gzipped)')
    required.add_argument('--fasta', required=True,
                         help='Path to reference genome FASTA file')
    required.add_argument('--gene', required=True,
                         help='Single gene symbol or comma-separated gene symbols '
                              '(e.g., BRCA1 or BRCA1,TP53)')
    required.add_argument('--transcript', required=True,
                         help='Single transcript or comma-separated transcripts, '
                              'one-to-one with --gene '
                              '(e.g., NM_007294 or NM_007294,NM_000546)')
    required.add_argument('--database', required=True,
                         help='Path to annotation database directory (humandb)')
    required.add_argument('--output_dir', required=True,
                         help='Output directory for all results')
    required.add_argument('--merge_output', nargs='?', const=True, default=False,
                          type=_parse_bool,
                          help="True = merge all genes into one VCF/annotation; "
                               "False = each gene separate (default). "
                               "Usage: --merge_output or --merge_output TRUE/FALSE")

    # ── Optional arguments ───────────────────────────────────────────────────
    optional = parser.add_argument_group('Optional Arguments')
    optional.add_argument('--variant-types', default='SNV,insertion,deletion',
                          help='Comma-separated variant types (default: SNV,insertion,deletion)')
    optional.add_argument('--max-indel-length', type=int, default=15,
                          help='Maximum indel length in bp (default: 15)')
    optional.add_argument('--min-indel-length', type=int, default=1,
                          help='Minimum indel length in bp (default: 1)')
    optional.add_argument('--no-synonymous', action='store_false',
                          dest='synonymous', default=True,
                          help='Exclude synonymous SNVs')
    optional.add_argument('--no-stop-codon', action='store_false',
                          dest='include_stop_codon', default=True,
                          help='Exclude stop-codon variants')
    optional.add_argument('--max-splice-offset', type=int, default=20,
                          help='Maximum intron offset for splice site variants (default: 20)')
    optional.add_argument('--min-splice-offset', type=int, default=1,
                          help='Minimum intron offset for splice site variants (default: 1)')
    optional.add_argument('--no-classic-splice-sites', action='store_false',
                          dest='include_classic_splice_sites', default=True,
                          help='Exclude classic splice sites (+/-1, +/-2)')
    optional.add_argument('--max-variants', type=int, default=None,
                          help='Maximum total variants to generate (default: no limit)')
    optional.add_argument('--protocols', default='refGene',
                          help='Comma-separated annotation protocols (default: refGene)')
    optional.add_argument('--operations', default='g',
                          help='Comma-separated operations: g=gene, f=filter (default: g)')
    optional.add_argument('--buildver', default='hg19',
                          choices=['hg19', 'hg38'],
                          help='Genome version (default: hg19)')
    optional.add_argument('--threads', type=int, default=4,
                          help='Number of threads (default: 4)')
    optional.add_argument('--no-visualization', action='store_true',
                          help='Skip generating figures (faster run)')
    optional.add_argument('--keep-temp', action='store_true',
                          help='Keep temporary files for debugging')

    # ── Logging options ──────────────────────────────────────────────────────
    logging_group = parser.add_argument_group('Logging Options')
    logging_group.add_argument('--verbose', '-v', action='store_true',
                               help='Enable verbose DEBUG output')
    logging_group.add_argument('--log-file',
                               help='Save logs to file (default: auto-generated)')

    # ── Info options ─────────────────────────────────────────────────────────
    info_group = parser.add_argument_group('Information')
    info_group.add_argument('--version', action='version',
                            version=f'%(prog)s {__version__}')
    info_group.add_argument('--list-protocols', action='store_true',
                            help='List available annotation protocols')
    info_group.add_argument('--list-variant-types', action='store_true',
                            help='List supported variant types')

    args = parser.parse_args()

    # ── Info-only requests ────────────────────────────────────────────────────
    if args.list_protocols:
        print("Available annotation protocols:")
        print(" Gene annotation: refGene, ensGene, knownGene")
        print(" Region annotation: cytoBand, gwasCatalog")
        print(" Filter operations: exac03, esp6500si, gnomad, avsift, dbnsfp42a")
        return 0

    if args.list_variant_types:
        print("Supported variant types:")
        print(" SNV           - Single nucleotide variants")
        print(" insertion     - Frameshift insertions")
        print(" deletion      - Frameshift deletions")
        print(" splice_site   - Splice site variants")
        print(" inframe       - In-frame indels (3-6bp)")
        return 0

    # ── Logging setup ─────────────────────────────────────────────────────────
    log_level = logging.DEBUG if args.verbose else logging.INFO
    log_handlers = [logging.StreamHandler(sys.stdout)]

    if args.log_file:
        log_handlers.append(
            logging.FileHandler(args.log_file, encoding='utf-8'))

    logging.basicConfig(
        level=log_level,
        format='%(asctime)s - %(levelname)s - %(message)s',
        handlers=log_handlers,
    )

    # ── Input validation ─────────────────────────────────────────────────────
    try:
        if not os.path.exists(args.gtf):
            logger.error(f"GTF file not found: {args.gtf}")
            return 1
        if not os.path.exists(args.fasta):
            logger.error(f"FASTA file not found: {args.fasta}")
            return 1
        if not os.path.isdir(args.database):
            logger.error(f" Database directory not found: {args.database}")
            return 1
    except Exception as exc:
        logger.error(f" Input validation failed: {exc}")
        return 1

    try:
        logger.info("=" * 70)
        logger.info(f" MATCHVAR Pipeline v{__version__}")
        logger.info("=" * 70)

        # Log configuration summary
        logger.info("Pipeline configuration:")
        logger.info(f" GTF file:    {args.gtf}")
        logger.info(f" FASTA file:  {args.fasta}")
        logger.info(f" Gene:        {args.gene}")
        logger.info(f" Transcript:  {args.transcript}")
        logger.info(f" Database:    {args.database}")
        logger.info(f" Output dir:  {args.output_dir}")
        logger.info(f" Variant types: {args.variant_types}")
        logger.info(f" Protocols:   {args.protocols}")
        logger.info(f" Operations:  {args.operations}")
        logger.info(f" Buildver:    {args.buildver}")
        logger.info(f" Threads:     {args.threads}")
        logger.info(f" Max indel:   {args.max_indel_length}")
        logger.info(f" Min indel:   {args.min_indel_length}")
        logger.info(f" Max splice:  {args.max_splice_offset}")
        logger.info(f" Min splice:  {args.min_splice_offset}")
        logger.info(f" Merge output: {args.merge_output}")

        logger.info("=" * 70)
        logger.info(" Starting pipeline execution...")
        logger.info("=" * 70)

        results = run_pipeline(
            gtf_file=args.gtf,
            fasta_file=args.fasta,
            gene_name=args.gene,
            transcript_id=args.transcript,
            database_dir=args.database,
            output_dir=args.output_dir,
            merge_output=args.merge_output,
            variant_types=[v.strip() for v in args.variant_types.split(',')],
            protocols=[p.strip() for p in args.protocols.split(',')],
            operations=[o.strip() for o in args.operations.split(',')],
            buildver=args.buildver,
            threads=args.threads,
            max_indel_length=args.max_indel_length,
            min_indel_length=args.min_indel_length,
            synonymous=args.synonymous,
            include_stop_codon=args.include_stop_codon,
            max_splice_offset=args.max_splice_offset,
            min_splice_offset=args.min_splice_offset,
            include_classic_splice_sites=args.include_classic_splice_sites,
            max_variants=args.max_variants,
        )

        # ── Summary ───────────────────────────────────────────────────────────
        logger.info("=" * 70)
        logger.info(" Pipeline Execution Summary")
        logger.info("=" * 70)

        if results.get('mode') == 'separate_output' or results.get('all_gene_results'):
            all_gene_results = results.get('all_gene_results', {})
            logger.info(f" Mode: multi-gene ({len(all_gene_results)} gene(s))")
            for gene, gene_res in all_gene_results.items():
                logger.info(f" Gene: {gene}")
                logger.info(f" Total variants: {gene_res.get('total_variants', 'N/A')}")
                logger.info(f" Simulated VCF: {gene_res.get('simulated_vcf', 'N/A')}")
                logger.info(f" Annotated TSV: {gene_res.get('annotated_tsv', 'N/A')}")
                if gene_res.get('auroc_scores'):
                    for tool, metrics in gene_res['auroc_scores'].items():
                        logger.info(f"    {tool} auROC: {metrics['auroc']:.4f}")
            if results.get('simulated_vcf'):
                logger.info(f" Merged VCF: {results['simulated_vcf']}")
            if results.get('annotated_tsv'):
                logger.info(f" Merged annotation: {results['annotated_tsv']}")
            if results.get('auroc_scores'):
                logger.info(" Merged auROC Scores:")
                for tool, metrics in results['auroc_scores'].items():
                    logger.info(f" {tool}: {metrics['auroc']:.4f}")
        else:
            logger.info(
                f" Gene: {results.get('gene_name', 'N/A')} "
                f"({results.get('transcript_id', 'N/A')})")
            logger.info(
                f" Total variants simulated: {results.get('total_variants', 'N/A')}")
            logger.info(f" Simulated VCF: {results.get('simulated_vcf', 'N/A')}")
            logger.info(f" Annotated TSV: {results.get('annotated_tsv', 'N/A')}")
            logger.info(f" Output directory: {args.output_dir}")

            if results.get('auroc_scores'):
                logger.info("\nauROC Scores:")
                for tool, metrics in results['auroc_scores'].items():
                    logger.info(f"  {tool}: {metrics['auroc']:.4f}")

        logger.info(" Generated files:")
        for file_type, file_path in results.get('figures', {}).items():
            logger.info(f" {file_type}: {file_path}")

        logger.info(" Pipeline completed successfully!")
        logger.info("=" * 70)

        return 0

    except KeyboardInterrupt:
        logger.info(" Operation cancelled by user")
        return 130
    except Exception as exc:
        logger.error(f" Pipeline failed: {exc}")
        if args.verbose:
            import traceback
            logger.error(traceback.format_exc())
        return 1


if __name__ == '__main__':
    sys.exit(main())
