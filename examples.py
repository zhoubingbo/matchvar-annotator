#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
MATCHVAR Pipeline - Quick Start Example

This script demonstrates how to use the integrated matchvar pipeline
in Python code.
"""

import sys
import os

# Add package to path if running from source
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from matchvar_annotator import (
    MatchingPipeline,
    simulate_variants,
    VariantMetricCalculator,
    create_summary_figure
)


def example_basic_pipeline():
    """
    Basic pipeline usage example

    This example shows the end-to-end workflow with minimal configuration.
    """
    print("\n" + "=" * 70)
    print("EXAMPLE 1: Basic Pipeline")
    print("=" * 70 + "\n")

    # Configuration (replace with your actual paths)
    config = {
        'gtf_file': '/path/to/annotation.gtf',
        'fasta_file': '/path/to/genome.fa',
        'gene_name': 'BRCA1',
        'transcript_id': 'NM_007294.4',
        'database_dir': '/path/to/humandb',
        'output_dir': './results/brca1_basic'
    }

    print("Configuration:")
    for k, v in config.items():
        print(f"  {k}: {v}")

    print("\nNote: Replace paths with actual file locations to run.")
    print("The pipeline will:")
    print("  1. Simulate variants from GTF → VCF")
    print("  2. Annotate VCF with external tools")
    print("  3. Calculate auROC scores")
    print("  4. Generate publication-quality figures")


def example_custom_protocols():
    """
    Example with custom annotation protocols

    Shows how to specify which annotation tools to use.
    """
    print("\n" + "=" * 70)
    print("EXAMPLE 2: Custom Protocols")
    print("=" * 70 + "\n")

    pipeline = MatchingPipeline(
        gtf_file='/data/annotation.gtf',
        fasta_file='/data/genome.fa',
        gene_name='TP53',
        transcript_id='NM_000546.6',
        database_dir='/data/humandb',
        output_dir='./results/tp53_custom',
        protocols=['refGene', 'exac03', 'avsift', 'dbnsfp42a'],
        operations=['g', 'f', 'f', 'f'],
        variant_types=['SNV', 'insertion', 'deletion'],
        buildver='hg38',
        threads=8
    )

    print("Pipeline configured with custom protocols:")
    print("  - refGene: Gene annotation")
    print("  - exac03: Population frequency")
    print("  - avsift: SIFT prediction")
    print("  - dbnsfp42a: Multiple prediction scores")


def example_step_by_step():
    """
    Example showing step-by-step control

    Run each component separately for more flexibility.
    """
    print("\n" + "=" * 70)
    print("EXAMPLE 3: Step-by-Step Execution")
    print("=" * 70 + "\n")

    from matchvar_annotator import GeneTranscript, TableAnnotator

    # Step 1: Simulate variants
    print("[Step 1] Simulating variants...")
    transcript = GeneTranscript.from_gtf(
        gene_name='CFTR',
        transcript_id='NM_000492.4',
        gtf_file='/data/gencode.v44.annotation.gtf',
        fasta_file='/data/GRCh38.fa'
    )

    variants = transcript.generate_all_variants()
    vcf_path = './results/cftr_variants.vcf'
    transcript.export_to_vcf(variants, vcf_path)
    print(f"  Saved {sum(len(v) for v in variants.values())} variants to {vcf_path}")

    # Step 2: Annotate
    print("[Step 2] Running annotation...")
    annotator = TableAnnotator(
        queryfile=vcf_path,
        dbloc='/data/humandb',
        outfile='./results/cftr_annotated',
        buildver='hg38',
        protocol='refGene,dbnsfp42a',
        operation='g,f',
        vcfinput=True,
        otherinfo=True
    )
    annotator.run_annotation()
    print(f"  Annotation complete: ./results/cftr_annotated.hg38_multianno.tsv")

    # Step 3: Calculate metrics
    print("[Step 3] Calculating performance metrics...")
    calculator = VariantMetricCalculator()
    df = pd.read_csv('./results/cftr_annotated.hg38_multianno.tsv', sep='\t')

    y_true = calculator.extract_labels_from_annotation(df)
    print(f"  Found {y_true.sum()} pathogenic variants")

    # Step 4: Visualize
    print("[Step 4] Generating figures...")
    scores = {}
    if 'dbnsfp42a_SIFT_score' in df.columns:
        y_scores = pd.to_numeric(df['dbnsfp42a_SIFT_score'], errors='coerce').values
        metrics = calculator.evaluate_tool(y_scores, 'SIFT', y_true)
        scores['SIFT'] = metrics

    figures = create_summary_figure(
        scores,
        gene_name='CFTR',
        output_dir='./results/figures'
    )
    print(f"  Figures saved to ./results/figures/")
    for name, path in figures.items():
        print(f"    {name}: {path}")


def example_batch_processing():
    """
    Example for batch processing multiple genes
    """
    print("\n" + "=" * 70)
    print("EXAMPLE 4: Batch Processing")
    print("=" * 70 + "\n")

    genes = [
        {'name': 'BRCA1', 'transcript': 'NM_007294.4'},
        {'name': 'TP53', 'transcript': 'NM_000546.6'},
        {'name': 'CFTR', 'transcript': 'NM_000492.4'}
    ]

    print("Batch configuration:")
    for gene in genes:
        print(f"  {gene['name']} ({gene['transcript']})")

    print("\nRun loop:")
    print("  for gene in genes:")
    print("      pipeline = MatchingPipeline(")
    print("          gene_name=gene['name'],")
    print("          transcript_id=gene['transcript'],")
    print("          ...")
    print("      pipeline.run()")


def example_evaluation_only():
    """
    Example for evaluating existing annotated file
    """
    print("\n" + "=" * 70)
    print("EXAMPLE 5: Evaluation-Only Mode")
    print("=" * 70 + "\n")

    from matchvar_annotator import VariantMetricCalculator

    annotated_file = './results/existing_annotation.tsv'
    output_prefix = './results/evaluation'

    print(f"Loading annotated file: {annotated_file}")
    df = pd.read_csv(annotated_file, sep='\t')

    calculator = VariantMetricCalculator()

    # Extract ground truth
    y_true = calculator.extract_labels_from_annotation(df)
    print(f"Extracted {len(y_true)} labels ({y_true.sum()} positive)")

    # Identify score columns
    score_columns = [col for col in df.columns if 'score' in col.lower() or 'Score' in col]
    print(f"Found {len(score_columns)} score columns: {score_columns}")

    # Evaluate each tool
    for col in score_columns:
        if col != 'Total_Score':
            scores = pd.to_numeric(df[col], errors='coerce').values
            metrics = calculator.evaluate_tool(scores, col, y_true)

    # Save results
    calculator.save_results(output_prefix)

    print(f"\nResults saved:")
    print(f"  - {output_prefix}_auroc_summary.tsv")
    print(f"  - {output_prefix}_detailed_metrics.json")


def print_cli_usage():
    """Print command-line usage examples"""
    print("\n" + "=" * 70)
    print("CLI USAGE")
    print("=" * 70 + "\n")

    print("Basic command:")
    print("  $ matchvar-pipeline \\")
    print("    --gtf annotation.gtf \\")
    print("    --fasta genome.fa \\")
    print("    --gene BRCA1 \\")
    print("    --transcript NM_007294.4 \\")
    print("    --database humandb \\")
    print("    --output results\n")

    print("With custom protocols:")
    print("  $ matchvar-pipeline \\")
    print("    --gtf gencode.gtf \\")
    print("    --fasta GRCh38.fa \\")
    print("    --gene TP53 \\")
    print("    --transcript NM_000546.6 \\")
    print("    --database humandb \\")
    print("    --output tp53_results \\")
    print("    --protocols refGene,exac03,avsift \\")
    print("    --operations g,f,f \\")
    print("    --buildver hg38 \\")
    print("    --threads 8\n")

    print("List available options:")
    print("  $ matchvar-pipeline --help")
    print("  $ matchvar-pipeline --list-protocols")
    print("  $ matchvar-pipeline --list-variant-types")


def main():
    """Run all examples"""
    print("\n" + "=" * 70)
    print("MATCHVAR PIPELINE - USAGE EXAMPLES")
    print("=" * 70)

    examples = {
        '1': ('Basic pipeline', example_basic_pipeline),
        '2': ('Custom protocols', example_custom_protocols),
        '3': ('Step-by-step execution', example_step_by_step),
        '4': ('Batch processing', example_batch_processing),
        '5': ('Evaluation-only mode', example_evaluation_only),
        'c': ('CLI usage', print_cli_usage)
    }

    print("\nAvailable examples:")
    for key, (desc, _) in examples.items():
        print(f"  {key}. {desc}")
    print("\nRun: python examples.py <number>")
    print("  or: python examples.py all  # to run all examples\n")

    if len(sys.argv) > 1:
        arg = sys.argv[1]
        if arg == 'all':
            for key, (_, func) in examples.items():
                if key != 'c':
                    func()
        elif arg in examples:
            examples[arg][1]()
        else:
            print(f"Unknown example: {arg}")
    else:
        print("No example specified. Use 'python examples.py <number>' or 'all'.")


if __name__ == '__main__':
    main()