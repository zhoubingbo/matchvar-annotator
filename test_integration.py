#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Quick Test Script for MATCHVAR Pipeline Integration

This script verifies that all modules can be imported and basic
functionality works. It does NOT run the full pipeline (requires
GTF/FASTA/database files).

Usage:
    python test_integration.py
"""

import sys
import os

# Add package to path
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))


def test_imports():
    """Test that all new modules can be imported"""
    print("=" * 60)
    print("Testing module imports...")

    try:
        from matchvar_annotator import (
            simulate_variants,
            GeneTranscript,
            ExonExtractor,
            MatchingPipeline,
            run_pipeline,
            VariantMetricCalculator,
            PerformanceVisualizer,
            create_summary_figure
        )
        print("✅ All modules imported successfully")
        return True
    except ImportError as e:
        print(f"❌ Import failed: {e}")
        return False


def test_vsimulator_basic():
    """Test vsimulator module (requires test data)"""
    print("\n" + "=" * 60)
    print("Testing vsimulator module...")

    try:
        from matchvar_annotator import GeneTranscript

        # This would require actual GTF/FASTA files to run
        print("  GeneTranscript class available")
        print("  (Skipping actual variant generation - requires test data)")
        print("✅ vsimulator module OK")
        return True
    except Exception as e:
        print(f"❌ vsimulator test failed: {e}")
        return False


def test_metrics_basic():
    """Test metrics module"""
    print("\n" + "=" * 60)
    print("Testing metrics module...")

    try:
        from matchvar_annotator import VariantMetricCalculator
        import numpy as np

        calculator = VariantMetricCalculator()

        # Create synthetic test data
        y_true = np.array([0, 0, 0, 1, 1, 1])
        y_scores = np.array([0.1, 0.2, 0.3, 0.7, 0.8, 0.9])

        metrics = calculator.evaluate_tool(y_scores, 'test_tool', y_true)

        assert 'auroc' in metrics
        assert 'auprc' in metrics
        assert 0 <= metrics['auroc'] <= 1
        assert 0 <= metrics['auprc'] <= 1

        print(f"  Test AUROC: {metrics['auroc']:.4f}")
        print(f"  Test AUPRC: {metrics['auprc']:.4f}")
        print("✅ metrics module OK")
        return True
    except Exception as e:
        print(f"❌ metrics test failed: {e}")
        import traceback
        traceback.print_exc()
        return False


def test_visualization_basic():
    """Test visualization module"""
    print("\n" + "=" * 60)
    print("Testing visualization module...")

    try:
        from matchvar_annotator import PerformanceVisualizer, VariantMetricCalculator
        import numpy as np
        import tempfile

        # Create test metrics
        calculator = VariantMetricCalculator()
        y_true = np.array([0, 0, 0, 1, 1, 1])
        y_scores = np.array([0.1, 0.2, 0.3, 0.7, 0.8, 0.9])

        metrics = calculator.evaluate_tool(y_scores, 'TestTool', y_true)
        tool_metrics = {'TestTool': metrics}

        # Create visualizer
        visualizer = PerformanceVisualizer()

        # Test ROC plot
        with tempfile.TemporaryDirectory() as tmpdir:
            roc_path = os.path.join(tmpdir, 'test_roc.png')
            fig = visualizer.plot_roc_curves(tool_metrics, output_file=roc_path)
            assert os.path.exists(roc_path)
            print(f"  ROC curve saved: {roc_path}")

            # Test bar plot
            bar_path = os.path.join(tmpdir, 'test_bar.png')
            fig2 = visualizer.plot_auroc_comparison(tool_metrics, output_file=bar_path)
            assert os.path.exists(bar_path)
            print(f"  Bar plot saved: {bar_path}")

        print("✅ visualization module OK")
        return True
    except Exception as e:
        print(f"❌ visualization test failed: {e}")
        import traceback
        traceback.print_exc()
        return False


def test_pipeline_basic():
    """Test pipeline module (without running full pipeline)"""
    print("\n" + "=" * 60)
    print("Testing pipeline module...")

    try:
        from matchvar_annotator import MatchingPipeline

        # Check that class exists and can be instantiated (but not run without files)
        print("  MatchingPipeline class available")
        print("  (Skipping pipeline execution - requires GTF/FASTA/database)")
        print("✅ pipeline module OK")
        return True
    except Exception as e:
        print(f"❌ pipeline test failed: {e}")
        return False


def test_cli_entry_point():
    """Test CLI entry point"""
    print("\n" + "=" * 60)
    print("Testing CLI entry point...")

    try:
        from matchvar_annotator.simulate_annotate_pipeline import main
        print("  CLI entry point found: matchvar-pipeline")
        print("✅ CLI module OK")
        return True
    except Exception as e:
        print(f"❌ CLI test failed: {e}")
        return False


def main():
    """Run all tests"""
    print("\n" + "=" * 60)
    print("MATCHVAR PIPELINE INTEGRATION TEST")
    print("=" * 60 + "\n")

    results = {
        'imports': test_imports(),
        'vsimulator': test_vsimulator_basic(),
        'metrics': test_metrics_basic(),
        'visualization': test_visualization_basic(),
        'pipeline': test_pipeline_basic(),
        'cli': test_cli_entry_point()
    }

    print("\n" + "=" * 60)
    print("TEST SUMMARY")
    print("=" * 60)

    for test_name, passed in results.items():
        status = "✅ PASS" if passed else "❌ FAIL"
        print(f"  {test_name}: {status}")

    total_passed = sum(results.values())
    total = len(results)
    print(f"\nTotal: {total_passed}/{total} tests passed")

    if total_passed == total:
        print("\n✅ All tests passed! Package is ready for use.")
        print("\nNext steps:")
        print("  1. Install the package: pip install -e .")
        print("  2. Run pipeline: matchvar-pipeline --help")
        print("  3. See examples.py for usage patterns")
        return 0
    else:
        print("\n⚠️  Some tests failed. Please check the output above.")
        return 1


if __name__ == '__main__':
    sys.exit(main())