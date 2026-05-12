# MATCHVAR Annotator

[![Python Version](https://img.shields.io/badge/python-3.8%2B-blue.svg)](https://python.org)
[![License](https://img.shields.io/badge/license-MIT-green.svg)](LICENSE)
[![PyPI](https://img.shields.io/pypi/v/matchvar-annotator.svg)](https://pypi.org/project/matchvar-annotator/)
[![Web Tool](https://img.shields.io/badge/Web%20Tool-Online-brightgreen.svg)](https://matchvar.intelligene.cn/)

**Version 1.2.0** - Now with integrated variant simulation and automated performance evaluation pipeline!

MATCHVAR Annotator is a comprehensive Python package for functional annotation and analysis of genomic variants. It provides complete MATCHVAR annotation functionality, including variant annotation, format conversion, coding change analysis, and an **all-new integrated pipeline** for variant simulation → annotation → evaluation → visualization.

## 🌐 Web Interface

**Try MATCHVAR Annotator online at: [https://matchvar.intelligene.cn/](https://matchvar.intelligene.cn/)**

### Web Tool Features

- **User-Friendly Interface**: Intuitive web-based interface for easy variant annotation
- **No Installation Required**: Run annotations directly in your browser
- **Multiple Input Formats**: Support for VCF, BED, and other standard formats
- **Real-time Processing**: Fast annotation with progress tracking
- **Interactive Results**: Browse and filter annotation results online
- **Download Options**: Export results in various formats (TSV, CSV, Excel)
- **Batch Processing**: Handle multiple files simultaneously
- **Cloud Storage**: Secure file handling with automatic cleanup

### Web Tool Usage

1. **Upload Your Data**: Drag and drop your variant files (VCF, BED, etc.)
2. **Select Annotation Options**: Choose from available protocols and operations
3. **Configure Parameters**: Set genome version, quality filters, and other options
4. **Run Annotation**: Click to start the annotation process
5. **View Results**: Browse annotated variants with interactive tables
6. **Download Results**: Export your results in your preferred format

### Web Tool Advantages

- **Accessibility**: No technical setup required
- **Scalability**: Handles large datasets efficiently
- **Collaboration**: Share results easily with team members
- **Updates**: Always uses the latest annotation databases
- **Support**: Built-in help and documentation

## Features

### Core Annotation Tools
- **Variant Annotation**: Supports multiple annotation protocols, including gene annotation, region annotation, and filtering operations
- **Format Conversion**: Supports conversion between VCF, BED, MATCHVAR, and other formats
- **Coding Change Analysis**: Analyzes the impact of DNA-level variations on protein sequences
- **Table Annotation**: Provides complete table annotation functionality
- **Database Management**: Supports compression, indexing, and validation of large databases

### 🆕 Integrated Pipeline (New in v1.2.0)
- **Variant Simulation**: Generate synthetic variants from GTF annotations with biologically-grounded impact scores
- **Automated Workflow**: One-command execution from GTF → VCF → annotation → evaluation → visualization
- **Performance Benchmarking**: Automatically calculate auROC scores comparing annotation tools
- **Publication Figures**: Generate ROC curves, auROC comparison plots (PNG + PDF, 300 DPI)
- **End-to-End Automation**: No manual file transfers; automatic path propagation and logging

## Installation

### Install from PyPI (Recommended)

```bash
pip install matchvar-annotator
```

For full functionality including visualization support:

```bash
pip install "matchvar-annotator[visualization]"
```

### Install from Source (Development)

```bash
git clone https://github.com/zhoubingbo/matchvar-annotator.git
cd matchvar-annotator
pip install -e .
```

### Verify Installation

```bash
python -c "import matchvar_annotator; print('Installation successful!')"
matchvar-pipeline --help  # Should display usage information
```

## Quick Start

### 🆕 One-Command Pipeline (Recommended)

The new `matchvar-pipeline` command runs the complete workflow:

```bash
# Basic usage
matchvar-pipeline \
  --gtf /path/to/annotation.gtf \
  --fasta /path/to/genome.fa \
  --gene BRCA1 \
  --transcript NM_007294.4 \
  --database /path/to/humandb \
  --output ./results

# With custom protocols and multi-threading
matchvar-pipeline \
  --gtf gencode.v44.annotation.gtf \
  --fasta GRCh38.primary_assembly.genome.fa \
  --gene TP53 \
  --transcript NM_000546.6 \
  --database humandb_hg38 \
  --output ./tp53_benchmark \
  --protocols refGene,exac03,avsift,dbnsfp42a \
  --operations g,f,f,f \
  --buildver hg38 \
  --threads 16
```

**What it does automatically:**
1. Parses GTF to extract gene structure
2. Generates all variant types (SNV, insertion, deletion, splice sites) with biological scores
3. Exports VCF with HGVS notation and Total_Score
4. Runs annotation via `table_matchvar` using specified protocols
5. Extracts ground truth labels from functional annotations
6. Calculates auROC for each prediction tool
7. Saves publication-ready figures (ROC curves, bar plots, heatmaps) and statistics table

**Output structure:**
```
results/
├── BRCA1_simulated.vcf
├── BRCA1_annotated.hg19_multianno.tsv
├── BRCA1_auroc_statistics.tsv
├── BRCA1_pipeline_summary.json
└── figures/
    ├── BRCA1_performance_roc.png
    ├── BRCA1_performance_roc.pdf
    ├── BRCA1_performance_auroc_bar.png
    └── BRCA1_performance_heatmap.png
```

## Traditional Usage (v1.x)

### Python API Usage

#### Traditional Annotation (v1.x)

```python
from matchvar_annotator import MatchvarRunner

# Create annotator instance
runner = MatchvarRunner(
    resources_dir="/path/to/resources",
    genome_version="hg19",
    thread_count=4
)

# Run annotation
result_df = runner.run_matchvar(
    input_file="variants.vcf",
    protocols=["refGene", "exac03", "avsift"],
    buildver="hg19",
    output_prefix="result"
)

print(f"Annotation completed with {len(result_df)} rows of data")

# Database Management
from matchvar_annotator import DatabaseManager

# Create database manager
db_manager = DatabaseManager(
    humandb_dir="/path/to/humandb",
    genome_version="hg19"
)

# View database status
db_manager.print_status_report()

# Build indexes
results = db_manager.build_indexes(
    min_size_gb=5.0,
    threads=8,
    force_rebuild=False
)

# Verify indexes
verify_results = db_manager.verify_indexes()
```

#### 🆕 Integrated Pipeline API (v1.2.0)

```python
from matchvar_annotator import MatchingPipeline, run_pipeline

# Method 1: Use MatchingPipeline class
pipeline = MatchingPipeline(
    gtf_file="annotation.gtf",
    fasta_file="genome.fa",
    gene_name="BRCA1",
    transcript_id="NM_007294.4",
    database_dir="humandb",
    output_dir="./results",
    protocols=["refGene", "exac03", "avsift"],
    operations=["g", "f", "f"],
    variant_types=["SNV", "insertion", "deletion", "splice_site"],
    buildver="hg19",
    threads=8
)

results = pipeline.run()
print(f"Generated {results['total_variants']} variants")
print(f"Simulated VCF: {results['simulated_vcf']}")
print(f"Annotated TSV: {results['annotated_tsv']}")
print(f"AUROC scores: {results['auroc_scores']}")

# Method 2: Use convenience function
results = run_pipeline(
    gtf_file="annotation.gtf",
    fasta_file="genome.fa",
    gene_name="TP53",
    transcript_id="NM_000546.6",
    database_dir="humandb",
    output_dir="./tp53_results",
    buildver="hg38",
    threads=16
)
```

#### Step-by-Step Control

```python
from matchvar_annotator import GeneTranscript, TableAnnotator, VariantMetricCalculator

# Step 1: Simulate variants
transcript = GeneTranscript.from_gtf(
    gene_name="CFTR",
    transcript_id="NM_000492.4",
    gtf_file="gencode.v44.annotation.gtf",
    fasta_file="GRCh38.fa"
)

variants = transcript.generate_all_variants()
print(f"Generated {sum(len(v) for v in variants.values())} variants")

transcript.export_to_vcf(variants, "cftr_variants.vcf")

# Step 2: Annotate
annotator = TableAnnotator(
    queryfile="cftr_variants.vcf",
    dbloc="humandb",
    outfile="cftr_annotated",
    buildver="hg38",
    protocol="refGene,dbnsfp42a",
    operation="g,f",
    vcfinput=True,
    otherinfo=True
)
annotator.run_annotation()

# Step 3: Evaluate
calculator = VariantMetricCalculator()
df = pd.read_csv("cftr_annotated.hg38_multianno.tsv", sep="\t")

y_true = calculator.extract_labels_from_annotation(df)
print(f"Found {y_true.sum()} pathogenic variants")

# Evaluate each tool's scores
score_columns = [col for col in df.columns if "score" in col.lower() or "Score" in col]
for col in score_columns:
    if col != "Total_Score":
        y_scores = pd.to_numeric(df[col], errors="coerce").values
        metrics = calculator.evaluate_tool(y_scores, col, y_true)
        print(f"{col}: AUROC={metrics['auroc']:.4f}")
```

### Command Line Usage

#### 🆕 Integrated Pipeline (New)

```bash
# Run complete pipeline with one command
matchvar-pipeline \
  --gtf annotation.gtf \
  --fasta genome.fa \
  --gene BRCA1 \
  --transcript NM_007294.4 \
  --database /path/to/humandb \
  --output results

# Customize variant types and protocols
matchvar-pipeline \
  --gtf annotation.gtf \
  --fasta genome.fa \
  --gene TP53 \
  --transcript NM_000546.6 \
  --database humandb \
  --output tp53_results \
  --variant-types SNV,insertion,deletion,splice_site \
  --protocols refGene,exac03,avsift,dbnsfp42a \
  --operations g,f,f,f \
  --buildver hg38 \
  --threads 8 \
  --verbose

# List available options
matchvar-pipeline --help
matchvar-pipeline --list-protocols
matchvar-pipeline --list-variant-types
```

**Key options:**
- `--gtf`: GTF annotation file (required)
- `--fasta`: Reference genome FASTA (required)
- `--gene`: Gene symbol (e.g., BRCA1) (required)
- `--transcript`: Transcript ID (e.g., NM_007294.4) (required)
- `--database`: Path to humandb directory (required)
- `--output`: Output directory (required)
- `--protocols`: Comma-separated annotation tools (default: refGene)
- `--operations`: Comma-separated operations: g=gene, f=filter (default: g)
- `--variant-types`: Types of variants to simulate (default: SNV,insertion,deletion)
- `--buildver`: Genome version, hg19 or hg38 (default: hg19)
- `--threads`: Number of threads (default: 4)
- `--no-visualization`: Skip generating figures (faster)
- `--keep-temp`: Keep intermediate files for debugging

#### Traditional Commands (v1.x)

##### Basic Annotation

```bash
# Use default protocols for annotation
matchvar-annotator input.vcf

# Specify protocols and operations
matchvar-annotator input.vcf --protocol refGene,exac03,avsift --operation g,f,f

# Specify output file
matchvar-annotator input.vcf --output result --protocol refGene,cytoBand,exac03

# Use custom resource directory
matchvar-annotator input.vcf --resources-dir /path/to/resources --genome-version hg38

# Multi-threaded processing
matchvar-annotator input.vcf --threads 8 --protocol refGene,ensGene,knownGene
```

#### Table Annotation Tool

```bash
# Basic table annotation
matchvar-table input.mvinput /path/to/humandb --protocol refGene,cytoBand --operation g,r

# VCF input
matchvar-table input.vcf /path/to/humandb --vcfinput --protocol refGene,exac03,avsift --operation g,f,f

# Specify threads and output format
matchvar-table input.mvinput /path/to/humandb --protocol refGene,exac03 --operation g,f --thread 8 --csvout

# Include additional information
matchvar-table input.mvinput /path/to/humandb --protocol refGene,cytoBand --operation g,r --otherinfo

# Polish gene names
matchvar-table input.mvinput /path/to/humandb --protocol refGene --operation g --polishgene
```

#### Format Conversion Tool

```bash
# VCF4 format conversion
matchvar-convert input.vcf --format vcf4 --outfile output

# Include additional information
matchvar-convert input.vcf --format vcf4 --includeinfo --outfile output

# Set quality thresholds
matchvar-convert input.vcf --format vcf4 --snpqual 20 --snppvalue 0.05 --outfile output

# Coverage filtering
matchvar-convert input.vcf --format vcf4 --coverage 10 --maxcoverage 1000 --outfile output

# Chromosome filtering
matchvar-convert input.vcf --format vcf4 --chr 1,2,3 --outfile output

# Allelic fraction filtering
matchvar-convert input.vcf --format vcf4 --allelicfrac --fraction 0.1 --outfile output
```

#### Coding Change Analysis Tool

```bash
# Basic coding change analysis
matchvar-coding evf.txt gene.txt mrna.fa --alltranscript

# Output to specified file
matchvar-coding evf.txt gene.txt mrna.fa --outfile result --alltranscript
```

#### Database Management Tools

```bash
# View database status
matchvar-db status --humandb /path/to/humandb --buildver hg19

# Build database indexes
matchvar-db index --humandb /path/to/humandb --buildver hg19 --threads 8 --min-size-gb 5.0

# Verify indexes
matchvar-db verify --humandb /path/to/humandb --buildver hg19

# Diagnose index issues
matchvar-db diagnose --humandb /path/to/humandb --buildver hg19

# Show compression statistics
matchvar-db stats --humandb /path/to/humandb --buildver hg19

# Direct use of index building tool
matchvar-index --humandb /path/to/humandb --buildver hg19 --threads 8 --pre-sort --verify
```

## Advanced Usage

### Custom Protocol Configuration

```python
from matchvar_annotator import MatchvarRunner

# Define custom protocol
custom_protocols = {
    "my_clinvar": {
        "file": "/path/to/my_clinvar.txt",
        "operation": "f",
        "description": "Custom ClinVar database"
    }
}

# Create annotator with custom protocols
runner = MatchvarRunner(custom_protocols=custom_protocols)

# Annotate with custom database
result_df = runner.run_matchvar(
    input_file="variants.vcf",
    protocols=["refGene", "my_clinvar"],
    buildver="hg19"
)
```

### Error Handling

```python
from matchvar_annotator import MatchvarRunner
import logging

# Set up logging
logging.basicConfig(level=logging.INFO)

try:
    runner = MatchvarRunner()
    result_df = runner.run_matchvar(
        input_file="variants.vcf",
        protocols=["refGene", "exac03"]
    )
    print("Annotation successful!")
    
except FileNotFoundError as e:
    print(f"File not found: {e}")
    
except ValueError as e:
    print(f"Invalid parameter: {e}")
    
except Exception as e:
    print(f"Unexpected error: {e}")
```

### Performance Optimization

```python
# Use multiple threads for annotation
runner = MatchvarRunner(thread_count=8)
result_df = runner.run_matchvar(
    input_file="variants.vcf",
    protocols=["refGene", "exac03"]
)

# Process large files in chunks
runner = MatchvarRunner()
result_df = runner.run_matchvar(
    input_file="large_variants.vcf",
    protocols=["refGene"],
    chunk_size=10000
)

# Build indexes for better performance
db_manager = DatabaseManager("/path/to/humandb", "hg19")
db_manager.build_indexes(
    min_size_gb=5.0,
    threads=8,
    pre_sort=True
)
```

## Supported Annotation Protocols

### Gene Information Related
- `refGene`: RefSeq gene annotation
- `ensGene`: Ensembl gene annotation  
- `knownGene`: UCSC known gene annotation

### Frequency Databases
- `exac03`: ExAC exome frequency data
- `gnomad211_genome`: gnomAD genome frequency data
- `esp6500siv2_all`: ESP6500 population frequency data
- `1000g2015aug_all`: 1000 Genomes frequency data
- `clinvar`: ClinVar clinical variant database
- `cosmic70`: COSMIC cancer variant database

### Functional Prediction Software
- `avsift`: SIFT functional prediction
- `dbnsfp42a`: dbNSFP functional prediction
- `revel`: REVEL pathogenicity prediction
- `cadd13gt10`: CADD deleteriousness prediction
- `AlphaMissense`: AlphaMissense pathogenicity prediction

### Region Annotation
- `cytoBand`: Cytogenetic band information
- `dgvMerged`: Database of Genomic Variants
- `wgRna`: miRNA and other regulatory elements

## Operation Types

- `g`: Gene annotation
- `r`: Region annotation
- `f`: Filter annotation

## Supported Formats

### Input Formats
- **VCF**: Variant Call Format
- **BED**: Browser Extensible Data
- **MATCHVAR**: MATCHVAR input format
- **Pileup**: Pileup format
- **GFF3-SOLID**: GFF3-SOLID format

### Output Formats
- **TSV**: Tab-separated values
- **CSV**: Comma-separated values
- **MATCHVAR**: Standard MATCHVAR output format

## Configuration

### Environment Variables

```bash
# Specify Python interpreter
export PYTHON_EXECUTABLE=/path/to/python

# Set resource directory
export MATCHVAR_RESOURCES_DIR=/path/to/resources

# Set genome version
export MATCHVAR_GENOME_VERSION=hg19
```

### Configuration File

Create a `config.yaml` file:

```yaml
# Resource configuration
resources_dir: "/path/to/resources"
genome_version: "hg19"
thread_count: 4

# Protocol configuration
protocols:
  - "refGene"
  - "exac03"
  - "avsift"

# Operation configuration
operations:
  - "g"
  - "f"
  - "f"

# Database configuration
database:
  humandb_dir: "/path/to/humandb"
  min_size_gb: 5.0
  threads: 8
```

## Examples

### Example 1: Basic Variant Annotation

```python
from matchvar_annotator import MatchvarRunner

# Create annotator
runner = MatchvarRunner()

# Annotate variants
result_df = runner.run_matchvar(
    input_file="variants.vcf",
    protocols=["refGene", "exac03", "avsift"],
    buildver="hg19"
)

# Save results
result_df.to_csv("annotated_variants.csv", index=False)
print(f"Annotated {len(result_df)} variants")
```

### Example 2: Database Management

```python
from matchvar_annotator import DatabaseManager

# Create database manager
db_manager = DatabaseManager("/path/to/humandb", "hg19")

# Check database status
db_manager.print_status_report()

# Build indexes for large databases
results = db_manager.build_indexes(
    min_size_gb=5.0,
    threads=8,
    force_rebuild=False
)

# Verify indexes
verify_results = db_manager.verify_indexes()
for filename, success in verify_results.items():
    print(f"{filename}: {'OK' if success else 'FAILED'}")
```

## Troubleshooting

### Common Issues

1. **File not found error**
   - Check if the input file exists
   - Verify the file path is correct
   - Ensure proper file permissions

2. **Protocol not found error**
   - Check if the protocol is supported
   - Verify the protocol file exists in the resources directory
   - Check the protocol configuration

3. **Memory error**
   - Reduce the number of threads
   - Process files in smaller chunks
   - Increase system memory

4. **Database indexing error**
   - Check if pysam is installed
   - Verify database file integrity
   - Use pre-sort option for large files

### Debug Mode

```python
import logging

# Enable debug logging
logging.basicConfig(level=logging.DEBUG)

# Run with debug information
runner = MatchvarRunner()
result_df = runner.run_matchvar(
    input_file="variants.vcf",
    protocols=["refGene"],
    debug=True
)
```

## Development

### Install Development Dependencies

```bash
pip install -e ".[dev]"
```

### Run Tests

```bash
pytest tests/
```

### Comprehensive Testing Guide

This section provides detailed testing methods to verify the complete functionality of MATCHVAR Annotator.

#### 1. Package Installation Test

Test the package installation from source:

```bash
# Install the package in development mode
pip install -e .

# Verify installation
python -c "import matchvar_annotator; print(f'Version: {matchvar_annotator.__version__}')"
```

Expected output:
```
Version: 1.1.2
```

#### 2. Package Import Test

Test importing all main classes:

```bash
python -c "
from matchvar_annotator import MatchvarRunner, TableAnnotator, Convert2Matchvar, CodingChange, DatabaseManager
print('All main classes imported successfully')
"
```

Expected output:
```
All main classes imported successfully
```

#### 3. Command Line Tools Test

Test all command line tools:

```bash
# Test main annotator
matchvar-annotator --help

# Test table annotation tool
matchvar-table --help

# Test format conversion tool
matchvar-convert --help

# Test coding change analysis tool
matchvar-coding --help

# Test database management tools
matchvar-db --help
matchvar-index --help
```

#### 4. End-to-End Annotation Test

Perform a complete annotation test using the provided test data:

```bash
# Create result directory
mkdir -p result

# Run comprehensive annotation test
matchvar-table tests/test.vcf resources/humandb \
    -outfile result/test_result \
    -buildver hg19 \
    -protocol refGene,cytoBand \
    -operation g,r \
    -vcfinput \
    -otherinfo \
    -thread 4
```

Expected results:
- Input: `tests/test.vcf` (856KB, ~2854 variants)
- Output: `result/test_result.hg19_multianno.tsv` (~1.2MB, ~2908 rows)
- Processing time: ~8 seconds
- Annotation types: Gene annotation + Region annotation

#### 5. Test Data Verification

Verify the test results:

```bash
# Check output file exists and has correct size
ls -lh result/test_result.hg19_multianno.tsv

# Count output rows
wc -l result/test_result.hg19_multianno.tsv

# Check column structure
head -1 result/test_result.hg19_multianno.tsv | tr '\t' '\n' | nl

# View sample results
head -5 result/test_result.hg19_multianno.tsv
```

Expected output:
- File size: ~1.2MB
- Row count: ~2908 (1 header + 2907 data rows)
- Column count: 27 columns
- Columns include: Chr, Start, End, Ref, Alt, Func.refGene, Gene.refGene, etc.

#### 6. Test File Structure

The test generates the following files in the `result/` directory:

```
result/
├── test_result.hg19_multianno.tsv          # Main annotation results (1.2MB)
├── test_result.mvinput                      # VCF to MV format conversion (993KB)
├── test_result.refGene.variant_function     # Gene variant function (1.1MB)
├── test_result.refGene.exonic_variant_function # Exonic variant function (208KB)
├── test_result.hg19_cytoBand                # Cytogenetic band annotation (1.0MB)
├── test_result.refGene.fa                   # Gene sequences (281KB)
└── test_result.refGene.mane                 # MANE transcript info (372KB)
```

#### 7. Annotation Quality Verification

Verify annotation quality by checking specific examples:

```bash
# Check for exonic variants
grep "exonic" result/test_result.hg19_multianno.tsv | head -3

# Check for missense mutations
grep "missense" result/test_result.hg19_multianno.tsv | head -3

# Check for synonymous variants
grep "synonymous" result/test_result.hg19_multianno.tsv | head -3
```

Expected annotations:
- Exonic variants with HGVS nomenclature (e.g., `p.Pro283Ser`)
- Synonymous variants marked as `p.(=)`
- Intronic variants with proper gene context
- Cytogenetic band information (e.g., `p13.33`)

#### 8. Performance Test

Test performance with different thread counts:

```bash
# Test with 1 thread
time matchvar-table tests/test.vcf resources/humandb \
    -outfile result/test_1thread \
    -buildver hg19 \
    -protocol refGene \
    -operation g \
    -vcfinput \
    -thread 1

# Test with 4 threads
time matchvar-table tests/test.vcf resources/humandb \
    -outfile result/test_4thread \
    -buildver hg19 \
    -protocol refGene \
    -operation g \
    -vcfinput \
    -thread 4
```

#### 9. Database Management Test

Test database management functionality:

```bash
# Check database status
matchvar-db status --humandb resources/humandb --buildver hg19

# Test database indexing (if large files exist)
matchvar-index --humandb resources/humandb --buildver hg19 --threads 4

# Verify indexes
matchvar-db verify --humandb resources/humandb --buildver hg19
```

#### 10. Error Handling Test

Test error handling with invalid inputs:

```bash
# Test with non-existent file
matchvar-table nonexistent.vcf resources/humandb -outfile test_error

# Test with invalid parameters
matchvar-table tests/test.vcf resources/humandb \
    -outfile test_error \
    -buildver invalid_build \
    -protocol nonexistent
```

#### Test Success Criteria

A successful test should meet the following criteria:

✅ **Installation**: Package installs without errors  
✅ **Import**: All main classes import successfully  
✅ **CLI Tools**: All command line tools show help without errors  
✅ **Annotation**: Complete annotation runs without errors  
✅ **Output**: Results file is generated with correct structure  
✅ **Quality**: Annotations contain expected gene and region information  
✅ **Performance**: Processing completes in reasonable time  
✅ **Error Handling**: Invalid inputs are handled gracefully  

#### Troubleshooting Test Issues

If tests fail, check:

1. **Dependencies**: Ensure all required packages are installed
2. **Database**: Verify humandb directory exists and contains required files
3. **Permissions**: Check file read/write permissions
4. **Memory**: Ensure sufficient memory for large file processing
5. **Threads**: Adjust thread count based on system capabilities

For detailed troubleshooting, see the [Troubleshooting](#troubleshooting) section.

### Code Formatting

```bash
black matchvar_annotator/
flake8 matchvar_annotator/
```

## Contributing

Contributions are welcome! Please follow these steps:

1. Fork the repository
2. Create a feature branch (`git checkout -b feature/amazing-feature`)
3. Commit your changes (`git commit -m 'Add some amazing feature'`)
4. Push to the branch (`git push origin feature/amazing-feature`)
5. Open a Pull Request

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Project Structure

```
matchvar-annotator/
├── matchvar_annotator/          # Main package
│   ├── __init__.py              # Package init (v1.2.0 API exports)
│   ├── cli.py                   # Original CLI: matchvar-annotator
│   ├── table_matchvar.py        # Table annotation (unchanged)
│   ├── matchvar_annotator.py    # Core annotator (unchanged)
│   ├── vsimulator.py            # 🆕 Variant simulation wrapper
│   ├── pipeline.py              # 🆕 End-to-end pipeline
│   ├── metrics.py               # 🆕 auROC calculation
│   ├── visualization.py         # 🆕 Figure generation
│   ├── simulate_annotate_pipeline.py  # 🆕 CLI: matchvar-pipeline
│   ├── convert2matchvar.py      # Format conversion (unchanged)
│   ├── coding_change.py         # Coding analysis (unchanged)
│   ├── database_manager.py      # Database utilities (unchanged)
│   └── build_tabix_indexes.py   # Index building (unchanged)
├── setup.py                     # Package setup (updated v1.2.0)
├── requirements.txt             # Dependencies (updated)
├── README.md                    # This file (updated)
├── examples.py                  # 🆕 Usage examples
├── test_integration.py          # 🆕 Integration test
└── LICENSE
```

## Examples

See `examples.py` for complete usage patterns:

```bash
# Run all examples
python examples.py all

# Run specific example
python examples.py 1   # Basic pipeline
python examples.py 2   # Custom protocols
python examples.py 3   # Step-by-step execution
python examples.py 4   # Batch processing
python examples.py 5   # Evaluation-only mode
```

## Testing

```bash
# Run integration tests (no data files required)
python test_integration.py

# Expected output:
# ✅ All modules imported successfully
# ✅ metrics module OK
# ✅ visualization module OK
# ✅ pipeline module OK
# ✅ All tests passed!
```

## Troubleshooting

### "GTF file not found" or "FASTA file not found"
- Verify file paths are correct and accessible
- Gzipped GTF files (`.gz`) are supported
- Check file permissions

### "Database directory missing"
- Ensure `humandb` directory contains annotation files (e.g., `hg19_refGene.txt`)
- Database files should be uncompressed or have `.gz` + `.tbi` indexes

### "No ground truth labels extracted"
- Check that annotation produced `Func.refGene` or `ExonicFunc.refGene` columns
- Variants may be all intergenic; try different gene/transcript
- Ensure simulation generated exonic variants

### Import errors after installation
```bash
# Reinstall in development mode
pip install -e . --force-reinstall

# Check package is importable
python -c "import matchvar_annotator; print(matchvar_annotator.__version__)"
```

### Matplotlib backend errors (headless systems)
```bash
# Use non-interactive backend
export MPLBACKEND=Agg
# Or set in Python before plotting:
import matplotlib
matplotlib.use('Agg')
```

### Memory issues with large genes
- Use `--threads` to control parallelism
- Consider limiting variant types: `--variant-types SNV`
- Pipeline automatically uses temporary files; ensure disk space

## FAQ

**Q: Do I need to modify `data_simulation.py`?**  
A: No. The original file remains untouched. The new `vsimulator` module is a clean wrapper.

**Q: Can I use just the simulation part?**  
A: Yes. Import `GeneTranscript` or `simulate_variants` directly.

**Q: Are old commands still available?**  
A: Yes. `matchvar-annotator`, `matchvar-table`, etc. work exactly as before.

**Q: Does the pipeline require all protocols?**  
A: No. Specify any subset via `--protocols`. Minimal: just `refGene` for gene annotation.

**Q: Can I add my own evaluation metrics?**  
A: Yes. `VariantMetricCalculator` accepts any score array. Extend for precision, F1, etc.

**Q: Where are figures saved?**  
A: In `<output_dir>/figures/` with both `.png` (300 DPI) and `.pdf` formats.

## Contributing

Contributions are welcome! Please follow these steps:

1. Fork the repository
2. Create a feature branch (`git checkout -b feature/amazing-feature`)
3. Commit your changes (`git commit -m 'Add some amazing feature'`)
4. Push to the branch (`git push origin feature/amazing-feature`)
5. Open a Pull Request

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Acknowledgments

- Thanks to all contributors for their support
- Based on the MATCHVAR tool development
- Uses many excellent open source libraries
- Variant scoring algorithm from the original `data_simulation.py`

## Contact

- Author: Bingbo Zhou
- Email: zhoubingbo@hotmail.com
- Project Link: [https://github.com/zhoubingbo/matchvar-annotator](https://github.com/zhoubingbo/matchvar-annotator)
- Web Tool: [https://matchvar.intelligene.cn/](https://matchvar.intelligene.cn/)