# MATCHVAR Annotator

[![Python Version](https://img.shields.io/badge/python-3.8%2B-blue.svg)](https://python.org)
[![License](https://img.shields.io/badge/license-MIT-green.svg)](LICENSE)
[![PyPI](https://img.shields.io/pypi/v/matchvar-annotator.svg)](https://pypi.org/project/matchvar-annotator/)
[![Web Tool](https://img.shields.io/badge/Web%20Tool-Online-brightgreen.svg)](https://matchvar.intelligene.cn/)

**Version 1.2.0** — gene/region/filter annotation, HGVS polish, MANE filtering, and an optional simulate → annotate → AUROC pipeline.

Online UI (no install): [https://matchvar.intelligene.cn/](https://matchvar.intelligene.cn/)

## What this package does

1. Annotate VCF (classic 4 columns `CHROM POS REF ALT`, or a full VCF) and MATCHVAR 5-column `.mvinput` against a local `humandb`.
2. Optionally convert pileup / GFF3 / VCF → MATCHVAR `.mvinput` (`matchvar-convert`) when you still want that table. Annotation itself does **not** require conversion.
3. Polish coding variants to c.HGVS / p.HGVS (including HGVS 3′ indel normalization).
4. Optionally simulate variants from GTF + FASTA, annotate them, and score tools with AUROC.

Gene-annotation columns follow the MATCHVAR contract (aligned with the web tool):

| Column | Meaning |
|--------|---------|
| `Chr Start End Ref Alt` | Variant coordinates (1-based closed) |
| `gHGVS` | Genomic HGVS, e.g. `chr17:g.41276045A>G` |
| `Function.{db}` | Region: exonic, intronic, UTR5, splicing, … |
| `Gene.{db}` | Gene symbol |
| `cHGVS.{db}` | Transcript / c. notation |
| `MANE Select.{db}` | `yes` / `no` per transcript fragment |
| `ExonicEffect.{db}` / `VarType.{db}` | Unified type: missense, frameshift, RV, SNV, INDEL, … |
| `pHGVS.{db}` | Protein HGVS |

`{db}` is `refGene`, `ensGene`, or `knownGene`. Older names such as `Func.refGene` / `GeneDetail.refGene` are still recognized when reading results.

## Installation

```bash
pip install matchvar-annotator
# or from source
git clone https://github.com/zhoubingbo/matchvar-annotator.git
cd matchvar-annotator
pip install -e .
```

Visualization, AUROC, FASTA/GTF simulation, and tabix indexing are included in the main install (`matplotlib`, `scikit-learn`, `biopython`, `pyfaidx`, `pysam`).

```bash
python -c "import matchvar_annotator; print(matchvar_annotator.__version__)"
matchvar-annotator --help
matchvar-pipeline --help
```

You need a local annotation database directory (`humandb`). Classic ANNOVAR-style files still work (`hg19_refGene.txt` / `refGene.txt`, `hg19_refGeneMrna.fa`, a genome FASTA such as `Homo_sapiens_assembly37.fasta`). The same directory can also hold UCSC/IGV files and they will be picked up automatically:

| File | Role |
|------|------|
| `ncbiRefSeq.bb` | NCBI RefSeq gene models (protocol `ncbiRefSeq`, and fallback for `refGene` if the txt is missing) |
| `gencodeV49lift37.bb` / `gencode*.bb` | GENCODE models (protocol `gencode`; fallback for `ensGene`) |
| `mane.bb` | MANE Select / Plus Clinical mapping (needs `pip install pyBigWig`) |
| `{buildver}.2bit` | Genome sequence when FASTA is absent (used to splice mRNA for polish and to fill omitted g.HGVS alleles) |
| `{buildver}.chrom.sizes` | Chromosome lengths for bigBed / coordinate checks |
| `cytoBand.txt.gz` | Cytoband region annotation (`cytoBand`) |

Gzipped files with `.tbi` indexes are supported. Reading `*.bb` requires `pyBigWig` (`pip install matchvar-annotator[bb]` or `pip install pyBigWig`).

## Command-line tools

| Command | Role |
|---------|------|
| `matchvar-annotator` | High-level annotation (VCF path → `MatchvarRunner`) |
| `matchvar-table` | Direct table annotation against `humandb` |
| `matchvar-convert` | VCF / pileup / GFF3 → `.mvinput` |
| `matchvar-coding` | EVF → c./p. HGVS polish |
| `matchvar-ghgvs` | g.HGVS file or `--hgvs` → annotation TSV |
| `matchvar-db` / `matchvar-index` | Database status, bgzip + tabix |

### 1. Annotate a VCF (`matchvar-annotator`)

Pass a **real file path** (absolute or relative). `--operation` must match `--protocol` one-to-one (`g` = gene, `r` = region, `f` = filter).

```bash
# Gene + region + a few filter databases
matchvar-annotator /data/sample.vcf \
  --resources-dir /path/to/resources \
  --genome-version hg38 \
  --protocol refGene,cytoBand,exac03,avsift,dbnsfp42a \
  --operation g,r,f,f,f \
  --output result \
  --threads 8

# MANE Select / Plus Clinical filtering
matchvar-annotator /data/sample.vcf \
  --resources-dir /path/to/resources \
  --genome-version hg38 \
  --protocol refGene \
  --operation g \
  --use-mane-transcript \
  --output result_mane
```

`--resources-dir` should contain `humandb/` (and optionally `custom_databases/`, `mane/`). If omitted, the runner looks at `MATCHVAR_RESOURCES_DIR`, then `./resources` next to the package or the current working directory.

Default for `.vcf` is native `CHROM POS REF ALT` annotation. Add `--convert-vcf` to run `convert2matchvar` first (same as `matchvar-table -convertvcf`).

Output: `{prefix}.{hg19|hg38}_multianno.tsv`.

### 2. Table annotation (`matchvar-table`)

Pass a VCF, a 4-column `CHROM POS REF ALT` TSV, or a 5-column `.mvinput`. Both annotation paths are supported:

| Path | When | What happens |
|------|------|----------------|
| Native (default) | `.vcf` / `.vcf.gz` / 4-column TSV | Parse `CHROM POS REF ALT` in the annotator. No conversion. |
| Convert then annotate | `-convertvcf` | `convert2matchvar` writes `{outfile}.mvinput`, then annotate that 5-column file (ANNOVAR-style, keeps INFO/samples via `-includeinfo`). |
| Already `.mvinput` | 5-column `Chr Start End Ref Alt` | Annotate as-is. `-convertvcf` is ignored. |

VCF alleles on the native path are mapped internally (insertion `A AAT` → `start==end`, Ref `-`; deletion `ATG A` → deleted bases with Alt `-`).

```bash
# Native VCF (or 4-column CHROM POS REF ALT) — default
matchvar-table sample.vcf /path/to/humandb \
  -outfile result \
  -buildver hg19 \
  -protocol refGene,cytoBand \
  -operation g,r \
  -thread 8

# Legacy: convert VCF → .mvinput, then annotate
matchvar-table sample.vcf /path/to/humandb \
  -convertvcf \
  -outfile result \
  -buildver hg19 \
  -protocol refGene -operation g

# MATCHVAR 5-column input
matchvar-table input.mvinput /path/to/humandb \
  -outfile result \
  -buildver hg38 \
  -protocol refGene,cytoBand,exac03 \
  -operation g,r,f \
  -otherinfo \
  -thread 8

# Force native VCF parsing when the file has no .vcf suffix
matchvar-table input.txt /path/to/humandb \
  -vcfinput \
  -outfile result -buildver hg38 -protocol refGene -operation g

# Disable coding-change polish (on by default)
matchvar-table input.mvinput /path/to/humandb \
  -outfile result -protocol refGene -operation g -nopolish
```

`-polishgene` is accepted for compatibility; polish is **on by default**. Use `-nopolish` to turn it off.

### 3. Format conversion (`matchvar-convert`)

Use this when you want a standalone 5-column `.mvinput` (or pileup / GFF3). For annotation you can either convert first, or pass `-convertvcf` to `matchvar-table` so conversion and annotation run together.

```bash
matchvar-convert input.vcf --format vcf4 --outfile sample
matchvar-convert input.vcf --format vcf4 --includeinfo --outfile sample
matchvar-convert input.vcf --format vcf4 --chr 1,2,17 --outfile sample
```

`--snpqual` is only for `pileup` / `vcf4old`. Do not pass it for `vcf4`.

Python API (vcf4 no longer requires quality filters):

```python
from matchvar_annotator import Convert2Matchvar

conv = Convert2Matchvar("input.vcf", format="vcf4", outfile="sample")
conv.convert()  # writes sample.mvinput depending on CLI/outfile usage
```

### 4. Annotate g.HGVS (`matchvar-ghgvs`)

Input is one g.HGVS per line (or `--hgvs` on the command line). The tool converts to MATCHVAR `.mvinput`, then writes the usual `*.{hg19|hg38}_multianno.tsv`.

```bash
# File of g.HGVS strings
matchvar-ghgvs variants.ghgvs --resources-dir /path/to/resources --output result

# Single variant
matchvar-ghgvs --hgvs 'chr5:g.131008083_131008084insTGACAGTTGTTTGCACAG' \
  --humandb /path/to/humandb --output fnip1 --buildver hg19

# Choose protocols (otherwise auto-detect refGene / ncbiRefSeq / gencode / cytoBand)
matchvar-ghgvs variants.ghgvs --humandb humandb -p refGene,cytoBand -op g,r -o out

# Only MANE Select / MANE Plus Clinical transcripts (cHGVS / pHGVS)
matchvar-ghgvs --hgvs 'chr10:g.114925413_114925426del' \
  --humandb /path/to/humandb --buildver hg19 \
  --use-mane-transcript -o tcf7l2_mane

# Explicit mapping file (GTF, mane.bb, or gene/refseq/ensembl TSV)
matchvar-ghgvs variants.ghgvs --humandb /path/to/humandb \
  --mane-file /path/to/mane.bb --use-mane-transcript -o result
```

`--use-mane-transcript` keeps MANE Select and MANE Plus Clinical in `cHGVS` / `pHGVS`. Mapping is auto-discovered from `humandb/mane.bb`, `resources/mane/*.gtf`, or `mane_transcript.txt`. Genes with no MANE hit keep the original transcripts. `MANE Select.{db}` still marks each fragment `yes`/`no`.

Python:

```python
from matchvar_annotator import annotate_g_hgvs, parse_g_hgvs

parse_g_hgvs("chr5:g.131008083_131008084delCC")
# {'chrom': 'chr5', 'start': 131008083, 'end': 131008084, 'ref': 'CC', 'alt': '-', ...}

df = annotate_g_hgvs(
    ["chr5:g.131008083_131008084insTGACAGTTGTTTGCACAG"],
    dbloc="resources/humandb",
    buildver="hg19",
    outfile="fnip1",
    use_mane_transcript=True,
)
```

Omitted deletion/duplication sequence (`g.123_125del`) is filled from FASTA or `*.2bit` when present.

### 5. End-to-end pipeline (`matchvar-pipeline`)

The output flag is **`--output_dir`**, not `--output`.

```bash
matchvar-pipeline \
  --gtf /data/gencode.v44.annotation.gtf \
  --fasta /data/GRCh38.fa \
  --gene BRCA1 \
  --transcript NM_007294.4 \
  --database /path/to/humandb \
  --output_dir ./results \
  --buildver hg38 \
  --protocols refGene,exac03,avsift,dbnsfp42a \
  --operations g,f,f,f \
  --variant-types SNV,insertion,deletion,splice_site \
  --threads 8

# Faster: skip figures; keep temp files for debugging
matchvar-pipeline \
  --gtf annotation.gtf --fasta genome.fa \
  --gene TP53 --transcript NM_000546.6 \
  --database humandb --output_dir ./tp53 \
  --no-visualization --keep-temp

# Optional ClinVar ROC (CSV with Chr/Start/Ref/Alt or VCF allele columns)
matchvar-pipeline \
  --gtf annotation.gtf --fasta genome.fa \
  --gene BRCA1 --transcript NM_007294.4 \
  --database humandb --output_dir ./brca1 \
  --clinvar /path/to/clinvar.csv

matchvar-pipeline --list-protocols
matchvar-pipeline --list-variant-types
```

`--simulation-strategy spectrum` draws SNVs from a trinucleotide mutation-rate model (default Alexandrov-style rates; optional gnomAD constraint scores under `resources/gnomad/`). `--simulation-strategy hybrid` mixes those SNVs with traditional SNV/splice enumeration (`--spectrum-ratio`, default `0.7`). Both use the bundled `EnhancedGeneTranscript`. If generation yields no variants, the pipeline logs a warning and falls back to `traditional`.

Typical output:

```
results/
├── BRCA1_simulated.vcf
├── BRCA1_annotated.hg38_multianno.tsv
├── BRCA1_auroc_statistics.tsv
├── BRCA1_pipeline_summary.json
└── figures/
    ├── BRCA1_performance_roc.png
    └── BRCA1_performance_roc.pdf
```

Simulated deletions are written as left-aligned VCF (`REF` = anchor + deleted bases, `ALT` = anchor). Frameshift INFO is `FRAMESHIFT=true` or `FRAMESHIFT=false`.

### 6. Database indexes

```bash
matchvar-db status  --humandb /path/to/humandb --buildver hg19
matchvar-db index   --humandb /path/to/humandb --buildver hg19 --threads 8 --min-size-gb 5.0
matchvar-db verify  --humandb /path/to/humandb --buildver hg19
matchvar-index --humandb /path/to/humandb --buildver hg19 --threads 8 --pre-sort --verify
```

## Python API

### Annotate

```python
from matchvar_annotator import MatchvarRunner

runner = MatchvarRunner(
    resources_dir="/path/to/resources",  # contains humandb/
    genome_version="hg38",
    thread_count=8,
)

df = runner.run_matchvar(
    input_file="/data/sample.vcf",  # real path, not basename
    protocols=["refGene", "cytoBand", "exac03"],
    buildver="hg38",
    output_prefix="result",
    additional_args={
        "operations": ["g", "r", "f"],
        "use_mane_transcript": True,
    },
)
print(df.columns.tolist())
# ..., gHGVS, Function.refGene, Gene.refGene, cHGVS.refGene, MANE Select.refGene, ...
```

Custom filter tables named `{buildver}_{dbname}.txt` under `resources/custom_databases/` (or `humandb/`) are picked up automatically as `operation=f`.

### Table annotator

```python
from matchvar_annotator import TableAnnotator

annotator = TableAnnotator(
    queryfile="sample.vcf",
    dbloc="/path/to/humandb",
    outfile="sample_annotated",
    buildver="hg38",
    protocol="refGene,dbnsfp42a",
    operation="g,f",
    otherinfo=True,
    use_mane_transcript=True,
    thread=8,
)
annotator.run_annotation()
# → sample_annotated.hg38_multianno.tsv

# Same VCF, but convert to .mvinput first (legacy path)
TableAnnotator(
    queryfile="sample.vcf",
    dbloc="/path/to/humandb",
    outfile="sample_via_mv",
    buildver="hg38",
    protocol="refGene",
    operation="g",
    convertvcf=True,
).run_annotation()
```

### Simulate variants, then annotate

```python
from matchvar_annotator import GeneTranscript, TableAnnotator, VariantMetricCalculator
import pandas as pd

tx = GeneTranscript.from_gtf(
    gene_name="CFTR",
    transcript_id="NM_000492.4",
    gtf_file="gencode.v44.annotation.gtf",
    fasta_file="GRCh38.fa",
)
variants = tx.generate_all_variants()
tx.export_to_vcf(variants, "cftr_variants.vcf")

from matchvar_annotator import EnhancedGeneTranscript
enh = EnhancedGeneTranscript.from_gtf(
    gene_name="CFTR",
    transcript_id="NM_000492.4",
    gtf_file="gencode.v44.annotation.gtf",
    fasta_file="GRCh38.fa",
    genome_version="hg38",
)
spectrum = enh.generate_spectrum_aware_variants(max_variants=200, use_constraint_filter=False)


TableAnnotator(
    queryfile="cftr_variants.vcf",
    dbloc="humandb",
    outfile="cftr_annotated",
    buildver="hg38",
    protocol="refGene,dbnsfp42a",
    operation="g,f",
    otherinfo=True,
).run_annotation()

calc = VariantMetricCalculator()
df = pd.read_csv("cftr_annotated.hg38_multianno.tsv", sep="\t")
y_true = calc.extract_labels_from_annotation(df)
```

Equivalent one-shot:

```python
from matchvar_annotator import MatchingPipeline, run_pipeline

pipeline = MatchingPipeline(
    gtf_file="annotation.gtf",
    fasta_file="genome.fa",
    gene_name="BRCA1",
    transcript_id="NM_007294.4",
    database_dir="humandb",
    output_dir="./results",
    protocols=["refGene", "exac03", "avsift"],
    operations=["g", "f", "f"],
    variant_types=["SNV", "insertion", "deletion"],
    buildver="hg38",
    threads=8,
    simulation_strategy="hybrid",
    spectrum_ratio=0.7,
    enable_visualization=True,
    keep_temp=False,
)
results = pipeline.run()
```

### Helpers

```python
from matchvar_annotator import format_g_hgvs, parse_g_hgvs, classify_row, annotate_g_hgvs

format_g_hgvs("17", 41276045, 41276045, "A", "G")
# 'chr17:g.41276045A>G'

effect, vartype = classify_row(
    function="UTR5",
    gene="BRCA1",
    c_hgvs="NM_007294:c.-12A>G",
    ref="A",
    alt="G",
)
# ('RV', 'RV')  — regulatory variation
```

## Protocols and operations

| Operation | Meaning | Typical protocols |
|-----------|---------|-------------------|
| `g` | Gene annotation | `refGene`, `ensGene`, `knownGene`, `ncbiRefSeq`, `gencode` |
| `r` | Region overlap | `cytoBand`, `dgvMerged`, `wgRna` |
| `f` | Point/filter DB | `exac03`, `gnomad211_genome`, `esp6500siv2_all`, `1000g2015aug_all`, `clinvar`, `cosmic102`, `avsift`, `dbnsfp42a`, `revel`, `cadd13gt10`, `AlphaMissense`, `dbscsnv11` |

`dbscsnv11` is registered but **not** in the default protocol list; pass it explicitly if the file exists in `humandb`.

Default `MatchvarRunner` protocols when `protocols=None`:

`refGene, ensGene, knownGene, cytoBand, exac03, avsift, dbnsfp42a, gnomad211_genome, esp6500siv2_all, revel, cadd13gt10, AlphaMissense`

## Environment variables

```bash
export MATCHVAR_RESOURCES_DIR=/path/to/resources   # humandb lives under this directory
export PYTHON_EXECUTABLE=/path/to/python           # interpreter for subprocesses
export MPLBACKEND=Agg                              # headless plotting
```

## Tests

```bash
pip install -e ".[dev]"
pytest tests/
```

## Troubleshooting

| Symptom | What to check |
|---------|----------------|
| File not found on `matchvar-annotator /abs/path.vcf` | Fixed in 1.2.0: the CLI uses the real path. Confirm `--resources-dir` / `humandb` exists. |
| `Unknown protocol: dbscsnv11` | It is no longer a default. Add it only if the DB file is present. |
| `Convert2Matchvar(..., format="vcf4")` ValueError | Do not set `snpqual`/`snppvalue` for vcf4; defaults are now `None`. |
| Simulated deletions rejected by annotators | Deletions are left-anchored VCF; regenerate the VCF with the current package. |
| Empty ClinVar ROC | Matching uses `Ref`/`Alt` (and VCF allele aliases) plus chromosome without requiring a `chr` prefix. |
| spectrum/hybrid still looks traditional | Check logs: empty spectrum output falls back. Ensure GTF+FASTA yield a CDS, and try a lower `--min-probability`. |
| `No MANE transcript file` | Place a MANE GTF under `resources/mane/`, `humandb/mane_transcript.txt`, or `humandb/mane.bb` (needs pyBigWig). |

## Project layout

```
matchvar-annotator/
├── matchvar_annotator/
│   ├── matchvar_annotator.py   # MatchvarRunner
│   ├── table_matchvar.py       # multi-protocol table annotation
│   ├── annotate_variation.py   # gene / region / filter engine
│   ├── convert2matchvar.py
│   ├── coding_change.py        # c./p. HGVS + 3′ indel rule
│   ├── column_names.py / variant_typing.py / mane_transcripts.py
│   ├── ghgvs.py / resource_files.py / twobit_reader.py
│   ├── variant_simulation.py / pipeline.py
│   ├── mutation_spectrum_config.py / mutation_spectrum_generator.py / enhanced_data_simulation.py
│   └── cli.py, db_cli.py, ...
├── tests/
├── examples.py
└── setup.py
```

```bash
python examples.py all    # prints usage patterns (edit paths first)
python examples.py 3      # step-by-step simulation + annotation sketch
```

## License

MIT. Author: Bingbo Zhou \<zhoubingbo@hotmail.com\>

- Source: [https://github.com/zhoubingbo/matchvar-annotator](https://github.com/zhoubingbo/matchvar-annotator)
- Web tool: [https://matchvar.intelligene.cn/](https://matchvar.intelligene.cn/)
