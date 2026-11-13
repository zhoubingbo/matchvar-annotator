# Changelog

All notable changes to `matchvar-annotator` are documented here.
The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.2.0]

### Added

- VCF (classic 4-column and full VCF) and MATCHVAR `.mvinput` annotation against a local `humandb`.
- Pileup / GFF3 / VCF → `.mvinput` conversion via `matchvar-convert`.
- c.HGVS / p.HGVS polishing, including HGVS 3′ indel normalization.
- MANE Select filtering and gene/region/filter annotation for `refGene`, `ensGene`, `knownGene`, `ncbiRefSeq`, `gencode`, and `cytoBand`.
- Optional simulate → annotate → AUROC pipeline (`matchvar-pipeline`).
- Console entry points: `matchvar-annotator`, `matchvar-table`, `matchvar-convert`, `matchvar-coding`, `matchvar-index`, `matchvar-db`, `matchvar-clinvardb`, `matchvar-dbnsfpdb`, `matchvar-pipeline`, `matchvar-ghgvs`.
- Optional `pyBigWig` extra (`pip install "matchvar-annotator[bb]"`) for reading `*.bb` resources.
- Packaging metadata for PyPI: `pyproject.toml`, `MANIFEST.in`, license metadata, and CI publishing workflow.

### Fixed

- CLI resolves absolute VCF paths correctly.
- `dbscsnv11` is no longer a default protocol.
- `Convert2Matchvar(format="vcf4")` no longer requires `snpqual` / `snppvalue`.
- Simulated deletions are left-anchored VCF alleles.
- ClinVar ROC joins on `Ref` / `Alt` plus chromosome without requiring a `chr` prefix.
- Empty spectrum output falls back to the traditional spectrum.
- `*.bb` protocols are skipped when `pyBigWig` is unavailable instead of emitting empty columns.

[1.2.0]: https://github.com/zhoubingbo/matchvar-annotator/releases/tag/v1.2.0
