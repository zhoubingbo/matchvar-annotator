#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""VCF CHROM/POS/REF/ALT can be annotated without convert2matchvar."""

from pathlib import Path

import pytest

from matchvar_annotator.query_format import (
    iter_query_records,
    normalize_vcf_alleles,
    parse_query_line,
    sniff_query_format,
)


class TestNormalizeVcfAlleles:
    def test_snv(self):
        assert normalize_vcf_alleles(100, "A", "G") == (100, 100, "A", "G")

    def test_insertion_left_anchor(self):
        # VCF POS A AAT → insertion of AT after POS
        assert normalize_vcf_alleles(103234292, "A", "AAT") == (
            103234292,
            103234292,
            "-",
            "AT",
        )

    def test_deletion_left_anchor(self):
        # VCF POS ATG A → delete TG spanning POS+1..POS+2
        assert normalize_vcf_alleles(103248997, "ATG", "A") == (
            103248998,
            103248999,
            "TG",
            "-",
        )

    def test_mnp(self):
        assert normalize_vcf_alleles(10, "AT", "GC") == (10, 11, "AT", "GC")


class TestParseQueryLine:
    def test_four_column_vcf_matches_mvinput(self):
        vcf = parse_query_line("chr12\t103234292\tA\tAAT", fmt="vcf")
        mv = parse_query_line("chr12\t103234292\t103234292\t-\tAT", fmt="mvinput")
        assert len(vcf) == 1 and len(mv) == 1
        for key in ("chrom", "start", "end", "ref", "alt"):
            assert vcf[0][key] == mv[0][key]
        assert vcf[0]["query_line"] == mv[0]["query_line"]

    def test_deletion_four_column_equals_mvinput(self):
        vcf = parse_query_line("chr12\t103248997\tATG\tA", fmt="vcf")
        mv = parse_query_line("chr12\t103248998\t103248999\tTG\t-", fmt="mvinput")
        assert vcf[0]["query_line"] == mv[0]["query_line"]

    def test_full_vcf_skips_id(self):
        recs = parse_query_line(
            "chr1\t100\t.\tA\tG\t.\tPASS\tDP=10", fmt="vcf"
        )
        assert recs[0]["start"] == 100
        assert recs[0]["ref"] == "A"
        assert recs[0]["alt"] == "G"
        assert recs[0]["extra_fields"][0] == "."

    def test_multiallelic_splits(self):
        recs = parse_query_line("chr1\t100\tA\tG,T", fmt="vcf")
        assert [r["alt"] for r in recs] == ["G", "T"]

    def test_symbolic_alt_skipped(self):
        assert parse_query_line("chr1\t100\tA\t<DEL>", fmt="vcf") == []


class TestSniffQueryFormat:
    def test_vcf_extension(self, tmp_path):
        p = tmp_path / "x.vcf"
        p.write_text("chr1\t100\tA\tG\n")
        assert sniff_query_format(str(p)) == "vcf"

    def test_four_column_without_extension(self, tmp_path):
        p = tmp_path / "sites.tsv"
        p.write_text("chr1\t100\tA\tG\n")
        assert sniff_query_format(str(p)) == "vcf"

    def test_mvinput_five_column(self, tmp_path):
        p = tmp_path / "sites.txt"
        p.write_text("1\t100\t100\tA\tG\n")
        assert sniff_query_format(str(p)) == "mvinput"

    def test_vcf_header(self, tmp_path):
        p = tmp_path / "x.txt"
        p.write_text(
            "##fileformat=VCFv4.2\n"
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
            "chr1\t100\t.\tA\tG\t.\tPASS\t.\n"
        )
        assert sniff_query_format(str(p)) == "vcf"


def _humandb():
    root = Path(__file__).resolve().parents[1]
    for dbloc in (
        root / "resources" / "humandb",
        root / "resources" / "humandb" / "hg19",
    ):
        if (dbloc / "hg19_refGene.txt").is_file() or (dbloc / "refGene.txt").is_file():
            return dbloc
    return None


@pytest.mark.skipif(_humandb() is None, reason="resources/humandb with refGene not present")
def test_vcf_and_mvinput_gene_annotation_agree(tmp_path):
    from matchvar_annotator.table_matchvar import TableAnnotator

    dbloc = str(_humandb())
    vcf = tmp_path / "q.vcf"
    mv = tmp_path / "q.mvinput"
    # 4-column VCF vs equivalent MATCHVAR alleles
    vcf.write_text("chr12\t103234292\tA\tAAT\nchr12\t103248997\tATG\tA\nchr1\t100\tA\tG\n")
    mv.write_text(
        "chr12\t103234292\t103234292\t-\tAT\n"
        "chr12\t103248998\t103248999\tTG\t-\n"
        "chr1\t100\t100\tA\tG\n"
    )
    out_vcf = tmp_path / "ann_vcf"
    out_mv = tmp_path / "ann_mv"
    TableAnnotator(
        queryfile=str(vcf),
        dbloc=dbloc,
        outfile=str(out_vcf),
        buildver="hg19",
        protocol="refGene",
        operation="g",
        nopolish=True,
        remove=True,
    ).run_annotation()
    TableAnnotator(
        queryfile=str(mv),
        dbloc=dbloc,
        outfile=str(out_mv),
        buildver="hg19",
        protocol="refGene",
        operation="g",
        nopolish=True,
        remove=True,
    ).run_annotation()
    tsv_vcf = Path(f"{out_vcf}.hg19_multianno.tsv")
    tsv_mv = Path(f"{out_mv}.hg19_multianno.tsv")
    assert tsv_vcf.is_file() and tsv_mv.is_file()
    lines_v = tsv_vcf.read_text().splitlines()
    lines_m = tsv_mv.read_text().splitlines()
    assert lines_v[0] == lines_m[0]
    assert len(lines_v) == len(lines_m)
    # Chr Start End Ref Alt and gene function columns should match
    for a, b in zip(lines_v[1:], lines_m[1:]):
        pa, pb = a.split("\t"), b.split("\t")
        assert pa[:5] == pb[:5]
        assert pa[6:8] == pb[6:8]


def test_iter_query_records_skips_header(tmp_path):
    p = tmp_path / "x.vcf"
    p.write_text(
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\n"
        "1\t100\t.\tA\tG\n"
    )
    recs = list(iter_query_records(str(p)))
    assert len(recs) == 1
    assert recs[0]["query_line"] == "1\t100\t100\tA\tG"


_FULL_VCF = (
    "##fileformat=VCFv4.2\n"
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
    "chr1\t100\t.\tA\tG\t.\tPASS\tDP=9\n"
    "chr12\t103234292\t.\tA\tAAT\t.\tPASS\t.\n"
    "chr12\t103248997\t.\tATG\tA\t.\tPASS\t.\n"
)


def test_convertvcf_writes_mvinput(tmp_path):
    from matchvar_annotator.table_matchvar import TableAnnotator

    vcf = tmp_path / "x.vcf"
    vcf.write_text(_FULL_VCF)
    ann = TableAnnotator(
        queryfile=str(vcf),
        dbloc=str(tmp_path),
        outfile=str(tmp_path / "out"),
        protocol="refGene",
        operation="g",
        convertvcf=True,
        nopolish=True,
    )
    ann._convert_vcf_to_mvinput()
    assert ann.query_format == "mvinput"
    mv = Path(ann.queryfile)
    assert mv.name.endswith(".mvinput")
    rows = [ln for ln in mv.read_text().splitlines() if ln.strip()]
    assert len(rows) == 3
    first = rows[0].split("\t")
    assert first[:5] == ["chr1", "100", "100", "A", "G"]
    ins = rows[1].split("\t")
    assert ins[:5] == ["chr12", "103234292", "103234292", "-", "AT"]
    dele = rows[2].split("\t")
    assert dele[:5] == ["chr12", "103248998", "103248999", "TG", "-"]


def test_convertvcf_skipped_when_already_mvinput(tmp_path):
    from matchvar_annotator.table_matchvar import TableAnnotator

    mv = tmp_path / "q.mvinput"
    mv.write_text("1\t100\t100\tA\tG\n")
    ann = TableAnnotator(
        queryfile=str(mv),
        dbloc=str(tmp_path),
        outfile=str(tmp_path / "out"),
        protocol="refGene",
        operation="g",
        convertvcf=True,
        nopolish=True,
    )
    ann._convert_vcf_to_mvinput()
    assert Path(ann.queryfile) == mv
    assert ann.query_format == "mvinput"


@pytest.mark.skipif(_humandb() is None, reason="resources/humandb with refGene not present")
def test_native_vcf_and_convertvcf_annotation_agree(tmp_path):
    from matchvar_annotator.table_matchvar import TableAnnotator

    dbloc = str(_humandb())
    vcf = tmp_path / "q.vcf"
    vcf.write_text(_FULL_VCF)
    out_native = tmp_path / "ann_native"
    out_conv = tmp_path / "ann_conv"
    TableAnnotator(
        queryfile=str(vcf),
        dbloc=dbloc,
        outfile=str(out_native),
        buildver="hg19",
        protocol="refGene",
        operation="g",
        nopolish=True,
        remove=True,
    ).run_annotation()
    TableAnnotator(
        queryfile=str(vcf),
        dbloc=dbloc,
        outfile=str(out_conv),
        buildver="hg19",
        protocol="refGene",
        operation="g",
        convertvcf=True,
        nopolish=True,
        remove=True,
    ).run_annotation()
    tsv_n = Path(f"{out_native}.hg19_multianno.tsv").read_text().splitlines()
    tsv_c = Path(f"{out_conv}.hg19_multianno.tsv").read_text().splitlines()
    assert len(tsv_n) == len(tsv_c)
    for a, b in zip(tsv_n[1:], tsv_c[1:]):
        pa, pb = a.split("\t"), b.split("\t")
        assert pa[:5] == pb[:5]
        assert pa[6:8] == pb[6:8]

