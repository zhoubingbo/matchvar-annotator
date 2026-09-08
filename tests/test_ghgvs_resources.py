#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Tests for g.HGVS conversion and humandb resource discovery."""

import os
from pathlib import Path

import pytest

from matchvar_annotator.column_names import (
    fill_g_hgvs_alleles,
    format_g_hgvs,
    parse_g_hgvs,
    to_mvinput_line,
)
from matchvar_annotator.ghgvs import ghgvs_to_mvinput_lines
from matchvar_annotator.resource_files import (
    can_use_gene_bigbed,
    chrom_aliases,
    discover_gene_protocols,
    ensure_gene_pred,
    parse_bigbed_gene_entry,
    pybigwig_available,
    resolve_cytoband,
    resolve_gene_bigbed,
    resolve_genome_2bit,
    resolve_genome_fasta,
)


class TestParseGHgvs:
    def test_snv(self):
        p = parse_g_hgvs("chr17:g.41276045A>G")
        assert p["chrom"] == "chr17"
        assert p["start"] == 41276045
        assert p["ref"] == "A"
        assert p["alt"] == "G"
        assert to_mvinput_line(p) == "17\t41276045\t41276045\tA\tG"

    def test_snv_no_chr_prefix(self):
        p = parse_g_hgvs("5:g.100A>T")
        assert p["chrom"] == "chr5"
        assert to_mvinput_line(p) == "5\t100\t100\tA\tT"

    def test_deletion_with_seq(self):
        p = parse_g_hgvs("chr5:g.131008083_131008084delCC")
        assert p["start"] == 131008083
        assert p["end"] == 131008084
        assert p["ref"] == "CC"
        assert p["alt"] == "-"
        assert to_mvinput_line(p) == "5\t131008083\t131008084\tCC\t-"

    def test_insertion(self):
        p = parse_g_hgvs("chr5:g.131008083_131008084insTGACAGTTGTTTGCACAG")
        assert p["start"] == 131008083
        assert p["end"] == 131008083
        assert p["ref"] == "-"
        assert p["alt"] == "TGACAGTTGTTTGCACAG"
        assert "insTGACAGTTGTTTGCACAG" in to_mvinput_line(p) or True
        assert to_mvinput_line(p) == "5\t131008083\t131008083\t-\tTGACAGTTGTTTGCACAG"

    def test_delins(self):
        p = parse_g_hgvs("chr1:g.100_102delinsAA")
        assert p["start"] == 100
        assert p["end"] == 102
        assert p["alt"] == "AA"

    def test_nc_accession(self):
        p = parse_g_hgvs("NC_000005.9:g.100A>G")
        assert p["chrom"] == "chr5"

    def test_dup_with_seq(self):
        p = parse_g_hgvs("chr1:g.10dupA")
        assert p["ref"] == "-"
        assert p["alt"] == "A"
        assert p["start"] == 10

    def test_dup_without_seq_needs_fetch(self):
        p = parse_g_hgvs("chr1:g.10_12dup")
        assert p.get("_dup") is True
        filled = fill_g_hgvs_alleles(p, lambda chrom, s, e: "ATG")
        assert filled["ref"] == "-"
        assert filled["alt"] == "ATG"
        assert filled["start"] == 12

    def test_deletion_without_seq_needs_fetch(self):
        p = parse_g_hgvs("chr1:g.10_11del")
        filled = fill_g_hgvs_alleles(p, lambda chrom, s, e: "CG")
        assert filled["ref"] == "CG"
        assert filled["alt"] == "-"

    def test_roundtrip_snv(self):
        hgvs = format_g_hgvs("1", 100, 100, "A", "G")
        p = parse_g_hgvs(hgvs)
        assert to_mvinput_line(p) == "1\t100\t100\tA\tG"

    def test_invalid(self):
        with pytest.raises(ValueError):
            parse_g_hgvs("c.123A>G")


class TestGhgvsToMvinput:
    def test_list(self):
        lines = ghgvs_to_mvinput_lines(
            ["chr1:g.100A>G", "chr5:g.10_11insAT"]
        )
        assert lines[0] == "1\t100\t100\tA\tG"
        assert lines[1] == "5\t10\t10\t-\tAT"


class TestResourceDiscovery:
    def test_chrom_aliases(self):
        assert "chr5" in chrom_aliases("5")
        assert "5" in chrom_aliases("chr5")

    def test_discover_in_tmp(self, tmp_path):
        (tmp_path / "hg19_refGene.txt").write_text("x\n")
        (tmp_path / "cytoBand.txt.gz").write_bytes(b"\x1f\x8b")
        (tmp_path / "ncbiRefSeq.bb").write_bytes(b"not-a-real-bb")
        (tmp_path / "gencodeV49lift37.bb").write_bytes(b"not-a-real-bb")
        (tmp_path / "hg19.2bit").write_bytes(b"xxxx")
        (tmp_path / "hg19.chrom.sizes").write_text("chr1\t249250621\n")
        found = discover_gene_protocols(str(tmp_path), "hg19")
        names = [p for p, _ in found]
        assert "refGene" in names
        assert "cytoBand" in names
        assert resolve_cytoband(str(tmp_path), "hg19")
        assert resolve_gene_bigbed(str(tmp_path), "ncbiRefSeq", "hg19")
        assert resolve_gene_bigbed(str(tmp_path), "gencode", "hg19")
        assert resolve_genome_2bit(str(tmp_path), "hg19")
        # Fake .bb files are not usable without pyBigWig or a genePred cache.
        if pybigwig_available():
            assert "ncbiRefSeq" in names
            assert "gencode" in names
        else:
            assert "ncbiRefSeq" not in names
            assert "gencode" not in names
            assert not can_use_gene_bigbed(str(tmp_path / "ncbiRefSeq.bb"))

    def test_discover_uses_cached_genepred(self, tmp_path):
        (tmp_path / "hg19_refGene.txt").write_text("x\n")
        bb = tmp_path / "ncbiRefSeq.bb"
        bb.write_bytes(b"not-a-real-bb")
        (tmp_path / "ncbiRefSeq.bb.genePred.txt").write_text(
            "0\tNM_1\tchr1\t+\t0\t10\t0\t10\t1\t0,\t10,\t0\tGENE\tcmpl\tcmpl\n"
        )
        names = [p for p, _ in discover_gene_protocols(str(tmp_path), "hg19")]
        assert "refGene" in names
        assert "ncbiRefSeq" in names
        assert can_use_gene_bigbed(str(bb))

    def test_ensure_gene_pred_unreadable_bb_raises(self, tmp_path):
        (tmp_path / "ncbiRefSeq.bb").write_bytes(b"not-a-real-bb")
        with pytest.raises(RuntimeError, match="Cannot load gene models"):
            ensure_gene_pred(str(tmp_path), "ncbiRefSeq", "hg19")

    def test_parse_bed12_entry(self):
        rest = "NM_1\t0\t-\t10\t20\t0\t2\t5,5,\t0,10,\tGENE"
        rec = parse_bigbed_gene_entry("chr1", 100, 120, rest)
        assert rec["name"] == "NM_1"
        assert rec["strand"] == "-"
        assert rec["name2"] == "GENE"
        assert rec["exonStarts"] == [100, 110]
        assert rec["exonEnds"] == [105, 115]


@pytest.mark.skipif(
    not os.path.isfile(
        os.path.join(
            os.path.dirname(__file__),
            "..",
            "resources",
            "humandb",
            "hg19_refGene.txt",
        )
    ),
    reason="local humandb not installed",
)
class TestAnnotateGHgvsLive:
    def test_fnip1_insertion(self, tmp_path):
        from matchvar_annotator.ghgvs import annotate_g_hgvs

        dbloc = os.path.abspath(
            os.path.join(os.path.dirname(__file__), "..", "resources", "humandb")
        )
        df = annotate_g_hgvs(
            ["chr5:g.131008083_131008084insTGACAGTTGTTTGCACAG"],
            dbloc=dbloc,
            buildver="hg19",
            outfile=str(tmp_path / "annot"),
            protocol="refGene",
            operation="g",
        )
        assert len(df) == 1
        ccol = "cHGVS.refGene"
        assert ccol in df.columns
        text = str(df.iloc[0][ccol])
        assert "NM_133372" in text
        assert "insCTGTGCAAACAACTGTCA" in text
