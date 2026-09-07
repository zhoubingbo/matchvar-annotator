#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Near-gene c.-N/c.*N HGVS, mutation-spectrum simulation, and annotation golden tests."""

import os
from pathlib import Path

import pytest

from matchvar_annotator.annotate_variation import AnnotateVariation
from matchvar_annotator.enhanced_data_simulation import EnhancedGeneTranscript
from matchvar_annotator.table_matchvar import TableAnnotator


def _annotator(tmp_path, neargene=1000, use_mane=False):
    query = tmp_path / "q.mvinput"
    query.write_text("1\t1\t1\tA\tG\n")
    dbloc = tmp_path / "humandb"
    dbloc.mkdir()
    return AnnotateVariation(
        str(query),
        str(dbloc),
        geneanno=True,
        dbtype="refGene",
        neargene=neargene,
        use_mane_transcript=use_mane,
    )


def _plus_strand_gene(**overrides):
    gene = {
        "name": "NM_TEST.1",
        "chrom": "1",
        "strand": "+",
        "txStart": 10001,
        "txEnd": 11000,
        "cdsStart": 10101,
        "cdsEnd": 10900,
        "name2": "TESTGENE",
    }
    gene.update(overrides)
    return gene


def _tiny_transcript(strand="+"):
    cds = "ATGAAACCCGGGTTTAAA"  # 18 bp
    return EnhancedGeneTranscript(
        gene_name="TESTGENE",
        transcript_id="NM_TEST.1",
        exons=[{
            "genomic_start": 1000,
            "genomic_end": 1018,
            "cds_sequence": cds,
        }],
        chromosome="chr1",
        strand=strand,
        genome_version="hg19",
    )


def _write_refgene(path: Path, buildver="hg19"):
    # UCSC 0-based genePred; loader converts txStart/cdsStart to 1-based.
    line = "\t".join([
        "0",
        "NM_TEST.1",
        "chr1",
        "+",
        "10000",
        "11000",
        "10100",
        "10900",
        "1",
        "10000,",
        "11000,",
        "0",
        "TESTGENE",
        "cmpl",
        "cmpl",
    ]) + "\n"
    (path / f"{buildver}_refGene.txt").write_text(line)


def _find_humandb():
    candidates = []
    env = os.environ.get("MATCHVAR_RESOURCES_DIR")
    if env:
        candidates.append(os.path.join(env, "humandb"))
        candidates.append(env)
    here = Path(__file__).resolve().parent
    pkg_root = here.parent
    candidates.extend([
        os.path.join(os.getcwd(), "resources", "humandb"),
        str(pkg_root / "resources" / "humandb"),
        str(pkg_root.parent / "gene_annotation" / "resources" / "humandb"),
        os.path.expanduser("~/humandb"),
    ])
    for dbloc in candidates:
        if not dbloc or not os.path.isdir(dbloc):
            continue
        for buildver in ("hg19", "hg38"):
            for name in (f"{buildver}_refGene.txt", "refGene.txt"):
                if os.path.isfile(os.path.join(dbloc, name)):
                    return dbloc, buildver
    return None, None


class TestNearGeneHgvs:
    def test_upstream_plus_strand_c_minus(self, tmp_path):
        ann = _annotator(tmp_path)
        gene_db = {"1": [_plus_strand_gene()]}
        variant = {
            "chrom": "1",
            "start": 9950,
            "end": 9950,
            "ref": "A",
            "alt": "G",
        }
        out = ann._annotate_variant_by_gene(variant, gene_db)
        assert out["function"] == "utr5"
        assert out["gene"] == "TESTGENE"
        assert "c.-151A>G" in out["gene_detail"]
        assert out["gene_detail"].startswith("NM_TEST.1:UTR5:")

    def test_downstream_plus_strand_c_star(self, tmp_path):
        ann = _annotator(tmp_path)
        gene_db = {"1": [_plus_strand_gene()]}
        variant = {
            "chrom": "1",
            "start": 11050,
            "end": 11050,
            "ref": "C",
            "alt": "T",
        }
        out = ann._annotate_variant_by_gene(variant, gene_db)
        assert out["function"] == "utr3"
        assert "c.*150C>T" in out["gene_detail"]
        assert ":UTR3:" in out["gene_detail"]

    def test_far_away_stays_intergenic(self, tmp_path):
        ann = _annotator(tmp_path, neargene=1000)
        gene_db = {"1": [_plus_strand_gene()]}
        variant = {
            "chrom": "1",
            "start": 20000,
            "end": 20000,
            "ref": "A",
            "alt": "G",
        }
        out = ann._annotate_variant_by_gene(variant, gene_db)
        assert out["function"] == "intergenic"
        assert out["gene_detail"] == "NA"

    def test_minus_strand_txstart_is_utr3(self, tmp_path):
        ann = _annotator(tmp_path)
        gene = _plus_strand_gene(strand="-", name="NM_MINUS.1")
        gene_db = {"1": [gene]}
        variant = {
            "chrom": "1",
            "start": 9950,
            "end": 9950,
            "ref": "A",
            "alt": "G",
        }
        out = ann._annotate_variant_by_gene(variant, gene_db)
        assert out["function"] == "utr3"
        assert "c.*" in out["gene_detail"]
        assert "T>C" in out["gene_detail"]  # reverse complement of A>G

    def test_near_gene_deletion_allele(self, tmp_path):
        ann = _annotator(tmp_path)
        out = ann._format_c_allele_change("AT", "-", strand="+")
        assert out == "delAT"


class TestSyntheticGenePredGolden:
    def test_annotate_variation_writes_near_gene_hgvs(self, tmp_path):
        dbloc = tmp_path / "humandb"
        dbloc.mkdir()
        _write_refgene(dbloc)
        query = tmp_path / "query.mvinput"
        query.write_text(
            "1\t9950\t9950\tA\tG\n"
            "1\t10500\t10500\tA\tG\n"
            "1\t11050\t11050\tA\tG\n"
            "1\t20000\t20000\tA\tG\n"
        )
        outfile = tmp_path / "out"
        ann = AnnotateVariation(
            str(query),
            str(dbloc),
            outfile=str(outfile),
            geneanno=True,
            dbtype="refGene",
            buildver="hg19",
            neargene=1000,
        )
        ann.run_annotation()
        vf = (tmp_path / "out.variant_function").read_text().splitlines()
        assert len(vf) == 4
        assert vf[0].startswith("utr5\t")
        assert "c.-151A>G" in vf[0]
        assert vf[1].startswith("exonic\t") or vf[1].startswith("UTR")
        assert vf[2].startswith("utr3\t")
        assert "c.*150A>G" in vf[2]
        assert vf[3].startswith("intergenic\t")

    def test_table_annotator_gene_columns(self, tmp_path):
        dbloc = tmp_path / "humandb"
        dbloc.mkdir()
        _write_refgene(dbloc)
        query = tmp_path / "query.mvinput"
        query.write_text("1\t9950\t9950\tA\tG\n1\t11050\t11050\tC\tT\n")
        outfile = tmp_path / "tbl"
        TableAnnotator(
            str(query),
            str(dbloc),
            outfile=str(outfile),
            buildver="hg19",
            protocol="refGene",
            operation="g",
            remove=True,
        ).run_annotation()
        tsv = tmp_path / "tbl.hg19_multianno.tsv"
        assert tsv.exists()
        header, *rows = tsv.read_text().splitlines()
        cols = header.split("\t")
        assert "Function.refGene" in cols
        assert "cHGVS.refGene" in cols
        func_i = cols.index("Function.refGene")
        chgvs_i = cols.index("cHGVS.refGene")
        funcs = [r.split("\t")[func_i] for r in rows]
        chgvs = [r.split("\t")[chgvs_i] for r in rows]
        assert "UTR5" in funcs[0] or "utr5" in funcs[0]
        assert "c.-" in chgvs[0]
        assert "UTR3" in funcs[1] or "utr3" in funcs[1]
        assert "c.*" in chgvs[1]


class TestSpectrumSimulation:
    def test_spectrum_generates_snvs_without_fasta(self):
        tx = _tiny_transcript()
        result = tx.generate_spectrum_aware_variants(
            max_variants=20,
            min_probability=0.001,
            use_constraint_filter=False,
        )
        assert result["total"] > 0
        assert result["spectrum_aware"] is True
        snv = result["variants"][0]
        assert snv["type"] == "SNV"
        assert snv["hgvs"].startswith("c.")
        assert snv["ref_seq"] in "ACGT"
        assert snv["alt_seq"] in "ACGT"
        assert snv["mutation_probability"] >= 0.001
        assert len(snv["trinucleotide_context"]) == 3

    def test_hybrid_includes_spectrum_and_traditional(self):
        tx = _tiny_transcript()
        result = tx.generate_hybrid_variants(
            spectrum_ratio=0.5,
            max_variants=10,
            min_probability=0.001,
            use_constraint_filter=False,
        )
        assert result["hybrid"] is True
        assert result["total"] > 0
        assert any(v.get("type") == "SNV" for v in result["all_variants"])


def _find_unoccupied_near_gene_site(gene_db, neargene=1000):
    """Pick a 1-based position within neargene of some transcript and overlapping none."""
    for chrom, genes in gene_db.items():
        occupied = [(int(g["txStart"]), int(g["txEnd"])) for g in genes]
        for gene in genes:
            tx_start = int(gene["txStart"])
            tx_end = int(gene["txEnd"])
            for pos in (tx_start - 50, tx_end + 50):
                if pos < 1:
                    continue
                if any(start <= pos <= end for start, end in occupied):
                    continue
                dist = (tx_start - pos) if pos < tx_start else (pos - tx_end)
                if 0 < dist <= neargene:
                    return chrom, pos
    return None, None


class TestHumandbGolden:
    def test_real_humandb_near_gene_if_present(self, tmp_path):
        dbloc, buildver = _find_humandb()
        if not dbloc:
            pytest.skip("humandb refGene not found")

        probe = _annotator(tmp_path)
        probe.dbloc = dbloc
        probe.buildver = buildver
        probe.dbtype1 = "refGene"
        gene_db = probe._load_gene_database()
        chrom, pos = _find_unoccupied_near_gene_site(gene_db)
        if not chrom:
            pytest.skip("no unoccupied near-gene site in this refGene")
        query = tmp_path / "real.mvinput"
        query.write_text(f"{chrom}\t{pos}\t{pos}\tA\tG\n")
        outfile = tmp_path / "real_out"
        AnnotateVariation(
            str(query),
            dbloc,
            outfile=str(outfile),
            geneanno=True,
            dbtype="refGene",
            buildver=buildver,
            neargene=1000,
        ).run_annotation()
        vf = (tmp_path / "real_out.variant_function").read_text()
        assert vf.strip()
        assert ("c.-" in vf) or ("c.*" in vf) or vf.startswith("utr5") or vf.startswith("utr3")
        assert not vf.startswith("intergenic")

    def test_ndufa10_intron_del_matches_vep_span(self, tmp_path):
        """VEP: NM_004544.4 c.1000-12385_1000-12365del — both intron ends."""
        dbloc, buildver = _find_humandb()
        if not dbloc or buildver != "hg19":
            pytest.skip("hg19 refGene not found")
        query = tmp_path / "ndufa10.mvinput"
        query.write_text("2\t240912968\t240912988\tACGATCAGCTTCAAGGAGGGA\t-\n")
        outfile = tmp_path / "ndufa10_out"
        AnnotateVariation(
            str(query),
            dbloc,
            outfile=str(outfile),
            geneanno=True,
            dbtype="refGene",
            buildver="hg19",
        ).run_annotation()
        vf = (tmp_path / "ndufa10_out.variant_function").read_text()
        assert "NM_004544.4" in vf
        assert "c.1000-12385_1000-12365del" in vf
        assert "NDUFA10" in vf
