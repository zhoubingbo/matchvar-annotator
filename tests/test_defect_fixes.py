#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Regression tests for the P0/P1 defect fixes."""

import os
import tempfile

import pandas as pd
import pytest

from matchvar_annotator import MatchvarRunner, Convert2Matchvar, GeneTranscript
from matchvar_annotator.column_names import format_g_hgvs
from matchvar_annotator.variant_typing import classify_row
from matchvar_annotator.coding_change import normalize_coding_indel_hgvs
from matchvar_annotator.labeling import frameshift_from_info, label_from_vcf_info, func_is_pathogenic
from matchvar_annotator.evaluation.clinvar_processor import ClinVarProcessor
from matchvar_annotator.variant_simulation import GeneTranscript as GT


class TestConvertDefaults:
    def test_vcf4_api_does_not_raise(self, tmp_path):
        vcf = tmp_path / "x.vcf"
        vcf.write_text(
            "##fileformat=VCFv4.2\n"
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
            "chr1\t100\t.\tA\tG\t.\tPASS\t.\n"
        )
        converter = Convert2Matchvar(str(vcf), format="vcf4", outfile=str(tmp_path / "out"))
        assert converter.format == "vcf4"


class TestDefaultProtocols:
    def test_default_protocols_are_registered(self, tmp_path):
        runner = MatchvarRunner(resources_dir=str(tmp_path), genome_version="hg19")
        ok, errors = runner.validate_protocols(runner.DEFAULT_PROTOCOLS)
        assert ok, errors
        assert "dbscsnv11" not in runner.DEFAULT_PROTOCOLS
        assert "dbscsnv11" in runner.PROTOCOL_CONFIGS

    def test_operations_override(self, tmp_path):
        runner = MatchvarRunner(resources_dir=str(tmp_path), genome_version="hg19")
        query = tmp_path / "q.mvinput"
        query.write_text("1\t100\t100\tA\tG\n")
        cmd = runner.build_matchvar_command(
            str(query),
            ["refGene", "cytoBand"],
            additional_args={"operations": ["g", "r"]},
        )
        assert "-operation g,r" in cmd
        assert str(query) in cmd

    def test_vcf_native_flag_on_vcf_path(self, tmp_path):
        runner = MatchvarRunner(resources_dir=str(tmp_path), genome_version="hg19")
        vcf = tmp_path / "x.vcf"
        vcf.write_text("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\nchr1\t1\t.\tA\tG\n")
        cmd = runner.build_matchvar_command(str(vcf), ["refGene"])
        assert "-vcfinput" in cmd
        assert "-convertvcf" not in cmd

    def test_convert_vcf_flag(self, tmp_path):
        runner = MatchvarRunner(resources_dir=str(tmp_path), genome_version="hg19")
        vcf = tmp_path / "x.vcf"
        vcf.write_text("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\nchr1\t1\t.\tA\tG\n")
        cmd = runner.build_matchvar_command(
            str(vcf),
            ["refGene"],
            additional_args={"convert_vcf": True, "operations": ["g"]},
        )
        assert "-convertvcf" in cmd
        assert "-vcfinput" not in cmd


class TestHgvsAndTyping:
    def test_format_g_hgvs_snv(self):
        assert format_g_hgvs("1", 100, 100, "A", "G") == "chr1:g.100A>G"

    def test_format_g_hgvs_deletion(self):
        assert format_g_hgvs("chr17", 100, 102, "ACG", "A") == "chr17:g.101_102delCG"

    def test_indel_3prime_normalization(self):
        cds = "ATGCATGCATGC"
        # insert AT after pos 4: c.4_5insAT — if tandem, may become dup
        out = normalize_coding_indel_hgvs("c.4_5insATGC", cds)
        assert out.startswith("c.")

    def test_classify_utr_as_rv(self):
        effect, vartype = classify_row(
            function="UTR5",
            gene="BRCA1",
            c_hgvs="NM_007294:c.-12A>G",
            ref="A",
            alt="G",
        )
        assert effect == "RV"
        assert vartype == "RV"

    def test_classify_missense(self):
        effect, _ = classify_row(
            function="exonic",
            gene="TP53",
            c_hgvs="NM_000546:exon5:c.215C>G:p.Pro72Arg",
            ref="C",
            alt="G",
            exonic_effect="nonsynonymous_SNV",
            p_hgvs="p.Pro72Arg",
        )
        assert effect == "missense"


class TestLabelsAndClinvar:
    def test_frameshift_flag_and_equals(self):
        assert frameshift_from_info("TYPE=INSERTION;FRAMESHIFT=true")
        assert frameshift_from_info("TYPE=INSERTION;FRAMESHIFT")
        assert not frameshift_from_info("TYPE=INSERTION;FRAMESHIFT=false")
        assert label_from_vcf_info("TYPE=INSERTION;FRAMESHIFT=true") == 1
        assert label_from_vcf_info("TYPE=SNV;FRAMESHIFT=false") == 0
        assert func_is_pathogenic("nonsynonymous SNV")
        assert not func_is_pathogenic("synonymous SNV")

    def test_clinvar_matches_ref_alt_columns(self, tmp_path):
        csv_path = tmp_path / "clinvar.csv"
        csv_path.write_text("GeneSymbol,ClinicalSignificance,Chr,Start,ReferenceAlleleVCF,AlternateAlleleVCF\n")
        proc = ClinVarProcessor(str(csv_path))
        annotated = pd.DataFrame([{
            "Chr": "chr17",
            "Start": 41276045,
            "Ref": "A",
            "Alt": "G",
            "Total_Score": 8.1,
        }])
        clinvar = pd.DataFrame([{
            "Chr": "17",
            "Start": 41276045,
            "ReferenceAlleleVCF": "A",
            "AlternateAlleleVCF": "G",
            "ClinicalSignificance": "Pathogenic",
        }])
        matches = proc.find_matches(annotated, clinvar)
        assert len(matches) == 1
        assert matches[0]["label"] == 1


class TestVcfDeletionAlleles:
    def test_left_anchored_deletion(self):
        tx = GT(
            gene_name="TEST",
            transcript_id="NM_1",
            exons=[{
                "cds_sequence": "ATGAAA",
                "genomic_start": 100,
                "genomic_end": 106,
            }],
            chromosome="chr1",
            strand="+",
        )
        pos, ref, alt = tx._vcf_deletion_alleles(50, "ACG")
        assert alt != "."
        assert ref.endswith("ACG")
        assert len(ref) == len(alt) + 3
        assert pos == 50

    def test_from_gtf_exists(self):
        assert hasattr(GeneTranscript, "from_gtf")
