#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
MATCHVAR table annotation tool
"""

import os
import sys
import argparse
import subprocess
import logging
import re
import locale
import io
from typing import Dict, List, Tuple, Optional

def _detect_python_executable() -> str:
    """Detect a suitable Python interpreter to invoke child scripts.
    Priority:
    1) Environment variable PYTHON_EXECUTABLE
    2) Project local .venv (Windows/Linux/macOS)
    3) Current interpreter (sys.executable)
    """
    # 1) Explicit env var
    env_py = os.environ.get('PYTHON_EXECUTABLE')
    if env_py and os.path.exists(env_py):
        return env_py

    # 2) Local .venv (package dir parent, then cwd)
    pkg_root = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
    for project_root in (pkg_root, os.getcwd()):
        venv_unix = os.path.join(project_root, '.venv', 'bin', 'python')
        venv_win = os.path.join(project_root, '.venv', 'Scripts', 'python.exe')
        if os.name == 'nt' and os.path.exists(venv_win):
            return venv_win
        if os.path.exists(venv_unix):
            return venv_unix

    # 3) Fallback to current interpreter
    return sys.executable


def _resolve_humandb_file(dbloc: str, dbname: str, buildver: str, kind: str = "txt") -> Optional[str]:
    """优先裸库名，兼容 {buildver}_* 旧文件名；基因协议还可回退 bigBed。"""
    try:
        from resources_layout import resolve_humandb_db_file  # optional platform helper
        p = resolve_humandb_db_file(dbloc, dbname, buildver, kind=kind)
        if p:
            return str(p)
    except Exception:
        pass
    try:
        from . import resource_files as rf
    except ImportError:
        here = os.path.dirname(os.path.abspath(__file__))
        if here not in sys.path:
            sys.path.insert(0, here)
        import resource_files as rf  # type: ignore
    if kind == "mrna":
        found = rf.resolve_mrna_fasta(dbloc, dbname, buildver)
        if found:
            return found
    elif kind == "gene":
        return rf.ensure_gene_pred(dbloc, dbname, buildver)
    elif dbname == "cytoBand" or kind == "region":
        found = rf.resolve_region_file(dbloc, dbname, buildver)
        if found:
            return found
    for stem in (dbname, f"{buildver}_{dbname}"):
        if kind == "mrna":
            for suf in (".fa", ".fasta"):
                cand = os.path.join(dbloc, f"{stem}Mrna{suf}")
                if os.path.exists(cand):
                    return cand
        else:
            if kind == "plain":
                sufs = (".txt",)
            elif kind in ("gz", "any"):
                sufs = (".txt.gz", ".txt")
            else:
                sufs = (".txt", ".txt.gz")
            for suf in sufs:
                cand = os.path.join(dbloc, f"{stem}{suf}")
                if os.path.exists(cand):
                    return cand
    if kind in ("plain", "txt", "any") and dbname in ("refGene", "ensGene", "ncbiRefSeq", "gencode"):
        return rf.ensure_gene_pred(dbloc, dbname, buildver)
    return None

PYTHON_EXECUTABLE = _detect_python_executable()

def _query_fmt_mod():
    try:
        from . import query_format as _m
        return _m
    except ImportError:
        pass
    here = os.path.dirname(os.path.abspath(__file__))
    if here not in sys.path:
        sys.path.insert(0, here)
    import query_format as _m  # type: ignore
    return _m


def _mane_mod():
    """加载 mane_transcripts：兼容包导入与脚本子进程两种运行方式。"""
    try:
        from . import mane_transcripts as _m
        return _m
    except ImportError:
        pass
    here = os.path.dirname(os.path.abspath(__file__))
    if here not in sys.path:
        sys.path.insert(0, here)
    import mane_transcripts as _m  # type: ignore
    return _m

logger = logging.getLogger(__name__)

def get_system_encoding():
    """Get the system encoding, ensuring cross-platform compatibility"""
    try:
        # Try to get the system default encoding
        system_encoding = locale.getpreferredencoding()
        # If the encoding is not available, use UTF-8
        if system_encoding.lower() in ['cp1252', 'gbk', 'gb2312']:
            return 'utf-8'
        return system_encoding
    except:
        return 'utf-8'

def run_subprocess_safe(command: str, **kwargs):
    """Run subprocess safely, handle encoding issues, and support real-time output"""
    # Set the default encoding parameters
    default_kwargs = {
        'shell': True,
        'text': True,
        'encoding': get_system_encoding(),
        'errors': 'replace'  # Replace the undecodable characters with placeholders
    }
    
    # Check if real-time output is needed
    real_time_output = kwargs.pop('real_time_output', True)
    
    if real_time_output:
        # Real-time output mode: do not capture output, directly display
        default_kwargs.update({
            'stdout': None,  # Directly output to the console
            'stderr': None,  # Directly output to the console
        })
    else:
        # Capture output mode: used for scenarios that require a return result
        default_kwargs.update({
            'capture_output': True,
        })
    
    # Merge other parameters
    default_kwargs.update(kwargs)
    
    try:
        return subprocess.run(command, **default_kwargs)
    except UnicodeDecodeError:
        # If the default encoding fails, try UTF-8
        default_kwargs['encoding'] = 'utf-8'
        return subprocess.run(command, **default_kwargs)

# Annotation header definition
# table_matchvar.py 既可作为包模块导入，也会被 CLI 直接 python 执行
try:
    from .column_names import (
        annotation_headers_for_gene_dbs,
        format_g_hgvs,
        gene_column_keys,
        gene_column_names,
    )
except ImportError:  # pragma: no cover - script entry
    from column_names import (
        annotation_headers_for_gene_dbs,
        format_g_hgvs,
        gene_column_keys,
        gene_column_names,
    )

ANNOTATION_HEADERS = {
    "ljb_all": ["LJB_PhyloP", "LJB_PhyloP_Pred", "LJB_SIFT", "LJB_SIFT_Pred", "LJB_PolyPhen2", "LJB_PolyPhen2_Pred", "LJB_LRT", "LJB_LRT_Pred", "LJB_MutationTaster", "LJB_MutationTaster_Pred", "LJB_GERP++"],
    "ljb2_all": ["LJB2_SIFT", "LJB2_PolyPhen2_HDIV", "LJB2_PP2_HDIV_Pred", "LJB2_PolyPhen2_HVAR", "LJB2_PolyPhen2_HVAR_Pred", "LJB2_LRT", "LJB2_LRT_Pred", "LJB2_MutationTaster", "LJB2_MutationTaster_Pred", "LJB_MutationAssessor", "LJB_MutationAssessor_Pred", "LJB2_FATHMM", "LJB2_GERP++", "LJB2_PhyloP", "LJB2_SiPhy"],
    "popfreq_all": ["PopFreqMax", "1000G2012APR_ALL", "1000G2012APR_AFR", "1000G2012APR_AMR", "1000G2012APR_ASN", "1000G2012APR_EUR", "ESP6500si_ALL", "ESP6500si_AA", "ESP6500si_EA", "CG46"],
    # MATCHVAR gene annotation columns (Function/Gene/cHGVS/ExonicEffect/VarType/pHGVS)
    **annotation_headers_for_gene_dbs(),
    # Add AlphaMissense annotation fields
    "AlphaMissense": ["am_pathogenicity", "am_class"]
}

class TableAnnotator:
    """Table annotator - complete implementation"""
    
    def __init__(self, queryfile: str, dbloc: str, **kwargs):
        self.queryfile = queryfile
        self.dbloc = dbloc
        self.outfile = kwargs.get('outfile')
        self.buildver = kwargs.get('buildver', 'hg19')
        self.remove = kwargs.get('remove', False)
        self.protocol = kwargs.get('protocol')
        self.operation = kwargs.get('operation')
        self.otherinfo = kwargs.get('otherinfo', False)
        self.nastring = kwargs.get('nastring')
        self.argument = kwargs.get('argument')
        self.thread = kwargs.get('thread')
        # Protein polish is on by default; -nopolish disables it. -polishgene is kept for CLI compat.
        if kwargs.get('nopolish'):
            self.polish = False
        else:
            self.polish = kwargs.get('polish', True)
        self.mane_file = kwargs.get('mane_file')
        self.use_mane_transcript = kwargs.get('use_mane_transcript', False)
        # PLATFORM_UNUSED kwargs — 平台 MatchvarRunner 命令行未传入（独立 CLI 仍可用）
        self.csvout = kwargs.get('csvout', False)
        self.vcfinput = kwargs.get('vcfinput', False)
        # Legacy ANNOVAR-style path: convert2matchvar → .mvinput → annotate
        self.convertvcf = kwargs.get('convertvcf', False)
        self.dot2underline = kwargs.get('dot2underline', False)
        self.intronhgvs = kwargs.get('intronhgvs', False)
        self.verbose = kwargs.get('verbose', False)
        self.man = kwargs.get('man', False)
        self.checkfile = kwargs.get('checkfile', False)
        self.onetranscript = kwargs.get('onetranscript', False)
        self.genericdbfile = kwargs.get('genericdbfile')
        self.gff3dbfile = kwargs.get('gff3dbfile')
        self.bedfile = kwargs.get('bedfile')
        self.vcfdbfile = kwargs.get('vcfdbfile')
        self.tempdir = kwargs.get('tempdir')
        self.xreffile = kwargs.get('xreffile')
        self.convertarg = kwargs.get('convertarg')
        self.codingarg = kwargs.get('codingarg')
        self.filter = kwargs.get('filter', False)
        self.regionanno = kwargs.get('regionanno', False)
        self.geneanno = kwargs.get('geneanno', False)
        
        # Internal variables
        self.unlink_files = []
        self.header = []
        self.varanno = {}
        self.protocols = []
        self.operations = []
        self.arguments = []
        self.dbtype1 = []
        
        # Preload MANE transcript mapping (only load when use_mane_transcript is True)
        if self.use_mane_transcript:
            self.mane_transcripts = self._load_mane_transcripts()
            self._mane_all_ids = _mane_mod().all_mane_base_ids(self.mane_transcripts)
            if not self.mane_transcripts:
                logger.warning(
                    "已启用 MANE 过滤但未加载到任何映射；将无法按 MANE 转录本过滤"
                )
                self._mane_all_ids = set()
            else:
                logger.info(
                    "MANE 全局转录本 ID 数: %s", len(self._mane_all_ids),
                )
        else:
            self.mane_transcripts = {}
            self._mane_all_ids = set()
        
        # Ensure outfile and tempfile use the same directory as the input file
        if self.outfile and not os.path.dirname(self.outfile):
            # If outfile is only a file name, use the directory of the input file
            input_dir = os.path.dirname(os.path.abspath(self.queryfile))
            self.outfile = os.path.join(input_dir, self.outfile)
            self.tempfile = self.outfile
        else:
            self.tempfile = self.outfile
        
        # Process arguments
        self._process_arguments()
    
    def _process_arguments(self):
        """Process arguments"""
        # Set default values
        if not self.outfile:
            self.outfile = self.queryfile
        
        # Verify required parameters
        if not self.protocol:
            raise ValueError("Error: --protocol is required")
        if not self.operation:
            raise ValueError("Error: --operation is required")
        
        # Verify VCF input related parameters
        if (self.vcfinput or self.convertvcf) and self.csvout:
            raise ValueError("Error in argument: -csvout is not compatible with -vcfinput / -convertvcf")
        
        if self.convertvcf and self.vcfinput:
            logger.info(
                "NOTICE: both -convertvcf and -vcfinput given; "
                "using convert2matchvar → .mvinput, then annotating"
            )
        
        # Verify file related parameters
        if self.genericdbfile and not (self.filter or self.regionanno):
            raise ValueError("Error in argument: the --genericdbfile argument is supported only for the --filter and --region operation")
        
        if self.gff3dbfile and not (self.geneanno or self.regionanno):
            raise ValueError("Error in argument: the --gff3dbfile argument is supported only for the --geneanno or --regionanno operation")
        
        if self.bedfile and not self.regionanno:
            raise ValueError("Error in argument: the --bedfile argument is supported only for the --regionanno operation")
        
        if self.vcfdbfile and not self.filter:
            raise ValueError("Error in argument: the --vcfdbfile argument is supported only for the --filter operation")
        
        # VCF input (native or convert): force -nastring '.' and enable -otherinfo
        if self.vcfinput or self.convertvcf:
            if self.nastring is not None and self.nastring != '.':
                raise ValueError("Error in argument: -nastring must be '.' when '-vcfinput' or '-convertvcf' is specified")
            self.nastring = '.'
            self.otherinfo = True
        else:
            if self.nastring is None:
                self.nastring = '.'

        # # Verify thread parameters
        # if self.thread and self.thread > self.maxgenethread:
        #     logger.info(f"NOTICE: number of threads is reduced to {self.maxgenethread}")
        #     self.thread = self.maxgenethread
        
        # Process protocols and operations
        self.protocols = [p.strip() for p in self.protocol.split(',')] if self.protocol else []
        self.operations = [o.strip() for o in self.operation.split(',')] if self.operation else []
        self.arguments = [a.strip() for a in self.argument.split(',')] if self.argument else []
        
        # Ensure the length of the parameter list matches
        while len(self.arguments) < len(self.protocols):
            self.arguments.append('')
        
        # Verify operation types
        valid_operations = ['g', 'r', 'f', 'gx']
        for op in self.operations:
            if op not in valid_operations:
                raise ValueError(f"Error: invalid operation '{op}'. Valid operations are: {', '.join(valid_operations)}")
        
        # Check if gx operation requires xreffile
        if 'gx' in self.operations and not self.xreffile:
            logger.warning("WARNING: the 'g' rather than 'gx' operation will be used due to lack of -xreffile argument")
            # Replace gx with g
            self.operations = ['g' if op == 'gx' else op for op in self.operations]
    
    def _check_file_existence(self, dbtype1_list: List[str]):
        """Check if the database files exist"""
        if not self.checkfile:
            return
        
        logger.info("Checking database file existence...")
        missing_files = []
        
        for dbtype1 in dbtype1_list:
            db_file = _resolve_humandb_file(self.dbloc, dbtype1, self.buildver, kind="any")
            if not db_file:
                missing_files.append(f"{dbtype1}.txt (or {self.buildver}_{dbtype1}.txt) in {self.dbloc}")
                logger.warning(f"Database file does not exist for {dbtype1} in {self.dbloc}")
        
        if missing_files:
            logger.error(f"Found {len(missing_files)} missing database files:")
            for file in missing_files:
                logger.error(f"  - {file}")
            raise FileNotFoundError(f"Missing required database files, please check the database directory: {self.dbloc}")
        else:
            logger.info("All database files checked")
    
    def _proxy_db_type(self, protocols: List[str]) -> List[str]:
        """Convert protocols to database types"""
        dbtype1 = []
        for protocol in protocols:
            if protocol in ['gene', 'refgene']:
                dbtype1.append('refGene')
            elif protocol == 'knowngene':
                dbtype1.append('knownGene')
            elif protocol == 'ensgene':
                dbtype1.append('ensGene')
            elif protocol.lower() == 'ncbirefseq':
                dbtype1.append('ncbiRefSeq')
            elif protocol.lower() == 'gencode':
                dbtype1.append('gencode')
            else:
                dbtype1.append(protocol)
        return dbtype1
    
    def run_annotation(self):
        """Run annotation"""
        logger.info("Starting table annotation...")
        
        # Convert database types
        self.dbtype1 = self._proxy_db_type(self.protocols)
        
        # Check file existence
        self._check_file_existence(self.dbtype1)
        
        # Two query paths:
        #   native VCF / 4-col  → parse CHROM POS REF ALT in annotate_variation
        #   -convertvcf         → convert2matchvar → .mvinput → annotate
        #   existing .mvinput   → annotate as-is
        if self.convertvcf:
            self._convert_vcf_to_mvinput()
        else:
            self._resolve_query_format()
        self._run_standard_annotation()
        
        # Print original output
        self._print_original_output()
        
        # Clean up temporary files
        if self.remove:
            self._cleanup_temp_files()
    
    def _resolve_query_format(self):
        """Detect VCF vs MATCHVAR 5-column input. ``-vcfinput`` forces VCF (native)."""
        if self.vcfinput:
            self.query_format = "vcf"
        else:
            try:
                self.query_format = _query_fmt_mod().sniff_query_format(self.queryfile)
            except Exception:
                self.query_format = "mvinput"
        if self.query_format == "vcf":
            logger.info(
                "NOTICE: annotating VCF CHROM/POS/REF/ALT directly "
                "(use -convertvcf to go through convert2matchvar → .mvinput)"
            )

    def _convert_vcf_to_mvinput(self):
        """Legacy path: VCF → 5-column .mvinput via convert2matchvar, then annotate that file."""
        if self.csvout:
            raise ValueError("Error: -csvout is not compatible with -convertvcf")

        sniffed = None
        try:
            sniffed = _query_fmt_mod().sniff_query_format(self.queryfile)
        except Exception:
            sniffed = None
        if sniffed == "mvinput":
            logger.info(
                "NOTICE: -convertvcf ignored because query is already MATCHVAR 5-column input"
            )
            self.query_format = "mvinput"
            return

        current_dir = os.path.dirname(os.path.abspath(__file__))
        convert2matchvar_script = os.path.join(current_dir, "convert2matchvar.py")
        if not os.path.exists(convert2matchvar_script):
            raise FileNotFoundError(f"convert2matchvar.py not found: {convert2matchvar_script}")

        convertarg_str = f"{self.convertarg} " if self.convertarg else ""
        mv_path = f"{self.tempfile}.mvinput"
        sc = (
            f"{PYTHON_EXECUTABLE} {convert2matchvar_script} {convertarg_str}"
            f"-includeinfo -allsample -withfreq -format vcf4 "
            f"{self.queryfile} > {mv_path}"
        )
        logger.info(f"NOTICE: converting VCF → .mvinput with <{sc}>")

        result = run_subprocess_safe(sc)
        if result.returncode != 0:
            raise RuntimeError(f"Error running system command: <{sc}>")
        if not os.path.exists(mv_path) or os.path.getsize(mv_path) == 0:
            raise RuntimeError(f"convert2matchvar produced no output: {mv_path}")

        self.queryfile = mv_path
        self.query_format = "mvinput"
        self.unlink_files.append(mv_path)
        logger.info(f"NOTICE: converted query is {mv_path}; annotating as MATCHVAR 5-column input")
    
    def _annotate_variation_query_flag(self) -> str:
        if getattr(self, "query_format", None) == "vcf":
            return " -vcfinput"
        return ""
    
    def _run_standard_annotation(self):
        """Run standard annotation"""
        for i, (protocol, operation) in enumerate(zip(self.protocols, self.operations)):
            logger.info("-----------------------------------------------------------------")
            logger.info(f"NOTICE: Processing operation={operation} protocol={protocol}")
            
            try:
                if operation in ['g', 'gx']:
                    self._gene_operation(protocol, self.dbtype1[i], self.arguments[i] if i < len(self.arguments) else None, operation)
                elif operation == 'r':
                    self._region_operation(protocol, self.dbtype1[i], self.arguments[i] if i < len(self.arguments) else None)
                elif operation == 'f':
                    self._filter_operation(protocol, self.dbtype1[i], self.arguments[i] if i < len(self.arguments) else None)
                
                logger.info(f"Successfully completed protocol {protocol}")
                
            except Exception as e:
                logger.error(f"Failed to process protocol {protocol}: {e}")
                logger.warning(f"Continuing with next protocol...")
                # Continue with the next protocol instead of stopping
    
    def _gene_operation(self, protocol: str, dbtype1: str, argument: str, operation: str):
        """Gene annotation operation"""
        # Process protocol names
        genetype = {
            'gene': 'refGene',
            'refgene': 'refGene',
            'knowngene': 'knownGene',
            'ensgene': 'ensGene',
            'ncbirefseq': 'ncbiRefSeq',
            'gencode': 'gencode',
        }
        if protocol in genetype:
            protocol = genetype[protocol]
        
        # Get the absolute path of the current script
        current_dir = os.path.dirname(os.path.abspath(__file__))
        annotate_variation_script = os.path.join(current_dir, 'annotate_variation.py')
        
        # Build command, using the Python interpreter in the virtual environment
        sc = f"{PYTHON_EXECUTABLE} {annotate_variation_script} -geneanno -buildver {self.buildver} -dbtype {protocol} -outfile {self.tempfile}.{protocol} -exonsort -nofirstcodondel {self.queryfile} {self.dbloc}{self._annotate_variation_query_flag()}"
        
        # Add splicing_threshold parameter
        if self.intronhgvs:
            sc += f" -splicing_threshold {self.intronhgvs}"
        
        # Add MANE transcript mapping parameter (only when use_mane_transcript is True)
        if self.use_mane_transcript and self.mane_transcripts:
            # Write MANE transcript mapping to temporary file
            mane_file = f"{self.tempfile}.{protocol}.mane"
            try:
                _mane_mod().write_mane_tsv(self.mane_transcripts, mane_file)
                sc += f" -mane_file {mane_file}"
                sc += f" -use_mane_transcript"
                self.unlink_files.append(mane_file)
                logger.debug(f"Added MANE transcript filtering for protocol {protocol}")
            except Exception as e:
                logger.warning(f"Failed to write MANE transcript mapping file: {e}")
        
        if argument and argument.strip():
            sc += f" {argument}"
        
        if self.thread:
            sc += f" -thread {self.thread}"
        
        logger.info(f"NOTICE: Running with system command <{sc}>")
        result = run_subprocess_safe(sc)
        if result.returncode != 0:
            raise RuntimeError(f"Error running system command: <{sc}>")

        # Protein annotation optimization based on mRNA (using {buildver}_{protocol}Mrna.fa)
        anno_outfile = f"{self.tempfile}.{protocol}.variant_function"
        e_anno_outfile = f"{self.tempfile}.{protocol}.exonic_variant_function"

        if self.polish:
            try:
                # Prepare the required file paths (裸库名优先，兼容旧 {buildver}_*，以及 bigBed / 2bit)
                gene_file = _resolve_humandb_file(self.dbloc, protocol, self.buildver, kind="plain")
                if not gene_file:
                    gene_file = _resolve_humandb_file(self.dbloc, protocol, self.buildver, kind="gene")
                mrna_fa = _resolve_humandb_file(self.dbloc, protocol, self.buildver, kind="mrna")
                coding_change_script_py = os.path.join(current_dir, 'coding_change.py')

                # Only execute polish step when the key files exist
                if not gene_file:
                    logger.warning(f"Gene definition file does not exist, skipping polish step: {protocol} in {self.dbloc}")
                else:
                    if not mrna_fa:
                        # Splice transcripts from genome 2bit/FASTA for IDs in this EVF
                        try:
                            from .resource_files import (
                                GenomeSequence,
                                extract_mrna_fasta,
                                parse_evf_transcript_ids,
                            )
                        except ImportError:
                            from resource_files import (  # type: ignore
                                GenomeSequence,
                                extract_mrna_fasta,
                                parse_evf_transcript_ids,
                            )
                        try:
                            genome = GenomeSequence.from_dbloc(self.dbloc, self.buildver)
                            mrna_fa = f"{self.tempfile}.{protocol}.spliced.fa"
                            extract_mrna_fasta(
                                gene_file,
                                genome,
                                mrna_fa,
                                transcript_ids=parse_evf_transcript_ids(e_anno_outfile),
                            )
                            genome.close()
                            self.unlink_files.append(mrna_fa)
                            logger.info(f"Built spliced mRNA FASTA from genome sequence: {mrna_fa}")
                        except Exception as exc:
                            logger.warning(
                                f"mRNA FASTA missing and genome splice failed ({exc}); skipping polish"
                            )
                            mrna_fa = None
                    if not mrna_fa:
                        logger.warning(
                            f"mRNA FASTA file does not exist, skipping polish step: {protocol}Mrna.fa in {self.dbloc}"
                        )
                    else:
                        # Rename the original exonic_variant_function to .orig
                        e_anno_outfile_orig = f"{e_anno_outfile}.orig"
                        try:
                            os.rename(e_anno_outfile, e_anno_outfile_orig)
                        except Exception as e:
                            logger.warning(f"Failed to rename {e_anno_outfile} -> {e_anno_outfile_orig}, skipping polish step: {e}")
                            e_anno_outfile_orig = None

                        if e_anno_outfile_orig:
                            sc_cc = (
                                f"{PYTHON_EXECUTABLE} {coding_change_script_py} "
                                f"{self.codingarg + ' ' if self.codingarg else ''}"
                                f"{e_anno_outfile_orig} {gene_file} {mrna_fa} "
                                f"-includesnp -alltranscript "
                                f"-outfile {self.tempfile}.{protocol}.fa "
                                f"-newevf {e_anno_outfile}"
                            )
                            logger.info(f"NOTICE: Running with system command <{sc_cc}>")
                            result_cc = run_subprocess_safe(sc_cc, real_time_output=False)
                            if result_cc.returncode != 0:
                                logger.error(f"coding_change failed, return code: {result_cc.returncode}")
                                logger.error(f"stdout: {getattr(result_cc, 'stdout', '')}")
                                logger.error(f"stderr: {getattr(result_cc, 'stderr', '')}")
                                logger.warning(f"Polish step failed due to coding_change error, using unmodified exonic annotation")

                            self.unlink_files.append(f"{self.tempfile}.{protocol}.fa")
                            self.unlink_files.append(e_anno_outfile_orig)
            except Exception as e:
                logger.warning(f"Polish step failed, using unmodified exonic annotation: {e}")
        
        # Set the header — MATCHVAR gene columns
        self.header.extend(gene_column_names(protocol, underline=self.dot2underline))
        
        # Read the annotation results
        self._read_gene_annotation(anno_outfile, e_anno_outfile, protocol)
    
    def _read_gene_annotation(self, anno_outfile: str, e_anno_outfile: str, protocol: str):
        """Read the gene annotation results"""
        try:
            with open(anno_outfile, 'r', encoding='utf-8') as f:
                for line in f:
                    line = line.strip()
                    if not line:
                        continue
                    
                    parts = line.split('\t')
                    if len(parts) >= 3:
                        function = parts[0]
                        gene = parts[1]
                        varstring = '\t'.join(parts[2:7])  # Chr Start End Ref Alt
                        
                        # Process GeneDetail information
                        gene_detail = ''
                        aa_change = 'p.?'
                        gene_name = gene
                        
                        # Extract the transcript information in parentheses
                        if '(' in gene and ')' in gene:
                            import re
                            # Extract all the content in parentheses
                            transcript_matches = re.findall(r'\(([^)]+)\)', gene)
                            # Remove the content in parentheses, keep the gene name
                            gene_name = re.sub(r'\([^)]+\)', '', gene).strip()
                            
                            if transcript_matches:
                                # Always use all transcript information first
                                gene_detail = ','.join(transcript_matches)
                                # upstream/downstream 的 dist=N 不是 cHGVS，不要写入该列
                                if re.fullmatch(
                                    r'(?:dist=\d+)(?:,(?:dist=\d+))*',
                                    gene_detail.replace(' ', ''),
                                ):
                                    gene_detail = ''
                                
                                # MANE 模式：按「全局 MANE 转录本 ID」过滤（覆盖多基因重叠行）
                                if gene_detail and self.use_mane_transcript and self._mane_all_ids:
                                    mt = _mane_mod()
                                    # 多基因名时无法用单一 gene 查表，改用全局 ID 集
                                    if ',' in gene_name or gene_name not in self.mane_transcripts:
                                        filtered = mt.filter_parts_by_global_mane_ids(
                                            [p for p in transcript_matches if p.strip()],
                                            self._mane_all_ids,
                                            joiner=',',
                                        )
                                    else:
                                        filtered = mt.filter_annotation_parts(
                                            transcript_matches,
                                            self.mane_transcripts[gene_name],
                                            joiner=',',
                                        )
                                    # 无 MANE 匹配时回退保留原转录本
                                    gene_detail = filtered or gene_detail

                        # 基因名去重（upstream 等多异构体残留：PLEKHN1,PLEKHN1,PLEKHN1）
                        if ',' in (gene_name or ''):
                            gparts = [p.strip() for p in gene_name.split(',') if p.strip()]
                            gene_name = ','.join(dict.fromkeys(gparts))
                        if ',' in (function or ''):
                            fparts = [p.strip() for p in function.split(',') if p.strip()]
                            if len(fparts) > 1 and len(set(fparts)) == 1:
                                function = fparts[0]

                        # Process GeneDetail and AAChange for intronic variants
                        if function == 'intronic' and gene_detail:
                            # For intronic variants, gene_detail contains the full annotation
                            # Add :p.? if not present
                            if ':p.' not in gene_detail:
                                gene_detail = gene_detail + ':p.?'
                            # AAChange should be p.? for intronic variants
                            aa_change = 'p.?'
                        
                        # Store the annotation - MATCHVAR column names
                        if varstring not in self.varanno:
                            self.varanno[varstring] = {}
                        
                        cols = gene_column_keys(protocol, underline=self.dot2underline)
                        self.varanno[varstring].update({
                            cols["function"]: function,
                            cols["gene"]: gene_name,
                            cols["c_hgvs"]: gene_detail,
                            cols["mane_select"]: '',
                            cols["exonic_effect"]: 'NA',
                            cols["vartype"]: 'NA',
                            cols["p_hgvs"]: aa_change,
                        })
        
        except Exception as e:
            logger.error(f"Error reading gene annotation file: {e}")
        
        # Read the exonic annotation file
        try:
            with open(e_anno_outfile, 'r', encoding='utf-8') as f:
                for line in f:
                    line = line.strip()
                    if not line:
                        continue
                    
                    parts = line.split('\t')
                    if len(parts) >= 3:
                        if parts[0].startswith('line'):
                            if len(parts) >= 8:
                                exonic_function = parts[1]
                                aa_change = parts[2]
                                varstring = '\t'.join(parts[3:8])
                            else:
                                continue
                        else:
                            exonic_function = parts[0]
                            aa_change = parts[1]
                            varstring = '\t'.join(parts[2:7])  # Chr Start End Ref Alt
                        
                        # Update the exonic annotation information + variant type (VarType)
                        if varstring in self.varanno:
                            cols = gene_column_keys(protocol, underline=self.dot2underline)
                            # Extract p.HGVS from cHGVS for pHGVS column
                            gene_detail = self.varanno[varstring].get(cols["c_hgvs"], '')
                            
                            # Process pHGVS: extract p.HGVS from aa_change variable
                            aa_change_final = 'p.?'
                            if aa_change and self.use_mane_transcript and self._mane_all_ids:
                                # 外显子注释也可能含多转录本 / 多基因，按全局 MANE ID 过滤
                                aa_parts = [p.strip() for p in aa_change.split(',') if p.strip()]
                                filtered_aa = _mane_mod().filter_parts_by_global_mane_ids(
                                    aa_parts, self._mane_all_ids, joiner=',',
                                )
                                # 无 MANE 匹配时回退保留原 aa_change
                                aa_change = filtered_aa or aa_change
                            if aa_change and ':p.' in aa_change:
                                # Extract p.HGVS from all transcripts in aa_change
                                import re
                                p_matches = re.findall(r':p\.([^:\s,]+)', aa_change)
                                if p_matches:
                                    # Join all p.HGVS parts with comma
                                    aa_change_final = ','.join([f"p.{p}" for p in p_matches])
                                else:
                                    # If no p.HGVS found, use p.?
                                    aa_change_final = 'p.?'
                            elif self.use_mane_transcript and not aa_change:
                                aa_change_final = 'p.?'
                            
                            # Write temporary ExonicEffect hint + pHGVS（最终类型在 _finalize_effect_columns 统一命名）
                            self.varanno[varstring][cols["exonic_effect"]] = exonic_function
                            self.varanno[varstring][cols["p_hgvs"]] = aa_change_final

                            # Refine cHGVS from EVF transcript detail when available
                            try:
                                detail_candidate = aa_change or ''
                                if detail_candidate and ':' in detail_candidate:
                                    # Normalize c.hgvs fragment
                                    try:
                                        import re
                                        # Extract and normalize c.hgvs
                                        m = re.search(r'(c\.[^:\s]+)', detail_candidate)
                                        if m:
                                            try:
                                                from .coding_change import normalize_c_hgvs
                                            except ImportError:
                                                from coding_change import normalize_c_hgvs
                                            norm_c = normalize_c_hgvs(m.group(1))
                                            detail_candidate = detail_candidate.replace(m.group(1), norm_c)
                                    except Exception:
                                        pass
                                    
                                    # Process cHGVS: keep full transcript detail format
                                    gene_detail_final = detail_candidate
                                    # If c. exists but lacks :p., add a placeholder p.?
                                    if ('c.' in gene_detail_final) and (':p.' not in gene_detail_final):
                                        gene_detail_final = gene_detail_final + ':p.?'
                                    # Only when c. or :p. is present is it considered a complete detail
                                    if ('c.' in gene_detail_final) or (':p.' in gene_detail_final):
                                        # 多转录本用逗号分隔
                                        gene_detail_final = gene_detail_final.replace(';', ',')
                                        self.varanno[varstring][cols["c_hgvs"]] = gene_detail_final
                            except Exception:
                                pass
        
        except Exception as e:
            logger.error(f"Error reading exonic annotation file: {e}")

        # 统一按区域规则重写 ExonicEffect / VarType（与 Gene、cHGVS 条数逗号对齐）
        self._finalize_effect_columns(protocol)

    def _finalize_effect_columns(self, protocol: str) -> None:
        """根据 Function / Gene / cHGVS 重写 ExonicEffect、VarType，并填充 MANE Select。"""
        try:
            try:
                from .variant_typing import classify_row, split_annotation_parts
            except ImportError:
                from variant_typing import classify_row, split_annotation_parts
        except Exception as e:
            logger.warning(f"variant_typing unavailable: {e}")
            return

        # MANE Select 列需要转录本目录（即使未启用过滤也尽量加载）
        if not getattr(self, '_mane_all_ids', None):
            try:
                if not getattr(self, 'mane_transcripts', None):
                    self.mane_transcripts = self._load_mane_transcripts()
                self._mane_all_ids = _mane_mod().all_mane_base_ids(self.mane_transcripts or {})
            except Exception:
                self._mane_all_ids = set()

        cols = gene_column_keys(protocol, underline=self.dot2underline)
        for varstring, anno in list(self.varanno.items()):
            try:
                toks = varstring.split('\t')
                ref = toks[3] if len(toks) >= 5 else ''
                alt = toks[4] if len(toks) >= 5 else ''
                function = anno.get(cols["function"], '') or ''
                gene = anno.get(cols["gene"], '') or ''
                c_hgvs = (anno.get(cols["c_hgvs"], '') or '').replace(';', ',')
                exonic_hint = anno.get(cols["exonic_effect"], '') or ''
                p_hgvs = anno.get(cols["p_hgvs"], '') or ''

                effect, vartype = classify_row(
                    function=function,
                    gene=gene,
                    c_hgvs=c_hgvs,
                    ref=ref,
                    alt=alt,
                    exonic_effect=exonic_hint,
                    p_hgvs=p_hgvs,
                )
                anno[cols["c_hgvs"]] = c_hgvs
                anno[cols["exonic_effect"]] = effect
                anno[cols["vartype"]] = vartype
                anno[cols["mane_select"]] = self._mane_select_flags(c_hgvs, split_annotation_parts)
            except Exception as e:
                logger.debug(f"finalize effect columns failed for {varstring}: {e}")

    def _mane_select_flags(self, c_hgvs: str, split_fn) -> str:
        """按 cHGVS 各转录本片段生成 yes/no（逗号对齐）。"""
        parts = split_fn(c_hgvs)
        if not parts:
            return ''
        allowed = getattr(self, '_mane_all_ids', None) or set()
        flags = []
        for part in parts:
            tid = ''
            try:
                m = re.search(r'\b((?:NM|NR|XM|XR|ENST)_[0-9]+(?:\.[0-9]+)?)\b', part or '')
                if m:
                    tid = m.group(1)
            except Exception:
                tid = ''
            if not tid or not allowed:
                flags.append('no')
                continue
            try:
                bid = _mane_mod().base_transcript_id(tid)
                flags.append('yes' if bid and bid in allowed else 'no')
            except Exception:
                flags.append('no')
        return ','.join(flags)
    
    def _region_operation(self, protocol: str, dbtype1: str, argument: str):
        """Region annotation operation"""
        # Get the absolute path of the current script
        current_dir = os.path.dirname(os.path.abspath(__file__))
        annotate_variation_script = os.path.join(current_dir, 'annotate_variation.py')
        
        sc = f"{PYTHON_EXECUTABLE} {annotate_variation_script} -regionanno -dbtype {protocol} -buildver {self.buildver} -outfile {self.tempfile} {self.queryfile} {self.dbloc}{self._annotate_variation_query_flag()}"
        
        if argument and argument.strip():
            sc += f" {argument}"
        
        if self.thread:
            sc += f" -thread {self.thread}"
        
        logger.info(f"NOTICE: Running with system command <{sc}>")
        result = run_subprocess_safe(sc)
        if result.returncode != 0:
            raise RuntimeError(f"Error running system command: <{sc}>")
        
        # Set the header
        header = protocol
        self.header.append(header)
        
        # Read the annotation results
        region_file = f"{self.tempfile}.{self.buildver}_{dbtype1}"
        try:
            with open(region_file, 'r', encoding='utf-8') as f:
                for line in f:
                    line = line.strip()
                    if not line:
                        continue
                    
                    parts = line.split('\t')
                    if len(parts) >= 3:
                        db = parts[0]
                        anno = parts[1]
                        varstring = '\t'.join(parts[2:7])
                        
                        if varstring not in self.varanno:
                            self.varanno[varstring] = {}
                        self.varanno[varstring][header] = anno
        
        except Exception as e:
            logger.error(f"Error reading region annotation file: {e}")
        
        self.unlink_files.append(region_file)
    
    def _filter_operation(self, protocol: str, dbtype1: str, argument: str):
        """Filter operation"""
        # Get the absolute path of the current script
        current_dir = os.path.dirname(os.path.abspath(__file__))
        annotate_variation_script = os.path.join(current_dir, 'annotate_variation.py')
        
        sc = f"{PYTHON_EXECUTABLE} {annotate_variation_script} -filter -dbtype {protocol} -buildver {self.buildver} -outfile {self.tempfile} {self.queryfile} {self.dbloc}{self._annotate_variation_query_flag()}"
        
        if argument and argument.strip():
            sc += f" {argument}"
        
        if self.thread:
            sc += f" -thread {self.thread}"
        
        # Check for both uncompressed and compressed versions (裸库名优先)
        dbfile = _resolve_humandb_file(self.dbloc, protocol, self.buildver, kind="any")
        if dbfile and dbfile.endswith(".gz"):
            tbi_file = dbfile + ".tbi"
            if os.path.exists(tbi_file):
                logger.info(f"Using compressed database with Tabix index: {dbfile}")
            else:
                logger.info(f"Using compressed database (no index): {dbfile}")
        elif dbfile:
            logger.info(f"Using uncompressed database: {dbfile}")
        
        if dbfile and os.path.exists(dbfile):
            try:
                # Handle compressed files
                if dbfile.endswith('.gz'):
                    import gzip
                    with gzip.open(dbfile, 'rt', encoding='utf-8') as f:
                        first_line = f.readline().strip()
                        # Compatible with UTF-8 BOM
                        if first_line.startswith('\ufeff'):
                            first_line = first_line.lstrip('\ufeff')
                        if first_line.startswith('#'):
                            # Read the header information: compatible with TAB or whitespace
                            header_raw = first_line[1:].strip()
                            # Support TAB, comma, or any whitespace
                            if '\t' in header_raw:
                                fields_all = header_raw.split('\t')
                            elif ',' in header_raw:
                                fields_all = header_raw.split(',')
                            else:
                                fields_all = re.split(r'\s+', header_raw)
                            # Remove the first 5 columns (Chr, Start, End, Ref, Alt)
                            fields = fields_all[5:] if len(fields_all) > 5 else []
                            if not fields:
                                logger.warning(
                                    f"WARNING: Parsed zero annotation columns for {protocol}; header_raw='{header_raw}'")
                            ANNOTATION_HEADERS[protocol] = fields
                            logger.info(
                                f"NOTICE: Finished reading {len(fields)} column headers for '-dbtype {protocol}': {fields}")
                            sc += " -otherinfo"
                            self.otherinfo = True  # Set otherinfo flag
                        else:
                            # If there is no header information, record the warning but continue processing
                            logger.warning(f"WARNING: No header found in {protocol} database, but expected to have one. Using default field name.")
                            ANNOTATION_HEADERS[protocol] = [protocol]
                            sc += " -otherinfo"
                            self.otherinfo = True  # Set otherinfo flag
                else:
                    # Handle uncompressed files
                    with open(dbfile, 'r', encoding='utf-8') as f:
                        first_line = f.readline().strip()
                        # Compatible with UTF-8 BOM
                        if first_line.startswith('\ufeff'):
                            first_line = first_line.lstrip('\ufeff')
                        if first_line.startswith('#'):
                            # Read the header information: compatible with TAB or whitespace
                            header_raw = first_line[1:].strip()
                            # Support TAB, comma, or any whitespace
                            if '\t' in header_raw:
                                fields_all = header_raw.split('\t')
                            elif ',' in header_raw:
                                fields_all = header_raw.split(',')
                            else:
                                fields_all = re.split(r'\s+', header_raw)
                            # Remove the first 5 columns (Chr, Start, End, Ref, Alt)
                            fields = fields_all[5:] if len(fields_all) > 5 else []
                            if not fields:
                                logger.warning(
                                    f"WARNING: Parsed zero annotation columns for {protocol}; header_raw='{header_raw}'")
                            ANNOTATION_HEADERS[protocol] = fields
                            logger.info(
                                f"NOTICE: Finished reading {len(fields)} column headers for '-dbtype {protocol}': {fields}")
                            sc += " -otherinfo"
                            self.otherinfo = True  # Set otherinfo flag
                        else:
                            # If there is no header information, record the warning but continue processing
                            logger.warning(f"WARNING: No header found in {protocol} database, but expected to have one. Using default field name.")
                            ANNOTATION_HEADERS[protocol] = [protocol]
                            sc += " -otherinfo"
                            self.otherinfo = True  # Set otherinfo flag
            except Exception as e:
                logger.error(f"Error reading database file header: {e}")
                ANNOTATION_HEADERS[protocol] = [protocol]
                sc += " -otherinfo"
                self.otherinfo = True  # Set otherinfo flag
        else:
            logger.error(f"Database file not found: {dbfile_txt} (or compressed version with index)")
            ANNOTATION_HEADERS[protocol] = [protocol]
            sc += " -otherinfo"
            self.otherinfo = True  # Set otherinfo flag
        # Special handling (keep avsift's threshold parameters)
        if protocol == 'avsift':
            sc += " -sift_threshold 0"
        logger.info(f"NOTICE: Running system command <{sc}>")
        result = run_subprocess_safe(sc)
        if result.returncode != 0:
            raise RuntimeError(f"Error running system command: <{sc}>")
        # Set the header
        header = protocol
        self.header.append(header)
        # Read the filter results - need to read both filtered and dropped files
        filtered_file = f"{self.tempfile}.{self.buildver}_{dbtype1}_filtered"
        dropped_file = f"{self.tempfile}.{self.buildver}_{dbtype1}_dropped"
        
        # Function to process filter result lines
        def process_filter_line(line, is_filtered=True):
            line = line.strip()
            if not line:
                return None, None
            parts = line.split('\t')
            if len(parts) >= 3:
                db = parts[0]
                # Uniformly process all database annotation fields (assuming all have header information)
                if self.otherinfo or (protocol in ANNOTATION_HEADERS):
                    # More reliable parsing method: determine the number of annotation columns based on the database header field
                    expected_anno_fields = len(ANNOTATION_HEADERS.get(protocol, []))
                    # At least needed: dbtype(1) + annotation columns(expected_anno_fields) + variant 5 columns
                    min_needed = 1 + expected_anno_fields + 5
                    if len(parts) >= min_needed and expected_anno_fields > 0:
                        anno = '\t'.join(parts[1:1 + expected_anno_fields])
                        # Variant key takes the first 5 columns after the variant (Chr, Start, End, Ref, Alt), not the last 5 columns of the line
                        varstring = '\t'.join(parts[1 + expected_anno_fields:1 + expected_anno_fields + 5])
                    else:
                        # Degraded processing: keep compatible with old format (no header or cannot be reliably inferred)
                        if len(parts) >= 6:
                            annotation_count = len(parts) - 6
                            anno = '\t'.join(parts[1:1 + annotation_count]) if annotation_count > 0 else parts[1]
                            # Still try to assemble varstring from the 5 columns immediately after the annotation, avoiding the last 5 columns of the line
                            start_idx = 1 + max(annotation_count, 0)
                            if len(parts) >= start_idx + 5:
                                varstring = '\t'.join(parts[start_idx:start_idx + 5])
                            else:
                                # If it is impossible to determine, fall back
                                varstring = '\t'.join(parts[-5:])
                        else:
                            return None, None
                else:
                    anno = parts[1]
                    varstring = '\t'.join(parts[2:7])
                
                # For dropped variants (not matched), use nastring
                # For filtered variants (matched), keep the actual annotation data
                if is_filtered:
                    anno = str(self.nastring)
                
                return varstring, anno
            return None, None
        
        # Read filtered file (variants that match the database, use actual annotation data)
        try:
            with open(filtered_file, 'r', encoding='utf-8') as f:
                for line in f:
                    varstring, anno = process_filter_line(line, is_filtered=False)
                    if varstring and anno:
                        if varstring not in self.varanno:
                            self.varanno[varstring] = {}
                        # For filtered variants (matched variants), use the actual annotation data
                        self.varanno[varstring][header] = anno
        except Exception as e:
            logger.error(f"Error reading filtered file: {e}")
        
        # Read dropped file (variants that match the database, use actual annotation data)
        try:
            with open(dropped_file, 'r', encoding='utf-8') as f:
                for line in f:
                    varstring, anno = process_filter_line(line, is_filtered=False)
                    if varstring and anno:
                        if varstring not in self.varanno:
                            self.varanno[varstring] = {}
                        # For dropped variants (matched), use the actual annotation data
                        self.varanno[varstring][header] = anno
        except Exception as e:
            logger.error(f"Error reading dropped file: {e}")
        
        # Add both files to cleanup list
        self.unlink_files.append(filtered_file)
        self.unlink_files.append(dropped_file)
    
    def _lookup_varanno(self, rec: dict):
        """Match annotation keys for VCF or 5-column query records."""
        chrom = str(rec.get("chrom", "")).strip()
        start = str(rec.get("start", "")).strip()
        end = str(rec.get("end", "")).strip()
        ref = str(rec.get("ref", "")).strip()
        alt = str(rec.get("alt", "")).strip()
        keys = []
        qline = rec.get("query_line")
        if qline:
            keys.append(str(qline))
        keys.append("\t".join([chrom, start, end, ref, alt]))
        disp = chrom
        if disp and not disp.lower().startswith("chr") and (
            disp.isdigit() or disp.upper() in ("X", "Y", "M", "MT")
        ):
            disp = f"chr{disp}"
            keys.append("\t".join([disp, start, end, ref, alt]))
        elif chrom.lower().startswith("chr"):
            keys.append("\t".join([chrom[3:], start, end, ref, alt]))
        for key in keys:
            if key in self.varanno:
                return self.varanno[key]
        return None

    def _print_original_output(self):
        """Print original output"""
        # Ensure the output file uses the same directory as the input file
        # Fixed output as TSV, file name *_multianno.tsv
        final_out = os.path.abspath(f"{self.outfile}.{self.buildver}_multianno.tsv")
        
        logger.info("-----------------------------------------------------------------")
        logger.info(f"NOTICE: Multianno output file is written to {final_out}")
        
        # Debug information: check varanno dictionary (avoid screen spam, only output the number, change to DEBUG and limit to the first 20)
        logger.info(f"varanno dictionary contains {len(self.varanno)} entries")
        for idx, (varstring, annotations) in enumerate(self.varanno.items()):
            if idx < 20:
                logger.debug(f"  {varstring}: {list(annotations.keys())}")
            elif idx == 20:
                logger.debug("  ... (more entries omitted)")
                break
        
        # Expand the header and remove duplicates
        expanded_header = []
        seen = set()
        for item in self.header:
            if item in ANNOTATION_HEADERS:
                for h in ANNOTATION_HEADERS[item]:
                    if h not in seen:
                        expanded_header.append(h)
                        seen.add(h)
            else:
                if item not in seen:
                    expanded_header.append(item)
                    seen.add(item)
        
        logger.info(f"Expanded header: {expanded_header}")
        
        try:
            base_cols = ['Chr', 'Start', 'End', 'Ref', 'Alt', 'gHGVS']
            header_cols = list(base_cols) + list(expanded_header)
            qfmt = _query_fmt_mod()
            fmt = getattr(self, "query_format", None) or qfmt.sniff_query_format(self.queryfile)

            with open(final_out, 'w', encoding='utf-8', newline='') as f:
                linecount = 0
                for rec in qfmt.iter_query_records(self.queryfile, fmt=fmt):
                    chrom = str(rec.get("chrom", "")).strip()
                    start = str(rec.get("start", "")).strip()
                    end = str(rec.get("end", "")).strip()
                    ref = str(rec.get("ref", "")).strip()
                    alt = str(rec.get("alt", "")).strip()
                    if chrom and not chrom.lower().startswith('chr') and (
                        chrom.isdigit() or chrom.upper() in ('X', 'Y', 'M', 'MT')
                    ):
                        chrom = f"chr{chrom}"

                    extra = [str(x) for x in (rec.get("extra_fields") or [])]
                    info_cols = [re.sub(r'[\t\r\n]+', ' ', x) for x in extra]

                    if linecount == 0:
                        header_cols_write = list(header_cols)
                        if self.otherinfo:
                            for i in range(1, len(info_cols) + 1):
                                header_cols_write.append(f"Otherinfo{i}")
                        f.write("\t".join(header_cols_write) + "\n")
                        expected_cols = len(header_cols_write)

                    oneline = []
                    anno_map = self._lookup_varanno(rec)
                    for item in self.header:
                        if item in ANNOTATION_HEADERS:
                            expanded_field = len(ANNOTATION_HEADERS[item])
                            if anno_map is not None and item in anno_map:
                                values = str(anno_map[item]).split('\t')
                                values = [v.replace('\\x2c', ',').replace('\\x23', '#') for v in values]
                                values = ['.' if v in ['-1', '', 'NA', 'N/A', 'null', 'NULL'] else v for v in values]
                                values = [re.sub(r'[\t\r\n]+', ' ', v) for v in values]
                                if len(values) < expanded_field:
                                    values.extend([str(self.nastring)] * (expanded_field - len(values)))
                                elif len(values) > expanded_field:
                                    values = values[:expanded_field]
                                oneline.extend(values)
                            else:
                                oneline.extend([str(self.nastring)] * expanded_field)
                        else:
                            if anno_map is not None and item in anno_map:
                                cell = re.sub(r'[\t\r\n]+', ' ', str(anno_map[item]))
                                oneline.append(cell)
                            else:
                                oneline.append(str(self.nastring))

                    g_hgvs = format_g_hgvs(chrom, start, end, ref, alt)
                    row_cols = [chrom, start, end, ref, alt, g_hgvs] + oneline
                    if self.otherinfo:
                        n_other = max(0, expected_cols - len(base_cols) - len(expanded_header))
                        if len(info_cols) < n_other:
                            info_cols.extend([''] * (n_other - len(info_cols)))
                        elif len(info_cols) > n_other:
                            info_cols = info_cols[:n_other]
                        row_cols.extend(info_cols)

                    if len(row_cols) < expected_cols:
                        row_cols.extend([''] * (expected_cols - len(row_cols)))
                    elif len(row_cols) > expected_cols:
                        row_cols = row_cols[:expected_cols]

                    f.write("\t".join(row_cols) + "\n")
                    linecount += 1
                
                logger.info(f"Successfully wrote {linecount} rows to {final_out}")
        
        except Exception as e:
            logger.error(f"Error writing output file: {e}")
            import traceback
            logger.error(f"Detailed error information: {traceback.format_exc()}")

    def _cleanup_temp_files(self):
        """Clean up temporary files"""
        for temp_file in self.unlink_files:
            try:
                if os.path.exists(temp_file):
                    os.remove(temp_file)
                    logger.info(f"Deleted temporary file: {temp_file}")
            except Exception as e:
                logger.warning(f"Failed to delete temporary file {temp_file}: {e}")

    def _load_mane_transcripts(self) -> Dict[str, str]:
        """Load MANE transcript information（优先 resources/mane GTF，其次 -mane_file / humandb）。"""
        try:
            mt = _mane_mod()
            resources_dir = mt.project_resources_dir()
            # dbloc 通常为 .../resources/humandb → 上一级即 resources
            if self.dbloc:
                parent = os.path.dirname(os.path.abspath(self.dbloc))
                if os.path.basename(parent) == 'resources' or os.path.isdir(os.path.join(parent, 'mane')):
                    resources_dir = parent
            explicit = self.mane_file if getattr(self, 'mane_file', None) else None
            if not explicit:
                for name in ("mane.bb", "mane_transcript.txt"):
                    cand = os.path.join(self.dbloc, name) if self.dbloc else ""
                    if cand and os.path.isfile(cand):
                        explicit = cand
                        break
            # 若 CLI 给了不存在的路径，仍回退到 resources/mane
            if explicit and not os.path.isfile(explicit):
                logger.warning(f"指定的 MANE 文件不存在，回退自动发现: {explicit}")
                explicit = None
            resolved = mt.resolve_mane_file(explicit, resources_dir=resources_dir)
            logger.debug(f"MANE 映射文件: {resolved}")
            return mt.load_mane_transcripts(resolved, resources_dir=resources_dir)
        except Exception as e:
            logger.error(f"Failed to load MANE transcript file: {e}")
            return {}
    
    def _is_mane_transcript_match(self, transcript_part: str, mane_info: Dict) -> bool:
        """Check if a transcript part matches MANE transcript information"""
        try:
            return _mane_mod().is_mane_transcript_token(transcript_part, mane_info)
        except Exception as e:
            logger.warning(f"Failed to check MANE transcript match: {e}")
            return False

def main():
    """Main function"""
    examples = (
        "Examples:\n"
        "1) Gene annotation + region annotation based on MV input file (TSV output)\n"
        "   python utils/matchvar/table_matchvar.py \\\n+        /Users/James/PycharmProjects/Variant_Data_Simulation_2.0/resources/202511.family.mvinput \\\n+        /Users/James/PycharmProjects/Variant_Data_Simulation_2.0/resources/humandb \\\n+        -outfile /Users/James/PycharmProjects/Variant_Data_Simulation_2.0/resources/matchvar \\\n+        -buildver hg19 -protocol refGene,cytoBand -operation g,r\n\n"
        "2) Start directly from VCF or 4-column CHROM POS REF ALT (native, no convert2matchvar)\n"
        "   python -m matchvar_annotator.table_matchvar sample.vcf resources/humandb \\\n"
        "        -outfile result -buildver hg19 -protocol refGene -operation g\n\n"
        "3) Legacy path: convert VCF → .mvinput (convert2matchvar), then annotate\n"
        "   python -m matchvar_annotator.table_matchvar sample.vcf resources/humandb \\\n"
        "        -convertvcf -otherinfo -outfile result -buildver hg19 \\\n"
        "        -protocol refGene,exac03 -operation g,f\n\n"
        "4) Specify threads and NA placeholder, output to current directory\n"
        "   python utils/matchvar/table_matchvar.py input.mvinput resources/humandb \\\n+        -outfile out -thread 8 -nastring . -buildver hg19 \\\n+        -protocol refGene,clinvar -operation g,f\n\n"
        "Tip: You can specify the subprocess interpreter by setting the PYTHON_EXECUTABLE environment variable; otherwise, automatically find the project .venv or fall back to the current interpreter."
    )
    parser = argparse.ArgumentParser(
        description='MATCHVAR table annotation tool',
        epilog=examples,
        formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument('queryfile', help='Input file')
    parser.add_argument('dbloc', help='Database location')
    parser.add_argument('-outfile', help='Output file prefix')
    parser.add_argument('-buildver', default='hg19', help='Genome version')
    parser.add_argument('-remove', action='store_true', help='Remove temporary files')
    parser.add_argument('-protocol', help='Protocol list, separated by commas')
    parser.add_argument('-operation', help='Operation list, separated by commas')
    parser.add_argument('-otherinfo', action='store_true', help='Include other information')
    parser.add_argument('-nastring', help='NA string')
    parser.add_argument('-argument', help='Parameter list, separated by commas')
    parser.add_argument('-thread', type=int, help='Thread number')
    parser.add_argument('-mane_file', type=str, help='MANE transcript mapping file')
    parser.add_argument('-use_mane_transcript', action='store_true', help='Use MANE transcript filtering')
    # PLATFORM_UNUSED CLI — MatchvarRunner 不传下列参数（独立 CLI 仍可用）
    parser.add_argument('-csvout', action='store_true', help='Output CSV format')
    parser.add_argument('-vcfinput', action='store_true', help='Force native VCF CHROM/POS/REF/ALT parsing (auto-detected for .vcf and 4-column files)')
    parser.add_argument('-convertvcf', action='store_true', help='Convert VCF to .mvinput via convert2matchvar, then annotate (legacy path)')
    parser.add_argument('-dot2underline', action='store_true', help='Replace dots with underscores')
    parser.add_argument('-polishgene', action='store_true', help='Optimize gene annotation (enabled by default)')
    parser.add_argument('-nopolish', action='store_true', help='Disable coding-change polish')
    parser.add_argument('-intronhgvs', action='store_true', help='Output intronic HGVSp')
    parser.add_argument('-verbose', '-v', action='store_true', help='Verbose output')
    parser.add_argument('-man', '-m', action='store_true', help='Display manual')
    parser.add_argument('-checkfile', action='store_true', help='Check file existence')
    parser.add_argument('-onetranscript', action='store_true', help='Single transcript mode')
    parser.add_argument('-genericdbfile', type=str, help='Generic database file')
    parser.add_argument('-gff3dbfile', type=str, help='GFF3 database file')
    parser.add_argument('-bedfile', type=str, help='BED file')
    parser.add_argument('-vcfdbfile', type=str, help='VCF database file')
    parser.add_argument('-tempdir', type=str, help='Temporary directory')
    parser.add_argument('-xreffile', type=str, help='Cross-reference file')
    parser.add_argument('-convertarg', type=str, help='Conversion parameter')
    parser.add_argument('-codingarg', type=str, help='Coding parameter')
    
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(0)

    args = parser.parse_args()
    
    # Create the annotator
    annotator = TableAnnotator(
        queryfile=args.queryfile,
        dbloc=args.dbloc,
        outfile=args.outfile,
        buildver=args.buildver,
        remove=args.remove,
        protocol=args.protocol,
        operation=args.operation,
        otherinfo=args.otherinfo,
        nastring=args.nastring,
        csvout=args.csvout,
        argument=args.argument,
        vcfinput=args.vcfinput,
        convertvcf=args.convertvcf,
        dot2underline=args.dot2underline,
        thread=args.thread,
        nopolish=args.nopolish,
        intronhgvs=args.intronhgvs,
        verbose=args.verbose,
        man=args.man,
        checkfile=args.checkfile,
        onetranscript=args.onetranscript,
        genericdbfile=args.genericdbfile,
        gff3dbfile=args.gff3dbfile,
        bedfile=args.bedfile,
        vcfdbfile=args.vcfdbfile,
        tempdir=args.tempdir,
        # maxgenethread=args.maxgenethread,
        xreffile=args.xreffile,
        convertarg=args.convertarg,
        codingarg=args.codingarg,
        mane_file=args.mane_file,
        use_mane_transcript=args.use_mane_transcript
    )
    
    # Run annotation
    try:
        annotator.run_annotation()
        logger.info("Annotation completed successfully")
    except Exception as e:
        logger.error(f"Error during annotation: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main() 