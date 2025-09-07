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

    # 2) Local .venv
    project_root = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))
    venv_unix = os.path.join(project_root, '.venv', 'bin', 'python')
    venv_win = os.path.join(project_root, '.venv', 'Scripts', 'python.exe')
    if os.name == 'nt' and os.path.exists(venv_win):
        return venv_win
    if os.path.exists(venv_unix):
        return venv_unix

    # 3) Fallback to current interpreter
    return sys.executable

PYTHON_EXECUTABLE = _detect_python_executable()

# Set the encoding of standard output and error output to UTF-8
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
sys.stderr = io.TextIOWrapper(sys.stderr.buffer, encoding='utf-8', errors='replace')

# Set the log
logging.basicConfig(
    level=logging.INFO, 
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.StreamHandler(sys.stdout)
    ]
)
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
ANNOTATION_HEADERS = {
    "ljb_all": ["LJB_PhyloP", "LJB_PhyloP_Pred", "LJB_SIFT", "LJB_SIFT_Pred", "LJB_PolyPhen2", "LJB_PolyPhen2_Pred", "LJB_LRT", "LJB_LRT_Pred", "LJB_MutationTaster", "LJB_MutationTaster_Pred", "LJB_GERP++"],
    "ljb2_all": ["LJB2_SIFT", "LJB2_PolyPhen2_HDIV", "LJB2_PP2_HDIV_Pred", "LJB2_PolyPhen2_HVAR", "LJB2_PolyPhen2_HVAR_Pred", "LJB2_LRT", "LJB2_LRT_Pred", "LJB2_MutationTaster", "LJB2_MutationTaster_Pred", "LJB_MutationAssessor", "LJB_MutationAssessor_Pred", "LJB2_FATHMM", "LJB2_GERP++", "LJB2_PhyloP", "LJB2_SiPhy"],
    "popfreq_all": ["PopFreqMax", "1000G2012APR_ALL", "1000G2012APR_AFR", "1000G2012APR_AMR", "1000G2012APR_ASN", "1000G2012APR_EUR", "ESP6500si_ALL", "ESP6500si_AA", "ESP6500si_EA", "CG46"],
    # Add gene annotation fields
    "refGene": ["Func.refGene", "Gene.refGene", "GeneDetail.refGene", "ExonicFunc.refGene", "AAChange.refGene"],
    "ensGene": ["Func.ensGene", "Gene.ensGene", "GeneDetail.ensGene", "ExonicFunc.ensGene", "AAChange.ensGene"],
    "knownGene": ["Func.knownGene", "Gene.knownGene", "GeneDetail.knownGene", "ExonicFunc.knownGene", "AAChange.knownGene"],
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
        self.csvout = kwargs.get('csvout', False)
        self.argument = kwargs.get('argument')
        self.vcfinput = kwargs.get('vcfinput', False)
        self.dot2underline = kwargs.get('dot2underline', False)
        self.thread = kwargs.get('thread')
        # Default to enable protein annotation optimization based on mRNA
        self.polish = True
        self.intronhgvs = kwargs.get('intronhgvs', False)
        
        # Add new important parameters
        self.verbose = kwargs.get('verbose', False)  # Verbose output
        self.man = kwargs.get('man', False)  # Manual
        self.checkfile = kwargs.get('checkfile', False)  # File check
        self.onetranscript = kwargs.get('onetranscript', False)  # Single transcript
        self.genericdbfile = kwargs.get('genericdbfile')  # Generic database file
        self.gff3dbfile = kwargs.get('gff3dbfile')  # GFF3 database file
        self.bedfile = kwargs.get('bedfile')  # BED file
        self.vcfdbfile = kwargs.get('vcfdbfile')  # VCF database file
        self.tempdir = kwargs.get('tempdir')  # Temporary directory
        self.maxgenethread = kwargs.get('maxgenethread', 16)  # Maximum gene thread number
        self.xreffile = kwargs.get('xreffile')  # Cross-reference file
        self.convertarg = kwargs.get('convertarg')  # Conversion parameters
        self.codingarg = kwargs.get('codingarg')  # Coding parameters
        self.mane_file = kwargs.get('mane_file')  # MANE transcript mapping file
        self.use_mane_transcript = kwargs.get('use_mane_transcript', False)  # Use MANE transcript filtering
        
        # Add missing attributes that are referenced in _process_arguments
        self.filter = kwargs.get('filter', False)  # Filter operation
        self.regionanno = kwargs.get('regionanno', False)  # Region annotation
        self.geneanno = kwargs.get('geneanno', False)  # Gene annotation
        
        # Internal variables
        self.unlink_files = []
        self.header = []
        self.varanno = {}
        self.protocols = []
        self.operations = []
        self.arguments = []
        self.dbtype1 = []
        
        # Preload MANE transcript mapping (only load once)
        self.mane_transcripts = self._load_mane_transcripts()
        
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
        if self.vcfinput and self.csvout:
            raise ValueError("Error in argument: -csvout is not compatible with -vcfinput")
        
        # Verify file related parameters
        if self.genericdbfile and not (self.filter or self.regionanno):
            raise ValueError("Error in argument: the --genericdbfile argument is supported only for the --filter and --region operation")
        
        if self.gff3dbfile and not (self.geneanno or self.regionanno):
            raise ValueError("Error in argument: the --gff3dbfile argument is supported only for the --geneanno or --regionanno operation")
        
        if self.bedfile and not self.regionanno:
            raise ValueError("Error in argument: the --bedfile argument is supported only for the --regionanno operation")
        
        if self.vcfdbfile and not self.filter:
            raise ValueError("Error in argument: the --vcfdbfile argument is supported only for the --filter operation")
        
        # VCF input: force -nastring '.' and automatically enable -otherinfo
        if self.vcfinput:
            if self.nastring is not None and self.nastring != '.':
                raise ValueError("Error in argument: -nastring must be '.' when '-vcfinput' is specified")
            self.nastring = '.'
            self.otherinfo = True
        else:
            if self.nastring is None:
                self.nastring = '.'

        # Verify thread parameters
        if self.thread and self.thread > self.maxgenethread:
            logger.info(f"NOTICE: number of threads is reduced to {self.maxgenethread}")
            self.thread = self.maxgenethread
        
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
            db_file = os.path.join(self.dbloc, f"{self.buildver}_{dbtype1}.txt")
            if not os.path.exists(db_file):
                missing_files.append(db_file)
                logger.warning(f"Database file does not exist: {db_file}")
        
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
        
        # Process VCF input
        if self.vcfinput:
            self._handle_vcf_input()
        else:
            self._run_standard_annotation()
        
        # Print original output
        self._print_original_output()
        
        # Clean up temporary files
        if self.remove:
            self._cleanup_temp_files()
    
    def _handle_vcf_input(self):
        """Process VCF input"""
        if self.csvout:
            raise ValueError("Error: -csvout is not compatible with -vcfinput")
        
        # Get the absolute path of the current script
        current_dir = os.path.dirname(os.path.abspath(__file__))
        # Get the path of the convert2matchvar.py script
        convert2matchvar_script = os.path.join(current_dir, 'convert2matchvar.py')
        
        # Convert VCF to MATCHVAR input format
        convertarg_str = f"{self.convertarg} " if self.convertarg else ""
        if self.queryfile.endswith('.vcf') or self.queryfile.endswith('.vcf.gz'):
            sc = f"{PYTHON_EXECUTABLE} {convert2matchvar_script} {convertarg_str}-includeinfo -allsample -withfreq -format vcf4 {self.queryfile} > {self.tempfile}.mvinput"
        else:
            sc = f"{PYTHON_EXECUTABLE} {convert2matchvar_script} {convertarg_str}-includeinfo -allsample -withfreq -format vcf4 {self.queryfile} > {self.tempfile}.mvinput"
        logger.info(f"NOTICE: Running with system command <{sc}>")
        
        result = run_subprocess_safe(sc)
        if result.returncode != 0:
            raise RuntimeError(f"Error running system command: <{sc}>")
        
        # Use standard gene annotation process (supports -dbtype, -exonsort, and polishgene)
        # Switch queryfile to the converted .mvinput, then reuse _run_standard_annotation
        self.queryfile = f"{self.tempfile}.mvinput"
        self._run_standard_annotation()
    
    def _run_standard_annotation(self):
        """Run standard annotation"""
        for i, (protocol, operation) in enumerate(zip(self.protocols, self.operations)):
            logger.info("-----------------------------------------------------------------")
            logger.info(f"NOTICE: Processing operation={operation} protocol={protocol}")
            
            if operation in ['g', 'gx']:
                self._gene_operation(protocol, self.dbtype1[i], self.arguments[i] if i < len(self.arguments) else None, operation)
            elif operation == 'r':
                self._region_operation(protocol, self.dbtype1[i], self.arguments[i] if i < len(self.arguments) else None)
            elif operation == 'f':
                self._filter_operation(protocol, self.dbtype1[i], self.arguments[i] if i < len(self.arguments) else None)
    
    def _gene_operation(self, protocol: str, dbtype1: str, argument: str, operation: str):
        """Gene annotation operation"""
        # Process protocol names
        genetype = {'gene': 'refGene', 'refgene': 'refGene', 'knowngene': 'knownGene', 'ensgene': 'ensGene'}
        if protocol in genetype:
            protocol = genetype[protocol]
        
        # Get the absolute path of the current script
        current_dir = os.path.dirname(os.path.abspath(__file__))
        annotate_variation_script = os.path.join(current_dir, 'annotate_variation.py')
        
        # Build command, using the Python interpreter in the virtual environment
        sc = f"{PYTHON_EXECUTABLE} {annotate_variation_script} -geneanno -buildver {self.buildver} -dbtype {protocol} -outfile {self.tempfile}.{protocol} -exonsort -nofirstcodondel {self.queryfile} {self.dbloc}"
        
        # Add splicing_threshold parameter
        if self.intronhgvs:
            sc += f" -splicing_threshold {self.intronhgvs}"
        
        # Add MANE transcript mapping parameter
        if self.mane_transcripts:
            # Write MANE transcript mapping to temporary file
            mane_file = f"{self.tempfile}.{protocol}.mane"
            try:
                with open(mane_file, 'w', encoding='utf-8') as f:
                    for gene_id, transcript_id in self.mane_transcripts.items():
                        f.write(f"{gene_id}\t{transcript_id}\n")
                sc += f" -mane_file {mane_file}"
                self.unlink_files.append(mane_file)
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
                # Prepare the required file paths
                gene_file = os.path.join(self.dbloc, f"{self.buildver}_{protocol}.txt")
                mrna_fa = os.path.join(self.dbloc, f"{self.buildver}_{protocol}Mrna.fa")
                coding_change_script_py = os.path.join(current_dir, 'coding_change.py')

                # Only execute polish step when the key files exist
                if not os.path.exists(gene_file):
                    logger.warning(f"Gene definition file does not exist, skipping polish step: {gene_file}")
                elif not os.path.exists(mrna_fa):
                    logger.warning(f"mRNA FASTA file does not exist, skipping polish step: {mrna_fa}")
                else:
                    # Rename the original exonic_variant_function to .orig
                    e_anno_outfile_orig = f"{e_anno_outfile}.orig"
                    try:
                        os.rename(e_anno_outfile, e_anno_outfile_orig)
                    except Exception as e:
                        logger.warning(f"Failed to rename {e_anno_outfile} -> {e_anno_outfile_orig}, skipping polish step: {e}")
                        e_anno_outfile_orig = None

                    if e_anno_outfile_orig:
                        # Only use the coding_change, using the current interpreter
                        sc_cc = (
                            f"{PYTHON_EXECUTABLE} {coding_change_script_py} "
                            f"{self.codingarg + ' ' if self.codingarg else ''}"
                            f"{e_anno_outfile_orig} {gene_file} {mrna_fa} "
                            f"-includesnp -alltranscript -out {self.tempfile}.{protocol}.fa -newevf {e_anno_outfile}"
                        )
                        logger.info(f"NOTICE: Running with system command <{sc_cc}>")
                        # Run with blocking and capturing output, if failed, output stderr, avoid silent
                        result_cc = run_subprocess_safe(sc_cc, real_time_output=False)
                        if result_cc.returncode != 0:
                            logger.error(f"coding_change failed, return code: {result_cc.returncode}")
                            logger.error(f"stdout: {getattr(result_cc, 'stdout', '')}")
                            logger.error(f"stderr: {getattr(result_cc, 'stderr', '')}")
                        if result_cc.returncode != 0:
                            raise RuntimeError(f"Error running system command: <{sc_cc}>")

                        # Record the temporary files to be cleaned up
                        self.unlink_files.append(f"{self.tempfile}.{protocol}.fa")
                        self.unlink_files.append(e_anno_outfile_orig)
            except Exception as e:
                logger.warning(f"Polish step failed, using unmodified exonic annotation: {e}")
        
        # Set the header
        if self.dot2underline:
            # Add VarType column after ExonicFunc
            self.header.extend([f"Func_{protocol}", f"Gene_{protocol}", f"GeneDetail_{protocol}", f"ExonicFunc_{protocol}", f"VarType_{protocol}", f"AAChange_{protocol}"])
        else:
            self.header.extend([f"Func.{protocol}", f"Gene.{protocol}", f"GeneDetail.{protocol}", f"ExonicFunc.{protocol}", f"VarType.{protocol}", f"AAChange.{protocol}"])
        
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
                        gene_name = gene
                        aa_change = 'p.?'
                        
                        # Extract the transcript information in parentheses
                        if '(' in gene and ')' in gene:
                            import re
                            # Extract all the content in parentheses
                            transcript_matches = re.findall(r'\(([^)]+)\)', gene)
                            if transcript_matches:
                                # Use the transcript information as GeneDetail
                                gene_detail = ';'.join(transcript_matches)
                                # Remove the content in parentheses, keep the gene name
                                gene_name = re.sub(r'\([^)]+\)', '', gene).strip()
                        
                        # Process GeneDetail and AAChange for intronic variants
                        if function == 'intronic' and gene_detail:
                            # For intronic variants, gene_detail contains the full annotation
                            # Add :p.? if not present
                            if ':p.' not in gene_detail:
                                gene_detail = gene_detail + ':p.?'
                            # AAChange should be p.? for intronic variants
                            aa_change = 'p.?'
                        
                        # Store the annotation - ensure not to overwrite the previous annotation
                        if varstring not in self.varanno:
                            self.varanno[varstring] = {}
                        
                        if self.dot2underline:
                            self.varanno[varstring].update({
                                f"Func_{protocol}": function,
                                f"Gene_{protocol}": gene_name,
                                f"GeneDetail_{protocol}": gene_detail,
                                f"ExonicFunc_{protocol}": 'NA',
                                f"VarType_{protocol}": 'NA',
                                f"AAChange_{protocol}": aa_change
                            })
                        else:
                            self.varanno[varstring].update({
                                f"Func.{protocol}": function,
                                f"Gene.{protocol}": gene_name,
                                f"GeneDetail.{protocol}": gene_detail,
                                f"ExonicFunc.{protocol}": 'NA',
                                f"VarType.{protocol}": 'NA',
                                f"AAChange.{protocol}": aa_change
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
                            # Write ExonicFunc and AAChange
                            if self.dot2underline:
                                self.varanno[varstring][f"ExonicFunc_{protocol}"] = exonic_function
                                self.varanno[varstring][f"AAChange_{protocol}"] = aa_change
                            else:
                                self.varanno[varstring][f"ExonicFunc.{protocol}"] = exonic_function
                                self.varanno[varstring][f"AAChange.{protocol}"] = aa_change

                            # Calculate VarType (uniform standard: non-frameshift = inframe)
                            try:
                                func_key = f"Func_{protocol}" if self.dot2underline else f"Func.{protocol}"
                                vartype_key = f"VarType_{protocol}" if self.dot2underline else f"VarType.{protocol}"
                                func_val = self.varanno[varstring].get(func_key, '').lower()
                                exonic_val = exonic_function.lower() if exonic_function else ''

                                vartype = 'unknown'
                                # 1) Highest priority: splicing/terminator
                                if func_val.find('splicing') != -1 or exonic_val.find('splicing') != -1:
                                    vartype = 'splicing'
                                elif 'stopgain' in exonic_val:
                                    vartype = 'nonsense'
                                elif 'stoploss' in exonic_val:
                                    vartype = 'stoploss'
                                else:
                                    # 2) First, classify by effect text
                                    if 'frameshift' in exonic_val:
                                        vartype = 'frameshift'
                                    elif (
                                        'nonframeshift' in exonic_val
                                        or 'inframe' in exonic_val
                                        or ('duplication' in exonic_val and 'frameshift' not in exonic_val)
                                    ):
                                        vartype = 'inframe'
                                    elif 'nonsynonymous' in exonic_val or 'missense' in exonic_val:
                                        vartype = 'missense'
                                    elif 'synonymous' in exonic_val:
                                        vartype = 'synonymous'

                                    # 3) Cross-validation (strong validation):
                                    #    a) If p.HGVS contains fs*, it is frameshift; otherwise continue
                                    #    b) Determine if the change in Ref/Alt nucleotides (|len(alt)-len(ref)|) is a multiple of 3
                                    try:
                                        toks = varstring.split('\t')
                                        if len(toks) >= 5:
                                            ref_nt = toks[3].replace('-', '').replace('*', '')
                                            alt_nt = toks[4].replace('-', '').replace('*', '')
                                            # a) Protein layer priority: fs*
                                            if aa_change and 'p.' in aa_change and 'fs*' in aa_change:
                                                vartype = 'frameshift'
                                            else:
                                                # b) Change in nucleotides
                                                delta = abs(len(alt_nt) - len(ref_nt))
                                                if delta == 0:
                                                    # Equal length and length>1: inframe replacement (if not defined by missense/synonymous)
                                                    if len(ref_nt) > 1 and vartype not in ('missense', 'synonymous'):
                                                        vartype = 'inframe'
                                                else:
                                                    vartype = 'inframe' if (delta % 3 == 0) else 'frameshift'
                                    except Exception:
                                        pass

                                    # 4) Backtracking for dup: If the text does not explicitly indicate the frame shift property, parse the length of the dup sequence
                                    if 'duplication' in exonic_val and ('nonframeshift' not in exonic_val and 'frameshift' not in exonic_val):
                                        try:
                                            mdup = re.search(r'c\.[^\s:]*dup([ACGTN]+)', aa_change or '', re.IGNORECASE)
                                            if mdup:
                                                dup_nt = mdup.group(1).upper()
                                                vartype = 'inframe' if (len(dup_nt) % 3 == 0) else 'frameshift'
                                        except Exception:
                                            pass

                                self.varanno[varstring][vartype_key] = vartype
                            except Exception as _:
                                # Safe downgrade: keep NA
                                pass

                            # Process GeneDetail and AAChange separately
                            try:
                                gene_detail_key = f"GeneDetail_{protocol}" if self.dot2underline else f"GeneDetail.{protocol}"
                                aa_change_key = f"AAChange_{protocol}" if self.dot2underline else f"AAChange.{protocol}"
                                
                                detail_candidate = aa_change or ''
                                if detail_candidate and ':' in detail_candidate:
                                    # Normalize c.hgvs fragment
                                    try:
                                        import re
                                        # Extract and normalize c.hgvs
                                        m = re.search(r'(c\.[^:\s]+)', detail_candidate)
                                        if m:
                                            from .coding_change import normalize_c_hgvs  # Reuse the same normalization method
                                            norm_c = normalize_c_hgvs(m.group(1))
                                            detail_candidate = detail_candidate.replace(m.group(1), norm_c)
                                    except Exception:
                                        pass
                                    
                                    # Process GeneDetail: keep full format
                                    gene_detail_final = detail_candidate
                                    # If c. exists but lacks :p., add a placeholder p.?
                                    if ('c.' in gene_detail_final) and (':p.' not in gene_detail_final):
                                        gene_detail_final = gene_detail_final + ':p.?'
                                    # Only when c. or :p. is present is it considered a complete detail
                                    if ('c.' in gene_detail_final) or (':p.' in gene_detail_final):
                                        self.varanno[varstring][gene_detail_key] = gene_detail_final
                                    
                                    # Process AAChange: extract only p.HGVS part
                                    aa_change_final = 'p.?'
                                    if ':p.' in detail_candidate:
                                        # Extract p.HGVS part
                                        p_match = re.search(r':p\.([^:\s]+)', detail_candidate)
                                        if p_match:
                                            aa_change_final = f"p.{p_match.group(1)}"
                                    self.varanno[varstring][aa_change_key] = aa_change_final
                            except Exception:
                                pass
        
        except Exception as e:
            logger.error(f"Error reading exonic annotation file: {e}")
    
    def _region_operation(self, protocol: str, dbtype1: str, argument: str):
        """Region annotation operation"""
        # Get the absolute path of the current script
        current_dir = os.path.dirname(os.path.abspath(__file__))
        annotate_variation_script = os.path.join(current_dir, 'annotate_variation.py')
        
        sc = f"{PYTHON_EXECUTABLE} {annotate_variation_script} -regionanno -dbtype {protocol} -buildver {self.buildver} -outfile {self.tempfile} {self.queryfile} {self.dbloc}"
        
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
        
        sc = f"{PYTHON_EXECUTABLE} {annotate_variation_script} -filter -dbtype {protocol} -buildver {self.buildver} -outfile {self.tempfile} {self.queryfile} {self.dbloc}"
        
        if argument and argument.strip():
            sc += f" {argument}"
        
        if self.thread:
            sc += f" -thread {self.thread}"
        
        # Check for both uncompressed and compressed versions
        dbfile_txt = os.path.join(self.dbloc, f"{self.buildver}_{protocol}.txt")
        dbfile_gz = dbfile_txt + ".gz"
        tbi_file = dbfile_gz + ".tbi"
        
        # Use compressed version if available with index, otherwise use uncompressed
        dbfile = None
        if os.path.exists(dbfile_gz) and os.path.exists(tbi_file):
            dbfile = dbfile_gz
            logger.info(f"Using compressed database with Tabix index: {dbfile}")
        elif os.path.exists(dbfile_txt):
            dbfile = dbfile_txt
            logger.info(f"Using uncompressed database: {dbfile}")
        elif os.path.exists(dbfile_gz):
            dbfile = dbfile_gz
            logger.info(f"Using compressed database (no index): {dbfile}")
        
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
            with open(final_out, 'w', encoding='utf-8') as f:
                linecount = 0
                with open(self.queryfile, 'r', encoding='utf-8') as input_f:
                    for line in input_f:
                        line = line.strip()
                        if not line or line.startswith('#'):
                            continue
                        
                        parts = line.split('\t')
                        if len(parts) >= 5:
                            varstring = '\t'.join(parts[:5])
                            info = '\t'.join(parts[5:]) if len(parts) > 5 else ''
                            
                            logger.debug(f"Processing line: {varstring}")
                            
                            if linecount == 0:
                                # Write the header (fixed as TSV)
                                header_line = "\t".join(['Chr', 'Start', 'End', 'Ref', 'Alt'] + expanded_header)
                                if self.otherinfo:
                                    num_info = len(parts) - 5
                                    for i in range(1, num_info + 1):
                                        header_line += f"\tOtherinfo{i}"
                                f.write(header_line + "\n")
                            
                            # Write the data line
                            oneline = []
                            for item in self.header:
                                if item in ANNOTATION_HEADERS:
                                    # Process the extended field
                                    expanded_field = len(ANNOTATION_HEADERS[item])
                                    if varstring in self.varanno and item in self.varanno[varstring]:
                                        # Uniformly split the annotation values by tab
                                        values = str(self.varanno[varstring][item]).split('\t')
                                        # Replace escape characters and standardize empty values
                                        values = [v.replace('\\x2c', ',').replace('\\x23', '#') for v in values]
                                        # Standardize empty values: convert -1, empty strings, and other common empty representations to "."
                                        values = ['.' if v in ['-1', '', 'NA', 'N/A', 'null', 'NULL'] else v for v in values]
                                        # Align by expected number of columns: if less, fill with NA, if more, truncate
                                        if len(values) < expanded_field:
                                            values.extend([str(self.nastring)] * (expanded_field - len(values)))
                                        elif len(values) > expanded_field:
                                            values = values[:expanded_field]
                                        oneline.extend(values)
                                    else:
                                        for _ in range(expanded_field):
                                            oneline.append(str(self.nastring))
                                else:
                                    # Process the normal field
                                    if varstring in self.varanno and item in self.varanno[varstring]:
                                        oneline.append(str(self.varanno[varstring][item]))
                                    else:
                                        oneline.append(str(self.nastring))

                            # Fixed as TSV output, ensure the first 5 columns are always written
                            output_line = "\t".join(parts[:5] + oneline)
                            if self.otherinfo:
                                output_line += f"\t{info}"
                            f.write(output_line + "\n")

                            # Verify that the number of columns in the output line matches the header
                            output_cols = output_line.split('\t')
                            if linecount == 0:
                                expected_cols = len(output_cols)
                            elif len(output_cols) != expected_cols:
                                logger.warning(
                                    f"Row {linecount + 1} has inconsistent number of columns: expected {expected_cols} columns, actual {len(output_cols)} columns")
                                # Fill or truncate to the correct number of columns
                                if len(output_cols) < expected_cols:
                                    output_cols.extend([''] * (expected_cols - len(output_cols)))
                                else:
                                    output_cols = output_cols[:expected_cols]
                                output_line = '\t'.join(output_cols)
                                f.seek(f.tell() - len(output_line + '\n'))
                                f.write(output_line + '\n')
                            
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
        """Load MANE transcript information"""
        mane_transcripts = {}
        mane_file = os.path.join(self.dbloc, 'mane_transcript.txt')
        
        if not os.path.exists(mane_file):
            logger.warning(f"MANE transcript file does not exist: {mane_file}")
            return mane_transcripts
        
        try:
            with open(mane_file, 'r', encoding='utf-8') as f:
                for line in f:
                    line = line.strip()
                    if not line or line.startswith('#'):
                        continue
                    
                    parts = line.split('\t')
                    if len(parts) >= 9:
                        # GTF format: chr, source, feature, start, end, score, strand, frame, attributes
                        attributes = parts[8]
                        
                        # Parse the attribute field
                        gene_id = None
                        transcript_id = None
                        
                        for attr in attributes.split(';'):
                            attr = attr.strip()
                            if attr.startswith('gene_id'):
                                gene_id = attr.split('"')[1] if '"' in attr else attr.split()[1]
                            elif attr.startswith('transcript_id'):
                                transcript_id = attr.split('"')[1] if '"' in attr else attr.split()[1]
                        
                        if gene_id and transcript_id:
                            mane_transcripts[gene_id] = transcript_id
            
            logger.info(f"Loaded {len(mane_transcripts)} MANE transcript mappings")
        except Exception as e:
            logger.error(f"Failed to load MANE transcript file: {e}")
        
        return mane_transcripts

def main():
    """Main function"""
    examples = (
        "Examples:\n"
        "1) Gene annotation + region annotation based on MV input file (TSV output)\n"
        "   python utils/matchvar/table_matchvar.py \\\n+        /Users/James/PycharmProjects/Variant_Data_Simulation_2.0/resources/202511.family.mvinput \\\n+        /Users/James/PycharmProjects/Variant_Data_Simulation_2.0/resources/humandb \\\n+        -outfile /Users/James/PycharmProjects/Variant_Data_Simulation_2.0/resources/matchvar \\\n+        -buildver hg19 -protocol refGene,cytoBand -operation g,r\n\n"
        "2) Start directly from VCF (automatically converts to MV input internally) and retain original info columns\n"
        "   python utils/matchvar/table_matchvar.py \\\n+        /Users/James/PycharmProjects/Variant_Data_Simulation_2.0/resources/202511.family.vcf \\\n+        /Users/James/PycharmProjects/Variant_Data_Simulation_2.0/resources/humandb \\\n+        -vcfinput -otherinfo -outfile result -buildver hg19 \\\n+        -protocol refGene,exac03,avsift -operation g,f,f\n\n"
        "3) Specify threads and NA placeholder, output to current directory\n"
        "   python utils/matchvar/table_matchvar.py input.mvinput resources/humandb \\\n+        -outfile out -thread 8 -nastring . -buildver hg19 \\\n+        -protocol refGene,clinvar -operation g,f\n\n"
        "Tip: You can specify the subprocess interpreter by setting the environment variable PYTHON_EXECUTABLE; otherwise, it will automatically look for the project .venv or fall back to the current interpreter."
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
    parser.add_argument('-csvout', action='store_true', help='Output CSV format')
    parser.add_argument('-argument', help='Parameter list, separated by commas')
    parser.add_argument('-vcfinput', action='store_true', help='VCF input')
    parser.add_argument('-dot2underline', action='store_true', help='Replace dots with underscores')
    parser.add_argument('-thread', type=int, help='Thread number')
    parser.add_argument('-polishgene', action='store_true', help='Optimize gene annotation')
    parser.add_argument('-intronhgvs', action='store_true', help='Output intronic HGVSp')
    
    # Add new important parameters
    parser.add_argument('-verbose', '-v', action='store_true', help='Verbose output')
    parser.add_argument('-man', '-m', action='store_true', help='Display manual')
    parser.add_argument('-checkfile', action='store_true', help='Check file existence')
    parser.add_argument('-onetranscript', action='store_true', help='Single transcript mode')
    parser.add_argument('-genericdbfile', type=str, help='Generic database file')
    parser.add_argument('-gff3dbfile', type=str, help='GFF3 database file')
    parser.add_argument('-bedfile', type=str, help='BED file')
    parser.add_argument('-vcfdbfile', type=str, help='VCF database file')
    parser.add_argument('-tempdir', type=str, help='Temporary directory')
    parser.add_argument('-maxgenethread', type=int, default=16, help='Maximum gene thread number')
    parser.add_argument('-xreffile', type=str, help='Cross-reference file')
    parser.add_argument('-convertarg', type=str, help='Conversion parameter')
    parser.add_argument('-codingarg', type=str, help='Coding parameter')
    parser.add_argument('-mane_file', type=str, help='MANE transcript mapping file')
    parser.add_argument('-use_mane_transcript', action='store_true', help='Use MANE transcript filtering')
    
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
        dot2underline=args.dot2underline,
        thread=args.thread,
        polishgene=args.polishgene,
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
        maxgenethread=args.maxgenethread,
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