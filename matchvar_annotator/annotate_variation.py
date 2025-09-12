#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
MATCHVAR variation annotation tool
"""

import os
import sys
import argparse
import logging
import re
import threading
import io
from typing import Dict, List, Tuple, Optional
from datetime import datetime
from urllib import request, error
try:
    import pysam  # Optional: used for Tabix-accelerated access to huge databases
except Exception:
    pysam = None

# Set the encoding of standard output and error output to UTF-8
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
sys.stderr = io.TextIOWrapper(sys.stderr.buffer, encoding='utf-8', errors='replace')

# Set log
logging.basicConfig(
    level=logging.INFO, 
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.StreamHandler(sys.stdout)
    ]
)
logger = logging.getLogger(__name__)

# Codon table
CODON_TABLE = {
    'TTT': 'F', 'TTC': 'F', 'TCT': 'S', 'TCC': 'S', 'TAT': 'Y', 'TAC': 'Y',
    'TGT': 'C', 'TGC': 'C', 'TTA': 'L', 'TCA': 'S', 'TAA': '*', 'TGA': '*',
    'TTG': 'L', 'TCG': 'S', 'TAG': '*', 'TGG': 'W', 'CTT': 'L', 'CTC': 'L',
    'CCT': 'P', 'CCC': 'P', 'CAT': 'H', 'CAC': 'H', 'CGT': 'R', 'CGC': 'R',
    'CTA': 'L', 'CTG': 'L', 'CCA': 'P', 'CCG': 'P', 'CAA': 'Q', 'CAG': 'Q',
    'CGA': 'R', 'CGG': 'R', 'ATT': 'I', 'ATC': 'I', 'ACT': 'T', 'ACC': 'T',
    'AAT': 'N', 'AAC': 'N', 'AGT': 'S', 'AGC': 'S', 'ATA': 'I', 'ACA': 'T',
    'AAA': 'K', 'AGA': 'R', 'ATG': 'M', 'ACG': 'T', 'AAG': 'K', 'AGG': 'R',
    'GTT': 'V', 'GTC': 'V', 'GCT': 'A', 'GCC': 'A', 'GAT': 'D', 'GAC': 'D',
    'GGT': 'G', 'GGC': 'G', 'GTA': 'V', 'GTG': 'V', 'GCA': 'A', 'GCG': 'A',
    'GAA': 'E', 'GAG': 'E', 'GGA': 'G', 'GGG': 'G'
}

# Mitochondrial codon table
CODON_TABLE_MT = {
    'TTT': 'F', 'TTC': 'F', 'TCT': 'S', 'TCC': 'S', 'TAT': 'Y', 'TAC': 'Y',
    'TGT': 'C', 'TGC': 'C', 'TTA': 'L', 'TCA': 'S', 'TAA': '*', 'TGA': 'W',
    'TTG': 'L', 'TCG': 'S', 'TAG': '*', 'TGG': 'W', 'CTT': 'L', 'CTC': 'L',
    'CCT': 'P', 'CCC': 'P', 'CAT': 'H', 'CAC': 'H', 'CGT': 'R', 'CGC': 'R',
    'CTA': 'L', 'CTG': 'L', 'CCA': 'P', 'CCG': 'P', 'CAA': 'Q', 'CAG': 'Q',
    'CGA': 'R', 'CGG': 'R', 'ATT': 'I', 'ATC': 'I', 'ACT': 'T', 'ACC': 'T',
    'AAT': 'N', 'AAC': 'N', 'AGT': 'S', 'AGC': 'S', 'ATA': 'M', 'ACA': 'T',
    'AAA': 'K', 'AGA': '*', 'ATG': 'M', 'ACG': 'T', 'AAG': 'K', 'AGG': '*',
    'GTT': 'V', 'GTC': 'V', 'GCT': 'A', 'GCC': 'A', 'GAT': 'D', 'GAC': 'D',
    'GGT': 'G', 'GGC': 'G', 'GTA': 'V', 'GTG': 'V', 'GCA': 'A', 'GCG': 'A',
    'GAA': 'E', 'GAG': 'E', 'GGA': 'G', 'GGG': 'G'
}

# IUPAC codon table
IUPAC = {
    'R': 'AG', 'Y': 'CT', 'S': 'GC', 'W': 'AT', 'K': 'GT', 'M': 'AC',
    'A': 'AA', 'C': 'CC', 'G': 'GG', 'T': 'TT', 'B': 'CGT', 'D': 'AGT',
    'H': 'ACT', 'V': 'ACG', 'N': 'ACGT', '.': '-', '-': '-'
}

class AnnotateVariation:
    """Variation annotator"""
    
    def __init__(self, queryfile: str, dbloc: str, **kwargs):
        self.queryfile = queryfile
        self.dbloc = dbloc
        self.outfile = kwargs.get('outfile')
        self.dbtype = kwargs.get('dbtype')
        self.geneanno = kwargs.get('geneanno', False)
        self.regionanno = kwargs.get('regionanno', False)
        self.filter = kwargs.get('filter', False)
        self.buildver = kwargs.get('buildver', 'hg19')
        self.thread = kwargs.get('thread')
        self.maxgenethread = kwargs.get('maxgenethread', 4)
        self.mingenelinecount = kwargs.get('mingenelinecount', 1000000)
        self.splicing_threshold = kwargs.get('splicing_threshold', 2)  # Add splicing_threshold parameter
        self.indel_splicing_threshold = kwargs.get('indel_splicing_threshold')  # Add indel_splicing_threshold parameter
        self.otherinfo = kwargs.get('otherinfo', False)  # Add otherinfo parameter
        self.sift_threshold = kwargs.get('sift_threshold', 0.05)  # Add sift_threshold parameter
        self.score_threshold = kwargs.get('score_threshold')  # Add score_threshold parameter
        self.reverse = kwargs.get('reverse', False)  # Add reverse parameter
        self.rawscore = kwargs.get('rawscore', False)  # Add rawscore parameter
        # Intronic dup recognition (remote Ensembl REST)
        self.intronic_dup_remote = kwargs.get('intronic_dup_remote', False)
        # MANE transcript filtering
        self.use_mane_transcript = kwargs.get('use_mane_transcript', False)
        self.intronic_dup_window = int(kwargs.get('intronic_dup_window', 50) or 50)
        
        # Add new important parameters
        self.verbose = kwargs.get('verbose', False)  # Detailed output
        self.batchsize = kwargs.get('batchsize', 5000000)  # Batch size
        self.genomebinsize = kwargs.get('genomebinsize')  # Genome bin size
        self.neargene = kwargs.get('neargene', 1000)  # Gene distance threshold
        self.expandbin = kwargs.get('expandbin')  # Expand bin number
        self.maf_threshold = kwargs.get('maf_threshold', 0)  # MAF threshold
        self.normscore_threshold = kwargs.get('normscore_threshold')  # Normalized score threshold
        self.indexfilter_threshold = kwargs.get('indexfilter_threshold', 0.9)  # Index filter threshold
        self.chromosome = kwargs.get('chromosome')  # Chromosome limit
        self.gff3dbfile = kwargs.get('gff3dbfile')  # GFF3 database file
        self.genericdbfile = kwargs.get('genericdbfile')  # Generic database file
        self.vcfdbfile = kwargs.get('vcfdbfile')  # VCF database file
        self.precedence = kwargs.get('precedence')  # Precedence
        self.hgvs = kwargs.get('hgvs', False)  # HGVS format
        self.exonicsplicing = kwargs.get('exonicsplicing', False)  # Exonic splicing
        self.separate = kwargs.get('separate', False)  # Separate output
        self.downdb = kwargs.get('downdb', False)  # Download database
        self.time = kwargs.get('time', False)  # Time information
        self.wget = kwargs.get('wget', True)  # Use wget
        self.comment = kwargs.get('comment', False)  # Annotation information
        self.transcript_function = kwargs.get('transcript_function', False)  # Transcript function
        self.avcolumn = kwargs.get('avcolumn')  # AV column
        self.bedfile = kwargs.get('bedfile')  # BED file
        self.seq_padding = kwargs.get('seq_padding')  # Sequence padding
        self.infoasscore = kwargs.get('infoasscore', False)  # Information as score
        self.firstcodondel = kwargs.get('firstcodondel', True)  # Delete first codon
        self.aamatrixfile = kwargs.get('aamatrixfile')  # Amino acid matrix file
        self.gff3attribute = kwargs.get('gff3attribute', False)  # GFF3 attribute
        self.infosep = kwargs.get('infosep', False)  # Information separator
        self.dbm = kwargs.get('dbm', False)  # DBM mode
        self.idasscore = kwargs.get('idasscore', False)  # ID as score
        self.minqueryfrac = kwargs.get('minqueryfrac', 0)  # Minimum query fraction
        self.scorecolumn = kwargs.get('scorecolumn')  # Score column
        self.poscolumn = kwargs.get('poscolumn')  # Position column
        self.webfrom = kwargs.get('webfrom')  # Web source
        self.colsWanted = kwargs.get('colsWanted')  # Wanted columns
        self.zerostart = kwargs.get('zerostart', False)  # Zero start
        self.memfree = kwargs.get('memfree')  # Available memory
        self.memtotal = kwargs.get('memtotal')  # Total memory
        self.mane_file = kwargs.get('mane_file')  # MANE transcript mapping file
        # Record the line number of gene annotation output, used to generate exonic_variant_function compatible with Perl (starting with lineN)
        self._gene_exonic_line_counter = 0
        
        # Internal variables
        self.valichr = {}
        self.dbtype1 = None
        self._sift_threshold_explicitly_set = kwargs.get('_sift_threshold_explicitly_set', False) # Add flag to track whether sift_threshold is set through command line parameters
        self._indexfilter_threshold_explicitly_set = kwargs.get('_indexfilter_threshold_explicitly_set', False) # Add flag to track whether indexfilter_threshold is set through command line parameters
        self._maf_threshold_explicitly_set = kwargs.get('_maf_threshold_explicitly_set', False) # Add flag to track whether maf_threshold is set through command line parameters
        self._minqueryfrac_explicitly_set = kwargs.get('_minqueryfrac_explicitly_set', False) # Add flag to track whether minqueryfrac is set through command line parameters
        self._wget_explicitly_set = kwargs.get('_wget_explicitly_set', False) # Add flag to track whether wget is set through command line parameters
        self._precedence_explicitly_set = kwargs.get('_precedence_explicitly_set', False) # Add flag to track whether precedence is set through command line parameters
        
        # Preload MANE transcript mapping (only load when use_mane_transcript is True)
        if self.use_mane_transcript:
            if self.mane_file and os.path.exists(self.mane_file):
                # Load MANE transcript mapping from file
                self.mane_transcripts = self._load_mane_transcripts_from_file(self.mane_file)
            else:
                # Load MANE transcript mapping from default location
                self.mane_transcripts = self._load_mane_transcripts()
        else:
            self.mane_transcripts = {}
        
        # Process arguments
        self._process_arguments()
    
    def _process_arguments(self):
        """Process arguments"""
        # Set default values
        if not self.outfile:
            self.outfile = self.queryfile
        
        # Validate arguments
        if not any([self.geneanno, self.regionanno, self.filter]):
            raise ValueError("Error: must specify one of --geneanno, --regionanno, or --filter")
        
        if self.geneanno and self.regionanno:
            raise ValueError("Error: cannot specify both --geneanno and --regionanno")
        
        if self.geneanno and self.filter:
            raise ValueError("Error: cannot specify both --geneanno and --filter")
        
        if self.regionanno and self.filter:
            raise ValueError("Error: cannot specify both --regionanno and --filter")
        
        # Process database type
        if self.dbtype:
            self.dbtype1 = self._process_dbtype(self.dbtype)
        
        # Validate sift_threshold parameter
        if hasattr(self, 'sift_threshold') and self.sift_threshold is not None:
            # Only validate when sift_threshold is explicitly set
            if hasattr(self, '_sift_threshold_explicitly_set') and self._sift_threshold_explicitly_set:
                if not self.filter:
                    raise ValueError("Error in argument: the --sift_threshold is supported only for the --filter operation")
                if self.dbtype1 != 'avsift':
                    raise ValueError("Error in argument: the --sift_threshold argument can be used only if '--dbtype avsift' is used")
                if not (0 <= self.sift_threshold <= 1):
                    raise ValueError("Error in argument: the --sift_threshold must be between 0 and 1 inclusive")
        else:
            if self.dbtype1 == 'avsift':
                logger.info("NOTICE: The --sift_threshold is set as 0.05 by default")
                self.sift_threshold = 0.05
        
        # Validate score_threshold parameter
        if hasattr(self, 'score_threshold') and self.score_threshold is not None:
            if self.geneanno:
                raise ValueError("Error in argument: the --score_threshold is not useful for --geneanno operations")
        
        # Validate normscore_threshold parameter
        if hasattr(self, 'normscore_threshold') and self.normscore_threshold is not None:
            if not self.regionanno:
                raise ValueError("Error in argument: the --normscore_threshold is supported only for the --regionanno operation")
            if not (0 <= self.normscore_threshold <= 1000):
                raise ValueError("Error in argument: the --normscore_threshold must be between 0 and 1000")
        
        # Validate indexfilter_threshold parameter
        if hasattr(self, 'indexfilter_threshold') and self.indexfilter_threshold is not None:
            # Only validate when indexfilter_threshold is explicitly set
            if hasattr(self, '_indexfilter_threshold_explicitly_set') and self._indexfilter_threshold_explicitly_set:
                if not self.filter:
                    raise ValueError("Error in argument: the --indexfilter_threshold is supported only for the --filter operation")
                if not (0 <= self.indexfilter_threshold <= 1):
                    raise ValueError("Error in argument: the --indexfilter_threshold must be between 0 and 1 inclusive")
        
        # Validate maf_threshold parameter
        if hasattr(self, 'maf_threshold') and self.maf_threshold is not None:
            # Only validate when maf_threshold is explicitly set
            if hasattr(self, '_maf_threshold_explicitly_set') and self._maf_threshold_explicitly_set:
                if not self.filter:
                    raise ValueError("Error in argument: the --maf_threshold is supported only for the --filter operation")
                if self.dbtype and not self.dbtype.startswith('1000g'):
                    raise ValueError("Error in argument: the --maf_threshold is supported only for 1000 Genomes Project data set")
        
        # Validate minqueryfrac parameter
        if hasattr(self, 'minqueryfrac') and self.minqueryfrac is not None:
            # Only validate when minqueryfrac is explicitly set
            if hasattr(self, '_minqueryfrac_explicitly_set') and self._minqueryfrac_explicitly_set:
                if not self.regionanno:
                    raise ValueError("Error in argument: the --minqueryfrac is supported only for the --regionanno operation")
        
        # Validate gff3dbfile parameter
        if hasattr(self, 'gff3dbfile') and self.gff3dbfile is not None:
            if self.dbtype1 != 'gff3':
                raise ValueError("Error in argument: the --gff3dbfile argument can be used only if '--dbtype gff3' is used")
            if not (self.geneanno or self.regionanno):
                raise ValueError("Error in argument: the --gff3dbfile argument is supported only for the --geneanno or --regionanno operation")
        
        # Validate bedfile parameter
        if hasattr(self, 'bedfile') and self.bedfile is not None:
            if self.dbtype1 != 'bed':
                raise ValueError("Error in argument: the --bedfile argument can be used only if '--dbtype bed' is used")
            if not self.regionanno:
                raise ValueError("Error in argument: the --bedfile argument is supported only for the --regionanno operation")
        
        # Validate genericdbfile parameter
        if hasattr(self, 'genericdbfile') and self.genericdbfile is not None:
            if not (self.filter or self.regionanno):
                raise ValueError("Error in argument: the --genericdbfile argument is supported only for the --filter and --region operation")
        
        # Validate vcfdbfile parameter
        if hasattr(self, 'vcfdbfile') and self.vcfdbfile is not None:
            if self.dbtype != 'vcf':
                raise ValueError("Error in argument: the --vcfdbfile argument can be used only if '--dbtype vcf' is used")
        
        # Validate wget parameter
        if hasattr(self, 'wget') and self.wget is not None:
            # Only validate when wget is explicitly set
            if hasattr(self, '_wget_explicitly_set') and self._wget_explicitly_set:
                if not self.downdb:
                    raise ValueError("Error in argument: the --wget argument is supported only for the --downdb operation")
        
        # Validate precedence parameter
        if hasattr(self, 'precedence') and self.precedence is not None:
            # Only validate when precedence is explicitly set
            if hasattr(self, '_precedence_explicitly_set') and self._precedence_explicitly_set:
                if not self.geneanno:
                    raise ValueError("Error in argument: the --precedence argument is supported only for the --geneanno operation")
        
        # Set default values
        if not self.genomebinsize:
            if self.geneanno:
                self.genomebinsize = 100000  # Genes usually span large genomic regions
            else:
                self.genomebinsize = 10000   # MCE, TFBS, miRNA, etc. are small genomic regions
        
        if not self.expandbin:
            self.expandbin = int(2000000 / self.genomebinsize)
        
        # Process zerostart parameter (deprecated)
        if hasattr(self, 'zerostart') and self.zerostart:
            raise ValueError("Error: the -zerostart argument is now obsolete and will no longer be supported in MATCHVAR")

        # Remote dup recognition window limit
        if self.intronic_dup_remote:
            if self.intronic_dup_window < 1 or self.intronic_dup_window > 2000:
                raise ValueError("Error: --intronic_dup_window must be between 1 and 2000")
    
    def _process_dbtype(self, dbtype: str) -> str:
        """Process database type"""
        # Process database type aliases
        dbtype_aliases = {
            'gene': 'refGene',
            'refgene': 'refGene',
            'knowngene': 'knownGene',
            'ensgene': 'ensGene'
        }
        
        if dbtype in dbtype_aliases:
            return dbtype_aliases[dbtype]
        else:
            return dbtype
    
    def run_annotation(self):
        """Run annotation"""
        logger.info(f"Starting annotation, output file: {self.outfile}")
        
        # Print output file name
        if self.geneanno:
            logger.info(f"NOTICE: Output files are written to {self.outfile}.variant_function, {self.outfile}.exonic_variant_function")
        elif self.regionanno:
            logger.info(f"NOTICE: Output file is written to {self.outfile}.{self.buildver}_{self.dbtype1}")
        elif self.filter:
            logger.info(f"NOTICE: Output file with variants matching filtering criteria is written to {self.outfile}.{self.buildver}_{self.dbtype1}_dropped, and output file with other variants is written to {self.outfile}.{self.buildver}_{self.dbtype1}_filtered")
        
        # Check input file line count and adjust thread count
        if self.thread:
            queryfile_line_count, chunk_line_count = self._calculate_chunk_line(self.queryfile, self.thread)
            
            if self.geneanno:
                if queryfile_line_count < self.mingenelinecount:
                    logger.info(f"NOTICE: threading is disabled for gene-based annotation on file with less than {self.mingenelinecount} input lines")
                    self.thread = None
                
                if self.thread and self.thread > self.maxgenethread:
                    logger.info(f"NOTICE: number of threads is reduced to {self.maxgenethread}")
                    self.thread = self.maxgenethread
                    queryfile_line_count, chunk_line_count = self._calculate_chunk_line(self.queryfile, self.thread)
        
        # Start main program
        if self.thread:
            # 使用多线程实现
            self._run_multi_threaded_annotation(queryfile_line_count, chunk_line_count)
        else:
            self._run_single_threaded_annotation()
    
    def _calculate_chunk_line(self, queryfile: str, thread: int) -> Tuple[int, int]:
        """Calculate chunk line count with optimized chunking strategy"""
        try:
            with open(queryfile, 'r', encoding='utf-8') as f:
                line_count = sum(1 for line in f if line.strip() and not line.startswith('#'))
            
            # 优化分块策略：对于大文件使用更大的块，减少线程间切换开销
            if line_count < 1000:
                chunk_line_count = max(1, line_count // thread)
            elif line_count < 10000:
                chunk_line_count = max(100, line_count // thread)
            else:
                # 对于大文件，使用更大的块大小，但不超过5000行
                chunk_line_count = min(5000, max(1000, line_count // thread))
            
            return line_count, chunk_line_count
        except Exception as e:
            logger.error(f"Error calculating chunk line: {e}")
            return 1000, 1000  # Default value
    
    def _run_multi_threaded_annotation(self, queryfile_line_count: int, chunk_line_count: int):
        """Run multi-threaded annotation"""
        threads = []
        
        # 预加载数据库，避免每个线程重复加载
        if self.geneanno:
            logger.info("Pre-loading gene database for multi-threading...")
            gene_db = self._load_gene_database()
        elif self.regionanno:
            logger.info("Pre-loading region database for multi-threading...")
            region_db = self._load_region_database()
        elif self.filter:
            logger.info("Pre-loading filter database for multi-threading...")
            filter_db = self._load_filter_database()
        
        for i in range(self.thread):
            start_line = i * chunk_line_count + 1
            end_line = start_line + chunk_line_count - 1
            if end_line > queryfile_line_count:
                end_line = queryfile_line_count
            
            logger.info(f"NOTICE: Creating new threads for query line {start_line} to {end_line}")
            
            if self.geneanno:
                thread = threading.Thread(
                    target=self._annotate_query_by_gene_thread,
                    args=(f"{self.outfile}.variant_function.{i}", 
                          f"{self.outfile}.exonic_variant_function.{i}",
                          f"{self.outfile}.invalid_input.{i}", 
                          start_line, end_line, i, gene_db)
                )
            elif self.regionanno:
                thread = threading.Thread(
                    target=self._annotate_query_by_region_thread,
                    args=(f"{self.outfile}.{self.buildver}_{self.dbtype1}.{i}",
                          f"{self.outfile}.invalid_input.{i}",
                          start_line, end_line, i, region_db)
                )
            elif self.filter:
                thread = threading.Thread(
                    target=self._filter_query_thread,
                    args=(f"{self.outfile}.{self.buildver}_{self.dbtype1}_filtered.{i}",
                          f"{self.outfile}.{self.buildver}_{self.dbtype1}_dropped.{i}",
                          f"{self.outfile}.invalid_input.{i}",
                          start_line, end_line, i, filter_db)
                )
            
            threads.append(thread)
            thread.start()
        
        # Wait for all threads to complete
        for thread in threads:
            thread.join()
        
        # Merge results
        self._merge_thread_results()
    
    def _run_single_threaded_annotation(self):
        """Run single-threaded annotation"""
        if self.geneanno:
            self._annotate_query_by_gene()
        elif self.regionanno:
            self._annotate_query_by_region()
        elif self.filter:
            self._filter_query()
    
    def _annotate_query_by_gene(self):
        """Gene annotation"""
        logger.info("Starting gene annotation...")
        
        try:
            # Read gene database
            gene_db = self._load_gene_database()
            
            # Process query file
            with open(self.queryfile, 'r', encoding='utf-8') as query_f, \
                 open(f"{self.outfile}.variant_function", 'w', encoding='utf-8') as var_f, \
                 open(f"{self.outfile}.exonic_variant_function", 'w', encoding='utf-8') as exonic_f:
                
                for line in query_f:
                    line = line.strip()
                    if not line or line.startswith('#'):
                        continue
                    
                    # Parse variant
                    variant = self._parse_variant_line(line)
                    if not variant:
                        continue
                    
                    # Annotate variant
                    annotation = self._annotate_variant_by_gene(variant, gene_db)
                    
                    # Write result
                    self._write_gene_annotation(variant, annotation, var_f, exonic_f)
        
        except Exception as e:
            logger.error(f"Gene annotation failed: {e}")
            raise
    
    def _load_gene_database(self) -> Dict:
        """Load gene database"""
        gene_db = {}
        
        try:
            gene_file = os.path.join(self.dbloc, f"{self.buildver}_{self.dbtype1}.txt")
            if not os.path.exists(gene_file):
                logger.warning(f"Gene database file does not exist: {gene_file}")
                return gene_db
            
            # Load kgXref file (if it exists and is knownGene database)
            kgxref = {}
            if self.dbtype1 == 'knownGene':
                kgxreffile = os.path.join(self.dbloc, f"{self.buildver}_kgXref.txt")
                if os.path.exists(kgxreffile):
                    logger.info(f"Loaded knownGene cross-reference file: {kgxreffile}")
                    kgxref = self._load_kgxref(kgxreffile)
                else:
                    logger.warning(f"knownGene cross-reference file does not exist: {kgxreffile}")
            
            with open(gene_file, 'r', encoding='utf-8') as f:
                for line in f:
                    line = line.strip()
                    if not line or line.startswith('#'):
                        continue
                    
                    parts = line.split('\t')
                    
                    # Parse gene information based on database type
                    if self.dbtype1 == 'refGene':
                        if len(parts) < 15:
                            continue
                        # refGene格式: bin, name, chrom, strand, txStart, txEnd, cdsStart, cdsEnd, exonCount, exonStarts, exonEnds, score, name2, cdsStartStat, cdsEndStat
                        gene_info = {
                            'name': parts[1],
                            'chrom': parts[2],
                            'strand': parts[3],
                            'txStart': int(parts[4]),
                            'txEnd': int(parts[5]),
                            'cdsStart': int(parts[6]),
                            'cdsEnd': int(parts[7]),
                            'exonCount': int(parts[8]),
                            'exonStarts': parts[9],
                            'exonEnds': parts[10],
                            'score': parts[11],
                            'name2': parts[12] if len(parts) > 12 else parts[1],
                            'cdsStartStat': parts[13] if len(parts) > 13 else '',
                            'cdsEndStat': parts[14] if len(parts) > 14 else ''
                        }
                    elif self.dbtype1 == 'knownGene':
                        if len(parts) < 10:
                            continue
                        # knownGene format: name, chrom, strand, txStart, txEnd, cdsStart, cdsEnd, exonCount, exonStarts, exonEnds
                        name = parts[0]
                        name2 = kgxref.get(name, name)  # Use kgXref to get gene name
                        gene_info = {
                            'name': name,
                            'chrom': parts[1],
                            'strand': parts[2],
                            'txStart': int(parts[3]),
                            'txEnd': int(parts[4]),
                            'cdsStart': int(parts[5]),
                            'cdsEnd': int(parts[6]),
                            'exonCount': int(parts[7]),
                            'exonStarts': parts[8],
                            'exonEnds': parts[9],
                            'score': '',
                            'name2': name2,
                            'cdsStartStat': '',
                            'cdsEndStat': ''
                        }
                    elif self.dbtype1 == 'ensGene':
                        if len(parts) < 16:
                            continue
                        # ensGene格式: bin, name, chrom, strand, txStart, txEnd, cdsStart, cdsEnd, exonCount, exonStarts, exonEnds, score, name2, cdsStartStat, cdsEndStat, exonFrames
                        gene_info = {
                            'name': parts[1],
                            'chrom': parts[2],
                            'strand': parts[3],
                            'txStart': int(parts[4]),
                            'txEnd': int(parts[5]),
                            'cdsStart': int(parts[6]),
                            'cdsEnd': int(parts[7]),
                            'exonCount': int(parts[8]),
                            'exonStarts': parts[9],
                            'exonEnds': parts[10],
                            'score': parts[11],
                            'name2': parts[12],
                            'cdsStartStat': parts[13],
                            'cdsEndStat': parts[14]
                        }
                    else:
                        # Generic format processing
                        if len(parts) < 11:
                            continue
                        gene_info = {
                            'name': parts[1] if len(parts) > 1 else parts[0],
                            'chrom': parts[2] if len(parts) > 2 else parts[1],
                            'strand': parts[3] if len(parts) > 3 else '+',
                            'txStart': int(parts[4]) if len(parts) > 4 else 0,
                            'txEnd': int(parts[5]) if len(parts) > 5 else 0,
                            'cdsStart': int(parts[6]) if len(parts) > 6 else 0,
                            'cdsEnd': int(parts[7]) if len(parts) > 7 else 0,
                            'exonCount': int(parts[8]) if len(parts) > 8 else 0,
                            'exonStarts': parts[9] if len(parts) > 9 else '',
                            'exonEnds': parts[10] if len(parts) > 10 else '',
                            'score': parts[11] if len(parts) > 11 else '',
                            'name2': parts[12] if len(parts) > 12 else parts[1] if len(parts) > 1 else parts[0],
                            'cdsStartStat': parts[13] if len(parts) > 13 else '',
                            'cdsEndStat': parts[14] if len(parts) > 14 else ''
                        }
                    
                    # Process chromosome name (remove chr prefix)
                    chrom = gene_info['chrom']
                    if chrom.startswith('chr'):
                        chrom = chrom[3:]
                    gene_info['chrom'] = chrom
                    
                    # Convert to 1-based coordinates
                    gene_info['txStart'] += 1
                    gene_info['cdsStart'] += 1
                    if gene_info['exonStarts']:
                        exon_starts = gene_info['exonStarts'].rstrip(',').split(',')
                        gene_info['exonStarts'] = [int(x) + 1 for x in exon_starts]
                    if gene_info['exonEnds']:
                        exon_ends = gene_info['exonEnds'].rstrip(',').split(',')
                        # Note: UCSC genePred's exonEnds is 0-based end (exclusive),
                        # When converted to 1-based closed interval, it should remain int(x) (not +1)
                        gene_info['exonEnds'] = [int(x) for x in exon_ends]
                    
                    # Organize by chromosome
                    if chrom not in gene_db:
                        gene_db[chrom] = []
                    gene_db[chrom].append(gene_info)
        
        except Exception as e:
            logger.error(f"Failed to load gene database: {e}")
        
        return gene_db
    
    def _load_kgxref(self, kgxreffile: str) -> Dict[str, str]:
        """Load knownGene cross-reference file"""
        kgxref = {}
        try:
            with open(kgxreffile, 'r', encoding='utf-8') as f:
                for line in f:
                    line = line.strip()
                    if line.startswith('#'):
                        continue
                    
                    parts = line.split('\t')
                    if len(parts) >= 5:
                        transcript_id = parts[0]
                        gene_name = parts[4]
                        # Remove "Em:" prefix from gene name
                        gene_name = gene_name.replace('Em:', 'Em.')
                        gene_name = gene_name.replace(' ', '')
                        
                        # Process duplicate transcript IDs
                        if transcript_id in kgxref:
                            if kgxref[transcript_id].startswith(('BC', 'AK')):
                                kgxref[transcript_id] = gene_name
                        else:
                            kgxref[transcript_id] = gene_name
            
            logger.info(f"Loaded {len(kgxref)} knownGene cross-references")
        except Exception as e:
            logger.error(f"Failed to load kgXref file: {e}")
        
        return kgxref
    
    def _parse_variant_line(self, line: str) -> Optional[Dict]:
        """Parse variant line"""
        parts = line.split('\t')
        if len(parts) < 5:
            return None
        
        return {
            'chrom': parts[0],
            'start': int(parts[1]),
            'end': int(parts[2]),
            'ref': parts[3],
            'alt': parts[4],
            'original_line': line
        }
    
    def _annotate_variant_by_gene(self, variant: Dict, gene_db: Dict) -> Dict:
        """Annotate single variant"""
        chrom = variant['chrom']
        start = variant['start']
        end = variant['end']
        
        # Process chromosome name (remove chr prefix, consistent with gene database)
        if chrom.startswith('chr'):
            chrom = chrom[3:]
        
        # Find overlapping genes
        overlapping_genes = []
        if chrom in gene_db:
            for gene in gene_db[chrom]:
                if gene['txStart'] <= end and gene['txEnd'] >= start:
                    overlapping_genes.append(gene)
        
        if not overlapping_genes:
            return {
                'function': 'intergenic',
                'gene': 'NA',
                'gene_detail': 'NA',
                'exonic_function': 'NA',
                'aa_change': 'NA'
            }
        
        # Analyze variant type
        annotations = []
        for gene in overlapping_genes:
            annotation = self._analyze_variant_in_gene(variant, gene)
            annotations.append(annotation)
        
        # Select best annotation
        return self._select_best_annotation(annotations)
    
    def _analyze_variant_in_gene(self, variant: Dict, gene: Dict) -> Dict:
        """Analyze variant location in gene"""
        start = variant['start']
        end = variant['end']
        
        # Get gene information
        gene_name = gene.get('name', 'Unknown')
        gene_name2 = gene.get('name2', gene_name)
        strand = gene.get('strand', '+')
        tx_start = gene.get('txStart', 0)
        tx_end = gene.get('txEnd', 0)
        cds_start = gene.get('cdsStart', 0)
        cds_end = gene.get('cdsEnd', 0)
        
        # Parse exon information
        exon_starts = gene.get('exonStarts', [])
        exon_ends = gene.get('exonEnds', [])
        
        # Generate transcript information
        transcript_info = self._generate_transcript_info(variant, gene)

        # Check if it truly falls within any exon interval
        in_exon = False
        for exon_start, exon_end in zip(exon_starts, exon_ends):
            if not (end < exon_start or start > exon_end):
                in_exon = True
                break
        
        # First detect splice site (within |splicing_threshold| range of exon boundary)
        is_splicing, splicing_detail = self._detect_splicing(variant, gene)

        # Determine variant type
        if in_exon:
            # Exon (possibly UTR or CDS)
            if strand == '+':
                # Positive strand gene: transcript 5'→3' direction is the same as the genome direction
                if start < cds_start:
                    function = 'utr5'
                    exonic_function = 'UTR5'
                    aa_change = 'NA'
                elif start > cds_end:
                    function = 'utr3'
                    exonic_function = 'UTR3'
                    aa_change = 'NA'
                else:
                    function = 'exonic'
                    exonic_function = self._determine_exonic_function(variant, gene)
                    aa_change = self._calculate_aa_change(variant, gene)
            else:
                # Negative strand gene: transcript 5'→3' direction is opposite to the genome direction
                if start < cds_start:
                    function = 'utr3'
                    exonic_function = 'UTR3'
                    aa_change = 'NA'
                elif start > cds_end:
                    function = 'utr5'
                    exonic_function = 'UTR5'
                    aa_change = 'NA'
                else:
                    function = 'exonic'
                    exonic_function = self._determine_exonic_function(variant, gene)
                    aa_change = self._calculate_aa_change(variant, gene)

            # exonic;splicing switch (align with Perl -exonicsplicing behavior)
            if self.exonicsplicing and is_splicing:
                function = 'exonic;splicing'
                # Append splice expression to gene_detail
                if splicing_detail:
                    transcript_info = splicing_detail
        elif strand == '+':
            # Positive strand gene: non-exonic region
            if start < cds_start:
                # 5' UTR
                function = 'utr5'
                exonic_function = 'UTR5'
                aa_change = 'NA'
            elif start > cds_end:
                # 3' UTR
                function = 'utr3'
                exonic_function = 'UTR3'
                aa_change = 'NA'
            else:
                # Intron
                if is_splicing:
                    function = 'splicing'
                    exonic_function = 'splicing'
                    # splicing_detail contains c.* or r.spl/UTR5/UTR3 expression
                    if splicing_detail:
                        transcript_info = splicing_detail
                    aa_change = 'NA'
                else:
                    function = 'intronic'
                    exonic_function = 'intronic'
                    # Generate c.HGVS for intron (c.N±offset, only for location expression, not for protein calculation)
                    aa_change = 'NA'
        else:
            # Negative strand gene: non-exonic region
            if start < cds_start:
                # 3' UTR (in negative strand gene, the region with smaller genome coordinates corresponds to the 3' end of the transcript)
                function = 'utr3'
                exonic_function = 'UTR3'
                aa_change = 'NA'
            elif start > cds_end:
                # 5' UTR (in negative strand gene, the region with larger genome coordinates corresponds to the 5' end of the transcript)
                function = 'utr5'
                exonic_function = 'UTR5'
                aa_change = 'NA'
            else:
                # Intron
                if is_splicing:
                    function = 'splicing'
                    exonic_function = 'splicing'
                    # splicing_detail contains c.* or r.spl/UTR5/UTR3 expression
                    if splicing_detail:
                        transcript_info = splicing_detail
                    aa_change = 'NA'
                else:
                    function = 'intronic'
                    exonic_function = 'intronic'
                    # 为内含子生成 c.HGVS（c.N±offset，仅作定位表达，不计算蛋白）
                    aa_change = 'NA'
        
        return {
            'function': function,
            'gene': gene_name2,
            'gene_detail': transcript_info,
            'exonic_function': exonic_function,
            'aa_change': aa_change
        }

    def _detect_splicing(self, variant: Dict, gene: Dict) -> Tuple[bool, str]:
        """Detect splice variant (near exon boundary).
        Return (whether it is splicing, detail string).
        Rules:
        - If the variant site is within ±splicing_threshold range of the exon boundary or crosses the boundary, it is considered splicing
        - If it is SNV and in the adjacent CDS region, output transcript:exonN:c.<cpos>+/-<offset>Ref>Alt
        - If it is indel, output transcript:exonN:r.spl
        - If it is near UTR, output transcript:exonN:UTR5/UTR3
        """
        try:
            start = variant['start']
            end = variant['end']
            ref = variant['ref']
            alt = variant['alt']
            strand = gene.get('strand', '+')
            cds_start = gene.get('cdsStart', 0)
            cds_end = gene.get('cdsEnd', 0)
            exon_starts = gene.get('exonStarts', [])
            exon_ends = gene.get('exonEnds', [])
            threshold = self.indel_splicing_threshold if (len(ref) != len(alt)) and (self.indel_splicing_threshold is not None) else self.splicing_threshold

            # Standard transcript ID
            # gene.get('name') is transcript ID, gene.get('name2') is gene name
            transcript_id = self._get_standard_transcript_id(gene.get('name2', 'Unknown'), gene.get('name', 'Unknown'), self.dbtype1)

            # Build CDS segments and calculate the cumulative cDNA intervals for cpos
            coding_segments = []  # (seg_start, seg_end, exon_index)
            for idx, (es, ee) in enumerate(zip(exon_starts, exon_ends)):
                if strand == '-':
                    # For negative strand genes: cds_start > cds_end
                    # We need to find the intersection of exon with CDS range [cds_end, cds_start]
                    seg_s = max(es, cds_end)  # Use cds_end as the lower bound
                    seg_e = min(ee, cds_start)  # Use cds_start as the upper bound
                else:
                    # For positive strand genes: cds_start < cds_end
                    seg_s = max(es, cds_start)
                    seg_e = min(ee, cds_end)
                if seg_s <= seg_e:
                    coding_segments.append((seg_s, seg_e, idx))

            # Sort transcript direction
            if strand == '+':
                coding_segments.sort(key=lambda x: x[0])
            else:
                coding_segments.sort(key=lambda x: x[0], reverse=True)

            cumulative = []  # (exon_index, seg_start, seg_end, c_start, c_end)
            acc = 0
            for seg_s, seg_e, idx in coding_segments:
                length = seg_e - seg_s + 1
                c_start = acc + 1
                c_end = acc + length
                cumulative.append((idx, seg_s, seg_e, c_start, c_end))
                acc += length

            # Check each exon boundary's adjacent region
            for i, (exon_start, exon_end) in enumerate(zip(exon_starts, exon_ends)):
                # Adjacent range
                left_region = (exon_start - threshold, exon_start - 1)
                right_region = (exon_end + 1, exon_end + threshold)
                overlaps_boundary = not (end < exon_start or start > exon_end) and (start <= exon_end and end >= exon_start)

                near_left = (start >= left_region[0] and start <= left_region[1]) or (end >= left_region[0] and end <= left_region[1])
                near_right = (start >= right_region[0] and start <= right_region[1]) or (end >= right_region[0] and end <= right_region[1])

                if overlaps_boundary or near_left or near_right:
                    # Calculate exon number, considering the positive and negative strand direction of the gene
                    # Create exon coordinate list for sorting
                    exon_coords = list(zip(exon_starts, exon_ends))
                    
                    # Sort exons based on the gene strand direction
                    if strand == '-':
                        # Negative strand gene: sort by transcript 5'→3' direction (genome coordinates from large to small)
                        exon_coords.sort(key=lambda x: x[0], reverse=True)
                    else:
                        # Positive strand gene: sort by transcript 5'→3' direction (genome coordinates from small to large)
                        exon_coords.sort(key=lambda x: x[0])
                    
                    # Find the current exon in the sorted position
                    current_exon = (exon_start, exon_end)
                    exon_num = 0
                    for idx, (es, ee) in enumerate(exon_coords):
                        if es == exon_start and ee == exon_end:
                            exon_num = idx + 1
                            break

                    # UTR near: if the site is outside the CDS, mark as UTR
                    if strand == '+':
                        # Positive strand gene
                        if start < cds_start:
                            return True, f"{transcript_id}:exon{exon_num}:UTR5"
                        if start > cds_end:
                            return True, f"{transcript_id}:exon{exon_num}:UTR3"
                    else:
                        # Negative strand gene
                        if start < cds_start:
                            return True, f"{transcript_id}:exon{exon_num}:UTR3"
                        if start > cds_end:
                            return True, f"{transcript_id}:exon{exon_num}:UTR5"

                    # indel: use r.spl
                    if len(ref) != len(alt):
                        return True, f"{transcript_id}:exon{exon_num}:r.spl"

                    # SNV: construct cpos and +/- offset
                    # Determine the cDNA start and end of the exon in CDS
                    cpos_exon_start = cpos_exon_end = None
                    # Create a mapping from exon coordinates to indices
                    exon_to_idx = {}
                    for idx, (es, ee) in enumerate(zip(exon_starts, exon_ends)):
                        exon_to_idx[(es, ee)] = idx
                    
                    # Use the current exon coordinates to find the corresponding index
                    current_exon_idx = exon_to_idx.get((exon_start, exon_end))
                    if current_exon_idx is not None:
                        for idx, seg_s, seg_e, c_s, c_e in cumulative:
                            if idx == current_exon_idx:
                                cpos_exon_start, cpos_exon_end = c_s, c_e
                                break
                    if cpos_exon_start is None:
                        # Non-coding exon (no CDS), process as UTR
                        if strand == '+':
                            # Positive strand gene
                            if start < cds_start:
                                return True, f"{transcript_id}:exon{exon_num}:UTR5"
                            if start > cds_end:
                                return True, f"{transcript_id}:exon{exon_num}:UTR3"
                        else:
                            # Negative strand gene
                            if start < cds_start:
                                return True, f"{transcript_id}:exon{exon_num}:UTR3"
                            if start > cds_end:
                                return True, f"{transcript_id}:exon{exon_num}:UTR5"
                        return True, f"{transcript_id}:exon{exon_num}:r.spl"

                    # Calculate the offset to the boundary and decide to use + or -
                    # Consider the transcript direction:
                    # Positive strand: use '-' near exon_start, use '+' near exon_end
                    # Negative strand: opposite (because the transcript 5'→3' direction is opposite to the genome direction)
                    if strand == '+':
                        if near_left:
                            boundary_cdna = cpos_exon_start
                            offset = exon_start - start  # >=1
                            sign = '-'
                        elif near_right or overlaps_boundary:
                            boundary_cdna = cpos_exon_end
                            offset = start - exon_end if start > exon_end else 0
                            offset = max(1, offset)
                            sign = '+'
                        else:
                            continue
                    else:
                        if near_left:
                            # Negative strand: the left adjacent region of the genome corresponds to the '+' end of the transcript
                            boundary_cdna = cpos_exon_end
                            offset = exon_start - start  # >=1
                            sign = '+'
                        elif near_right or overlaps_boundary:
                            boundary_cdna = cpos_exon_start
                            offset = start - exon_end if start > exon_end else 0
                            offset = max(1, offset)
                            sign = '-'
                        else:
                            continue

                    # Process bases based on strand
                    if strand == '+':
                        ref_base, alt_base = ref, alt
                    else:
                        ref_base = self._reverse_complement(ref)
                        alt_base = self._reverse_complement(alt)

                    return True, f"{transcript_id}:exon{exon_num}:c.{boundary_cdna}{sign}{offset}{ref_base}>{alt_base}"

            return False, ''
        except Exception:
            return False, ''
    
    def _generate_transcript_info(self, variant: Dict, gene: Dict) -> str:
        """Generate transcript information, format: transcript ID:exon number:c.position reference base>alternative base"""
        try:
            gene_name = gene.get('name', 'Unknown')
            strand = gene.get('strand', '+')
            tx_start = gene.get('txStart', 0)
            cds_start = gene.get('cdsStart', 0)
            cds_end = gene.get('cdsEnd', 0)
            
            # Parse exon information
            exon_starts = gene.get('exonStarts', [])
            exon_ends = gene.get('exonEnds', [])
            
            if not exon_starts or not exon_ends:
                return f"{gene_name}:NA"
            
            # Find the exon number of the variant
            variant_start = variant['start']
            variant_end = variant['end']
            ref = variant['ref']
            alt = variant['alt']
            # Normalize variant alleles, process '.'、'-'、'*' placeholders
            ref, alt = self._normalize_variant_alleles(ref, alt)
            
            # Calculate exon number, considering the positive and negative strand direction of the gene
            exon_num = 0
            # Create exon coordinate list for sorting
            exon_coords = list(zip(exon_starts, exon_ends))
            
            # Sort exons based on the gene strand direction
            if strand == '-':
                # Negative strand gene: sort by transcript 5'→3' direction (genome coordinates from large to small)
                exon_coords.sort(key=lambda x: x[0], reverse=True)
            else:
                # Positive strand gene: sort by transcript 5'→3' direction (genome coordinates from small to large)
                exon_coords.sort(key=lambda x: x[0])
            
            # Find the exon number of the variant
            for i, (exon_start, exon_end) in enumerate(exon_coords):
                if exon_start <= variant_start <= exon_end:
                    exon_num = i + 1
                    break
            
            # Get standard transcript ID
            # gene_name is actually the transcript ID, gene name is in gene.get('name2')
            original_transcript_id = gene_name  # This is actually transcript ID
            gene_name_actual = gene.get('name2', gene_name)  # This is the actual gene name
            standard_transcript_id = self._get_standard_transcript_id(gene_name_actual, original_transcript_id, self.dbtype1)
            
            # First determine if it is in any exon
            in_exon = False
            upstream_exon_idx = -1
            downstream_exon_idx = -1
            # Create exon coordinate list for sorting
            exon_coords = list(zip(exon_starts, exon_ends))
            
            # Sort exons based on the gene strand direction
            if strand == '-':
                # Negative strand gene: sort by transcript 5'→3' direction (genome coordinates from large to small)
                exon_coords.sort(key=lambda x: x[0], reverse=True)
            else:
                # Positive strand gene: sort by transcript 5'→3' direction (genome coordinates from small to large)
                exon_coords.sort(key=lambda x: x[0])
            
            for i, (es, ee) in enumerate(exon_coords):
                if es <= variant_start <= ee:
                    in_exon = True
                    exon_num = i + 1
                    break
                if ee < variant_start:
                    upstream_exon_idx = i
                if downstream_exon_idx == -1 and es > variant_start:
                    downstream_exon_idx = i

            # Non-coding region (not in CDS range or not in any exon)
            if not in_exon:
                # Check if it is in UTR region
                # In refGene.txt, cds_start and cds_end are always cds_start < cds_end regardless of strand
                # For negative strand genes, we need to interpret them correctly for CDS position calculation
                if strand == '-':
                    # Negative strand: CDS range is [cds_start, cds_end] in genomic coordinates
                    # But 5' end corresponds to cds_end (larger coordinate) and 3' end to cds_start (smaller coordinate)
                    in_cds_range = cds_start <= variant_start <= cds_end
                else:
                    # Positive strand: CDS range is [cds_start, cds_end]
                    in_cds_range = cds_start <= variant_start <= cds_end
                
                if not in_cds_range:
                    # UTR intronic: need to generate c.HGVS format
                    # This logic will be handled in the later code
                    # Directly call UTR intronic processing logic
                    return self._generate_utr_intronic_hgvs(variant_start, ref, alt, strand, 
                                                          exon_starts, exon_ends, cds_start, cds_end, 
                                                          standard_transcript_id)
                else:
                    # Non-UTR intronic: select the nearest boundary (the smaller the number after ±, the higher the priority), output c.N+/-offset expression
                    # First build CDS segments and calculate the CDS boundaries of each exon in cDNA
                    coding_segments = []  # (seg_start, seg_end, exon_index)
                for idx, (es, ee) in enumerate(zip(exon_starts, exon_ends)):
                    if strand == '-':
                        # For negative strand genes: CDS range is [cds_start, cds_end] in genomic coordinates
                        # We need to find the intersection of exon with CDS range
                        seg_s = max(es, cds_start)  # Use cds_start as the lower bound
                        seg_e = min(ee, cds_end)  # Use cds_end as the upper bound
                    else:
                        # For positive strand genes: cds_start < cds_end
                        seg_s = max(es, cds_start)
                        seg_e = min(ee, cds_end)
                    if seg_s <= seg_e:
                        coding_segments.append((seg_s, seg_e, idx))
                if not coding_segments:
                    return f"{standard_transcript_id}:intronic"

                # Sort by transcript 5'→3'
                if strand == '+':
                    coding_segments.sort(key=lambda x: x[0])
                else:
                    coding_segments.sort(key=lambda x: x[0], reverse=True)

                # Calculate the cumulative length of each cDNA segment
                cumulative = []  # (exon_index, seg_start, seg_end, c_start, c_end)
                acc = 0
                for seg_s, seg_e, idx in coding_segments:
                    length = seg_e - seg_s + 1
                    c_start = acc + 1
                    c_end = acc + length
                    cumulative.append((idx, seg_s, seg_e, c_start, c_end))
                    acc += length
                # Traverse adjacent CDS segments, locate the intronic interval where the variant is located, and compare the distance between the two boundaries
                best = None  # (offset, sign, N_base, intron_info)
                for i in range(len(cumulative) - 1):
                    idx_p, s_p, e_p, c_s_p, c_e_p = cumulative[i]
                    idx_n, s_n, e_n, c_s_n, c_e_n = cumulative[i + 1]
                    # Whether the variant is between e_p and s_n
                    low = min(e_p, s_n)
                    high = max(e_p, s_n)
                    if not (low <= variant_start <= high):
                        continue
                    candidates = []
                    if strand == '+':
                        # Positive strand: the previous exon end is '+', the next exon start is '-'
                        candidates.append((abs(variant_start - e_p), '+', c_e_p, f"intron{i+1}"))
                        candidates.append((abs(s_n - variant_start), '-', c_s_n, f"intron{i+1}"))
                    else:
                        # Negative strand: the previous exon end is '+', the next exon start is '-'
                        # Note: For negative strand, we need to use the correct coordinate variables
                        candidates.append((abs(variant_start - s_p), '+', c_e_p, f"intron{i+1}"))
                        candidates.append((abs(e_n - variant_start), '-', c_s_n, f"intron{i+1}"))
                    # Select offset smaller; if equal, prioritize '+'
                    candidates.sort(key=lambda x: (x[0], 0 if x[1] == '+' else 1))
                    best = candidates[0]
                    break

                if not best:
                    return f"{standard_transcript_id}:intronic"

                offset, sign, N_base, intron_info = best
                # Process bases based on strand
                ref_b, alt_b = ref, alt
                if strand == '-':
                    ref_b = self._reverse_complement(ref_b) if ref_b else ''
                    alt_b = self._reverse_complement(alt_b) if alt_b else ''

                # Organize intronic c.HGVS based on variant type (SNV/deletion/insertion/substitution)
                if len(ref_b) == 1 and len(alt_b) == 1:
                    return f"{standard_transcript_id}:{intron_info}:c.{N_base}{sign}{int(offset)}{ref_b}>{alt_b}:p.?"
                elif len(ref_b) > len(alt_b):
                    if len(alt_b) == 0:
                        return f"{standard_transcript_id}:{intron_info}:c.{N_base}{sign}{int(offset)}del{ref_b}:p.?"
                    else:
                        return f"{standard_transcript_id}:{intron_info}:c.{N_base}{sign}{int(offset)}del{ref_b}ins{alt_b}:p.?"
                elif len(ref_b) < len(alt_b):
                    if len(ref_b) == 0:
                        # Insertion: try to identify dup remotely, otherwise output ins
                        a = int(offset)
                        b = max(1, a - 1)
                        if sign == '+':
                            o1, o2 = sorted([a, b])
                        else:
                            o1, o2 = sorted([a, b], reverse=True)
                        if self.intronic_dup_remote:
                            try:
                                if self._is_intronic_dup_by_remote(variant['chrom'], variant['start'], alt_b, strand):
                                    return f"{standard_transcript_id}:{intron_info}:c.{N_base}{sign}{o1}_{N_base}{sign}{o2}dup{alt_b}:p.?"
                            except Exception:
                                pass
                        return f"{standard_transcript_id}:{intron_info}:c.{N_base}{sign}{o1}_{N_base}{sign}{o2}ins{alt_b}:p.?"
                    else:
                        return f"{standard_transcript_id}:{intron_info}:c.{N_base}{sign}{int(offset)}del{ref_b}ins{alt_b}:p.?"
                else:
                    return f"{standard_transcript_id}:intronic:p.?"
            
            # Check if it is in UTR region (exon or intron), if so, use special c.HGVS format
            # In refGene.txt, cds_start and cds_end are always cds_start < cds_end regardless of strand
            # For negative strand genes, we need to interpret them correctly for CDS position calculation
            if strand == '-':
                # Negative strand: CDS range is [cds_start, cds_end] in genomic coordinates
                # But 5' end corresponds to cds_end (larger coordinate) and 3' end to cds_start (smaller coordinate)
                in_cds_range = cds_start <= variant_start <= cds_end
            else:
                # Positive strand: CDS range is [cds_start, cds_end]
                in_cds_range = cds_start <= variant_start <= cds_end
            
            if not in_cds_range:
                # UTR region: need to determine if it is in UTR exon or intron
                
                # First check if it is in any exon
                in_utr_exon = False
                for exon_start, exon_end in zip(exon_starts, exon_ends):
                    if exon_start <= variant_start <= exon_end:
                        in_utr_exon = True
                        break
                
                # If it is not in any exon, but close to the exon, also process as UTR intronic
                if not in_utr_exon:
                    # Check if it is in UTR intronic (close to the exon)
                    min_distance_to_exon = float('inf')
                    for exon_start, exon_end in zip(exon_starts, exon_ends):
                        dist_to_start = abs(variant_start - exon_start)
                        dist_to_end = abs(variant_start - exon_end)
                        min_distance = min(dist_to_start, dist_to_end)
                        if min_distance < min_distance_to_exon:
                            min_distance_to_exon = min_distance
                    
                    # If the distance to the nearest exon is less than a certain threshold, process as UTR intronic
                    if min_distance_to_exon <= 2000:  # 2kb threshold
                        in_utr_exon = False  # Ensure processing as intronic
                    else:
                        # Too far away, possibly not a variant of this gene
                        in_utr_exon = True  # Process as exon, but generate default format
                
                if in_utr_exon:
                    # UTR exonic: according to UTR region naming rules
                    if strand == '+':
                        # Positive strand gene
                        if variant_start < cds_start:
                            # UTR5: the first nucleotide of the start codon is +1, 5'UTR is negative
                            # Calculate the distance to the start codon (negative)
                            utr_position = variant_start - cds_start
                            if len(ref) == 1 and len(alt) == 1:
                                hgvs = f"{standard_transcript_id}:UTR5:c.-{abs(utr_position)}{ref}>{alt}:p.?"
                            else:
                                hgvs = f"{standard_transcript_id}:UTR5:c.-{abs(utr_position)}delins{alt}:p.?"
                        else:
                            # UTR3: the last nucleotide of the stop codon is the boundary, 3'UTR is *1, *2, *3...
                            # Calculate the distance to the stop codon (positive)
                            utr_position = variant_start - cds_end
                            if len(ref) == 1 and len(alt) == 1:
                                hgvs = f"{standard_transcript_id}:UTR3:c.*{abs(utr_position)}{ref}>{alt}:p.?"
                            else:
                                hgvs = f"{standard_transcript_id}:UTR3:c.*{abs(utr_position)}delins{alt}:p.?"
                    else:
                        # Negative strand gene
                        if variant_start > cds_end:
                            # Negative strand gene: the region with larger genome coordinates corresponds to the 5' end of the transcript
                            # So this is UTR5, the start codon is +1, 5'UTR is negative
                            utr_position = variant_start - cds_end
                            ref_seq = self._reverse_complement(ref) if ref else ''
                            alt_seq = self._reverse_complement(alt) if alt else ''
                            if len(ref_seq) == 1 and len(alt_seq) == 1:
                                hgvs = f"{standard_transcript_id}:UTR5:c.-{abs(utr_position)}{ref_seq}>{alt_seq}:p.?"
                            else:
                                hgvs = f"{standard_transcript_id}:UTR5:c.-{abs(utr_position)}delins{alt_seq}:p.?"
                        else:
                            # Negative strand gene: the region with smaller genome coordinates corresponds to the 3' end of the transcript
                            # So this is UTR3, the stop codon is the boundary, 3'UTR is *1, *2, *3...
                            utr_position = variant_start - cds_start
                            ref_seq = self._reverse_complement(ref) if ref else ''
                            alt_seq = self._reverse_complement(alt) if alt else ''
                            if len(ref_seq) == 1 and len(alt_seq) == 1:
                                hgvs = f"{standard_transcript_id}:UTR3:c.*{abs(utr_position)}{ref_seq}>{alt_seq}:p.?"
                            else:
                                hgvs = f"{standard_transcript_id}:UTR3:c.*{abs(utr_position)}delins{alt_seq}:p.?"
                else:
                    # UTR intronic: based on the neighboring exon number naming
                    # This logic has been moved to the _generate_utr_intronic_hgvs method
                    return self._generate_utr_intronic_hgvs(variant_start, ref, alt, strand, 
                                                          exon_starts, exon_ends, cds_start, cds_end, 
                                                          standard_transcript_id)
                
                return hgvs
            else:
                # CDS region: accurately calculate the cDNA position
                coding_segments = []  # (seg_start, seg_end)
                for es, ee in zip(exon_starts, exon_ends):
                    if strand == '-':
                        # For negative strand genes: CDS range is [cds_start, cds_end] in genomic coordinates
                        seg_s = max(es, cds_start)  # Use cds_start as the lower bound
                        seg_e = min(ee, cds_end)  # Use cds_end as the upper bound
                    else:
                        # For positive strand genes: cds_start < cds_end
                        seg_s = max(es, cds_start)
                        seg_e = min(ee, cds_end)
                    if seg_s <= seg_e:
                        coding_segments.append((seg_s, seg_e))

                if not coding_segments:
                    return f"{standard_transcript_id}:NA"

                # Sort by transcript 5'→3'
                if strand == '+':
                    coding_segments.sort(key=lambda x: x[0])
                else:
                    coding_segments.sort(key=lambda x: x[0], reverse=True)

                # Calculate the cumulative cDNA length and locate the segment
                acc = 0
                cds_pos = None
                for seg_s, seg_e in coding_segments:
                    seg_len = seg_e - seg_s + 1
                    if seg_s <= variant_start <= seg_e:
                        # Offset within the segment
                        if strand == '+':
                            offset = variant_start - seg_s
                            cds_pos = acc + offset + 1
                        else:
                            # For negative strand genes, calculate position from the 5' end of CDS
                            # For negative strand, 5' end corresponds to cds_end (larger genomic coordinate)
                            # We need to calculate how many bases from the 5' end
                            cds_pos = cds_end - variant_start + 1
                        break
                    acc += seg_len

                if cds_pos is None:
                    # Not hit CDS segment (theoretically should not happen, because it is in in_exon and in CDS range), fallback
                    cds_pos = 1

                # Process bases based on strand
                if strand == '+':
                    ref_seq = ref
                    alt_seq = alt
                else:
                    ref_seq = self._reverse_complement(ref) if ref else ''
                    alt_seq = self._reverse_complement(alt) if alt else ''
                
                # Generate HGVS format - consistent with _calculate_aa_change
                if len(ref_seq) == 1 and len(alt_seq) == 1:
                    # Single nucleotide substitution
                    hgvs = f"{standard_transcript_id}:exon{exon_num}:c.{cds_pos}{ref_seq}>{alt_seq}"
                elif len(ref_seq) > len(alt_seq):
                    # Deletion or deletion+insertion (delins, alt is not empty)
                    if len(alt_seq) == 0:
                        if strand == '+':
                            start_pos = cds_pos
                            end_pos = cds_pos + len(ref_seq) - 1
                        else:
                            start_pos = cds_pos - len(ref_seq) + 1
                            end_pos = cds_pos
                        hgvs = f"{standard_transcript_id}:exon{exon_num}:c.{start_pos}_{end_pos}del{ref_seq}"
                    else:
                        if strand == '+':
                            start_pos = cds_pos
                            end_pos = cds_pos + len(ref_seq) - 1
                        else:
                            start_pos = cds_pos - len(ref_seq) + 1
                            end_pos = cds_pos
                        hgvs = f"{standard_transcript_id}:exon{exon_num}:c.{start_pos}_{end_pos}del{ref_seq}ins{alt_seq}"
                elif len(ref_seq) < len(alt_seq):
                    # Insertion or delins (ref is not empty)
                    if len(ref_seq) == 0:
                        # Pure insertion: check the original variant information
                        original_ref = variant.get('ref', '')
                        original_alt = variant.get('alt', '')
                        
                        # If the original REF is "-", prioritize using delins format
                        if original_ref == '-' or original_ref == '.' or original_ref == '*':
                            # When REF is "-", use a single position instead of a position range
                            hgvs = f"{standard_transcript_id}:exon{exon_num}:c.{cds_pos}delins{alt_seq}"
                        else:
                            # Regular insertion format
                            if strand == '+':
                                hgvs = f"{standard_transcript_id}:exon{exon_num}:c.{cds_pos}_{cds_pos+1}ins{alt_seq}"
                            else:
                                hgvs = f"{standard_transcript_id}:exon{exon_num}:c.{cds_pos-1}_{cds_pos}ins{alt_seq}"
                    else:
                        if strand == '+':
                            start_pos = cds_pos
                            end_pos = cds_pos + len(ref_seq) - 1
                        else:
                            start_pos = cds_pos - len(ref_seq) + 1
                            end_pos = cds_pos
                        hgvs = f"{standard_transcript_id}:exon{exon_num}:c.{start_pos}_{end_pos}del{ref_seq}ins{alt_seq}"
                else:
                    # Other cases
                    hgvs = f"{standard_transcript_id}:exon{exon_num}:c.{cds_pos}{ref_seq}>{alt_seq}"
            
            return hgvs
            
        except Exception as e:
            logger.warning(f"Failed to generate transcript information: {e}")
            return f"{gene.get('name', 'Unknown')}:NA"
    
    def _calculate_cds_cdna_position(self, genomic_pos: int, gene: Dict) -> int:
        """Calculate the position of the specified genomic position in cDNA"""
        try:
            strand = gene.get('strand', '+')
            exon_starts = gene.get('exonStarts', [])
            exon_ends = gene.get('exonEnds', [])
            
            # Create exon coordinate list
            exon_coords = list(zip(exon_starts, exon_ends))
            
            # Sort exons based on the gene strand direction
            if strand == '-':
                exon_coords.sort(key=lambda x: x[0], reverse=True)
            else:
                exon_coords.sort(key=lambda x: x[0])
            
            # Calculate the cumulative cDNA length
            acc = 0
            for exon_start, exon_end in exon_coords:
                if exon_start <= genomic_pos <= exon_end:
                    # Find the position, calculate the offset within the exon
                    if strand == '+':
                        offset = genomic_pos - exon_start
                    else:
                        offset = exon_end - genomic_pos
                    return acc + offset + 1
                acc += exon_end - exon_start + 1
            
            # If not found, return 0
            return 0
        except Exception:
            return 0
    
    def _calculate_variant_cdna_position(self, variant_pos: int, gene: Dict) -> int:
        """Calculate the position of the variant in cDNA"""
        return self._calculate_cds_cdna_position(variant_pos, gene)
    
    def _determine_exonic_function(self, variant: Dict, gene: Dict) -> str:
        """Determine the functional impact of the exonic variant"""
        # Normalize variant alleles, process '.'、'-'、'*' placeholders
        ref, alt = self._normalize_variant_alleles(variant['ref'], variant['alt'])
        
        if len(ref) == len(alt):
            if len(ref) == 1:
                return 'nonsynonymous_SNV'
            else:
                return 'nonsynonymous_indel'
        elif len(ref) > len(alt):
            return 'frameshift_deletion'
        else:
            return 'frameshift_insertion'
    
    def _calculate_aa_change(self, variant: Dict, gene: Dict) -> str:
        """Calculate amino acid changes, generate HGVS format, consider strand direction"""
        try:
            # Get gene information
            gene_name = gene.get('name2', 'Unknown')
            strand = gene.get('strand', '+')  # Get strand direction
            
            # Calculate CDS position
            cds_start = gene.get('cdsStart', 0)
            cds_end = gene.get('cdsEnd', 0)
            
            # Calculate the position of the variant in CDS
            # In refGene.txt, cds_start and cds_end are always cds_start < cds_end regardless of strand
            # For negative strand genes, we need to interpret them correctly for CDS position calculation
            if strand == '-':
                # Negative strand: CDS range is [cds_start, cds_end] in genomic coordinates
                # But 5' end corresponds to cds_end (larger coordinate) and 3' end to cds_start (smaller coordinate)
                if variant['start'] < cds_start or variant['start'] > cds_end:
                    return f"{gene_name}:p.?"
            else:
                # Positive strand: CDS range is [cds_start, cds_end]
                if variant['start'] < cds_start or variant['start'] > cds_end:
                    return f"{gene_name}:p.?"
            
            # Accurately calculate the cDNA position: based on the cumulative length of the CDS after splicing
            exon_starts = gene.get('exonStarts', [])
            exon_ends = gene.get('exonEnds', [])
            vpos = variant['start']
            coding_segments = []
            for es, ee in zip(exon_starts, exon_ends):
                if strand == '-':
                    # For negative strand genes: CDS range is [cds_start, cds_end] in genomic coordinates
                    seg_s = max(es, cds_start)  # Use cds_start as the lower bound
                    seg_e = min(ee, cds_end)  # Use cds_end as the upper bound
                else:
                    # For positive strand genes: cds_start < cds_end
                    seg_s = max(es, cds_start)
                    seg_e = min(ee, cds_end)
                if seg_s <= seg_e:
                    coding_segments.append((seg_s, seg_e))
            if strand == '+':
                coding_segments.sort(key=lambda x: x[0])
            else:
                coding_segments.sort(key=lambda x: x[0], reverse=True)
            acc = 0
            cds_pos = None
            for seg_s, seg_e in coding_segments:
                seg_len = seg_e - seg_s + 1
                if seg_s <= vpos <= seg_e:
                    if strand == '+':
                        offset = vpos - seg_s
                        cds_pos = acc + offset + 1
                    else:
                        # For negative strand genes, calculate position from the 5' end of CDS
                        # For negative strand, 5' end corresponds to cds_end (larger genomic coordinate)
                        cds_pos = cds_end - vpos + 1
                    break
                acc += seg_len
            if cds_pos is None:
                cds_pos = 1
            # Normalize variant alleles, process '.'、'-'、'*' placeholders
            ref, alt = self._normalize_variant_alleles(variant['ref'], variant['alt'])
            
            # Process bases based on strand
            if strand == '+':
                ref_seq = ref
                alt_seq = alt
            else:
                ref_seq = self._reverse_complement(ref)
                alt_seq = self._reverse_complement(alt)
            
            # Validate HGVS format
            if not self._validate_hgvs_format(ref_seq, alt_seq):
                logger.warning(f"HGVS格式校验失败: ref={ref_seq}, alt={alt_seq}")
                return f"{gene_name}:p.?"
            
            # Generate HGVS format based on variant type
            if ref_seq == alt_seq:
                return f"{gene_name}:p.?"
            elif len(ref_seq) == 1 and len(alt_seq) == 1:
                # Single nucleotide substitution - standard HGVS format
                return f"{gene_name}:c.{cds_pos}{ref_seq}>{alt_seq}"
            elif len(ref_seq) > len(alt_seq):
                # Deletion
                if len(alt_seq) == 0:
                    # Pure deletion
                    end_pos = cds_pos + len(ref_seq) - 1
                    return f"{gene_name}:c.{cds_pos}_{end_pos}del{ref_seq}"
                else:
                    # Substitution deletion
                    end_pos = cds_pos + len(ref_seq) - 1
                    return f"{gene_name}:c.{cds_pos}_{end_pos}del{ref_seq}ins{alt_seq}"
            elif len(ref_seq) < len(alt_seq):
                # Insertion
                if len(ref_seq) == 0:
                    # Pure insertion
                    return f"{gene_name}:c.{cds_pos}_{cds_pos+1}ins{alt_seq}"
                else:
                    # Substitution insertion
                    end_pos = cds_pos + len(ref_seq) - 1
                    return f"{gene_name}:c.{cds_pos}_{end_pos}del{ref_seq}ins{alt_seq}"
            else:
                # Other cases
                return f"{gene_name}:c.{cds_pos}{ref_seq}>{alt_seq}"
                
        except Exception as e:
            logger.warning(f"Failed to calculate amino acid changes: {e}")
            return f"{gene.get('name2', 'Unknown')}:p.?"
    
    def _reverse_complement(self, sequence: str) -> str:
        """Calculate the reverse complementary sequence"""
        complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N'}
        return ''.join(complement.get(base, base) for base in reversed(sequence))
    
    def _normalize_variant_alleles(self, ref: str, alt: str) -> tuple:
        """Normalize variant alleles, process '.'、'-'、'*' placeholders"""
        # Normalize '.'、'-'、'*'
        normalized_ref = ref.replace('.', '-').replace('*', '-') if ref else ''
        normalized_alt = alt.replace('.', '-').replace('*', '-') if alt else ''
        
        # Process special cases
        if normalized_ref == '-' and normalized_alt != '-':
            # ref='.'、'-'、'*' alt='A' -> insertion
            return ('', normalized_alt)
        elif normalized_alt == '-' and normalized_ref != '-':
            # ref='A' alt='.'、'-'、'*' -> deletion
            return (normalized_ref, '')
        elif normalized_ref == '-' and normalized_alt == '-':
            # ref='.'、'-'、'*' alt='.'、'-'、'*' -> invalid variant
            return ('', '')
        else:
            # Other cases remain the same
            return (normalized_ref, normalized_alt)
    
    def _validate_hgvs_format(self, ref_seq: str, alt_seq: str) -> bool:
        """Validate the validity of the HGVS format"""
        # Allow deletion placeholder'-', gVCF placeholder'*' and VCF placeholder'.'
        valid_bases = {'A', 'T', 'C', 'G', 'N', '-', '*', '.'}
        
        for seq in [ref_seq, alt_seq]:
            if not all(base in valid_bases for base in seq):
                return False
        
        # Check sequence length
        if len(ref_seq) == 0 and len(alt_seq) == 0:
            return False
        
        return True
    
    def _select_best_annotation(self, annotations: List[Dict]) -> Dict:
        """Select the best annotation"""
        if not annotations:
            return {
                'function': 'intergenic',
                'gene': 'NA',
                'gene_detail': 'NA',
                'exonic_function': 'NA',
                'aa_change': 'NA'
            }
        
        # Always combine all transcript annotations first
        combined_annotation = self._combine_transcript_annotations(annotations)
        
        # If MANE transcript filtering is enabled, filter the combined annotation
        if self.use_mane_transcript:
            combined_annotation = self._filter_mane_annotation(combined_annotation)
        
        return combined_annotation
    
    def _filter_mane_annotation(self, annotation: Dict) -> Dict:
        """Filter annotation to only include MANE transcript information"""
        if not self.use_mane_transcript or not annotation:
            return annotation
        
        gene_name = annotation.get('gene', '')
        gene_detail = annotation.get('gene_detail', '')
        
        # Check if this gene has MANE transcript information
        if gene_name not in self.mane_transcripts:
            return annotation
        
        mane_info = self.mane_transcripts[gene_name]
        filtered_annotation = annotation.copy()
        
        # Filter gene_detail to only include MANE transcript information
        if gene_detail and gene_detail != 'NA':
            # Split by comma and find parts containing MANE transcript IDs
            parts = gene_detail.split(',')
            mane_parts = []
            
            for part in parts:
                part = part.strip()
                # Check if this part contains MANE transcript ID (by base ID)
                if self._is_mane_transcript_part(part, mane_info):
                    mane_parts.append(part)
            
            if mane_parts:
                filtered_annotation['gene_detail'] = ','.join(mane_parts)
                logger.info(f"Filtered gene_detail for {gene_name}: {len(mane_parts)} MANE transcript parts")
            else:
                # If no MANE transcript parts found, keep original
                logger.warning(f"No MANE transcript parts found in gene_detail for {gene_name}")
        
        # Filter AAChange based on filtered gene_detail
        if filtered_annotation.get('gene_detail') != annotation.get('gene_detail'):
            # Extract p.HGVS from the filtered gene_detail
            filtered_annotation['aa_change'] = self._extract_aa_change_from_gene_detail(
                filtered_annotation.get('gene_detail', '')
            )
        
        return filtered_annotation
    
    def _is_mane_transcript_part(self, part: str, mane_info: Dict) -> bool:
        """Check if a gene_detail part contains MANE transcript information"""
        try:
            # Extract transcript ID from the part (format: transcript_id:exon:c.position)
            if ':' in part:
                transcript_id = part.split(':')[0]
                base_transcript_id = transcript_id.split('.')[0] if '.' in transcript_id else transcript_id
                
                # Check if it matches MANE RefSeq or Ensembl transcript
                if (mane_info.get('base_refseq') == base_transcript_id or 
                    mane_info.get('base_ensembl') == base_transcript_id):
                    return True
        except Exception as e:
            logger.warning(f"Failed to check MANE transcript part: {e}")
        
        return False
    
    def _extract_aa_change_from_gene_detail(self, gene_detail: str) -> str:
        """Extract p.HGVS from gene_detail"""
        try:
            if not gene_detail or gene_detail == 'NA':
                return 'NA'
            
            # Look for p.HGVS pattern in gene_detail
            import re
            p_hgvs_pattern = r':p\.[^,]+'
            matches = re.findall(p_hgvs_pattern, gene_detail)
            
            if matches:
                # Return the first p.HGVS found
                return matches[0][1:]  # Remove the leading ':'
            else:
                return 'NA'
        except Exception as e:
            logger.warning(f"Failed to extract p.HGVS from gene_detail: {e}")
            return 'NA'
    
    def _filter_mane_transcript_details(self, annotation: Dict, mane_transcript_id: str) -> Dict:
        """Filter annotation details to only include MANE transcript information"""
        filtered_annotation = annotation.copy()
        
        # Filter gene_detail to only include MANE transcript
        gene_detail = annotation['gene_detail']
        if gene_detail != 'NA':
            # Split by comma and find parts containing the MANE transcript ID
            parts = gene_detail.split(',')
            mane_parts = [part for part in parts if mane_transcript_id in part]
            filtered_annotation['gene_detail'] = ','.join(mane_parts) if mane_parts else 'NA'
        
        # Filter exonic_function and aa_change similarly
        exonic_function = annotation['exonic_function']
        if exonic_function != 'NA':
            # For exonic_function, we need to be more careful as it might not directly contain transcript IDs
            # For now, we'll keep the original logic but this might need refinement
            pass
        
        aa_change = annotation['aa_change']
        if aa_change != 'NA':
            # For aa_change, we need to be more careful as it might not directly contain transcript IDs
            # For now, we'll keep the original logic but this might need refinement
            pass
        
        return filtered_annotation
    
    def _combine_transcript_annotations(self, annotations: List[Dict]) -> Dict:
        """Combine all transcript annotations into a single annotation"""
        if not annotations:
            return {
                'function': 'intergenic',
                'gene': 'NA',
                'gene_detail': 'NA',
                'exonic_function': 'NA',
                'aa_change': 'NA'
            }
        
        # Group annotations by gene name
        gene_groups = {}
        for annotation in annotations:
            gene_name = annotation['gene']
            if gene_name not in gene_groups:
                gene_groups[gene_name] = []
            gene_groups[gene_name].append(annotation)
        
        # For each gene, combine all transcript annotations
        combined_annotations = []
        for gene_name, gene_annotations in gene_groups.items():
            # Sort by function priority
            priority_order = ['exonic', 'intronic', 'utr5', 'utr3', 'upstream', 'downstream']
            gene_annotations.sort(key=lambda x: priority_order.index(x['function']) if x['function'] in priority_order else len(priority_order))
            
            # Combine gene details (transcript information) - use comma separation for coding_change.py
            gene_details = []
            exonic_functions = []
            aa_changes = []
            
            for annotation in gene_annotations:
                if annotation['gene_detail'] != 'NA':
                    gene_details.append(annotation['gene_detail'])
                    # Extract p.HGVS from gene_detail for AAChange column
                    p_hgvs = self._extract_p_hgvs_from_gene_detail(annotation['gene_detail'])
                    if p_hgvs:
                        aa_changes.append(p_hgvs)
                if annotation['exonic_function'] != 'NA':
                    exonic_functions.append(annotation['exonic_function'])
            
            # Create combined annotation - use comma separation for coding_change.py to process
            combined_annotation = {
                'function': gene_annotations[0]['function'],  # Use the highest priority function
                'gene': gene_name,
                'gene_detail': ','.join(gene_details) if gene_details else 'NA',
                'exonic_function': ','.join(exonic_functions) if exonic_functions else 'NA',
                'aa_change': ','.join(aa_changes) if aa_changes else 'NA'
            }
            
            # Apply MANE filtering if enabled
            if self.use_mane_transcript:
                combined_annotation = self._filter_mane_annotation(combined_annotation)
            
            combined_annotations.append(combined_annotation)
        
        # If multiple genes, return the first one (this shouldn't happen in normal cases)
        return combined_annotations[0] if combined_annotations else annotations[0]
    
    def _extract_p_hgvs_from_gene_detail(self, gene_detail: str) -> str:
        """Extract p.HGVS from gene_detail string"""
        try:
            # Split by comma to get individual transcript annotations
            transcript_annotations = [x.strip() for x in gene_detail.split(',') if x.strip()]
            
            p_hgvs_list = []
            for annotation in transcript_annotations:
                # Look for p.HGVS pattern: :p.xxx
                if ':p.' in annotation:
                    # Extract the p.HGVS part after the last ':p.'
                    p_hgvs = annotation.split(':p.')[-1]
                    if p_hgvs:
                        p_hgvs_list.append(f"p.{p_hgvs}")
                else:
                    # If no p.HGVS found, add p.?
                    p_hgvs_list.append("p.?")
            
            return ','.join(p_hgvs_list) if p_hgvs_list else "p.?"
            
        except Exception as e:
            logger.warning(f"Failed to extract p.HGVS from gene_detail: {e}")
            return "p.?"
    
    def _write_gene_annotation(self, variant: Dict, annotation: Dict, var_f, exonic_f):
        """Write the gene annotation result"""
        # Build gene name and transcript information
        gene_name = annotation['gene']
        gene_detail = annotation['gene_detail']
        
        # If GeneDetail contains transcript information, add it to the gene name
        if gene_detail and gene_detail != 'NA':
            # Format: gene_name(transcript_info)
            gene_with_detail = f"{gene_name}({gene_detail})"
        else:
            gene_with_detail = gene_name
        
        # Write to variant_function file
        var_line = f"{annotation['function']}\t{gene_with_detail}\t{variant['original_line']}\n"
        var_f.write(var_line)
        
        # If it is an exonic variant, write to the exonic_variant_function file compatible with Perl
        if annotation['function'] == 'exonic':
            # Increment line counter
            self._gene_exonic_line_counter += 1
            # Second column uses the annotation containing transcript information, which facilitates the subsequent coding_change based on mRNA modification
            annotation_with_transcript = annotation.get('gene_detail') or 'NA'
            # Concatenate standard five column variant information
            var_cols = f"{variant['chrom']}\t{variant['start']}\t{variant['end']}\t{variant['ref']}\t{variant['alt']}"
            # Perl format: lineN \t ExonicFunc \t AAChange-like(含转录本) \t Chr Start End Ref Alt
            exonic_line = f"line{self._gene_exonic_line_counter}\t{annotation['exonic_function']}\t{annotation_with_transcript}\t{var_cols}\n"
            exonic_f.write(exonic_line)
    
    def _annotate_query_by_region(self):
        """Region annotation"""
        logger.info("Start region annotation...")
        
        try:
            # Read the region database
            region_db = self._load_region_database()
            
            # Process the query file
            with open(self.queryfile, 'r', encoding='utf-8') as query_f, \
                 open(f"{self.outfile}.{self.buildver}_{self.dbtype1}", 'w', encoding='utf-8') as output_f:
                
                for line in query_f:
                    line = line.strip()
                    if not line or line.startswith('#'):
                        continue
                    
                    # Parse the variant
                    variant = self._parse_variant_line(line)
                    if not variant:
                        continue
                    
                    # Annotate the variant
                    annotation = self._annotate_variant_by_region(variant, region_db)
                    
                    # Write the result
                    output_line = f"{self.dbtype1}\t{annotation}\t{variant['original_line']}\n"
                    output_f.write(output_line)
        
        except Exception as e:
            logger.error(f"Region annotation failed: {e}")
            raise
    
    def _load_region_database(self) -> Dict:
        """Load the region database"""
        region_db = {}
        
        try:
            region_file = os.path.join(self.dbloc, f"{self.buildver}_{self.dbtype1}.txt")
            if not os.path.exists(region_file):
                logger.warning(f"Region database file does not exist: {region_file}")
                return region_db
            
            with open(region_file, 'r', encoding='utf-8') as f:
                for line in f:
                    line = line.strip()
                    if not line or line.startswith('#'):
                        continue
                    
                    parts = line.split('\t')
                    if len(parts) < 3:
                        continue
                    
                    chrom = parts[0]
                    start = int(parts[1])
                    end = int(parts[2])
                    info = parts[3] if len(parts) > 3 else ''
                    
                    if chrom not in region_db:
                        region_db[chrom] = []
                    region_db[chrom].append({
                        'start': start,
                        'end': end,
                        'info': info
                    })
        
        except Exception as e:
            logger.error(f"Failed to load the region database: {e}")
        
        return region_db
    
    def _annotate_variant_by_region(self, variant: Dict, region_db: Dict) -> str:
        """Annotate a single region variant"""
        chrom = variant['chrom']
        start = variant['start']
        end = variant['end']
        
        if chrom not in region_db:
            return 'Unknown'
        
        # Find overlapping regions
        for region in region_db[chrom]:
            if region['start'] <= end and region['end'] >= start:
                return region['info']
        
        return 'Unknown'
    
    def _filter_query(self):
        """Filter the query"""
        logger.info("Start filtering the query...")
        
        try:
            # Read the filter database
            filter_db = self._load_filter_database()
            
            # Process the query file
            with open(self.queryfile, 'r', encoding='utf-8') as query_f, \
                 open(f"{self.outfile}.{self.buildver}_{self.dbtype1}_filtered", 'w', encoding='utf-8') as filtered_f, \
                 open(f"{self.outfile}.{self.buildver}_{self.dbtype1}_dropped", 'w', encoding='utf-8') as dropped_f:
                
                for line in query_f:
                    line = line.strip()
                    if not line or line.startswith('#'):
                        continue
                    
                    # Parse the variant
                    variant = self._parse_variant_line(line)
                    if not variant:
                        continue
                    
                    # Filter the variant
                    is_filtered = self._filter_variant(variant, filter_db)
                    
                    # Write the result
                    if is_filtered:
                        # Get the annotation information for matched variants数据库
                        annotation_info = self._get_filter_annotation_info(variant, filter_db)
                        
                        if self.otherinfo:
                            # When using otherinfo, output the complete annotation information
                            filtered_f.write(f"{self.dbtype1}\t{annotation_info}\t{variant['original_line']}\n")
                        else:
                            # When not using otherinfo, only output the matching status
                            filtered_f.write(f"{self.dbtype1}\tFiltered\t{variant['original_line']}\n")
                    else:
                        # For variants that don't match the database, use nastring (default ".")
                        if self.otherinfo:
                            dropped_f.write(f"{self.dbtype1}\t.\t{variant['original_line']}\n")
                        else:
                            dropped_f.write(f"{self.dbtype1}\t.\t{variant['original_line']}\n")
        
        except Exception as e:
            logger.error(f"Filtering the query failed: {e}")
            raise
    
    def _load_filter_database(self) -> Dict:
        """Load the filter database with optimized memory usage.
        For huge files (e.g., 30G gnomad), prefer Tabix index if available to avoid full load.
        Returns either:
          - dict for small files (legacy behavior)
          - a special handler dict with key '__tabix__' carrying a TabixFile and mode='tabix'
        """
        filter_file_txt = os.path.join(self.dbloc, f"{self.buildver}_{self.dbtype1}.txt")
        filter_file_gz = filter_file_txt + ".gz"
        tabix_tbi = filter_file_gz + ".tbi"
        
        # 检查文件大小，决定加载策略
        file_size_mb = 0
        if os.path.exists(filter_file_txt):
            file_size_mb = os.path.getsize(filter_file_txt) / (1024 * 1024)
        elif os.path.exists(filter_file_gz):
            file_size_mb = os.path.getsize(filter_file_gz) / (1024 * 1024)
        
        # 对于大于100MB的文件，优先使用Tabix索引
        if file_size_mb > 100 and pysam and os.path.exists(filter_file_gz) and os.path.exists(tabix_tbi):
            try:
                logger.info(f"Using Tabix-indexed database for fast access: {filter_file_gz} ({file_size_mb:.1f}MB)")
                tbx = pysam.TabixFile(filter_file_gz)
                # Return a lightweight handler
                return {
                    '__mode__': 'tabix',
                    '__tabix__': tbx,
                    '__path__': filter_file_gz
                }
            except Exception as e:
                logger.warning(f"Failed to open Tabix database, falling back to full load: {e}")
                # fallthrough to full load
        elif file_size_mb > 100:
            logger.info(f"Large database detected ({file_size_mb:.1f}MB), but no Tabix index found. Consider building index for better performance.")
        
        # Fall back to loading into memory (for small/medium files)
        filter_db: Dict[str, str] = {}
        try:
            # Check if compressed file exists when txt file doesn't exist
            if not os.path.exists(filter_file_txt):
                if os.path.exists(filter_file_gz):
                    logger.info(f"Using compressed database file: {filter_file_gz}")
                    # Use gzip to read compressed file
                    import gzip
                    with gzip.open(filter_file_gz, 'rt', encoding='utf-8') as f:
                        for line in f:
                            line = line.strip()
                            if not line or line.startswith('#'):
                                continue
                            parts = line.split('\t')
                            if len(parts) < 5:
                                continue
                            chrom = parts[0]
                            if parts[1] == parts[2]:
                                pos = parts[1]
                                ref = parts[3]
                                alt = parts[4]
                            else:
                                pos = parts[1]
                                ref = parts[3]
                                alt = parts[4]
                            key = f"{chrom}:{pos}:{ref}:{alt}"
                            value = '\t'.join(parts[5:]) if len(parts) > 5 else ''
                            filter_db[key] = value
                    return filter_db
                else:
                    logger.warning(f"Filter database file does not exist: {filter_file_txt}")
                    return filter_db
            
            with open(filter_file_txt, 'r', encoding='utf-8') as f:
                for line in f:
                    line = line.strip()
                    if not line or line.startswith('#'):
                        continue
                    parts = line.split('\t')
                    if len(parts) < 5:
                        continue
                    chrom = parts[0]
                    if parts[1] == parts[2]:
                        pos = parts[1]
                        ref = parts[3]
                        alt = parts[4]
                    else:
                        pos = parts[1]
                        ref = parts[3]
                        alt = parts[4]
                    key = f"{chrom}:{pos}:{ref}:{alt}"
                    value = '\t'.join(parts[5:]) if len(parts) > 5 else ''
                    filter_db[key] = value
        except Exception as e:
            logger.error(f"Failed to load the filter database: {e}")
        return filter_db
    
    def _filter_variant(self, variant: Dict, filter_db: Dict) -> bool:
        """Filter a single variant.
        Supports two modes:
          - dict exact-lookup (legacy)
          - Tabix on-the-fly region query (for huge .gz+.tbi)
        """
        if filter_db and filter_db.get('__mode__') == 'tabix' and '__tabix__' in filter_db:
            # Tabix expects 1-based inclusive coordinates, query by exact position window
            chrom = str(variant['chrom'])
            pos = int(variant['start'])
            ref = variant['ref']
            alt = variant['alt']
            tbx = filter_db['__tabix__']
            try:
                # Tabix chromosome naming sometimes lacks 'chr'
                chroms = [chrom]
                if chrom.startswith('chr'):
                    chroms.append(chrom[3:])
                else:
                    chroms.append('chr' + chrom)
                for c in chroms:
                    try:
                        for rec in tbx.fetch(c, max(0, pos-1), pos):
                            parts = rec.strip().split('\t')
                            if len(parts) < 5:
                                continue
                            db_chrom = parts[0]
                            # Handle both interval and POS-POS line
                            db_pos = parts[1]
                            db_ref = parts[3]
                            db_alt = parts[4]
                            # Match by exact key
                            if (db_chrom == chrom or db_chrom == chrom.replace('chr','') or ('chr'+db_chrom)==chrom) \
                               and str(db_pos) == str(pos) and db_ref == ref and db_alt == alt:
                                return True
                    except ValueError:
                        continue
            except Exception as e:
                logger.warning(f"Tabix query failed for {chrom}:{pos}: {e}")
            return False
        # Fallback: dict lookup
        key = f"{variant['chrom']}:{variant['start']}:{variant['ref']}:{variant['alt']}"
        return key in filter_db

    def _get_filter_annotation_info(self, variant: Dict, filter_db: Dict) -> str:
        """Fetch annotation info string for a filtered variant.
        Works for both dict mode and Tabix mode.
        """
        if filter_db and filter_db.get('__mode__') == 'tabix' and '__tabix__' in filter_db:
            chrom = str(variant['chrom'])
            pos = int(variant['start'])
            ref = variant['ref']
            alt = variant['alt']
            tbx = filter_db['__tabix__']
            try:
                chroms = [chrom]
                if chrom.startswith('chr'):
                    chroms.append(chrom[3:])
                else:
                    chroms.append('chr' + chrom)
                for c in chroms:
                    try:
                        for rec in tbx.fetch(c, max(0, pos-1), pos):
                            parts = rec.strip().split('\t')
                            if len(parts) < 5:
                                continue
                            db_chrom = parts[0]
                            db_pos = parts[1]
                            db_ref = parts[3]
                            db_alt = parts[4]
                            if (db_chrom == chrom or db_chrom == chrom.replace('chr','') or ('chr'+db_chrom)==chrom) \
                               and str(db_pos) == str(pos) and db_ref == ref and db_alt == alt:
                                return '\t'.join(parts[5:]) if len(parts) > 5 else ''
                    except ValueError:
                        continue
            except Exception as e:
                logger.warning(f"Tabix query for annotation failed for {chrom}:{pos}: {e}")
            return ''
        # dict mode
        key = f"{variant['chrom']}:{variant['start']}:{variant['ref']}:{variant['alt']}"
        return filter_db.get(key, '')
    
    def _annotate_query_by_gene_thread(self, var_output: str, exonic_output: str, invalid_output: str, 
                                      start_line: int, end_line: int, thread_id: int, gene_db):
        """Gene annotation thread function"""
        logger.info(f"Thread {thread_id}: processing lines {start_line} to {end_line}")
        
        try:
            # Process the specified line range
            with open(self.queryfile, 'r', encoding='utf-8') as query_f, \
                 open(var_output, 'w', encoding='utf-8') as var_f, \
                 open(exonic_output, 'w', encoding='utf-8') as exonic_f, \
                 open(invalid_output, 'w', encoding='utf-8') as invalid_f:
                
                line_count = 0
                for line in query_f:
                    line_count += 1
                    
                    if line_count < start_line:
                        continue
                    if line_count > end_line:
                        break
                    
                    line = line.strip()
                    if not line or line.startswith('#'):
                        continue
                    
                    # Parse the variant
                    variant = self._parse_variant_line(line)
                    if not variant:
                        invalid_f.write(f"line{line_count}\tInvalid format\t{line}\n")
                        continue
                    
                    # Annotate the variant
                    annotation = self._annotate_variant_by_gene(variant, gene_db)
                    
                    # Write the result
                    self._write_gene_annotation(variant, annotation, var_f, exonic_f)
        
        except Exception as e:
            logger.error(f"Thread {thread_id} failed: {e}")
    
    def _annotate_query_by_region_thread(self, output_file: str, invalid_output: str, 
                                        start_line: int, end_line: int, thread_id: int, region_db):
        """Region annotation thread function"""
        logger.info(f"Thread {thread_id}: processing lines {start_line} to {end_line}")
        
        try:
            # Process the specified line range
            with open(self.queryfile, 'r', encoding='utf-8') as query_f, \
                 open(output_file, 'w', encoding='utf-8') as output_f, \
                 open(invalid_output, 'w', encoding='utf-8') as invalid_f:
                
                line_count = 0
                for line in query_f:
                    line_count += 1
                    
                    if line_count < start_line:
                        continue
                    if line_count > end_line:
                        break
                    
                    line = line.strip()
                    if not line or line.startswith('#'):
                        continue
                    
                    # Parse the variant
                    variant = self._parse_variant_line(line)
                    if not variant:
                        invalid_f.write(f"line{line_count}\tInvalid format\t{line}\n")
                        continue
                    
                    # Annotate the variant
                    annotation = self._annotate_variant_by_region(variant, region_db)
                    
                    # Write the result
                    output_line = f"{self.dbtype1}\t{annotation}\t{variant['original_line']}\n"
                    output_f.write(output_line)
        
        except Exception as e:
            logger.error(f"Thread {thread_id} failed: {e}")
    
    def _filter_query_thread(self, filtered_output: str, dropped_output: str, invalid_output: str,
                           start_line: int, end_line: int, thread_id: int, filter_db):
        """Filter query thread function"""
        logger.info(f"Thread {thread_id}: processing lines {start_line} to {end_line}")
        
        try:
            # For Tabix mode, create a separate TabixFile for this thread to avoid concurrency issues
            if filter_db and filter_db.get('__mode__') == 'tabix':
                try:
                    import pysam
                    filter_db = {
                        '__mode__': 'tabix',
                        '__tabix__': pysam.TabixFile(filter_db['__path__']),
                        '__path__': filter_db['__path__']
                    }
                except Exception as e:
                    logger.error(f"Thread {thread_id} failed to create TabixFile: {e}")
                    return
            
            # Process the specified line range
            with open(self.queryfile, 'r', encoding='utf-8') as query_f, \
                 open(filtered_output, 'w', encoding='utf-8') as filtered_f, \
                 open(dropped_output, 'w', encoding='utf-8') as dropped_f, \
                 open(invalid_output, 'w', encoding='utf-8') as invalid_f:
                
                line_count = 0
                for line in query_f:
                    line_count += 1
                    
                    if line_count < start_line:
                        continue
                    if line_count > end_line:
                        break
                    
                    line = line.strip()
                    if not line or line.startswith('#'):
                        continue
                    
                    # Parse the variant
                    variant = self._parse_variant_line(line)
                    if not variant:
                        invalid_f.write(f"line{line_count}\tInvalid format\t{line}\n")
                        continue
                    
                    # Filter the variant
                    is_filtered = self._filter_variant(variant, filter_db)
                    
                    # Write the result
                    if is_filtered:
                        # Get the annotation information
                        annotation_info = self._get_filter_annotation_info(variant, filter_db)
                        
                        if self.otherinfo:
                            # When using otherinfo, output the complete annotation information
                            filtered_f.write(f"{self.dbtype1}\t{annotation_info}\t{variant['original_line']}\n")
                        else:
                            # When not using otherinfo, only output the matching status
                            filtered_f.write(f"{self.dbtype1}\tFiltered\t{variant['original_line']}\n")
                    else:
                        # For variants that don't match the database, use nastring (default ".")
                        if self.otherinfo:
                            dropped_f.write(f"{self.dbtype1}\t.\t{variant['original_line']}\n")
                        else:
                            dropped_f.write(f"{self.dbtype1}\t.\t{variant['original_line']}\n")
        
        except Exception as e:
            logger.error(f"Thread {thread_id} failed: {e}")
    
    def _merge_thread_results(self):
        """Merge thread results"""
        logger.info("Merge thread results...")
        
        try:
            if self.geneanno:
                # Merge the variant_function file
                with open(f"{self.outfile}.variant_function", 'w', encoding='utf-8') as merged_f:
                    for i in range(self.thread):
                        thread_file = f"{self.outfile}.variant_function.{i}"
                        if os.path.exists(thread_file):
                            with open(thread_file, 'r', encoding='utf-8') as f:
                                merged_f.write(f.read())
                            os.remove(thread_file)
                
                # Merge the exonic_variant_function file
                with open(f"{self.outfile}.exonic_variant_function", 'w', encoding='utf-8') as merged_f:
                    for i in range(self.thread):
                        thread_file = f"{self.outfile}.exonic_variant_function.{i}"
                        if os.path.exists(thread_file):
                            with open(thread_file, 'r', encoding='utf-8') as f:
                                merged_f.write(f.read())
                            os.remove(thread_file)
            
            elif self.regionanno:
                # Merge the region annotation file
                with open(f"{self.outfile}.{self.buildver}_{self.dbtype1}", 'w', encoding='utf-8') as merged_f:
                    for i in range(self.thread):
                        thread_file = f"{self.outfile}.{self.buildver}_{self.dbtype1}.{i}"
                        if os.path.exists(thread_file):
                            with open(thread_file, 'r', encoding='utf-8') as f:
                                merged_f.write(f.read())
                            os.remove(thread_file)
            
            elif self.filter:
                # Merge the filter file
                with open(f"{self.outfile}.{self.buildver}_{self.dbtype1}_filtered", 'w', encoding='utf-8') as merged_f:
                    for i in range(self.thread):
                        thread_file = f"{self.outfile}.{self.buildver}_{self.dbtype1}_filtered.{i}"
                        if os.path.exists(thread_file):
                            with open(thread_file, 'r', encoding='utf-8') as f:
                                merged_f.write(f.read())
                            os.remove(thread_file)
                
                with open(f"{self.outfile}.{self.buildver}_{self.dbtype1}_dropped", 'w', encoding='utf-8') as merged_f:
                    for i in range(self.thread):
                        thread_file = f"{self.outfile}.{self.buildver}_{self.dbtype1}_dropped.{i}"
                        if os.path.exists(thread_file):
                            with open(thread_file, 'r', encoding='utf-8') as f:
                                merged_f.write(f.read())
                            os.remove(thread_file)
        
        except Exception as e:
            logger.error(f"Failed to merge thread results: {e}")
    
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
    
    def _get_standard_transcript_id(self, gene_name: str, transcript_id: str, protocol: str) -> str:
        """Get the standard transcript ID"""
        try:
            # Only use MANE transcript filtering if use_mane_transcript is True
            if not self.use_mane_transcript:
                return transcript_id
            
            # Process the gene name based on the protocol type
            if protocol == 'refGene':
                # refGene: use the gene name directly
                gene_id = gene_name
            elif protocol == 'ensGene':
                # ensGene: add 'Ensembl:' prefix
                gene_id = f"Ensembl:{gene_name}"
            elif protocol == 'knownGene':
                # knownGene: use the gene name directly
                gene_id = gene_name
            else:
                # Other protocols: use the gene name directly
                gene_id = gene_name
            
            # Find MANE transcript (using the preloaded mapping)
            if gene_id in self.mane_transcripts:
                mane_info = self.mane_transcripts[gene_id]
                
                # Get the base transcript ID (without version) for matching
                base_transcript_id = transcript_id.split('.')[0] if '.' in transcript_id else transcript_id
                
                # Check if current transcript matches MANE transcript (by base ID)
                if protocol == 'refGene' and mane_info.get('base_refseq') == base_transcript_id:
                    logger.info(f"Using MANE RefSeq transcript {mane_info['refseq']} for gene {gene_id}")
                    return mane_info['refseq']
                elif protocol == 'ensGene' and mane_info.get('base_ensembl') == base_transcript_id:
                    logger.info(f"Using MANE Ensembl transcript {mane_info['ensembl']} for gene {gene_id}")
                    return mane_info['ensembl']
                else:
                    # If no exact match, return the original transcript ID
                    logger.info(f"No matching MANE transcript found for gene {gene_id}, using original transcript {transcript_id}")
                    return transcript_id
            else:
                # If the MANE transcript is not found, return the original transcript ID
                logger.info(f"No MANE transcript found for gene {gene_id}, using original transcript {transcript_id}")
                return transcript_id
                
        except Exception as e:
            logger.warning(f"Failed to get the standard transcript ID: {e}")
            return transcript_id
    
    def _load_mane_transcripts_from_file(self, mane_file: str) -> Dict[str, str]:
        """Load MANE transcript mappings from file"""
        mane_transcripts = {}
        
        try:
            with open(mane_file, 'r', encoding='utf-8') as f:
                for line in f:
                    line = line.strip()
                    if not line or line.startswith('#'):
                        continue
                    
                    parts = line.split('\t')
                    if len(parts) >= 3:
                        # 三列格式：gene_id, refseq_id, ensembl_id
                        gene_id = parts[0]
                        refseq_id = parts[1]
                        ensembl_id = parts[2]
                        
                        if gene_id and refseq_id:
                            # Store both RefSeq and Ensembl transcript IDs
                            base_refseq = refseq_id.split('.')[0] if '.' in refseq_id else refseq_id
                            base_ensembl = ensembl_id.split('.')[0] if ensembl_id and '.' in ensembl_id else ensembl_id
                            
                            mane_transcripts[gene_id] = {
                                'refseq': refseq_id,
                                'ensembl': ensembl_id,
                                'base_refseq': base_refseq,
                                'base_ensembl': base_ensembl
                            }
                    elif len(parts) >= 9:
                        # GTF格式：chr, source, feature, start, end, score, strand, frame, attributes
                        attributes = parts[8]
                        
                        # Parse the attribute field
                        gene_id = None
                        transcript_id = None
                        ensembl_transcript_id = None
                        
                        for attr in attributes.split(';'):
                            attr = attr.strip()
                            if attr.startswith('gene_id'):
                                gene_id = attr.split('"')[1] if '"' in attr else attr.split()[1]
                            elif attr.startswith('transcript_id'):
                                transcript_id = attr.split('"')[1] if '"' in attr else attr.split()[1]
                            elif attr.startswith('db_xref') and 'Ensembl:' in attr:
                                # Extract Ensembl transcript ID
                                ensembl_part = attr.split('Ensembl:')[1].split('"')[0] if '"' in attr else attr.split('Ensembl:')[1]
                                ensembl_transcript_id = ensembl_part
                        
                        if gene_id and transcript_id:
                            # Store both RefSeq and Ensembl transcript IDs
                            # Use the base transcript ID (without version) as key for matching
                            base_transcript_id = transcript_id.split('.')[0] if '.' in transcript_id else transcript_id
                            base_ensembl_id = ensembl_transcript_id.split('.')[0] if ensembl_transcript_id and '.' in ensembl_transcript_id else ensembl_transcript_id
                            
                            mane_transcripts[gene_id] = {
                                'refseq': transcript_id,
                                'ensembl': ensembl_transcript_id,
                                'base_refseq': base_transcript_id,
                                'base_ensembl': base_ensembl_id
                            }
            
            logger.info(f"Loaded {len(mane_transcripts)} MANE transcript mappings")
        except Exception as e:
            logger.error(f"Failed to load MANE transcript mappings from file: {e}")
        
        return mane_transcripts

    def _ensembl_rest_base(self) -> str:
        # Select the REST base URL based on buildver
        if str(self.buildver).lower() in ("hg19", "grch37"):
            return "https://grch37.rest.ensembl.org"
        return "https://rest.ensembl.org"

    def _fetch_sequence_region(self, chrom: str, start: int, end: int) -> Optional[str]:
        base = self._ensembl_rest_base()
        # Remove the 'chr' prefix
        if chrom.lower().startswith('chr'):
            chrom = chrom[3:]
        url = f"{base}/sequence/region/human/{chrom}:{start}..{end}:1?content-type=text/plain"
        req = request.Request(url, headers={"Content-Type": "text/plain", "User-Agent": "matchvar-python"})
        try:
            with request.urlopen(req, timeout=8) as resp:
                data = resp.read().decode('utf-8').strip().upper()
                return data if data else None
        except error.URLError:
            return None
        except Exception:
            return None

    def _is_intronic_dup_by_remote(self, chrom: str, pos: int, alt_seq: str, strand: str) -> bool:
        """
        Identify intronic insertions as left-duplicates (dup) through Ensembl REST.
        Take the sequence in the window [pos - window, pos - 1], compare the last k=len(alt) bases with alt;
        Use reverse complement for negative strand.
        """
        try:
            k = len(alt_seq)
            if k == 0:
                return False
            if strand == '+':
                # Positive strand: take the left window, compare the last k bases with alt
                start = max(1, pos - self.intronic_dup_window)
                end = pos - 1
                if end < start:
                    return False
                seq = self._fetch_sequence_region(chrom, start, end)
                if not seq or len(seq) < k:
                    return False
                return seq[-k:] == alt_seq
            else:
                # Negative strand: duplicate sequence on the right side of the genome (upstream of the transcript)
                start = pos + 1
                end = pos + self.intronic_dup_window
                seq = self._fetch_sequence_region(chrom, start, end)
                if not seq or len(seq) < k:
                    return False
                # Negative strand Alt has been complemented in the upstream logic, compare the genome direction sequence directly
                return seq[:k] == self._reverse_complement(alt_seq)
        except Exception:
            return False

    def _generate_utr_intronic_hgvs(self, variant_start, ref, alt, strand, exon_starts, exon_ends, 
                                   cds_start, cds_end, standard_transcript_id):
        """Generate the HGVS format for UTR intronic variants""" 
        # For UTR intronic variants, select the reference exon based on the position of the variant relative to the CDS
        # Find the CDS start and end exons
        cds_start_exon_idx = None
        cds_end_exon_idx = None
        cds_start_exon_num = None
        cds_end_exon_num = None
        
        for i, (exon_start, exon_end) in enumerate(zip(exon_starts, exon_ends)):
            if strand == '-':
                transcript_exon_num = len(exon_starts) - i
            else:
                transcript_exon_num = i + 1
            
            # Check if the CDS start position is included
            if exon_start <= cds_start <= exon_end:
                cds_start_exon_idx = i
                cds_start_exon_num = transcript_exon_num
            
            # Check if the CDS end position is included
            if exon_start <= cds_end <= exon_end:
                cds_end_exon_idx = i
                cds_end_exon_num = transcript_exon_num
        
        # For negative strand genes, swap the exon coordinates (start and end positions)
        if strand == '-':
            # Negative strand gene: exon coordinates from largest to smallest
            # Swap the start and end positions in the refGene database
            real_exon_starts = exon_ends  # The end position in the database is the real start position
            real_exon_ends = exon_starts   # The start position in the database is the real end position
        else:
            # Positive strand gene: use the original coordinates
            real_exon_starts = exon_starts
            real_exon_ends = exon_ends
        
        # Find the exon closest to the variant site as the reference
        min_distance = float('inf')
        nearest_exon_idx = 0
        nearest_exon_num = 1
        nearest_exon_start = None
        nearest_exon_end = None
        
        for i, (exon_start, exon_end) in enumerate(zip(real_exon_starts, real_exon_ends)):
            if strand == '-':
                transcript_exon_num = len(exon_starts) - i
            else:
                transcript_exon_num = i + 1
            
            # Calculate the distance to the exon
            if variant_start < exon_start:
                distance = exon_start - variant_start
            elif variant_start > exon_end:
                distance = variant_start - exon_end
            else:
                distance = 0
            
            if distance < min_distance:
                min_distance = distance
                nearest_exon_idx = i
                nearest_exon_num = transcript_exon_num
                nearest_exon_start = exon_start
                nearest_exon_end = exon_end
        
        # Determine the reference point: select the start or end position of the nearest exon based on the distance
        reference_point = None
        reference_exon_num = None
        
        # Calculate the distance to the start and end positions of the nearest exon
        distance_to_start = abs(variant_start - nearest_exon_start)
        distance_to_end = abs(variant_start - nearest_exon_end)
        
        # For negative strand genes, we need to consider the transcript direction
        # The transcript direction is 5' to 3', so we need to determine which boundary
        # is the "upstream" boundary in transcript coordinates
        if strand == '-':
            # For negative strand genes:
            # - nearest_exon_start (in real coordinates) corresponds to the 3' end of the exon in transcript
            # - nearest_exon_end (in real coordinates) corresponds to the 5' end of the exon in transcript
            # We should use the boundary that is upstream in transcript coordinates
            if distance_to_start <= distance_to_end:
                # Use the 3' end (nearest_exon_start in real coordinates) as reference
                reference_point = nearest_exon_start
                reference_exon_num = nearest_exon_num
                # In transcript coordinates, if variant is after the 3' end, it's downstream (-)
                if variant_start > reference_point:
                    sign = '-'
                else:
                    sign = '+'
            else:
                # Use the 5' end (nearest_exon_end in real coordinates) as reference
                reference_point = nearest_exon_end
                reference_exon_num = nearest_exon_num
                # In transcript coordinates, if variant is before the 5' end, it's upstream (+)
                if variant_start < reference_point:
                    sign = '+'
                else:
                    sign = '-'
        else:
            # For positive strand genes: use the closer boundary
            if distance_to_start <= distance_to_end:
                # Use the start position of the exon as the reference point
                reference_point = nearest_exon_start
                reference_exon_num = nearest_exon_num
                # If the site is after the reference point, the sign is "-"
                if variant_start > reference_point:
                    sign = '-'
                else:
                    sign = '+'
            else:
                # Use the end position of the exon as the reference point
                reference_point = nearest_exon_end
                reference_exon_num = nearest_exon_num
                # If the site is before the reference point, the sign is "+"
                if variant_start < reference_point:
                    sign = '+'
                else:
                    sign = '-'
        

        
        # Use the determined reference point to calculate the offset
        if reference_point is not None:
            
            # Calculate the offset to the reference point
            offset = abs(variant_start - reference_point)
            
            # The sign has already been determined above
            # For negative strand genes, the sign needs to be adjusted based on the transcript direction
            if strand == '-':
                # Special handling for negative strand genes: adjust the sign based on the transcript direction
                if sign == '+':
                    # On the transcript direction, it is positive
                    sign = '+'
                else:
                    # On the transcript direction, it is negative
                    sign = '-'
            
            reference_exon = reference_exon_num
            
            # Process the base sequence
            ref_seq = ref
            alt_seq = alt
            if strand == '-':
                ref_seq = self._reverse_complement(ref) if ref else ''
                alt_seq = self._reverse_complement(alt) if alt else ''
            
            # Generate the HGVS format
            if len(ref_seq) == 1 and len(alt_seq) == 1:
                hgvs = f"{standard_transcript_id}:c.{reference_exon}{sign}{offset}{ref_seq}>{alt_seq}:p.?"
            else:
                hgvs = f"{standard_transcript_id}:c.{reference_exon}{sign}{offset}delins{alt_seq}:p.?"
        else:
            hgvs = f"{standard_transcript_id}:intronic:p.?"
        
        return hgvs

def main():
    """Main function"""
    examples = (
        "示例:\n"
        "1) 基因注释（RefSeq）：\n"
        "   python utils/matchvar/annotate_variation.py -geneanno -buildver hg19 \\\n+        -dbtype refGene -outfile out input.mvinput resources/humandb\n\n"
        "2) 区域注释（cytoBand）：\n"
        "   python utils/matchvar/annotate_variation.py -regionanno -buildver hg19 \\\n+        -dbtype cytoBand -outfile out input.mvinput resources/humandb\n\n"
        "3) 过滤（ClinVar）：\n"
        "   python utils/matchvar/annotate_variation.py -filter -buildver hg19 \\\n+        -dbtype clinvar -outfile out input.mvinput resources/humandb -otherinfo\n\n"
        "4) 使用MANE转录本辅助（基因注释）：\n"
        "   python utils/matchvar/annotate_variation.py -geneanno -buildver hg19 \\\n+        -dbtype refGene -outfile out input.mvinput resources/humandb -mane_file resources/humandb/mane_transcript.txt\n"
    )
    parser = argparse.ArgumentParser(
        description='MATCHVAR variant annotation tool',
        epilog=examples,
        formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument('queryfile', help='Query file')
    parser.add_argument('dbloc', help='Database location')
    parser.add_argument('-outfile', help='Output file')
    parser.add_argument('-dbtype', help='Database type')
    parser.add_argument('-geneanno', action='store_true', help='Gene annotation')
    parser.add_argument('-regionanno', action='store_true', help='Region annotation')
    parser.add_argument('-filter', action='store_true', help='Filter')
    parser.add_argument('-buildver', default='hg19', help='Genome version')
    parser.add_argument('-thread', type=int, help='Thread number')
    parser.add_argument('-maxgenethread', type=int, default=4, help='Maximum gene thread number')
    parser.add_argument('-mingenelinecount', type=int, default=1000000, help='Minimum gene line count')
    parser.add_argument('-exonsort', action='store_true', help='Exon sorting')
    parser.add_argument('-nofirstcodondel', action='store_true', help='Do not delete the first codon')
    parser.add_argument('-otherinfo', action='store_true', help='Other information')
    parser.add_argument('-nastring', default='.', help='NA string')
    parser.add_argument('-splicing_threshold', type=int, default=2, help='Splicing variant and exon/intron boundary distance threshold')
    parser.add_argument('-indel_splicing_threshold', type=int, help='Indel splicing variant threshold (default equals splicing_threshold)')
    parser.add_argument('-sift_threshold', type=float, default=0.05, help='SIFT score threshold')
    parser.add_argument('-score_threshold', type=float, help='SIFT score threshold')
    parser.add_argument('-reverse', action='store_true', help='Reverse strand')
    parser.add_argument('-rawscore', action='store_true', help='Output raw SIFT score')
    
    # Add new important parameters
    parser.add_argument('-verbose', '-v', action='store_true', help='Verbose output')
    parser.add_argument('-man', '-m', action='store_true', help='Display manual')
    parser.add_argument('-separate', action='store_true', help='Separate output')
    parser.add_argument('-batchsize', type=str, default='5000000', help='Batch size')
    parser.add_argument('-genomebinsize', type=str, help='Genome bin size')
    parser.add_argument('-neargene', type=int, default=1000, help='Gene distance threshold')
    parser.add_argument('-expandbin', type=int, help='Expand bin number')
    parser.add_argument('-maf_threshold', type=float, default=0, help='MAF threshold')
    parser.add_argument('-normscore_threshold', type=int, help='Normalized score threshold')
    parser.add_argument('-indexfilter_threshold', type=float, default=0.9, help='Index filter threshold')
    parser.add_argument('-chromosome', type=str, help='Chromosome restriction')
    parser.add_argument('-gff3dbfile', type=str, help='GFF3 database file')
    parser.add_argument('-genericdbfile', type=str, help='Generic database file')
    parser.add_argument('-vcfdbfile', type=str, help='VCF database file')
    parser.add_argument('-precedence', type=str, help='Precedence')
    parser.add_argument('-hgvs', action='store_true', help='HGVS format')
    parser.add_argument('-exonicsplicing', action='store_true', help='Exonic splicing')
    parser.add_argument('-downdb', action='store_true', help='Download database')
    parser.add_argument('-time', action='store_true', help='Time information')
    parser.add_argument('-wget', action='store_true', default=True, help='Use wget')
    parser.add_argument('-comment', action='store_true', help='Comment information')
    parser.add_argument('-transcript_function', action='store_true', help='Transcript function')
    parser.add_argument('-avcolumn', type=str, help='AV column')
    parser.add_argument('-bedfile', type=str, help='BED file')
    parser.add_argument('-seq_padding', type=int, help='Sequence padding')
    parser.add_argument('-infoasscore', action='store_true', help='Information as score')
    parser.add_argument('-firstcodondel', action='store_true', default=True, help='Delete the first codon')
    parser.add_argument('-aamatrixfile', type=str, help='Amino acid matrix file')
    parser.add_argument('-gff3attribute', action='store_true', help='GFF3 attribute')
    parser.add_argument('-infosep', action='store_true', help='Information separator')
    parser.add_argument('-dbm', action='store_true', help='DBM mode')
    parser.add_argument('-idasscore', action='store_true', help='ID as score')
    parser.add_argument('-minqueryfrac', type=float, default=0, help='Minimum query fraction')
    parser.add_argument('-scorecolumn', type=int, help='Score column')
    parser.add_argument('-poscolumn', type=str, help='Position column')
    parser.add_argument('-webfrom', type=str, help='Web source')
    parser.add_argument('-colsWanted', type=str, help='Desired columns')
    parser.add_argument('-zerostart', action='store_true', help='Zero start')
    parser.add_argument('-memfree', type=int, help='Available memory')
    parser.add_argument('-memtotal', type=int, help='Total memory')
    parser.add_argument('-mane_file', type=str, help='MANE transcript mapping file')
    parser.add_argument('-use_mane_transcript', action='store_true', help='Use MANE transcript filtering')
    # Intronic dup online recognition
    parser.add_argument('--intronic_dup_remote', action='store_true', help='Use Ensembl REST to recognize intronic dup')
    parser.add_argument('--intronic_dup_window', type=int, default=50, help='Intronic dup recognition left window size')
    
    args = parser.parse_args()
    
    # Check which parameters are explicitly set
    sift_threshold_explicitly_set = 'sift_threshold' in sys.argv
    indexfilter_threshold_explicitly_set = 'indexfilter_threshold' in sys.argv
    maf_threshold_explicitly_set = 'maf_threshold' in sys.argv
    minqueryfrac_explicitly_set = 'minqueryfrac' in sys.argv
    wget_explicitly_set = 'wget' in sys.argv
    precedence_explicitly_set = 'precedence' in sys.argv
    
    # Create the annotator
    annotator = AnnotateVariation(
        queryfile=args.queryfile,
        dbloc=args.dbloc,
        outfile=args.outfile,
        dbtype=args.dbtype,
        geneanno=args.geneanno,
        regionanno=args.regionanno,
        filter=args.filter,
        buildver=args.buildver,
        thread=args.thread,
        maxgenethread=args.maxgenethread,
        mingenelinecount=args.mingenelinecount,
        splicing_threshold=args.splicing_threshold,
        indel_splicing_threshold=args.indel_splicing_threshold,
        otherinfo=args.otherinfo,
        sift_threshold=args.sift_threshold,
        score_threshold=args.score_threshold,
        reverse=args.reverse,
        rawscore=args.rawscore,
        verbose=args.verbose,
        batchsize=args.batchsize,
        genomebinsize=args.genomebinsize,
        neargene=args.neargene,
        expandbin=args.expandbin,
        maf_threshold=args.maf_threshold,
        normscore_threshold=args.normscore_threshold,
        indexfilter_threshold=args.indexfilter_threshold,
        chromosome=args.chromosome,
        gff3dbfile=args.gff3dbfile,
        genericdbfile=args.genericdbfile,
        vcfdbfile=args.vcfdbfile,
        precedence=args.precedence,
        hgvs=args.hgvs,
        exonicsplicing=args.exonicsplicing,
        separate=args.separate,
        downdb=args.downdb,
        time=args.time,
        wget=args.wget,
        comment=args.comment,
        transcript_function=args.transcript_function,
        avcolumn=args.avcolumn,
        bedfile=args.bedfile,
        seq_padding=args.seq_padding,
        infoasscore=args.infoasscore,
        firstcodondel=args.firstcodondel,
        aamatrixfile=args.aamatrixfile,
        gff3attribute=args.gff3attribute,
        infosep=args.infosep,
        dbm=args.dbm,
        idasscore=args.idasscore,
        minqueryfrac=args.minqueryfrac,
        scorecolumn=args.scorecolumn,
        poscolumn=args.poscolumn,
        webfrom=args.webfrom,
        colsWanted=args.colsWanted,
        zerostart=args.zerostart,
        memfree=args.memfree,
        memtotal=args.memtotal,
        mane_file=args.mane_file,
        use_mane_transcript=args.use_mane_transcript,
        intronic_dup_remote=args.intronic_dup_remote,
        intronic_dup_window=args.intronic_dup_window,
        _sift_threshold_explicitly_set=sift_threshold_explicitly_set,
        _indexfilter_threshold_explicitly_set=indexfilter_threshold_explicitly_set,
        _maf_threshold_explicitly_set=maf_threshold_explicitly_set,
        _minqueryfrac_explicitly_set=minqueryfrac_explicitly_set,
        _wget_explicitly_set=wget_explicitly_set,
        _precedence_explicitly_set=precedence_explicitly_set
    )
    
    # When no arguments are provided, show help with examples
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(0)
    
    # Run the annotation
    try:
        annotator.run_annotation()
        logger.info("Annotation completed successfully")
    except Exception as e:
        logger.error(f"Error during annotation: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main() 