import os, re, itertools, gzip, argparse
from typing import List, Optional, Dict, Tuple
from Bio import SeqIO
from Bio.Seq import Seq, reverse_complement
from pyfaidx import Fasta

class ExonExtractor:
    """
    A class for extracting exon CDS information of specified genes and transcripts
    from GTF annotation files and reference genome FASTA files.
    """

    def __init__(self, gtf_file: str, fasta_file: str):
        """
        Initialize the ExonExtractor with GTF and FASTA file paths.

        Parameters:
            gtf_file: Path to the GTF annotation file (human gene annotation,
                      such as Ensembl or Gencode GTF)
            fasta_file: Path to the reference genome FASTA file (used for extracting CDS sequences)
        """
        self.gtf_file = gtf_file
        self.fasta_file = fasta_file
        # Load genome data using pyfaidx for better performance
        self.genome = Fasta(fasta_file)

def extract_exons(self, gene_name: str, transcript_id: str) -> Tuple[List[Dict], str, str]:
    """
    Extract exon CDS information for a specified gene and transcript,
    generating a list of exon dictionaries.

    Parameters:
        gene_name: Target gene name (e.g., "BRCA1")
        transcript_id: Target transcript ID (e.g., "NM_007294.4")

    Returns:
        Tuple containing:
        - List[Dict]: A list in the standard format, where each element contains:
            - 'cds_sequence': CDS sequence of the exon (uppercase)
            - 'genomic_start': 0-based genomic start coordinate
            - 'genomic_end': 0-based genomic end coordinate (half-open interval)
        - str: Chromosome name
        - str: Strand direction (+/-)
    """
        # ------------------------------
        # Step 1: Parse GTF file, extract exon and CDS information for target transcript
        # ------------------------------
        exons = {}  # Store exon information (key: exon_id)
        cds_regions = {}  # Store CDS information (key: exon_id)

        # 支持.gz压缩文件
        if self.gtf_file.endswith('.gz'):
            file_handle = gzip.open(self.gtf_file, 'rt', encoding='utf-8')
        else:
            file_handle = open(self.gtf_file, 'r')
        
        with file_handle as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith('#'):
                    continue  # Skip comments and empty lines

                # GTF format: chromosome source feature_type start end score strand phase attributes
                parts = line.split('\t')
                if len(parts) != 9:
                    continue  # Skip invalid lines

                chrom, _, feature, start, end, _, strand, frame, attrs = parts

                # Parse attribute field
                attr_dict = self._parse_gtf_attributes(attrs)

                # Filter for target gene and transcript
                if (attr_dict.get('gene_name') != gene_name or
                        attr_dict.get('transcript_id') != transcript_id):
                    continue

                # Process exons
                if feature == 'exon':
                    exon_number = attr_dict.get('exon_number')
                    # GTF coordinates are 1-based closed intervals, convert to 0-based half-open intervals
                    exons[exon_number] = {
                        'chrom': chrom,
                        'start': int(start),  # 1-based start
                        'end': int(end),  # 1-based end (closed interval)
                        'strand': strand,
                        'exon_number': exon_number
                    }

                # Process CDS (only keep CDS belonging to exons)
                elif feature == 'CDS':
                    cds_start = int(start)
                    cds_end = int(end)
                    matched_exon = None

                    # Find the matching exon for this CDS region
                    for exon_id, exon in exons.items():
                        # Check if CDS is completely within the exon
                        if (exon['start'] <= cds_start and
                                cds_end <= exon['end'] and
                                exon['chrom'] == chrom):
                            matched_exon = exon_id
                            break

                    if matched_exon:
                        cds_regions[matched_exon] = (cds_start, cds_end)

        # Check if exons and CDS were found
        if not exons:
            raise ValueError(f"No exon information found for transcript {transcript_id}")
        if not cds_regions:
            raise ValueError(f"No CDS information found for transcript {transcript_id}")

        # ------------------------------
        # Step 2: Sort exons by genomic position (considering strand direction)
        # ------------------------------
        exon_list = sorted(exons.values(),
                           key=lambda x: x['start'] if x['strand'] == '+' else -x['start'])

        # Verify all exons belong to the same strand and chromosome
        chrom = exon_list[0]['chrom']
        strand = exon_list[0]['strand']
        for exon in exon_list:
            if exon['chrom'] != chrom:
                raise ValueError("Exons are distributed on different chromosomes, possibly incorrect annotation")
            if exon['strand'] != strand:
                raise ValueError("Exon strand directions are inconsistent, possibly incorrect annotation")

        # ------------------------------
        # Step 3: Extract CDS sequences from reference genome
        # ------------------------------
        chrom_seq = self.genome[chrom]  # Chromosome sequence using pyfaidx

        # Build result structure
        result_exons = []
        for exon in exon_list:
            exon_number = exon['exon_number']

            # Get CDS coordinates for this exon
            if exon_number not in cds_regions:
                continue  # Skip exons without CDS (non-coding exons)

            cds_start, cds_end = cds_regions[exon_number]

            # Extract CDS sequence (convert GTF 1-based to 0-based for Biopython)
            cds_seq = Seq(chrom_seq[cds_start - 1:cds_end].seq.upper())

            # If negative strand, take reverse complement
            if strand == '-':
                cds_seq = cds_seq.reverse_complement()

            # Build exon dictionary with 0-based coordinates
            result_exons.append({
                'cds_sequence': str(cds_seq),
                'genomic_start': cds_start,  # Convert to 0-based
                'genomic_end': cds_end  # 0-based (half-open interval)
            })

        # Check if valid CDS was generated
        if not result_exons:
            raise ValueError("No valid CDS sequences extracted, possibly a non-coding transcript")

        # 新增：返回染色体和链信息，与PostgreSQL格式对齐
        return result_exons, chrom, strand

    def _parse_gtf_attributes(self, attrs_str: str) -> Dict[str, str]:
        """
        Parse GTF attribute field (e.g., 'gene_name "BRCA1"; transcript_id "NM_007294.4";')

        Parameters:
            attrs_str: Attribute string from GTF file

        Returns:
            Dictionary of parsed attributes
        """
        attrs = {}
        # Regular expression to match key-value pairs (handles quoted and unquoted values)
        pattern = re.compile(r'(\w+)\s*["\']?([^; "\']+)["\']?;?')
        for match in pattern.findall(attrs_str):
            key, value = match
            attrs[key] = value.strip()
        return attrs

class GeneTranscript:
    # Codon to amino acid mapping table
    CODON_TABLE = {
        'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
        'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
        'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
        'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
        'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
        'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
        'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
        'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
        'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
        'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
        'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
        'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
        'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
        'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
        'TAC':'Y', 'TAT':'Y', 'TAA':'*', 'TAG':'*',
        'TGC':'C', 'TGT':'C', 'TGA':'*', 'TGG':'W'
    }

    def __init__(self, 
                 gene_name: str, 
                 transcript_id: str, 
                 exons: List[Dict],  # Exon structure with genomic coordinates
                 chromosome: str = "chr1",
                 strand: str = "+",
                 utr5: Optional[Dict] = None,  # 5' UTR region
                 utr3: Optional[Dict] = None):  # 3' UTR region
        """
        Initialize gene transcript with explicit exon coordinates (handles intron spacing)
        :param gene_name: Name of the gene
        :param transcript_id: Transcript ID (e.g., NM_000000.1)
        :param exons: List of exon dictionaries with structure:
                      {
                          'cds_sequence': 'ATG...',  # CDS portion of this exon
                          'genomic_start': 1000000,  # 0-based start coordinate in genome
                          'genomic_end': 1000020     # 0-based end coordinate in genome (exclusive)
                      }
        :param chromosome: Chromosome identifier
        :param strand: Strand direction ('+' or '-')
        :param utr5: 5' UTR region dict with 'genomic_start' and 'genomic_end'
        :param utr3: 3' UTR region dict with 'genomic_start' and 'genomic_end'
        """
        self.gene_name = gene_name
        self.transcript_id = transcript_id
        self.exons = exons
        self.chromosome = chromosome
        self.strand = strand if strand in ['+', '-'] else '+'
        self.utr5 = utr5
        self.utr3 = utr3
        self.valid_nucleotides = {'A', 'T', 'C', 'G'}
        
        # Validate exons and build CDS-genomic coordinate mapping
        self._validate_exons()
        self._build_cds_genomic_mapping()
    
    @staticmethod
    def _reverse_complement_static(sequence: str) -> str:
        """Static method to reverse complement a DNA sequence"""
        complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N'}
        return ''.join(complement.get(base, base) for base in sequence[::-1])

    def _validate_exons(self):
        """Validate exon structure and CDS sequences"""
        # Check exon order and non-overlapping
        for i in range(1, len(self.exons)):
            prev_end = self.exons[i-1]
            curr_start = self.exons[i]
            if self.strand == '+':
                if curr_start['genomic_start'] <= prev_end['genomic_end']:
                    raise ValueError(f"Sense strand Exons {i} and {i+1} overlap or are out of order")
            else:
                if curr_start['genomic_end'] >= prev_end['genomic_start']:
                    raise ValueError(f"Antisense strand Exons {i} and {i+1} overlap or are out of order")

        # Validate each exon's CDS
        for i, exon in enumerate(self.exons):
            cds = exon['cds_sequence'].upper()
            if not set(cds).issubset(self.valid_nucleotides):
                invalid = set(cds) - self.valid_nucleotides
                raise ValueError(f"Exon {i+1} contains invalid nucleotides: {invalid}")
     
            # Check genomic coordinates are valid
            if exon['genomic_start'] >= exon['genomic_end']:
                raise ValueError(f"Exon {i+1} has invalid coordinates (start >= end)")
            
            # Check CDS length is reasonable (should be <= coordinate span)
            coord_length = exon['genomic_end'] - exon['genomic_start'] + 1
            if len(cds) > coord_length:
                raise ValueError(f"Exon {i+1} CDS length ({len(cds)}) exceeds genomic span ({coord_length})")
            # Note: CDS length can be less than genomic span if exon contains UTR regions

        # Check total CDS length is multiple of 3
        total_cds_length = sum(len(exon['cds_sequence']) for exon in self.exons)
        if total_cds_length % 3 != 0:
            raise ValueError(f"Total CDS length must be multiple of 3, got {total_cds_length}")

    def _build_cds_genomic_mapping(self):
        """
        Create mapping from CDS position (0-based) to genomic coordinate
        Result: dict where key = cds_position, value = genomic_coordinate
        """
        self.cds_genomic_map = {}
        self.cds_sequence = ""
        current_cds_pos = 0
    
        for exon in self.exons:
            exon_cds = exon['cds_sequence'].upper()
            self.cds_sequence += exon_cds
    
            # Increment genomic position (handles intron gaps automatically)
            if self.strand == '+':
                # Sense strand: Add coordinates for each base in this exon's CDS
                genomic_pos = exon['genomic_start'] - 1  # Convert to 0-based
                step = 1
            else:
                # Antisense strand: Add coordinates for each base in this exon's CDS
                # For negative strand, we start from the end and go backwards
                # This ensures that the first CDS position maps to the highest genomic position
                genomic_pos = exon['genomic_end'] - 1  # Start from the end position (1-based)
                step = -1
    
            # Map each CDS base to genomic coordinates
            for base in exon_cds:
                self.cds_genomic_map[current_cds_pos] = genomic_pos
                current_cds_pos += 1
                genomic_pos += step  # Sense strand +1, antisense strand -1
    
        self.cds_length = current_cds_pos

    def get_genomic_coord(self, cds_position: int) -> int:
        """
        Convert CDS position (0-based) to genomic coordinate
        :param cds_position: Position in CDS sequence (0-based)
        :return: Corresponding genomic coordinate (0-based)
        """
        if cds_position < 0 or cds_position >= self.cds_length:
            raise ValueError(f"CDS position {cds_position} out of range (0-{self.cds_length-1})")
        return self.cds_genomic_map[cds_position]

    def is_synonymous(self, cds_position: int, mutant_nucleotide: str) -> bool:
        """
        Determine if a point mutation is synonymous
        :param cds_position: Position in CDS (0-based)
        :param mutant_nucleotide: Mutated nucleotide
        :return: True if synonymous
        """
        codon_start = (cds_position // 3) * 3
        codon_pos_in_triplet = cds_position % 3
        
        # Get original codon
        original_codon = self.cds_sequence[codon_start:codon_start+3]
        if len(original_codon) != 3:
            return False
        
        # Create mutant codon
        mutant_codon_list = list(original_codon)
        mutant_codon_list[codon_pos_in_triplet] = mutant_nucleotide
        mutant_codon = ''.join(mutant_codon_list)
        
        # Compare amino acids
        original_aa = self.CODON_TABLE.get(original_codon)
        mutant_aa = self.CODON_TABLE.get(mutant_codon)
        
        return original_aa == mutant_aa if (original_aa and mutant_aa) else False

    def _calculate_snv_score(self, cds_pos: int, is_syn: bool, is_stop: bool, codon_pos: int) -> float:
        """Calculate score for SNV variants (0-10 scale)"""
        score = 0.0
        
        # Base score based on synonymity (non-synonymous has higher score)
        score += 4.0 if not is_syn else 1.0
        
        # Add score if introduces stop codon
        if is_stop:
            score += 3.0
            
        # Adjust score based on codon position (1st and 2nd positions more impactful)
        if codon_pos in [1, 2]:
            score += 1.0
            
        # Normalize to 0-10 range
        return min(score, 10.0)

    def _calculate_insertion_score(self, insert_length: int, is_frameshift: bool, insert_pos: int) -> float:
        """Calculate score for insertion variants (0-10 scale)"""
        score = 0.0
        
        # Base score based on length (longer insertions have higher score)
        length_factor = min(insert_length / 6, 1.0)  # Normalize to 0-1
        score += 5.0 * length_factor
        
        # Add score for frameshift mutations
        if is_frameshift:
            score += 3.0
            
        # Adjust based on position (earlier in CDS has higher impact)
        position_factor = 1.0 - min(insert_pos / self.cds_length, 1.0)
        score += 2.0 * position_factor
        
        # Normalize to 0-10 range
        return min(score, 10.0)

    def _calculate_deletion_score(self, del_length: int, is_frameshift: bool, start_pos: int, end_pos: int) -> float:
        """Calculate score for deletion variants (0-10 scale)"""
        score = 0.0
        
        # Base score based on length (longer deletions have higher score)
        length_factor = min(del_length / 6, 1.0)  # Normalize to 0-1
        score += 5.0 * length_factor
        
        # Add score for frameshift mutations
        if is_frameshift:
            score += 3.0
            
        # Adjust based on position (earlier in CDS has higher impact)
        position_factor = 1.0 - min(start_pos / self.cds_length, 1.0)
        score += 2.0 * position_factor
        
        # Normalize to 0-10 range
        return min(score, 10.0)

    def _calculate_splice_score(self, offset: int, site_type: str, is_classic: bool = False, is_canonical: bool = False) -> float:
        """Calculate score for splice site variants (0-10 scale)
        
        Score decreases with distance from splice site:
        - Classic sites (±1, ±2): highest score
        - Canonical sites: high score  
        - Distant sites: lower score
        """
        # Use absolute offset for distance calculation
        abs_offset = abs(offset)
        
        # Score decreases with distance from splice site
        # Classic sites (offset ±1, ±2) get highest score
        if abs_offset <= 2:
            distance_factor = 1.0
        elif abs_offset <= 5:
            distance_factor = 0.8
        elif abs_offset <= 10:
            distance_factor = 0.6
        elif abs_offset <= 15:
            distance_factor = 0.4
        else:
            distance_factor = 0.2
        
        # Donor sites generally have higher impact than acceptor sites
        site_factor = 1.0 if site_type == 'donor' else 0.9
        
        # Calculate base score
        score = 8.0 * distance_factor * site_factor
        
        # Add bonus for classic splice sites (GT-AG at ±1, ±2)
        if is_classic:
            score += 1.5
            
        # Add bonus for canonical splice sites
        if is_canonical:
            score += 1.0
            
        # Normalize to 0-10 range
        return min(score, 10.0)

    def _calculate_inframe_score(self, 
                                length_factor: float, 
                                conservation_factor: float, 
                                domain_factor: float, 
                                aa_property_factor: float,
                                max_length: int = 6) -> float:
        """
        Calculate pathogenicity score for in-frame variants
        
        Dynamic Formula: S_inframe = min(2×L_f + 3×C_f + 3×D_f + 2×A_f, 10)
        Where:
        - L_f: Length factor (0.2-1.0) - dynamically adjusted based on max_length
        - C_f: Conservation factor (0.0-1.0) - based on domain conservation
        - D_f: Domain factor (0.1-1.0) - expanded range for better discrimination
        - A_f: Amino acid property factor (0.1-1.0) - expanded range
        - max_length: Maximum indel length set by user (default 6bp)
        
        This formula provides better score distribution and more balanced weighting
        with dynamic length factor adjustment based on user settings
        """
        # Apply non-linear transformation to length factor for better distribution
        # Use square root to reduce the dominance of length
        adjusted_length_factor = 0.2 + 0.8 * (length_factor ** 0.5)
        
        # Adjust length factor weight based on max_length to maintain score distribution
        # For larger max_length, reduce the impact of length factor
        if max_length <= 6:
            length_weight = 2.0  # Standard weight for small indels
        elif max_length <= 12:
            length_weight = 1.5  # Reduced weight for medium indels
        else:
            length_weight = 1.0  # Minimal weight for large indels
        
        # Apply logarithmic transformation to conservation factor for better spread
        # Higher conservation (lower factor) should have more impact
        adjusted_conservation_factor = conservation_factor ** 0.7
        
        # Expand domain factor range for better discrimination
        adjusted_domain_factor = 0.1 + 0.9 * domain_factor
        
        # Expand amino acid property factor range
        adjusted_aa_factor = 0.1 + 0.9 * aa_property_factor
        
        score = (length_weight * adjusted_length_factor) + \
                (3 * adjusted_conservation_factor) + \
                (3 * adjusted_domain_factor) + \
                (2 * adjusted_aa_factor)
        
        return min(score, 10.0)

    def _get_region_conservation(self, start: int, end: int, conserved_regions: List[Dict]) -> int:
        """
        Get conservation score (0-100) for a region based on protein domain data:
        - Domain/Family: 90 (highly conserved)
        - Motif: 85 (highly conserved)
        - Repeat: 70 (moderately conserved)
        - Coiled-coil: 60 (moderately conserved)
        - Disordered: 30 (low conservation)
        - Non-domain regions: 40 (low conservation)
        """
        # Get protein domains for this gene
        domains = self._get_protein_domains()
        
        if not domains:
            # Fallback to original method if no domain data available
            for region in conserved_regions:
                if region['start'] <= start and region['end'] >= end:
                    return region.get('conservation', 50)
            return 50  # Default to moderate conservation
        
        # Convert CDS positions to protein positions (1-based)
        protein_start = (start // 3) + 1
        protein_end = (end // 3) + 1
        
        # Find overlapping domains and calculate weighted conservation score
        conservation_scores = []
        total_overlap = 0
        
        for domain in domains:
            domain_start = domain['seq_start']
            domain_end = domain['seq_end']
            domain_type = domain['protein_domain']
            
            # Calculate overlap between variant and domain
            overlap_start = max(protein_start, domain_start)
            overlap_end = min(protein_end, domain_end)
            
            if overlap_start <= overlap_end:
                overlap_length = overlap_end - overlap_start + 1
                total_overlap += overlap_length
                
                # Assign conservation score based on domain type
                if domain_type in ['Domain', 'Family']:
                    conservation_score = 90  # Highly conserved
                elif domain_type == 'Motif':
                    conservation_score = 85  # Highly conserved
                elif domain_type == 'Repeat':
                    conservation_score = 70  # Moderately conserved
                elif domain_type == 'Coiled-coil':
                    conservation_score = 60  # Moderately conserved
                elif domain_type == 'Disordered':
                    conservation_score = 30  # Low conservation
                else:
                    conservation_score = 50  # Default
                
                conservation_scores.append((conservation_score, overlap_length))
        
        # Calculate weighted average conservation score
        if conservation_scores:
            total_weighted_score = sum(score * weight for score, weight in conservation_scores)
            weighted_avg = total_weighted_score / total_overlap
            
            # If there's partial overlap, blend with non-domain score
            variant_length = protein_end - protein_start + 1
            overlap_ratio = total_overlap / variant_length
            
            if overlap_ratio >= 0.8:
                return int(weighted_avg)  # Mostly overlapping with domains
            elif overlap_ratio >= 0.3:
                # Blend domain and non-domain scores
                non_domain_score = 40  # Low conservation for non-domain regions
                blended_score = weighted_avg * overlap_ratio + non_domain_score * (1 - overlap_ratio)
                return int(blended_score)
            else:
                # Mostly non-domain, but some domain influence
                non_domain_score = 40
                blended_score = weighted_avg * overlap_ratio + non_domain_score * (1 - overlap_ratio)
                return int(blended_score)
        else:
            # No domain overlap, use non-domain score
            return 40  # Low conservation for non-domain regions

    def _get_domain_factor(self, start: int, end: int) -> float:
        """
        Calculate domain factor using real protein domain data:
        - Within functional domain: 1.0
        - Domain boundary: 0.7
        - Non-functional region: 0.3
        """
        # Get protein domains for this gene
        domains = self._get_protein_domains()
        
        if not domains:
            # Fallback to simplified method if no domain data available
            mid_point = self.cds_length // 2
            if start < mid_point and end < mid_point:
                return 1.0  # Functional domain
            elif (start < mid_point and end >= mid_point) or (start < mid_point and end >= mid_point):
                return 0.7  # Boundary region
            else:
                return 0.3  # Non-functional region
        
        # Check if the variant region overlaps with any functional domains
        variant_start = start
        variant_end = end
        
        # Convert CDS positions to protein positions (1-based)
        protein_start = (variant_start // 3) + 1
        protein_end = (variant_end // 3) + 1
        variant_length = protein_end - protein_start + 1
        
        # Find overlapping domains
        overlapping_domains = []
        total_overlap = 0
        
        for domain in domains:
            domain_start = domain['seq_start']
            domain_end = domain['seq_end']
            
            # Calculate overlap between variant and domain
            overlap_start = max(protein_start, domain_start)
            overlap_end = min(protein_end, domain_end)
            
            if overlap_start <= overlap_end:
                overlap_length = overlap_end - overlap_start + 1
                total_overlap += overlap_length
                overlapping_domains.append({
                    'domain': domain,
                    'overlap_start': overlap_start,
                    'overlap_end': overlap_end,
                    'overlap_length': overlap_length
                })
        
        if not overlapping_domains:
            return 0.3  # No domain overlap
        
        # Calculate overlap ratio
        overlap_ratio = total_overlap / variant_length
        
        # Check for boundary conditions (variant near domain edges)
        is_boundary = False
        for overlap_info in overlapping_domains:
            domain = overlap_info['domain']
            domain_start = domain['seq_start']
            domain_end = domain['seq_end']
            
            # Check if variant is near domain boundary (within 10 amino acids)
            boundary_threshold = 10
            if (protein_start <= domain_start + boundary_threshold and protein_end >= domain_start) or \
               (protein_start <= domain_end and protein_end >= domain_end - boundary_threshold):
                is_boundary = True
                break
        
        # Determine factor based on overlap ratio and boundary conditions
        if overlap_ratio >= 0.8:
            return 1.0  # Mostly within functional domain
        elif overlap_ratio >= 0.5:
            return 0.8 if not is_boundary else 0.7  # Significant overlap
        elif overlap_ratio >= 0.2:
            return 0.7 if not is_boundary else 0.6  # Moderate overlap
        elif overlap_ratio > 0:
            return 0.6 if not is_boundary else 0.5  # Minimal overlap
        else:
            return 0.3  # No overlap
    
    def _get_protein_domains(self):
        """
        Get protein domains for the current gene from the protein domain database
        """
        # Cache domains to avoid repeated file reading
        if not hasattr(self, '_cached_domains'):
            self._cached_domains = self._load_protein_domains()
        
        return self._cached_domains.get(self.gene_name, [])
    
    def _load_protein_domains(self):
        """
        Load protein domain data from file and create gene-based index
        """
        domain_file = os.path.join(os.path.dirname(__file__), '..', '..', 'resources', 'proteininfo', 'protein_domain.txt')
        
        if not os.path.exists(domain_file):
            return {}
        
        gene_domains = {}
        
        try:
            with open(domain_file, 'r', encoding='utf-8') as f:
                # Skip header
                next(f)
                
                for line in f:
                    line = line.strip()
                    if not line:
                        continue
                    
                    parts = line.split('\t')
                    if len(parts) >= 5:
                        pfamseq_acc = parts[0]
                        seq_start = int(parts[1])
                        seq_end = int(parts[2])
                        protein_domain = parts[3]
                        gene_name = parts[4]
                        
                        if gene_name not in gene_domains:
                            gene_domains[gene_name] = []
                        
                        gene_domains[gene_name].append({
                            'pfamseq_acc': pfamseq_acc,
                            'seq_start': seq_start,
                            'seq_end': seq_end,
                            'protein_domain': protein_domain,
                            'gene_name': gene_name
                        })
        
        except Exception as e:
            print(f"Warning: Could not load protein domain data: {e}")
            return {}
        
        return gene_domains

    def _calculate_aa_property_factor(self, insert_pos: int, inserted_seq: str) -> float:
        """
        Calculate amino acid property factor:
        - Significant property change: 1.0
        - Similar property change: 0.3
        - Uncertain: 0.5
        """
        # 1. Translate inserted sequence to amino acids
        try:
            inserted_aa = self._translate_cds_to_aa(inserted_seq)
        except:
            return 0.5  # Return default for translation failures
        
        # 2. Get flanking amino acids around insertion position
        flanking_aa = self._get_flanking_amino_acids(insert_pos)
        if not flanking_aa:
            return 0.5
        
        # 3. Compare amino acid properties
        return self._compare_amino_acid_properties(flanking_aa, inserted_aa)

    def _translate_cds_to_aa(self, cds_seq: str) -> str:
        """Translate CDS sequence to amino acid sequence"""
        if len(cds_seq) % 3 != 0:
            raise ValueError("CDS sequence length must be multiple of 3")
        
        aa_seq = []
        for i in range(0, len(cds_seq), 3):
            codon = cds_seq[i:i+3]
            aa = self.CODON_TABLE.get(codon, 'X')  # 'X' represents unknown amino acid
            aa_seq.append(aa)
        return ''.join(aa_seq)

    def _get_flanking_amino_acids(self, cds_pos: int) -> Optional[str]:
        """Get flanking amino acids around insertion position (one on each side)"""
        if cds_pos == 0:
            # Insertion at start position, get only next amino acid
            if self.cds_length < 3:
                return None
            return self._translate_cds_to_aa(self.cds_sequence[0:3])[0]
        elif cds_pos >= self.cds_length:
            # Insertion at end position, get only previous amino acid
            if self.cds_length < 3:
                return None
            return self._translate_cds_to_aa(self.cds_sequence[-3:])[-1]
        else:
            # Insertion at middle position, get both flanking amino acids
            prev_codon_start = (cds_pos // 3) * 3
            next_codon_start = prev_codon_start + 3
            if next_codon_start >= self.cds_length:
                return None
            prev_aa = self._translate_cds_to_aa(self.cds_sequence[prev_codon_start:prev_codon_start+3])[-1]
            next_aa = self._translate_cds_to_aa(self.cds_sequence[next_codon_start:next_codon_start+3])[0]
            return prev_aa + next_aa

    def _compare_amino_acid_properties(self, flanking_aa: str, inserted_aa: str) -> float:
        """
        Compare physicochemical properties using real AAindex data
        """
        # AAindex data directly embedded (KYTJ820101: hydrophobicity, BIGC670101: volume, KLEP840101: charge)
        aa_properties = {
            'A': [-1.8, 52.6, 0.0],   # Ala: hydrophobicity, volume, charge
            'R': [4.5, 109.1, 1.0],   # Arg: hydrophobicity, volume, charge  
            'N': [3.5, 75.7, 0.0],    # Asn: hydrophobicity, volume, charge
            'D': [3.5, 68.4, -1.0],   # Asp: hydrophobicity, volume, charge
            'C': [-2.5, 68.3, 0.0],   # Cys: hydrophobicity, volume, charge
            'Q': [3.5, 89.7, 0.0],    # Gln: hydrophobicity, volume, charge
            'E': [3.5, 84.7, -1.0],   # Glu: hydrophobicity, volume, charge
            'G': [0.4, 36.3, 0.0],    # Gly: hydrophobicity, volume, charge
            'H': [3.2, 91.9, 0.0],    # His: hydrophobicity, volume, charge
            'I': [-4.5, 102.0, 0.0],  # Ile: hydrophobicity, volume, charge
            'L': [-3.9, 102.0, 0.0],  # Leu: hydrophobicity, volume, charge
            'K': [1.9, 105.1, 1.0],   # Lys: hydrophobicity, volume, charge
            'M': [2.8, 97.7, 0.0],    # Met: hydrophobicity, volume, charge
            'F': [1.6, 113.9, 0.0],   # Phe: hydrophobicity, volume, charge
            'P': [-1.6, 73.6, 0.0],   # Pro: hydrophobicity, volume, charge
            'S': [0.8, 54.9, 0.0],    # Ser: hydrophobicity, volume, charge
            'T': [0.7, 71.2, 0.0],    # Thr: hydrophobicity, volume, charge
            'W': [0.9, 135.4, 0.0],   # Trp: hydrophobicity, volume, charge
            'Y': [1.3, 116.2, 0.0],   # Tyr: hydrophobicity, volume, charge
            'V': [-4.2, 85.1, 0.0],  # Val: hydrophobicity, volume, charge
            'X': [0.0, 50.0, 0.0]     # Unknown: neutral values
        }
        
        # Normalize properties to 0-1 range
        def normalize_property(value, prop_type):
            if prop_type == 'hydrophobicity':
                return (value + 4.5) / 9.0  # Range: -4.5 to 4.5
            elif prop_type == 'volume':
                return (value - 36.3) / 99.1  # Range: 36.3 to 135.4
            elif prop_type == 'charge':
                return (value + 1.0) / 2.0  # Range: -1.0 to 1.0
            return 0.5
        
        # Get property vectors for flanking and inserted amino acids
        def get_property_vector(aa_seq):
            if not aa_seq:
                return [0.5, 0.5, 0.5]  # Default neutral values
            
            vectors = []
            for aa in aa_seq.upper():
                if aa in aa_properties:
                    props = aa_properties[aa]
                    norm_h = normalize_property(props[0], 'hydrophobicity')
                    norm_v = normalize_property(props[1], 'volume')
                    norm_c = normalize_property(props[2], 'charge')
                    vectors.append([norm_h, norm_v, norm_c])
                else:
                    vectors.append([0.5, 0.5, 0.5])  # Unknown amino acid
            
            # Return average vector
            if not vectors:
                return [0.5, 0.5, 0.5]
            return [sum(v[i] for v in vectors) / len(vectors) for i in range(3)]
        
        flanking_vec = get_property_vector(flanking_aa)
        inserted_vec = get_property_vector(inserted_aa)
        
        # Calculate Euclidean distance between property vectors
        import math
        distance = math.sqrt(sum((a - b) ** 2 for a, b in zip(flanking_vec, inserted_vec)))
        max_distance = math.sqrt(3)  # Maximum possible distance in 3D space
        
        # Convert distance to factor (0.3 = similar, 1.0 = very different)
        normalized_diff = distance / max_distance
        
        if normalized_diff < 0.3:
            return 0.3
        elif normalized_diff > 0.7:
            return 1.0
        else:
            return round(0.3 + (normalized_diff - 0.3) * (0.7 / 0.4), 1)

    def generate_snv_variants(self, 
                             synonymous: bool = True,
                             codon_positions: Optional[List[int]] = None,
                             include_stop_codon: bool = True) -> List[dict]:
        """Generate SNV variants with accurate genomic coordinates and HGVS validation"""
        snv_variants = []
        
        for cds_pos in range(self.cds_length):
            # CDS sequence (already in CDS direction)
            cds_nucleotide = self.cds_sequence[cds_pos]
            codon_pos = (cds_pos % 3) + 1  # 1-based position in codon
            
            # Filter by codon position
            if codon_positions and codon_pos not in codon_positions:
                continue
                
            for mutant_cds_nucleotide in self.valid_nucleotides:
                if mutant_cds_nucleotide == cds_nucleotide:
                    continue
                    
                is_syn = self.is_synonymous(cds_pos, mutant_cds_nucleotide)
                if not synonymous and is_syn:
                    continue
                    
                # Check for stop codon creation
                codon_start = (cds_pos // 3) * 3
                original_codon = self.cds_sequence[codon_start:codon_start+3]
                mutant_codon = original_codon[:codon_pos-1] + mutant_cds_nucleotide + original_codon[codon_pos:]
                is_stop = self.CODON_TABLE.get(mutant_codon) == '*'
                
                if not include_stop_codon and is_stop:
                    continue
                    
                # Get accurate genomic position
                genomic_pos = self.get_genomic_coord(cds_pos)
                
                # Get the actual genomic reference base at this position
                try:
                    genomic_ref_base = self._get_reference_base(genomic_pos + 1)  # Convert to 1-based
                    if not genomic_ref_base:
                        genomic_ref_base = cds_nucleotide  # Fallback
                except:
                    genomic_ref_base = cds_nucleotide  # Fallback
                
                # Determine genomic and HGVS sequences based on strand direction
                if self.strand == '+':
                    # Sense strand: CDS sequence = genomic sequence = HGVS sequence
                    genomic_ref = genomic_ref_base
                    genomic_alt = mutant_cds_nucleotide
                    hgvs_ref = cds_nucleotide
                    hgvs_alt = mutant_cds_nucleotide
                else:
                    # Antisense strand:                   
                    genomic_ref = genomic_ref_base
                    genomic_alt = self._reverse_complement(mutant_cds_nucleotide)
                    hgvs_ref = cds_nucleotide  
                    hgvs_alt = mutant_cds_nucleotide  
                
                # HGVS format validation
                if not self._validate_hgvs_format(hgvs_ref, hgvs_alt):
                    print(f"Warning: HGVS format validation failed: ref={hgvs_ref}, alt={hgvs_alt}")
                    continue
                
                # Generate HGVS format
                hgvs = f"c.{cds_pos+1}{hgvs_ref}>{hgvs_alt}"
                
                # Validate HGVS format
                if not self._validate_hgvs_notation(hgvs, 'SNV', cds_pos+1, hgvs_ref, hgvs_alt):
                    print(f"Warning: Incorrect HGVS format: {hgvs}")
                    continue
                
                # Calculate variant score
                score = self._calculate_snv_score(
                    cds_pos=cds_pos,
                    is_syn=is_syn,
                    is_stop=is_stop,
                    codon_pos=codon_pos
                )
                
                snv_variants.append({
                    'type': 'SNV',
                    'hgvs': hgvs,
                    'cds_position': cds_pos + 1,  # 1-based for HGVS
                    'original': genomic_ref,  # Genomic sequence (for VCF/BED REF column)
                    'mutant': genomic_alt,    # Genomic sequence (for VCF/BED ALT column)
                    'synonymous': is_syn,
                    'codon_position': codon_pos,
                    'is_stop_codon': is_stop,
                    'genomic_pos': genomic_pos,  # 0-based coordinate for VCF/BED
                    'ref_seq': genomic_ref,      # Genomic reference sequence
                    'alt_seq': genomic_alt,      # Genomic variant sequence
                    'hgvs_ref': hgvs_ref,        # HGVS reference sequence (CDS direction)
                    'hgvs_alt': hgvs_alt,        # HGVS variant sequence (CDS direction)
                    'Total_Score': round(score, 2)  # Added Total_Score field
                })
        
        return snv_variants

    def generate_insertion_variants(self, 
                                   max_length: int = 6, 
                                   min_length: int = 1,
                                   frameshift_only: bool = False) -> List[dict]:
        """Generate insertion variants with accurate genomic coordinates and HGVS validation"""
        insertion_variants = []
        
        # Insertions can happen between any two CDS bases (including start/end)
        for insert_pos in range(self.cds_length + 1):
            # Get genomic position for insertion (between CDS bases)
            if insert_pos == 0:
                # Insert before first CDS base
                genomic_pos = self.get_genomic_coord(0) if self.cds_length > 0 else self.exons[0]['genomic_start']
            elif insert_pos == self.cds_length:
                # Insert after last CDS base
                genomic_pos = self.get_genomic_coord(self.cds_length - 1) + 1
            else:
                # Insert between two CDS bases
                genomic_pos = self.get_genomic_coord(insert_pos - 1) + 1

            # Get reference base at insertion position
            if insert_pos == 0:
                ref_base = self.cds_sequence[0] if self.cds_length > 0 else 'N'
            else:
                ref_base = self.cds_sequence[insert_pos - 1]

            # Determine genomic reference base based on strand direction
            if self.strand == '+':
                genomic_ref_base = ref_base
            else:
                # Antisense strand: CDS sequence is already reverse complemented, need to reverse complement again to get genomic sequence
                genomic_ref_base = self._reverse_complement(ref_base)
            
            for insert_length in range(min_length, max_length + 1):
                # Filter frameshift insertions
                if frameshift_only and (insert_length % 3 == 0):
                    continue
                    
                for inserted_sequence in itertools.product(self.valid_nucleotides, repeat=insert_length):
                    inserted_sequence = ''.join(inserted_sequence)
                    
                    # Determine genomic and HGVS sequences based on strand direction
                    if self.strand == '+':
                        # Sense strand: Inserted sequence directly used for genomic and HGVS
                        genomic_alt_seq = inserted_sequence
                        hgvs_alt_seq = inserted_sequence
                    else:
                        # Antisense strand:
                        # - Genomic sequence needs reverse complement
                        # - HGVS sequence should show CDS sequence (not reverse complement)
                        genomic_alt_seq = self._reverse_complement(inserted_sequence)
                        hgvs_alt_seq = inserted_sequence  # CDS direction, already correct
                    
                    # HGVS format validation
                    if not self._validate_hgvs_format('', hgvs_alt_seq):
                        print(f"Warning: Insertion HGVS format validation failed: alt={hgvs_alt_seq}")
                        continue
                    
                    # HGVS notation
                    if insert_pos == 0:
                        hgvs = f"c.1_{2}ins{hgvs_alt_seq}"
                    elif insert_pos == self.cds_length:
                        hgvs = f"c.{insert_pos}_{insert_pos+1}ins{hgvs_alt_seq}"
                    else:
                        hgvs = f"c.{insert_pos}_{insert_pos+1}ins{hgvs_alt_seq}"
                    
                    # Validate HGVS format
                    if not self._validate_hgvs_notation(hgvs, 'insertion', insert_pos+1, '', hgvs_alt_seq):
                        print(f"Warning: Incorrect insertion HGVS format: {hgvs}")
                        continue
                    
                    # Calculate variant score
                    is_frameshift = (insert_length % 3 != 0)
                    score = self._calculate_insertion_score(
                        insert_length=insert_length,
                        is_frameshift=is_frameshift,
                        insert_pos=insert_pos
                    )
                    
                    insertion_variants.append({
                        'type': 'insertion',
                        'hgvs': hgvs,
                        'cds_position': insert_pos + 1,
                        'ref_base': genomic_ref_base,      # Genomic reference base
                        'inserted_sequence': inserted_sequence,  # Original inserted sequence (CDS direction)
                        'length': insert_length,
                        'is_frameshift': (insert_length % 3 != 0),
                        'genomic_pos': genomic_pos,
                        'genomic_alt_seq': genomic_alt_seq,  # Genomic variant sequence
                        'hgvs_alt_seq': hgvs_alt_seq,        # HGVS variant sequence (CDS direction)
                        'Total_Score': round(score, 2)        # Added Total_Score field
                    })
        
        return insertion_variants

    def generate_deletion_variants(self, 
                                  max_length: int = 6, 
                                  min_length: int = 1,
                                  frameshift_only: bool = False) -> List[dict]:
        """Generate deletion variants with accurate genomic coordinates and HGVS validation"""
        deletion_variants = []
        
        for start_pos in range(self.cds_length):
            for del_length in range(min_length, max_length + 1):
                end_pos = start_pos + del_length - 1
                if end_pos >= self.cds_length:
                    continue
                    
                # Filter frameshift deletions
                if frameshift_only and (del_length % 3 == 0):
                    continue
                    
                # Get genomic coordinates
                genomic_start = self.get_genomic_coord(start_pos)
                genomic_end = self.get_genomic_coord(end_pos) + 1  # End is exclusive in BED/VCF
                
                # Get deleted sequence
                deleted_cds_sequence = self.cds_sequence[start_pos:end_pos+1]
                
                # Determine genomic and HGVS sequences based on strand direction
                if self.strand == '+':
                    # Sense strand: CDS sequence = genomic sequence = HGVS sequence
                    genomic_deleted_sequence = deleted_cds_sequence
                    hgvs_deleted_sequence = deleted_cds_sequence
                else:
                    # Antisense strand:
                    # - CDS sequence is already reverse complemented (CDS direction)
                    # - Genomic sequence needs to reverse complement CDS sequence
                    # - HGVS sequence should show CDS sequence (not reverse complement)
                    genomic_deleted_sequence = self._reverse_complement(deleted_cds_sequence)
                    hgvs_deleted_sequence = deleted_cds_sequence  # CDS direction, already correct
                
                # HGVS format validation
                if not self._validate_hgvs_format(hgvs_deleted_sequence, ''):
                    print(f"Warning: Deletion HGVS format validation failed: ref={hgvs_deleted_sequence}")
                    continue
                
                # Generate HGVS format
                hgvs = f"c.{start_pos+1}_{end_pos+1}del{hgvs_deleted_sequence}"
                
                # Validate HGVS format
                if not self._validate_hgvs_notation(hgvs, 'deletion', start_pos+1, hgvs_deleted_sequence, ''):
                    print(f"Warning: Incorrect deletion HGVS format: {hgvs}")
                    continue
                
                # Calculate variant score
                is_frameshift = (del_length % 3 != 0)
                score = self._calculate_deletion_score(
                    del_length=del_length,
                    is_frameshift=is_frameshift,
                    start_pos=start_pos,
                    end_pos=end_pos
                )
                
                deletion_variants.append({
                    'type': 'deletion',
                    'hgvs': hgvs,
                    'cds_start': start_pos + 1,
                    'cds_end': end_pos + 1,
                    'deleted_sequence': genomic_deleted_sequence,  # Genomic deleted sequence (for VCF/BED)
                    'length': del_length,
                    'is_frameshift': (del_length % 3 != 0),
                    'genomic_start': genomic_start,
                    'genomic_end': genomic_end,
                    'ref_seq': genomic_deleted_sequence,  # Genomic reference sequence
                    'hgvs_ref_seq': hgvs_deleted_sequence,  # HGVS reference sequence (CDS direction)
                    'Total_Score': round(score, 2)          # Added Total_Score field
                })
        
        return deletion_variants

    def generate_splice_variants(self, 
                                max_intron_offset: int = 20,
                                min_intron_offset: int = 1,
                                include_classic_sites: bool = False) -> List[dict]:
        """Generate splice site variants with accurate genomic coordinates and HGVS validation"""
        splice_variants = []
        
        # Get all exon boundaries (splice sites)
        # For splice sites, we need to consider the actual intron-exon boundaries
        # Also include UTR boundaries if available
        exon_boundaries = []
        
        for i, exon in enumerate(self.exons):
            # For each exon, we need to identify the splice sites
            # Donor site: where intron starts (exon end for + strand, exon start for - strand)
            # Acceptor site: where intron ends (exon start for + strand, exon end for - strand)
            
            if self.strand == '+':
                # Positive strand: donor site is at exon end, acceptor site is at exon start
                # Donor site (exon end -> intron start) - only for non-last exons
                if i < len(self.exons) - 1:  # Not the last exon
                    exon_boundaries.append({
                        'position': exon['genomic_end'],
                        'type': 'donor',
                        'exon_index': i,
                        'is_exon_end': True,
                        'boundary_type': 'exon_end'
                    })
                # Acceptor site (intron end -> exon start) - only for non-first exons
                if i > 0:  # Not the first exon
                    exon_boundaries.append({
                        'position': exon['genomic_start'],
                        'type': 'acceptor',
                        'exon_index': i,
                        'is_exon_start': True,
                        'boundary_type': 'exon_start'
                    })
            else:
                # Negative strand: donor site is at exon start, acceptor site is at exon end
                # Donor site (exon start -> intron start) - only for non-last exons
                if i < len(self.exons) - 1:  # Not the last exon
                    exon_boundaries.append({
                        'position': exon['genomic_start'],
                        'type': 'donor',
                        'exon_index': i,
                        'is_exon_start': True,
                        'boundary_type': 'exon_start'
                    })
                # Acceptor site (intron end -> exon end) - only for non-first exons
                if i > 0:  # Not the first exon
                    exon_boundaries.append({
                        'position': exon['genomic_end'],
                        'type': 'acceptor',
                        'exon_index': i,
                        'is_exon_end': True,
                        'boundary_type': 'exon_end'
                    })
        
        # Add UTR boundaries if available
        if self.utr5:
            # 5' UTR boundaries
            if self.strand == '+':
                # 5' UTR end -> first exon start (acceptor site)
                exon_boundaries.append({
                    'position': self.utr5['genomic_end'],
                    'type': 'acceptor',
                    'exon_index': -1,  # Special index for UTR
                    'is_exon_start': True,
                    'boundary_type': 'utr5_end',
                    'region': 'utr5'
                })
            else:
                # 5' UTR start -> first exon start (acceptor site)
                exon_boundaries.append({
                    'position': self.utr5['genomic_start'],
                    'type': 'acceptor',
                    'exon_index': -1,  # Special index for UTR
                    'is_exon_start': True,
                    'boundary_type': 'utr5_start',
                    'region': 'utr5'
                })
        
        if self.utr3:
            # 3' UTR boundaries
            if self.strand == '+':
                # Last exon end -> 3' UTR start (donor site)
                exon_boundaries.append({
                    'position': self.utr3['genomic_start'],
                    'type': 'donor',
                    'exon_index': -2,  # Special index for UTR
                    'is_exon_end': True,
                    'boundary_type': 'utr3_start',
                    'region': 'utr3'
                })
            else:
                # Last exon end -> 3' UTR end (donor site)
                exon_boundaries.append({
                    'position': self.utr3['genomic_end'],
                    'type': 'donor',
                    'exon_index': -2,  # Special index for UTR
                    'is_exon_end': True,
                    'boundary_type': 'utr3_end',
                    'region': 'utr3'
                })
        
        # Sort boundaries by genomic position
        exon_boundaries.sort(key=lambda x: x['position'])
        
        for boundary in exon_boundaries:
            genomic_pos = boundary['position']
            site_type = boundary['type']
            exon_idx = boundary['exon_index']
            boundary_type = boundary['boundary_type']
            
            # Determine the range of intron positions to test
            # For donor sites: test positions in the intron (after exon end for + strand, before exon start for - strand)
            # For acceptor sites: test positions in the intron (before exon start for + strand, after exon end for - strand)
            
            if site_type == 'donor':
                if self.strand == '+':
                    # Donor site: test positions after exon end (in intron)
                    start_offset = min_intron_offset
                    end_offset = max_intron_offset
                else:
                    # Donor site: test positions before exon start (in intron)
                    start_offset = -max_intron_offset
                    end_offset = -min_intron_offset
            else:  # acceptor site
                if self.strand == '+':
                    # Acceptor site: test positions before exon start (in intron)
                    start_offset = -max_intron_offset
                    end_offset = -min_intron_offset
                else:
                    # Acceptor site: test positions after exon end (in intron)
                    start_offset = min_intron_offset
                    end_offset = max_intron_offset
            
            # Generate variants for this splice site
            for offset in range(start_offset, end_offset + 1):
                # Skip if offset is 0 (no change)
                if offset == 0:
                    continue
                
                # Check if this is a classic splice site
                is_classic = abs(offset) in [1, 2]
                
                # Skip non-classic sites if only classic sites are requested
                if include_classic_sites and not is_classic:
                    continue
                
                # Calculate the actual genomic position for the variant
                variant_genomic_pos = genomic_pos + offset
                
                # Skip if position would be negative
                if variant_genomic_pos < 0:
                    continue
                
                # Get reference sequence at the variant position
                # For PostgreSQL data, we need to extract from actual genome
                try:
                    ref_base = self._get_reference_base(variant_genomic_pos)
                    if not ref_base:
                        continue
                except (KeyError, IndexError):
                    # Skip if position is outside chromosome bounds
                    continue
                
                # Generate all possible nucleotide substitutions
                for alt_base in self.valid_nucleotides:
                    if alt_base == ref_base:
                        continue
                    
                    # Convert genomic position to CDS position for HGVS notation
                    cds_pos = self._genomic_to_cds_position(genomic_pos, site_type, exon_idx)
                    
                    # Determine HGVS notation based on strand and splice site type
                    if self.strand == '+':
                        if site_type == 'donor':
                            # Donor site: c.ExonEnd+offsetRef>Alt
                            hgvs = f"c.{cds_pos}+{offset}{ref_base}>{alt_base}"
                        else:  # acceptor
                            # Acceptor site: c.ExonStart-offsetRef>Alt
                            hgvs = f"c.{cds_pos}-{abs(offset)}{ref_base}>{alt_base}"
                    else:
                        # For negative strand, need to reverse complement the bases for HGVS notation
                        hgvs_ref = self._reverse_complement(ref_base)
                        hgvs_alt = self._reverse_complement(alt_base)
                        
                        if site_type == 'donor':
                            # Donor site: c.ExonStart+offsetRef>Alt (positive offset for intron 5' end)
                            hgvs = f"c.{cds_pos}+{abs(offset)}{hgvs_ref}>{hgvs_alt}"
                        else:  # acceptor
                            # Acceptor site: c.ExonEnd-offsetRef>Alt (negative offset for intron 3' end)
                            hgvs = f"c.{cds_pos}-{abs(offset)}{hgvs_ref}>{hgvs_alt}"
                    
                    # Validate HGVS format
                    if not self._validate_splice_hgvs_format(hgvs, site_type, offset):
                        continue
                    
                    # Determine if this is a canonical splice site variant
                    is_canonical = self._is_canonical_splice_site(ref_base, alt_base, site_type, offset)
                    
                    # Calculate variant score
                    score = self._calculate_splice_score(
                        offset=offset,
                        site_type=site_type,
                        is_classic=is_classic,
                        is_canonical=is_canonical
                    )
                    
                    splice_variants.append({
                        'type': 'splice_site',
                        'hgvs': hgvs,
                        'genomic_pos': variant_genomic_pos,
                        'ref_base': ref_base,
                        'alt_base': alt_base,
                        'site_type': site_type,
                        'exon_index': exon_idx,
                        'boundary_type': boundary_type,
                        'intron_offset': offset,
                        'is_classic': is_classic,
                        'is_canonical': is_canonical,
                        'severity': self._get_splice_severity(offset, is_classic, is_canonical),
                        'strand': self.strand,
                        'region': boundary.get('region', 'cds'),  # 'cds', 'utr5', or 'utr3'
                        'Total_Score': round(score, 2)           # Added Total_Score field
                    })
        
        return splice_variants

    def generate_inframe_variants(self, 
                                 min_length: int = 3, 
                                 max_length: int = 6,
                                 step: int = 3,  # In-frame variants must be multiples of 3
                                 conserved_regions: Optional[List[Dict]] = None) -> List[dict]:
        """
        Generate in-frame variants (In-frame Insertion/Deletion) and their scores
        In-frame variants refer to insertions/deletions with lengths that are multiples of 3,
        which do not alter the codon reading frame
        
        Parameters:
            min_length: Minimum length (default 3bp)
            max_length: Maximum length (default 6bp)
            step: Length step (fixed at 3 to ensure multiples of 3)
            conserved_regions: List of conserved regions, each containing
                              {'start': 0-based start position, 'end': 0-based end position, 'conservation': 0-100}
        
        Returns:
            List of in-frame variants with type, HGVS notation, score, etc.
        """
        inframe_variants = []
        conserved_regions = conserved_regions or []
        
        # Generate in-frame deletions
        for del_length in range(min_length, max_length + 1, step):
            for start_pos in range(self.cds_length):
                end_pos = start_pos + del_length - 1
                if end_pos >= self.cds_length:
                    continue
                    
                # Get deleted sequence and coordinates
                deleted_cds_sequence = self.cds_sequence[start_pos:end_pos+1]
                genomic_start = self.get_genomic_coord(start_pos)
                genomic_end = self.get_genomic_coord(end_pos) + 1  # Half-open interval
                
                # Strand direction handling
                if self.strand == '+':
                    genomic_deleted = deleted_cds_sequence
                    hgvs_deleted = deleted_cds_sequence
                else:
                    genomic_deleted = self._reverse_complement(deleted_cds_sequence)
                    hgvs_deleted = deleted_cds_sequence
                
                # Calculate conservation factor
                conservation = self._get_region_conservation(
                    start_pos, end_pos, conserved_regions)
                conservation_factor = min((100 - conservation) / 100, 1.0)
                
                # Calculate domain factor
                domain_factor = self._get_domain_factor(start_pos, end_pos)
                
                # Amino acid property factor (0.5 for in-frame deletions with no inserted amino acids)
                aa_factor = 0.5
                
                # Length factor
                length_factor = min(del_length / max_length, 1.0)
                
                # Calculate score
                score = self._calculate_inframe_score(
                    length_factor=length_factor,
                    conservation_factor=conservation_factor,
                    domain_factor=domain_factor,
                    aa_property_factor=aa_factor,
                    max_length=max_length
                )
                
                # HGVS notation
                hgvs = f"c.{start_pos+1}_{end_pos+1}del{hgvs_deleted}"
                
                if self._validate_hgvs_notation(hgvs, 'inframe_deletion', start_pos+1, hgvs_deleted, ''):
                    inframe_variants.append({
                        'type': 'inframe_deletion',
                        'hgvs': hgvs,
                        'cds_start': start_pos + 1,
                        'cds_end': end_pos + 1,
                        'length': del_length,
                        'deleted_sequence': genomic_deleted,
                        'genomic_start': genomic_start,
                        'genomic_end': genomic_end,
                        'conservation': conservation,
                        'Total_Score': round(score, 2)
                    })
        
        # Generate in-frame insertions
        for insert_length in range(min_length, max_length + 1, step):
            for insert_pos in range(self.cds_length + 1):
                # Determine genomic coordinate for insertion
                if insert_pos == 0:
                    genomic_pos = self.get_genomic_coord(0) if self.cds_length > 0 else self.exons[0]['genomic_start']
                elif insert_pos == self.cds_length:
                    genomic_pos = self.get_genomic_coord(self.cds_length - 1) + 1
                else:
                    genomic_pos = self.get_genomic_coord(insert_pos - 1) + 1

                # Generate all possible inserted sequences (lengths are multiples of 3)
                for inserted_seq in itertools.product(self.valid_nucleotides, repeat=insert_length):
                    inserted_seq = ''.join(inserted_seq)
                    
                    # Strand direction handling
                    if self.strand == '+':
                        genomic_inserted = inserted_seq
                        hgvs_inserted = inserted_seq
                    else:
                        genomic_inserted = self._reverse_complement(inserted_seq)
                        hgvs_inserted = inserted_seq
                    
                    # Calculate conservation factor (conservation at insertion position)
                    conservation = self._get_region_conservation(
                        insert_pos, insert_pos, conserved_regions)  # Conservation at point position
                    conservation_factor = min((100 - conservation) / 100, 1.0)
                    
                    # Calculate domain factor
                    domain_factor = self._get_domain_factor(insert_pos, insert_pos)
                    
                    # Calculate amino acid property factor
                    aa_factor = self._calculate_aa_property_factor(insert_pos, inserted_seq)
                    
                    # Length factor
                    length_factor = min(insert_length / max_length, 1.0)
                    
                    # Calculate score
                    score = self._calculate_inframe_score(
                        length_factor=length_factor,
                        conservation_factor=conservation_factor,
                        domain_factor=domain_factor,
                        aa_property_factor=aa_factor,
                        max_length=max_length
                    )
                    
                    # HGVS notation
                    if insert_pos == 0:
                        hgvs = f"c.1_{2}ins{hgvs_inserted}"
                    elif insert_pos == self.cds_length:
                        hgvs = f"c.{insert_pos}_{insert_pos+1}ins{hgvs_inserted}"
                    else:
                        hgvs = f"c.{insert_pos}_{insert_pos+1}ins{hgvs_inserted}"
                    
                    if self._validate_hgvs_notation(hgvs, 'inframe_insertion', insert_pos+1, '', hgvs_inserted):
                        inframe_variants.append({
                            'type': 'inframe_insertion',
                            'hgvs': hgvs,
                            'cds_position': insert_pos + 1,
                            'length': insert_length,
                            'inserted_sequence': genomic_inserted,
                            'genomic_pos': genomic_pos,
                            'conservation': conservation,
                            'Total_Score': round(score, 2)
                        })
        
        return inframe_variants

    def _get_reference_base(self, genomic_pos: int) -> str:
        """Get reference base at genomic position"""
        try:
            # Try to get from genome if available
            if hasattr(self, 'genome') and self.genome:
                chrom_seq = self.genome[self.chromosome]
                return str(chrom_seq[genomic_pos - 1].seq).upper()  # Convert to 1-based indexing
            else:
                # For PostgreSQL data without genome, use FASTA file
                return self._get_reference_base_from_fasta(genomic_pos)
        except (KeyError, IndexError, AttributeError):
            return None

    def _get_reference_base_from_fasta(self, genomic_pos: int) -> str:
        """Get reference base from FASTA file"""
        try:
            import pysam
            # Try to open the FASTA file
            fasta_path = f"resources/{self.chromosome}.fa"
            if not os.path.exists(fasta_path):
                # Fallback to Homo_sapiens_assembly37.fasta
                fasta_path = "resources/Homo_sapiens_assembly37.fasta"
            
            if os.path.exists(fasta_path):
                with pysam.FastaFile(fasta_path) as fasta:
                    # Handle different chromosome naming conventions
                    chrom_name = self.chromosome
                    if not chrom_name.startswith('chr'):
                        chrom_name = f"chr{chrom_name}"
                    
                    # Try different chromosome name formats
                    for chrom_variant in [chrom_name, self.chromosome, chrom_name.replace('chr', '')]:
                        try:
                            base = fasta.fetch(chrom_variant, genomic_pos - 1, genomic_pos).upper()
                            if base and base in 'ATCGN':
                                return base
                        except Exception as e:
                            print(f"Failed to fetch base for {chrom_variant} at {genomic_pos}: {e}")
                            continue
                    
                    # If all variants fail, try to get from available chromosomes
                    try:
                        available_chroms = fasta.references
                        print(f"Available chromosomes: {available_chroms}")
                        # Try to find a matching chromosome
                        for available_chrom in available_chroms:
                            if (self.chromosome in available_chrom or 
                                available_chrom in self.chromosome or
                                available_chrom.replace('chr', '') == self.chromosome.replace('chr', '')):
                                try:
                                    base = fasta.fetch(available_chrom, genomic_pos - 1, genomic_pos).upper()
                                    if base and base in 'ATCGN':
                                        print(f"Successfully fetched base {base} from {available_chrom} at {genomic_pos}")
                                        return base
                                except:
                                    continue
                    except:
                        pass
                    
                    # If all attempts fail, return N instead of A
                    print(f"Warning: Could not fetch reference base for {self.chromosome} at {genomic_pos}")
                    return 'N'
            else:
                print(f"FASTA file not found: {fasta_path}")
                return 'N'  # Return N instead of A when file not found
        except Exception as e:
            print(f"Error getting reference base from FASTA: {e}")
            return 'N'  # Return N instead of A on error

    @staticmethod
    def _extract_cds_sequence_from_genome(start_pos: int, end_pos: int, chromosome: str) -> str:
        """Extract CDS sequence from genome FASTA file using pyfaidx"""
        try:
            from pyfaidx import Fasta
            
            # Use the unified fasta file
            fasta_path = "resources/Homo_sapiens_assembly37.fasta"
            
            if not os.path.exists(fasta_path):
                print(f"FASTA file not found: {fasta_path}")
                return 'N' * (end_pos - start_pos + 1)
            
            with Fasta(fasta_path) as fasta:
                # Handle different chromosome naming conventions
                chrom_variants = [
                    chromosome,
                    f"chr{chromosome}" if not chromosome.startswith('chr') else chromosome,
                    chromosome.replace('chr', '') if chromosome.startswith('chr') else chromosome
                ]
                
                # Try different chromosome name formats
                for chrom_name in chrom_variants:
                    try:
                        if chrom_name in fasta:
                            sequence = str(fasta[chrom_name][start_pos-1:end_pos].seq).upper()
                            if sequence and len(sequence) > 0:
                                print(f"Successfully extracted sequence from {chrom_name}: {sequence[:20]}...")
                                return sequence
                    except Exception as e:
                        print(f"Failed to extract sequence from {chrom_name}: {e}")
                        continue
                
                # Try to find a matching chromosome
                try:
                    available_chroms = list(fasta.keys())
                    print(f"Available chromosomes: {available_chroms[:10]}...")
                    for available_chrom in available_chroms:
                        if (chromosome in available_chrom or 
                            available_chrom in chromosome or
                            available_chrom.replace('chr', '') == chromosome.replace('chr', '')):
                            try:
                                sequence = str(fasta[available_chrom][start_pos-1:end_pos].seq).upper()
                                if sequence and len(sequence) > 0:
                                    print(f"Successfully extracted sequence from {available_chrom}: {sequence[:20]}...")
                                    return sequence
                            except:
                                continue
                except:
                    pass
                
                # If all attempts fail, return N
                print(f"Warning: Could not extract sequence for {chromosome} from {start_pos} to {end_pos}")
                return 'N' * (end_pos - start_pos + 1)
                
        except Exception as e:
            print(f"Error extracting sequence from genome: {e}")
            return 'N' * (end_pos - start_pos + 1)

    def _is_canonical_splice_site(self, ref_base: str, alt_base: str, site_type: str, offset: int) -> bool:
        """Check if this is a canonical splice site variant"""
        # Canonical donor sites: GT (positions +1, +2)
        # Canonical acceptor sites: AG (positions -2, -1)
        
        if site_type == 'donor':
            if offset in [1, 2]:
                return (ref_base == 'G' and alt_base != 'G') or (ref_base == 'T' and alt_base != 'T')
        elif site_type == 'acceptor':
            if offset in [-2, -1]:
                return (ref_base == 'A' and alt_base != 'A') or (ref_base == 'G' and alt_base != 'G')
        
        return False

    def _get_splice_severity(self, offset: int, is_classic: bool, is_canonical: bool) -> str:
        """Determine the severity of the splice site variant"""
        if is_canonical:
            return 'high'
        elif is_classic:
            return 'moderate'
        else:
            return 'low'

    def _genomic_to_cds_position(self, genomic_pos: int, site_type: str, exon_idx: int) -> int:
        """Convert genomic position to CDS position for HGVS notation
        
        This function handles cases where CDS start position differs from 
        the genomic start of the first exon (e.g., when ATG is not at exon start)
        """
        if site_type == 'donor':
            if self.strand == '+':
                # Donor site: at exon end, CDS position is the end of the exon
                return self._get_exon_cds_end_position(exon_idx)
            else:
                # Donor site: at exon start, CDS position is the start of the exon
                return self._get_exon_cds_start_position(exon_idx)
        else:  # acceptor
            if self.strand == '+':
                # Acceptor site: at exon start, CDS position is the start of the exon
                return self._get_exon_cds_start_position(exon_idx)
            else:
                # Acceptor site: at exon end, CDS position is the start of the exon
                return self._get_exon_cds_start_position(exon_idx)
    
    def _get_exon_cds_start_position(self, exon_idx: int) -> int:
        """Get CDS start position for an exon (1-based)
        
        This function handles cases where CDS start position differs from 
        the genomic start of the first exon (e.g., when ATG is not at exon start)
        """
        if self.strand == '+':
            # Positive strand: count from start
            cds_pos = 1
            for i in range(exon_idx):
                cds_pos += len(self.exons[i]['cds_sequence'])
            
            # If this is the first exon, adjust for CDS start offset
            if exon_idx == 0:
                cds_start_offset = self._get_cds_start_offset()
                cds_pos += cds_start_offset
                
            return cds_pos
        else:
            # Negative strand: exons are ordered by genomic position (descending)
            # CDS position should be calculated from the end
            total_cds_length = sum(len(exon['cds_sequence']) for exon in self.exons)
            cds_pos = total_cds_length
            for i in range(exon_idx + 1, len(self.exons)):
                cds_pos -= len(self.exons[i]['cds_sequence'])
            
            # If this is the first exon, adjust for CDS start offset
            if exon_idx == 0:
                cds_start_offset = self._get_cds_start_offset()
                cds_pos -= cds_start_offset
                
            return cds_pos - len(self.exons[exon_idx]['cds_sequence']) + 1
    
    def _get_cds_start_offset(self) -> int:
        """Get the offset from genomic start to CDS start for the first exon
        
        This handles cases where the CDS start (ATG) is not at the beginning 
        of the first exon. Returns the number of bases from exon start to CDS start.
        """
        if not self.exons:
            return 0
            
        if self.strand == '+':
            # For positive strand: first exon is 5' end, ATG is at the beginning
            first_exon = self.exons[0]
            first_exon_seq = first_exon['cds_sequence']
            
            # Find the position of ATG in the first exon
            atg_pos = first_exon_seq.find('ATG')
            if atg_pos == -1:
                # If no ATG found, assume CDS starts at the beginning
                return 0
                
            return atg_pos
        else:
            # For negative strand: first exon is 3' end, need to find ATG in the last exon
            last_exon = self.exons[-1]
            last_exon_seq = last_exon['cds_sequence']
            
            # Find the position of ATG in the last exon (from the end)
            atg_pos = last_exon_seq.rfind('ATG')
            if atg_pos == -1:
                # If no ATG found, assume CDS starts at the beginning
                return 0
                
            # For negative strand, ATG position is from the end of the exon
            return len(last_exon_seq) - atg_pos - 3  # -3 because ATG is 3 bases long
    
    def _get_exon_cds_end_position(self, exon_idx: int) -> int:
        """Get CDS end position for an exon (1-based)"""
        if self.strand == '+':
            # Positive strand: count from start
            cds_pos = 0
            for i in range(exon_idx + 1):
                cds_pos += len(self.exons[i]['cds_sequence'])
            return cds_pos
        else:
            # Negative strand: exons are ordered by genomic position (descending)
            # CDS position should be calculated from the end
            total_cds_length = sum(len(exon['cds_sequence']) for exon in self.exons)
            cds_pos = total_cds_length
            for i in range(exon_idx + 1, len(self.exons)):
                cds_pos -= len(self.exons[i]['cds_sequence'])
            return cds_pos

    def _validate_splice_hgvs_format(self, hgvs: str, site_type: str, offset: int) -> bool:
        """Validate HGVS format for splice site variants"""
        try:
            # Basic format check
            if not hgvs.startswith('c.'):
                return False
            
            # Check for proper splice site notation
            if '+' in hgvs or '-' in hgvs:
                # Check if the offset matches the notation
                if '+' in hgvs:
                    parts = hgvs.split('+')
                    if len(parts) != 2:
                        return False
                    offset_part = parts[1].split('>')[0]
                    if not offset_part[:-1].isdigit():
                        return False
                    if int(offset_part[:-1]) != abs(offset):
                        return False
                elif '-' in hgvs:
                    parts = hgvs.split('-')
                    if len(parts) != 2:
                        return False
                    offset_part = parts[1].split('>')[0]
                    if not offset_part[:-1].isdigit():
                        return False
                    if int(offset_part[:-1]) != abs(offset):
                        return False
                
                # Check if bases are valid
                if '>' in hgvs:
                    base_part = hgvs.split('>')
                    if len(base_part) != 2:
                        return False
                    ref_base = base_part[0][-1]
                    alt_base = base_part[1]
                    if not all(base in 'ATCGN' for base in [ref_base, alt_base]):
                        return False
                
                return True
            
            return False
            
        except Exception:
            return False

    def generate_all_variants(self, 
                             variant_types: Optional[List[str]] = None,
                             max_indel_length: int = 15,
                             min_indel_length: int = 1,
                             synonymous: bool = True,
                             codon_positions: Optional[List[int]] = None,
                             include_stop_codon: bool = True,
                             frameshift_only: bool = False,
                             max_splice_offset: int = 20,
                             min_splice_offset: int = 1,
                             include_classic_splice_sites: bool = True,
                             max_variants: Optional[int] = None) -> dict:
        """Generate all variants with accurate genomic coordinates"""
        variant_types = variant_types or ['SNV', 'insertion', 'deletion']
        variants = {'snvs': [], 'insertions': [], 'deletions': [], 'splice_sites': [], 'inframe_variants': []}
        
        print(f"Generating variants for {self.gene_name} ({self.transcript_id})...")
        print(f"Total CDS length: {self.cds_length}bp across {len(self.exons)} exons")
        
        if 'SNV' in variant_types:
            variants['snvs'] = self.generate_snv_variants(
                synonymous=synonymous,
                codon_positions=codon_positions,
                include_stop_codon=include_stop_codon
            )
        
        if 'insertion' in variant_types:
            # Generate frameshift insertions only (exclude in-frame insertions)
            variants['insertions'] = self.generate_insertion_variants(
                max_length=max_indel_length,
                min_length=min_indel_length,
                frameshift_only=True  # Only generate frameshift insertions
            )
        
        if 'deletion' in variant_types:
            # Generate frameshift deletions only (exclude in-frame deletions)
            variants['deletions'] = self.generate_deletion_variants(
                max_length=max_indel_length,
                min_length=min_indel_length,
                frameshift_only=True  # Only generate frameshift deletions
            )
        
        if 'splice_site' in variant_types:
            variants['splice_sites'] = self.generate_splice_variants(
                max_intron_offset=max_splice_offset,
                min_intron_offset=min_splice_offset,
                include_classic_sites=include_classic_splice_sites
            )
        
        if 'inframe' in variant_types:
            # Generate in-frame variants (insertions/deletions with lengths multiples of 3)
            variants['inframe_variants'] = self.generate_inframe_variants(
                min_length=3,
                max_length=min(max_indel_length, 6),  # Limit to 6bp for in-frame variants
                step=3,
                conserved_regions=None  # Can be extended with real conservation data
            )
        
        # Apply variant limit
        total = 0
        if max_variants:
            for key in variants:
                remaining = max_variants - total
                if remaining <= 0:
                    variants[key] = []
                else:
                    variants[key] = variants[key][:remaining]
                    total += len(variants[key])
        else:
            total = sum(len(v) for v in variants.values())
        
        print(f"Generated {total} total variants")
        return {**variants, 'total': total}

    def export_to_vcf(self, variants, output_file):
        """Export variants to VCF with accurate genomic coordinates and HGVS validation"""
        # Statistics
        total_variants = 0
        valid_hgvs_count = 0
        invalid_hgvs_count = 0
        invalid_hgvs_list = []
        
        with open(output_file, 'w') as f:
            f.write("##fileformat=VCFv4.2\n")
            f.write(f"##source=GeneVariantGenerator (with genomic mapping)\n")
            f.write(f"##contig=<ID={self.chromosome},length=250000000>\n")  # Approx human chromosome length
            f.write(f"##INFO=<ID=TYPE,Number=1,Type=String,Description=\"Variant type\">\n")
            f.write(f"##INFO=<ID=HGVS,Number=1,Type=String,Description=\"HGVS notation\">\n")
            f.write(f"##INFO=<ID=SYNONYMOUS,Number=1,Type=Flag,Description=\"Synonymous mutation\">\n")
            f.write(f"##INFO=<ID=FRAMESHIFT,Number=1,Type=Flag,Description=\"Frameshift mutation\">\n")
            f.write(f"##INFO=<ID=SPLICE_TYPE,Number=1,Type=String,Description=\"Splice site type (donor/acceptor)\">\n")
            f.write(f"##INFO=<ID=SPLICE_SEVERITY,Number=1,Type=String,Description=\"Splice variant severity (high/moderate/low)\">\n")
            f.write(f"##INFO=<ID=IS_CANONICAL,Number=1,Type=Flag,Description=\"Canonical splice site variant\">\n")
            f.write(f"##INFO=<ID=MUTATION_PROB,Number=1,Type=Float,Description=\"Mutation probability from MuRaL model\">\n")
            f.write(f"##INFO=<ID=TRINUCLEOTIDE_CONTEXT,Number=1,Type=String,Description=\"Trinucleotide context\">\n")
            f.write(f"##INFO=<ID=IS_CONSTRAINED,Number=1,Type=Flag,Description=\"Gene is constrained (gnomAD LOEUF < 0.35)\">\n")
            f.write(f"##INFO=<ID=GNOAMD_CONSTRAINT,Number=1,Type=Float,Description=\"gnomAD constraint score (LOEUF)\">\n")
            f.write(f"##INFO=<ID=TOTAL_SCORE,Number=1,Type=Float,Description=\"Total variant impact score (0-10 scale)\">\n")
            f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
            
            var_id = 1
            # Export SNVs
            for var in variants['snvs']:
                total_variants += 1
                
                # Validate HGVS format
                hgvs = var['hgvs']
                if self._validate_vcf_hgvs_format(hgvs, 'SNV', var):
                    valid_hgvs_count += 1
                else:
                    invalid_hgvs_count += 1
                    invalid_hgvs_list.append(f"SNV {var_id}: {hgvs}")
                
                info = [
                    f"TYPE=SNV",
                    f"HGVS={hgvs}"
                ]
                if var['synonymous']:
                    info.append("SYNONYMOUS")
                info.append("FRAMESHIFT=false")
                
                # Add mutation spectrum information
                if 'mutation_probability' in var:
                    info.append(f"MUTATION_PROB={var['mutation_probability']}")
                if 'trinucleotide_context' in var:
                    info.append(f"TRINUCLEOTIDE_CONTEXT={var['trinucleotide_context']}")
                if 'is_constrained_gene' in var:
                    if var['is_constrained_gene']:
                        info.append("IS_CONSTRAINED")
                if 'gnomad_constraint_score' in var:
                    info.append(f"GNOAMD_CONSTRAINT={var['gnomad_constraint_score']}")
                
                # Add Total_Score information
                if 'Total_Score' in var:
                    info.append(f"TOTAL_SCORE={var['Total_Score']}")
                
                f.write(f"{self.chromosome}\t{var['genomic_pos']+1}\t{var_id}\t{var['original']}\t{var['mutant']}\t.\tPASS\t{';'.join(info)}\n")
                var_id += 1
            
            # Export insertions
            for var in variants['insertions']:
                total_variants += 1
                
                # Validate HGVS format
                hgvs = var['hgvs']
                if self._validate_vcf_hgvs_format(hgvs, 'INSERTION', var):
                    valid_hgvs_count += 1
                else:
                    invalid_hgvs_count += 1
                    invalid_hgvs_list.append(f"INSERTION {var_id}: {hgvs}")
                
                info = [
                    f"TYPE=INSERTION",
                    f"HGVS={hgvs}"
                ]
                if var['is_frameshift']:
                    info.append("FRAMESHIFT")
                
                # Add mutation spectrum information (if available)
                if 'mutation_probability' in var:
                    info.append(f"MUTATION_PROB={var['mutation_probability']}")
                if 'trinucleotide_context' in var:
                    info.append(f"TRINUCLEOTIDE_CONTEXT={var['trinucleotide_context']}")
                if 'is_constrained_gene' in var:
                    if var['is_constrained_gene']:
                        info.append("IS_CONSTRAINED")
                if 'gnomad_constraint_score' in var:
                    info.append(f"GNOAMD_CONSTRAINT={var['gnomad_constraint_score']}")
                
                # Add Total_Score information
                if 'Total_Score' in var:
                    info.append(f"TOTAL_SCORE={var['Total_Score']}")
                
                f.write(f"{self.chromosome}\t{var['genomic_pos']+1}\t{var_id}\t{var['ref_base']}\t{var['ref_base']+var['inserted_sequence']}\t.\tPASS\t{';'.join(info)}\n")
                var_id += 1
            
            # Export deletions
            for var in variants['deletions']:
                total_variants += 1
                
                # Validate HGVS format
                hgvs = var['hgvs']
                if self._validate_vcf_hgvs_format(hgvs, 'DELETION', var):
                    valid_hgvs_count += 1
                else:
                    invalid_hgvs_count += 1
                    invalid_hgvs_list.append(f"DELETION {var_id}: {hgvs}")
                
                info = [
                    f"TYPE=DELETION",
                    f"HGVS={hgvs}"
                ]
                if var['is_frameshift']:
                    info.append("FRAMESHIFT")
                
                # Add mutation spectrum information (if available)
                if 'mutation_probability' in var:
                    info.append(f"MUTATION_PROB={var['mutation_probability']}")
                if 'trinucleotide_context' in var:
                    info.append(f"TRINUCLEOTIDE_CONTEXT={var['trinucleotide_context']}")
                if 'is_constrained_gene' in var:
                    if var['is_constrained_gene']:
                        info.append("IS_CONSTRAINED")
                if 'gnomad_constraint_score' in var:
                    info.append(f"GNOAMD_CONSTRAINT={var['gnomad_constraint_score']}")
                
                # Add Total_Score information
                if 'Total_Score' in var:
                    info.append(f"TOTAL_SCORE={var['Total_Score']}")
                
                f.write(f"{self.chromosome}\t{var['genomic_start']+1}\t{var_id}\t{var['deleted_sequence']}\t.\t.\tPASS\t{';'.join(info)}\n")
                var_id += 1
            
            # Export splice site variants
            for var in variants['splice_sites']:
                total_variants += 1
                
                # Validate HGVS format
                hgvs = var['hgvs']
                if self._validate_vcf_hgvs_format(hgvs, 'SPLICE_SITE', var):
                    valid_hgvs_count += 1
                else:
                    invalid_hgvs_count += 1
                    invalid_hgvs_list.append(f"SPLICE_SITE {var_id}: {hgvs}")
                
                info = [
                    f"TYPE=SPLICE_SITE",
                    f"HGVS={hgvs}",
                    f"SPLICE_TYPE={var['site_type']}",
                    f"SPLICE_SEVERITY={var['severity']}",
                    f"REGION={var['region']}"
                ]
                if var['is_canonical']:
                    info.append("IS_CANONICAL")
                if var['is_classic']:
                    info.append("IS_CLASSIC")
                
                # Add mutation spectrum information (if available)
                if 'mutation_probability' in var:
                    info.append(f"MUTATION_PROB={var['mutation_probability']}")
                if 'trinucleotide_context' in var:
                    info.append(f"TRINUCLEOTIDE_CONTEXT={var['trinucleotide_context']}")
                if 'is_constrained_gene' in var:
                    if var['is_constrained_gene']:
                        info.append("IS_CONSTRAINED")
                if 'gnomad_constraint_score' in var:
                    info.append(f"GNOAMD_CONSTRAINT={var['gnomad_constraint_score']}")
                
                # Add Total_Score information
                if 'Total_Score' in var:
                    info.append(f"TOTAL_SCORE={var['Total_Score']}")
                
                f.write(f"{self.chromosome}\t{var['genomic_pos']+1}\t{var_id}\t{var['ref_base']}\t{var['alt_base']}\t.\tPASS\t{';'.join(info)}\n")
                var_id += 1
            
            # Export in-frame variants
            for var in variants['inframe_variants']:
                total_variants += 1
                
                # Validate HGVS format
                hgvs = var['hgvs']
                variant_type = var['type']
                if self._validate_vcf_hgvs_format(hgvs, variant_type, var):
                    valid_hgvs_count += 1
                else:
                    invalid_hgvs_count += 1
                    invalid_hgvs_list.append(f"{variant_type.upper()} {var_id}: {hgvs}")
                
                info = [
                    f"TYPE={variant_type.upper()}",
                    f"HGVS={hgvs}",
                    f"LENGTH={var['length']}",
                    f"CONSERVATION={var['conservation']}"
                ]
                
                # Add frameshift information (always false for in-frame variants)
                info.append("FRAMESHIFT=false")
                
                # Add Total_Score information
                if 'Total_Score' in var:
                    info.append(f"TOTAL_SCORE={var['Total_Score']}")
                
                # Write variant line based on type
                if variant_type == 'inframe_deletion':
                    f.write(f"{self.chromosome}\t{var['genomic_start']+1}\t{var_id}\t{var['deleted_sequence']}\t.\t.\tPASS\t{';'.join(info)}\n")
                elif variant_type == 'inframe_insertion':
                    f.write(f"{self.chromosome}\t{var['genomic_pos']+1}\t{var_id}\t.\t{var['inserted_sequence']}\t.\tPASS\t{';'.join(info)}\n")
                
                var_id += 1
        
        # Print validation results
        print(f"VCF exported to: {os.path.abspath(output_file)}")
        print(f"HGVS validation results: Total {total_variants} variants, valid HGVS: {valid_hgvs_count}, invalid HGVS: {invalid_hgvs_count}")
        if invalid_hgvs_list:
            print(f"Invalid HGVS list (first 10):")
            for invalid_hgvs in invalid_hgvs_list[:10]:
                print(f"  {invalid_hgvs}")
            if len(invalid_hgvs_list) > 10:
                print(f"  ... and {len(invalid_hgvs_list)-10} more invalid HGVS")

    def export_to_bed(self, variants, output_file):
        """Export variants to BED with accurate genomic coordinates and HGVS validation"""
        # Statistics
        total_variants = 0
        valid_hgvs_count = 0
        invalid_hgvs_count = 0
        invalid_hgvs_list = []
        
        with open(output_file, 'w') as f:
            f.write(f"# BED file for {self.gene_name} ({self.transcript_id})\n")
            f.write(f"# chrom\tstart\tend\tname\ttype\tattributes\n")
            
            # Export SNVs
            for var in variants['snvs']:
                total_variants += 1
                
                # Validate HGVS format
                hgvs = var['hgvs']
                if self._validate_bed_hgvs_format(hgvs, 'SNV', var):
                    valid_hgvs_count += 1
                else:
                    invalid_hgvs_count += 1
                    invalid_hgvs_list.append(f"SNV: {hgvs}")
                
                attrs = f"syn={var['synonymous']};stop={var['is_stop_codon']}"
                
                # Add mutation spectrum information
                if 'mutation_probability' in var:
                    attrs += f";mut_prob={var['mutation_probability']}"
                if 'trinucleotide_context' in var:
                    attrs += f";context={var['trinucleotide_context']}"
                if 'is_constrained_gene' in var:
                    attrs += f";constrained={var['is_constrained_gene']}"
                if 'gnomad_constraint_score' in var:
                    attrs += f";gnomad_constraint={var['gnomad_constraint_score']}"
                
                # Add Total_Score information
                if 'Total_Score' in var:
                    attrs += f";total_score={var['Total_Score']}"
                
                f.write(f"{self.chromosome}\t{var['genomic_pos']}\t{var['genomic_pos'] + 1}\t{hgvs}\tSNV\t{attrs}\n")
            
            # Export insertions
            for var in variants['insertions']:
                total_variants += 1
                
                # Validate HGVS format
                hgvs = var['hgvs']
                if self._validate_bed_hgvs_format(hgvs, 'INSERTION', var):
                    valid_hgvs_count += 1
                else:
                    invalid_hgvs_count += 1
                    invalid_hgvs_list.append(f"INSERTION: {hgvs}")
                
                attrs = f"length={var['length']};fs={var['is_frameshift']}"
                
                # Add mutation spectrum information (if available)
                if 'mutation_probability' in var:
                    attrs += f";mut_prob={var['mutation_probability']}"
                if 'trinucleotide_context' in var:
                    attrs += f";context={var['trinucleotide_context']}"
                if 'is_constrained_gene' in var:
                    if var['is_constrained_gene']:
                        attrs += f";constrained={var['is_constrained_gene']}"
                if 'gnomad_constraint_score' in var:
                    attrs += f";gnomad_constraint={var['gnomad_constraint_score']}"
                
                # Add Total_Score information
                if 'Total_Score' in var:
                    attrs += f";total_score={var['Total_Score']}"
                
                f.write(f"{self.chromosome}\t{var['genomic_pos']}\t{var['genomic_pos']+1}\t{hgvs}\tINSERTION\t{attrs}\n")
            
            # Export deletions
            for var in variants['deletions']:
                total_variants += 1
                
                # Validate HGVS format
                hgvs = var['hgvs']
                if self._validate_bed_hgvs_format(hgvs, 'DELETION', var):
                    valid_hgvs_count += 1
                else:
                    invalid_hgvs_count += 1
                    invalid_hgvs_list.append(f"DELETION: {hgvs}")
                
                attrs = f"length={var['length']};fs={var['is_frameshift']}"
                
                # Add mutation spectrum information (if available)
                if 'mutation_probability' in var:
                    attrs += f";mut_prob={var['mutation_probability']}"
                if 'trinucleotide_context' in var:
                    attrs += f";context={var['trinucleotide_context']}"
                if 'is_constrained_gene' in var:
                    if var['is_constrained_gene']:
                        attrs += f";constrained={var['is_constrained_gene']}"
                if 'gnomad_constraint_score' in var:
                    attrs += f";gnomad_constraint={var['gnomad_constraint_score']}"
                
                # Add Total_Score information
                if 'Total_Score' in var:
                    attrs += f";total_score={var['Total_Score']}"
                
                f.write(f"{self.chromosome}\t{var['genomic_start']}\t{var['genomic_end']}\t{hgvs}\tDELETION\t{attrs}\n")
            
            # Export splice site variants
            for var in variants['splice_sites']:
                total_variants += 1
                
                # Validate HGVS format
                hgvs = var['hgvs']
                if self._validate_bed_hgvs_format(hgvs, 'SPLICE_SITE', var):
                    valid_hgvs_count += 1
                else:
                    invalid_hgvs_count += 1
                    invalid_hgvs_list.append(f"SPLICE_SITE: {hgvs}")
                
                attrs = f"type={var['site_type']};severity={var['severity']};canonical={var['is_canonical']};classic={var['is_classic']};offset={var['intron_offset']};region={var['region']}"
                
                # Add mutation spectrum information (if available)
                if 'mutation_probability' in var:
                    attrs += f";mut_prob={var['mutation_probability']}"
                if 'trinucleotide_context' in var:
                    attrs += f";context={var['trinucleotide_context']}"
                if 'is_constrained_gene' in var:
                    attrs += f";constrained={var['is_constrained_gene']}"
                if 'gnomad_constraint_score' in var:
                    attrs += f";gnomad_constraint={var['gnomad_constraint_score']}"
                
                # Add Total_Score information
                if 'Total_Score' in var:
                    attrs += f";total_score={var['Total_Score']}"
                
                f.write(f"{self.chromosome}\t{var['genomic_pos']}\t{var['genomic_pos']+1}\t{hgvs}\tSPLICE_SITE\t{attrs}\n")
            
            # Export in-frame variants
            for var in variants['inframe_variants']:
                total_variants += 1
                
                # Validate HGVS format
                hgvs = var['hgvs']
                variant_type = var['type']
                if self._validate_bed_hgvs_format(hgvs, variant_type, var):
                    valid_hgvs_count += 1
                else:
                    invalid_hgvs_count += 1
                    invalid_hgvs_list.append(f"{variant_type.upper()}: {hgvs}")
                
                attrs = f"length={var['length']};conservation={var['conservation']};frameshift=false"
                
                # Add Total_Score information
                if 'Total_Score' in var:
                    attrs += f";total_score={var['Total_Score']}"
                
                # Write variant line based on type
                if variant_type == 'inframe_deletion':
                    f.write(f"{self.chromosome}\t{var['genomic_start']}\t{var['genomic_end']}\t{hgvs}\tINFRAME_DELETION\t{attrs}\n")
                elif variant_type == 'inframe_insertion':
                    f.write(f"{self.chromosome}\t{var['genomic_pos']}\t{var['genomic_pos']+1}\t{hgvs}\tINFRAME_INSERTION\t{attrs}\n")
        
        # Print validation results
        print(f"BED exported to: {os.path.abspath(output_file)}")
        print(f"HGVS validation results: Total {total_variants} variants, valid HGVS: {valid_hgvs_count}, invalid HGVS: {invalid_hgvs_count}")
        if invalid_hgvs_list:
            print(f"Invalid HGVS list (first 10):")
            for invalid_hgvs in invalid_hgvs_list[:10]:
                print(f"  {invalid_hgvs}")
            if len(invalid_hgvs_list) > 10:
                print(f"  ... and {len(invalid_hgvs_list)-10} more invalid HGVS")

    def export_variants(self, variants, output_file, file_format='vcf'):
        if file_format.lower() == 'bed':
            self.export_to_bed(variants, output_file)
        elif file_format.lower() == 'vcf':
            self.export_to_vcf(variants, output_file)
        else:
            raise ValueError(f"Unsupported format: {file_format}. Use 'bed' or 'vcf'.")

    def _validate_hgvs_format(self, ref_seq: str, alt_seq: str) -> bool:
        """Validate HGVS format validity"""
        # Check if sequences contain only valid bases
        valid_bases = {'A', 'T', 'C', 'G', 'N'}
        
        for seq in [ref_seq, alt_seq]:
            if not all(base in valid_bases for base in seq):
                return False
        
        # Check sequence length
        if len(ref_seq) == 0 and len(alt_seq) == 0:
            return False
        
        return True
    
    def _reverse_complement(self, sequence: str) -> str:
        """Calculate reverse complement sequence"""
        complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N'}
        return ''.join(complement.get(base, base) for base in reversed(sequence))
    
    def _validate_hgvs_notation(self, hgvs: str, variant_type: str, cds_pos: int, ref_seq: str, alt_seq: str) -> bool:
        """Validate HGVS format correctness"""
        try:
            # Basic format check
            if not hgvs.startswith('c.'):
                return False
            
            # Check format based on variant type
            if variant_type == 'SNV':
                # Format: c.123A>G
                expected_format = f"c.{cds_pos}{ref_seq}>{alt_seq}"
                return hgvs == expected_format
            
            elif variant_type == 'insertion':
                # Format: c.123_124insABC
                if 'ins' not in hgvs:
                    return False
                # Check if inserted sequence is correct
                if not hgvs.endswith(alt_seq):
                    return False
                return True
            
            elif variant_type == 'deletion':
                # Format: c.123_125delABC
                if 'del' not in hgvs:
                    return False
                # Check if deleted sequence is correct
                if not hgvs.endswith(ref_seq):
                    return False
                return True
            
            elif variant_type == 'inframe_insertion':
                # Format: c.123_124insABC (3-6bp insertions)
                if 'ins' not in hgvs:
                    return False
                # Check if inserted sequence is correct and length is multiple of 3
                if not hgvs.endswith(alt_seq):
                    return False
                if len(alt_seq) % 3 != 0 or len(alt_seq) < 3 or len(alt_seq) > 6:
                    return False
                return True
            
            elif variant_type == 'inframe_deletion':
                # Format: c.123_125delABC (3-6bp deletions)
                if 'del' not in hgvs:
                    return False
                # Check if deleted sequence is correct and length is multiple of 3
                if not hgvs.endswith(ref_seq):
                    return False
                if len(ref_seq) % 3 != 0 or len(ref_seq) < 3 or len(ref_seq) > 6:
                    return False
                return True
            
            return True
            
        except Exception:
            return False

    def _validate_vcf_hgvs_format(self, hgvs: str, variant_type: str, variant_data: dict) -> bool:
        """Validate HGVS format in VCF file INFO field"""
        try:
            # Basic format check
            if not hgvs.startswith('c.'):
                return False
            
            # Check if HGVS format is included in INFO field
            if 'HGVS=' not in f"HGVS={hgvs}":
                return False
            
            # Detailed validation based on variant type
            if variant_type == 'SNV':
                # Check SNV format: c.123A>G
                if '>' not in hgvs:
                    return False
                # Check if position is numeric
                parts = hgvs.split('>')
                if len(parts) != 2:
                    return False
                position_part = parts[0].replace('c.', '')
                if not position_part[:-1].isdigit():
                    return False
                # Check if bases are valid
                ref_base = position_part[-1]
                alt_base = parts[1]
                if not all(base in 'ATCGN' for base in [ref_base, alt_base]):
                    return False
                return True
            
            elif variant_type == 'INSERTION':
                # Check insertion format: c.123_124insABC
                if 'ins' not in hgvs:
                    return False
                # Check position range
                if '_' not in hgvs:
                    return False
                position_part = hgvs.split('ins')[0].replace('c.', '')
                if not position_part.replace('_', '').isdigit():
                    return False
                # Check inserted sequence
                inserted_seq = hgvs.split('ins')[1]
                if not all(base in 'ATCGN' for base in inserted_seq):
                    return False
                return True
            
            elif variant_type == 'DELETION':
                # Check deletion format: c.123_125delABC
                if 'del' not in hgvs:
                    return False
                # Check position range
                if '_' not in hgvs:
                    return False
                position_part = hgvs.split('del')[0].replace('c.', '')
                if not position_part.replace('_', '').isdigit():
                    return False
                # Check deleted sequence
                deleted_seq = hgvs.split('del')[1]
                if not all(base in 'ATCGN' for base in deleted_seq):
                    return False
                return True
            
            elif variant_type == 'SPLICE_SITE':
                # Check splice site format: c.123+1A>G or c.123-2A>G
                if '+' in hgvs or '-' in hgvs:
                    if '>' not in hgvs:
                        return False
                    # Check if bases are valid
                    base_part = hgvs.split('>')
                    if len(base_part) != 2:
                        return False
                    ref_base = base_part[0][-1]
                    alt_base = base_part[1]
                    if not all(base in 'ATCGN' for base in [ref_base, alt_base]):
                        return False
                    return True
                return False
            
            elif variant_type == 'inframe_insertion':
                # Check in-frame insertion format: c.123_124insABC (3-6bp, multiple of 3)
                if 'ins' not in hgvs:
                    return False
                # Check position range
                if '_' not in hgvs:
                    return False
                position_part = hgvs.split('ins')[0].replace('c.', '')
                if not position_part.replace('_', '').isdigit():
                    return False
                # Check inserted sequence
                inserted_seq = hgvs.split('ins')[1]
                if not all(base in 'ATCGN' for base in inserted_seq):
                    return False
                # Check length is multiple of 3 and within range
                if len(inserted_seq) % 3 != 0 or len(inserted_seq) < 3 or len(inserted_seq) > 6:
                    return False
                return True
            
            elif variant_type == 'inframe_deletion':
                # Check in-frame deletion format: c.123_125delABC (3-6bp, multiple of 3)
                if 'del' not in hgvs:
                    return False
                # Check position range
                if '_' not in hgvs:
                    return False
                position_part = hgvs.split('del')[0].replace('c.', '')
                if not position_part.replace('_', '').isdigit():
                    return False
                # Check deleted sequence
                deleted_seq = hgvs.split('del')[1]
                if not all(base in 'ATCGN' for base in deleted_seq):
                    return False
                # Check length is multiple of 3 and within range
                if len(deleted_seq) % 3 != 0 or len(deleted_seq) < 3 or len(deleted_seq) > 6:
                    return False
                return True
            
            return True
            
        except Exception:
            return False
    
    def _validate_bed_hgvs_format(self, hgvs: str, variant_type: str, variant_data: dict) -> bool:
        """Validate HGVS format in BED file name column"""
        try:
            # Basic format check
            if not hgvs.startswith('c.'):
                return False
            
            # HGVS format in BED file should be the same as VCF, but as name column
            # Check if it contains tab or other characters not allowed in BED format
            if '\t' in hgvs or '\n' in hgvs:
                return False
            
            # Detailed validation based on variant type
            if variant_type == 'SNV':
                # Check SNV format: c.123A>G
                if '>' not in hgvs:
                    return False
                # Check if position is numeric
                parts = hgvs.split('>')
                if len(parts) != 2:
                    return False
                position_part = parts[0].replace('c.', '')
                if not position_part[:-1].isdigit():
                    return False
                # Check if bases are valid
                ref_base = position_part[-1]
                alt_base = parts[1]
                if not all(base in 'ATCGN' for base in [ref_base, alt_base]):
                    return False
                return True
            
            elif variant_type == 'INSERTION':
                # Check insertion format: c.123_124insABC
                if 'ins' not in hgvs:
                    return False
                # Check position range
                if '_' not in hgvs:
                    return False
                position_part = hgvs.split('ins')[0].replace('c.', '')
                if not position_part.replace('_', '').isdigit():
                    return False
                # Check inserted sequence
                inserted_seq = hgvs.split('ins')[1]
                if not all(base in 'ATCGN' for base in inserted_seq):
                    return False
                return True
            
            elif variant_type == 'DELETION':
                # Check deletion format: c.123_125delABC
                if 'del' not in hgvs:
                    return False
                # Check position range
                if '_' not in hgvs:
                    return False
                position_part = hgvs.split('del')[0].replace('c.', '')
                if not position_part.replace('_', '').isdigit():
                    return False
                # Check deleted sequence
                deleted_seq = hgvs.split('del')[1]
                if not all(base in 'ATCGN' for base in deleted_seq):
                    return False
                return True
            
            elif variant_type == 'SPLICE_SITE':
                # Check splice site format: c.123+1A>G or c.123-2A>G
                if '+' in hgvs or '-' in hgvs:
                    if '>' not in hgvs:
                        return False
                    # Check if bases are valid
                    base_part = hgvs.split('>')
                    if len(base_part) != 2:
                        return False
                    ref_base = base_part[0][-1]
                    alt_base = base_part[1]
                    if not all(base in 'ATCGN' for base in [ref_base, alt_base]):
                        return False
                    return True
                return False
            
            elif variant_type == 'inframe_insertion':
                # Check in-frame insertion format: c.123_124insABC (3-6bp, multiple of 3)
                if 'ins' not in hgvs:
                    return False
                # Check position range
                if '_' not in hgvs:
                    return False
                position_part = hgvs.split('ins')[0].replace('c.', '')
                if not position_part.replace('_', '').isdigit():
                    return False
                # Check inserted sequence
                inserted_seq = hgvs.split('ins')[1]
                if not all(base in 'ATCGN' for base in inserted_seq):
                    return False
                # Check length is multiple of 3 and within range
                if len(inserted_seq) % 3 != 0 or len(inserted_seq) < 3 or len(inserted_seq) > 6:
                    return False
                return True
            
            elif variant_type == 'inframe_deletion':
                # Check in-frame deletion format: c.123_125delABC (3-6bp, multiple of 3)
                if 'del' not in hgvs:
                    return False
                # Check position range
                if '_' not in hgvs:
                    return False
                position_part = hgvs.split('del')[0].replace('c.', '')
                if not position_part.replace('_', '').isdigit():
                    return False
                # Check deleted sequence
                deleted_seq = hgvs.split('del')[1]
                if not all(base in 'ATCGN' for base in deleted_seq):
                    return False
                # Check length is multiple of 3 and within range
                if len(deleted_seq) % 3 != 0 or len(deleted_seq) < 3 or len(deleted_seq) > 6:
                    return False
                return True
            
            return True
            
        except Exception:
            return False

def main():
    """Command - line main entry point: Parse arguments → Pre - process GTF/FASTA → Generate variants → Export files"""
    parser = argparse.ArgumentParser(
        description="Gene Variant Generation Tool: Generate SNV/InDel/splice - site variants from GTF + reference genome FASTA, and output in VCF/BED format",
        formatter_class=argparse.RawTextHelpFormatter
    )

    # ==================== Core input parameters (required) ====================
    required = parser.add_argument_group('Required Arguments')
    required.add_argument('--gtf', required=True, type=str, help='Path to the gene annotation GTF file (supports.gz compression)')
    required.add_argument('--fasta', required=True, type=str, help='Path to the reference genome FASTA file')
    required.add_argument('--gene-name', required=True, type=str, help='Target gene name (e.g., BRCA1)')
    required.add_argument('--transcript-id', required=True, type=str, help='Target transcript ID (e.g., NM_007294.4)')

    # ==================== Output parameters ====================
    output = parser.add_argument_group('Output Arguments')
    output.add_argument('--output-vcf', type=str, default='variants.vcf', help='VCF output path (default: variants.vcf)')
    output.add_argument('--output-bed', type=str, default='', help='BED output path (default: not generated)')
    output.add_argument('--output-format', choices=['vcf', 'bed', 'both'], default='vcf', help='Output format (default: vcf)')

    # ==================== Variant generation configuration parameters ====================
    variant = parser.add_argument_group('Variant Generation Arguments')
    variant.add_argument('--variant-types', nargs='+', default=['SNV', 'insertion', 'deletion'],
                         choices=['SNV', 'insertion', 'deletion','splice_site', 'inframe'],
                         help='Variant types to be generated (default: SNV insertion deletion)')
    variant.add_argument('--max-indel-length', type=int, default=15, help='Maximum insertion/deletion length (default: 15)')
    variant.add_argument('--min-indel-length', type=int, default=1, help='Minimum insertion/deletion length (default: 1)')
    variant.add_argument('--synonymous', action='store_true', default=True, help='Generate synonymous mutations (default: on)')
    variant.add_argument('--no-synonymous', action='store_false', dest='synonymous', help='Do not generate synonymous mutations')
    variant.add_argument('--include-stop-codon', action='store_true', default=True, help='Include stop - codon mutations (default: on)')
    variant.add_argument('--no-stop-codon', action='store_false', dest='include_stop_codon', help='Do not generate stop - codon mutations')
    variant.add_argument('--max-splice-offset', type=int, default=20, help='Maximum intron offset for splice sites (default: 20)')
    variant.add_argument('--min-splice-offset', type=int, default=1, help='Minimum intron offset for splice sites (default: 1)')
    variant.add_argument('--include-classic-splice', action='store_true', default=True, help='Include classic splice sites (±1/±2, default: on)')
    variant.add_argument('--max-variants', type=int, default=None, help='Maximum number of variants to generate (default: no limit)')

    # Parse arguments
    args = parser.parse_args()

    # ==================== GeneTranscript ====================
    print(f"Extracting gene structure from GTF/FASTA: {args.gene_name} | {args.transcript_id}")
    extractor = ExonExtractor(args.gtf, args.fasta)
    exons, chromosome, strand = extractor.extract_exons(args.gene_name, args.transcript_id)

    # Instantiate GeneTranscript
    transcript = GeneTranscript(
        gene_name=args.gene_name,
        transcript_id=args.transcript_id,
        exons=exons,
        chromosome=chromosome,
        strand=strand
    )
    print(f"Gene structure loaded: Chromosome {chromosome} | Strand {strand} | CDS length {transcript.cds_length}bp")

    # ==================== Generate all variants ====================
    variants = transcript.generate_all_variants(
        variant_types=args.variant_types,
        max_indel_length=args.max_indel_length,
        min_indel_length=args.min_indel_length,
        synonymous=args.synonymous,
        include_stop_codon=args.include_stop_codon,
        max_splice_offset=args.max_splice_offset,
        min_splice_offset=args.min_splice_offset,
        include_classic_splice_sites=args.include_classic_splice,
        max_variants=args.max_variants
    )

    # ==================== Export files ====================
    if args.output_format in ['vcf', 'both']:
        transcript.export_to_vcf(variants, args.output_vcf)
    if args.output_format in ['bed', 'both'] or args.output_bed:
        bed_path = args.output_bed if args.output_bed else 'variants.bed'
        transcript.export_to_bed(variants, bed_path)

    print("All tasks completed!")

if __name__ == '__main__':
    main()