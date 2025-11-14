"""
DNA Sequence Utilities

This module provides utilities for DNA sequence validation, processing, and analysis.
Supports FASTA format parsing, sequence statistics, and basic bioinformatics operations.
"""

import re
from collections import Counter


class DNASequence:
    """
    Represents a DNA sequence with validation and utility methods.

    Attributes:
        sequence (str): The DNA sequence (uppercase)
        name (str): Sequence name/identifier
        description (str): Optional description
    """

    VALID_NUCLEOTIDES = set('ATCGN')
    COMPLEMENT_MAP = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N'}

    def __init__(self, sequence, name="", description=""):
        """
        Initialize DNA sequence.

        Args:
            sequence (str): DNA sequence string
            name (str): Sequence identifier
            description (str): Optional description

        Raises:
            ValueError: If sequence contains invalid characters
        """
        self.sequence = self._validate_and_clean(sequence)
        self.name = name
        self.description = description

    def _validate_and_clean(self, sequence):
        """Validate and clean the DNA sequence."""
        cleaned = re.sub(r'\s+', '', sequence.upper())

        for char in cleaned:
            if char not in self.VALID_NUCLEOTIDES:
                raise ValueError(f"Invalid nucleotide '{char}' in sequence. "
                               f"Only A, T, C, G, N allowed.")

        return cleaned

    @property
    def length(self):
        """Get sequence length."""
        return len(self.sequence)

    @property
    def gc_content(self):
        """Calculate GC content as percentage."""
        if self.length == 0:
            return 0.0

        gc_count = self.sequence.count('G') + self.sequence.count('C')
        return (gc_count / self.length) * 100

    @property
    def nucleotide_counts(self):
        """Get counts of each nucleotide."""
        return Counter(self.sequence)

    @property
    def composition(self):
        """Get nucleotide composition as percentages."""
        if self.length == 0:
            return {'A': 0, 'T': 0, 'C': 0, 'G': 0, 'N': 0}

        counts = self.nucleotide_counts
        return {nuc: (counts.get(nuc, 0) / self.length) * 100
                for nuc in 'ATCGN'}

    def complement(self):
        """Get the complement sequence."""
        return ''.join(self.COMPLEMENT_MAP.get(nuc, nuc) for nuc in self.sequence)

    def reverse_complement(self):
        """Get the reverse complement sequence."""
        return self.complement()[::-1]

    def translate(self, frame=0):
        """
        Translate DNA sequence to amino acid sequence.

        Args:
            frame (int): Reading frame (0, 1, or 2)

        Returns:
            str: Amino acid sequence
        """
        if frame not in (0, 1, 2):
            raise ValueError("Frame must be 0, 1, or 2")

        # Genetic code (simplified, standard code)
        genetic_code = {
            'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
            'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
            'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
            'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
            'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
            'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
            'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
            'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
            'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*',
            'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
            'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
            'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
            'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W',
            'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
            'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
            'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G'
        }

        protein = []
        seq = self.sequence[frame:]

        for i in range(0, len(seq) - 2, 3):
            codon = seq[i:i+3]
            if len(codon) == 3:
                aa = genetic_code.get(codon, 'X')  # X for unknown codons
                protein.append(aa)

        return ''.join(protein)

    def find_orfs(self, min_length=30):
        """
        Find Open Reading Frames (ORFs) in the sequence.

        Args:
            min_length (int): Minimum ORF length in nucleotides

        Returns:
            list: List of tuples (start, end, frame, protein_sequence)
        """
        orfs = []

        # Check all three reading frames
        for frame in range(3):
            protein = self.translate(frame)
            start_pos = None

            for i, aa in enumerate(protein):
                if aa == 'M' and start_pos is None:
                    # Start of ORF
                    start_pos = i
                elif aa == '*' and start_pos is not None:
                    # End of ORF
                    orf_length = (i - start_pos + 1) * 3
                    if orf_length >= min_length:
                        dna_start = frame + start_pos * 3
                        dna_end = frame + i * 3 + 3
                        orf_seq = protein[start_pos:i+1]
                        orfs.append((dna_start, dna_end, frame, orf_seq))
                    start_pos = None

            # Handle ORF that goes to end of sequence
            if start_pos is not None:
                orf_length = (len(protein) - start_pos) * 3
                if orf_length >= min_length:
                    dna_start = frame + start_pos * 3
                    dna_end = len(self.sequence)
                    orf_seq = protein[start_pos:]
                    orfs.append((dna_start, dna_end, frame, orf_seq))

        return orfs

    def get_subsequence(self, start, end):
        """
        Get a subsequence from start to end positions.

        Args:
            start (int): Start position (0-based)
            end (int): End position (exclusive)

        Returns:
            str: Subsequence
        """
        return self.sequence[start:end]

    def __str__(self):
        """String representation of the sequence."""
        return f"DNASequence({self.name}): {self.sequence[:50]}{'...' if len(self.sequence) > 50 else ''}"

    def __len__(self):
        """Length of the sequence."""
        return self.length


def parse_fasta(content):
    """
    Parse FASTA format content into DNASequence objects.

    Args:
        content (str): FASTA format string

    Returns:
        list: List of DNASequence objects

    Raises:
        ValueError: If FASTA format is invalid
    """
    sequences = []
    lines = content.strip().split('\n')

    current_name = ""
    current_description = ""
    current_sequence = []

    for line in lines:
        line = line.strip()
        if not line:
            continue

        if line.startswith('>'):
            # Save previous sequence if exists
            if current_sequence:
                seq_str = ''.join(current_sequence)
                sequences.append(DNASequence(seq_str, current_name, current_description))
                current_sequence = []

            # Parse header
            header = line[1:].split(None, 1)
            current_name = header[0]
            current_description = header[1] if len(header) > 1 else ""
        else:
            current_sequence.append(line)

    # Save last sequence
    if current_sequence:
        seq_str = ''.join(current_sequence)
        sequences.append(DNASequence(seq_str, current_name, current_description))

    if not sequences:
        raise ValueError("No valid FASTA sequences found")

    return sequences


def read_fasta_file(filepath):
    """
    Read and parse a FASTA file.

    Args:
        filepath (str): Path to FASTA file

    Returns:
        list: List of DNASequence objects
    """
    with open(filepath, 'r') as f:
        content = f.read()

    return parse_fasta(content)


def validate_sequence(sequence):
    """
    Validate a DNA sequence string.

    Args:
        sequence (str): Sequence to validate

    Returns:
        bool: True if valid, False otherwise
    """
    try:
        DNASequence(sequence)
        return True
    except ValueError:
        return False
