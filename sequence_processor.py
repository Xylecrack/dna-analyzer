"""
Optimized DNA Sequence Processor for Large Sequences

Provides memory-efficient processing of large DNA sequences with streaming support,
progress tracking, and optimized algorithms for bioinformatics analysis.
"""

import os
import mmap
import logging
from typing import Iterator, Tuple, List, Optional, Callable
from dataclasses import dataclass
from collections import defaultdict
import time


@dataclass
class SequenceStats:
    """Statistics for a DNA sequence."""
    length: int
    gc_content: float
    composition: dict
    valid_chars: int
    invalid_chars: int


@dataclass
class MatchResult:
    """Result of a pattern match."""
    pattern: str
    position: int
    context: str = ""


@dataclass
class ProcessingProgress:
    """Progress information for sequence processing."""
    processed_bytes: int
    total_bytes: int
    matches_found: int
    time_elapsed: float
    current_operation: str


class LargeSequenceProcessor:
    """
    Memory-efficient processor for large DNA sequences.

    Supports streaming processing, chunked analysis, and progress tracking.
    """

    CHUNK_SIZE = 1024 * 1024  # 1MB chunks
    OVERLAP_SIZE = 1000  # Overlap between chunks to handle cross-chunk matches

    def __init__(self, chunk_size: int = CHUNK_SIZE, overlap_size: int = OVERLAP_SIZE):
        """
        Initialize the processor.

        Args:
            chunk_size: Size of chunks to process (bytes)
            overlap_size: Overlap between chunks (nucleotides)
        """
        self.chunk_size = chunk_size
        self.overlap_size = overlap_size
        self.valid_nucleotides = set('ATCGNatcgn')

    def process_file_streaming(self, file_path: str, patterns: List[str] = None,
                             progress_callback: Callable[[ProcessingProgress], None] = None) -> Iterator[Tuple[str, List[MatchResult]]]:
        """
        Process a large sequence file in streaming fashion.

        Args:
            file_path: Path to the sequence file
            patterns: List of patterns to search for
            progress_callback: Callback for progress updates

        Yields:
            Tuples of (chunk_sequence, matches_found)
        """
        file_size = os.path.getsize(file_path)
        processed_bytes = 0
        start_time = time.time()

        # Compile patterns into DFAs
        dfas = {}
        if patterns:
            from automata import DFA
            for pattern in patterns:
                if ':' in pattern:
                    name, seq = pattern.split(':', 1)
                    dfas[name] = DFA(seq.strip())
                else:
                    dfas[pattern] = DFA(pattern)

        with open(file_path, 'r', encoding='utf-8', errors='ignore') as f:
            buffer = ""
            chunk_start_pos = 0

            while True:
                chunk = f.read(self.chunk_size)
                if not chunk:
                    break

                # Clean and process the chunk
                cleaned_chunk = self._clean_sequence(chunk)
                buffer += cleaned_chunk

                # Process complete chunks
                if len(buffer) >= self.chunk_size + self.overlap_size:
                    # Find the last complete sequence position
                    process_end = len(buffer) - self.overlap_size

                    # Extract sequence to process
                    sequence_chunk = buffer[:process_end]

                    # Find matches in this chunk
                    matches = []
                    if patterns:
                        for pattern_name, dfa in dfas.items():
                            chunk_matches = dfa.search(sequence_chunk)
                            for pos in chunk_matches:
                                # Adjust position for global file position
                                global_pos = chunk_start_pos + pos
                                context_start = max(0, pos - 20)
                                context_end = min(len(sequence_chunk), pos + len(pattern_name) + 20)
                                context = sequence_chunk[context_start:context_end]

                                matches.append(MatchResult(
                                    pattern=pattern_name,
                                    position=global_pos,
                                    context=context
                                ))

                    # Update progress
                    processed_bytes += len(chunk)
                    if progress_callback:
                        progress = ProcessingProgress(
                            processed_bytes=processed_bytes,
                            total_bytes=file_size,
                            matches_found=len(matches),
                            time_elapsed=time.time() - start_time,
                            current_operation=f"Processing chunk at position {chunk_start_pos}"
                        )
                        progress_callback(progress)

                    yield sequence_chunk, matches

                    # Keep overlap for next iteration
                    buffer = buffer[process_end:]
                    chunk_start_pos += process_end

            # Process remaining buffer
            if buffer:
                matches = []
                if patterns:
                    for pattern_name, dfa in dfas.items():
                        chunk_matches = dfa.search(buffer)
                        for pos in chunk_matches:
                            global_pos = chunk_start_pos + pos
                            context_start = max(0, pos - 20)
                            context_end = min(len(buffer), pos + len(pattern_name) + 20)
                            context = buffer[context_start:context_end]

                            matches.append(MatchResult(
                                pattern=pattern_name,
                                position=global_pos,
                                context=context
                            ))

                yield buffer, matches

    def calculate_stats_streaming(self, file_path: str,
                                progress_callback: Callable[[ProcessingProgress], None] = None) -> SequenceStats:
        """
        Calculate sequence statistics in a streaming fashion.

        Args:
            file_path: Path to the sequence file
            progress_callback: Callback for progress updates

        Returns:
            SequenceStats object
        """
        file_size = os.path.getsize(file_path)
        processed_bytes = 0
        start_time = time.time()

        composition = defaultdict(int)
        total_length = 0
        valid_chars = 0
        invalid_chars = 0

        with open(file_path, 'r', encoding='utf-8', errors='ignore') as f:
            while True:
                chunk = f.read(self.chunk_size)
                if not chunk:
                    break

                # Process chunk
                cleaned_chunk = self._clean_sequence(chunk)
                total_length += len(cleaned_chunk)

                for char in cleaned_chunk:
                    if char in self.valid_nucleotides:
                        composition[char.upper()] += 1
                        valid_chars += 1
                    else:
                        invalid_chars += 1

                # Update progress
                processed_bytes += len(chunk)
                if progress_callback:
                    progress = ProcessingProgress(
                        processed_bytes=processed_bytes,
                        total_bytes=file_size,
                        matches_found=0,
                        time_elapsed=time.time() - start_time,
                        current_operation="Calculating statistics"
                    )
                    progress_callback(progress)

        # Calculate GC content
        gc_count = composition.get('G', 0) + composition.get('C', 0)
        gc_content = (gc_count / max(total_length, 1)) * 100

        return SequenceStats(
            length=total_length,
            gc_content=round(gc_content, 2),
            composition=dict(composition),
            valid_chars=valid_chars,
            invalid_chars=invalid_chars
        )

    def find_orfs_streaming(self, file_path: str, min_length: int = 30,
                           progress_callback: Callable[[ProcessingProgress], None] = None) -> List[Tuple[int, int, int, str]]:
        """
        Find Open Reading Frames in a large sequence file.

        Args:
            file_path: Path to the sequence file
            min_length: Minimum ORF length
            progress_callback: Callback for progress updates

        Returns:
            List of (start, end, frame, protein_sequence) tuples
        """
        file_size = os.path.getsize(file_path)
        processed_bytes = 0
        start_time = time.time()

        orfs = []
        sequence_buffer = ""
        global_position = 0

        with open(file_path, 'r', encoding='utf-8', errors='ignore') as f:
            while True:
                chunk = f.read(self.chunk_size)
                if not chunk:
                    break

                # Clean and add to buffer
                cleaned_chunk = self._clean_sequence(chunk)
                sequence_buffer += cleaned_chunk

                # Process complete codons (chunks of 3)
                process_length = (len(sequence_buffer) // 3) * 3

                if process_length >= 3:
                    chunk_sequence = sequence_buffer[:process_length]

                    # Check all three reading frames
                    for frame in range(3):
                        frame_sequence = chunk_sequence[frame:]
                        frame_orfs = self._find_orfs_in_frame(frame_sequence, min_length, global_position + frame)
                        orfs.extend([(start, end, frame, protein_seq) for start, end, protein_seq in frame_orfs])

                    # Keep remainder for next iteration
                    sequence_buffer = sequence_buffer[process_length:]
                    global_position += process_length

                # Update progress
                processed_bytes += len(chunk)
                if progress_callback:
                    progress = ProcessingProgress(
                        processed_bytes=processed_bytes,
                        total_bytes=file_size,
                        matches_found=len(orfs),
                        time_elapsed=time.time() - start_time,
                        current_operation="Finding ORFs"
                    )
                    progress_callback(progress)

        # Process remaining buffer
        if len(sequence_buffer) >= 3:
            for frame in range(3):
                frame_sequence = sequence_buffer[frame:]
                frame_orfs = self._find_orfs_in_frame(frame_sequence, min_length, global_position + frame)
                orfs.extend([(start, end, frame, protein_seq) for start, end, protein_seq in frame_orfs])

        return orfs

    def _clean_sequence(self, sequence: str) -> str:
        """Clean sequence by removing whitespace and converting to uppercase."""
        return ''.join(c.upper() for c in sequence if c.isalpha())

    def _find_orfs_in_frame(self, sequence: str, min_length: int, global_offset: int) -> List[Tuple[int, int, str]]:
        """Find ORFs in a single reading frame."""
        orfs = []
        i = 0

        # Genetic code for translation
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

        while i < len(sequence) - 2:
            codon = sequence[i:i+3]

            # Look for start codon
            if codon == 'ATG':
                start_pos = i
                j = i + 3

                # Look for stop codon
                while j < len(sequence) - 2:
                    stop_codon = sequence[j:j+3]
                    if stop_codon in ['TAA', 'TAG', 'TGA']:
                        orf_length = j + 3 - start_pos
                        if orf_length >= min_length:
                            # Translate the ORF to protein
                            orf_sequence = sequence[start_pos:j+3]
                            protein = []
                            for k in range(0, len(orf_sequence) - 2, 3):
                                codon = orf_sequence[k:k+3]
                                aa = genetic_code.get(codon, 'X')
                                protein.append(aa)
                            protein_seq = ''.join(protein)
                            orfs.append((global_offset + start_pos, global_offset + j + 3, protein_seq))
                        i = j + 2  # Skip to after stop codon
                        break
                    j += 3
                else:
                    # No stop codon found, check if ORF goes to end
                    orf_length = len(sequence) - start_pos
                    if orf_length >= min_length:
                        # Translate the ORF to protein
                        orf_sequence = sequence[start_pos:]
                        protein = []
                        for k in range(0, len(orf_sequence) - 2, 3):
                            codon = orf_sequence[k:k+3]
                            aa = genetic_code.get(codon, 'X')
                            protein.append(aa)
                        protein_seq = ''.join(protein)
                        orfs.append((global_offset + start_pos, global_offset + len(sequence), protein_seq))
                    i = len(sequence)  # End of sequence
            else:
                i += 3

        return orfs


class MemoryMappedSequence:
    """
    Memory-mapped sequence for extremely large files.

    Uses mmap for efficient random access to large sequence files.
    """

    def __init__(self, file_path: str):
        """
        Initialize memory-mapped sequence.

        Args:
            file_path: Path to the sequence file
        """
        self.file_path = file_path
        self.file_size = os.path.getsize(file_path)
        self._mmap = None
        self._sequence_cache = {}

    def __enter__(self):
        """Context manager entry."""
        self._mmap = mmap.mmap(os.open(self.file_path, os.O_RDONLY), 0, access=mmap.ACCESS_READ)
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        """Context manager exit."""
        if self._mmap:
            self._mmap.close()

    def get_sequence_chunk(self, start: int, length: int) -> str:
        """
        Get a chunk of sequence from the memory-mapped file.

        Args:
            start: Starting position
            length: Length of chunk to retrieve

        Returns:
            Sequence chunk as string
        """
        if not self._mmap:
            raise RuntimeError("MemoryMappedSequence must be used as context manager")

        # Read chunk from file
        self._mmap.seek(start)
        chunk = self._mmap.read(length).decode('utf-8', errors='ignore')

        # Clean the sequence
        return ''.join(c.upper() for c in chunk if c.isalpha())

    def find_pattern_efficient(self, pattern: str, start_pos: int = 0) -> Iterator[int]:
        """
        Efficiently find pattern occurrences using memory mapping.

        Args:
            pattern: Pattern to search for
            start_pos: Starting position in file

        Yields:
            Positions where pattern is found
        """
        if not self._mmap:
            raise RuntimeError("MemoryMappedSequence must be used as context manager")

        from automata import DFA
        dfa = DFA(pattern)

        chunk_size = 1024 * 1024  # 1MB chunks
        overlap = len(pattern) - 1

        pos = start_pos

        while pos < self.file_size:
            # Read chunk with overlap
            chunk_start = max(0, pos - overlap)
            chunk = self.get_sequence_chunk(chunk_start, chunk_size + overlap)

            if not chunk:
                break

            # Find matches in chunk
            matches = dfa.search(chunk)

            # Adjust positions and yield
            for match_pos in matches:
                global_pos = chunk_start + match_pos
                if global_pos >= start_pos:  # Ensure we don't yield positions before start_pos
                    yield global_pos

            # Move to next chunk
            pos += chunk_size

            # Break if we've processed the whole file
            if pos >= self.file_size:
                break
