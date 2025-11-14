#!/usr/bin/env python3
"""
DNA Sequence Analyzer - Command Line Interface

A tool for analyzing DNA sequences using finite automata.
Supports pattern matching, sequence statistics, and ORF finding.
"""

import argparse
import sys
import json
import os
from automata import DFA, create_dfa_for_pattern, DNA_PATTERNS
from sequence_utils import DNASequence, parse_fasta, read_fasta_file
from sequence_processor import LargeSequenceProcessor, SequenceStats, ProcessingProgress


def analyze_sequence(sequence, patterns=None, find_orfs=False, output_format='text'):
    """
    Analyze a DNA sequence with specified patterns.

    Args:
        sequence (DNASequence): The DNA sequence to analyze
        patterns (list): List of pattern names or custom patterns
        find_orfs (bool): Whether to find ORFs
        output_format (str): Output format ('text', 'json')

    Returns:
        dict: Analysis results
    """
    results = {
        'sequence_info': {
            'name': sequence.name,
            'length': sequence.length,
            'gc_content': round(sequence.gc_content, 2),
            'composition': {k: round(v, 2) for k, v in sequence.composition.items()}
        },
        'pattern_matches': {},
        'orfs': []
    }

    # Pattern matching
    if patterns:
        for pattern_spec in patterns:
            if ':' in pattern_spec:
                # Custom pattern: name:sequence
                name, pat_seq = pattern_spec.split(':', 1)
                try:
                    dfa = DFA(pat_seq)
                    matches = dfa.search(sequence.sequence)
                    results['pattern_matches'][name] = {
                        'pattern': pat_seq,
                        'matches': matches,
                        'count': len(matches)
                    }
                except ValueError as e:
                    results['pattern_matches'][name] = {'error': str(e)}
            else:
                # Predefined pattern
                try:
                    dfa = create_dfa_for_pattern(pattern_spec)
                    matches = dfa.search(sequence.sequence)
                    results['pattern_matches'][pattern_spec] = {
                        'pattern': DNA_PATTERNS[pattern_spec],
                        'matches': matches,
                        'count': len(matches)
                    }
                except ValueError as e:
                    results['pattern_matches'][pattern_spec] = {'error': str(e)}

    # ORF finding
    if find_orfs:
        orfs = sequence.find_orfs()
        results['orfs'] = [
            {
                'start': start,
                'end': end,
                'frame': frame,
                'length': end - start,
                'protein_length': len(protein_seq),
                'protein_sequence': protein_seq
            }
            for start, end, frame, protein_seq in orfs
        ]

    return results


def print_text_results(results):
    """Print results in human-readable text format."""
    seq_info = results['sequence_info']
    print(f"\n=== DNA Sequence Analysis ===")
    print(f"Name: {seq_info['name']}")
    print(f"Length: {seq_info['length']} bp")
    print(f"GC Content: {seq_info['gc_content']}%")
    print("Nucleotide Composition:")
    for nuc, pct in seq_info['composition'].items():
        print(f"  {nuc}: {pct}%")

    if results['pattern_matches']:
        print(f"\n=== Pattern Matches ===")
        for pattern_name, match_info in results['pattern_matches'].items():
            if 'error' in match_info:
                print(f"{pattern_name}: ERROR - {match_info['error']}")
            else:
                print(f"{pattern_name} ({match_info['pattern']}): {match_info['count']} matches")
                if match_info['matches']:
                    print(f"  Positions: {match_info['matches'][:10]}{'...' if len(match_info['matches']) > 10 else ''}")

    if results['orfs']:
        print(f"\n=== Open Reading Frames ===")
        print(f"Found {len(results['orfs'])} ORFs")
        for i, orf in enumerate(results['orfs'][:5], 1):  # Show first 5
            print(f"ORF {i}: {orf['start']}-{orf['end']} (frame {orf['frame']}), "
                  f"length {orf['length']}bp, protein {orf['protein_length']}aa")
        if len(results['orfs']) > 5:
            print(f"... and {len(results['orfs']) - 5} more")


def print_json_results(results):
    """Print results in JSON format."""
    print(json.dumps(results, indent=2))


def progress_callback(progress: ProcessingProgress):
    """Progress callback for streaming operations."""
    percent = (progress.processed_bytes / progress.total_bytes) * 100 if progress.total_bytes > 0 else 0
    print(f"\r[{progress.current_operation}] {progress.processed_bytes}/{progress.total_bytes} bytes "
          f"({percent:.1f}%) - {progress.matches_found} matches - {progress.time_elapsed:.1f}s", end='', flush=True)


def main():
    parser = argparse.ArgumentParser(
        description="DNA Sequence Analyzer using Finite Automata",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Analyze a sequence for start codons
  python dna_analyzer.py -s "ATGGATTACA" -p start_codon

  # Find multiple patterns
  python dna_analyzer.py -s "GAATTCATGAAGCTT" -p ecori_site hindiii_site

  # Custom pattern
  python dna_analyzer.py -s "ATCGATCG" -p "my_pattern:ATCG"

  # Analyze FASTA file
  python dna_analyzer.py -f sequence.fasta -p tata_box -o json

  # Find ORFs
  python dna_analyzer.py -s "ATGGAAATAG" --find-orfs

  # Large file processing
  python dna_analyzer.py -f large_sequence.fasta --large-file --progress -p start_codon

Available patterns: """ + ', '.join(DNA_PATTERNS.keys())
    )

    parser.add_argument('-s', '--sequence',
                       help='DNA sequence string')
    parser.add_argument('-f', '--fasta-file',
                       help='FASTA file path')
    parser.add_argument('-p', '--patterns', nargs='+',
                       help='Patterns to search for (predefined or custom:name:sequence)')
    parser.add_argument('--find-orfs', action='store_true',
                       help='Find Open Reading Frames')
    parser.add_argument('-o', '--output', choices=['text', 'json'],
                       default='text', help='Output format (default: text)')
    parser.add_argument('--list-patterns', action='store_true',
                       help='List available predefined patterns')
    parser.add_argument('--large-file', action='store_true',
                       help='Use streaming processing for large files (>10MB)')
    parser.add_argument('--chunk-size', type=int, default=1024*1024,
                       help='Chunk size for streaming processing (default: 1MB)')
    parser.add_argument('--progress', action='store_true',
                       help='Show progress during processing')
    parser.add_argument('--stats-only', action='store_true',
                       help='Calculate only sequence statistics (faster for large files)')

    args = parser.parse_args()

    if args.list_patterns:
        print("Available predefined patterns:")
        for name, pattern in DNA_PATTERNS.items():
            print(f"  {name}: {pattern}")
        return

    if not args.sequence and not args.fasta_file:
        parser.error("Either --sequence or --fasta-file must be provided")

    if args.sequence and args.fasta_file:
        parser.error("Cannot specify both --sequence and --fasta-file")

    try:
        # Check if we should use large file processing
        use_large_processor = args.large_file
        if args.fasta_file and not use_large_processor:
            # Auto-detect large files (>10MB)
            file_size = os.path.getsize(args.fasta_file)
            if file_size > 10 * 1024 * 1024:  # 10MB
                print(f"Large file detected ({file_size/1024/1024:.1f}MB). Using streaming processing.")
                use_large_processor = True

        if args.sequence:
            # Single sequence - use original method
            seq = DNASequence(args.sequence, "input_sequence")
            results = analyze_sequence(seq, args.patterns, args.find_orfs, args.output)
            if args.output == 'json':
                print_json_results(results)
            else:
                print_text_results(results)

        elif args.fasta_file:
            if use_large_processor:
                # Use streaming processor for large files
                processor = LargeSequenceProcessor(chunk_size=args.chunk_size)
                callback = progress_callback if args.progress else None

                if args.stats_only:
                    # Calculate only statistics
                    stats = processor.calculate_stats_streaming(args.fasta_file, callback)
                    results = {
                        'sequence_info': {
                            'name': os.path.basename(args.fasta_file),
                            'length': stats.length,
                            'gc_content': stats.gc_content,
                            'composition': {k: round((v/stats.length)*100, 2) for k, v in stats.composition.items()},
                            'valid_chars': stats.valid_chars,
                            'invalid_chars': stats.invalid_chars
                        }
                    }
                    if callback:
                        print()  # New line after progress
                    if args.output == 'json':
                        print_json_results(results)
                    else:
                        print_text_results(results)

                elif args.find_orfs and not args.patterns:
                    # Find ORFs only
                    orfs = processor.find_orfs_streaming(args.fasta_file, progress_callback=callback)
                    results = {
                        'sequence_info': {'name': os.path.basename(args.fasta_file)},
                        'orfs': [
                            {
                                'start': start,
                                'end': end,
                                'frame': frame,
                                'length': end - start
                            }
                            for start, end, frame in orfs
                        ]
                    }
                    if callback:
                        print()  # New line after progress
                    print(f"\nFound {len(orfs)} ORFs in large file")
                    if args.output == 'json':
                        print_json_results(results)

                else:
                    # Pattern matching with streaming
                    all_matches = {}
                    total_chunks = 0

                    for chunk_seq, matches in processor.process_file_streaming(args.fasta_file, args.patterns, callback):
                        total_chunks += 1
                        for match in matches:
                            pattern = match.pattern
                            if pattern not in all_matches:
                                all_matches[pattern] = []
                            all_matches[pattern].append({
                                'position': match.position,
                                'context': match.context
                            })

                    if callback:
                        print()  # New line after progress

                    # Summarize results
                    results = {
                        'sequence_info': {
                            'name': os.path.basename(args.fasta_file),
                            'processing_method': 'streaming',
                            'chunks_processed': total_chunks
                        },
                        'pattern_matches': {}
                    }

                    for pattern, matches in all_matches.items():
                        results['pattern_matches'][pattern] = {
                            'count': len(matches),
                            'matches': matches[:100] if len(matches) > 100 else matches,  # Limit for display
                            'truncated': len(matches) > 100
                        }

                    if args.output == 'json':
                        print_json_results(results)
                    else:
                        print(f"\n=== Large File Analysis ({os.path.basename(args.fasta_file)}) ===")
                        print(f"Processed in {total_chunks} chunks")
                        for pattern, match_info in results['pattern_matches'].items():
                            print(f"{pattern}: {match_info['count']} matches")
                            if match_info.get('truncated'):
                                print("  (showing first 100 matches)")

            else:
                # Use original method for smaller files
                sequences = read_fasta_file(args.fasta_file)
                all_results = []

                for seq in sequences:
                    result = analyze_sequence(seq, args.patterns, args.find_orfs, args.output)
                    all_results.append(result)

                if args.output == 'json':
                    print_json_results(all_results)
                else:
                    for i, results in enumerate(all_results):
                        if len(all_results) > 1:
                            print(f"\n--- Sequence {i+1} ---")
                        print_text_results(results)

    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == '__main__':
    main()
