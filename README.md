# DNA Sequence Analyzer using DFA

## Finite Automata in Bioinformatics

This project demonstrates the application of **Deterministic Finite Automata (DFA)** from Theory of Computation in bioinformatics, specifically for DNA sequence analysis. DFA provides efficient exact pattern matching with guaranteed linear-time performance, making it ideal for searching biological motifs in genomic data.

### DFA Implementation Highlights

- **Pattern Matching via Finite Automata**: Each DNA pattern (e.g., restriction sites, start codons) is compiled into a DFA using the prefix table approach, ensuring O(n + m) time complexity where n is sequence length and m is pattern length.
- **KMP Algorithm Principles**: The DFA construction employs failure function computation similar to Knuth-Morris-Pratt algorithm, handling overlaps and preventing backtracking.
- **Optimization Techniques**: The graph-based DFA stores only necessary transitions, reducing space complexity while maintaining correctness.
- **Educational Visualization**: Interactive state diagrams show how automata process input symbols, providing insight into computational efficiency.

## Features

### Core Functionality
- **Exact Pattern Matching**: Search for predefined and custom DNA patterns using DFA (deterministic finite automaton)
- **Open Reading Frame (ORF) Detection**: Identify potential protein-coding regions in DNA sequences
- **Sequence Statistics**: GC content, nucleotide composition, sequence length analysis

### Predefined Patterns (TOC Examples)
- Start/Stop Codons: ATG, TAA, TAG, TGA
- Regulatory Elements: TATA box, Polyadenylation signals
- Restriction Enzymes: EcoRI (GAATTC), HindIII (AAGCTT)
- Splicing Signals: GT donor, AG acceptor sites

### Interfaces
- **Command-Line Tool**: `dna_analyzer.py` - Comprehensive analysis with multiple output formats
- **Web Application**: Interactive web interface with DFA state diagram visualization and JSON API endpoints

## Installation

### Prerequisites
- Python 3.7+
- Flask
- NetworkX

### Setup
```bash
# Clone the repository
git clone https://github.com/Xylecrack/dna-analyzer.git
cd dna-analyzer

# Install dependencies (ensure pip, flask, networkx are installed)
pip install flask networkx
```

## Usage

### Command Line Interface

#### Basic DNA Pattern Search
```bash
# Search for start codons in a sequence
python dna_analyzer.py -s ATGGATTACA -p start_codon

# Multiple patterns
python dna_analyzer.py -s GAATTCATGAAGCTT -p ecori_site hindiii_site

# Custom pattern
python dna_analyzer.py -s ATCGATCG -p "my_pattern:ATCG"
```

#### ORF Detection
```bash
# Find Open Reading Frames
python dna_analyzer.py -s ATGGAAATAG --find-orfs
```

#### Large File Analysis
```bash
# Process FASTA files >10MB with streaming
python dna_analyzer.py -f large_sequence.fasta --progress -p start_codon
```

#### Output Formats
```bash
# JSON output for programmatic use
python dna_analyzer.py -f sequence.fasta -p tata_box -o json
```

### Web Application

Run the web server:
```bash
python app.py
```
Navigate to `http://localhost:5000` for the interactive interface featuring:
- DFA state diagram visualization
- Real-time pattern matching
- File upload for large sequences
- Educational insights into automaton behavior

## Theory of Computation Educators

This implementation serves as a practical example of:

1. **Language Recognition**: DFA as recognizers for biological languages (DNA motifs)
2. **Complexity Analysis**: Demonstrating optimal algorithms in bioinformatics
3. **Automata Construction**: Prefix table computation and transition function building
4. **Pattern Matching Efficiency**: Linear time guarantees in practice
5. **State Diagram Interpretation**: Understanding computation through visualization

## Algorithm Details

The DFA construction algorithm:

1. **Validate Alphabet**: Ensure pattern contains only A, C, G, T
2. **Compute Failure Function**: KMP-style prefix table for state transitions
3. **Build Transition Graph**: Only store necessary edges for space efficiency
4. **Search Implementation**: Single pass through text with no backtracking

For pattern "ATG":
```
States: 0 → 1 → 2 → 3 (accept)
Transitions:
  0 -A→ 1
  1 -T→ 2
  2 -G→ 3
  Any invalid char from intermediate states → appropriate fallback via failure function
```

## Performance Characteristics

- **Time Complexity**: O(n) for search after O(m) preprocessing
- **Space Complexity**: O(m) for DFA storage
- **Streaming Support**: Processes arbitrarily large files in chunks
- **Memory Efficiency**: Minimal overhead for pattern matching vs regex alternatives

## Molecular Biology Context

- **Motif Detection**: Identify regulatory elements, binding sites, repeats
- **Gene Prediction**: ORF finding for potential protein-coding regions
- **Restriction Mapping**: Locate enzyme cutting sites
- **Sequence Annotation**: Automated classification and marking
- **Comparative Genomics**: Pattern conservation analysis
