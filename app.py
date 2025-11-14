"""
DNA Sequence Analyzer - Web Application

A Flask web application for analyzing DNA sequences using finite automata.
Provides a user-friendly interface for pattern matching and sequence analysis.
"""

from flask import Flask, render_template, request, jsonify
import json
import os
import tempfile
from werkzeug.utils import secure_filename
from automata import DFA, create_dfa_for_pattern, DNA_PATTERNS
from sequence_utils import DNASequence, parse_fasta, read_fasta_file
from sequence_processor import LargeSequenceProcessor, ProcessingProgress

app = Flask(__name__)
app.config['MAX_CONTENT_LENGTH'] = 500 * 1024 * 1024  # 500MB max file size
app.config['UPLOAD_FOLDER'] = tempfile.gettempdir()

# Logging file path
LOG_FILE = 'dfa_analysis_log.txt'

def log_analysis(user_input, pattern, states, transitions):
    """
    Log user input and DFA transition details.

    Args:
        user_input (str): User's DNA sequence input
        pattern (str): Pattern being analyzed
        states (list): List of DFA states
        transitions (list): List of transitions with labels
    """
    try:
        with open(LOG_FILE, 'a', encoding='utf-8') as f:
            f.write('='*80 + '\n')
            f.write(f'USER INPUT: {user_input}\n')
            f.write(f'PATTERN: {pattern}\n')
            f.write(f'STATES: {states}\n')
            f.write('\nTRANSITIONS:\n')

            # Categorize transitions
            forward_transitions = []
            self_loops = []
            backward_transitions = []

            for trans in transitions:
                from_state = trans['from']
                to_state = trans['to']
                label = trans['label']

                if from_state == to_state:
                    self_loops.append(trans)
                    f.write(f"  [SELF-LOOP] State {from_state} --{label}--> State {to_state}\n")
                elif to_state > from_state:
                    forward_transitions.append(trans)
                    if from_state == 0:
                        f.write(f"  [START] State {from_state} --{label}--> State {to_state}\n")
                    elif to_state == len(states) - 1:
                        f.write(f"  [END] State {from_state} --{label}--> State {to_state}\n")
                    else:
                        f.write(f"  [FORWARD] State {from_state} --{label}--> State {to_state}\n")
                else:
                    backward_transitions.append(trans)
                    f.write(f"  [BACKWARD] State {from_state} --{label}--> State {to_state}\n")

            f.write(f'\nSUMMARY:\n')
            f.write(f'  Total States: {len(states)}\n')
            f.write(f'  Total Transitions: {len(transitions)}\n')
            f.write(f'  Forward Transitions: {len(forward_transitions)}\n')
            f.write(f'  Self-Loops: {len(self_loops)}\n')
            f.write(f'  Backward Transitions: {len(backward_transitions)}\n')
            f.write('='*80 + '\n\n')
    except Exception as e:
        print(f"Logging error: {e}")

@app.route('/')
def index():
    """Render the main page."""
    return render_template('index.html', patterns=DNA_PATTERNS)

@app.route('/upload', methods=['POST'])
def upload_file():
    """
    Handle file upload for large sequence files.

    Returns:
        JSON response with upload status
    """
    try:
        if 'file' not in request.files:
            return jsonify({'success': False, 'error': 'No file provided'})

        file = request.files['file']
        if file.filename == '':
            return jsonify({'success': False, 'error': 'No file selected'})

        if file:
            filename = secure_filename(file.filename)
            file_path = os.path.join(app.config['UPLOAD_FOLDER'], filename)
            file.save(file_path)

            # Check file size
            file_size = os.path.getsize(file_path)
            return jsonify({
                'success': True,
                'filename': filename,
                'file_path': file_path,
                'file_size': file_size,
                'file_size_mb': round(file_size / (1024 * 1024), 2)
            })

    except Exception as e:
        return jsonify({'success': False, 'error': str(e)})


@app.route('/analyze', methods=['POST'])
def analyze():
    """
    Analyze DNA sequence based on form submission.

    Returns:
        JSON response with analysis results
    """
    try:
        # Check if this is a file upload analysis
        uploaded_file = request.form.get('uploaded_file', '')
        use_large_processor = request.form.get('use_large_processor') == 'true'

        if uploaded_file and os.path.exists(uploaded_file):
            # Handle large file processing
            file_size = os.path.getsize(uploaded_file)
            if file_size > 10 * 1024 * 1024 or use_large_processor:  # >10MB or explicitly requested
                return analyze_large_file(uploaded_file)

        # Get form data for regular processing
        sequence_input = request.form.get('sequence', '').strip()
        input_type = request.form.get('input_type', 'sequence')
        selected_patterns = request.form.getlist('patterns')
        custom_patterns = request.form.get('custom_patterns', '').strip()
        find_orfs = request.form.get('find_orfs') == 'on'

        # Parse sequences
        if input_type == 'fasta':
            sequences = parse_fasta(sequence_input)
        else:
            sequences = [DNASequence(sequence_input, "input_sequence")]

        results = []

        for seq in sequences:
            result = {
                'sequence_info': {
                    'name': seq.name,
                    'length': seq.length,
                    'gc_content': round(seq.gc_content, 2),
                    'composition': {k: round(v, 2) for k, v in seq.composition.items()}
                },
                'pattern_matches': {},
                'orfs': []
            }

            # Process selected patterns
            for pattern_name in selected_patterns:
                try:
                    dfa = DFA(DNA_PATTERNS[pattern_name])
                    matches = dfa.search(seq.sequence)
                    result['pattern_matches'][pattern_name] = {
                        'pattern': DNA_PATTERNS[pattern_name],
                        'matches': matches,
                        'count': len(matches)
                    }
                except ValueError as e:
                    result['pattern_matches'][pattern_name] = {'error': str(e)}

            # Process custom patterns
            if custom_patterns:
                for line in custom_patterns.split('\n'):
                    line = line.strip()
                    if line and ':' in line:
                        name, pat_seq = line.split(':', 1)
                        name = name.strip()
                        pat_seq = pat_seq.strip()
                        try:
                            dfa = DFA(pat_seq)
                            matches = dfa.search(seq.sequence)
                            result['pattern_matches'][name] = {
                                'pattern': pat_seq,
                                'matches': matches,
                                'count': len(matches)
                            }
                        except ValueError as e:
                            result['pattern_matches'][name] = {'error': str(e)}

            # Find ORFs if requested
            if find_orfs:
                orfs = seq.find_orfs(min_length=6)  # Minimum 6bp (2 codons) for ORFs
                result['orfs'] = [
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

            results.append(result)

        return jsonify({'success': True, 'results': results})

    except Exception as e:
        return jsonify({'success': False, 'error': str(e)})


def analyze_large_file(file_path):
    """
    Analyze a large file using streaming processor.

    Args:
        file_path (str): Path to the uploaded file

    Returns:
        JSON response with analysis results
    """
    try:
        # Get analysis parameters from form
        selected_patterns = request.form.getlist('patterns')
        custom_patterns = request.form.get('custom_patterns', '').strip()
        find_orfs = request.form.get('find_orfs') == 'on'
        stats_only = request.form.get('stats_only') == 'true'

        # Combine patterns
        all_patterns = selected_patterns[:]
        if custom_patterns:
            for line in custom_patterns.split('\n'):
                line = line.strip()
                if line and ':' in line:
                    name, _ = line.split(':', 1)
                    all_patterns.append(name.strip())

        processor = LargeSequenceProcessor()

        if stats_only:
            # Calculate only statistics
            stats = processor.calculate_stats_streaming(file_path)
            results = [{
                'sequence_info': {
                    'name': os.path.basename(file_path),
                    'length': stats.length,
                    'gc_content': stats.gc_content,
                    'composition': {k: round((v/stats.length)*100, 2) for k, v in stats.composition.items()},
                    'valid_chars': stats.valid_chars,
                    'invalid_chars': stats.invalid_chars,
                    'processing_method': 'streaming_stats'
                }
            }]
            return jsonify({'success': True, 'results': results})

        elif find_orfs and not all_patterns:
            # Find ORFs only
            orfs = processor.find_orfs_streaming(file_path, min_length=6)
            results = [{
                'sequence_info': {'name': os.path.basename(file_path)},
                'orfs': [
                    {
                        'start': start,
                        'end': end,
                        'frame': frame,
                        'length': end - start,
                        'protein_length': len(protein_seq),
                        'protein_sequence': protein_seq
                    }
                    for start, end, frame, protein_seq in orfs[:100]  # Limit for web display
                ],
                'total_orfs': len(orfs),
                'truncated': len(orfs) > 100
            }]
            return jsonify({'success': True, 'results': results})

        else:
            # Pattern matching with streaming
            all_matches = {}
            total_chunks = 0

            for chunk_seq, matches in processor.process_file_streaming(file_path, all_patterns if all_patterns else None):
                total_chunks += 1
                for match in matches:
                    pattern = match.pattern
                    if pattern not in all_matches:
                        all_matches[pattern] = []
                    all_matches[pattern].append({
                        'position': match.position,
                        'context': match.context
                    })

            # Build results
            results = [{
                'sequence_info': {
                    'name': os.path.basename(file_path),
                    'processing_method': 'streaming',
                    'chunks_processed': total_chunks
                },
                'pattern_matches': {}
            }]

            for pattern, matches in all_matches.items():
                results[0]['pattern_matches'][pattern] = {
                    'count': len(matches),
                    'matches': matches[:50] if len(matches) > 50 else matches,  # Limit for web display
                    'truncated': len(matches) > 50
                }

            return jsonify({'success': True, 'results': results})

    except Exception as e:
        return jsonify({'success': False, 'error': str(e)})

@app.route('/state_diagram/<pattern_name>')
def state_diagram(pattern_name):
    """
    Get state diagram data for a pattern.

    Args:
        pattern_name (str): Name of the predefined pattern

    Returns:
        JSON response with state diagram data
    """
    try:
        dfa = create_dfa_for_pattern(pattern_name)
        diagram_data = dfa.get_state_diagram_data()

        # Log the diagram generation
        log_analysis(
            user_input=f"Selected Pattern: {pattern_name}",
            pattern=diagram_data['pattern'],
            states=diagram_data['states'],
            transitions=diagram_data['transitions']
        )

        return jsonify({'success': True, 'diagram': diagram_data})
    except ValueError as e:
        return jsonify({'success': False, 'error': str(e)})

@app.route('/custom_state_diagram', methods=['POST'])
def custom_state_diagram():
    """
    Get state diagram data for a custom pattern.

    Returns:
        JSON response with state diagram data
    """
    try:
        pattern = request.json.get('pattern', '')
        dfa = DFA(pattern)
        diagram_data = dfa.get_state_diagram_data()

        # Log the custom pattern analysis
        log_analysis(
            user_input=f"Custom Pattern Input: {pattern}",
            pattern=diagram_data['pattern'],
            states=diagram_data['states'],
            transitions=diagram_data['transitions']
        )

        return jsonify({'success': True, 'diagram': diagram_data})
    except ValueError as e:
        return jsonify({'success': False, 'error': str(e)})

if __name__ == '__main__':
    app.run(debug=True, host='0.0.0.0', port=5000)
