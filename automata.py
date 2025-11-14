"""
Finite Automata for DNA Sequence Analysis

This module implements Deterministic Finite Automata (DFA) for exact pattern matching
in DNA sequences. The DFA is constructed using the Knuth-Morris-Pratt algorithm
principles for efficient pattern matching, and represented using NetworkX graphs.
"""

import networkx as nx

class DFA:
    """
    Deterministic Finite Automaton for exact pattern matching.

    Attributes:
        pattern (str): The pattern to search for
        alphabet (set): Set of valid characters (A, C, G, T)
        states (int): Number of states in the DFA
        graph (nx.DiGraph): NetworkX directed graph representing the DFA
        failure_function (list): Failure function for pattern matching
    """

    def __init__(self, pattern):
        """
        Initialize DFA for the given pattern.

        Args:
            pattern (str): The DNA pattern to search for (uppercase A,C,G,T)
        """
        self.pattern = pattern.upper()
        self.alphabet = {'A', 'C', 'G', 'T'}
        self.states = len(pattern) + 1
        self.graph = nx.DiGraph()
        self.failure_function = [0] * (len(pattern) + 1)

        self._validate_pattern()
        self._build_failure_function()
        self._build_graph()

    def _validate_pattern(self):
        """Validate that pattern contains only valid DNA nucleotides."""
        for char in self.pattern:
            if char not in self.alphabet:
                raise ValueError(f"Invalid character '{char}' in pattern. Only A, C, G, T allowed.")

    def _build_failure_function(self):
        """Build the failure function (prefix table) for the pattern."""
        length = 0
        i = 1

        while i < len(self.pattern):
            if self.pattern[i] == self.pattern[length]:
                length += 1
                self.failure_function[i] = length
                i += 1
            else:
                if length != 0:
                    length = self.failure_function[length - 1]
                else:
                    self.failure_function[i] = 0
                    i += 1

    def _build_graph(self):
        """Build the DFA as a NetworkX graph with only necessary transitions."""
        # Add all states as nodes
        for state in range(self.states):
            self.graph.add_node(state, is_accept=(state == len(self.pattern)))

        # Build only meaningful transitions
        # Track which transitions we've already added to avoid duplicates
        added_transitions = set()

        for state in range(self.states):
            for char in self.alphabet:
                next_state = self._get_next_state(state, char)

                # Only add transition if:
                # 1. It's a forward progress (matches pattern)
                # 2. It's a necessary failure transition (goes to state > 0)
                # 3. Not a self-loop to state 0 (reject state)
                edge_key = (state, next_state, char)

                if edge_key not in added_transitions:
                    # Add forward progress transitions (matching pattern characters)
                    if state < len(self.pattern) and char == self.pattern[state]:
                        self.graph.add_edge(state, next_state, label=char)
                        added_transitions.add(edge_key)
                    # Add necessary failure transitions (non-zero target states)
                    elif next_state > 0:
                        self.graph.add_edge(state, next_state, label=char)
                        added_transitions.add(edge_key)

    def _get_next_state(self, state, char):
        """
        Get the next state from current state on input character.

        Args:
            state (int): Current state
            char (str): Input character

        Returns:
            int: Next state
        """
        if state < len(self.pattern) and char == self.pattern[state]:
            return state + 1
        else:
            # Follow failure function
            if state == 0:
                return 0
            else:
                return self._get_next_state(self.failure_function[state - 1], char)

    def search(self, text):
        """
        Search for pattern occurrences in the text.

        Args:
            text (str): The DNA sequence to search in

        Returns:
            list: List of starting positions where pattern was found
        """
        text = text.upper()
        current_state = 0
        matches = []

        for i, char in enumerate(text):
            if char not in self.alphabet:
                # Skip invalid characters or reset state
                current_state = 0
                continue

            # Find next state using graph
            next_state = None
            for neighbor in self.graph.neighbors(current_state):
                edge_data = self.graph.get_edge_data(current_state, neighbor)
                if edge_data and edge_data.get('label') == char:
                    next_state = neighbor
                    break

            if next_state is not None:
                current_state = next_state
            else:
                current_state = 0  # Should not happen in a complete DFA

            if current_state == len(self.pattern):
                # Pattern found, record starting position
                matches.append(i - len(self.pattern) + 1)
                # Continue searching (don't reset to failure function for overlapping matches)
                current_state = self.failure_function[current_state - 1]

        return matches

    def get_state_diagram_data(self):
        """
        Get data for visualizing the state diagram.

        Returns:
            dict: Dictionary containing states, transitions, and accept states
        """
        states = list(self.graph.nodes())
        accept_states = [node for node in self.graph.nodes() if self.graph.nodes[node].get('is_accept', False)]
        transitions = []

        # Include all edges from the graph (including self-loops if present)
        for u, v, data in self.graph.edges(data=True):
            label = data.get('label', '')
            transitions.append({
                'from': u,
                'to': v,
                'label': label
            })

        return {
            'states': states,
            'accept_states': accept_states,
            'transitions': transitions,
            'pattern': self.pattern
        }

    def get_graph_properties(self):
        """
        Get various graph properties using NetworkX algorithms.

        Returns:
            dict: Dictionary containing graph analysis results
        """
        return {
            'is_strongly_connected': nx.is_strongly_connected(self.graph),
            'is_weakly_connected': nx.is_weakly_connected(self.graph),
            'has_cycles': len(list(nx.simple_cycles(self.graph))) > 0,
            'diameter': nx.diameter(self.graph) if nx.is_connected(self.graph.to_undirected()) else None,
            'average_shortest_path': nx.average_shortest_path_length(self.graph) if nx.is_connected(self.graph.to_undirected()) else None,
            'density': nx.density(self.graph),
            'number_of_edges': self.graph.number_of_edges(),
            'number_of_nodes': self.graph.number_of_nodes()
        }

    def get_shortest_path(self, start_state=0, end_state=None):
        """
        Find shortest path from start_state to end_state (or any accept state).

        Args:
            start_state (int): Starting state
            end_state (int, optional): Target state, defaults to accept state

        Returns:
            list: Shortest path as list of states, or None if no path
        """
        if end_state is None:
            accept_states = [node for node in self.graph.nodes() if self.graph.nodes[node].get('is_accept', False)]
            if not accept_states:
                return None
            end_state = accept_states[0]

        try:
            return nx.shortest_path(self.graph, start_state, end_state)
        except nx.NetworkXNoPath:
            return None

    def export_graphml(self, filename):
        """
        Export the DFA graph to GraphML format.

        Args:
            filename (str): Output filename
        """
        nx.write_graphml(self.graph, filename)


# Common DNA patterns
DNA_PATTERNS = {
    'start_codon': 'ATG',
    'stop_codon_taa': 'TAA',
    'stop_codon_tag': 'TAG',
    'stop_codon_tga': 'TGA',
    'tata_box': 'TATAAA',
    'ecori_site': 'GAATTC',
    'hindiii_site': 'AAGCTT',
    'polyadenylation_signal': 'AATAAA',
    'splice_donor': 'GT',
    'splice_acceptor': 'AG'
}


def create_dfa_for_pattern(pattern_name):
    """
    Create a DFA for a predefined DNA pattern.

    Args:
        pattern_name (str): Name of the pattern from DNA_PATTERNS

    Returns:
        DFA: Configured DFA for the pattern

    Raises:
        ValueError: If pattern_name is not recognized
    """
    if pattern_name not in DNA_PATTERNS:
        available = ', '.join(DNA_PATTERNS.keys())
        raise ValueError(f"Unknown pattern '{pattern_name}'. Available: {available}")

    return DFA(DNA_PATTERNS[pattern_name])
