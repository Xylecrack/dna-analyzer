/**
 * DFA State Diagram Visualization
 * Professional, clean, and interactive graph visualization for Deterministic Finite Automata
 */

class DFAVisualizer {
    constructor(containerId, width = 800, height = 600) {
        this.containerId = containerId;
        this.width = width;
        this.height = height;
        this.svg = null;
        this.zoom = null;
        this.data = null;
        this.nodes = [];
        this.links = [];
        this.nodeRadius = 30;
        this.selfLoopRadius = 25;
        this.centerX = width / 2;
        this.centerY = height / 2;

        this.init();
    }

    init() {
        // Create SVG container with better sizing
        this.svg = d3.select(`#${this.containerId}`)
            .append('svg')
            .attr('width', this.width)
            .attr('height', this.height)
            .attr('viewBox', [0, 0, this.width, this.height])
            .attr('style', 'max-width: 100%; height: auto; border: 1px solid #e9ecef; border-radius: 8px;');

        // Add zoom behavior with better controls
        this.zoom = d3.zoom()
            .scaleExtent([0.3, 3])
            .on('zoom', (event) => {
                this.svg.select('.graph-group').attr('transform', event.transform);
            });

        this.svg.call(this.zoom);

        // Create main group for graph elements
        this.svg.append('g').attr('class', 'graph-group');

        // Add arrow markers for transitions
        const defs = this.svg.append('defs');

        // End arrow (for forward direction) - larger and more visible
        defs.append('marker')
            .attr('id', 'arrowhead-end')
            .attr('viewBox', '0 -5 10 10')
            .attr('refX', 28)
            .attr('refY', 0)
            .attr('markerWidth', 10)
            .attr('markerHeight', 10)
            .attr('orient', 'auto')
            .append('path')
            .attr('d', 'M0,-5L10,0L0,5')
            .attr('fill', '#495057');

        // Start arrow (for reverse direction on bidirectional edges) - larger
        defs.append('marker')
            .attr('id', 'arrowhead-start')
            .attr('viewBox', '0 -5 10 10')
            .attr('refX', -18)
            .attr('refY', 0)
            .attr('markerWidth', 10)
            .attr('markerHeight', 10)
            .attr('orient', 'auto')
            .append('path')
            .attr('d', 'M10,-5L0,0L10,5')
            .attr('fill', '#495057');

        // Highlighted end arrow
        defs.append('marker')
            .attr('id', 'arrowhead-end-highlight')
            .attr('viewBox', '0 -5 10 10')
            .attr('refX', 28)
            .attr('refY', 0)
            .attr('markerWidth', 10)
            .attr('markerHeight', 10)
            .attr('orient', 'auto')
            .append('path')
            .attr('d', 'M0,-5L10,0L0,5')
            .attr('fill', '#ff6b35');

        // Start state arrow
        defs.append('marker')
            .attr('id', 'start-arrow')
            .attr('viewBox', '0 -5 10 10')
            .attr('refX', 5)
            .attr('refY', 0)
            .attr('markerWidth', 6)
            .attr('markerHeight', 6)
            .attr('orient', 'auto')
            .append('path')
            .attr('d', 'M10,-5L0,0L10,5')
            .attr('fill', '#007bff');
    }

    render(data) {
        this.data = data;
        this.prepareData(data);
        this.layoutCircular();
        this.renderEdges();
        this.renderNodes();
        this.renderLabels();
        this.addPatternInfo();
        this.setupInteractions();

        return this;
    }

    prepareData(data) {
        // Prepare nodes
        this.nodes = data.states.map(state => ({
            id: state,
            isAccept: data.accept_states.includes(state),
            label: state.toString(),
            x: 0,
            y: 0
        }));

        // Group transitions by source-target pairs
        const linkGroups = new Map();
        data.transitions.forEach((transition) => {
            const key = `${transition.from}-${transition.to}`;
            if (!linkGroups.has(key)) {
                linkGroups.set(key, []);
            }
            linkGroups.get(key).push(transition);
        });

        // Create links with positioning information
        this.links = [];
        linkGroups.forEach(group => {
            group.forEach((transition, groupIndex) => {
                this.links.push({
                    source: transition.from,
                    target: transition.to,
                    label: transition.label,
                    isSelfLoop: transition.from === transition.to,
                    groupIndex: groupIndex,
                    totalInGroup: group.length,
                    originalTransition: transition
                });
            });
        });

        // Detect and merge bidirectional edges
        this.mergeBidirectionalEdges();
    }

    mergeBidirectionalEdges() {
        const mergedLinks = [];
        const processedPairs = new Set();

        this.links.forEach(link => {
            if (link.isSelfLoop) {
                mergedLinks.push(link);
                return;
            }

            const pairKey = link.source < link.target ?
                `${link.source}-${link.target}` :
                `${link.target}-${link.source}`;

            if (processedPairs.has(pairKey)) {
                return; // Already processed this bidirectional pair
            }

            // Check if there's a reverse edge
            const reverseLink = this.links.find(l =>
                l.source === link.target &&
                l.target === link.source &&
                !l.isSelfLoop
            );

            if (reverseLink) {
                // Create a bidirectional edge
                mergedLinks.push({
                    source: link.source,
                    target: link.target,
                    reverseSource: reverseLink.source,
                    reverseTarget: reverseLink.target,
                    label: link.label,
                    reverseLabel: reverseLink.label,
                    isBidirectional: true,
                    isSelfLoop: false,
                    groupIndex: 0,
                    totalInGroup: 1,
                    originalTransition: link.originalTransition,
                    reverseTransition: reverseLink.originalTransition
                });
                processedPairs.add(pairKey);
            } else {
                // Single direction edge
                mergedLinks.push(link);
            }
        });

        this.links = mergedLinks;
    }

    layoutCircular() {
        const numStates = this.nodes.length;
        // Auto-adjust radius based on number of states to prevent crowding
        let baseRadius = Math.min(this.width, this.height) * 0.35;

        // Increase radius for more states
        if (numStates > 6) {
            baseRadius = Math.min(this.width, this.height) * 0.4;
        }
        if (numStates > 10) {
            baseRadius = Math.min(this.width, this.height) * 0.42;
        }

        this.nodes.forEach((node, i) => {
            const angle = (2 * Math.PI * i) / numStates - Math.PI / 2; // Start from top
            node.x = this.centerX + baseRadius * Math.cos(angle);
            node.y = this.centerY + baseRadius * Math.sin(angle);
        });
    }

    renderEdges() {
        const graphGroup = this.svg.select('.graph-group');

        // For bidirectional edges, we need two path elements
        const allPaths = [];
        this.links.forEach(link => {
            if (link.isBidirectional) {
                // Create two separate paths for bidirectional edges
                allPaths.push({
                    ...link,
                    direction: 'forward',
                    pathId: `${link.source}-${link.target}-forward`
                });
                allPaths.push({
                    ...link,
                    direction: 'reverse',
                    pathId: `${link.source}-${link.target}-reverse`
                });
            } else {
                allPaths.push({
                    ...link,
                    pathId: `${link.source}-${link.target}-${link.label}`
                });
            }
        });

        const linkElements = graphGroup.selectAll('.dfa-link')
            .data(allPaths, d => d.pathId)
            .join('path')
            .attr('class', 'dfa-link')
            .attr('marker-end', d => {
                if (d.isBidirectional) {
                    return d.direction === 'forward' ? 'url(#arrowhead-end)' : null;
                }
                return 'url(#arrowhead-end)';
            })
            .attr('marker-start', d => {
                if (d.isBidirectional) {
                    return d.direction === 'reverse' ? 'url(#arrowhead-start)' : null;
                }
                return null;
            })
            .attr('stroke', '#495057')
            .attr('stroke-width', 2.5)
            .attr('fill', 'none')
            .style('opacity', 0.85)
            .attr('d', d => this.getEdgePath(d));
    }

    getEdgePath(link) {
        const sourceNode = this.nodes.find(n => n.id === link.source);
        const targetNode = this.nodes.find(n => n.id === link.target);

        if (!sourceNode || !targetNode) return '';

        if (link.isSelfLoop) {
            // Self-loop tangent to the node at the right side
            const x = sourceNode.x;
            const y = sourceNode.y;
            const r = 30; // Radius for the loop
            // Start and end at node boundary, curve out and back
            return `M ${x + this.nodeRadius} ${y} C ${x + this.nodeRadius + r} ${y - r}, ${x + this.nodeRadius + r} ${y + r}, ${x + this.nodeRadius} ${y}`;
        } else {
            const dx = targetNode.x - sourceNode.x;
            const dy = targetNode.y - sourceNode.y;
            const distance = Math.sqrt(dx * dx + dy * dy);

            if (distance === 0) return '';

            // For bidirectional edges, use gentle curves instead of parallel lines
            if (link.isBidirectional) {
                // Calculate perpendicular direction for offset
                const perpX = -dy / distance;
                const perpY = dx / distance;
                const offset = 18; // Increased separation between curves

                // Calculate the start and end points at node boundaries
                const startX = sourceNode.x + (dx / distance) * this.nodeRadius;
                const startY = sourceNode.y + (dy / distance) * this.nodeRadius;
                const endX = targetNode.x - (dx / distance) * this.nodeRadius;
                const endY = targetNode.y - (dy / distance) * this.nodeRadius;

                // Midpoint for curve control
                const midX = (startX + endX) / 2;
                const midY = (startY + endY) / 2;

                if (link.direction === 'forward') {
                    // Upper curve with offset control point
                    const ctrlX = midX + perpX * offset;
                    const ctrlY = midY + perpY * offset;
                    return `M ${startX} ${startY} Q ${ctrlX} ${ctrlY} ${endX} ${endY}`;
                } else {
                    // Lower curve with opposite offset - draw from target to source
                    const ctrlX = midX - perpX * offset;
                    const ctrlY = midY - perpY * offset;
                    return `M ${endX} ${endY} Q ${ctrlX} ${ctrlY} ${startX} ${startY}`;
                }
            } else {
                // Regular curved edge for single-direction transitions
                const controlOffset = Math.min(distance * 0.2, 50);
                const midX = (sourceNode.x + targetNode.x) / 2;
                const midY = (sourceNode.y + targetNode.y) / 2;

                // Perpendicular offset for multiple edges between same nodes
                if (link.totalInGroup > 1) {
                    const perpX = -dy / distance;
                    const perpY = dx / distance;
                    const offset = (link.groupIndex - (link.totalInGroup - 1) / 2) * 30;
                    const offsetX = perpX * offset;
                    const offsetY = perpY * offset;

                    return `M ${sourceNode.x} ${sourceNode.y} Q ${midX + offsetX} ${midY + offsetY} ${targetNode.x} ${targetNode.y}`;
                } else {
                    // Slight curve for better visibility
                    return `M ${sourceNode.x} ${sourceNode.y} Q ${midX} ${midY} ${targetNode.x} ${targetNode.y}`;
                }
            }
        }
    }

    renderNodes() {
        const graphGroup = this.svg.select('.graph-group');

        // Create node groups
        const nodeGroups = graphGroup.selectAll('.node-group')
            .data(this.nodes)
            .enter().append('g')
            .attr('class', 'node-group')
            .attr('transform', d => `translate(${d.x}, ${d.y})`);

        // Main node circles
        nodeGroups.append('circle')
            .attr('class', 'dfa-node')
            .attr('r', this.nodeRadius)
            .attr('fill', d => d.isAccept ? '#28a745' : '#007bff')
            .attr('stroke', '#ffffff')
            .attr('stroke-width', 3)
            .style('filter', 'drop-shadow(0 2px 4px rgba(0,0,0,0.1))');

        // Double ring for accept states - more visible
        nodeGroups.filter(d => d.isAccept)
            .append('circle')
            .attr('r', this.nodeRadius - 5)
            .attr('fill', 'none')
            .attr('stroke', '#ffffff')
            .attr('stroke-width', 2.5);

        // Add state number labels inside each node
        nodeGroups.append('text')
            .attr('class', 'node-label')
            .attr('text-anchor', 'middle')
            .attr('dy', 6)
            .attr('font-size', '18px')
            .attr('font-weight', 'bold')
            .attr('fill', 'white')
            .attr('pointer-events', 'none')
            .text(d => d.id);

        // Add "START" indicator for start node
        nodeGroups.filter(d => d.id === 0)
            .append('text')
            .attr('class', 'start-indicator')
            .attr('text-anchor', 'middle')
            .attr('dy', -this.nodeRadius - 12)
            .attr('font-size', '11px')
            .attr('font-weight', 'bold')
            .attr('fill', '#007bff')
            .attr('pointer-events', 'none')
            .text('START');

        // Add tooltips
        nodeGroups.append('title')
            .text(d => {
                if (d.id === 0 && d.isAccept) return `Start & Accept State ${d.id}`;
                if (d.id === 0) return `Start State ${d.id}`;
                if (d.isAccept) return `Accept State ${d.id}`;
                return `State ${d.id}`;
            });
    }

    renderLabels() {
        const graphGroup = this.svg.select('.graph-group');

        // Prepare label data (including bidirectional labels)
        const labelData = this.prepareLabelData();

        // Add white background rectangles for labels (for better readability)
        const labelBgs = graphGroup.selectAll('.edge-label-bg')
            .data(labelData)
            .enter().append('rect')
            .attr('class', 'edge-label-bg')
            .attr('fill', 'rgba(255, 255, 255, 0.95)')
            .attr('stroke', 'rgba(0, 0, 0, 0.1)')
            .attr('stroke-width', 1)
            .attr('rx', 4)
            .attr('ry', 4);

        // Edge labels with collision detection
        const edgeLabels = graphGroup.selectAll('.edge-label')
            .data(labelData)
            .enter().append('text')
            .attr('class', 'edge-label')
            .attr('text-anchor', 'middle')
            .attr('font-size', '14px')
            .attr('font-weight', '600')
            .attr('fill', '#212529')
            .attr('pointer-events', 'none')
            .text(d => d.label);

        // Position labels with collision avoidance
        this.positionLabelsWithCollisionAvoidance(edgeLabels);

        // Update background rectangles to match text size
        edgeLabels.each(function(d, i) {
            const bbox = this.getBBox();
            const bg = d3.select(labelBgs.nodes()[i]);
            bg.attr('x', bbox.x - 4)
              .attr('y', bbox.y - 2)
              .attr('width', bbox.width + 8)
              .attr('height', bbox.height + 4);
        });
    }

    prepareLabelData() {
        const labelData = [];

        this.links.forEach((link, linkIndex) => {
            if (link.isBidirectional) {
                // For bidirectional edges, create two labels
                const pos1 = this.getBidirectionalLabelPosition(link, true);
                const pos2 = this.getBidirectionalLabelPosition(link, false);

                labelData.push({
                    link: link,
                    label: link.label,
                    x: pos1.x,
                    y: pos1.y,
                    isForward: true,
                    index: labelData.length
                });

                labelData.push({
                    link: link,
                    label: link.reverseLabel,
                    x: pos2.x,
                    y: pos2.y,
                    isForward: false,
                    index: labelData.length
                });
            } else {
                // For regular edges, create one label
                const pos = this.getLabelPosition(link);
                labelData.push({
                    link: link,
                    label: link.label,
                    x: pos.x,
                    y: pos.y,
                    index: labelData.length
                });
            }
        });

        return labelData;
    }

    getBidirectionalLabelPosition(link, isForward) {
        // For bidirectional edges, both directions use the same source/target
        const sourceNode = this.nodes.find(n => n.id === link.source);
        const targetNode = this.nodes.find(n => n.id === link.target);

        if (!sourceNode || !targetNode) return { x: 0, y: 0 };

        // Position along the curved edge
        const dx = targetNode.x - sourceNode.x;
        const dy = targetNode.y - sourceNode.y;
        const distance = Math.sqrt(dx * dx + dy * dy);

        if (distance === 0) return { x: sourceNode.x, y: sourceNode.y };

        const midX = (sourceNode.x + targetNode.x) / 2;
        const midY = (sourceNode.y + targetNode.y) / 2;

        // Position labels along the curves with increased separation
        const offset = 22; // Increased distance from center line
        const perpX = -dy / distance;
        const perpY = dx / distance;

        if (isForward) {
            return {
                x: midX + perpX * offset,
                y: midY + perpY * offset
            };
        } else {
            return {
                x: midX - perpX * offset,
                y: midY - perpY * offset
            };
        }
    }

    getLabelPosition(link) {
        const sourceNode = this.nodes.find(n => n.id === link.source);
        const targetNode = this.nodes.find(n => n.id === link.target);

        if (!sourceNode || !targetNode) return { x: 0, y: 0 };

        if (link.isSelfLoop) {
            // Position to the right of the node to avoid collision with other edges
            return {
                x: sourceNode.x + this.nodeRadius + 25,
                y: sourceNode.y
            };
        } else {
            // Position along the curved edge
            const dx = targetNode.x - sourceNode.x;
            const dy = targetNode.y - sourceNode.y;
            const distance = Math.sqrt(dx * dx + dy * dy);

            // Control point for curve
            const midX = (sourceNode.x + targetNode.x) / 2;
            const midY = (sourceNode.y + targetNode.y) / 2;

            // Perpendicular offset for multiple edges
            if (link.totalInGroup > 1) {
                const perpX = -dy / distance;
                const perpY = dx / distance;
                const offset = (link.groupIndex - (link.totalInGroup - 1) / 2) * 30;
                return {
                    x: midX + perpX * offset,
                    y: midY + perpY * offset
                };
            } else {
                return { x: midX, y: midY };
            }
        }
    }

    addStartArrow() {
        const graphGroup = this.svg.select('.graph-group');
        const startNode = this.nodes.find(n => n.id === 0);

        if (startNode) {
            // Calculate position for start arrow (left of start state)
            const arrowLength = 60;
            const arrowX = startNode.x - this.nodeRadius - arrowLength;
            const arrowY = startNode.y;

            // Draw arrow line
            graphGroup.append('line')
                .attr('x1', arrowX)
                .attr('y1', arrowY)
                .attr('x2', startNode.x - this.nodeRadius - 5)
                .attr('y2', arrowY)
                .attr('stroke', '#007bff')
                .attr('stroke-width', 3)
                .attr('marker-end', 'url(#start-arrow)');

            // Add "Start" label
            graphGroup.append('text')
                .attr('x', arrowX - 10)
                .attr('y', arrowY - 5)
                .attr('text-anchor', 'end')
                .attr('font-size', '12px')
                .attr('font-weight', 'bold')
                .attr('fill', '#007bff')
                .text('Start');
        }
    }

    addPatternInfo() {
        const graphGroup = this.svg.select('.graph-group');

        if (this.data.pattern) {
            // Position pattern info at the top center
            graphGroup.append('text')
                .attr('x', this.centerX)
                .attr('y', 30)
                .attr('text-anchor', 'middle')
                .attr('font-size', '18px')
                .attr('font-weight', 'bold')
                .attr('fill', '#495057')
                .text(`Pattern: ${this.data.pattern}`);
        }
    }

    setupInteractions() {
        const graphGroup = this.svg.select('.graph-group');

        // Make nodes draggable with smooth transitions
        const nodeGroups = graphGroup.selectAll('.node-group');

        nodeGroups
            .call(d3.drag()
                .on('start', (event, d) => {
                    d3.select(event.sourceEvent.target.parentNode)
                        .raise()
                        .style('cursor', 'grabbing');
                })
                .on('drag', (event, d) => {
                    d.x = event.x;
                    d.y = event.y;
                    this.updatePositions();
                })
                .on('end', (event, d) => {
                    d3.select(event.sourceEvent.target.parentNode)
                        .style('cursor', 'grab');
                })
            )
            .style('cursor', 'grab');

        // Add hover effects to show connections
        nodeGroups.on('mouseenter', (event, d) => {
            // Highlight the hovered node
            d3.select(event.currentTarget).select('.dfa-node')
                .transition()
                .duration(200)
                .attr('stroke-width', 5)
                .style('filter', 'drop-shadow(0 4px 8px rgba(0,0,0,0.2))');

            // Highlight connected edges
            graphGroup.selectAll('.dfa-link')
                .transition()
                .duration(200)
                .style('opacity', link => {
                    return (link.source === d.id || link.target === d.id) ? 1 : 0.3;
                })
                .attr('stroke-width', link => {
                    return (link.source === d.id || link.target === d.id) ? 3.5 : 2.5;
                });

            // Highlight connected edge labels
            graphGroup.selectAll('.edge-label')
                .transition()
                .duration(200)
                .style('opacity', label => {
                    return (label.link.source === d.id || label.link.target === d.id) ? 1 : 0.4;
                })
                .style('font-weight', label => {
                    return (label.link.source === d.id || label.link.target === d.id) ? '700' : '600';
                });
        })
        .on('mouseleave', (event, d) => {
            // Reset node appearance
            d3.select(event.currentTarget).select('.dfa-node')
                .transition()
                .duration(200)
                .attr('stroke-width', 3)
                .style('filter', 'drop-shadow(0 2px 4px rgba(0,0,0,0.1))');

            // Reset all edges
            graphGroup.selectAll('.dfa-link')
                .transition()
                .duration(200)
                .style('opacity', 0.85)
                .attr('stroke-width', 2.5);

            // Reset all labels
            graphGroup.selectAll('.edge-label')
                .transition()
                .duration(200)
                .style('opacity', 1)
                .style('font-weight', '600');
        });

        // Add smooth transitions to edges
        graphGroup.selectAll('.dfa-link')
            .on('mouseenter', function() {
                d3.select(this)
                    .transition()
                    .duration(200)
                    .attr('stroke-width', 4)
                    .style('opacity', 1);
            })
            .on('mouseleave', function() {
                d3.select(this)
                    .transition()
                    .duration(200)
                    .attr('stroke-width', 2.5)
                    .style('opacity', 0.85);
            });
    }

    positionLabelsWithCollisionAvoidance(edgeLabels) {
        // Get the current label data from the edgeLabels selection
        const labelData = [];
        edgeLabels.each(function(d, i) {
            labelData.push({
                link: d.link,
                label: d.label,
                x: d.x,
                y: d.y,
                width: 40, // Approximate width
                height: 16, // Approximate height
                index: i
            });
        });

        // Handle perpendicular intersections (collision avoidance)
        this.handlePerpendicularIntersections(labelData);

        // Apply final positions
        edgeLabels
            .attr('x', (d, i) => labelData[i].x)
            .attr('y', (d, i) => labelData[i].y);
    }

    handleBidirectionalLabels(labelData) {
        // Group labels by node pairs (bidirectional edges)
        const pairGroups = new Map();

        labelData.forEach(label => {
            const link = label.link;
            if (!link.isSelfLoop) {
                const pairKey = link.source < link.target ?
                    `${link.source}-${link.target}` :
                    `${link.target}-${link.source}`;
                if (!pairGroups.has(pairKey)) {
                    pairGroups.set(pairKey, []);
                }
                pairGroups.get(pairKey).push(label);
            }
        });

        // For bidirectional pairs, position labels on opposite sides
        pairGroups.forEach(group => {
            if (group.length === 2) {
                const [label1, label2] = group;
                const link1 = label1.link;
                const link2 = label2.link;

                // Determine which direction each label represents
                const isForward1 = link1.source === Math.min(link1.source, link1.target);
                const isForward2 = link2.source === Math.min(link2.source, link2.target);

                // Position labels on opposite sides of the edge
                const offset = 15; // Distance from center line

                // Get edge direction
                const source1 = this.nodes.find(n => n.id === link1.source);
                const target1 = this.nodes.find(n => n.id === link1.target);
                const dx = target1.x - source1.x;
                const dy = target1.y - source1.y;
                const length = Math.sqrt(dx * dx + dy * dy);

                if (length > 0) {
                    const perpX = -dy / length;
                    const perpY = dx / length;

                    // Position first label
                    label1.x += perpX * offset;
                    label1.y += perpY * offset;

                    // Position second label on opposite side
                    label2.x -= perpX * offset;
                    label2.y -= perpY * offset;
                }
            }
        });
    }

    handlePerpendicularIntersections(labelData) {
        const minDistance = 30; // Increased minimum distance between labels
        const maxIterations = 5; // Multiple passes for better separation

        // Multiple iterations for better collision resolution
        for (let iteration = 0; iteration < maxIterations; iteration++) {
            let collisionDetected = false;

            // Check all pairs and separate if too close (exclude self-loop labels)
            for (let i = 0; i < labelData.length; i++) {
                for (let j = i + 1; j < labelData.length; j++) {
                    const label1 = labelData[i];
                    const label2 = labelData[j];

                    // Skip collision avoidance for self-loop labels
                    if (label1.link.isSelfLoop || label2.link.isSelfLoop) continue;

                    // Calculate distance between label centers
                    const dx = label2.x - label1.x;
                    const dy = label2.y - label1.y;
                    const distance = Math.sqrt(dx * dx + dy * dy);

                    if (distance < minDistance && distance > 0.1) {
                        collisionDetected = true;

                        // Calculate separation needed
                        const separation = (minDistance - distance) / 2 + 2;
                        const unitX = dx / distance;
                        const unitY = dy / distance;

                        // Move labels apart along their connection line
                        label1.x -= unitX * separation;
                        label1.y -= unitY * separation;
                        label2.x += unitX * separation;
                        label2.y += unitY * separation;
                    }
                }
            }

            // If no collisions detected, we can stop early
            if (!collisionDetected) break;
        }

        // Additional pass: push labels away from node positions (exclude self-loops)
        this.nodes.forEach(node => {
            labelData.forEach(label => {
                if (label.link.isSelfLoop) return; // Skip self-loop labels

                const dx = label.x - node.x;
                const dy = label.y - node.y;
                const distance = Math.sqrt(dx * dx + dy * dy);
                const minNodeDistance = this.nodeRadius + 15;

                if (distance < minNodeDistance && distance > 0.1) {
                    const pushDistance = minNodeDistance - distance;
                    const unitX = dx / distance;
                    const unitY = dy / distance;

                    label.x += unitX * pushDistance;
                    label.y += unitY * pushDistance;
                }
            });
        });
    }

    updatePositions() {
        const graphGroup = this.svg.select('.graph-group');

        // Update node positions
        graphGroup.selectAll('.node-group')
            .attr('transform', d => `translate(${d.x}, ${d.y})`);

        // Node labels are inside the translated groups, so they move automatically
        // No need to update their absolute positions

        // Update edge paths
        graphGroup.selectAll('.dfa-link')
            .attr('d', d => this.getEdgePath(d));

        // Remove and recreate edge labels with updated positions
        graphGroup.selectAll('.edge-label').remove();
        graphGroup.selectAll('.edge-label-bg').remove();

        const labelData = this.prepareLabelData();

        // Add background rectangles
        const labelBgs = graphGroup.selectAll('.edge-label-bg')
            .data(labelData)
            .enter().append('rect')
            .attr('class', 'edge-label-bg')
            .attr('fill', 'rgba(255, 255, 255, 0.95)')
            .attr('stroke', 'rgba(0, 0, 0, 0.1)')
            .attr('stroke-width', 1)
            .attr('rx', 4)
            .attr('ry', 4);

        // Add labels
        const edgeLabels = graphGroup.selectAll('.edge-label')
            .data(labelData)
            .enter().append('text')
            .attr('class', 'edge-label')
            .attr('text-anchor', 'middle')
            .attr('font-size', '14px')
            .attr('font-weight', '600')
            .attr('fill', '#212529')
            .attr('pointer-events', 'none')
            .text(d => d.label)
            .attr('x', d => d.x)
            .attr('y', d => d.y);

        // Apply collision avoidance
        this.positionLabelsWithCollisionAvoidance(edgeLabels);

        // Update background rectangles
        edgeLabels.each(function(d, i) {
            const bbox = this.getBBox();
            const bg = d3.select(labelBgs.nodes()[i]);
            bg.attr('x', bbox.x - 4)
              .attr('y', bbox.y - 2)
              .attr('width', bbox.width + 8)
              .attr('height', bbox.height + 4);
        });
    }

    // Method to highlight a path through the DFA with improved visuals
    highlightPath(inputString) {
        if (!this.data || !inputString) {
            this.resetHighlight();
            return;
        }

        const graphGroup = this.svg.select('.graph-group');

        // Reset all highlights first
        this.resetHighlight();

        let currentState = 0;
        const pathStates = [0];
        const pathLinks = [];

        // Simulate DFA execution
        for (let i = 0; i < inputString.length; i++) {
            const char = inputString[i].toUpperCase();
            const transition = this.data.transitions.find(t =>
                t.from === currentState && t.label === char
            );

            if (transition) {
                currentState = transition.to;
                pathStates.push(currentState);
                pathLinks.push(transition);
            } else {
                break;
            }
        }

        // Highlight path states with smooth transitions
        graphGroup.selectAll('.dfa-node')
            .filter(d => pathStates.includes(d.id))
            .transition()
            .duration(400)
            .attr('stroke', '#ff6b35')
            .attr('stroke-width', 5)
            .style('filter', 'drop-shadow(0 4px 8px rgba(255,107,53,0.4))');

        // Highlight path links with smooth transitions
        graphGroup.selectAll('.dfa-link')
            .filter(d => {
                return pathLinks.some(t => {
                    if (d.isBidirectional) {
                        // Check both directions
                        return (t.from === d.source && t.to === d.target && d.direction === 'forward') ||
                               (t.from === d.reverseSource && t.to === d.reverseTarget && d.direction === 'reverse');
                    } else {
                        return t === d.originalTransition;
                    }
                });
            })
            .transition()
            .duration(400)
            .attr('stroke', '#ff6b35')
            .attr('stroke-width', 4)
            .style('opacity', 1)
            .attr('marker-end', d => {
                if (d.isBidirectional) {
                    return d.direction === 'forward' ? 'url(#arrowhead-end-highlight)' : null;
                }
                return 'url(#arrowhead-end-highlight)';
            });

        // Show accept/reject feedback
        const isAccepted = this.data.accept_states.includes(currentState) &&
                          pathStates.length === inputString.length + 1;

        this.showPathFeedback(isAccepted, currentState, inputString);
    }

    resetHighlight() {
        const graphGroup = this.svg.select('.graph-group');

        // Reset nodes
        graphGroup.selectAll('.dfa-node')
            .transition()
            .duration(300)
            .attr('stroke', '#ffffff')
            .attr('stroke-width', 3)
            .style('filter', 'drop-shadow(0 2px 4px rgba(0,0,0,0.1))');

        // Reset links
        graphGroup.selectAll('.dfa-link')
            .transition()
            .duration(300)
            .attr('stroke', '#495057')
            .attr('stroke-width', 2.5)
            .style('opacity', 0.85)
            .attr('marker-end', d => {
                if (d.isBidirectional) {
                    return d.direction === 'forward' ? 'url(#arrowhead-end)' : null;
                }
                return 'url(#arrowhead-end)';
            });

        // Remove feedback message
        d3.select('#path-feedback').remove();
    }

    showPathFeedback(isAccepted, finalState, inputString) {
        // Remove any existing feedback
        d3.select('#path-feedback').remove();

        // Create feedback element
        const feedbackColor = isAccepted ? '#28a745' : '#dc3545';
        const feedbackIcon = isAccepted ? '✓' : '✗';
        const feedbackText = isAccepted ?
            `Accepted! Ended in accept state ${finalState}` :
            `Rejected. Ended in state ${finalState}`;

        const feedback = d3.select('#dfa-viz-container')
            .append('div')
            .attr('id', 'path-feedback')
            .style('position', 'absolute')
            .style('top', '15px')
            .style('left', '50%')
            .style('transform', 'translateX(-50%)')
            .style('background-color', feedbackColor)
            .style('color', 'white')
            .style('padding', '12px 24px')
            .style('border-radius', '8px')
            .style('font-weight', '600')
            .style('font-size', '15px')
            .style('box-shadow', '0 4px 12px rgba(0,0,0,0.2)')
            .style('z-index', '1000')
            .style('opacity', '0')
            .html(`<span style="font-size: 18px; margin-right: 8px;">${feedbackIcon}</span>${feedbackText}`);

        // Animate in
        feedback.transition()
            .duration(400)
            .style('opacity', '1');
    }

    // Resize the visualization
    resize(width, height) {
        this.width = width;
        this.height = height;
        this.svg.attr('width', width).attr('height', height);

        if (this.simulation) {
            this.simulation.force('center', d3.forceCenter(width / 2, height / 2));
            this.simulation.restart();
        }
    }
}

// Global function to render state diagram (called from HTML)
function renderStateDiagram(diagramData) {
    const container = document.getElementById('stateDiagram');

    // Clear container
    container.innerHTML = '';

    // Add controls
    const controls = document.createElement('div');
    controls.className = 'dfa-controls mb-3';
    controls.innerHTML = `
        <div class="d-flex justify-content-between align-items-center mb-2">
            <small class="text-muted">
                <i class="fas fa-info-circle"></i>
                Drag nodes • Scroll to zoom • Hover for details
            </small>
            <div class="dfa-zoom-controls">
                <button onclick="zoomIn()" title="Zoom In" class="btn btn-sm btn-outline-secondary">+</button>
                <button onclick="zoomOut()" title="Zoom Out" class="btn btn-sm btn-outline-secondary">-</button>
                <button onclick="resetZoom()" title="Reset Zoom" class="btn btn-sm btn-outline-secondary">↻</button>
            </div>
        </div>
        <div class="row">
            <div class="col-md-8">
                <label for="inputString" class="form-label">Test Input String:</label>
                <div class="input-group input-group-sm">
                    <input type="text" class="form-control" id="inputString"
                           placeholder="Enter DNA sequence to highlight path" maxlength="20">
                    <button class="btn btn-outline-primary" onclick="highlightPath()">
                        <i class="fas fa-play"></i> Highlight Path
                    </button>
                </div>
            </div>
            <div class="col-md-4">
                <div class="dfa-legend">
                    <div class="dfa-legend-item">
                        <div class="dfa-legend-color" style="background-color: #007bff;"></div>
                        <span class="dfa-legend-text">Regular State</span>
                    </div>
                    <div class="dfa-legend-item">
                        <div class="dfa-legend-color" style="background-color: #28a745; border: 2px dashed #28a745;"></div>
                        <span class="dfa-legend-text">Accept State</span>
                    </div>
                    <div class="dfa-legend-item">
                        <div class="dfa-legend-color" style="background-color: #ff6b35;"></div>
                        <span class="dfa-legend-text">Active Path</span>
                    </div>
                </div>
            </div>
        </div>
    `;
    container.appendChild(controls);

    // Create visualization container
    const vizContainer = document.createElement('div');
    vizContainer.id = 'dfa-viz-container';
    vizContainer.className = 'dfa-visualization-container';
    vizContainer.style.width = '100%';
    vizContainer.style.height = '600px';
    vizContainer.style.position = 'relative';
    vizContainer.style.border = '1px solid #e9ecef';
    vizContainer.style.borderRadius = '8px';
    vizContainer.style.backgroundColor = '#f8f9fa';
    container.appendChild(vizContainer);

    // Create visualizer with larger size
    const visualizer = new DFAVisualizer('dfa-viz-container', 800, 600);
    visualizer.render(diagramData);

    // Store visualizer globally for interaction
    window.currentDFAVisualizer = visualizer;
}

// Global function to highlight path
function highlightPath() {
    const inputString = document.getElementById('inputString').value.toUpperCase();
    if (window.currentDFAVisualizer) {
        window.currentDFAVisualizer.highlightPath(inputString);
    }
}

// Zoom control functions - fixed implementation
function zoomIn() {
    if (window.currentDFAVisualizer && window.currentDFAVisualizer.svg && window.currentDFAVisualizer.zoom) {
        window.currentDFAVisualizer.svg
            .transition()
            .duration(300)
            .call(window.currentDFAVisualizer.zoom.scaleBy, 1.3);
    }
}

function zoomOut() {
    if (window.currentDFAVisualizer && window.currentDFAVisualizer.svg && window.currentDFAVisualizer.zoom) {
        window.currentDFAVisualizer.svg
            .transition()
            .duration(300)
            .call(window.currentDFAVisualizer.zoom.scaleBy, 0.77);
    }
}

function resetZoom() {
    if (window.currentDFAVisualizer && window.currentDFAVisualizer.svg && window.currentDFAVisualizer.zoom) {
        window.currentDFAVisualizer.svg
            .transition()
            .duration(500)
            .call(window.currentDFAVisualizer.zoom.transform, d3.zoomIdentity);
    }
}
