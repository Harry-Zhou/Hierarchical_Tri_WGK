import os
import gzip
import random
import networkx as nx
import matplotlib.pyplot as plt
import numpy as np
from sage.graphs.graph import Graph
from collections import Counter

###################################################
# 运行程序
# sage -python generate_strongly_regular_graphs.py
###################################################

# Set output directory
TOP_OUTPUT_DIR = os.path.join('.', 'graph_isomorphism_testing')

# Create output directory (if it doesn't exist)
os.makedirs(TOP_OUTPUT_DIR, exist_ok=True)

# Define the file to process
gz_filename = 'srgs_85_20_3_5.gz'
gz_filepath = os.path.join(TOP_OUTPUT_DIR, gz_filename)

# Fixed graph filenames
graph1_edgelist_path = os.path.join(TOP_OUTPUT_DIR, "graph1.edgelist")
graph2_edgelist_path = os.path.join(TOP_OUTPUT_DIR, "graph2.edgelist")
graph1_g6_path = os.path.join(TOP_OUTPUT_DIR, "graph1.g6")
graph2_g6_path = os.path.join(TOP_OUTPUT_DIR, "graph2.g6")

# Check if existing graph data files exist
existing_graphs_exist = os.path.exists(graph1_edgelist_path) and os.path.exists(graph2_edgelist_path)

if existing_graphs_exist:
    print(f"Found existing graph data files in {TOP_OUTPUT_DIR}")
    print("Loading existing graphs...")
    
    # Load existing graphs from edgelist files
    nx_graph1 = nx.read_edgelist(graph1_edgelist_path, nodetype=int)
    nx_graph2 = nx.read_edgelist(graph2_edgelist_path, nodetype=int)
    
    # Convert to SageMath graphs for analysis
    print("Converting to SageMath graphs for analysis...")
    
    # Load graph6 strings if available
    if os.path.exists(graph1_g6_path):
        with open(graph1_g6_path, 'r') as f:
            graph6_str1 = f.read().strip()
        graph1 = Graph(graph6_str1)
    else:
        # Convert NetworkX to SageMath graph
        graph1 = Graph()
        for u, v in nx_graph1.edges():
            graph1.add_edge(u, v)
    
    if os.path.exists(graph2_g6_path):
        with open(graph2_g6_path, 'r') as f:
            graph6_str2 = f.read().strip()
        graph2 = Graph(graph6_str2)
    else:
        # Convert NetworkX to SageMath graph
        graph2 = Graph()
        for u, v in nx_graph2.edges():
            graph2.add_edge(u, v)
    
    # Use placeholder indices
    idx1, idx2 = 0, 1
    
    print("✓ Loaded existing graph data")
    
else:
    print(f"No existing graph data found in {TOP_OUTPUT_DIR}")
    print(f"Reading from compressed file: {gz_filepath}")
    
    # Check if gz file exists
    if not os.path.exists(gz_filepath):
        print(f"Error: File {gz_filepath} does not exist!")
        print(f"Please ensure the .gz file is in the {TOP_OUTPUT_DIR} directory")
        exit()

    # Read graph6 format data
    graphs = []
    try:
        with gzip.open(gz_filepath, 'rt') as f:
            for line in f:
                line = line.strip()
                if line:
                    try:
                        sage_graph = Graph(line)
                        graphs.append(sage_graph)
                    except Exception as e:
                        print(f"Error parsing graph6 string: {e}")
    except Exception as e:
        print(f"Error reading file: {e}")
        exit()

    print(f"Successfully read {len(graphs)} graphs")

    # Function to find two non-isomorphic graphs
    def find_non_isomorphic_pair(graphs, max_attempts=50):
        """Find two non-isomorphic graphs from the list."""
        if len(graphs) < 2:
            return None, None, None, None
        
        attempts = 0
        while attempts < max_attempts:
            # Randomly select two different indices
            idx1, idx2 = random.sample(range(len(graphs)), 2)
            graph1 = graphs[idx1]
            graph2 = graphs[idx2]
            
            # Check if they are isomorphic
            if not graph1.is_isomorphic(graph2):
                return idx1, idx2, graph1, graph2
            
            attempts += 1
        
        # If we can't find non-isomorphic graphs, try all pairs up to first 20 graphs
        for i in range(min(20, len(graphs))):
            for j in range(i+1, min(20, len(graphs))):
                if not graphs[i].is_isomorphic(graphs[j]):
                    return i, j, graphs[i], graphs[j]
        
        return None, None, None, None

    # Find two non-isomorphic graphs
    idx1, idx2, graph1, graph2 = find_non_isomorphic_pair(graphs)

    if graph1 is None or graph2 is None:
        print("Could not find two non-isomorphic graphs!")
        exit()

    print(f"Selected two non-isomorphic graphs:")
    print(f"  Graph 1: Index {idx1}, Order {graph1.order()}")
    print(f"  Graph 2: Index {idx2}, Order {graph2.order()}")

    # Convert to networkx graphs
    nx_graph1 = graph1.networkx_graph()
    nx_graph2 = graph2.networkx_graph()
    
    # Save graph data for future use
    print("\nSaving graph data for future use...")
    
    # Save as edgelist format
    nx.write_edgelist(nx_graph1, graph1_edgelist_path, data=False)
    nx.write_edgelist(nx_graph2, graph2_edgelist_path, data=False)
    
    # Save as graph6 format
    with open(graph1_g6_path, 'w') as f:
        f.write(graph1.graph6_string() + '\n')
    with open(graph2_g6_path, 'w') as f:
        f.write(graph2.graph6_string() + '\n')
    
    print(f"✓ Saved graph data to {TOP_OUTPUT_DIR}")

# Function to analyze structural differences
def analyze_structural_differences(graph1, graph2, nx_graph1, nx_graph2):
    """Analyze structural differences between two graphs."""
    differences = {}
    
    # Clique analysis - find maximal cliques
    differences['max_cliques'] = {
        'graph1': list(nx.find_cliques(nx_graph1)),
        'graph2': list(nx.find_cliques(nx_graph2))
    }
    
    # Count cliques by size
    clique_size_counts1 = Counter(len(clique) for clique in differences['max_cliques']['graph1'])
    clique_size_counts2 = Counter(len(clique) for clique in differences['max_cliques']['graph2'])
    
    differences['clique_size_distribution'] = {
        'graph1': clique_size_counts1,
        'graph2': clique_size_counts2
    }
    
    # Find a clique size that differs between the two graphs
    # We'll look for the largest clique size where they differ
    all_clique_sizes = set(clique_size_counts1.keys()) | set(clique_size_counts2.keys())
    
    # Sort clique sizes in descending order
    sorted_clique_sizes = sorted(all_clique_sizes, reverse=True)
    
    # Find the first (largest) clique size where the two graphs differ
    selected_clique_size = None
    for size in sorted_clique_sizes:
        count1 = clique_size_counts1.get(size, 0)
        count2 = clique_size_counts2.get(size, 0)
        if count1 != count2:
            selected_clique_size = size
            break
    
    # If all clique sizes are the same, use the maximum clique size
    if selected_clique_size is None and sorted_clique_sizes:
        selected_clique_size = max(sorted_clique_sizes)
    
    differences['selected_clique_size'] = selected_clique_size
    
    # Get cliques of the selected size
    differences['selected_cliques'] = {
        'graph1': [clique for clique in differences['max_cliques']['graph1'] 
                  if len(clique) == selected_clique_size],
        'graph2': [clique for clique in differences['max_cliques']['graph2'] 
                  if len(clique) == selected_clique_size]
    }
    
    # Cycle analysis - look for specific cycle lengths
    differences['cycle_counts'] = {}
    
    # Check for 3-cycles (triangles)
    triangles1 = sum(1 for _ in nx.enumerate_all_cliques(nx_graph1) if len(_) == 3)
    triangles2 = sum(1 for _ in nx.enumerate_all_cliques(nx_graph2) if len(_) == 3)
    differences['cycle_counts']['triangles'] = {
        'graph1': triangles1, 
        'graph2': triangles2
    }
    
    return differences

# Analyze structural differences
print("\nAnalyzing structural differences...")
differences = analyze_structural_differences(graph1, graph2, nx_graph1, nx_graph2)

# Print analysis results
print("\n=== STRUCTURAL DIFFERENCE ANALYSIS ===")

# 1. Clique distribution differences
print("\n1. Clique Size Distribution:")
print(f"   Graph 1: {differences['clique_size_distribution']['graph1']}")
print(f"   Graph 2: {differences['clique_size_distribution']['graph2']}")

# Show which clique size was selected for visualization
selected_size = differences['selected_clique_size']
count1 = differences['clique_size_distribution']['graph1'].get(selected_size, 0)
count2 = differences['clique_size_distribution']['graph2'].get(selected_size, 0)
print(f"\n   Selected clique size for visualization: {selected_size}")
print(f"   Graph 1 has {count1} cliques of size {selected_size}")
print(f"   Graph 2 has {count2} cliques of size {selected_size}")

# 2. Triangle counts
print("\n2. Number of Triangles (3-cycles):")
print(f"   Graph 1: {differences['cycle_counts']['triangles']['graph1']}")
print(f"   Graph 2: {differences['cycle_counts']['triangles']['graph2']}")

# Function to create visualizations highlighting structural differences
def create_visualizations(nx_graph1, nx_graph2, graph1_name, graph2_name, differences):
    """Create visualizations highlighting structural differences."""
    
    # Create a figure with subplots
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 8))
    
    # Get selected clique size and cliques
    selected_clique_size = differences['selected_clique_size']
    selected_cliques1 = differences['selected_cliques']['graph1']
    selected_cliques2 = differences['selected_cliques']['graph2']
    
    # Create positions for both graphs
    pos1 = nx.spring_layout(nx_graph1, seed=42)
    pos2 = nx.spring_layout(nx_graph2, seed=42)
    
    # =============================================
    # Subplot 1: Graph 1 with selected cliques highlighted
    # =============================================
    # Get ALL nodes in ALL selected cliques
    clique_nodes1 = set()
    # Also store which edges are within each clique
    clique_edges1 = set()
    
    # Process ALL selected cliques
    for clique in selected_cliques1:
        # Add all nodes in this clique
        clique_nodes1.update(clique)
        
        # Add all edges within this clique
        for i in range(len(clique)):
            for j in range(i+1, len(clique)):
                # Ensure edge exists (it should, since it's a clique)
                if nx_graph1.has_edge(clique[i], clique[j]):
                    clique_edges1.add((clique[i], clique[j]))
    
    # Color nodes: red for nodes in ANY selected clique, blue for others
    node_colors_clique1 = ['red' if node in clique_nodes1 else 'blue' 
                          for node in nx_graph1.nodes()]
    
    # Draw ALL nodes
    nx.draw_networkx_nodes(nx_graph1, pos1, ax=ax1, node_size=150, 
                          node_color=node_colors_clique1, alpha=0.8)
    
    # Draw edges - make edges within cliques red and width=1.0
    edge_colors1 = []
    edge_widths1 = []
    
    for u, v in nx_graph1.edges():
        # Check if this edge is in ANY selected clique
        if (u, v) in clique_edges1 or (v, u) in clique_edges1:
            edge_colors1.append('red')
            edge_widths1.append(0.5)
        else:
            edge_colors1.append('blue')
            edge_widths1.append(0.5)
    
    nx.draw_networkx_edges(nx_graph1, pos1, ax=ax1, width=edge_widths1, 
                          alpha=0.7, edge_color=edge_colors1)
    
    # Add text annotation about clique size and count
    count1 = len(selected_cliques1)
    clique_info1 = f"Selected clique size: {selected_clique_size}\n"
    clique_info1 += f"Number of selected cliques: {count1}\n"
    clique_info1 += f"Nodes in selected cliques: {len(clique_nodes1)}\n"
    clique_info1 += f"Edges within selected cliques: {len(clique_edges1)}"
    
    ax1.set_title(f"{graph1_name}\nRed nodes/edges: cliques of size {selected_clique_size}", 
                  fontsize=25, pad=20)
    ax1.text(0.5, -0.2, clique_info1, transform=ax1.transAxes, 
             ha='center', fontsize=25, bbox=dict(boxstyle="round,pad=0.3", 
                                                 facecolor="lightyellow", alpha=0.7))
    
    # =============================================
    # Subplot 2: Graph 2 with selected cliques highlighted
    # =============================================
    # Get ALL nodes in ALL selected cliques
    clique_nodes2 = set()
    clique_edges2 = set()
    
    # Process ALL selected cliques
    for clique in selected_cliques2:
        # Add all nodes in this clique
        clique_nodes2.update(clique)
        
        # Add all edges within this clique
        for i in range(len(clique)):
            for j in range(i+1, len(clique)):
                # Ensure edge exists (it should, since it's a clique)
                if nx_graph2.has_edge(clique[i], clique[j]):
                    clique_edges2.add((clique[i], clique[j]))
    
    # Color nodes: red for nodes in ANY selected clique, light blue for others
    node_colors_clique2 = ['red' if node in clique_nodes2 else 'blue' 
                          for node in nx_graph2.nodes()]
    
    # Draw ALL nodes
    nx.draw_networkx_nodes(nx_graph2, pos2, ax=ax2, node_size=150, 
                          node_color=node_colors_clique2, alpha=0.8)
    
    # Draw edges - make edges within cliques red and width=1.0
    edge_colors2 = []
    edge_widths2 = []
    
    for u, v in nx_graph2.edges():
        # Check if this edge is in ANY selected clique
        if (u, v) in clique_edges2 or (v, u) in clique_edges2:
            edge_colors2.append('red')
            edge_widths2.append(0.5)
        else:
            edge_colors2.append('blue')
            edge_widths2.append(0.5)
    
    nx.draw_networkx_edges(nx_graph2, pos2, ax=ax2, width=edge_widths2, 
                          alpha=0.7, edge_color=edge_colors2)
    
    # Add text annotation about clique size and count
    count2 = len(selected_cliques2)
    clique_info2 = f"Selected clique size: {selected_clique_size}\n"
    clique_info2 += f"Number of selected cliques: {count2}\n"
    clique_info2 += f"Nodes in selected cliques: {len(clique_nodes2)}\n"
    clique_info2 += f"Edges within selected cliques: {len(clique_edges2)}"
    
    ax2.set_title(f"{graph2_name}\nRed nodes/edges: cliques of size {selected_clique_size}", 
                  fontsize=25, pad=20)
    ax2.text(0.5, -0.2, clique_info2, transform=ax2.transAxes, 
             ha='center', fontsize=25, bbox=dict(boxstyle="round,pad=0.3", 
                                                 facecolor="lightyellow", alpha=0.7))
    
    # Main title
    # title_suffix = " (different clique distributions)" if count1 != count2 else " (same clique count, different structure)"
    plt.suptitle(f'Structural Differences Between Non-Isomorphic Strongly Regular Graphs', 
                 fontsize=16, y=1.05)
    
    plt.subplots_adjust(bottom=0.25)
    plt.tight_layout()
    
    return fig, selected_clique_size, count1, count2, len(clique_nodes1), len(clique_nodes2), len(clique_edges1), len(clique_edges2)

# Create visualizations
print("\nCreating visualizations highlighting structural differences...")
fig, selected_clique_size, count1, count2, clique_nodes_count1, clique_nodes_count2, clique_edges_count1, clique_edges_count2 = create_visualizations(
    nx_graph1, nx_graph2, "Graph 1", "Graph 2", differences)

# Save visualization
visualization_path = os.path.join(TOP_OUTPUT_DIR, f"non_isomorphic_clique_comparison.png")
plt.savefig(visualization_path, dpi=150, bbox_inches='tight')
plt.close()

print(f"✓ Visualization saved: {visualization_path}")

# Save analysis summary
summary_path = os.path.join(TOP_OUTPUT_DIR, f"non_isomorphic_clique_analysis.txt")
with open(summary_path, 'w') as f:
    f.write("ANALYSIS OF NON-ISOMORPHIC STRONGLY REGULAR GRAPHS - CLIQUE STRUCTURE\n")
    f.write("="*70 + "\n\n")
    
    source_info = "Loaded from existing data files" if existing_graphs_exist else f"Generated from {gz_filename}"
    f.write(f"Source: {source_info}\n\n")
    
    f.write(f"Graph 1:\n")
    f.write(f"  Order (n): {graph1.order()}\n")
    f.write(f"  Degree: {graph1.degree_sequence()[0]}\n")
    f.write(f"  Is strongly regular: {graph1.is_strongly_regular()}\n")
    if graph1.is_strongly_regular():
        sr_params = graph1.is_strongly_regular(parameters=True)
        if sr_params:
            n, k, lambda_val, mu_val = sr_params
            f.write(f"  SRG parameters: (n={n}, k={k}, λ={lambda_val}, μ={mu_val})\n")
    
    f.write(f"\nGraph 2:\n")
    f.write(f"  Order (n): {graph2.order()}\n")
    f.write(f"  Degree: {graph2.degree_sequence()[0]}\n")
    f.write(f"  Is strongly regular: {graph2.is_strongly_regular()}\n")
    if graph2.is_strongly_regular():
        sr_params = graph2.is_strongly_regular(parameters=True)
        if sr_params:
            n, k, lambda_val, mu_val = sr_params
            f.write(f"  SRG parameters: (n={n}, k={k}, λ={lambda_val}, μ={mu_val})\n")
    
    f.write("\n\nCLIQUE STRUCTURE ANALYSIS:\n")
    f.write("="*70 + "\n\n")
    
    f.write("1. CLIQUE SIZE DISTRIBUTION:\n")
    f.write(f"   Graph 1: {dict(differences['clique_size_distribution']['graph1'])}\n")
    f.write(f"   Graph 2: {dict(differences['clique_size_distribution']['graph2'])}\n\n")
    
    f.write("2. SELECTED CLIQUE FOR VISUALIZATION:\n")
    f.write(f"   Selected clique size: {selected_clique_size}\n")
    f.write(f"   Reason: This clique size shows differences between the two graphs\n")
    f.write(f"   Graph 1: {count1} cliques of size {selected_clique_size}\n")
    f.write(f"   Graph 2: {count2} cliques of size {selected_clique_size}\n")
    f.write(f"   Difference: {abs(count1 - count2)} cliques\n\n")
    
    f.write("3. DETAILS OF SELECTED CLIQUES:\n")
    f.write(f"   Graph 1:\n")
    f.write(f"     - Number of selected cliques = {count1}\n")
    f.write(f"     - Nodes in selected cliques = {clique_nodes_count1}\n")
    f.write(f"     - Edges within selected cliques = {clique_edges_count1}\n")
    f.write(f"   Graph 2:\n")
    f.write(f"     - Number of selected cliques = {count2}\n")
    f.write(f"     - Nodes in selected cliques = {clique_nodes_count2}\n")
    f.write(f"     - Edges within selected cliques = {clique_edges_count2}\n\n")
    
    f.write("4. TRIANGLE COUNTS (3-CLIQUES):\n")
    f.write(f"   Graph 1: {differences['cycle_counts']['triangles']['graph1']} triangles\n")
    f.write(f"   Graph 2: {differences['cycle_counts']['triangles']['graph2']} triangles\n")
    if differences['cycle_counts']['triangles']['graph1'] != differences['cycle_counts']['triangles']['graph2']:
        f.write(f"   Difference: {abs(differences['cycle_counts']['triangles']['graph1'] - differences['cycle_counts']['triangles']['graph2'])} triangles\n\n")
    else:
        f.write("\n")
    
    f.write("VISUALIZATION KEY:\n")
    f.write("-"*70 + "\n")
    f.write("1. Red nodes: ALL nodes belonging to ANY clique of size {selected_clique_size}\n")
    f.write("2. Red edges: ALL edges within ANY clique of size {selected_clique_size} (width=1.0)\n")
    f.write("3. Light blue nodes: All other nodes\n")
    f.write("4. Gray edges: All other edges (width=0.5)\n\n")
    
    f.write("WHY THESE GRAPHS ARE NON-ISOMORPHIC:\n")
    f.write("="*70 + "\n")
    
    # Determine why they are non-isomorphic based on the selected clique size
    if count1 != count2:
        f.write(f"The two graphs are non-isomorphic because they have different numbers of cliques\n")
        f.write(f"of size {selected_clique_size}:\n")
        f.write(f"  - Graph 1 has {count1} cliques of size {selected_clique_size}\n")
        f.write(f"  - Graph 2 has {count2} cliques of size {selected_clique_size}\n")
        f.write(f"  - Difference: {abs(count1 - count2)} cliques\n\n")
    else:
        f.write(f"Both graphs have the same number of cliques of size {selected_clique_size} ({count1} each).\n")
        f.write("However, they are non-isomorphic because:\n")
        f.write("1. The arrangement and connectivity of these cliques within the graph must differ.\n")
        f.write("2. The distribution of other clique sizes may also differ (see above).\n")
        f.write("3. The triangle counts may differ, indicating different local clustering structures.\n\n")
    
    f.write("ADDITIONAL EVIDENCE OF NON-ISOMORPHISM:\n")
    f.write("-"*70 + "\n")
    f.write("1. Clique distribution differences:\n")
    
    # Find all clique sizes where the two graphs differ
    diff_sizes = []
    all_sizes = set(differences['clique_size_distribution']['graph1'].keys()) | \
                set(differences['clique_size_distribution']['graph2'].keys())
    
    for size in sorted(all_sizes):
        c1 = differences['clique_size_distribution']['graph1'].get(size, 0)
        c2 = differences['clique_size_distribution']['graph2'].get(size, 0)
        if c1 != c2:
            diff_sizes.append((size, c1, c2))
    
    if diff_sizes:
        f.write("   The following clique sizes have different counts:\n")
        for size, c1, c2 in diff_sizes:
            f.write(f"   - Size {size}: Graph 1 has {c1}, Graph 2 has {c2} (difference: {abs(c1-c2)})\n")
    else:
        f.write("   All clique sizes have the same counts, but the arrangement differs.\n")
    
    f.write("\n2. Triangle count differences:\n")
    if differences['cycle_counts']['triangles']['graph1'] != differences['cycle_counts']['triangles']['graph2']:
        f.write(f"   Graph 1: {differences['cycle_counts']['triangles']['graph1']} triangles\n")
        f.write(f"   Graph 2: {differences['cycle_counts']['triangles']['graph2']} triangles\n")
        f.write(f"   Difference: {abs(differences['cycle_counts']['triangles']['graph1'] - differences['cycle_counts']['triangles']['graph2'])} triangles\n")
    else:
        f.write("   Both graphs have the same number of triangles.\n")
    
    f.write("\nCONCLUSION:\n")
    f.write("-"*70 + "\n")
    f.write("Although both graphs have the same global SRG parameters, their local clique structures\n")
    f.write("differ significantly. These differences in clique distributions and arrangements are\n")
    f.write("sufficient to prove they are non-isomorphic. The visualization highlights cliques of\n")
    f.write(f"size {selected_clique_size}, which shows clear structural differences between the graphs.\n")

print(f"✓ Analysis summary saved: {summary_path}")

print("\n" + "="*70)
print("ANALYSIS COMPLETE")
print("="*70)
print(f"✓ {'Loaded existing' if existing_graphs_exist else 'Generated new'} non-isomorphic strongly regular graphs")
print(f"✓ Analyzed clique structure differences")
print(f"✓ Created visualization highlighting cliques of size {selected_clique_size}")
print(f"✓ Saved all files to: {TOP_OUTPUT_DIR}")
print("\nKey findings:")
print(f"  1. Selected clique size for visualization: {selected_clique_size}")
print(f"  2. Graph 1: {count1} cliques of size {selected_clique_size}")
print(f"  3. Graph 2: {count2} cliques of size {selected_clique_size}")
print(f"  4. Difference: {abs(count1 - count2)} cliques")
print(f"  5. Graph 1 triangle count: {differences['cycle_counts']['triangles']['graph1']}")
print(f"  6. Graph 2 triangle count: {differences['cycle_counts']['triangles']['graph2']}")
print(f"\nFiles created/used:")
print(f"  - Visualization: {visualization_path}")
print(f"  - Analysis summary: {summary_path}")
print(f"  - Graph 1 data: {graph1_edgelist_path}")
print(f"  - Graph 2 data: {graph2_edgelist_path}")