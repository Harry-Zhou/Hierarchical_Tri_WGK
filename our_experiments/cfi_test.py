import networkx as nx
import matplotlib
matplotlib.use('Agg')  # 非交互后端
import matplotlib.pyplot as plt
import numpy as np

# ========================
# 1. Base graph K4 definition
# ========================
vertices = [1, 2, 3, 4]
edges = [(1,2), (2,3), (1,3), (1,4), (2,4), (3,4)]

# Triangle (1,2),(2,3),(1,3) as twisted edge set
twisted_edges = {(1,2), (2,3), (1,3)}

# ========================
# 2. Build unlabeled CFI graph
# ========================
def build_cfi(twist_set):
    G = nx.Graph()
    # Add vertex clones: v0, v1 (all nodes have no color label)
    for v in vertices:
        G.add_node(f"v{v}_0")
        G.add_node(f"v{v}_1")
    # Add edge nodes: each edge produces e_a, e_b
    for (u,v) in edges:
        G.add_node(f"e{u}{v}_a")
        G.add_node(f"e{u}{v}_b")
        # Check if this edge is twisted
        twist = (u,v) in twist_set or (v,u) in twist_set
        if not twist:   # parallel connection
            G.add_edge(f"e{u}{v}_a", f"v{u}_0")
            G.add_edge(f"e{u}{v}_a", f"v{v}_0")
            G.add_edge(f"e{u}{v}_b", f"v{u}_1")
            G.add_edge(f"e{u}{v}_b", f"v{v}_1")
        else:           # twisted connection
            G.add_edge(f"e{u}{v}_a", f"v{u}_0")
            G.add_edge(f"e{u}{v}_a", f"v{v}_1")
            G.add_edge(f"e{u}{v}_b", f"v{u}_1")
            G.add_edge(f"e{u}{v}_b", f"v{v}_0")
    return G

G1 = build_cfi(set())           # CFI(K4): all parallel
G2 = build_cfi(twisted_edges)   # CFI'(K4): triangle twisted

# ========================
# 3. Analyze minimum cycle basis
# ========================
from networkx.algorithms.cycles import minimum_cycle_basis

cycles_G1 = minimum_cycle_basis(G1)
cycles_G2 = minimum_cycle_basis(G2)

# Count cycles by length
def count_cycle_lengths(cycles):
    length_counts = {}
    for cycle in cycles:
        length = len(cycle)
        length_counts[length] = length_counts.get(length, 0) + 1
    return length_counts

length_dist_G1 = count_cycle_lengths(cycles_G1)
length_dist_G2 = count_cycle_lengths(cycles_G2)

print("\n" + "=" * 60)
print("MINIMUM CYCLE BASIS ANALYSIS")
print("=" * 60)
print(f"\nCFI(K4) - Number of cycles in minimum basis: {len(cycles_G1)}")
print(f"  Cycle length distribution: {length_dist_G1}")

print(f"\nCFI'(K4) - Number of cycles in minimum basis: {len(cycles_G2)}")
print(f"  Cycle length distribution: {length_dist_G2}")

print("\n" + "=" * 60)
print("FULL MINIMUM CYCLE BASIS")
print("=" * 60)

print("\nCFI(K4) - 6 cycles (all length 6):")
for i, cyc in enumerate(cycles_G1):
    print(f"  Cycle {i+1} (length {len(cyc)}): {list(cyc)}")

print("\nCFI'(K4) - 5 cycles (all length 8):")
for i, cyc in enumerate(cycles_G2):
    print(f"  Cycle {i+1} (length {len(cyc)}): {list(cyc)}")

# ========================
# 3. Find structural difference: a cycle with different parity
# ========================
# Cycle through triangle vertices: v1_0 -> e12_a -> v2_0 -> e23_a -> v3_0 -> e13_a -> v1_0
# In CFI(K4): 0 twisted edges (even parity)
# In CFI'(K4): 3 twisted edges (odd parity)
cycle_path_G1 = [('v1_0', 'e12_a'), ('e12_a', 'v2_0'), ('v2_0', 'e23_a'),
                 ('e23_a', 'v3_0'), ('v3_0', 'e13_a'), ('e13_a', 'v1_0')]

cycle_path_G2 = [('v1_0', 'e12_a'), ('e12_a', 'v2_0'), ('v2_0', 'e23_a'),
                 ('e23_a', 'v3_0'), ('v3_0', 'e13_a'), ('e13_a', 'v1_0')]

# Find twisted edges: edges connecting to edge-nodes from twisted base edges
# These are the edges that show "crossed" connection pattern
twisted_connection_edges_G1 = []  # No twisted edges in G1
twisted_connection_edges_G2 = []  # Edges from twisted edge-nodes in G2

for (u, v) in G2.edges():
    edge_node = u if u.startswith('e') else v
    if edge_node.startswith('e'):
        parts = edge_node[1:].split("_")[0]
        a = int(parts[0])
        b = int(parts[1]) if len(parts) > 1 else int(parts[0])
        if (a, b) in twisted_edges or (b, a) in twisted_edges:
            twisted_connection_edges_G2.append((u, v))

def count_twisted_edges(cycle_edges, twist_set):
    count = 0
    for u, v in cycle_edges:
        edge_node = u if u.startswith('e') else v
        parts = edge_node[1:].split("_")[0]
        a = int(parts[0])
        b = int(parts[1]) if len(parts)>1 else int(parts[0])
        if (a,b) in twist_set or (b,a) in twist_set:
            count += 1
    return count

twisted_count_G1 = count_twisted_edges(cycle_path_G1, set()) // 2
twisted_count_G2 = count_twisted_edges(cycle_path_G2, twisted_edges) // 2

print("Isomorphic?", nx.is_isomorphic(G1, G2))
print("=" * 60)
print("What are TWISTED EDGES in CFI construction?")
print("=" * 60)
print("For each original edge (u,v), CFI creates two edge-nodes (e_a, e_b)")
print("and connects them to vertex clones (v_u_0, v_u_1, v_v_0, v_v_1).")
print()
print("Parallel edge (non-twisted): e_a connects to v_u_0 AND v_v_0 (same color)")
print("                              e_b connects to v_u_1 AND v_v_1 (same color)")
print()
print("Twisted edge:                 e_a connects to v_u_0 AND v_v_1 (crossed)")
print("                              e_b connects to v_u_1 AND v_v_0 (crossed)")
print()
print(f"In this example, base edges (1,2), (2,3), (1,3) are TWISTED.")
print("Red dashed edges in right graph connect to these twisted edge-nodes.")
print("=" * 60)
print(f"Cycle v1->v2->v3->v1 through e_a nodes:")
print(f"  - CFI(K4):   {twisted_count_G1} twisted edge(s) in cycle -> parity EVEN")
print(f"  - CFI'(K4):  {twisted_count_G2} twisted edge(s) in cycle -> parity ODD")
print()
print("This parity difference makes the two graphs non-isomorphic!")
print("3-WL cannot detect this (needs 4-WL), but the cycle parity reveals it.")

# ========================
# 4. Layout and visualization
# ========================
def layout_cfi(G):
    """Manual node coordinates to ensure matching positions for both graphs"""
    pos = {}
    # Place 4 original vertices at the four corners of a square
    base_pos = {
        1: np.array([-1, 1]),
        2: np.array([1, 1]),
        3: np.array([1, -1]),
        4: np.array([-1, -1])
    }
    # Place v0 outside, v1 inside
    for v in vertices:
        center = base_pos[v]
        direction = center / np.linalg.norm(center) * 0.3
        pos[f"v{v}_0"] = center + direction
        pos[f"v{v}_1"] = center - direction

    # Edge nodes at midpoint, with offset to distinguish a and b
    for (u,v) in edges:
        mid = (base_pos[u] + base_pos[v]) / 2
        perp = np.array([-(base_pos[v][1]-base_pos[u][1]),
                          base_pos[v][0]-base_pos[u][0]])
        perp = perp / np.linalg.norm(perp) * 0.15
        pos[f"e{u}{v}_a"] = mid + perp
        pos[f"e{u}{v}_b"] = mid - perp
    return pos

pos = layout_cfi(G1)

# Prepare plotting
fig, axes = plt.subplots(1, 2, figsize=(14, 6))
fig.suptitle("CFI(K4) vs CFI'(K4) - Twisted edges connect different-colored vertex clones", fontsize=14)

# Node colors: light blue for vertex clones, light orange for edge nodes
def node_colors(G):
    colors = []
    for node in G.nodes():
        if node.startswith("v"):
            colors.append("#a6cee3")
        else:
            colors.append("#fdbf6f")
    return colors

# 公共参数
options = {
    "pos": pos,
    "node_color": node_colors(G1),
    "node_size": 400,
    "edgecolors": "black",
    "linewidths": 1.2,
    "with_labels": True,
    "font_size": 7
}

# Plot CFI(K4)
ax1 = axes[0]
ax1.set_title("CFI(K4) - All Parallel")
node_options = {k: v for k, v in options.items() if k not in ['with_labels', 'font_size']}
nx.draw_networkx_nodes(G1, **node_options, ax=ax1)
nx.draw_networkx_edges(G1, pos, edge_color='gray', style='dashed', ax=ax1)
# Highlight the cycle (green thick line)
nx.draw_networkx_edges(G1, pos, edgelist=cycle_path_G1, edge_color='#2ca02c', width=4, ax=ax1)
nx.draw_networkx_labels(G1, pos, font_size=7, ax=ax1)

# Add annotations
ax1.annotate(f"Parity: {twisted_count_G1} (even)", xy=(-0.5, 0.15), fontsize=11, color='#2ca02c',
             ha='center', fontweight='bold')
ax1.annotate("All parallel edges\n(dashed)", xy=(0.7, -0.8), fontsize=9, color='gray',
             ha='center', style='italic')

# Plot CFI'(K4) with twisted edges as red dashed
ax2 = axes[1]
ax2.set_title("CFI'(K4) - Triangle (1,2,3) Twisted")
node_options = {k: v for k, v in options.items() if k not in ['with_labels', 'font_size']}
nx.draw_networkx_nodes(G2, **node_options, ax=ax2)
# Separate twisted edges from normal edges
twisted_edge_list = []
normal_edge_list = []
for (u,v) in G2.edges():
    is_twisted = False
    for n in (u,v):
        if n.startswith("e"):
            parts = n[1:].split("_")[0]
            a = int(parts[0])
            b = int(parts[1]) if len(parts)>1 else int(parts[0])
            if (a,b) in twisted_edges or (b,a) in twisted_edges:
                is_twisted = True
                break
    if is_twisted:
        twisted_edge_list.append((u,v))
    else:
        normal_edge_list.append((u,v))

nx.draw_networkx_edges(G2, pos, edgelist=normal_edge_list, edge_color='red', style='dashed', width=1.5, ax=ax2)
nx.draw_networkx_edges(G2, pos, edgelist=twisted_edge_list, edge_color='black', width=2, ax=ax2)
nx.draw_networkx_labels(G2, pos, font_size=7, ax=ax2)

# Add annotations
ax2.annotate(f"Parity: {twisted_count_G2} (odd)\n(all 3 are twisted)", xy=(-0.5, 0.15), fontsize=11, color='black',
             ha='center', fontweight='bold')
ax2.annotate("Red dashed: parallel edges\nBlack solid: twisted edges", xy=(0.7, -0.8), fontsize=9, color='gray',
             ha='center', style='italic')

plt.tight_layout()
plt.savefig('/media/zhouyu/data/Workspace/Short_Long_Topology_Aware_WGK/cfi_comparison.png', dpi=150, bbox_inches='tight')
print("Figure saved to: /media/zhouyu/data/Workspace/Short_Long_Topology_Aware_WGK/cfi_comparison.png")