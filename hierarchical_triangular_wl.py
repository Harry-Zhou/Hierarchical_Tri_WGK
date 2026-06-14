"""
Hierarchical Triangulated Weisfeiler-Lehman Test (Step-2 & Step-3).

This module provides the main entry point for the hierarchical WL algorithm
that operates on multi-layer cyclic schematic graphs (CSG).

Step-1 (CSG construction and mappings) is delegated to :mod:`cyclic_schema.cyclic_schema`.
Step-2 (single-iteration message passing) and Step-3 (multi-layer iterative
message passing) are implemented by importing from :mod:`cyclic_schema.htn_wl`.

Usage
-----
    import networkx as nx
    import numpy as np
    from hierarchical_triangular_wl import hierarchical_triangular_wl

    G1 = nx.cycle_graph(5)
    G2 = nx.path_graph(5)
    labels = np.array([0, 1, 2, 3, 4], dtype=np.int32)

    wl1, wl2 = hierarchical_triangular_wl(G1, G2, labels, labels, L=1, I=3)
"""

from typing import Any, Dict, List, Optional, Tuple

import networkx as nx
import numpy as np

# ---------------------------------------------------------------------------
# Re-export from cyclic_schema.cyclic_schema  (Step-1: CSG construction)
# ---------------------------------------------------------------------------
from cyclic_schema.cyclic_schema import (
    cyclic_schematic_graph,
    build_csg_to_input_mapping,
    build_input_to_csg_mapping,
    build_multilayer_csg,
    build_multilayer_csg_with_mappings,
    get_cycle_edges,
    get_edge_in_cycle,
    get_node_type,
    get_original_nodes_for_csg_node,
    node_sort_key,
    draw_input_graph,
    draw_cyclic_schematic,
    draw_input_graph_with_cycles,
    draw_side_by_side,
    draw_side_by_side_with_cycles,
    print_analysis_info,
    build_example1_graph,
)

# ---------------------------------------------------------------------------
# Re-export from cyclic_schema.htn_wl  (Step-2 & Step-3: message passing)
# ---------------------------------------------------------------------------
from cyclic_schema.htn_wl import (
    canonicalize_cycle_label,
    compute_initial_label_tuples,
    compute_final_label_tuples,
    forward_aggregate,
    forward_message_passing_both,
    backward_message_passing_both,
    backward_message_passing_both_with_edges,
    hierarchical_triangular_wl,
    hierarchical_triangular_wl_with_edges,
    hierarchical_triangular_wl_unified,
)

# ---------------------------------------------------------------------------
# Re-export multi-layer CSG construction
# ---------------------------------------------------------------------------
from cyclic_schema.multilayer_csg import (
    build_multilayer_csg_single,
    build_multilayer_csg_pair,
    SingleGraphCSG,
    PairedGraphCSG,
)

# ---------------------------------------------------------------------------
# Re-export canonical MCB
# ---------------------------------------------------------------------------
from cyclic_schema.canonical_mcb import (
    canonical_mcb,
    compute_canonical_vertex_ids,
    canonical_mcb_invariant_summary,
)


__all__ = [
    # Step-1: CSG construction
    "cyclic_schematic_graph",
    "build_csg_to_input_mapping",
    "build_input_to_csg_mapping",
    "build_multilayer_csg",
    "build_multilayer_csg_with_mappings",
    "get_cycle_edges",
    "get_edge_in_cycle",
    "get_node_type",
    "get_original_nodes_for_csg_node",
    "node_sort_key",
    # Step-2 & Step-3: message passing
    "canonicalize_cycle_label",
    "compute_initial_label_tuples",
    "compute_final_label_tuples",
    "forward_aggregate",
    "forward_message_passing_both",
    "backward_message_passing_both",
    "backward_message_passing_both_with_edges",
    "hierarchical_triangular_wl",
    "hierarchical_triangular_wl_with_edges",
    "hierarchical_triangular_wl_unified",
    # Multi-layer CSG
    "build_multilayer_csg_single",
    "build_multilayer_csg_pair",
    "SingleGraphCSG",
    "PairedGraphCSG",
    # Canonical MCB
    "canonical_mcb",
    "compute_canonical_vertex_ids",
    "canonical_mcb_invariant_summary",
    # Utilities
    "build_example1_graph",
    "draw_input_graph",
    "draw_cyclic_schematic",
    "draw_input_graph_with_cycles",
    "draw_side_by_side",
    "draw_side_by_side_with_cycles",
    "print_analysis_info",
]
