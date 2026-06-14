"""
Cyclic Schema - Graph neural network tools for cyclic structure analysis.

Exports:
  Step-1 CSG construction (cyclic_schema):
    cyclic_schematic_graph, build_multilayer_csg, build_multilayer_csg_with_mappings,
    build_csg_to_input_mapping, build_input_to_csg_mapping, ...

  Multi-layer CSG construction (multilayer_csg):
    build_multilayer_csg_single, build_multilayer_csg_pair,
    SingleGraphCSG, PairedGraphCSG, ...

  Canonical minimum cycle basis (canonical_mcb):
    canonical_mcb, compute_canonical_vertex_ids, canonical_mcb_invariant_summary

  Step-2/Step-3 hierarchical WL (hierarchical_triangular_wl):
    hierarchical_triangular_wl, hierarchical_triangular_wl_with_edges,
    hierarchical_triangular_wl_unified, ...
"""

__version__ = "0.1.0"

from .cyclic_schema import (
    node_sort_key,
    cyclic_schematic_graph,
    build_example1_graph,
    build_multilayer_csg,
    build_multilayer_csg_with_mappings,
    get_node_type,
    get_original_nodes_for_csg_node,
    build_csg_to_input_mapping,
    build_input_to_csg_mapping,
    draw_input_graph,
    draw_cyclic_schematic,
    draw_input_graph_with_cycles,
    draw_side_by_side,
    draw_side_by_side_with_cycles,
    print_analysis_info,
    get_cycle_edges,
    get_edge_in_cycle,
)

from .multilayer_csg import (
    SingleGraphCSG,
    PairedGraphCSG,
    build_multilayer_csg_single,
    build_multilayer_csg_pair,
    get_layer_graph,
    get_cycle_basis,
    get_mapping,
    get_neighbor_components,
    extract_htn_wl_args,
)

from .canonical_mcb import (
    canonical_mcb,
    compute_canonical_vertex_ids,
    canonical_mcb_invariant_summary,
)

from .htn_wl import (
    hierarchical_triangular_wl,
    hierarchical_triangular_wl_with_edges,
    hierarchical_triangular_wl_unified,
    forward_aggregate,
    forward_message_passing_both,
    backward_message_passing_both,
    backward_message_passing_both_with_edges,
    canonicalize_cycle_label,
    compute_initial_label_tuples,
    compute_final_label_tuples,
)

from .csg_transformer import (
    CSGTransformer,
    CSGTransformerLayer,
    TNAAttention,
    SparseGlobalAttention,
    ForwardCrossAttention,
    BackwardCrossAttention,
    LaplacianPositionalEncoding,
    RandomWalkStructuralEncoding,
    StructuralEncoding,
    CompositePositionalEncoding,
    build_model_from_config,
    csg_transformer_unified,
    compute_graph_similarity,
)

# (htn_wl is the canonical implementation module; hierarchical_triangular_wl.py at project root re-exports)
