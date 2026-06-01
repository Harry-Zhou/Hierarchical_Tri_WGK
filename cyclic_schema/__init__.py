"""
Cyclic Schema - Graph neural network tools for cyclic structure analysis.
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
