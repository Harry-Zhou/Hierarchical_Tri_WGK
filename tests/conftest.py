"""Shared fixtures for the test suite."""

import networkx as nx
import numpy as np
import pytest


# ---------------------------------------------------------------------------
# Simple graph fixtures
# ---------------------------------------------------------------------------

@pytest.fixture
def triangle():
    """3-cycle (triangle)."""
    G = nx.Graph()
    G.add_edges_from([(0, 1), (1, 2), (0, 2)])
    return G


@pytest.fixture
def path3():
    """Path graph with 3 nodes."""
    G = nx.Graph()
    G.add_edges_from([(0, 1), (1, 2)])
    return G


@pytest.fixture
def path4():
    """Path graph with 4 nodes."""
    G = nx.Graph()
    G.add_edges_from((i, i + 1) for i in range(3))
    return G


@pytest.fixture
def path5():
    """Path graph with 5 nodes."""
    G = nx.Graph()
    G.add_edges_from((i, i + 1) for i in range(4))
    return G


@pytest.fixture
def cycle6():
    """6-cycle."""
    G = nx.Graph()
    G.add_edges_from((i, (i + 1) % 6) for i in range(6))
    return G


@pytest.fixture
def two_disjoint_triangles():
    """Two disjoint triangles (6 nodes, 6 edges)."""
    G = nx.Graph()
    G.add_edges_from([(0, 1), (1, 2), (0, 2), (3, 4), (4, 5), (3, 5)])
    return G


@pytest.fixture
def chord_cycle():
    """4-cycle with a chord."""
    G = nx.Graph()
    G.add_edges_from([(0, 1), (1, 2), (2, 3), (3, 0), (0, 2)])
    return G


@pytest.fixture
def star():
    """Star graph: center 0 with leaves 1,2."""
    G = nx.Graph()
    G.add_edges_from([(0, 1), (0, 2)])
    return G


@pytest.fixture
def single_node():
    """Single isolated node."""
    G = nx.Graph()
    G.add_node(0)
    return G


@pytest.fixture
def empty_graph():
    """Empty graph."""
    return nx.Graph()


@pytest.fixture
def non_consecutive_graph():
    """Graph with non-consecutive node IDs (0, 5, 10). Default test example."""
    G = nx.Graph()
    G.add_edges_from([(0, 5), (5, 10), (10, 0)])
    return G


# ---------------------------------------------------------------------------
# Label fixtures
# ---------------------------------------------------------------------------

@pytest.fixture
def labels_123():
    """Reusable integer label array [1, 2, 3]."""
    return np.array([1, 2, 3], dtype=np.int32)


@pytest.fixture
def labels_uniform():
    """Uniform label array [5, 5, 5]."""
    return np.array([5, 5, 5], dtype=np.int32)


@pytest.fixture
def seeded_random_labels():
    """Deterministic random labels, shape depends on graph."""
    np.random.seed(42)
    return np.random.randint(0, 10, 26)  # sized for example1


# ---------------------------------------------------------------------------
# Edge label fixtures
# ---------------------------------------------------------------------------

@pytest.fixture
def triangle_elabel_uniform():
    """All triangle edges have label 10."""
    return {(0, 1): 10, (1, 2): 10, (0, 2): 10}


@pytest.fixture
def triangle_elabel_distinct():
    """Triangle edges have distinct labels 10, 20, 30."""
    return {(0, 1): 10, (1, 2): 20, (0, 2): 30}


# ---------------------------------------------------------------------------
# Dataset info fixtures
# ---------------------------------------------------------------------------

@pytest.fixture
def g_info_node_only():
    """Dataset info with node labels only (no edge labels)."""
    return {'nl': True, 'el': False}


@pytest.fixture
def g_info_with_edges():
    """Dataset info with both node and edge labels."""
    return {'nl': True, 'el': True}
