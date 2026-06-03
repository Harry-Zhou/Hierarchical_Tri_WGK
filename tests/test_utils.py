"""Unit tests for utils.py."""

import numpy as np
import pytest

from utils import (
    get_start_end_indices,
    get_recommended_nproc,
    precomp_node_neighs,
)


class TestGetStartEndIndices:
    """Tests for multiprocessing chunking utility."""

    def test_single_sample(self):
        indices = get_start_end_indices(1, num_cells=4)
        assert indices == [(0, 1)]

    def test_fewer_samples_than_cells(self):
        indices = get_start_end_indices(3, num_cells=8)
        assert len(indices) == 3
        for start, end in indices:
            assert end - start == 1
        assert indices[0] == (0, 1)
        assert indices[-1] == (2, 3)

    def test_exact_multiple(self):
        indices = get_start_end_indices(8, num_cells=4)
        assert len(indices) == 4
        for start, end in indices:
            assert end - start == 2
        assert indices[0] == (0, 2)
        assert indices[-1] == (6, 8)

    def test_with_remainder(self):
        indices = get_start_end_indices(10, num_cells=4)
        assert len(indices) == 4
        total = sum(end - start for start, end in indices)
        assert total == 10

    def test_uneven_distribution(self):
        """When samples not evenly divisible, first cells get +1."""
        indices = get_start_end_indices(10, num_cells=3)
        assert len(indices) == 3
        # 10/3 = 3 remainder 1 → first 1 cell gets 4, rest get 3
        sizes = [end - start for start, end in indices]
        assert sizes == [4, 3, 3]

    def test_no_overlap(self):
        """All partitions must be disjoint and cover full range."""
        indices = get_start_end_indices(37, num_cells=8)
        seen = set()
        for start, end in indices:
            for i in range(start, end):
                assert i not in seen, f"Overlap at index {i}"
                seen.add(i)
        assert len(seen) == 37

    def test_single_cell(self):
        indices = get_start_end_indices(100, num_cells=1)
        assert indices == [(0, 100)]

    def test_many_cells(self):
        indices = get_start_end_indices(5, num_cells=100)
        assert len(indices) == 5
        assert indices[-1] == (4, 5)


class TestPrecompNodeNeighs:
    """Tests for neighbor precomputation."""

    def test_triangle(self, triangle):
        neighs = precomp_node_neighs(triangle)
        assert len(neighs) == 3
        assert set(neighs[0]) == {1, 2}

    def test_single_node(self, single_node):
        neighs = precomp_node_neighs(single_node)
        assert neighs[0] == []

    def test_non_consecutive_ids(self, non_consecutive_graph):
        neighs = precomp_node_neighs(non_consecutive_graph)
        assert set(neighs[0]) == {5, 10}
        assert set(neighs[5]) == {0, 10}

    def test_empty_graph(self, empty_graph):
        neighs = precomp_node_neighs(empty_graph)
        assert neighs == {}


class TestGetRecommendedNproc:
    """Smoke tests for resource-aware nproc recommendation."""

    def test_returns_positive_int(self):
        nproc = get_recommended_nproc()
        assert isinstance(nproc, int)
        assert nproc >= 1

    def test_respects_reserve_cores(self):
        nproc_default = get_recommended_nproc(reserve_cores=1)
        nproc_reserved = get_recommended_nproc(reserve_cores=4)
        # May be equal if system has few cores; should never increase
        assert nproc_reserved <= nproc_default

    def test_with_mem_constraint(self):
        nproc = get_recommended_nproc(est_mem_per_proc_gb=999999)
        assert nproc >= 1
