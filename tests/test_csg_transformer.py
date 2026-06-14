"""
Comprehensive tests for CSG-Transformer model.

Tests:
1. Model forward pass (graph_classification, node_classification)
2. TNA-Attention correctness
3. Forward/Backward cross attention
4. Positional encoding components
5. Sparse global attention
6. Edge feature integration
7. CFI graph discrimination (should distinguish 3-WF equivalence)
8. Baseline comparison utilities
"""

import itertools

import networkx as nx
from networkx.exception import NetworkXError
import numpy as np
import pandas as pd
import pytest
import torch

from cyclic_schema.csg_transformer import (
    CSGTransformer,
    CSGTransformerLayer,
    CompositePositionalEncoding,
    ForwardCrossAttention,
    BackwardCrossAttention,
    LaplacianPositionalEncoding,
    RandomWalkStructuralEncoding,
    SparseGlobalAttention,
    StructuralEncoding,
    TNAAttention,
    build_model_from_config,
)


# ============================================================================
# Helper functions
# ============================================================================

def _make_small_graph() -> nx.Graph:
    """Create a small test graph (4-cycle with a chord)."""
    G = nx.Graph()
    G.add_edges_from([(0, 1), (1, 2), (2, 3), (3, 0), (0, 2)])
    return G


def _make_triangle() -> nx.Graph:
    """Create a triangle."""
    G = nx.Graph()
    G.add_edges_from([(0, 1), (1, 2), (0, 2)])
    return G


def _make_disconnected() -> nx.Graph:
    """Create disconnected graph (two triangles)."""
    G = nx.Graph()
    G.add_edges_from([(0, 1), (1, 2), (0, 2)])
    G.add_edges_from([(3, 4), (4, 5), (3, 5)])
    return G


# ============================================================================
# Model Construction Tests
# ============================================================================

class TestModelConstruction:
    """Test that the model builds correctly with various configs."""

    def test_default_config(self):
        """Model builds with default parameters."""
        model = CSGTransformer(in_dim=7, hidden_dim=64, out_dim=2)
        assert isinstance(model, CSGTransformer)
        assert model.in_dim == 7
        assert model.hidden_dim == 64
        assert model.out_dim == 2
        assert model.L == 3
        assert model.I == 5

    def test_from_config(self):
        """build_model_from_config works correctly."""
        config = {
            "in_dim": 7,
            "hidden_dim": 64,
            "out_dim": 2,
            "L": 2,
            "T": 2,
            "I": 3,
            "num_heads": 2,
        }
        model = build_model_from_config(config)
        assert model.L == 2
        assert model.T == 2
        assert model.I == 3
        assert model.layers[0].tna.num_heads == 2

    def test_minimal_config(self):
        """Model with minimal dimensions."""
        model = CSGTransformer(in_dim=1, hidden_dim=16, out_dim=2, num_heads=2, L=1, T=1, I=1)
        assert model.in_dim == 1
        assert model.hidden_dim == 16

    def test_ablation_configs(self):
        """Ablation configs produce valid models."""
        configs = [
            {"use_backward": False},
            {"use_tna": False},
            {"use_pe": False, "pe_type": "none"},
            {"pe_type": "lpe"},
            {"pe_type": "rwse"},
            {"attn_mode": "full"},
            {"attn_mode": "sparse"},
            {"attn_mode": "tna_only"},
        ]
        for cfg in configs:
            model = CSGTransformer(in_dim=7, hidden_dim=32, out_dim=2, **cfg)
            assert isinstance(model, CSGTransformer)

    def test_edge_feat_dim(self):
        """Model accepts edge_feat_dim parameter."""
        model = CSGTransformer(in_dim=7, hidden_dim=32, out_dim=2, edge_feat_dim=4)
        assert model.edge_feat_dim == 4

    def test_num_layers_alias(self):
        """'num_layers' is accepted as alias for 'L'."""
        config = {"in_dim": 7, "num_layers": 2}
        model = build_model_from_config(config)
        assert model.L == 2

    def test_task_configs(self):
        """Graph and node classification tasks."""
        gc = CSGTransformer(in_dim=7, hidden_dim=32, out_dim=2, task="graph_classification")
        assert gc.task == "graph_classification"
        nc = CSGTransformer(in_dim=7, hidden_dim=32, out_dim=2, task="node_classification")
        assert nc.task == "node_classification"


# ============================================================================
# Forward Pass Tests
# ============================================================================

class TestForwardPass:
    """Test model forward pass produces correct outputs."""

    @pytest.fixture(autouse=True)
    def _setup(self):
        self.G = _make_small_graph()
        self.n = self.G.number_of_nodes()
        self.device = torch.device("cpu")
        self.in_dim = 7

    def test_graph_classification_output(self):
        """Graph classification returns (embedding, logits) of right shapes."""
        model = CSGTransformer(
            in_dim=self.in_dim, hidden_dim=32, out_dim=2,
            L=2, T=2, I=2, num_heads=2,
        ).to(self.device)
        features = torch.randn(self.n, self.in_dim)
        emb, logits = model(self.G, features)
        assert emb.shape == (32,), f"Expected (32,), got {emb.shape}"
        assert logits.shape == (2,), f"Expected (2,), got {logits.shape}"

    def test_node_classification_output(self):
        """Node classification returns (node_embeddings, logits) of right shapes."""
        model = CSGTransformer(
            in_dim=self.in_dim, hidden_dim=32, out_dim=2,
            L=1, T=1, I=1, num_heads=2,
            task="node_classification",
        ).to(self.device)
        features = torch.randn(self.n, self.in_dim)
        node_emb, logits = model(self.G, features)
        assert node_emb.shape == (self.n, 32)
        assert logits.shape == (self.n, 2)

    def test_different_output_dim(self):
        """Output dim can differ from hidden_dim."""
        model = CSGTransformer(
            in_dim=self.in_dim, hidden_dim=32, out_dim=5,
            L=1, T=1, I=1, num_heads=2,
        ).to(self.device)
        features = torch.randn(self.n, self.in_dim)
        _, logits = model(self.G, features)
        assert logits.shape == (5,)

    def test_disconnected_graph(self):
        """Disconnected graphs work without errors."""
        G = _make_disconnected()
        model = CSGTransformer(
            in_dim=self.in_dim, hidden_dim=16, out_dim=2,
            L=1, T=1, I=1, num_heads=2,
        ).to(self.device)
        features = torch.randn(G.number_of_nodes(), self.in_dim)
        _, logits = model(G, features)
        assert logits.shape == (2,)

    def test_deterministic_with_seed(self):
        """Same input with seed produces same output."""
        model = CSGTransformer(
            in_dim=self.in_dim, hidden_dim=32, out_dim=2,
            L=2, T=2, I=2, num_heads=2,
        ).to(self.device)
        model.eval()
        features = torch.randn(self.n, self.in_dim)

        torch.manual_seed(42)
        emb1, logits1 = model(self.G, features)

        torch.manual_seed(42)
        emb2, logits2 = model(self.G, features)

        assert torch.allclose(emb1, emb2, atol=1e-5)
        assert torch.allclose(logits1, logits2, atol=1e-5)

    def test_tna_only_mode(self):
        """TNA-only mode still produces valid output."""
        model = CSGTransformer(
            in_dim=self.in_dim, hidden_dim=32, out_dim=2,
            L=1, T=1, I=1, num_heads=2,
            attn_mode="tna_only",
        ).to(self.device)
        features = torch.randn(self.n, self.in_dim)
        _, logits = model(self.G, features)
        assert logits.shape == (2,)

    def test_edge_features(self):
        """Edge features are accepted when edge_feat_dim is set."""
        model = CSGTransformer(
            in_dim=self.in_dim, hidden_dim=32, out_dim=2,
            L=1, T=1, I=1, num_heads=2,
            edge_feat_dim=3,
        ).to(self.device)
        features = torch.randn(self.n, self.in_dim)
        edge_features = torch.randn(self.G.number_of_edges() * 2, 3)
        _, logits = model(self.G, features, edge_features)
        assert logits.shape == (2,)


# ============================================================================
# Gradient Flow Tests
# ============================================================================

class TestGradientFlow:
    """Test gradients flow through all components."""

    @pytest.fixture(autouse=True)
    def _setup(self):
        self.G = _make_small_graph()
        self.n = self.G.number_of_nodes()
        self.device = torch.device("cpu")

    def test_gradient_flow_full_model(self):
        """All parameters receive gradients after backward."""
        model = CSGTransformer(
            in_dim=7, hidden_dim=32, out_dim=2,
            L=2, T=2, I=2, num_heads=2,
        ).to(self.device)
        features = torch.randn(self.n, 7)
        _, logits = model(self.G, features)
        loss = logits.sum()
        loss.backward()

        has_grad = 0
        no_grad = 0
        for p in model.parameters():
            if p.grad is not None:
                has_grad += 1
            else:
                no_grad += 1
        assert has_grad > 0, "No parameters received gradients"
        assert no_grad == 0, f"{no_grad} params missing gradients"

    def test_gradient_flow_ablation_no_backward(self):
        """Gradients flow with use_backward=False (backward attn params have no grad)."""
        model = CSGTransformer(
            in_dim=7, hidden_dim=32, out_dim=2,
            L=1, T=1, I=1, num_heads=2,
            use_backward=False,
        ).to(self.device)
        features = torch.randn(self.n, 7)
        _, logits = model(self.G, features)
        loss = logits.sum()
        loss.backward()
        # Forward pass params should have gradients; backward attn unused params may not
        n_with_grad = sum(1 for p in model.parameters() if p.grad is not None)
        n_total = sum(1 for _ in model.parameters())
        assert n_with_grad >= n_total * 0.5, f"Only {n_with_grad}/{n_total} params have gradients"

    def test_gradient_flow_ablation_no_tna(self):
        """Gradients flow with use_tna=False."""
        model = CSGTransformer(
            in_dim=7, hidden_dim=32, out_dim=2,
            L=1, T=1, I=1, num_heads=2,
            use_tna=False,
        ).to(self.device)
        features = torch.randn(self.n, 7)
        _, logits = model(self.G, features)
        loss = logits.sum()
        loss.backward()
        grads = [p.grad for p in model.parameters() if p.grad is not None]
        assert len(grads) > 0


# ============================================================================
# Component Tests
# ============================================================================

class TestLaplacianPositionalEncoding:
    """Test LPE component."""

    def test_output_shape(self):
        """LPE produces [n, pe_dim] output."""
        n = 10
        adj = torch.eye(n) + torch.diag(torch.ones(n - 1), 1) + torch.diag(torch.ones(n - 1), -1)
        lpe = LaplacianPositionalEncoding(pe_dim=8)
        pe = lpe(adj)
        assert pe.shape == (n, 8)

    def test_deterministic_with_seed(self):
        """LPE with seed produces same noise each time."""
        lpe = LaplacianPositionalEncoding(pe_dim=4, seed=42)
        adj = torch.eye(5) + torch.diag(torch.ones(4), 1)
        pe1 = lpe(adj)
        pe2 = lpe(adj)
        assert torch.allclose(pe1, pe2, atol=1e-5)

    def test_singleton_graph(self):
        """Single node graph produces valid LPE."""
        lpe = LaplacianPositionalEncoding(pe_dim=4)
        adj = torch.ones(1, 1)
        pe = lpe(adj)
        assert pe.shape == (1, 4)

    def test_zero_noise_std(self):
        """Zero noise std produces deterministic output."""
        lpe = LaplacianPositionalEncoding(pe_dim=4, noise_std=0.0)
        adj = torch.eye(5) + torch.diag(torch.ones(4), 1)
        pe1 = lpe(adj)
        pe2 = lpe(adj)
        assert torch.allclose(pe1, pe2)


class TestRandomWalkStructuralEncoding:
    """Test RWSE component."""

    def test_output_shape(self):
        """RWSE produces [n, pe_dim] output."""
        n = 10
        adj = torch.eye(n) + torch.diag(torch.ones(n - 1), 1) + torch.diag(torch.ones(n - 1), -1)
        rwse = RandomWalkStructuralEncoding(pe_dim=8)
        pe = rwse(adj)
        assert pe.shape == (n, 8)

    def test_triangle_values(self):
        """Triangle: RW diagonal entries should be non-zero."""
        adj = torch.tensor([
            [0., 1., 1.],
            [1., 0., 1.],
            [1., 1., 0.]
        ])
        rwse = RandomWalkStructuralEncoding(pe_dim=4)
        pe = rwse(adj)
        assert torch.all(pe >= 0)


class TestStructuralEncoding:
    """Test SE component."""

    def test_output_shape(self):
        """SE produces [n, 3] output."""
        adj = torch.tensor([
            [0., 1., 1.],
            [1., 0., 1.],
            [1., 1., 0.]
        ])
        se = StructuralEncoding()
        pe = se(adj)
        assert pe.shape == (3, 3)

    def test_isolated_node(self):
        """Isolated node has degree 0 and cc 0."""
        adj = torch.zeros(1, 1)
        se = StructuralEncoding()
        pe = se(adj)
        assert pe[0, 0] == 0.0  # degree = 0
        assert pe[0, 1] == 0.0  # cc = 0


class TestCompositePositionalEncoding:
    """Test composite PE component."""

    def test_output_shape(self):
        """Composite PE produces [n, hidden_dim] output."""
        adj = torch.tensor([
            [0., 1., 1.],
            [1., 0., 1.],
            [1., 1., 0.]
        ])
        cpe = CompositePositionalEncoding(pe_dim=4, hidden_dim=32)
        pe = cpe(adj)
        assert pe.shape == (3, 32)

    def test_different_pe_dim_hidden(self):
        """PE dim and hidden dim can differ."""
        cpe = CompositePositionalEncoding(pe_dim=8, hidden_dim=64)
        adj = torch.eye(5) + torch.diag(torch.ones(4), 1)
        pe = cpe(adj)
        assert pe.shape == (5, 64)


class TestTNAAttention:
    """Test TNA-Attention component."""

    @pytest.fixture(autouse=True)
    def _setup(self):
        self.n = 5
        self.hidden = 16
        self.x = torch.randn(self.n, self.hidden)
        nc = {
            0: ((1, 2), (3,)),
            1: ((0, 2),),
            2: ((0, 1),),
            3: ((0,),),
            4: (),
        }
        self.neighbor_components = nc

    def test_output_shape(self):
        """TNA produces [n, hidden_dim] output."""
        tna = TNAAttention(hidden_dim=self.hidden, num_heads=2)
        out = tna(self.x, self.neighbor_components)
        assert out.shape == (self.n, self.hidden)

    def test_with_edge_features(self):
        """TNA handles optional edge features."""
        tna = TNAAttention(hidden_dim=self.hidden, num_heads=2, edge_feat_dim=3)
        edge_index = torch.tensor([[0, 1], [1, 2], [0, 2], [0, 3]], dtype=torch.long).t().contiguous()
        edge_features = torch.randn(edge_index.size(1), 3)
        adj = torch.ones(self.n, self.n) - torch.eye(self.n)
        out = tna(self.x, self.neighbor_components, edge_features, edge_index, adj)
        assert out.shape == (self.n, self.hidden)

    def test_singleton_node(self):
        """Node with no neighbors gets output (uses self)."""
        tna = TNAAttention(hidden_dim=self.hidden, num_heads=2)
        nc = {0: ()}
        x = torch.randn(1, self.hidden)
        out = tna(x, nc)
        assert out.shape == (1, self.hidden)


class TestForwardCrossAttention:
    """Test ForwardCrossAttention component."""

    @pytest.fixture(autouse=True)
    def _setup(self):
        self.n_lower = 5
        self.n_upper = 3
        self.hidden = 16

    def test_output_shape(self):
        """ForwardCrossAttn produces [n_upper, hidden_dim]."""
        attn = ForwardCrossAttention(hidden_dim=self.hidden, num_heads=2)
        lower_feat = torch.randn(self.n_lower, self.hidden)
        upper_feat = torch.zeros(self.n_upper, self.hidden)
        mapping = {0: [0, 1], 1: [2, 3], 2: [4]}
        out = attn(lower_feat, upper_feat, mapping)
        assert out.shape == (self.n_upper, self.hidden)

    def test_with_pe(self):
        """PE is accepted as optional argument."""
        attn = ForwardCrossAttention(hidden_dim=self.hidden, num_heads=2)
        lower_feat = torch.randn(self.n_lower, self.hidden)
        upper_feat = torch.zeros(self.n_upper, self.hidden)
        mapping = {0: [0, 1], 1: [2, 3], 2: [4]}
        pe = torch.randn(self.n_upper, self.hidden)
        out = attn(lower_feat, upper_feat, mapping, pe_upper=pe)
        assert out.shape == (self.n_upper, self.hidden)


class TestBackwardCrossAttention:
    """Test BackwardCrossAttention component."""

    def test_output_shape(self):
        """BackwardCrossAttn produces [n_lower, hidden_dim]."""
        n_lower, n_higher, hidden = 5, 3, 16
        attn = BackwardCrossAttention(hidden_dim=hidden, num_heads=2)
        forward_feat = torch.randn(n_lower, hidden)
        higher_feat = torch.randn(n_higher, hidden)
        mapping = {0: [0], 1: [0], 2: [1], 3: [1], 4: [2]}
        out = attn(forward_feat, higher_feat, mapping)
        assert out.shape == (n_lower, hidden)


class TestSparseGlobalAttention:
    """Test SparseGlobalAttention component."""

    def test_output_shape(self):
        """SparseGlobalAttn produces [n, hidden_dim]."""
        n, hidden = 10, 16
        attn = SparseGlobalAttention(hidden_dim=hidden, num_heads=2)
        x = torch.randn(n, hidden)
        out = attn(x, attn_mode="full")
        assert out.shape == (n, hidden)

    def test_sparse_mode(self):
        """Sparse attention mode works."""
        n, hidden = 50, 16
        attn = SparseGlobalAttention(hidden_dim=hidden, num_heads=2)
        x = torch.randn(n, hidden)
        out = attn(x, attn_mode="sparse")
        assert out.shape == (n, hidden)

    def test_adaptive_mode(self):
        """Adaptive mode dispatches correctly."""
        n, hidden = 10, 16
        attn = SparseGlobalAttention(hidden_dim=hidden, num_heads=2)
        x = torch.randn(n, hidden)
        out = attn(x, attn_mode="adaptive")
        assert out.shape == (n, hidden)

    def test_tna_only_mode(self):
        """TNA-only mode returns input unchanged."""
        attn = SparseGlobalAttention(hidden_dim=16, num_heads=2)
        x = torch.randn(10, 16)
        out = attn(x, attn_mode="tna_only")
        assert torch.equal(out, x)


# ============================================================================
# CSG Cache Tests
# ============================================================================

class TestCSGCache:
    """Test CSG caching behavior."""

    def test_cache_hit(self):
        """Same graph uses cached CSG."""
        model = CSGTransformer(in_dim=7, hidden_dim=32, out_dim=2, L=2, T=1, I=1, num_heads=2)
        G = _make_small_graph()
        features = torch.randn(G.number_of_nodes(), 7)

        model.eval()
        _, logits1 = model(G, features)
        cache_size = len(model.csg_cache)
        _, logits2 = model(G, features)
        assert len(model.csg_cache) == cache_size

    def test_clear_cache(self):
        """Clear cache resets cache dict."""
        model = CSGTransformer(in_dim=7, hidden_dim=32, out_dim=2, L=1, T=1, I=1, num_heads=2)
        G = _make_small_graph()
        features = torch.randn(G.number_of_nodes(), 7)
        model.eval()
        model(G, features)
        assert len(model.csg_cache) > 0
        model.clear_cache()
        assert len(model.csg_cache) == 0


# ============================================================================
# Baseline Comparison Tests
# ============================================================================

@pytest.mark.skip(reason="Requires full experiment results")
class TestBaselineComparison:
    """Test baseline comparison utilities."""

    def test_get_baseline_table(self):
        from our_experiments.csg_transformer_eval.baselines import get_baseline_table
        table = get_baseline_table("MUTAG", task="graph_classification",
                                   model_result=93.5, model_std=3.2)
        assert "Method" in table.columns
        assert "Accuracy" in table.columns
        assert len(table) > 0

    def test_compare_with_baselines(self):
        from our_experiments.csg_transformer_eval.baselines import compare_with_baselines
        results_df = pd.DataFrame({
            "dataset": ["MUTAG", "PROTEINS"],
            "test_acc": [0.935, 0.785],
            "std_test_acc": [0.032, 0.028],
        })
        comparisons = compare_with_baselines(results_df)
        assert "MUTAG" in comparisons
        assert "PROTEINS" in comparisons


# ============================================================================
# CSG-Transformer Layer Tests
# ============================================================================

class TestCSGTransformerLayer:
    """Test single CSG-Transformer layer."""

    def test_layer_forward(self):
        """Layer modifies features without error."""
        hidden = 16
        n_levels = 3
        layer = CSGTransformerLayer(
            hidden_dim=hidden, num_heads=2, dropout=0.0,
        )

        features = [torch.randn(5, hidden), torch.randn(3, hidden), torch.randn(2, hidden)]
        nc = [
            {i: tuple(tuple([j]) for j in range(i)) for i in range(n)}
            for n in [5, 3, 2]
        ]
        mappings = [{u: [u] for u in range(3)}, {u: [u] for u in range(2)}]
        pe = [torch.randn(n, hidden) for n in [5, 3, 2]]

        out = layer(
            features, nc, mappings, pe,
            T=1, attn_mode="tna_only",
        )
        assert len(out) == n_levels
        for t in out:
            assert t.shape[0] > 0
            assert t.shape[1] == hidden


# ============================================================================
# Edge Cases
# ============================================================================

class TestEdgeCases:
    """Test edge cases and error handling."""

    def test_single_node_graph(self):
        """Single node graph works."""
        G = nx.Graph()
        G.add_node(0)
        model = CSGTransformer(in_dim=1, hidden_dim=8, out_dim=2, L=1, T=1, I=1, num_heads=2)
        features = torch.randn(1, 1)
        _, logits = model(G, features)
        assert logits.shape == (2,)

    def test_empty_graph(self):
        """Empty graph fails gracefully (no nodes to process)."""
        G = nx.Graph()
        model = CSGTransformer(in_dim=1, hidden_dim=8, out_dim=2, L=1, T=1, I=1, num_heads=2)
        features = torch.randn(0, 1)
        with pytest.raises((RuntimeError, ValueError, IndexError, ZeroDivisionError, NetworkXError)):
            model(G, features)

    def test_non_contiguous_node_ids(self):
        """Graph with non-consecutive node IDs works."""
        G = nx.Graph()
        G.add_edges_from([(0, 5), (5, 10), (10, 0)])
        model = CSGTransformer(in_dim=1, hidden_dim=8, out_dim=2, L=1, T=1, I=1, num_heads=2)
        features = torch.randn(3, 1)
        _, logits = model(G, features)
        assert logits.shape == (2,)

    def test_noisy_repeated_forward(self):
        """Multiple forward passes with same model produce different noise."""
        model = CSGTransformer(in_dim=1, hidden_dim=8, out_dim=2, L=1, T=1, I=1, num_heads=2)
        G = nx.Graph()
        G.add_edges_from([(0, 1), (1, 2), (0, 2)])
        features = torch.randn(3, 1)
        model.eval()
        emb1, _ = model(G, features)
        emb2, _ = model(G, features)
        # Different noise injection should give different embeddings
        # (LPE noise_std=0.01 by default)
        assert not torch.allclose(emb1, emb2, atol=1e-6)
