"""
CSG-Transformer: Graph Representation Learning with Cyclic Schematic Graph and Transformer.

This module implements the CSG-Transformer model as described in the specification
(csg_transformer方案.md). It combines multi-layer Cyclic Schematic Graph (CSG)
abstraction with Transformer architecture for graph representation learning.

Key features:
1. Triangulated Neighborhood Attention (TNA) with connected component aggregation
   and edge feature integration
2. Forward/Backward cross attention for multi-scale fusion
3. Multi-layer composite positional encoding (LPE + RWSE + Structural)
4. Sparse global attention with adaptive strategy (full / sparse / TNA-only)
5. Edge label support in TNA-Attention
6. CSG caching for efficient repeated forward passes
7. Unified interface for graph pairs

References:
- Specification: cyclic_schema/csg_transformer方案.md
- Infrastructure: cyclic_schema/multilayer_csg.py
"""

import math
from typing import Any, Dict, Hashable, List, Optional, Tuple, Union

import networkx as nx
import numpy as np
import torch
import torch.nn as nn
import torch.nn.functional as F

from .multilayer_csg import (
    SingleGraphCSG,
    build_multilayer_csg_single,
    get_cycle_basis,
    get_layer_graph,
    get_mapping,
    get_neighbor_components,
)


# ============================================================================
# Positional Encoding Components
# ============================================================================

class LaplacianPositionalEncoding(nn.Module):
    """
    Laplacian Positional Encoding (LPE) using eigenvectors of normalized Laplacian.
    
    Computes LPE with noise injection for symmetry breaking as described in
    Section 3.3.1 of the specification.
    """
    
    def __init__(self, pe_dim: int = 16, noise_std: float = 0.01, seed: Optional[int] = None):
        """
        Args:
            pe_dim: Output dimension for positional encoding
            noise_std: Standard deviation for Gaussian noise (default: 0.01)
            seed: Optional seed for reproducible noise injection
        """
        super().__init__()
        self.pe_dim = pe_dim
        self.noise_std = noise_std
        self.seed = seed
        # Pre-create noise generator if seed is provided for reproducibility
        self._noise_generator: Optional[torch.Generator] = None
        if seed is not None:
            self._noise_generator = torch.Generator()
            self._noise_generator.manual_seed(seed)
    
    def forward(self, adj_matrix: torch.Tensor) -> torch.Tensor:
        """
        Compute LPE for all nodes in the graph.
        
        Args:
            adj_matrix: Adjacency matrix [n, n]
            
        Returns:
            pe: Positional encoding [n, pe_dim]
        """
        n = adj_matrix.size(0)
        device = adj_matrix.device
        
        # Compute degree matrix
        deg = adj_matrix.sum(dim=-1)
        deg_inv_sqrt = torch.where(deg > 0, deg.pow(-0.5), torch.zeros_like(deg))
        D_inv_sqrt = torch.diag(deg_inv_sqrt)
        
        # Normalized Laplacian: L = I - D^{-1/2} A D^{-1/2}
        I = torch.eye(n, device=device)
        L = I - D_inv_sqrt @ adj_matrix @ D_inv_sqrt
        
        # Compute eigenvectors
        try:
            eigenvalues, eigenvectors = torch.linalg.eigh(L)
        except Exception:
            # Fallback for singular matrices
            eigenvalues, eigenvectors = torch.linalg.eigh(L + 1e-6 * I)
        
        # Sort by eigenvalue magnitude
        sorted_indices = torch.argsort(eigenvalues)
        
        # Take pe_dim smallest non-zero eigenvalues
        # Skip the first eigenvector (corresponding to eigenvalue 0)
        pe_indices = sorted_indices[1:min(self.pe_dim + 1, n)]
        
        if len(pe_indices) < self.pe_dim:
            # Pad with zeros if not enough eigenvectors
            pe = torch.zeros(n, self.pe_dim, device=device)
            pe[:, :len(pe_indices)] = eigenvectors[:, pe_indices]
        else:
            pe = eigenvectors[:, pe_indices]
        
        # Add noise for symmetry breaking (deterministic if seed is set)
        if self._noise_generator is not None and self.seed is not None:
            noise = torch.randn(pe.shape, device=device, generator=self._noise_generator) * self.noise_std
            # Re-seed to produce same noise on next forward (reproducibility)
            self._noise_generator.manual_seed(self.seed)
        elif self.seed is not None:
            noise = torch.randn_like(pe) * self.noise_std
        else:
            noise = torch.randn_like(pe) * self.noise_std
        pe = pe + noise
        
        return pe


class RandomWalkStructuralEncoding(nn.Module):
    """
    Random Walk Structural Encoding (RWSE) as described in Section 3.3.2.
    
    Computes P^k(v,v) for k=1..p_RW, capturing local random walk features.
    """
    
    def __init__(self, pe_dim: int = 16):
        """
        Args:
            pe_dim: Output dimension for positional encoding
        """
        super().__init__()
        self.pe_dim = pe_dim
    
    def forward(self, adj_matrix: torch.Tensor) -> torch.Tensor:
        """
        Compute RWSE for all nodes.
        
        Args:
            adj_matrix: Adjacency matrix [n, n]
            
        Returns:
            pe: Positional encoding [n, pe_dim]
        """
        n = adj_matrix.size(0)
        device = adj_matrix.device
        
        # Compute degree matrix and transition probability matrix P = D^{-1} A
        deg = adj_matrix.sum(dim=-1)
        deg_inv = torch.where(deg > 0, deg.pow(-1.0), torch.zeros_like(deg))
        D_inv = torch.diag(deg_inv)
        P = D_inv @ adj_matrix
        
        # Compute P^k for k=1..pe_dim
        pe = torch.zeros(n, self.pe_dim, device=device)
        P_k = torch.eye(n, device=device)
        
        for k in range(1, self.pe_dim + 1):
            P_k = P_k @ P
            # Extract diagonal (P^k(v,v))
            pe[:, k - 1] = torch.diag(P_k)
        
        return pe


class StructuralEncoding(nn.Module):
    """
    Structural Encoding based on degree and clustering coefficient.
    
    Computes PE_SE(v) = (deg(v), cc(v), deg^2(v)) as described in Section 3.3.3.
    """
    
    def __init__(self):
        super().__init__()
    
    def forward(self, adj_matrix: torch.Tensor) -> torch.Tensor:
        """
        Compute structural encoding for all nodes.
        
        Args:
            adj_matrix: Adjacency matrix [n, n]
            
        Returns:
            pe: Positional encoding [n, 3]
        """
        n = adj_matrix.size(0)
        device = adj_matrix.device
        
        # Degree (normalized)
        deg = adj_matrix.sum(dim=-1)
        deg_norm = deg / (n - 1) if n > 1 else deg
        
        # Clustering coefficient
        cc = torch.zeros(n, device=device)
        adj_bool = adj_matrix > 0
        for v in range(n):
            neighbors = torch.where(adj_bool[v])[0]
            k = len(neighbors)
            if k < 2:
                cc[v] = 0.0
            else:
                # Count triangles using matrix ops for efficiency
                sub_adj = adj_bool[neighbors][:, neighbors]
                triangles = sub_adj.sum().item() // 2
                cc[v] = 2.0 * triangles / (k * (k - 1))
        
        # Stack features
        pe = torch.stack([deg_norm, cc, deg_norm ** 2], dim=-1)
        
        return pe


class CompositePositionalEncoding(nn.Module):
    """
    Multi-layer composite positional encoding combining LPE, RWSE, and SE.
    
    PE^(l)(v) = MLP_PE([PE_LPE^(l)(v) || PE_RW^(l)(v) || PE_SE^(l)(v)])
    """
    
    def __init__(self, pe_dim: int = 16, hidden_dim: int = 128):
        """
        Args:
            pe_dim: Dimension for each PE component
            hidden_dim: Output dimension
        """
        super().__init__()
        self.pe_dim = pe_dim
        self.lpe = LaplacianPositionalEncoding(pe_dim)
        self.rwse = RandomWalkStructuralEncoding(pe_dim)
        self.se = StructuralEncoding()
        
        # MLP to project concatenated PE to hidden_dim
        input_dim = pe_dim * 2 + 3  # LPE + RWSE + SE
        self.mlp = nn.Sequential(
            nn.Linear(input_dim, hidden_dim),
            nn.ReLU(),
            nn.Linear(hidden_dim, hidden_dim)
        )
    
    def forward(self, adj_matrix: torch.Tensor) -> torch.Tensor:
        """
        Compute composite positional encoding.
        
        Args:
            adj_matrix: Adjacency matrix [n, n]
            
        Returns:
            pe: Positional encoding [n, hidden_dim]
        """
        lpe_pe = self.lpe(adj_matrix)
        rwse_pe = self.rwse(adj_matrix)
        se_pe = self.se(adj_matrix)
        
        # Concatenate and project
        combined = torch.cat([lpe_pe, rwse_pe, se_pe], dim=-1)
        pe = self.mlp(combined)
        
        return pe


# ============================================================================
# TNA-Attention
# ============================================================================

class TNAAttention(nn.Module):
    """
    Triangulated Neighborhood Attention (TNA) as described in Section 3.5.1.
    
    For each node v:
    1. Compute connected components of induced subgraph G[N(v)]
    2. Compute component mean: m_R = (1/|R|) * sum_{u in R} h_u
    3. Attention over components using scaled dot-product:
       alpha_{v,R} = softmax(LeakyReLU( (W_Q h_v)^T (W_K m_R) / sqrt(d_k) ))
    4. Output: h_v' = LayerNorm(h_v + Dropout(sum_R alpha_{v,R} W_V m_R))
    
    Edge features (if provided) are integrated into component means:
    for each component R, the mean edge features of edges inside R are
    concatenated with the mean node features.
    """
    
    def __init__(self, hidden_dim: int, num_heads: int = 4, dropout: float = 0.1,
                 edge_feat_dim: Optional[int] = None):
        """
        Args:
            hidden_dim: Hidden dimension
            num_heads: Number of attention heads
            dropout: Dropout rate
            edge_feat_dim: Dimension of edge features (None if no edge features)
        """
        super().__init__()
        self.hidden_dim = hidden_dim
        self.num_heads = num_heads
        self.head_dim = hidden_dim // num_heads
        assert hidden_dim % num_heads == 0, "hidden_dim must be divisible by num_heads"
        
        # Projections for query, key, value
        self.W_Q = nn.Linear(hidden_dim, hidden_dim)
        self.W_K = nn.Linear(hidden_dim, hidden_dim)
        self.W_V = nn.Linear(hidden_dim, hidden_dim)
        
        # Edge feature integration (optional)
        self.edge_feat_dim = edge_feat_dim
        if edge_feat_dim is not None:
            self.edge_proj = nn.Linear(edge_feat_dim, hidden_dim // 4)
            # Adjust K/V projections to account for concatenated edge features
            self.comp_feat_dim = hidden_dim + hidden_dim // 4
            self.W_K_comp = nn.Linear(self.comp_feat_dim, hidden_dim)
            self.W_V_comp = nn.Linear(self.comp_feat_dim, hidden_dim)
        else:
            self.edge_proj = None
        
        # LayerNorm and Dropout
        self.norm = nn.LayerNorm(hidden_dim)
        self.dropout = nn.Dropout(dropout)
        
        # FFN
        self.ffn = nn.Sequential(
            nn.Linear(hidden_dim, hidden_dim * 4),
            nn.GELU(),
            nn.Dropout(dropout),
            nn.Linear(hidden_dim * 4, hidden_dim),
            nn.Dropout(dropout)
        )
        self.ffn_norm = nn.LayerNorm(hidden_dim)
    
    def forward(
        self,
        x: torch.Tensor,
        neighbor_components: Dict[int, Tuple[Tuple[int, ...], ...]],
        edge_features: Optional[torch.Tensor] = None,
        edge_index: Optional[torch.Tensor] = None,
        adjacency_matrix: Optional[torch.Tensor] = None,
    ) -> torch.Tensor:
        """
        Args:
            x: Node features [n, hidden_dim]
            neighbor_components: Dict mapping node -> tuple of component tuples
            edge_features: Optional edge features [m, edge_dim]
            edge_index: Optional edge index [2, m] for looking up edge features by nodes
            adjacency_matrix: Optional adjacency matrix [n, n] for edge lookup
            
        Returns:
            Updated node features [n, hidden_dim]
        """
        n = x.size(0)
        device = x.device
        
        # Compute queries for all nodes
        Q_full = self.W_Q(x)  # [n, hidden_dim]
        K_full = self.W_K(x)
        V_full = self.W_V(x)
        
        # Reshape for multi-head
        Q_full = Q_full.view(n, self.num_heads, self.head_dim)
        K_full = K_full.view(n, self.num_heads, self.head_dim)
        V_full = V_full.view(n, self.num_heads, self.head_dim)
        
        # Build edge lookup if edge features provided with edge_index
        edge_feat_map = None
        if edge_features is not None and edge_index is not None and adjacency_matrix is not None:
            # Build (u,v) -> edge_feat_idx lookup
            edge_feat_map = {}
            for e_idx in range(edge_index.size(1)):
                u = edge_index[0, e_idx].item()
                v = edge_index[1, e_idx].item()
                edge_feat_map[(u, v)] = e_idx
                edge_feat_map[(v, u)] = e_idx
        
        # Compute TNA attention
        attn_output = torch.zeros_like(x)
        
        for v in range(n):
            components = neighbor_components.get(v, ())
            
            if len(components) == 0:
                # No neighbors, just use self
                attn_output[v] = x[v]
                continue
            
            # Compute component representations
            comp_representations = []
            for comp in components:
                if len(comp) == 0:
                    continue
                comp_indices = list(comp)
                comp_tensor = torch.tensor(comp_indices, device=device, dtype=torch.long)
                
                # Start with mean of node features in component
                comp_node_mean = x[comp_tensor].mean(dim=0)  # [hidden_dim]
                comp_repr = comp_node_mean
                
                # Integrate edge features if available: concatenate mean edge features
                if (edge_features is not None and edge_feat_map is not None
                        and len(comp_indices) >= 2 and self.edge_proj is not None):
                    # Collect edge features for edges inside this component
                    comp_edge_feats = []
                    for i_idx in range(len(comp_indices)):
                        for j_idx in range(i_idx + 1, len(comp_indices)):
                            u_node = comp_indices[i_idx]
                            v_node = comp_indices[j_idx]
                            key = (u_node, v_node)
                            if key in edge_feat_map:
                                e_idx = edge_feat_map[key]
                                comp_edge_feats.append(edge_features[e_idx])
                    if comp_edge_feats:
                        comp_edge_mean = torch.stack(comp_edge_feats).mean(dim=0)
                        comp_edge_proj = self.edge_proj(comp_edge_mean)
                        comp_repr = torch.cat([comp_node_mean, comp_edge_proj], dim=-1)
                    else:
                        pad = torch.zeros(self.hidden_dim // 4, device=device)
                        comp_repr = torch.cat([comp_node_mean, pad], dim=-1)
                else:
                    if self.edge_proj is not None:
                        pad = torch.zeros(self.hidden_dim // 4, device=device)
                        comp_repr = torch.cat([comp_node_mean, pad], dim=-1)
                
                comp_representations.append(comp_repr)
            
            if len(comp_representations) == 0:
                attn_output[v] = x[v]
                continue
            
            # Stack component representations
            comp_stack = torch.stack(comp_representations)  # [num_comp, repr_dim]
            
            # Project component representations through K/V
            if self.edge_proj is not None:
                # Use component-specific projections
                comp_k = self.W_K_comp(comp_stack)  # [num_comp, hidden_dim]
                comp_v = self.W_V_comp(comp_stack)
            else:
                comp_k = self.W_K(comp_stack)
                comp_v = self.W_V(comp_stack)
            
            # Reshape for multi-head
            comp_k = comp_k.view(-1, self.num_heads, self.head_dim)  # [num_comp, num_heads, head_dim]
            comp_v = comp_v.view(-1, self.num_heads, self.head_dim)
            
            # Compute scaled dot-product attention scores
            # score_R = (W_Q h_v)^T (W_K m_R) / sqrt(d_k)
            q_v = Q_full[v].unsqueeze(0)  # [1, num_heads, head_dim]
            
            # [1, num_heads, head_dim] x [num_comp, num_heads, head_dim]^T
            # -> [num_heads, 1, head_dim] x [num_heads, head_dim, num_comp]
            scores = torch.matmul(q_v.permute(1, 0, 2), comp_k.permute(1, 2, 0))  # [num_heads, 1, num_comp]
            scores = scores.squeeze(1)  # [num_heads, num_comp]
            scores = scores / math.sqrt(self.head_dim)
            scores = F.leaky_relu(scores, negative_slope=0.2)
            
            # Softmax over components (dim=1 = num_comp)
            attn_weights = F.softmax(scores, dim=1)  # [num_heads, num_comp]
            attn_weights = self.dropout(attn_weights)
            
            # Weighted sum of component values
            # [num_heads, num_comp] x [num_comp, num_heads, head_dim] -> [num_heads, head_dim]
            attn_out = torch.einsum('hc,chd->hd', attn_weights, comp_v)
            attn_out = attn_out.contiguous().view(-1)  # [hidden_dim]
            
            attn_output[v] = attn_out
        
        # Residual connection and LayerNorm
        attn_output = self.norm(x + self.dropout(attn_output))
        
        # FFN
        ffn_out = self.ffn(attn_output)
        output = self.ffn_norm(attn_output + ffn_out)
        
        return output


# ============================================================================
# Cross Attention Layers
# ============================================================================

class ForwardCrossAttention(nn.Module):
    """
    Forward Cross Attention as described in Section 3.5.2.
    
    Aggregates information from lower level to higher level in CSG hierarchy.
    For each upper node u:
    - beta_v = softmax(q_u^T W_K h_v^{(l)}) for v in down(u)
    - x_u = sum_v beta_v W_V h_v^{(l)}
    - h_u^{(l+1)} = LayerNorm(h_u^{(l)} + Dropout(MLP_fwd([x_u || PE^{(l+1)}(u)])))
    """
    
    def __init__(self, hidden_dim: int, num_heads: int = 4, dropout: float = 0.1):
        """
        Args:
            hidden_dim: Hidden dimension
            num_heads: Number of attention heads
            dropout: Dropout rate
        """
        super().__init__()
        self.hidden_dim = hidden_dim
        self.num_heads = num_heads
        self.head_dim = hidden_dim // num_heads
        
        # Projections
        self.W_Q = nn.Linear(hidden_dim, hidden_dim)
        self.W_K = nn.Linear(hidden_dim, hidden_dim)
        self.W_V = nn.Linear(hidden_dim, hidden_dim)
        
        # Attention
        self.attention = nn.MultiheadAttention(
            hidden_dim, num_heads, dropout=dropout, batch_first=True
        )
        
        # LayerNorm and Dropout
        self.norm1 = nn.LayerNorm(hidden_dim)
        self.norm2 = nn.LayerNorm(hidden_dim)
        self.dropout = nn.Dropout(dropout)
        
        # FFN
        self.ffn = nn.Sequential(
            nn.Linear(hidden_dim, hidden_dim * 4),
            nn.GELU(),
            nn.Dropout(dropout),
            nn.Linear(hidden_dim * 4, hidden_dim),
            nn.Dropout(dropout)
        )
        
        # PE projection
        self.pe_proj = nn.Linear(hidden_dim, hidden_dim)
    
    def forward(
        self,
        lower_features: torch.Tensor,
        upper_features: torch.Tensor,
        down_mapping: Dict[int, List[int]],
        pe_upper: Optional[torch.Tensor] = None
    ) -> torch.Tensor:
        """
        Args:
            lower_features: Features from lower level [n_lower, hidden_dim]
            upper_features: Features from upper level [n_upper, hidden_dim]
            down_mapping: Mapping from upper node (index) to list of lower node indices
            pe_upper: Positional encoding for upper level [n_upper, hidden_dim]
            
        Returns:
            Updated upper features [n_upper, hidden_dim]
        """
        n_upper = upper_features.size(0)
        n_lower = lower_features.size(0)
        device = upper_features.device
        
        # Build attention mask for cross attention
        # Each upper node attends only to its mapped lower nodes
        attn_mask = torch.full((n_upper, n_lower), float('-inf'), device=device)
        
        for u in range(n_upper):
            lower_indices = down_mapping.get(u, [])
            if len(lower_indices) > 0:
                indices = torch.tensor(lower_indices, device=device, dtype=torch.long)
                attn_mask[u, indices] = 0.0
        
        # Cross attention: upper attends to lower
        Q = self.W_Q(upper_features).unsqueeze(0)  # [1, n_upper, hidden_dim]
        K = self.W_K(lower_features).unsqueeze(0)  # [1, n_lower, hidden_dim]
        V = self.W_V(lower_features).unsqueeze(0)  # [1, n_lower, hidden_dim]
        
        # Remove NaN/inf rows in mask (nodes with no lower connections)
        # These nodes will attend uniformly to all lower nodes
        valid_mask = ~torch.isinf(attn_mask).all(dim=1)
        if not valid_mask.all():
            # For rows where all entries are -inf, set to 0 (attend uniformly)
            uniform_rows = torch.isinf(attn_mask).all(dim=1)
            attn_mask[uniform_rows] = 0.0
        
        attn_out, _ = self.attention(Q, K, V, attn_mask=attn_mask)
        attn_out = attn_out.squeeze(0)  # [n_upper, hidden_dim]
        
        # Residual connection and LayerNorm
        upper_features = self.norm1(upper_features + self.dropout(attn_out))
        
        # Add positional encoding if provided
        if pe_upper is not None:
            upper_features = upper_features + self.pe_proj(pe_upper)
        
        # FFN
        ffn_out = self.ffn(upper_features)
        output = self.norm2(upper_features + ffn_out)
        
        return output


class BackwardCrossAttention(nn.Module):
    """
    Backward Cross Attention as described in Section 3.5.3.
    
    Aggregates information from higher level to lower level.
    For each lower node v:
    - gamma_u = softmax((W_Q' h_v^{(l-1)})^T (W_K' h_u^{(l)})) for u in up(v)
    - b_v = sum_u gamma_u W_V' h_u^{(l)}
    - h_v^{(l-1),new} = LayerNorm(h_v^{(l-1)} + Dropout(MLP_bwd([h_v^{(l-1)} || b_v])))
    """
    
    def __init__(self, hidden_dim: int, num_heads: int = 4, dropout: float = 0.1):
        """
        Args:
            hidden_dim: Hidden dimension
            num_heads: Number of attention heads
            dropout: Dropout rate
        """
        super().__init__()
        self.hidden_dim = hidden_dim
        self.num_heads = num_heads
        self.head_dim = hidden_dim // num_heads
        
        # Projections
        self.W_Q = nn.Linear(hidden_dim, hidden_dim)
        self.W_K = nn.Linear(hidden_dim, hidden_dim)
        self.W_V = nn.Linear(hidden_dim, hidden_dim)
        
        # Attention
        self.attention = nn.MultiheadAttention(
            hidden_dim, num_heads, dropout=dropout, batch_first=True
        )
        
        # Fusion MLP for combining forward and backward features
        self.fusion_mlp = nn.Sequential(
            nn.Linear(hidden_dim * 2, hidden_dim * 4),
            nn.GELU(),
            nn.Dropout(dropout),
            nn.Linear(hidden_dim * 4, hidden_dim),
            nn.Dropout(dropout)
        )
        
        # LayerNorm and Dropout
        self.norm1 = nn.LayerNorm(hidden_dim)
        self.norm2 = nn.LayerNorm(hidden_dim)
        self.dropout = nn.Dropout(dropout)
    
    def forward(
        self,
        forward_features: torch.Tensor,
        higher_features: torch.Tensor,
        up_mapping: Dict[int, List[int]]
    ) -> torch.Tensor:
        """
        Args:
            forward_features: Features from forward pass [n_lower, hidden_dim]
            higher_features: Features from higher level [n_higher, hidden_dim]
            up_mapping: Mapping from lower node (index) to list of higher node indices
            
        Returns:
            Updated features combining forward and backward information
            [n_lower, hidden_dim]
        """
        n_lower = forward_features.size(0)
        n_higher = higher_features.size(0)
        device = forward_features.device
        
        # Build attention mask: each lower node attends to its upper nodes
        attn_mask = torch.full((n_lower, n_higher), float('-inf'), device=device)
        
        for v in range(n_lower):
            higher_indices = up_mapping.get(v, [])
            if len(higher_indices) > 0:
                indices = torch.tensor(higher_indices, device=device, dtype=torch.long)
                attn_mask[v, indices] = 0.0
        
        # Handle rows where all entries are -inf (nodes with no upper connections)
        uniform_rows = torch.isinf(attn_mask).all(dim=1)
        attn_mask[uniform_rows] = 0.0
        
        # Cross attention: lower attends to higher
        Q = self.W_Q(forward_features).unsqueeze(0)  # [1, n_lower, hidden_dim]
        K = self.W_K(higher_features).unsqueeze(0)   # [1, n_higher, hidden_dim]
        V = self.W_V(higher_features).unsqueeze(0)   # [1, n_higher, hidden_dim]
        
        attn_out, _ = self.attention(Q, K, V, attn_mask=attn_mask)
        attn_out = attn_out.squeeze(0)  # [n_lower, hidden_dim]
        
        # Backward features with residual
        backward_features = self.norm1(forward_features + self.dropout(attn_out))
        
        # Concatenate forward and backward features as per spec:
        # h_v^{(l-1),new} = LayerNorm(h_v^{(l-1)} + Dropout(MLP([h_v^{(l-1)} || b_v])))
        combined = torch.cat([forward_features, backward_features], dim=-1)
        fused = self.norm2(self.fusion_mlp(combined))
        
        return fused


# ============================================================================
# Sparse Global Attention
# ============================================================================

class SparseGlobalAttention(nn.Module):
    """
    Sparse global attention with adaptive strategy based on graph size.
    
    Adaptive strategy (per spec §3.5.3):
    - n <= 500: full attention
    - 500 < n <= 5000: sparse attention (top-k, k = ceil(min(n, n * log(n))))
    - n > 5000: TNA-only (no global attention)
    """
    
    def __init__(self, hidden_dim: int, num_heads: int = 4, dropout: float = 0.1):
        super().__init__()
        self.hidden_dim = hidden_dim
        self.num_heads = num_heads
        self.head_dim = hidden_dim // num_heads
        
        self.attention = nn.MultiheadAttention(
            hidden_dim, num_heads, dropout=dropout, batch_first=True
        )
        
        self.norm = nn.LayerNorm(hidden_dim)
        self.dropout = nn.Dropout(dropout)
        
        self.ffn = nn.Sequential(
            nn.Linear(hidden_dim, hidden_dim * 4),
            nn.GELU(),
            nn.Dropout(dropout),
            nn.Linear(hidden_dim * 4, hidden_dim),
            nn.Dropout(dropout)
        )
        self.ffn_norm = nn.LayerNorm(hidden_dim)
    
    def _compute_top_k(self, n: int) -> int:
        """
        Compute top-k for sparse attention: k = ceil(min(n, max(2, n * log(n)))).
        
        Args:
            n: Number of nodes
            
        Returns:
            k: Number of top neighbors to attend to
        """
        if n <= 1:
            return 1
        k = int(math.ceil(min(n, max(2, n * math.log(n)))))
        return min(k, n)
    
    def forward(
        self,
        x: torch.Tensor,
        adj_matrix: Optional[torch.Tensor] = None,
        attn_mode: str = 'adaptive'
    ) -> torch.Tensor:
        """
        Args:
            x: Node features [n, hidden_dim]
            adj_matrix: Adjacency matrix [n, n] (optional, used for future enhancement)
            attn_mode: Attention mode ('full', 'sparse', 'adaptive', 'tna_only')
            
        Returns:
            Updated features [n, hidden_dim]
        """
        n = x.size(0)
        
        # Determine attention strategy
        if attn_mode == 'tna_only':
            return x
        elif attn_mode == 'full':
            use_full = True
        elif attn_mode == 'sparse':
            use_full = False
        elif attn_mode == 'adaptive':
            # Adaptive: n <= 500 → full, 500 < n ≤ 5000 → sparse, n > 5000 → skip
            if n <= 500:
                use_full = True
            elif n <= 5000:
                use_full = False
            else:
                return x
        else:
            # Default to full attention
            use_full = True
        
        if use_full:
            # Full attention
            attn_out, _ = self.attention(x.unsqueeze(0), x.unsqueeze(0), x.unsqueeze(0))
            attn_out = attn_out.squeeze(0)
        else:
            # Sparse attention with top-k
            k = self._compute_top_k(n)
            attn_out = self._sparse_attention(x, k)
        
        # Residual and FFN
        x = self.norm(x + self.dropout(attn_out))
        ffn_out = self.ffn(x)
        x = self.ffn_norm(x + ffn_out)
        
        return x
    
    def _sparse_attention(self, x: torch.Tensor, k: int) -> torch.Tensor:
        """
        Memory-efficient sparse attention using top-k selection.
        Uses chunked scatter_reduce to avoid building a dense n×n mask.
        
        Args:
            x: Node features [n, hidden_dim]
            k: Number of top neighbors
            
        Returns:
            Attended features [n, hidden_dim]
        """
        n = x.size(0)
        device = x.device
        head_dim = self.head_dim
        num_heads = self.num_heads
        k_actual = min(k, n)
        
        # Compute multi-head Q, K, V from x
        Q = x.unsqueeze(0)  # [1, n, hidden_dim]
        K = x.unsqueeze(0)
        V = x.unsqueeze(0)
        
        if k_actual >= n:
            # Fallback to full attention when k >= n
            attn_out, _ = self.attention(Q, K, V)
            return attn_out.squeeze(0)
        
        # Compute dot-product scores for top-k selection: [n, n]
        scores = torch.mm(x, x.t()) / math.sqrt(head_dim)
        
        # Select top-k indices per node
        _, topk_indices = torch.topk(scores, k_actual, dim=-1)  # [n, k_actual]
        
        # Gather K and V for only top-k indices to avoid dense n×n mask
        # K_selected: [n, k_actual, hidden_dim]
        K_selected = K[0, topk_indices]   # [n, k_actual, hidden_dim]
        V_selected = V[0, topk_indices]   # [n, k_actual, hidden_dim]
        
        # Compute attention scores efficiently: Q * K^T for selected pairs
        # Q: [1, n, hidden_dim] -> expand to [n, 1, hidden_dim]
        Q_expanded = Q.permute(1, 0, 2)  # [n, 1, hidden_dim]
        
        # Multi-head: reshape all to [n, num_heads, k_actual, head_dim]
        q_h = Q_expanded.view(n, 1, num_heads, head_dim).transpose(1, 2)   # [n, num_heads, 1, head_dim]
        k_h = K_selected.view(n, k_actual, num_heads, head_dim).permute(0, 2, 1, 3)  # [n, num_heads, k_actual, head_dim]
        v_h = V_selected.view(n, k_actual, num_heads, head_dim).permute(0, 2, 1, 3)  # [n, num_heads, k_actual, head_dim]
        
        # Scaled dot-product: [n, num_heads, 1, k_actual]
        attn_scores = torch.matmul(q_h, k_h.transpose(-2, -1)) / math.sqrt(head_dim)
        attn_weights = F.softmax(attn_scores, dim=-1)  # [n, num_heads, 1, k_actual]
        
        # Weighted sum: [n, num_heads, 1, head_dim] -> [n, num_heads, head_dim]
        attn_out = torch.matmul(attn_weights, v_h).squeeze(2)
        # [n, num_heads, head_dim] -> [n, hidden_dim]
        attn_out = attn_out.transpose(1, 0).contiguous().view(n, -1)
        
        return attn_out


# ============================================================================
# CSG-Transformer Layer
# ============================================================================

class CSGTransformerLayer(nn.Module):
    """
    Single CSG-Transformer layer combining TNA, forward/backward cross attention,
    and optional sparse global attention.
    
    Implements one iteration of Algorithm 2 from the specification.
    """
    
    def __init__(
        self,
        hidden_dim: int,
        num_heads: int = 4,
        dropout: float = 0.1,
        edge_feat_dim: Optional[int] = None,
        use_backward: bool = True,
        use_tna: bool = True,
    ):
        super().__init__()
        self.hidden_dim = hidden_dim
        self.use_backward = use_backward
        self.use_tna = use_tna
        
        # TNA Attention
        self.tna = TNAAttention(hidden_dim, num_heads, dropout, edge_feat_dim)
        
        # Forward cross attention
        self.forward_cross_attn = ForwardCrossAttention(hidden_dim, num_heads, dropout)
        
        # Backward cross attention
        self.backward_cross_attn = BackwardCrossAttention(hidden_dim, num_heads, dropout)
        
        # Sparse global attention
        self.global_attn = SparseGlobalAttention(hidden_dim, num_heads, dropout)
    
    def forward(
        self,
        layer_features: List[torch.Tensor],
        layer_neighbor_components: List[Dict],
        layer_mappings: List[Dict[int, List[int]]],
        layer_pe: List[torch.Tensor],
        T: int = 3,
        attn_mode: str = 'adaptive',
        edge_features: Optional[torch.Tensor] = None,
        edge_index: Optional[torch.Tensor] = None,
        adjacency_matrix: Optional[torch.Tensor] = None,
    ) -> List[torch.Tensor]:
        """
        Perform one iteration of CSG-Transformer (Algorithm 2).
        
        Forward phase (bottom-up):
        for l = 0..L-1:
            T_times: TNA on layer l
            ForwardCrossAttn: layer l -> layer l+1
            T_times: TNA on layer l+1
        
        Backward phase (top-down):
        for l = L..1:
            BackwardCrossAttn: layer l -> layer l-1
            T_times: TNA on layer l-1
        
        Global attention applied to layer 0 at the end.
        
        Args:
            layer_features: List of feature tensors for each level [L+1]
            layer_neighbor_components: List of neighbor component dicts [L+1]
            layer_mappings: List of down_mappings (lower -> upper) [L]
            layer_pe: List of positional encodings [L+1]
            T: Number of TNA rounds per stage
            attn_mode: Attention mode for global attention
            edge_features: Optional edge features
            edge_index: Optional edge index
            adjacency_matrix: Optional adjacency matrix
            
        Returns:
            Updated layer_features [L+1]
        """
        L = len(layer_features) - 1  # Number of CSG layers (0 = input graph)
        
        # ============================================================
        # Forward phase (自底向上)
        # ============================================================
        for l in range(L):
            # TNA on current level (T times)
            if self.use_tna:
                for t in range(T):
                    layer_features[l] = self.tna(
                        layer_features[l],
                        layer_neighbor_components[l],
                        edge_features=edge_features,
                        edge_index=edge_index,
                        adjacency_matrix=adjacency_matrix,
                    )
            
            # Forward cross attention to construct upper level initial features
            down_mapping = layer_mappings[l]
            pe_upper = layer_pe[l + 1] if l + 1 < len(layer_pe) else None
            
            # Initialize upper level features from lower (if zero/empty)
            # Uses a learned parameter per upper node via attention-weighted pooling
            upper_current = layer_features[l + 1]
            if torch.equal(upper_current, torch.zeros_like(upper_current)):
                # Use attention-weighted mean of mapped lower features per upper node
                down_mapping = layer_mappings[l]
                upper_init = []
                for u in range(upper_current.size(0)):
                    lower_indices = down_mapping.get(u, [])
                    if lower_indices:
                        lower_feats = layer_features[l][lower_indices]
                        upper_init.append(lower_feats.mean(dim=0))
                    else:
                        upper_init.append(torch.zeros(self.hidden_dim, device=upper_current.device))
                if upper_init:
                    upper_current = torch.stack(upper_init)
                else:
                    upper_current = torch.zeros_like(upper_current)
            
            layer_features[l + 1] = self.forward_cross_attn(
                layer_features[l],
                upper_current,
                down_mapping,
                pe_upper,
            )
            
            # TNA on next level (T times)
            if self.use_tna:
                for t in range(T):
                    layer_features[l + 1] = self.tna(
                        layer_features[l + 1],
                        layer_neighbor_components[l + 1],
                        edge_features=edge_features,
                        edge_index=edge_index,
                        adjacency_matrix=adjacency_matrix,
                    )
        
        # ============================================================
        # Backward phase (自顶向下)
        # ============================================================
        if self.use_backward:
            for l in range(L, 0, -1):
                # Build up_mapping: for each lower node, list of upper nodes it maps to
                down_mapping = layer_mappings[l - 1]
                up_mapping: Dict[int, List[int]] = {}
                for upper_node, lower_nodes in down_mapping.items():
                    if isinstance(upper_node, int) and isinstance(lower_nodes, list):
                        for lower_node in lower_nodes:
                            if lower_node not in up_mapping:
                                up_mapping[lower_node] = []
                            up_mapping[lower_node].append(upper_node)
                
                # Backward cross attention: inject higher-level info into lower level
                layer_features[l - 1] = self.backward_cross_attn(
                    layer_features[l - 1],
                    layer_features[l],
                    up_mapping,
                )
                
                # TNA on current level (T times)
                if self.use_tna:
                    for t in range(T):
                        layer_features[l - 1] = self.tna(
                            layer_features[l - 1],
                            layer_neighbor_components[l - 1],
                            edge_features=edge_features,
                            edge_index=edge_index,
                            adjacency_matrix=adjacency_matrix,
                        )
        
        # ============================================================
        # Optional: Apply global attention to level 0
        # ============================================================
        layer_features[0] = self.global_attn(
            layer_features[0],
            attn_mode=attn_mode,
        )
        
        return layer_features


# ============================================================================
# Main CSG-Transformer Model
# ============================================================================

class CSGTransformer(nn.Module):
    """
    CSG-Transformer Model.
    
    A graph neural network that combines cyclic schematic graph abstraction
    with transformer architecture for graph representation learning.
    
    Key features:
    - Triangulated Neighborhood Attention (TNA) with edge features
    - Forward/Backward cross attention for multi-scale fusion
    - Multi-layer composite positional encoding (LPE + RWSE + Structural)
    - Sparse global attention with adaptive strategy
    - CSG caching for efficient repeated forward passes
    - Edge label support
    """
    
    def __init__(
        self,
        in_dim: int,
        hidden_dim: int = 128,
        out_dim: Optional[int] = None,
        num_heads: int = 4,
        dropout: float = 0.1,
        L: int = 3,
        T: int = 3,
        I: int = 5,
        task: str = 'graph_classification',
        pe_type: str = 'composite',
        attn_mode: str = 'adaptive',
        edge_feat_dim: Optional[int] = None,
        use_backward: bool = True,
        use_tna: bool = True,
        use_pe: bool = True,
    ):
        """
        Args:
            in_dim: Input feature dimension
            hidden_dim: Hidden dimension
            out_dim: Output dimension (defaults to hidden_dim)
            num_heads: Number of attention heads
            dropout: Dropout rate
            L: Number of CSG layers
            T: Number of TNA rounds per layer per iteration
            I: Number of global iterations (Algorithm 2 repeated I times)
            task: 'graph_classification' or 'node_classification'
            pe_type: Positional encoding type ('lpe', 'rwse', 'composite')
            attn_mode: Attention mode ('full', 'sparse', 'adaptive', 'tna_only')
            edge_feat_dim: Edge feature dimension (None if no edge features)
            use_backward: Enable backward cross attention (for ablation)
            use_tna: Enable TNA attention (for ablation)
            use_pe: Enable positional encoding (for ablation)
        """
        super().__init__()
        
        self.in_dim = in_dim
        self.hidden_dim = hidden_dim
        self.out_dim = out_dim or hidden_dim
        self.L = L
        self.T = T
        self.I = I
        self.task = task
        self.pe_type = pe_type
        self.attn_mode = attn_mode
        self.edge_feat_dim = edge_feat_dim
        self.use_backward = use_backward
        self.use_tna = use_tna
        self.use_pe = use_pe
        
        # Input projection: projects raw features to hidden_dim
        self.input_proj = nn.Sequential(
            nn.Linear(in_dim, hidden_dim),
            nn.LayerNorm(hidden_dim),
            nn.Dropout(dropout)
        )
        
        # Positional encoding
        if pe_type == 'composite' and use_pe:
            self.pe_encoder = CompositePositionalEncoding(pe_dim=16, hidden_dim=hidden_dim)
        elif pe_type == 'lpe' and use_pe:
            self.pe_encoder = LaplacianPositionalEncoding(pe_dim=hidden_dim)
        elif pe_type == 'rwse' and use_pe:
            self.pe_encoder = RandomWalkStructuralEncoding(pe_dim=hidden_dim)
        else:
            self.pe_encoder = None
        
        # PE projection for single-type PE (lpe/rwse produce pe_dim dims, not hidden_dim)
        if pe_type in ('lpe', 'rwse') and use_pe:
            self.pe_proj = nn.Linear(hidden_dim if pe_type == 'lpe' else hidden_dim, hidden_dim)
        else:
            self.pe_proj = None
        
        # CSG-Transformer layers: I iterations of Algorithm 2
        self.layers = nn.ModuleList([
            CSGTransformerLayer(
                hidden_dim, num_heads, dropout, edge_feat_dim=edge_feat_dim,
                use_backward=use_backward, use_tna=use_tna,
            )
            for _ in range(I)
        ])
        
        # Classification head
        self.classifier = nn.Sequential(
            nn.Linear(hidden_dim, hidden_dim),
            nn.LayerNorm(hidden_dim),
            nn.GELU(),
            nn.Dropout(dropout),
            nn.Linear(hidden_dim, self.out_dim)
        )
        
        # CSG cache: stores precomputed CSG structures keyed by graph id
        self.csg_cache: Dict[int, SingleGraphCSG] = {}
        
        # Initialize weights
        self._init_weights()
    
    def _init_weights(self):
        """Initialize model weights using Xavier uniform initialization."""
        for m in self.modules():
            if isinstance(m, nn.Linear):
                nn.init.xavier_uniform_(m.weight)
                if m.bias is not None:
                    nn.init.zeros_(m.bias)
            elif isinstance(m, nn.LayerNorm):
                nn.init.ones_(m.weight)
                nn.init.zeros_(m.bias)
    
    def _get_or_build_csg(self, G: nx.Graph) -> SingleGraphCSG:
        """
        Get cached CSG or build and cache it.
        
        Args:
            G: NetworkX graph
            
        Returns:
            SingleGraphCSG instance
        """
        g_id = id(G)
        if g_id not in self.csg_cache:
            self.csg_cache[g_id] = build_multilayer_csg_single(G, L=self.L)
        return self.csg_cache[g_id]
    
    def clear_cache(self):
        """Clear the CSG cache."""
        self.csg_cache.clear()
    
    def _compute_adjacency_matrix(self, G: nx.Graph, device: torch.device) -> torch.Tensor:
        """Compute adjacency matrix from NetworkX graph."""
        return torch.tensor(
            nx.adjacency_matrix(G).todense(),
            dtype=torch.float32,
            device=device,
        )
    
    def _prepare_layer_features(
        self,
        G: nx.Graph,
        csg: SingleGraphCSG,
        node_features: torch.Tensor,
        edge_features: Optional[torch.Tensor],
        device: torch.device,
    ) -> Tuple[
        List[torch.Tensor],
        List[torch.Tensor],
        List[Dict],
        List[Dict[int, List[int]]],
        Optional[torch.Tensor],
        Optional[torch.Tensor],
    ]:
        """
        Prepare all layer features, PE, mappings, and neighbor components.
        
        Returns:
            layer_features, layer_pe, layer_neighbor_components,
            layer_mappings, edge_index_tensor, adj_matrix
        """
        # Build layer graphs list
        layer_graphs = [G]
        for k in range(self.L):
            H_k, _, _ = csg.layers[k]
            layer_graphs.append(H_k)
        
        # Input projection and PE for level 0
        x = self.input_proj(node_features)
        
        # Compute PE for level 0
        adj_0 = self._compute_adjacency_matrix(G, device)
        pe_0 = None
        if self.pe_encoder is not None and self.use_pe:
            pe_0 = self.pe_encoder(adj_0)
            x = x + pe_0
        
        # Initialize layer features
        layer_features = [x]
        for k in range(self.L):
            H_k, _, _ = csg.layers[k]
            n_k = H_k.number_of_nodes()
            # Initialize with zeros (will be populated by ForwardCrossAttn)
            layer_features.append(torch.zeros(n_k, self.hidden_dim, device=device))
        
        layer_node_to_idx = [
            {node: idx for idx, node in enumerate(H.nodes())}
            for H in layer_graphs
        ]
        
        layer_neighbor_components = []
        for k in range(len(layer_graphs)):
            nc_dict: Dict[int, Tuple[Tuple, ...]] = {}
            if k <= self.L:
                nc = get_neighbor_components(csg, k)
                node_to_idx = layer_node_to_idx[k]
                for node, comps in nc.items():
                    node_idx = node_to_idx.get(node)
                    if node_idx is None:
                        continue
                    nc_dict[node_idx] = tuple(
                        tuple(node_to_idx.get(n, n) for n in comp)
                        for comp in comps
                    )
            layer_neighbor_components.append(nc_dict)
        
        layer_mappings = []
        for k in range(self.L):
            mapping = get_mapping(csg, k + 1)
            lower_to_idx = layer_node_to_idx[k]
            upper_to_idx = layer_node_to_idx[k + 1]
            
            down_mapping: Dict[int, List[int]] = {}
            for lower_node, upper_nodes in mapping.items():
                l_idx = lower_to_idx.get(lower_node)
                if l_idx is None:
                    continue
                for upper_node in upper_nodes:
                    u_idx = upper_to_idx.get(upper_node)
                    if u_idx is None:
                        continue
                    if u_idx not in down_mapping:
                        down_mapping[u_idx] = []
                    down_mapping[u_idx].append(l_idx)
            layer_mappings.append(down_mapping)
        
        # Compute PE for all layers
        layer_pe = [pe_0] if pe_0 is not None else []
        for k in range(1, self.L + 1):
            H_k, _, _ = csg.layers[k - 1]
            adj_k = self._compute_adjacency_matrix(H_k, device)
            pe_k = self.pe_encoder(adj_k) if self.pe_encoder is not None else None
            if pe_k is not None:
                layer_pe.append(pe_k)
            else:
                layer_pe.append(torch.zeros(H_k.number_of_nodes(), self.hidden_dim, device=device))
        
        # Prepare edge features for TNA
        edge_index_tensor = None
        if edge_features is not None:
            # Build edge_index from G
            edges = list(G.edges())
            if edges:
                edge_index_tensor = torch.tensor(
                    [[u, v] for u, v in edges] + [[v, u] for u, v in edges],
                    dtype=torch.long, device=device
                ).t().contiguous()
        
        return (
            layer_features, layer_pe, layer_neighbor_components,
            layer_mappings, edge_index_tensor, adj_0,
        )
    
    def forward(
        self,
        G: nx.Graph,
        node_features: torch.Tensor,
        edge_features: Optional[torch.Tensor] = None,
    ) -> Tuple[torch.Tensor, torch.Tensor]:
        """
        Forward pass of CSG-Transformer.
        
        Args:
            G: NetworkX graph
            node_features: Node features [n, in_dim]
            edge_features: Optional edge features [m, edge_dim]
            
        Returns:
            If graph_classification:
                - graph_embedding: [hidden_dim]
                - logits: [out_dim]
            If node_classification:
                - node_embeddings: [n, hidden_dim]
                - logits: [n, out_dim]
        """
        device = node_features.device
        
        # Get or build CSG (with caching)
        csg = self._get_or_build_csg(G)
        
        # Prepare all layer structures
        (
            layer_features, layer_pe, layer_neighbor_components,
            layer_mappings, edge_index_tensor, adj_0,
        ) = self._prepare_layer_features(
            G, csg, node_features, edge_features, device
        )
        
        # Apply I iterations of CSG-Transformer layers (Algorithm 2)
        for layer in self.layers:
            layer_features = layer(
                layer_features,
                layer_neighbor_components,
                layer_mappings,
                layer_pe,
                T=self.T,
                attn_mode=self.attn_mode,
                edge_features=edge_features,
                edge_index=edge_index_tensor,
                adjacency_matrix=adj_0,
            )
        
        # Get final features for level 0 (input graph nodes)
        node_embeddings = layer_features[0]
        
        # Pooling for graph-level tasks
        if self.task == 'graph_classification':
            graph_embedding = node_embeddings.mean(dim=0)
            logits = self.classifier(graph_embedding)
            return graph_embedding, logits
        else:
            logits = self.classifier(node_embeddings)
            return node_embeddings, logits


# ============================================================================
# Unified Interface
# ============================================================================

def csg_transformer_unified(
    G1: nx.Graph,
    G2: nx.Graph,
    node_features1: np.ndarray,
    node_features2: np.ndarray,
    edge_features1: Optional[np.ndarray] = None,
    edge_features2: Optional[np.ndarray] = None,
    L: int = 3,
    T: int = 3,
    I: int = 5,
    hidden_dim: int = 128,
    num_heads: int = 4,
    dropout: float = 0.1,
    task: str = 'graph_classification',
    pe_type: str = 'composite',
    attn_mode: str = 'adaptive',
) -> Tuple[torch.Tensor, torch.Tensor]:
    """
    Unified interface for CSG-Transformer.
    
    Computes graph embeddings for a pair of graphs using CSG-Transformer.
    
    Args:
        G1: First NetworkX graph
        G2: Second NetworkX graph
        node_features1: Node features for G1 [n1, in_dim]
        node_features2: Node features for G2 [n2, in_dim]
        edge_features1: Optional edge features for G1 [m1, edge_dim]
        edge_features2: Optional edge features for G2 [m2, edge_dim]
        L: Number of CSG layers
        T: Number of TNA rounds per layer
        I: Number of global iterations
        hidden_dim: Hidden dimension
        num_heads: Number of attention heads
        dropout: Dropout rate
        task: 'graph_classification' or 'node_classification'
        pe_type: Positional encoding type
        attn_mode: Attention mode
        
    Returns:
        graph_embedding1: [hidden_dim]
        graph_embedding2: [hidden_dim]
    """
    # Determine input dimension
    in_dim = node_features1.shape[1]
    
    # Determine edge feature dimension
    edge_feat_dim = None
    if edge_features1 is not None:
        edge_feat_dim = edge_features1.shape[1]
    
    # Create model
    model = CSGTransformer(
        in_dim=in_dim,
        hidden_dim=hidden_dim,
        num_heads=num_heads,
        dropout=dropout,
        L=L,
        T=T,
        I=I,
        task=task,
        pe_type=pe_type,
        attn_mode=attn_mode,
        edge_feat_dim=edge_feat_dim,
    )
    
    # Convert to tensors
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    model = model.to(device)
    
    x1 = torch.tensor(node_features1, dtype=torch.float32, device=device)
    x2 = torch.tensor(node_features2, dtype=torch.float32, device=device)
    
    ef1 = torch.tensor(edge_features1, dtype=torch.float32, device=device) if edge_features1 is not None else None
    ef2 = torch.tensor(edge_features2, dtype=torch.float32, device=device) if edge_features2 is not None else None
    
    # Forward pass
    model.eval()
    with torch.no_grad():
        emb1, _ = model(G1, x1, ef1)
        emb2, _ = model(G2, x2, ef2)
    
    return emb1, emb2


# ============================================================================
# Utility Functions
# ============================================================================

def compute_graph_similarity(
    emb1: torch.Tensor,
    emb2: torch.Tensor,
    metric: str = 'cosine'
) -> float:
    """
    Compute similarity between two graph embeddings.
    
    Args:
        emb1: First graph embedding [hidden_dim]
        emb2: Second graph embedding [hidden_dim]
        metric: Similarity metric ('cosine', 'euclidean', 'dot')
        
    Returns:
        Similarity score
    """
    if metric == 'cosine':
        return F.cosine_similarity(emb1.unsqueeze(0), emb2.unsqueeze(0)).item()
    elif metric == 'euclidean':
        return -torch.norm(emb1 - emb2, p=2).item()
    elif metric == 'dot':
        return torch.dot(emb1, emb2).item()
    else:
        raise ValueError(f"Unknown metric: {metric}")


def build_model_from_config(config: Dict[str, Any]) -> CSGTransformer:
    """
    Build CSG-Transformer model from config dictionary.
    
    Supports both `L` and `num_layers` as the CSG layer count parameter
    (num_layers is an alias for L, used in config YAML files).
    
    Args:
        config: Dictionary with model hyperparameters
        
    Returns:
        CSGTransformer model instance
    """
    # Support num_layers as alias for L
    L_value = config.get('L', config.get('num_layers', 3))
    
    return CSGTransformer(
        in_dim=config['in_dim'],
        hidden_dim=config.get('hidden_dim', 128),
        out_dim=config.get('out_dim', None),
        num_heads=config.get('num_heads', 4),
        dropout=config.get('dropout', 0.1),
        L=L_value,
        T=config.get('T', 3),
        I=config.get('I', 5),
        task=config.get('task', 'graph_classification'),
        pe_type=config.get('pe_type', 'composite'),
        attn_mode=config.get('attn_mode', 'adaptive'),
        edge_feat_dim=config.get('edge_feat_dim', None),
        use_backward=config.get('use_backward', True),
        use_tna=config.get('use_tna', True),
        use_pe=config.get('use_pe', True),
    )
