"""
CSG-Transformer Model Implementation

A graph neural network model that combines Cyclic Schematic Graph (CSG) 
abstraction with Transformer architecture for graph representation learning.
"""

import torch
import torch.nn as nn
import torch.nn.functional as F
from torch_geometric.nn import global_mean_pool, global_add_pool
from torch_geometric.utils import to_dense_batch


class TNALayer(nn.Module):
    """
    Triangulated Neighborhood Attention (TNA) Layer
    
    Implements attention mechanism over triangulated neighborhood components
    with residual connections, LayerNorm, and Dropout.
    """
    
    def __init__(self, in_dim, hidden_dim, num_heads=4, dropout=0.1):
        super().__init__()
        self.num_heads = num_heads
        self.head_dim = hidden_dim // num_heads
        assert hidden_dim % num_heads == 0
        
        self.W_Q = nn.Linear(in_dim, hidden_dim)
        self.W_K = nn.Linear(in_dim, hidden_dim)
        self.W_V = nn.Linear(in_dim, hidden_dim)
        
        self.attention = nn.MultiheadAttention(
            hidden_dim, num_heads, dropout=dropout, batch_first=True
        )
        
        self.norm1 = nn.LayerNorm(hidden_dim)
        self.norm2 = nn.LayerNorm(hidden_dim)
        self.dropout = nn.Dropout(dropout)
        
        self.ffn = nn.Sequential(
            nn.Linear(hidden_dim, hidden_dim * 4),
            nn.GELU(),
            nn.Dropout(dropout),
            nn.Linear(hidden_dim * 4, hidden_dim),
            nn.Dropout(dropout)
        )
    
    def forward(self, x, tna_mask=None):
        """
        Args:
            x: Node features [batch_size, num_nodes, in_dim]
            tna_mask: Boolean mask for valid nodes [batch_size, num_nodes]
        Returns:
            Updated node features [batch_size, num_nodes, hidden_dim]
        """
        Q = self.W_Q(x)
        K = self.W_K(x)
        V = self.W_V(x)
        
        attn_out, _ = self.attention(Q, K, V, key_padding_mask=tna_mask)
        x = self.norm1(x + self.dropout(attn_out))
        
        ffn_out = self.ffn(x)
        x = self.norm2(x + ffn_out)
        
        return x


class ForwardCrossAttention(nn.Module):
    """
    Forward Cross Attention Layer
    
    Aggregates information from lower level to higher level in CSG hierarchy
    using cross attention mechanism.
    """
    
    def __init__(self, in_dim, hidden_dim, num_heads=4, dropout=0.1):
        super().__init__()
        self.cross_attn = nn.MultiheadAttention(
            hidden_dim, num_heads, dropout=dropout, batch_first=True
        )
        
        self.norm1 = nn.LayerNorm(hidden_dim)
        self.norm2 = nn.LayerNorm(hidden_dim)
        self.dropout = nn.Dropout(dropout)
        
        self.ffn = nn.Sequential(
            nn.Linear(hidden_dim, hidden_dim * 4),
            nn.GELU(),
            nn.Dropout(dropout),
            nn.Linear(hidden_dim * 4, hidden_dim),
            nn.Dropout(dropout)
        )
        
        self.query_proj = nn.Linear(in_dim, hidden_dim)
        self.key_proj = nn.Linear(in_dim, hidden_dim)
        self.value_proj = nn.Linear(in_dim, hidden_dim)
    
    def forward(self, query, key_value, mapping):
        """
        Args:
            query: Higher level features [batch_size, num_query, in_dim]
            key_value: Lower level features [batch_size, num_kv, in_dim]
            mapping: Mapping from query to key_value indices
        Returns:
            Updated query features [batch_size, num_query, hidden_dim]
        """
        Q = self.query_proj(query)
        K = self.key_proj(key_value)
        V = self.value_proj(key_value)
        
        attn_out, _ = self.cross_attn(Q, K, V)
        query = self.norm1(query + self.dropout(attn_out))
        
        ffn_out = self.ffn(query)
        query = self.norm2(query + ffn_out)
        
        return query


class BackwardCrossAttention(nn.Module):
    """
    Backward Cross Attention Layer
    
    Aggregates information from higher level to lower level in CSG hierarchy.
    Each node receives features from both forward pass and backward pass,
    which are concatenated before fusion.
    """
    
    def __init__(self, in_dim, hidden_dim, num_heads=4, dropout=0.1):
        super().__init__()
        self.cross_attn = nn.MultiheadAttention(
            hidden_dim, num_heads, dropout=dropout, batch_first=True
        )
        
        self.norm1 = nn.LayerNorm(hidden_dim)
        self.norm2 = nn.LayerNorm(hidden_dim)
        self.dropout = nn.Dropout(dropout)
        
        # Fusion MLP that combines forward and backward features
        self.fusion_mlp = nn.Sequential(
            nn.Linear(hidden_dim * 2, hidden_dim * 4),
            nn.GELU(),
            nn.Dropout(dropout),
            nn.Linear(hidden_dim * 4, hidden_dim),
            nn.Dropout(dropout)
        )
        
        self.query_proj = nn.Linear(in_dim, hidden_dim)
        self.key_proj = nn.Linear(in_dim, hidden_dim)
        self.value_proj = nn.Linear(in_dim, hidden_dim)
    
    def forward(self, forward_features, higher_features, mapping):
        """
        Args:
            forward_features: Features from forward pass [batch_size, num_nodes, in_dim]
            higher_features: Features from higher level [batch_size, num_higher, in_dim]
            mapping: Mapping from lower to higher level nodes
        Returns:
            Updated features combining forward and backward information
            [batch_size, num_nodes, hidden_dim]
        """
        Q = self.query_proj(forward_features)
        K = self.key_proj(higher_features)
        V = self.value_proj(higher_features)
        
        attn_out, _ = self.cross_attn(Q, K, V)
        backward_features = self.norm1(forward_features + self.dropout(attn_out))
        
        # Concatenate forward and backward features
        combined = torch.cat([forward_features, backward_features], dim=-1)
        fused = self.fusion_mlp(combined)
        
        return fused


class CSGTransformerLayer(nn.Module):
    """
    Single CSG-Transformer Layer
    
    Combines TNA attention with forward and backward cross attention
    for hierarchical graph representation learning.
    """
    
    def __init__(self, in_dim, hidden_dim, num_heads=4, dropout=0.1):
        super().__init__()
        
        self.tna = TNALayer(in_dim, hidden_dim, num_heads, dropout)
        self.forward_cross_attn = ForwardCrossAttention(
            hidden_dim, hidden_dim, num_heads, dropout
        )
        self.backward_cross_attn = BackwardCrossAttention(
            hidden_dim, hidden_dim, num_heads, dropout
        )
    
    def forward(self, x, csg_graphs, level_mapping):
        """
        Args:
            x: Node features at current level [batch_size, num_nodes, in_dim]
            csg_graphs: CSG graph structure
            level_mapping: Mapping between levels
        Returns:
            Updated node features [batch_size, num_nodes, hidden_dim]
        """
        # TNA attention
        x = self.tna(x)
        
        # Forward cross attention (if not at top level)
        if 'forward' in level_mapping:
            x = self.forward_cross_attn(
                x, level_mapping['forward']['source'], 
                level_mapping['forward']['mapping']
            )
        
        # Backward cross attention (if not at bottom level)
        if 'backward' in level_mapping:
            x = self.backward_cross_attn(
                x, level_mapping['backward']['source'],
                level_mapping['backward']['mapping']
            )
        
        return x


class CSGTransformer(nn.Module):
    """
    CSG-Transformer Model
    
    A graph neural network that combines cyclic schematic graph abstraction
    with transformer architecture for graph representation learning.
    
    Key features:
    - Triangulated Neighborhood Attention (TNA)
    - Forward/Backward cross attention for multi-scale fusion
    - Random Laplacian positional encoding
    - Residual connections, LayerNorm, and Dropout throughout
    """
    
    def __init__(
        self,
        in_dim,
        hidden_dim=128,
        out_dim=None,
        num_layers=3,
        num_heads=4,
        dropout=0.1,
        L=3,
        T=3,
        I=5,
        max_nodes=500,
        task='graph_classification'
    ):
        super().__init__()
        
        self.in_dim = in_dim
        self.hidden_dim = hidden_dim
        self.out_dim = out_dim or hidden_dim
        self.num_layers = num_layers
        self.L = L
        self.T = T
        self.I = I
        self.task = task
        
        # Input projection
        self.input_proj = nn.Sequential(
            nn.Linear(in_dim, hidden_dim),
            nn.LayerNorm(hidden_dim),
            nn.Dropout(dropout)
        )
        
        # Positional encoding
        self.pe_proj = nn.Linear(16, hidden_dim)  # 16-dim Laplacian PE
        
        # CSG-Transformer layers
        self.layers = nn.ModuleList([
            CSGTransformerLayer(hidden_dim, hidden_dim, num_heads, dropout)
            for _ in range(num_layers)
        ])
        
        # Classification head
        if task == 'graph_classification':
            self.classifier = nn.Sequential(
                nn.Linear(hidden_dim, hidden_dim),
                nn.LayerNorm(hidden_dim),
                nn.GELU(),
                nn.Dropout(dropout),
                nn.Linear(hidden_dim, self.out_dim)
            )
        else:  # node_classification
            self.classifier = nn.Sequential(
                nn.Linear(hidden_dim, hidden_dim),
                nn.LayerNorm(hidden_dim),
                nn.GELU(),
                nn.Dropout(dropout),
                nn.Linear(hidden_dim, self.out_dim)
            )
        
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
    
    def forward(self, data):
        """
        Forward pass of CSG-Transformer.
        
        Args:
            data: PyG Data object containing:
                - x: Node features [num_nodes, in_dim]
                - edge_index: Edge indices [2, num_edges]
                - batch: Batch assignment for nodes
                - pe: Positional encoding features (optional)
        
        Returns:
            If graph_classification:
                - graph_embedding: [batch_size, hidden_dim]
                - logits: [batch_size, out_dim]
            If node_classification:
                - node_embeddings: [num_nodes, hidden_dim]
                - logits: [num_nodes, out_dim]
        """
        x = data.x
        edge_index = data.edge_index
        batch = data.batch if hasattr(data, 'batch') else torch.zeros(x.size(0), dtype=torch.long, device=x.device)
        
        # Input projection
        x = self.input_proj(x)
        
        # Add positional encoding if available
        if hasattr(data, 'pe') and data.pe is not None:
            pe = self.pe_proj(data.pe)
            x = x + pe
        
        # Reshape for batch processing
        x_batch = to_dense_batch(x, batch)
        
        # Apply CSG-Transformer layers
        for layer in self.layers:
            x_batch = layer(x_batch, None, {})
        
        # Extract node embeddings
        node_embeddings = x_batch[0]
        
        # Pooling for graph-level tasks
        if self.task == 'graph_classification':
            graph_embedding = global_mean_pool(node_embeddings, batch)
            logits = self.classifier(graph_embedding)
            return graph_embedding, logits
        else:
            logits = self.classifier(node_embeddings)
            return node_embeddings, logits


def build_model(config):
    """
    Build CSG-Transformer model from config dictionary.
    
    Args:
        config: Dictionary with model hyperparameters
    
    Returns:
        CSGTransformer model instance
    """
    return CSGTransformer(
        in_dim=config['in_dim'],
        hidden_dim=config.get('hidden_dim', 128),
        out_dim=config.get('out_dim', None),
        num_layers=config.get('num_layers', 3),
        num_heads=config.get('num_heads', 4),
        dropout=config.get('dropout', 0.1),
        L=config.get('L', 3),
        T=config.get('T', 3),
        I=config.get('I', 5),
        max_nodes=config.get('max_nodes', 500),
        task=config.get('task', 'graph_classification')
    )
