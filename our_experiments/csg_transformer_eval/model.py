"""
CSG-Transformer Model - Thin wrapper around the core implementation.

This module re-exports the CSG-Transformer model from cyclic_schema.csg_transformer
for use in the evaluation pipeline. It provides the same interface with additional
convenience for loading from checkpoints.

Usage:
    from our_experiments.csg_transformer_eval.model import CSGTransformer
    model = CSGTransformer(in_dim=7, hidden_dim=128, out_dim=2)
    emb, logits = model(G, features)
"""

import sys
import os
from typing import Any, Dict, Optional

import torch
import torch.nn as nn

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))

from cyclic_schema.csg_transformer import CSGTransformer as _CSGTransformer
from cyclic_schema.csg_transformer import build_model_from_config as _build_model_from_config


class CSGTransformer(_CSGTransformer):
    """
    CSG-Transformer Model (evaluation wrapper).
    
    Inherits from cyclic_schema.csg_transformer.CSGTransformer.
    Adds checkpoint saving/loading for the evaluation pipeline.
    """
    
    def save_checkpoint(self, filepath: str) -> None:
        """Save model state dict and config to checkpoint.
        
        Compatible with load_checkpoint() - can be loaded back directly.
        """
        config = {
            'in_dim': self.in_dim,
            'hidden_dim': self.hidden_dim,
            'out_dim': self.out_dim,
            'num_heads': getattr(self, 'num_heads', 4),
            'dropout': getattr(self, '_dropout_rate', 0.1),
            'L': self.L,
            'T': self.T,
            'I': self.I,
            'task': self.task,
            'pe_type': self.pe_type,
            'attn_mode': self.attn_mode,
            'edge_feat_dim': getattr(self, 'edge_feat_dim', None),
            'use_backward': self.use_backward,
            'use_tna': self.use_tna,
            'use_pe': self.use_pe,
            'use_virtual_node': getattr(self, 'use_virtual_node', False),
            'use_rel_bias': getattr(self, 'use_rel_bias', False),
            'pooling': getattr(self, 'pooling', 'mean'),
        }
        torch.save({
            'config': config,
            'state_dict': self.state_dict(),
        }, filepath)
    
    @classmethod
    def load_checkpoint(cls, filepath: str, map_location: str = 'cpu') -> 'CSGTransformer':
        """Load model from checkpoint."""
        checkpoint = torch.load(filepath, map_location=map_location, weights_only=False)
        config = checkpoint['config']
        model = cls(
            in_dim=config['in_dim'],
            hidden_dim=config.get('hidden_dim', 128),
            out_dim=config.get('out_dim', None),
            num_heads=config.get('num_heads', 4),
            dropout=config.get('dropout', 0.1),
            L=config.get('L', 3),
            T=config.get('T', 3),
            I=config.get('I', 5),
            task=config.get('task', 'graph_classification'),
            pe_type=config.get('pe_type', 'composite'),
            attn_mode=config.get('attn_mode', 'adaptive'),
            edge_feat_dim=config.get('edge_feat_dim', None),
            use_backward=config.get('use_backward', True),
            use_tna=config.get('use_tna', True),
            use_pe=config.get('use_pe', True),
            use_virtual_node=config.get('use_virtual_node', False),
            use_rel_bias=config.get('use_rel_bias', False),
            pooling=config.get('pooling', 'mean'),
        )
        model.load_state_dict(checkpoint['state_dict'])
        return model


def build_model(config: Dict[str, Any]) -> CSGTransformer:
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
        num_heads=config.get('num_heads', 4),
        dropout=config.get('dropout', 0.1),
        L=config.get('L', 3),
        T=config.get('T', 3),
        I=config.get('I', 5),
        task=config.get('task', 'graph_classification'),
        pe_type=config.get('pe_type', 'composite'),
        attn_mode=config.get('attn_mode', 'adaptive'),
        edge_feat_dim=config.get('edge_feat_dim', None),
        use_virtual_node=config.get('use_virtual_node', False),
        use_rel_bias=config.get('use_rel_bias', False),
        pooling=config.get('pooling', 'mean'),
    )
