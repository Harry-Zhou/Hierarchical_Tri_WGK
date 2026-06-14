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
        """Save model state dict and config to checkpoint."""
        config = {
            'in_dim': self.in_dim,
            'hidden_dim': self.hidden_dim,
            'out_dim': self.out_dim,
            'num_heads': self.layers[0].tna.num_heads if self.layers else 4,
            'dropout': 0.1,  # Not stored, use default
            'L': self.L,
            'T': self.T,
            'I': self.I,
            'task': self.task,
            'pe_type': self.pe_type,
            'attn_mode': self.attn_mode,
        }
        torch.save({
            'config': config,
            'state_dict': self.state_dict(),
        }, filepath)
    
    @classmethod
    def load_checkpoint(cls, filepath: str, map_location: str = 'cpu') -> 'CSGTransformer':
        """Load model from checkpoint."""
        checkpoint = torch.load(filepath, map_location=map_location)
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
    )
