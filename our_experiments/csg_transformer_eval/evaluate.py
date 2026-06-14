"""
Evaluation script for CSG-Transformer model.

Provides evaluation utilities for graph and node classification tasks.
"""

import torch
import numpy as np
from sklearn.metrics import accuracy_score, roc_auc_score, f1_score, precision_score, recall_score
import json
import os


def evaluate_model(model, data_loader, device, task='graph_classification'):
    """
    Evaluate model on test set.
    
    Args:
        model: Trained model
        data_loader: Test data loader
        device: Device to use
        task: 'graph_classification' or 'node_classification'
    
    Returns:
        Dictionary with evaluation metrics
    """
    model.eval()
    all_preds = []
    all_targets = []
    all_probs = []
    
    with torch.no_grad():
        for batch in data_loader:
            batch = batch.to(device)
            _, logits = model(batch)
            
            probs = torch.softmax(logits, dim=-1)
            preds = logits.argmax(dim=-1)
            
            all_preds.extend(preds.cpu().numpy())
            all_targets.extend(batch.y.cpu().numpy())
            all_probs.extend(probs.cpu().numpy())
    
    all_preds = np.array(all_preds)
    all_targets = np.array(all_targets)
    all_probs = np.array(all_probs)
    
    metrics = {
        'accuracy': accuracy_score(all_targets, all_preds),
        'f1_macro': f1_score(all_targets, all_preds, average='macro'),
        'precision_macro': precision_score(all_targets, all_preds, average='macro'),
        'recall_macro': recall_score(all_targets, all_preds, average='macro'),
    }
    
    if task == 'graph_classification' and len(np.unique(all_targets)) == 2:
        metrics['auc_roc'] = roc_auc_score(all_targets, all_probs[:, 1])
    else:
        metrics['auc_roc'] = roc_auc_score(
            all_targets, all_probs, multi_class='ovr', average='macro'
        )
    
    return metrics


def evaluate_node_classification(model, data, device, mask_name='test_mask'):
    """
    Evaluate model on node classification task.
    
    Args:
        model: Trained model
        data: PyG Data object
        device: Device to use
        mask_name: Mask to use for evaluation ('test_mask', 'val_mask', 'train_mask')
    
    Returns:
        Dictionary with evaluation metrics
    """
    model.eval()
    data = data.to(device)
    
    with torch.no_grad():
        _, logits = model(data)
    
    mask = getattr(data, mask_name)
    logits = logits[mask]
    targets = data.y[mask]
    
    probs = torch.softmax(logits, dim=-1)
    preds = logits.argmax(dim=-1)
    
    all_preds = preds.cpu().numpy()
    all_targets = targets.cpu().numpy()
    all_probs = probs.cpu().numpy()
    
    metrics = {
        'accuracy': accuracy_score(all_targets, all_preds),
        'f1_macro': f1_score(all_targets, all_preds, average='macro'),
        'precision_macro': precision_score(all_targets, all_preds, average='macro'),
        'recall_macro': recall_score(all_targets, all_preds, average='macro'),
    }
    
    if len(np.unique(all_targets)) == 2:
        metrics['auc_roc'] = roc_auc_score(all_targets, all_probs[:, 1])
    else:
        metrics['auc_roc'] = roc_auc_score(
            all_targets, all_probs, multi_class='ovr', average='macro'
        )
    
    return metrics


def save_results(results, save_dir, experiment_name):
    """
    Save evaluation results to JSON file.
    
    Args:
        results: Dictionary with evaluation results
        save_dir: Directory to save results
        experiment_name: Name of the experiment
    """
    os.makedirs(save_dir, exist_ok=True)
    
    filename = f"{experiment_name}_results.json"
    filepath = os.path.join(save_dir, filename)
    
    with open(filepath, 'w') as f:
        json.dump(results, f, indent=2)
    
    print(f"Results saved to {filepath}")


def load_results(save_dir, experiment_name):
    """
    Load evaluation results from JSON file.
    
    Args:
        save_dir: Directory containing results
        experiment_name: Name of the experiment
    
    Returns:
        Dictionary with evaluation results
    """
    filename = f"{experiment_name}_results.json"
    filepath = os.path.join(save_dir, filename)
    
    with open(filepath, 'r') as f:
        results = json.load(f)
    
    return results


def print_metrics(metrics, task='graph_classification'):
    """
    Print evaluation metrics in formatted manner.
    
    Args:
        metrics: Dictionary with evaluation metrics
        task: 'graph_classification' or 'node_classification'
    """
    print(f"\n{'='*50}")
    print(f"{task.replace('_', ' ').title()} Results")
    print(f"{'='*50}")
    print(f"Accuracy:      {metrics['accuracy']:.4f}")
    print(f"AUC-ROC:       {metrics['auc_roc']:.4f}")
    print(f"F1 (macro):    {metrics['f1_macro']:.4f}")
    print(f"Precision:     {metrics['precision_macro']:.4f}")
    print(f"Recall:        {metrics['recall_macro']:.4f}")
    print(f"{'='*50}\n")
