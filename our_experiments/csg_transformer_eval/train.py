"""
Training pipeline for CSG-Transformer model.

Includes:
- Loss function with L2 regularization
- Model parameter initialization
- Validation set verification
- Early stopping mechanism
- 10-fold cross-validation for graph classification
- Results tracking and CSV output
"""

import os
import sys
import json
from typing import Any, Dict, List, Optional, Tuple

import networkx as nx
import numpy as np
import pandas as pd
import torch
import torch.nn as nn
import torch.nn.functional as F
from sklearn.metrics import (
    accuracy_score, balanced_accuracy_score,
    roc_auc_score, f1_score, precision_score, recall_score,
)
from sklearn.model_selection import StratifiedKFold

from torch.amp.grad_scaler import GradScaler
from torch.amp.autocast_mode import autocast

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))

from cyclic_schema.csg_transformer import CSGTransformer, build_model_from_config


class EarlyStopping:
    """Early stopping to prevent overfitting."""
    
    def __init__(self, patience: int = 20, min_delta: float = 1e-4, mode: str = 'max'):
        self.patience = patience
        self.min_delta = min_delta
        self.mode = mode
        self.counter = 0
        self.best_score = None
        self.early_stop = False
    
    def __call__(self, score: float) -> bool:
        if self.best_score is None:
            self.best_score = score
            return False
        
        if self.mode == 'max':
            improved = score > self.best_score + self.min_delta
        else:
            improved = score < self.best_score - self.min_delta
        
        if improved:
            self.best_score = score
            self.counter = 0
            return False
        else:
            self.counter += 1
            if self.counter >= self.patience:
                self.early_stop = True
                return True
            return False


class L2Regularization:
    """L2 regularization for model parameters.

    Computes lambda * ||theta||_2^2  (squared L2 norm), matching the
    specification formula: L = CE + lambda * ||theta||_2^2.

    Note: When the optimizer already has decoupled weight_decay (e.g. AdamW),
    this explicit term is skipped (lambda_reg=0) to avoid double-counting.
    """
    
    def __init__(self, model: nn.Module, lambda_reg: float = 1e-4):
        self.model = model
        self.lambda_reg = lambda_reg
    
    def __call__(self) -> torch.Tensor:
        if self.lambda_reg == 0.0:
            return torch.tensor(0.0, device=next(self.model.parameters()).device)
        device = next(self.model.parameters()).device
        l2_reg = torch.tensor(0.0, device=device)
        for param in self.model.parameters():
            l2_reg = l2_reg + param.pow(2).sum()
        return self.lambda_reg * l2_reg


def init_weights(model: nn.Module) -> None:
    """Initialize model weights using Xavier uniform initialization."""
    for m in model.modules():
        if isinstance(m, nn.Linear):
            nn.init.xavier_uniform_(m.weight)
            if m.bias is not None:
                nn.init.zeros_(m.bias)
        elif isinstance(m, nn.LayerNorm):
            nn.init.ones_(m.weight)
            nn.init.zeros_(m.bias)


class EMA(nn.Module):
    """Exponential Moving Average of model parameters for self-knowledge distillation."""
    def __init__(self, model: nn.Module, decay: float = 0.999):
        super().__init__()
        self.decay = decay
        self.shadow: Dict[str, torch.Tensor] = {}
        self.backup: Dict[str, torch.Tensor] = {}
        self.register(model)

    def register(self, model: nn.Module) -> None:
        for name, param in model.named_parameters():
            if param.requires_grad:
                self.shadow[name] = param.data.clone()

    def update(self, model: nn.Module) -> None:
        for name, param in model.named_parameters():
            if param.requires_grad and name in self.shadow:
                new_average = self.decay * self.shadow[name] + (1.0 - self.decay) * param.data
                self.shadow[name] = new_average.clone()

    def apply_shadow(self, model: nn.Module) -> None:
        for name, param in model.named_parameters():
            if param.requires_grad and name in self.shadow:
                self.backup[name] = param.data.clone()
                param.data = self.shadow[name]

    def restore(self, model: nn.Module) -> None:
        for name, param in model.named_parameters():
            if name in self.backup:
                param.data = self.backup[name]
        self.backup = {}


def _drop_edge(G: nx.Graph, drop_rate: float, rng: Optional[np.random.RandomState] = None) -> nx.Graph:
    """Randomly drop edges from a graph, ensuring connectivity is preserved."""
    if drop_rate <= 0.0:
        return G
    G_sub = G.copy()
    edges = list(G_sub.edges())
    if not edges:
        return G_sub
    if rng is None:
        rng = np.random.RandomState()
    n_keep = max(1, int(len(edges) * (1.0 - drop_rate)))
    keep = rng.choice(len(edges), size=n_keep, replace=False)
    drop_edges = [edges[i] for i in range(len(edges)) if i not in keep]
    G_sub.remove_edges_from(drop_edges)
    # Reconnect disconnected components
    if not nx.is_connected(G_sub) and G_sub.number_of_nodes() > 1:
        components = list(nx.connected_components(G_sub))
        for i in range(len(components) - 1):
            u = rng.choice(list(components[i]))
            v = rng.choice(list(components[i + 1]))
            G_sub.add_edge(u, v)
    return G_sub


def compute_loss(
    outputs: torch.Tensor,
    targets: torch.Tensor,
    model: nn.Module,
    lambda_reg: float = 0.0,
    teacher_logits: Optional[torch.Tensor] = None,
    kd_temperature: float = 4.0,
    kd_weight: float = 0.0,
) -> torch.Tensor:
    """Compute cross-entropy loss with L2 regularization and optional KD."""
    ce_loss = F.cross_entropy(outputs, targets)
    l2_reg = L2Regularization(model, lambda_reg)()
    total_loss = ce_loss + l2_reg
    if teacher_logits is not None and kd_weight > 0.0:
        student_log_sm = F.log_softmax(outputs / kd_temperature, dim=-1)
        teacher_sm = F.softmax(teacher_logits / kd_temperature, dim=-1)
        kd_loss = F.kl_div(student_log_sm, teacher_sm, reduction='batchmean')
        kd_loss = kd_loss * (kd_temperature ** 2)
        total_loss = total_loss + kd_weight * kd_loss
    return total_loss


def train_one_epoch(
    model: CSGTransformer,
    data_list: List[Tuple],
    optimizer: torch.optim.Optimizer,
    device: torch.device,
    lambda_reg: float = 0.0,
    grad_clip: Optional[float] = None,
    use_amp: bool = False,
    scaler: Optional[GradScaler] = None,
    drop_edge_rate: float = 0.0,
    ema: Optional[EMA] = None,
    kd_temperature: float = 4.0,
    kd_weight: float = 0.0,
    epoch: int = 0,
    ema_start_epoch: int = 0,
    rng: Optional[np.random.RandomState] = None,
) -> float:
    """Train model for one epoch with optional AMP, gradient clipping, DropEdge, and Self-KD."""
    model.train()
    total_loss = 0.0
    
    for G, features, label in data_list:
        features = features.to(device)
        label = label.to(device)
        
        G_forward = _drop_edge(G, drop_edge_rate, rng=rng) if drop_edge_rate > 0.0 else G
        
        optimizer.zero_grad()
        
        teacher_logits = None
        if ema is not None and epoch >= ema_start_epoch:
            ema.apply_shadow(model)
            with torch.no_grad():
                _, t_logits = model(G_forward, features)
                if t_logits.dim() == 1:
                    t_logits = t_logits.unsqueeze(0)
                teacher_logits = t_logits.detach()
            ema.restore(model)
        
        if use_amp and scaler is not None:
            with autocast('cuda'):
                _, logits = model(G_forward, features)
                if logits.dim() == 1:
                    logits = logits.unsqueeze(0)
                if label.dim() == 0:
                    label = label.unsqueeze(0)
                loss = compute_loss(logits, label, model, lambda_reg,
                                    teacher_logits=teacher_logits,
                                    kd_temperature=kd_temperature,
                                    kd_weight=kd_weight)
            scaler.scale(loss).backward()
            if grad_clip is not None:
                scaler.unscale_(optimizer)
                torch.nn.utils.clip_grad_norm_(model.parameters(), grad_clip)
            scaler.step(optimizer)
            scaler.update()
        else:
            _, logits = model(G_forward, features)
            if logits.dim() == 1:
                logits = logits.unsqueeze(0)
            if label.dim() == 0:
                label = label.unsqueeze(0)
            loss = compute_loss(logits, label, model, lambda_reg,
                                teacher_logits=teacher_logits,
                                kd_temperature=kd_temperature,
                                kd_weight=kd_weight)
            loss.backward()
            if grad_clip is not None:
                torch.nn.utils.clip_grad_norm_(model.parameters(), grad_clip)
            optimizer.step()
        
        if ema is not None and epoch >= ema_start_epoch:
            ema.update(model)
        
        total_loss += loss.item()
    
    if drop_edge_rate > 0.0 and hasattr(model, 'clear_cache'):
        model.clear_cache()
    
    return total_loss / len(data_list)


def evaluate(
    model: CSGTransformer,
    data_list: List[Tuple],
    device: torch.device,
    task: str = 'graph_classification',
) -> Dict[str, float]:
    """Evaluate model on validation or test set."""
    model.eval()
    all_preds = []
    all_targets = []
    all_probs = []
    
    with torch.no_grad():
        for G, features, label in data_list:
            features = features.to(device)
            label = label.to(device)
            
            _, logits = model(G, features)
            
            if logits.dim() == 1:
                logits = logits.unsqueeze(0)
            if label.dim() == 0:
                label = label.unsqueeze(0)
            
            probs = F.softmax(logits, dim=-1)
            preds = logits.argmax(dim=-1)
            
            all_preds.extend(preds.cpu().numpy())
            all_targets.extend(label.cpu().numpy())
            all_probs.extend(probs.cpu().numpy())
    
    all_preds = np.array(all_preds)
    all_targets = np.array(all_targets)
    all_probs = np.array(all_probs)
    
    num_classes = len(np.unique(all_targets))
    
    metrics = {
        'accuracy': accuracy_score(all_targets, all_preds),
        'balanced_accuracy': balanced_accuracy_score(all_targets, all_preds),
        'f1_macro': f1_score(all_targets, all_preds, average='macro'),
        'f1_weighted': f1_score(all_targets, all_preds, average='weighted'),
        'precision_macro': precision_score(all_targets, all_preds, average='macro', zero_division=0),
        'recall_macro': recall_score(all_targets, all_preds, average='macro', zero_division=0),
    }
    
    if num_classes == 2:
        metrics['auc_roc'] = roc_auc_score(all_targets, all_probs[:, 1])
    elif num_classes > 2:
        try:
            metrics['auc_roc'] = roc_auc_score(
                all_targets, all_probs, multi_class='ovr', average='macro'
            )
        except Exception:
            metrics['auc_roc'] = 0.0
    
    return metrics


def _get_lr_warmup_factor(epoch: int, warmup_epochs: int) -> float:
    """Compute warmup factor for linear LR warmup."""
    if epoch >= warmup_epochs:
        return 1.0
    return (epoch + 1) / warmup_epochs


def train_model(
    model: CSGTransformer,
    train_data: List[Tuple],
    val_data: List[Tuple],
    config: Dict,
    device: torch.device,
    verbose: bool = True,
) -> Dict:
    """
    Full training pipeline with validation, early stopping, gradient clipping, and AMP.
    
    Args:
        model: CSG-Transformer model
        train_data: List of (G, features, label)
        val_data: List of (G, features, label)
        config: Training configuration dict
        device: Device
        verbose: Print progress
        
    Returns:
        Dict with 'model', 'history', 'best_val_acc', 'best_model_state'
    """
    optimizer = torch.optim.AdamW(
        model.parameters(),
        lr=config.get('lr', 1e-3),
        weight_decay=config.get('weight_decay', 1e-4),
    )
    
    num_epochs = config.get('epochs', 200)
    warmup_epochs = min(config.get('warmup_epochs', 10), num_epochs)
    
    scheduler = torch.optim.lr_scheduler.CosineAnnealingLR(
        optimizer,
        T_max=max(1, num_epochs - warmup_epochs),
        eta_min=config.get('min_lr', 1e-6),
    )
    
    early_stopping = EarlyStopping(
        patience=config.get('patience', 30),
        min_delta=config.get('min_delta', 1e-4),
        mode='max',
    )
    
    use_amp = config.get('use_amp', False)
    grad_clip = config.get('grad_clip', 1.0)  # Default gradient clipping
    scaler = GradScaler('cuda') if use_amp else None
    
    # DropEdge
    drop_edge_rate = config.get('drop_edge_rate', 0.0)
    
    # Self-KD (EMA teacher)
    use_ema = config.get('use_ema', False)
    ema_decay = config.get('ema_decay', 0.999)
    ema_start_epoch = config.get('ema_start_epoch', 0)
    kd_temperature = config.get('kd_temperature', 4.0)
    kd_weight = config.get('kd_weight', 0.0)
    ema_model = EMA(model, decay=ema_decay) if use_ema else None
    _rng = np.random.RandomState(config.get('seed', 42))
    
    best_val_acc = 0.0
    best_model_state = None
    history = {
        'train_loss': [],
        'val_acc': [],
        'val_auc': [],
        'val_f1': [],
    }
    
    for epoch in range(num_epochs):
        # Linear LR warmup
        if epoch < warmup_epochs:
            factor = _get_lr_warmup_factor(epoch, warmup_epochs)
            for param_group in optimizer.param_groups:
                param_group['lr'] = config.get('lr', 1e-3) * factor
        else:
            scheduler.step()
        
        train_loss = train_one_epoch(
            model, train_data, optimizer, device,
            lambda_reg=config.get('lambda_reg', 1e-4),
            grad_clip=grad_clip,
            use_amp=use_amp,
            scaler=scaler,
            drop_edge_rate=drop_edge_rate,
            ema=ema_model,
            kd_temperature=kd_temperature,
            kd_weight=kd_weight,
            epoch=epoch,
            ema_start_epoch=ema_start_epoch,
            rng=_rng,
        )
        
        val_metrics = evaluate(
            model, val_data, device,
            config.get('task', 'graph_classification'),
        )
        
        history['train_loss'].append(train_loss)
        history['val_acc'].append(val_metrics['accuracy'])
        history['val_auc'].append(val_metrics.get('auc_roc', 0.0))
        history['val_f1'].append(val_metrics['f1_macro'])
        
        if verbose and (epoch + 1) % 10 == 0:
            print(f"Epoch {epoch+1}/{config.get('epochs', 200)}: "
                  f"Loss={train_loss:.4f}, Val ACC={val_metrics['accuracy']:.4f}, "
                  f"AUC={val_metrics.get('auc_roc', 0):.4f}")
        
        if val_metrics['accuracy'] > best_val_acc:
            best_val_acc = val_metrics['accuracy']
            best_model_state = model.state_dict().copy()
        
        if early_stopping(val_metrics['accuracy']):
            if verbose:
                print(f"Early stopping at epoch {epoch+1}")
            break
    
    if best_model_state is not None:
        model.load_state_dict(best_model_state)
    
    return {
        'model': model,
        'history': history,
        'best_val_acc': best_val_acc,
        'best_model_state': best_model_state,
    }


def cross_validate(
    dataset: List[Tuple],
    config: Dict,
    device: torch.device,
    n_folds: int = 10,
    verbose: bool = True,
    checkpoint_dir: Optional[str] = None,
) -> pd.DataFrame:
    """
    Perform N-fold cross-validation for graph classification.
    
    Args:
        dataset: List of (G, features, label) tuples
        config: Training configuration
        device: Device to use
        n_folds: Number of folds (default: 10)
        verbose: Whether to print progress
        checkpoint_dir: Optional dir to save best model per fold
        
    Returns:
        DataFrame with results for each fold
    """
    # Use stratified split to preserve label distribution across folds
    labels = [d[2].item() if torch.is_tensor(d[2]) else d[2] for d in dataset]
    skf = StratifiedKFold(n_splits=n_folds, shuffle=True, random_state=config.get('seed', 42))
    
    if checkpoint_dir:
        os.makedirs(checkpoint_dir, exist_ok=True)
    
    fold_results = []
    config = config.copy()
    
    for fold, (train_idx, val_idx) in enumerate(skf.split(dataset, labels)):
        if verbose:
            print(f"\n{'='*50}")
            print(f"Fold {fold + 1}/{n_folds}")
            print(f"{'='*50}")
        
        train_data = [dataset[i] for i in train_idx]
        val_data = [dataset[i] for i in val_idx]
        
        in_dim = train_data[0][1].shape[1]
        config['in_dim'] = in_dim
        
        model = build_model_from_config(config)
        model = model.to(device)
        
        result = train_model(
            model, train_data, val_data, config, device, verbose=verbose
        )
        
        # Save checkpoint to disk if requested
        if checkpoint_dir and result.get('best_model_state') is not None:
            ckpt_path = os.path.join(checkpoint_dir, f"fold_{fold+1}_best.pt")
            model_for_save = build_model_from_config(config)
            model_for_save.load_state_dict(result['best_model_state'])
            model_for_save.save_checkpoint(ckpt_path)
            if verbose:
                print(f"Checkpoint saved: {ckpt_path}")
        
        test_metrics = evaluate(
            model, val_data, device,
            config.get('task', 'graph_classification'),
        )
        
        fold_results.append({
            'fold': fold + 1,
            'train_loss': result['history']['train_loss'][-1] if result['history']['train_loss'] else 0,
            'val_acc': result['best_val_acc'],
            'test_acc': test_metrics['accuracy'],
            'test_balanced_accuracy': test_metrics.get('balanced_accuracy', 0.0),
            'test_auc': test_metrics.get('auc_roc', 0.0),
            'test_f1': test_metrics['f1_macro'],
            'test_precision_macro': test_metrics.get('precision_macro', 0.0),
            'test_recall_macro': test_metrics.get('recall_macro', 0.0),
            'epochs': len(result['history']['train_loss']),
        })
        
        if verbose:
            print(f"Fold {fold + 1} - Val Acc: {result['best_val_acc']:.4f}, "
                  f"Test Acc: {test_metrics['accuracy']:.4f}")
    
    df = pd.DataFrame(fold_results)
    
    if verbose:
        print(f"\n{'='*50}")
        print(f"Cross-validation Summary ({n_folds}-fold)")
        print(f"{'='*50}")
        print(f"Val Acc: {df['val_acc'].mean():.4f} ± {df['val_acc'].std():.4f}")
        print(f"Test Acc: {df['test_acc'].mean():.4f} ± {df['test_acc'].std():.4f}")
        print(f"Test Balanced Acc: {df['test_balanced_accuracy'].mean():.4f} ± {df['test_balanced_accuracy'].std():.4f}")
        print(f"Test AUC: {df['test_auc'].mean():.4f} ± {df['test_auc'].std():.4f}")
        print(f"Test F1: {df['test_f1'].mean():.4f} ± {df['test_f1'].std():.4f}")
        print(f"Test Precision: {df['test_precision_macro'].mean():.4f} ± {df['test_precision_macro'].std():.4f}")
        print(f"Test Recall: {df['test_recall_macro'].mean():.4f} ± {df['test_recall_macro'].std():.4f}")
    
    return df


def save_results_to_csv(
    results: pd.DataFrame,
    save_dir: str,
    experiment_name: str,
    config: Optional[Dict] = None,
) -> str:
    """
    Save results to CSV file with experiment metadata.
    
    Args:
        results: Results DataFrame
        save_dir: Directory to save to
        experiment_name: Name of the experiment
        config: Optional configuration dict to include
        
    Returns:
        Path to saved CSV file
    """
    os.makedirs(save_dir, exist_ok=True)
    
    result_df = results.copy()
    
    if config:
        for key, value in config.items():
            if key not in result_df.columns:
                result_df[key] = str(value)
    
    filepath = os.path.join(save_dir, f"{experiment_name}_results.csv")
    result_df.to_csv(filepath, index=False)
    
    # Compute statistics with the expanded set of metrics
    has_auc = 'test_auc' in results.columns
    has_f1 = 'test_f1' in results.columns
    has_balanced_acc = 'test_balanced_accuracy' in results.columns
    has_precision = 'test_precision_macro' in results.columns
    has_recall = 'test_recall_macro' in results.columns
    
    summary_dict = {
        'experiment': experiment_name,
        'mean_val_acc': results['val_acc'].mean(),
        'std_val_acc': results['val_acc'].std(),
        'mean_test_acc': results['test_acc'].mean(),
        'std_test_acc': results['test_acc'].std(),
    }
    if has_auc:
        summary_dict['mean_test_auc'] = results['test_auc'].mean()
        summary_dict['std_test_auc'] = results['test_auc'].std()
    if has_f1:
        summary_dict['mean_test_f1'] = results['test_f1'].mean()
        summary_dict['std_test_f1'] = results['test_f1'].std()
    if has_balanced_acc:
        summary_dict['mean_test_balanced_accuracy'] = results['test_balanced_accuracy'].mean()
        summary_dict['std_test_balanced_accuracy'] = results['test_balanced_accuracy'].std()
    if has_precision:
        summary_dict['mean_test_precision_macro'] = results['test_precision_macro'].mean()
        summary_dict['std_test_precision_macro'] = results['test_precision_macro'].std()
    if has_recall:
        summary_dict['mean_test_recall_macro'] = results['test_recall_macro'].mean()
        summary_dict['std_test_recall_macro'] = results['test_recall_macro'].std()
    
    summary_path = os.path.join(save_dir, f"{experiment_name}_summary.csv")
    summary = pd.DataFrame([summary_dict])
    summary.to_csv(summary_path, index=False)
    
    print(f"Results saved to {filepath}")
    print(f"Summary saved to {summary_path}")
    return filepath
