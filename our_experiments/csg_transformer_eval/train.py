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

import numpy as np
import pandas as pd
import torch
import torch.nn as nn
import torch.nn.functional as F
from sklearn.metrics import accuracy_score, roc_auc_score, f1_score
from sklearn.model_selection import StratifiedKFold

from torch.cuda.amp import GradScaler, autocast

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
    """L2 regularization for model parameters."""
    
    def __init__(self, model: nn.Module, lambda_reg: float = 1e-4):
        self.model = model
        self.lambda_reg = lambda_reg
    
    def __call__(self) -> torch.Tensor:
        device = next(self.model.parameters()).device
        l2_reg = torch.tensor(0.0, device=device)
        for param in self.model.parameters():
            l2_reg += torch.norm(param, p=2)
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


def compute_loss(
    outputs: torch.Tensor,
    targets: torch.Tensor,
    model: nn.Module,
    lambda_reg: float = 1e-4,
) -> torch.Tensor:
    """Compute cross-entropy loss with L2 regularization."""
    ce_loss = F.cross_entropy(outputs, targets)
    l2_reg = L2Regularization(model, lambda_reg)()
    return ce_loss + l2_reg


def train_one_epoch(
    model: CSGTransformer,
    data_list: List[Tuple],
    optimizer: torch.optim.Optimizer,
    device: torch.device,
    lambda_reg: float = 1e-4,
    grad_clip: Optional[float] = None,
    use_amp: bool = False,
    scaler: Optional[GradScaler] = None,
) -> float:
    """Train model for one epoch with optional AMP and gradient clipping."""
    model.train()
    total_loss = 0.0
    
    for G, features, label in data_list:
        features = features.to(device)
        label = label.to(device)
        
        optimizer.zero_grad()
        
        if use_amp and scaler is not None:
            with autocast():
                _, logits = model(G, features)
                if logits.dim() == 1:
                    logits = logits.unsqueeze(0)
                if label.dim() == 0:
                    label = label.unsqueeze(0)
                loss = compute_loss(logits, label, model, lambda_reg)
            scaler.scale(loss).backward()
            if grad_clip is not None:
                scaler.unscale_(optimizer)
                torch.nn.utils.clip_grad_norm_(model.parameters(), grad_clip)
            scaler.step(optimizer)
            scaler.update()
        else:
            _, logits = model(G, features)
            if logits.dim() == 1:
                logits = logits.unsqueeze(0)
            if label.dim() == 0:
                label = label.unsqueeze(0)
            loss = compute_loss(logits, label, model, lambda_reg)
            loss.backward()
            if grad_clip is not None:
                torch.nn.utils.clip_grad_norm_(model.parameters(), grad_clip)
            optimizer.step()
        
        total_loss += loss.item()
    
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
    
    metrics = {
        'accuracy': accuracy_score(all_targets, all_preds),
        'f1_macro': f1_score(all_targets, all_preds, average='macro'),
        'f1_weighted': f1_score(all_targets, all_preds, average='weighted'),
    }
    
    unique_classes = len(np.unique(all_targets))
    if unique_classes == 2:
        metrics['auc_roc'] = roc_auc_score(all_targets, all_probs[:, 1])
    elif unique_classes > 2:
        try:
            metrics['auc_roc'] = roc_auc_score(
                all_targets, all_probs, multi_class='ovr', average='macro'
            )
        except Exception:
            metrics['auc_roc'] = 0.0
    
    return metrics


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
    
    scheduler = torch.optim.lr_scheduler.CosineAnnealingLR(
        optimizer,
        T_max=config.get('epochs', 200),
        eta_min=config.get('min_lr', 1e-6),
    )
    
    early_stopping = EarlyStopping(
        patience=config.get('patience', 30),
        min_delta=config.get('min_delta', 1e-4),
        mode='max',
    )
    
    use_amp = config.get('use_amp', False)
    grad_clip = config.get('grad_clip', None)
    scaler = GradScaler() if use_amp else None
    
    best_val_acc = 0.0
    best_model_state = None
    history = {
        'train_loss': [],
        'val_acc': [],
        'val_auc': [],
        'val_f1': [],
    }
    
    for epoch in range(config.get('epochs', 200)):
        train_loss = train_one_epoch(
            model, train_data, optimizer, device,
            lambda_reg=config.get('lambda_reg', 1e-4),
            grad_clip=grad_clip,
            use_amp=use_amp,
            scaler=scaler,
        )
        
        val_metrics = evaluate(
            model, val_data, device,
            config.get('task', 'graph_classification'),
        )
        
        scheduler.step()
        
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
            torch.save({
                'fold': fold + 1,
                'model_state_dict': result['best_model_state'],
                'config': config,
                'val_acc': result['best_val_acc'],
                'history': result['history'],
            }, ckpt_path)
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
            'test_auc': test_metrics.get('auc_roc', 0.0),
            'test_f1': test_metrics['f1_macro'],
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
        print(f"Test AUC: {df['test_auc'].mean():.4f} ± {df['test_auc'].std():.4f}")
        print(f"Test F1: {df['test_f1'].mean():.4f} ± {df['test_f1'].std():.4f}")
    
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
    
    summary_path = os.path.join(save_dir, f"{experiment_name}_summary.csv")
    summary = pd.DataFrame([{
        'experiment': experiment_name,
        'mean_val_acc': results['val_acc'].mean(),
        'std_val_acc': results['val_acc'].std(),
        'mean_test_acc': results['test_acc'].mean(),
        'std_test_acc': results['test_acc'].std(),
        'mean_test_auc': results['test_auc'].mean(),
        'std_test_auc': results['test_auc'].std(),
        'mean_test_f1': results['test_f1'].mean(),
        'std_test_f1': results['test_f1'].std(),
    }])
    summary.to_csv(summary_path, index=False)
    
    print(f"Results saved to {filepath}")
    print(f"Summary saved to {summary_path}")
    return filepath
