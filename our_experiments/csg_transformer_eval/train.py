"""
Training pipeline for CSG-Transformer model.

Includes:
- Loss function with L2 regularization
- Model parameter initialization
- Validation set verification
- Early stopping mechanism
"""

import torch
import torch.nn as nn
import torch.nn.functional as F
from torch.optim import AdamW
from torch.optim.lr_scheduler import CosineAnnealingLR
import numpy as np
from sklearn.metrics import accuracy_score, roc_auc_score, f1_score


class EarlyStopping:
    """Early stopping to prevent overfitting."""
    
    def __init__(self, patience=20, min_delta=1e-4, mode='max'):
        self.patience = patience
        self.min_delta = min_delta
        self.mode = mode
        self.counter = 0
        self.best_score = None
        self.early_stop = False
    
    def __call__(self, score):
        if self.best_score is None:
            self.best_score = score
        elif self.mode == 'max':
            if score < self.best_score + self.min_delta:
                self.counter += 1
                if self.counter >= self.patience:
                    self.early_stop = True
            else:
                self.best_score = score
                self.counter = 0
        else:
            if score > self.best_score - self.min_delta:
                self.counter += 1
                if self.counter >= self.patience:
                    self.early_stop = True
            else:
                self.best_score = score
                self.counter = 0


class L2Regularization:
    """L2 regularization for model parameters."""
    
    def __init__(self, model, lambda_reg=1e-4):
        self.model = model
        self.lambda_reg = lambda_reg
    
    def __call__(self):
        l2_reg = torch.tensor(0.0, device=next(self.model.parameters()).device)
        for param in self.model.parameters():
            l2_reg += torch.norm(param, p=2)
        return self.lambda_reg * l2_reg


def init_weights(model):
    """Initialize model weights using Xavier uniform initialization."""
    for m in model.modules():
        if isinstance(m, nn.Linear):
            nn.init.xavier_uniform_(m.weight)
            if m.bias is not None:
                nn.init.zeros_(m.bias)
        elif isinstance(m, nn.LayerNorm):
            nn.init.ones_(m.weight)
            nn.init.zeros_(m.bias)


def compute_loss(outputs, targets, model, lambda_reg=1e-4, task='graph_classification'):
    """
    Compute loss with L2 regularization.
    
    Args:
        outputs: Model outputs (logits)
        targets: Ground truth labels
        model: Model for L2 regularization
        lambda_reg: L2 regularization coefficient
        task: 'graph_classification' or 'node_classification'
    
    Returns:
        Total loss (cross entropy + L2 regularization)
    """
    if task == 'graph_classification':
        ce_loss = F.cross_entropy(outputs, targets)
    else:
        ce_loss = F.cross_entropy(outputs, targets)
    
    l2_reg = L2Regularization(model, lambda_reg)()
    
    return ce_loss + l2_reg


def evaluate(model, data_loader, device, task='graph_classification'):
    """
    Evaluate model on validation/test set.
    
    Args:
        model: Trained model
        data_loader: Validation/test data loader
        device: Device to use
        task: 'graph_classification' or 'node_classification'
    
    Returns:
        Dictionary with metrics: loss, accuracy, auc, f1
    """
    model.eval()
    total_loss = 0
    all_preds = []
    all_targets = []
    all_probs = []
    
    with torch.no_grad():
        for batch in data_loader:
            batch = batch.to(device)
            _, logits = model(batch)
            
            loss = F.cross_entropy(logits, batch.y)
            total_loss += loss.item()
            
            probs = F.softmax(logits, dim=-1)
            preds = logits.argmax(dim=-1)
            
            all_preds.extend(preds.cpu().numpy())
            all_targets.extend(batch.y.cpu().numpy())
            all_probs.extend(probs.cpu().numpy())
    
    all_preds = np.array(all_preds)
    all_targets = np.array(all_targets)
    all_probs = np.array(all_probs)
    
    accuracy = accuracy_score(all_targets, all_preds)
    
    if task == 'graph_classification' and len(np.unique(all_targets)) == 2:
        auc = roc_auc_score(all_targets, all_probs[:, 1])
    else:
        auc = roc_auc_score(all_targets, all_probs, multi_class='ovr', average='macro')
    
    f1 = f1_score(all_targets, all_preds, average='macro')
    
    return {
        'loss': total_loss / len(data_loader),
        'accuracy': accuracy,
        'auc': auc,
        'f1': f1
    }


def train_one_epoch(model, data_loader, optimizer, device, lambda_reg=1e-4, task='graph_classification'):
    """
    Train model for one epoch.
    
    Args:
        model: Model to train
        data_loader: Training data loader
        optimizer: Optimizer
        device: Device to use
        lambda_reg: L2 regularization coefficient
        task: 'graph_classification' or 'node_classification'
    
    Returns:
        Average training loss
    """
    model.train()
    total_loss = 0
    
    for batch in data_loader:
        batch = batch.to(device)
        optimizer.zero_grad()
        
        _, logits = model(batch)
        loss = compute_loss(logits, batch.y, model, lambda_reg, task)
        
        loss.backward()
        optimizer.step()
        
        total_loss += loss.item()
    
    return total_loss / len(data_loader)


def train(model, train_loader, val_loader, config, device):
    """
    Full training pipeline with validation and early stopping.
    
    Args:
        model: Model to train
        train_loader: Training data loader
        val_loader: Validation data loader
        config: Training configuration dictionary
        device: Device to use
    
    Returns:
        Dictionary with training history and best model
    """
    optimizer = AdamW(
        model.parameters(),
        lr=config.get('lr', 1e-3),
        weight_decay=config.get('weight_decay', 1e-4)
    )
    
    scheduler = CosineAnnealingLR(
        optimizer,
        T_max=config.get('epochs', 200),
        eta_min=config.get('min_lr', 1e-6)
    )
    
    early_stopping = EarlyStopping(
        patience=config.get('patience', 30),
        min_delta=config.get('min_delta', 1e-4),
        mode='max'
    )
    
    best_val_acc = 0
    best_model_state = None
    history = {
        'train_loss': [],
        'val_loss': [],
        'val_acc': [],
        'val_auc': [],
        'val_f1': []
    }
    
    for epoch in range(config.get('epochs', 200)):
        train_loss = train_one_epoch(
            model, train_loader, optimizer, device,
            config.get('lambda_reg', 1e-4),
            config.get('task', 'graph_classification')
        )
        
        val_metrics = evaluate(model, val_loader, device, config.get('task', 'graph_classification'))
        
        scheduler.step()
        
        history['train_loss'].append(train_loss)
        history['val_loss'].append(val_metrics['loss'])
        history['val_acc'].append(val_metrics['accuracy'])
        history['val_auc'].append(val_metrics['auc'])
        history['val_f1'].append(val_metrics['f1'])
        
        print(f"Epoch {epoch+1}/{config.get('epochs', 200)}: "
              f"Train Loss: {train_loss:.4f}, "
              f"Val Acc: {val_metrics['accuracy']:.4f}, "
              f"Val AUC: {val_metrics['auc']:.4f}, "
              f"Val F1: {val_metrics['f1']:.4f}")
        
        if val_metrics['accuracy'] > best_val_acc:
            best_val_acc = val_metrics['accuracy']
            best_model_state = model.state_dict().copy()
        
        early_stopping(val_metrics['accuracy'])
        if early_stopping.early_stop:
            print(f"Early stopping at epoch {epoch+1}")
            break
    
    if best_model_state is not None:
        model.load_state_dict(best_model_state)
    
    return {
        'model': model,
        'history': history,
        'best_val_acc': best_val_acc
    }
