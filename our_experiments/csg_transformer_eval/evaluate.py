"""
Main evaluation script for CSG-Transformer model.

Provides:
- Graph classification evaluation with 10-fold cross-validation
- Node classification evaluation
- Results aggregation and CSV output
- Ablation study support matching the specification
"""

import argparse
import json
import os
import sys
from datetime import datetime
from typing import Any, Dict, List, Optional, Tuple

import networkx as nx
import numpy as np
import pandas as pd
import scipy.sparse as sp
import torch
import torch.nn.functional as F
import yaml

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))

# Default output directory: our_experiments/csg_transformer_eval/outputs/
_EVAL_DIR = os.path.dirname(os.path.abspath(__file__))
_DEFAULT_OUTPUT_DIR = os.path.join(_EVAL_DIR, 'outputs')

def _build_sparse_adj_norm(G: nx.Graph, device: torch.device, dtype: torch.dtype) -> torch.Tensor:
    """Build symmetrically normalized sparse adjacency tensor."""
    A = nx.to_scipy_sparse_array(G, nodelist=sorted(G.nodes()), format='csr')
    D = sp.diags(np.array(A.sum(axis=1)).flatten() ** -0.5, format='csr')
    A_norm = (D @ A @ D).tocoo()
    indices = torch.from_numpy(np.vstack([A_norm.row, A_norm.col])).long().to(device)
    values = torch.from_numpy(A_norm.data).to(device, dtype=dtype)
    return torch.sparse_coo_tensor(indices, values, A_norm.shape, device=device)


def _correct_predictions(
    G: nx.Graph,
    soft_labels: torch.Tensor,
    labels: torch.Tensor,
    train_mask: torch.Tensor,
    num_iterations: int = 50,
    alpha: float = 0.5,
) -> torch.Tensor:
    """Correct stage: spread residual errors from labeled to unlabeled nodes."""
    A_tensor = _build_sparse_adj_norm(G, soft_labels.device, soft_labels.dtype)
    
    error = torch.zeros_like(soft_labels)
    error[train_mask] = labels[train_mask] - soft_labels[train_mask]
    
    E = error.clone()
    for _ in range(num_iterations):
        E = (1 - alpha) * torch.sparse.mm(A_tensor, E) + alpha * error
    
    return soft_labels + E


def _smooth_predictions(
    G: nx.Graph,
    soft_labels: torch.Tensor,
    train_mask: torch.Tensor,
    num_iterations: int = 50,
    alpha: float = 0.5,
) -> torch.Tensor:
    """Smooth stage: spread predictions across the graph."""
    A_tensor = _build_sparse_adj_norm(G, soft_labels.device, soft_labels.dtype)
    
    Z = soft_labels.clone()
    Z_train = soft_labels[train_mask].clone()
    for _ in range(num_iterations):
        Z = (1 - alpha) * torch.sparse.mm(A_tensor, Z) + alpha * soft_labels
        Z[train_mask] = Z_train
    
    return Z


def apply_correct_and_smooth(
    G: nx.Graph,
    logits: torch.Tensor,
    labels: torch.Tensor,
    train_mask: torch.Tensor,
    correct_iter: int = 50,
    correct_alpha: float = 0.5,
    smooth_iter: int = 50,
    smooth_alpha: float = 0.5,
) -> torch.Tensor:
    """Apply Correct & Smooth post-processing to node classification logits."""
    soft_labels = F.softmax(logits, dim=-1)
    
    # One-hot encode labels for training nodes
    num_classes = soft_labels.shape[1]
    labels_onehot = torch.zeros_like(soft_labels)
    labels_onehot.scatter_(1, labels.unsqueeze(1), 1.0)
    
    # Stage 1: Correct
    corrected = _correct_predictions(
        G, soft_labels, labels_onehot, train_mask,
        num_iterations=correct_iter, alpha=correct_alpha,
    )
    
    # Stage 2: Smooth
    smoothed = _smooth_predictions(
        G, corrected, train_mask,
        num_iterations=smooth_iter, alpha=smooth_alpha,
    )
    
    return F.softmax(smoothed, dim=-1)


from cyclic_schema.csg_transformer import CSGTransformer, build_model_from_config
from our_experiments.csg_transformer_eval.train import (
    cross_validate, save_results_to_csv, evaluate, train_model,
)
from our_experiments.csg_transformer_eval.dataset import (
    load_graph_classification_dataset,
    load_node_classification_dataset,
    get_dataset_stats,
    GRAPH_CLASSIFICATION_DATASETS,
    NODE_CLASSIFICATION_DATASETS,
)
from our_experiments.csg_transformer_eval.utils import set_seed, get_device
from our_experiments.csg_transformer_eval.baselines import (
    compare_with_baselines, summarize_comparisons, get_baseline_table,
)


def prepare_graph_classification_data(
    dataset_name: str,
    root: str = './data',
    use_rewiring: bool = False,
    rewiring_alpha: float = 0.15,
) -> Tuple[List[Tuple], Dict]:
    """Prepare data for graph classification."""
    data_list, stats = load_graph_classification_dataset(
        dataset_name, root, use_rewiring=use_rewiring, rewiring_alpha=rewiring_alpha,
    )
    return data_list, stats


def run_graph_classification_experiment(
    dataset_name: str,
    config: Dict,
    device: torch.device,
    n_folds: int = 10,
    seed: int = 42,
    save_dir: str = '',
) -> pd.DataFrame:
    """Run graph classification experiment with cross-validation."""
    if not save_dir:
        save_dir = os.path.join(_DEFAULT_OUTPUT_DIR, 'graph_classification')
    set_seed(seed)
    
    print(f"\n{'='*60}")
    print(f"Graph Classification: {dataset_name}")
    print(f"{'='*60}")
    
    data_list, stats = prepare_graph_classification_data(
        dataset_name,
        use_rewiring=config.get('use_rewiring', False),
        rewiring_alpha=config.get('rewiring_alpha', 0.15),
    )
    
    config['in_dim'] = stats['in_dim']
    config['out_dim'] = stats['num_classes']
    config['task'] = 'graph_classification'
    
    print(f"Graphs: {stats['num_graphs']}, Classes: {stats['num_classes']}, "
          f"In dim: {stats['in_dim']}, Avg nodes: {stats['avg_nodes']:.1f}")
    
    checkpoint_dir = os.path.join(save_dir, "checkpoints", dataset_name)
    results_df = cross_validate(
        data_list, config, device, n_folds=n_folds, verbose=True,
        checkpoint_dir=checkpoint_dir,
    )
    
    results_df['dataset'] = dataset_name
    results_df['experiment'] = 'graph_classification'
    results_df['model_config'] = str({
        'L': config.get('L', 3), 'T': config.get('T', 3), 'I': config.get('I', 5),
        'hidden_dim': config.get('hidden_dim', 128),
    })
    results_df['timestamp'] = datetime.now().isoformat()
    
    os.makedirs(save_dir, exist_ok=True)
    save_results_to_csv(
        results_df, save_dir,
        f"{dataset_name}_graph_classification",
        config=config,
    )
    
    # Baseline comparison
    try:
        baseline_dir = os.path.join(save_dir, "baselines")
        compare_with_baselines(results_df, task="graph_classification", save_dir=baseline_dir)
    except Exception as e:
        print(f"  Baseline comparison skipped ({e})")
    
    return results_df


def run_node_classification_experiment(
    dataset_name: str,
    config: Dict,
    device: torch.device,
    seed: int = 42,
    save_dir: str = '',
    num_runs: int = 5,
) -> pd.DataFrame:
    """Run node classification experiment with multiple random runs."""
    if not save_dir:
        save_dir = os.path.join(_DEFAULT_OUTPUT_DIR, 'node_classification')
    set_seed(seed)
    
    print(f"\n{'='*60}")
    print(f"Node Classification: {dataset_name}")
    print(f"{'='*60}")
    
    G, features, labels, train_mask, val_mask, test_mask = \
        load_node_classification_dataset(
            dataset_name,
            use_rewiring=config.get('use_rewiring', False),
            rewiring_alpha=config.get('rewiring_alpha', 0.15),
        )
    
    in_dim = features.shape[1]
    num_classes = len(labels.unique())
    
    config['in_dim'] = in_dim
    config['out_dim'] = num_classes
    config['task'] = 'node_classification'
    
    print(f"Nodes: {G.number_of_nodes()}, Classes: {num_classes}, "
          f"In dim: {in_dim}")
    
    all_results = []
    
    for run in range(num_runs):
        set_seed(seed + run)
        
        model = build_model_from_config(config)
        model = model.to(device)
        
        optimizer = torch.optim.AdamW(
            model.parameters(),
            lr=config.get('lr', 0.01),
            weight_decay=config.get('weight_decay', 0.0005),
        )
        
        best_val_acc = 0.0
        best_state = None
        patience = config.get('patience', 50)
        patience_counter = 0
        
        for epoch in range(config.get('epochs', 200)):
            model.train()
            optimizer.zero_grad()
            
            _, logits = model(G, features)
            loss = F.cross_entropy(
                logits[train_mask], labels[train_mask]
            )
            loss.backward()
            optimizer.step()
            
            model.eval()
            with torch.no_grad():
                _, logits = model(G, features)
                val_preds = logits[val_mask].argmax(dim=1)
                val_acc = (val_preds == labels[val_mask]).float().mean().item()
            
            if val_acc > best_val_acc:
                best_val_acc = val_acc
                best_state = model.state_dict().copy()
                patience_counter = 0
            else:
                patience_counter += 1
            
            if patience_counter >= patience:
                if run == 0 and (epoch + 1) % 10 == 0:
                    print(f"  Early stopping at epoch {epoch+1}")
                break
            
            if run == 0 and (epoch + 1) % 10 == 0:
                print(f"  Epoch {epoch+1}: Loss={loss.item():.4f}, "
                      f"Val Acc={val_acc:.4f}")
        
        if best_state is not None:
            model.load_state_dict(best_state)
        
        model.eval()
        with torch.no_grad():
            _, logits = model(G, features)
        
        use_cs = config.get('use_correct_smooth', False)
        if use_cs:
            logits = apply_correct_and_smooth(
                G, logits, labels, train_mask,
                correct_iter=config.get('correct_iter', 50),
                correct_alpha=config.get('correct_alpha', 0.5),
                smooth_iter=config.get('smooth_iter', 50),
                smooth_alpha=config.get('smooth_alpha', 0.5),
            )
        
        test_preds = logits[test_mask].argmax(dim=1)
        test_acc = (test_preds == labels[test_mask]).float().mean().item()
        
        # Compute additional metrics for node classification
        all_preds = logits[test_mask].argmax(dim=1).cpu().numpy()
        all_targets = labels[test_mask].cpu().numpy()
        from sklearn.metrics import balanced_accuracy_score, precision_score, recall_score, f1_score
        bal_acc = balanced_accuracy_score(all_targets, all_preds)
        prec = precision_score(all_targets, all_preds, average='macro', zero_division=0)
        rec = recall_score(all_targets, all_preds, average='macro', zero_division=0)
        f1_macro = f1_score(all_targets, all_preds, average='macro', zero_division=0)
        
        all_results.append({
            'run': run + 1,
            'test_acc': test_acc,
            'test_balanced_accuracy': bal_acc,
            'test_precision_macro': prec,
            'test_recall_macro': rec,
            'test_f1_macro': f1_macro,
            'best_val_acc': best_val_acc,
        })
        
        print(f"Run {run+1}: Test Acc={test_acc:.4f}")
    
    results_df = pd.DataFrame(all_results)
    
    print(f"\n{'='*50}")
    print(f"Summary ({num_runs} runs)")
    print(f"{'='*50}")
    print(f"Test Acc: {results_df['test_acc'].mean():.4f} ± "
          f"{results_df['test_acc'].std():.4f}")
    
    results_df['dataset'] = dataset_name
    results_df['experiment'] = 'node_classification'
    results_df['timestamp'] = datetime.now().isoformat()
    
    os.makedirs(save_dir, exist_ok=True)
    summary_path = os.path.join(save_dir, f"{dataset_name}_node_classification.csv")
    summary = pd.DataFrame([{
        'dataset': dataset_name,
        'experiment': 'node_classification',
        'mean_test_acc': results_df['test_acc'].mean(),
        'std_test_acc': results_df['test_acc'].std(),
        'mean_test_balanced_accuracy': results_df['test_balanced_accuracy'].mean(),
        'std_test_balanced_accuracy': results_df['test_balanced_accuracy'].std(),
        'mean_test_precision_macro': results_df['test_precision_macro'].mean(),
        'std_test_precision_macro': results_df['test_precision_macro'].std(),
        'mean_test_recall_macro': results_df['test_recall_macro'].mean(),
        'std_test_recall_macro': results_df['test_recall_macro'].std(),
        'mean_test_f1_macro': results_df['test_f1_macro'].mean(),
        'mean_val_acc': results_df['best_val_acc'].mean(),
        'num_runs': num_runs,
    }])
    summary.to_csv(summary_path, index=False)
    print(f"Results saved to {summary_path}")
    
    return results_df


def run_ablation_study(
    dataset_name: str,
    base_config: Dict,
    device: torch.device,
    n_folds: int = 5,
    seed: int = 42,
    save_dir: str = '',
) -> pd.DataFrame:
    """
    Run ablation study matching the specification.
    
    Ablation variants (per spec §6.6.3):
    - full_model: Complete CSG-Transformer
    - w/o PE: No positional encoding
    - w/o Backward: No backward cross attention
    - w/o TNA: Replace TNA with simple neighbor mean
    - Single Layer: L=1
    - w/o RWSE: No random walk structural encoding
    - w/o SE: No structural encoding
    - Sparse Mode: Sparse attention only
    """
    if not save_dir:
        save_dir = os.path.join(_DEFAULT_OUTPUT_DIR, 'ablation')
    set_seed(seed)
    
    print(f"\n{'='*60}")
    print(f"Ablation Study: {dataset_name}")
    print(f"{'='*60}")
    
    data_list, stats = prepare_graph_classification_data(
        dataset_name,
        use_rewiring=base_config.get('use_rewiring', False),
        rewiring_alpha=base_config.get('rewiring_alpha', 0.15),
    )
    
    base_config = base_config.copy()
    base_config['in_dim'] = stats['in_dim']
    base_config['out_dim'] = stats['num_classes']
    base_config['task'] = 'graph_classification'
    
    ablation_variants = {
        'full_model': base_config.copy(),
        'w/o PE': {**base_config.copy(), 'pe_type': 'none', 'use_pe': False},
        'w/o Backward': {**base_config.copy(), 'use_backward': False},
        'w/o TNA': {**base_config.copy(), 'use_tna': False},
        'Single Layer': {**base_config.copy(), 'L': 1},
        'w/o RWSE': {**base_config.copy(), 'pe_type': 'lpe'},
        'w/o SE': {**base_config.copy(), 'pe_type': 'rwse'},
        'Sparse Mode': {**base_config.copy(), 'attn_mode': 'sparse'},
        'w/o Virtual Node': {**base_config.copy(), 'use_virtual_node': False},
        'w/o Rel Bias': {**base_config.copy(), 'use_rel_bias': False},
        'w/ DropEdge': {**base_config.copy(), 'drop_edge_rate': 0.1},
        'w/ Self-KD': {**base_config.copy(), 'use_ema': True, 'kd_weight': 0.5},
    }
    
    results = []
    
    for variant_name, config in ablation_variants.items():
        print(f"\nRunning: {variant_name}")
        
        fold_results = cross_validate(
            data_list, config, device, n_folds=n_folds, verbose=False
        )
        
        results.append({
            'variant': variant_name,
            'dataset': dataset_name,
            'mean_acc': fold_results['test_acc'].mean(),
            'std_acc': fold_results['test_acc'].std(),
            'mean_auc': fold_results['test_auc'].mean(),
            'std_auc': fold_results['test_auc'].std(),
            'mean_f1': fold_results['test_f1'].mean(),
            'std_f1': fold_results['test_f1'].std(),
        })
        
        print(f"  Acc: {results[-1]['mean_acc']:.4f} ± {results[-1]['std_acc']:.4f}")
    
    results_df = pd.DataFrame(results)
    
    os.makedirs(save_dir, exist_ok=True)
    filepath = os.path.join(save_dir, f"{dataset_name}_ablation.csv")
    results_df.to_csv(filepath, index=False)
    print(f"\nAblation results saved to {filepath}")
    
    return results_df


def load_config(config_path: str) -> Dict:
    """Load experiment configuration from YAML file."""
    with open(config_path, 'r') as f:
        config = yaml.safe_load(f)
    
    flat_config = {}
    
    for section in ['model', 'training', 'evaluation', 'dataset', 'hardware']:
        if section in config:
            flat_config.update(config[section])
    
    # Handle device
    if 'device' in flat_config:
        if flat_config['device'] == 'auto':
            flat_config['device'] = get_device()
        elif flat_config['device'] == 'cuda':
            flat_config['device'] = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
        else:
            flat_config['device'] = torch.device('cpu')
    
    return flat_config


def main():
    """Main evaluation entry point."""
    parser = argparse.ArgumentParser(description='CSG-Transformer Evaluation')
    
    parser.add_argument('--task', type=str, default='graph_classification',
                        choices=['graph_classification', 'node_classification',
                                 'ablation', 'all'],
                        help='Evaluation task')
    parser.add_argument('--dataset', type=str, default='MUTAG',
                        help='Dataset name')
    parser.add_argument('--config', type=str, default=None,
                        help='YAML config file path')
    parser.add_argument('--n_folds', type=int, default=10,
                        help='Cross-validation folds')
    parser.add_argument('--num_runs', type=int, default=5,
                        help='Number of runs for node classification')
    parser.add_argument('--seed', type=int, default=42,
                        help='Random seed')
    parser.add_argument('--save_dir', type=str, default='',
                        help='Output directory (default: our_experiments/csg_transformer_eval/outputs/)')
    parser.add_argument('--hidden_dim', type=int, default=128)
    parser.add_argument('--num_heads', type=int, default=4)
    parser.add_argument('--dropout', type=float, default=0.1)
    parser.add_argument('--L', type=int, default=3,
                        help='CSG layers')
    parser.add_argument('--T', type=int, default=3,
                        help='TNA rounds per layer')
    parser.add_argument('--I', type=int, default=5,
                        help='Global iterations')
    parser.add_argument('--epochs', type=int, default=200)
    parser.add_argument('--lr', type=float, default=1e-3)
    parser.add_argument('--patience', type=int, default=30)
    parser.add_argument('--pe_type', type=str, default='composite',
                        choices=['lpe', 'rwse', 'composite', 'none'])
    parser.add_argument('--attn_mode', type=str, default='adaptive',
                        choices=['full', 'sparse', 'adaptive', 'tna_only'])
    parser.add_argument('--pooling', type=str, default='mean',
                        choices=['mean', 'mean+max'],
                        help='Pooling strategy')
    parser.add_argument('--warmup_epochs', type=int, default=10,
                        help='Linear LR warmup epochs')
    parser.add_argument('--use_correct_smooth', action='store_true',
                        help='Enable Correct & Smooth post-processing')
    parser.add_argument('--correct_iter', type=int, default=50)
    parser.add_argument('--correct_alpha', type=float, default=0.5)
    parser.add_argument('--smooth_iter', type=int, default=50)
    parser.add_argument('--smooth_alpha', type=float, default=0.5)
    parser.add_argument('--drop_edge_rate', type=float, default=0.0,
                        help='DropEdge dropout rate')
    parser.add_argument('--use_ema', action='store_true',
                        help='Enable EMA teacher self-knowledge distillation')
    parser.add_argument('--ema_decay', type=float, default=0.999)
    parser.add_argument('--kd_temperature', type=float, default=4.0)
    parser.add_argument('--kd_weight', type=float, default=0.0)
    
    args = parser.parse_args()
    
    if not args.save_dir:
        args.save_dir = _DEFAULT_OUTPUT_DIR
    
    if args.config:
        config = load_config(args.config)
        device = config.get('device', get_device())
    else:
        device = get_device()
        config = {
            'hidden_dim': args.hidden_dim,
            'num_heads': args.num_heads,
            'dropout': args.dropout,
            'L': args.L,
            'T': args.T,
            'I': args.I,
            'epochs': args.epochs,
            'lr': args.lr,
            'patience': args.patience,
            'weight_decay': 1e-4,
            'lambda_reg': 1e-4,
            'min_lr': 1e-6,
            'min_delta': 1e-4,
            'seed': args.seed,
            'pe_type': args.pe_type,
            'attn_mode': args.attn_mode,
            'pooling': args.pooling,
            'warmup_epochs': args.warmup_epochs,
            'grad_clip': 1.0,
            'drop_edge_rate': args.drop_edge_rate,
            'use_ema': args.use_ema,
            'ema_decay': args.ema_decay,
            'kd_temperature': args.kd_temperature,
            'kd_weight': args.kd_weight,
            'use_correct_smooth': args.use_correct_smooth,
            'correct_iter': args.correct_iter,
            'correct_alpha': args.correct_alpha,
            'smooth_iter': args.smooth_iter,
            'smooth_alpha': args.smooth_alpha,
        }
    
    print(f"Device: {device}")
    print(f"Config: {config}")
    
    if args.task == 'graph_classification':
        datasets = [args.dataset] if args.dataset != 'all' else list(GRAPH_CLASSIFICATION_DATASETS.keys())
        for ds_name in datasets:
            run_graph_classification_experiment(
                ds_name, config, device, args.n_folds, args.seed,
                os.path.join(args.save_dir, 'graph_classification'),
            )
    
    elif args.task == 'node_classification':
        datasets = [args.dataset] if args.dataset != 'all' else list(NODE_CLASSIFICATION_DATASETS.keys())
        for ds_name in datasets:
            run_node_classification_experiment(
                ds_name, config, device, args.seed,
                os.path.join(args.save_dir, 'node_classification'),
                num_runs=args.num_runs,
            )
    
    elif args.task == 'ablation':
        run_ablation_study(
            args.dataset, config, device, args.n_folds, args.seed,
            os.path.join(args.save_dir, 'ablation'),
        )
    
    elif args.task == 'all':
        print("\n" + "="*60)
        print("RUNNING ALL EXPERIMENTS")
        print("="*60)
        
        for ds_name in list(GRAPH_CLASSIFICATION_DATASETS.keys())[:3]:
            run_graph_classification_experiment(
                ds_name, config, device, args.n_folds, args.seed,
                os.path.join(args.save_dir, 'graph_classification'),
            )
        
        run_ablation_study(
            'MUTAG', config, device, min(args.n_folds, 5), args.seed,
            os.path.join(args.save_dir, 'ablation'),
        )


if __name__ == '__main__':
    main()
