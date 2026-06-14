"""
Comprehensive experiment runner for CSG-Transformer.

Runs all experiments described in the specification (csg_transformer方案.md):
1. Graph classification on all datasets (10-fold CV)
2. Node classification on all datasets (5 runs)
3. Ablation study on MUTAG
4. CFI graph experiments

Outputs are saved to our_experiments/csg_transformer_eval/outputs/

Usage:
    python run_all_experiments.py                          # Run all experiments
    python run_all_experiments.py --mode graph             # Graph classification only
    python run_all_experiments.py --mode node              # Node classification only
    python run_all_experiments.py --mode ablation          # Ablation study only
    python run_all_experiments.py --mode cfi               # CFI experiments only
    python run_all_experiments.py --datasets MUTAG,NCI1    # Specific datasets
"""

import argparse
import os
import sys
from datetime import datetime

import torch

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))

from our_experiments.csg_transformer_eval.evaluate import (
    run_graph_classification_experiment,
    run_node_classification_experiment,
    run_ablation_study,
    GRAPH_CLASSIFICATION_DATASETS,
    NODE_CLASSIFICATION_DATASETS,
    _DEFAULT_OUTPUT_DIR,
)
from our_experiments.csg_transformer_eval.cfi_experiments import (
    run_all_cfi_experiments,
)
from our_experiments.csg_transformer_eval.utils import set_seed, get_device


def run_all_graph_classification(config, device, save_dir, n_folds=10, seed=42):
    """Run graph classification on all datasets."""
    print("\n" + "=" * 70)
    print("GRAPH CLASSIFICATION EXPERIMENTS (all datasets)")
    print("=" * 70)

    datasets = list(GRAPH_CLASSIFICATION_DATASETS.keys())
    all_results = {}

    for ds_name in datasets:
        try:
            result = run_graph_classification_experiment(
                ds_name, config, device, n_folds=n_folds,
                seed=seed, save_dir=save_dir,
            )
            all_results[ds_name] = {
                'mean_acc': result['test_acc'].mean() if 'test_acc' in result.columns else 0,
                'std_acc': result['test_acc'].std() if 'test_acc' in result.columns else 0,
                'mean_auc': result['test_auc'].mean() if 'test_auc' in result.columns else 0,
                'mean_f1': result['test_f1'].mean() if 'test_f1' in result.columns else 0,
            }
            print(f"  {ds_name}: {all_results[ds_name]['mean_acc']:.2f} +/- "
                  f"{all_results[ds_name]['std_acc']:.2f}")
        except Exception as e:
            print(f"  {ds_name}: FAILED - {e}")
            import traceback
            traceback.print_exc()

    return all_results


def run_all_node_classification(config, device, save_dir, num_runs=5, seed=42):
    """Run node classification on all datasets."""
    print("\n" + "=" * 70)
    print("NODE CLASSIFICATION EXPERIMENTS (all datasets)")
    print("=" * 70)

    datasets = list(NODE_CLASSIFICATION_DATASETS.keys())
    all_results = {}

    for ds_name in datasets:
        if ds_name == 'Pubmed':
            continue  # Skip large datasets for quick testing
        try:
            result = run_node_classification_experiment(
                ds_name, config, device, seed=seed,
                save_dir=save_dir, num_runs=num_runs,
            )
            all_results[ds_name] = {
                'mean_acc': result['test_acc'].mean(),
                'std_acc': result['test_acc'].std(),
            }
            print(f"  {ds_name}: {all_results[ds_name]['mean_acc']:.2f} +/- "
                  f"{all_results[ds_name]['std_acc']:.2f}")
        except Exception as e:
            print(f"  {ds_name}: FAILED - {e}")
            import traceback
            traceback.print_exc()

    return all_results


def main():
    parser = argparse.ArgumentParser(description='Run all CSG-Transformer experiments')

    parser.add_argument('--mode', type=str, default='all',
                        choices=['all', 'graph', 'node', 'ablation', 'cfi'],
                        help='Which experiments to run')
    parser.add_argument('--datasets', type=str, default='',
                        help='Comma-separated dataset names (default: all)')
    parser.add_argument('--n_folds', type=int, default=10,
                        help='Cross-validation folds')
    parser.add_argument('--num_runs', type=int, default=5,
                        help='Runs for node classification')
    parser.add_argument('--hidden_dim', type=int, default=128)
    parser.add_argument('--num_heads', type=int, default=4)
    parser.add_argument('--dropout', type=float, default=0.1)
    parser.add_argument('--L', type=int, default=3)
    parser.add_argument('--T', type=int, default=3)
    parser.add_argument('--I', type=int, default=5)
    parser.add_argument('--epochs', type=int, default=200)
    parser.add_argument('--lr', type=float, default=1e-3)
    parser.add_argument('--seed', type=int, default=42)
    parser.add_argument('--pe_type', type=str, default='composite')
    parser.add_argument('--attn_mode', type=str, default='adaptive')
    parser.add_argument('--cfi_pairs', type=int, default=100,
                        help='Number of CFI graph pairs')

    args = parser.parse_args()

    device = get_device()
    set_seed(args.seed)

    print(f"CSG-Transformer Comprehensive Evaluation")
    print(f"Device: {device}")
    print(f"Timestamp: {datetime.now().isoformat()}")
    print(f"Config: hidden_dim={args.hidden_dim}, L={args.L}, T={args.T}, I={args.I}")

    base_config = {
        'hidden_dim': args.hidden_dim,
        'num_heads': args.num_heads,
        'dropout': args.dropout,
        'L': args.L,
        'T': args.T,
        'I': args.I,
        'epochs': args.epochs,
        'lr': args.lr,
        'weight_decay': 1e-4,
        'lambda_reg': 1e-4,
        'patience': 30,
        'min_lr': 1e-6,
        'min_delta': 1e-4,
        'seed': args.seed,
        'pe_type': args.pe_type,
        'attn_mode': args.attn_mode,
    }

    save_dir = _DEFAULT_OUTPUT_DIR
    os.makedirs(save_dir, exist_ok=True)

    if args.mode in ('all', 'graph'):
        run_all_graph_classification(
            base_config, device,
            os.path.join(save_dir, 'graph_classification'),
            n_folds=args.n_folds, seed=args.seed,
        )

    if args.mode in ('all', 'node'):
        node_config = base_config.copy()
        node_config.update({'lr': 0.01, 'weight_decay': 5e-4, 'patience': 50})
        run_all_node_classification(
            node_config, device,
            os.path.join(save_dir, 'node_classification'),
            num_runs=args.num_runs, seed=args.seed,
        )

    if args.mode in ('all', 'ablation'):
        print("\n" + "=" * 70)
        print("ABLATION STUDY")
        print("=" * 70)
        run_ablation_study(
            'MUTAG', base_config, device,
            n_folds=min(args.n_folds, 5), seed=args.seed,
        )

    if args.mode in ('all', 'cfi'):
        print("\n" + "=" * 70)
        print("CFI EXPERIMENTS")
        print("=" * 70)
        cfi_config = base_config.copy()
        cfi_config.update({'in_dim': 2, 'epochs': 100})
        run_all_cfi_experiments(
            cfi_config, device,
            os.path.join(save_dir, 'cfi'),
        )

    print("\n" + "=" * 70)
    print("ALL EXPERIMENTS COMPLETE")
    print(f"Results saved to {save_dir}")
    print("=" * 70)


if __name__ == '__main__':
    main()
