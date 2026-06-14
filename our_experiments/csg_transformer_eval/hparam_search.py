"""
Hyperparameter search for CSG-Transformer.

Performs random/grid search over the parameter space defined in
csg_transformer方案.md §6.3.
"""

import itertools
import os
import random
from typing import Any, Callable, Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
import torch

from our_experiments.csg_transformer_eval.train import cross_validate


HPARAM_SPACE: Dict[str, List[Any]] = {
    "hidden_dim": [64, 128, 256],
    "num_heads": [2, 4, 8],
    "L": [1, 2, 3, 4],
    "T": [1, 3, 5],
    "I": [1, 3, 5],
    "dropout": [0.0, 0.1, 0.2],
    "lr": [5e-4, 1e-3, 5e-3],
    "weight_decay": [0, 1e-4, 5e-4],
    "pe_type": ["composite", "lpe", "rwse"],
    "attn_mode": ["adaptive", "full", "sparse"],
    "use_backward": [True, False],
    "use_tna": [True, False],
}


def sample_random_config(
    base_config: Dict[str, Any],
    search_space: Optional[Dict[str, List[Any]]] = None,
) -> Dict[str, Any]:
    """Sample a random hyperparameter configuration.
    
    Args:
        base_config: Base config to extend (model parameters preserved)
        search_space: Dict mapping param name to list of choices
    
    Returns:
        Extended config with sampled hyperparameters
    """
    if search_space is None:
        search_space = HPARAM_SPACE

    config = base_config.copy()
    for param, choices in search_space.items():
        config[param] = random.choice(choices)
    return config


def grid_configs(
    base_config: Dict[str, Any],
    fixed_params: Optional[Dict[str, Any]] = None,
    search_params: Optional[Dict[str, List[Any]]] = None,
) -> List[Dict[str, Any]]:
    """Generate all grid combinations over search parameters.
    
    Args:
        base_config: Base config to extend
        fixed_params: Params with fixed values (not searched)
        search_params: Params to grid-search over
    
    Returns:
        List of config dicts
    """
    if fixed_params is None:
        fixed_params = {}
    if search_params is None:
        search_params = {
            k: v for k, v in HPARAM_SPACE.items() if k not in fixed_params
        }

    keys = list(search_params.keys())
    value_lists = [search_params[k] for k in keys]

    configs = []
    for values in itertools.product(*value_lists):
        config = base_config.copy()
        config.update(fixed_params)
        for k, v in zip(keys, values):
            config[k] = v
        configs.append(config)

    return configs


def run_hparam_search(
    dataset: List[Tuple],
    base_config: Dict[str, Any],
    device: torch.device,
    n_folds: int = 5,
    n_trials: int = 20,
    search_mode: str = "random",
    search_space: Optional[Dict[str, List[Any]]] = None,
    fixed_params: Optional[Dict[str, Any]] = None,
    save_dir: str = "",
    verbose: bool = True,
) -> pd.DataFrame:
    """Run hyperparameter search.
    
    Args:
        dataset: List of (G, features, label) tuples
        base_config: Base training config
        device: Device
        n_folds: CV folds (lower for speed)
        n_trials: Number of random trials (used only in 'random' mode)
        search_mode: 'random' or 'grid'
        search_space: Param space to search
        fixed_params: Params held fixed
        save_dir: Optional save directory
        verbose: Print progress
    
    Returns:
        DataFrame with trial results sorted by val_acc
    """
    if search_space is None:
        search_space = HPARAM_SPACE
    if fixed_params is None:
        fixed_params = {}

    if search_mode == "grid":
        configs = grid_configs(base_config, fixed_params, search_space)
        if verbose:
            print(f"Grid search: {len(configs)} total configurations")
    else:
        configs = [
            sample_random_config(base_config, search_space)
            for _ in range(n_trials)
        ]
        if verbose:
            print(f"Random search: {n_trials} trials")

    trial_results = []
    for trial_idx, config in enumerate(configs):
        trial_config = {**base_config, **config}

        if verbose:
            print(f"\n--- Trial {trial_idx + 1}/{len(configs)} ---")
            print(f"  Config: {trial_config}")

        try:
            fold_results = cross_validate(
                dataset, trial_config, device,
                n_folds=n_folds, verbose=False,
            )

            trial_results.append({
                "trial": trial_idx + 1,
                "mean_val_acc": fold_results["val_acc"].mean(),
                "std_val_acc": fold_results["val_acc"].std(),
                "mean_test_acc": fold_results["test_acc"].mean(),
                "std_test_acc": fold_results["test_acc"].std(),
                "mean_test_auc": fold_results["test_auc"].mean(),
                "mean_test_f1": fold_results["test_f1"].mean(),
                **{k: str(config[k]) for k in search_space},
            })

            if verbose:
                print(
                    f"  Val Acc: {fold_results['val_acc'].mean():.4f} "
                    f"± {fold_results['val_acc'].std():.4f}"
                )
        except Exception as e:
            if verbose:
                print(f"  FAILED: {e}")
            continue

    if not trial_results:
        return pd.DataFrame()

    results_df = pd.DataFrame(trial_results)
    results_df = results_df.sort_values("mean_val_acc", ascending=False)

    if verbose:
        print(f"\n{'=' * 60}")
        print("Hyperparameter Search Results")
        print(f"{'=' * 60}")
        print("Top 5 configurations:")
        for i in range(min(5, len(results_df))):
            row = results_df.iloc[i]
            print(
                f"  #{i+1}: Val Acc={row['mean_val_acc']:.4f}, "
                f"Test Acc={row['mean_test_acc']:.4f}, "
                f"hidden_dim={row.get('hidden_dim', '?')}, "
                f"L={row.get('L', '?')}, T={row.get('T', '?')}, I={row.get('I', '?')}"
            )

    if save_dir:
        os.makedirs(save_dir, exist_ok=True)
        path = os.path.join(save_dir, "hparam_search_results.csv")
        results_df.to_csv(path, index=False)
        print(f"Results saved to {path}")

        best_path = os.path.join(save_dir, "best_config.json")
        best = results_df.iloc[0].to_dict()
        import json
        with open(best_path, "w") as f:
            json.dump(best, f, indent=2, default=str)
        print(f"Best config saved to {best_path}")

    return results_df
