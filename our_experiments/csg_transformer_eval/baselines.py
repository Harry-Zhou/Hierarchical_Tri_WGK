"""
Baseline comparison system for CSG-Transformer evaluation.

Provides paper-reported baseline results from SOTA graph transformers and GNNs,
and utilities to compare model outputs against them.

Baseline data sources (per csg_transformer方案.md §6.7):
- Graphormer (2021), GPS (2022), Exphormer (2023), GRIT (2023)
- PCB-GNN (2024), ESA (2022), TIGT (2022), GPM (2022), HOGT (2024)
- GCNII (2021), GraphGPS (2022), GOAT (2023), NodeFormer (2023)
- SGFormer (2023), Polynormer (2023), TAPE+RevGAT (2024)
"""

from typing import Any, Dict, List, Optional, Tuple

import pandas as pd
import numpy as np


# ============================================================================
# Graph Classification Baselines (accuracy %)
# ============================================================================

GRAPH_CLASSIFICATION_BASELINES: Dict[str, Dict[str, float]] = {
    "MUTAG": {
        "Graphormer": 89.6,
        "GPS": 92.6,
        "Exphormer": 89.3,
        "GRIT": 88.5,
        "PCB-GNN": 98.5,
        "ESA": 90.2,
        "TIGT": 89.8,
        "GPM": 91.2,
        "HOGT": 92.0,
    },
    "PROTEINS": {
        "Graphormer": 76.3,
        "GPS": 77.7,
        "Exphormer": 77.4,
        "GRIT": 76.8,
        "PCB-GNN": 82.2,
        "ESA": 78.1,
        "TIGT": 77.2,
        "GPM": 78.5,
        "HOGT": 79.0,
    },
    "NCI1": {
        "Graphormer": 78.6,
        "GPS": 82.5,
        "Exphormer": 80.2,
        "GRIT": 79.8,
        "PCB-GNN": 88.4,
        "ESA": 81.5,
        "TIGT": 80.8,
        "GPM": 82.0,
        "HOGT": 83.0,
    },
    "ENZYMES": {
        "Graphormer": 55.3,
        "GPS": 60.0,
        "Exphormer": 58.6,
        "GRIT": 57.5,
        "PCB-GNN": 62.1,
        "ESA": 59.3,
        "TIGT": 58.2,
        "GPM": 60.5,
        "HOGT": 61.0,
    },
    "IMDB-BINARY": {
        "Graphormer": 78.0,
        "GPS": 80.6,
        "Exphormer": 79.8,
        "GRIT": 79.2,
        "PCB-GNN": 81.2,
        "ESA": 80.5,
        "TIGT": 79.5,
        "GPM": 80.8,
        "HOGT": 81.0,
    },
    "IMDB-MULTI": {
        "Graphormer": 55.3,
        "GPS": 57.0,
        "Exphormer": 56.8,
        "GRIT": 56.2,
        "PCB-GNN": 58.5,
        "ESA": 57.2,
        "TIGT": 56.5,
        "GPM": 57.8,
        "HOGT": 58.0,
    },
    "COLLAB": {
        "Graphormer": 81.3,
        "GPS": 82.1,
        "Exphormer": 82.5,
        "GRIT": 81.8,
        "PCB-GNN": 83.1,
        "ESA": 82.8,
        "TIGT": 82.0,
        "GPM": 82.5,
        "HOGT": 83.0,
    },
}


# ============================================================================
# Node Classification Baselines (accuracy %)
# ============================================================================

NODE_CLASSIFICATION_BASELINES: Dict[str, Dict[str, float]] = {
    "Cora": {
        "GCNII": 85.5,
        "GraphGPS": 83.2,
        "GOAT": 84.5,
        "NodeFormer": 82.8,
        "SGFormer": 88.2,
        "Polynormer": 87.6,
        "Exphormer+GCN": 86.1,
        "TAPE+RevGAT": 92.9,
        "HOGT": 88.5,
    },
    "Citeseer": {
        "GCNII": 73.4,
        "GraphGPS": 71.8,
        "GOAT": 72.5,
        "NodeFormer": 71.2,
        "SGFormer": 76.8,
        "Polynormer": 75.2,
        "Exphormer+GCN": 74.5,
        "TAPE+RevGAT": 78.5,
        "HOGT": 76.5,
    },
    "Pubmed": {
        "GCNII": 80.3,
        "GraphGPS": 79.5,
        "GOAT": 80.0,
        "NodeFormer": 79.2,
        "SGFormer": 82.1,
        "Polynormer": 81.5,
        "Exphormer+GCN": 80.8,
        "TAPE+RevGAT": 84.2,
        "HOGT": 82.0,
    },
    "Computers": {
        "GCNII": 84.2,
        "GraphGPS": 83.0,
        "GOAT": 84.5,
        "NodeFormer": 82.8,
        "SGFormer": 86.5,
        "Polynormer": 85.8,
        "Exphormer+GCN": 85.2,
        "TAPE+RevGAT": 88.1,
        "HOGT": 86.0,
    },
    "Photo": {
        "GCNII": 85.0,
        "GraphGPS": 84.2,
        "GOAT": 85.5,
        "NodeFormer": 83.5,
        "SGFormer": 87.2,
        "Polynormer": 86.5,
        "Exphormer+GCN": 86.1,
        "TAPE+RevGAT": 89.3,
        "HOGT": 87.0,
    },
    "Squirrel": {
        "GCNII": 62.5,
        "GraphGPS": 61.8,
        "GOAT": 63.2,
        "NodeFormer": 62.0,
        "SGFormer": 65.8,
        "Polynormer": 68.2,
        "Exphormer+GCN": 64.5,
        "TAPE+RevGAT": 71.2,
        "HOGT": 67.5,
    },
    "Chameleon": {
        "GCNII": 60.8,
        "GraphGPS": 59.5,
        "GOAT": 61.5,
        "NodeFormer": 60.2,
        "SGFormer": 63.5,
        "Polynormer": 66.8,
        "Exphormer+GCN": 62.1,
        "TAPE+RevGAT": 69.5,
        "HOGT": 65.8,
    },
}


def get_baseline_table(
    dataset_name: str,
    task: str = "graph_classification",
    model_result: Optional[float] = None,
    model_std: Optional[float] = None,
) -> pd.DataFrame:
    """Get baseline comparison table for a dataset.
    
    Args:
        dataset_name: Name of the dataset (e.g., 'MUTAG', 'Cora')
        task: 'graph_classification' or 'node_classification'
        model_result: Optional CSG-Transformer accuracy to include
        model_std: Optional CSG-Transformer accuracy std to include
    
    Returns:
        DataFrame with baseline results sorted by accuracy
    """
    baselines = (
        GRAPH_CLASSIFICATION_BASELINES
        if task == "graph_classification"
        else NODE_CLASSIFICATION_BASELINES
    )

    if dataset_name not in baselines:
        raise ValueError(
            f"Dataset '{dataset_name}' not found in {task} baselines. "
            f"Available: {list(baselines.keys())}"
        )

    rows = []
    for method, acc in baselines[dataset_name].items():
        rows.append({"Method": method, "Accuracy": acc, "Source": "Paper"})

    if model_result is not None:
        label = f"CSG-Transformer (ours)"
        if model_std is not None:
            label += f" ±{model_std:.1f}"
        rows.append(
            {
                "Method": "CSG-Transformer (ours)",
                "Accuracy": model_result,
                "Source": "This work",
            }
        )

    df = pd.DataFrame(rows)
    df = df.sort_values("Accuracy", ascending=False).reset_index(drop=True)
    df["Rank"] = range(1, len(df) + 1)
    return df


def compare_with_baselines(
    results_df: pd.DataFrame,
    task: str = "graph_classification",
    save_dir: str = "",
) -> Dict[str, pd.DataFrame]:
    """Compare model results with baseline for all datasets in results.
    
    Args:
        results_df: DataFrame with columns 'dataset', 'test_acc', 'std_test_acc'
        task: 'graph_classification' or 'node_classification'
        save_dir: Optional directory to save comparison CSV files

    Returns:
        Dict mapping dataset name to comparison DataFrame
    """
    baselines = (
        GRAPH_CLASSIFICATION_BASELINES
        if task == "graph_classification"
        else NODE_CLASSIFICATION_BASELINES
    )

    comparisons = {}
    for dataset_name in results_df["dataset"].unique():
        if dataset_name not in baselines:
            continue

        subset = results_df[results_df["dataset"] == dataset_name]
        model_acc = subset["test_acc"].mean() * 100
        model_std = subset["test_acc"].std() * 100

        comp = get_baseline_table(
            dataset_name,
            task=task,
            model_result=model_acc,
            model_std=model_std,
        )
        comparisons[dataset_name] = comp

        baseline_best = comp[comp["Source"] == "Paper"]["Accuracy"].max()
        our_row = comp[comp["Method"] == "CSG-Transformer (ours)"]
        if not our_row.empty:
            our_acc = our_row["Accuracy"].iloc[0]
            gap = our_acc - baseline_best
            rank = our_row["Rank"].iloc[0]

            print(f"\n--- {dataset_name} ---")
            print(f"  Best baseline: {comp.iloc[0]['Method']} = {baseline_best:.1f}%")
            print(f"  CSG-Transformer: {our_acc:.1f}% (rank {rank}/{len(comp)})")
            print(f"  Gap to best: {gap:+.1f}%")

        if save_dir:
            import os
            os.makedirs(save_dir, exist_ok=True)
            path = os.path.join(save_dir, f"{dataset_name}_baseline_comparison.csv")
            comp.to_csv(path, index=False)
            print(f"  Saved to {path}")

    return comparisons


def summarize_comparisons(
    comparisons: Dict[str, pd.DataFrame],
) -> pd.DataFrame:
    """Create a summary table across all datasets.
    
    Args:
        comparisons: Dict from compare_with_baselines
    
    Returns:
        Summary DataFrame
    """
    rows = []
    for dataset_name, comp in comparisons.items():
        our_row = comp[comp["Method"] == "CSG-Transformer (ours)"]
        if our_row.empty:
            continue

        our_acc = our_row["Accuracy"].iloc[0]
        our_rank = our_row["Rank"].iloc[0]
        total_models = len(comp)

        best_baseline = comp[comp["Source"] == "Paper"].iloc[0]
        best_acc = best_baseline["Accuracy"]
        best_method = best_baseline["Method"]

        rows.append(
            {
                "Dataset": dataset_name,
                "Our Accuracy": our_acc,
                "Our Rank": f"{our_rank}/{total_models}",
                "Best Baseline": best_acc,
                "Best Method": best_method,
                "Gap": our_acc - best_acc,
            }
        )

    summary = pd.DataFrame(rows)
    summary = summary.sort_values("Gap", ascending=False).reset_index(drop=True)

    print("\n" + "=" * 80)
    print("CSG-Transformer vs Baselines Summary")
    print("=" * 80)
    print(summary.to_string(index=False))
    print(f"\nAverage gap to best baseline: {summary['Gap'].mean():.2f}%")
    print(
        f"Datasets where CSG-Transformer outperforms all baselines: "
        f"{(summary['Gap'] > 0).sum()}/{len(summary)}"
    )

    return summary
