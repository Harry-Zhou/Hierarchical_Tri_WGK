"""
Results aggregation and visualization for CSG-Transformer evaluation.

Loads experiment outputs from disk and produces:
1. Master summary tables (graph classification, node classification, ablation, CFI)
2. Comparison against paper-reported baselines
3. Structured data for external visualization

Usage:
    python results_aggregator.py                              # Aggregate all results
    python results_aggregator.py --output-dir ./outputs        # Custom output dir
    python results_aggregator.py --mode graph                  # Graph classification only
"""

import argparse
import json
import os
import sys
from datetime import datetime
from typing import Any, Dict, List, Optional, Tuple

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))

from our_experiments.csg_transformer_eval.baselines import (
    get_baseline_table, summarize_comparisons,
    GRAPH_CLASSIFICATION_BASELINES, NODE_CLASSIFICATION_BASELINES,
)
from our_experiments.csg_transformer_eval.evaluate import _DEFAULT_OUTPUT_DIR

# ============================================================================
# Result Loaders
# ============================================================================

def _find_files(root_dir: str, pattern: str) -> List[str]:
    """Find files matching pattern under root_dir."""
    matches = []
    for dirpath, _, filenames in os.walk(root_dir):
        for f in filenames:
            if f.endswith(pattern):
                matches.append(os.path.join(dirpath, f))
    return sorted(matches)


def load_graph_classification_results(
    results_dir: str,
) -> pd.DataFrame:
    """
    Load graph classification results from CSV files.

    Expects files named like: {dataset}_cv_results.csv or similar
    under results_dir/graph_classification/.

    Returns a DataFrame with columns: dataset, mean_acc, std_acc, mean_auc, ...
    """
    graph_dir = os.path.join(results_dir, 'graph_classification')
    if not os.path.isdir(graph_dir):
        print(f"Warning: graph classification directory not found: {graph_dir}")
        return pd.DataFrame()

    rows = []
    for dirpath, _, filenames in os.walk(graph_dir):
        for f in filenames:
            if not f.endswith('_cv_results.csv') and not f.endswith('_graph_classification.csv'):
                continue
            ds_name = f.replace('_cv_results.csv', '').replace('_graph_classification.csv', '')
            try:
                df = pd.read_csv(os.path.join(dirpath, f))
            except Exception as e:
                print(f"  Warning: failed to read {f}: {e}")
                continue

            row = {'dataset': ds_name}
            for col in df.columns:
                if df[col].dtype in (np.float64, np.float32, np.int64):
                    row[f'mean_{col}'] = df[col].mean()
                    row[f'std_{col}'] = df[col].std()
            rows.append(row)

    result = pd.DataFrame(rows)
    if not result.empty:
        result = result.sort_values('dataset').reset_index(drop=True)
    return result


def load_node_classification_results(
    results_dir: str,
) -> pd.DataFrame:
    """
    Load node classification results.

    Expects files named like: {dataset}_node_classification.csv
    under results_dir/node_classification/.
    """
    node_dir = os.path.join(results_dir, 'node_classification')
    if not os.path.isdir(node_dir):
        print(f"Warning: node classification directory not found: {node_dir}")
        return pd.DataFrame()

    rows = []
    for dirpath, _, filenames in os.walk(node_dir):
        for f in filenames:
            if not f.endswith('_node_classification.csv'):
                continue
            ds_name = f.replace('_node_classification.csv', '')
            try:
                df = pd.read_csv(os.path.join(dirpath, f))
            except Exception as e:
                print(f"  Warning: failed to read {f}: {e}")
                continue

            row = {'dataset': ds_name}
            for col in df.columns:
                if df[col].dtype in (np.float64, np.float32, np.int64):
                    row[f'mean_{col}'] = df[col].mean()
                    row[f'std_{col}'] = df[col].std()
            rows.append(row)

    result = pd.DataFrame(rows)
    if not result.empty:
        result = result.sort_values('dataset').reset_index(drop=True)
    return result


def load_ablation_results(
    results_dir: str,
) -> pd.DataFrame:
    """
    Load ablation study results.

    Expects files like: {dataset}_ablation.csv under results_dir/ablation/.
    """
    ablation_dir = os.path.join(results_dir, 'ablation')
    if not os.path.isdir(ablation_dir):
        print(f"Warning: ablation directory not found: {ablation_dir}")
        return pd.DataFrame()

    rows = []
    for dirpath, _, filenames in os.walk(ablation_dir):
        for f in filenames:
            if not f.endswith('_ablation.csv'):
                continue
            filepath = os.path.join(dirpath, f)
            ds_name = f.replace('_ablation.csv', '')
            try:
                df = pd.read_csv(filepath)
                df['dataset'] = ds_name
                rows.append(df)
            except Exception as e:
                print(f"  Warning: failed to read {f}: {e}")
                continue

    if not rows:
        return pd.DataFrame()
    result = pd.concat(rows, ignore_index=True)
    return result


def load_cfi_results(
    results_dir: str,
) -> Dict[str, Any]:
    """
    Load CFI experiment results.

    Expects files like: cfi_{base_graph}_results.json and
    cfi_layer_ablation.csv, cfi_ablation_study.csv under results_dir/cfi/.
    """
    cfi_dir = os.path.join(results_dir, 'cfi')
    if not os.path.isdir(cfi_dir):
        print(f"Warning: CFI directory not found: {cfi_dir}")
        return {}

    results = {}

    # Load per-base-graph accuracy JSONs
    base_accuracies = {}
    for f in sorted(os.listdir(cfi_dir)):
        if f.startswith('cfi_') and f.endswith('_results.json'):
            filepath = os.path.join(cfi_dir, f)
            try:
                with open(filepath) as fh:
                    data = json.load(fh)
                base_name = data.get('base_graph', f.replace('cfi_', '').replace('_results.json', ''))
                base_accuracies[base_name] = {
                    'accuracy': data.get('accuracy', 0),
                    'consistency': data.get('consistency', 0),
                    'correct': data.get('correct', 0),
                    'total': data.get('total', 0),
                }
            except Exception as e:
                print(f"  Warning: failed to load {f}: {e}")
    if base_accuracies:
        results['base_accuracy'] = base_accuracies

    # Load layer ablation
    layer_file = os.path.join(cfi_dir, 'cfi_layer_ablation.csv')
    if os.path.isfile(layer_file):
        try:
            results['layer_ablation'] = pd.read_csv(layer_file)
        except Exception as e:
            print(f"  Warning: failed to load layer ablation: {e}")

    # Load component ablation
    abl_file = os.path.join(cfi_dir, 'cfi_ablation_study.csv')
    if os.path.isfile(abl_file):
        try:
            results['component_ablation'] = pd.read_csv(abl_file)
        except Exception as e:
            print(f"  Warning: failed to load component ablation: {e}")

    return results


# ============================================================================
# Summary Generation
# ============================================================================

def generate_graph_classification_summary(
    results_df: pd.DataFrame,
) -> pd.DataFrame:
    """Add baseline comparison to graph classification results."""
    if results_df.empty:
        return results_df

    summary_rows = []
    for _, row in results_df.iterrows():
        ds = row['dataset']
        acc_col = next((c for c in results_df.columns if 'acc' in c.lower()), None)
        if acc_col is None:
            continue
        our_mean = row.get(f'mean_{acc_col}', row.get('mean_test_acc', None))
        our_std = row.get(f'std_{acc_col}', row.get('std_test_acc', None))
        if our_mean is None:
            continue

        entry = {
            'dataset': ds,
            'CSG-Transformer': f"{our_mean:.2f} +/- {our_std:.2f}" if our_std else f"{our_mean:.2f}",
        }

        # Add best baseline
        if ds in GRAPH_CLASSIFICATION_BASELINES:
            bsl = GRAPH_CLASSIFICATION_BASELINES[ds]
            best_name = max(bsl, key=bsl.get)
            best_val = bsl[best_name]
            entry['Best Baseline'] = f"{best_val:.1f} ({best_name})"
            entry['Gap to Best'] = f"{our_mean - best_val:.2f}"
            entry['Rank vs Baselines'] = sum(1 for v in bsl.values() if v < our_mean) + 1
        else:
            entry['Best Baseline'] = 'N/A'
            entry['Gap to Best'] = 'N/A'
            entry['Rank vs Baselines'] = 'N/A'

        summary_rows.append(entry)

    return pd.DataFrame(summary_rows)


def generate_node_classification_summary(
    results_df: pd.DataFrame,
) -> pd.DataFrame:
    """Add baseline comparison to node classification results."""
    if results_df.empty:
        return results_df

    summary_rows = []
    for _, row in results_df.iterrows():
        ds = row['dataset']
        acc_col = next((c for c in results_df.columns if 'acc' in c.lower()), None)
        if acc_col is None:
            continue
        our_mean = row.get(f'mean_{acc_col}', row.get('mean_test_acc', None))
        our_std = row.get(f'std_{acc_col}', row.get('std_test_acc', None))
        if our_mean is None:
            continue

        entry = {
            'dataset': ds,
            'CSG-Transformer': f"{our_mean:.2f} +/- {our_std:.2f}" if our_std else f"{our_mean:.2f}",
        }

        if ds in NODE_CLASSIFICATION_BASELINES:
            bsl = NODE_CLASSIFICATION_BASELINES[ds]
            best_name = max(bsl, key=bsl.get)
            best_val = bsl[best_name]
            entry['Best Baseline'] = f"{best_val:.1f} ({best_name})"
            entry['Gap to Best'] = f"{our_mean - best_val:.2f}"
            entry['Rank vs Baselines'] = sum(1 for v in bsl.values() if v < our_mean) + 1
        else:
            entry['Best Baseline'] = 'N/A'
            entry['Gap to Best'] = 'N/A'
            entry['Rank vs Baselines'] = 'N/A'

        summary_rows.append(entry)

    return pd.DataFrame(summary_rows)


def generate_cfi_summary(
    cfi_results: Dict[str, Any],
) -> Dict[str, Any]:
    """Generate human-readable summary from CFI results."""
    summary = {}

    if 'base_accuracy' in cfi_results:
        rows = []
        for base, metrics in cfi_results['base_accuracy'].items():
            acc = metrics['accuracy']
            cons = metrics.get('consistency', 0)
            correct = metrics['correct']
            total = metrics['total']
            verdict = "PASS" if acc > 0.9 else ("BORDERLINE" if acc > 0.5 else "FAIL")
            rows.append({
                'base_graph': base,
                'accuracy': f"{acc:.4f}",
                'consistency': f"{cons:.4f}",
                'correct/total': f"{correct}/{total}",
                'verdict': verdict,
            })
        summary['base_accuracy'] = pd.DataFrame(rows)

    if 'layer_ablation' in cfi_results:
        df = cfi_results['layer_ablation']
        summary['layer_ablation'] = df

    if 'component_ablation' in cfi_results:
        df = cfi_results['component_ablation']
        summary['component_ablation'] = df

    return summary


# ============================================================================
# Master Aggregation
# ============================================================================

def aggregate_all(output_dir: str) -> Tuple[Dict[str, Any], List[str]]:
    """Load all results and produce combined summary."""
    print(f"Aggregating results from: {output_dir}")
    print(f"Timestamp: {datetime.now().isoformat()}")
    print("=" * 60)

    output = {}
    missing = []

    # Graph classification
    graph_df = load_graph_classification_results(output_dir)
    if not graph_df.empty:
        graph_summary = generate_graph_classification_summary(graph_df)
        output['graph_classification'] = {
            'raw': graph_df,
            'summary': graph_summary,
        }
        print(f"\nGraph Classification: {len(graph_summary)} datasets loaded")
    else:
        missing.append('graph_classification')
        print("\nGraph Classification: no results found")

    # Node classification
    node_df = load_node_classification_results(output_dir)
    if not node_df.empty:
        node_summary = generate_node_classification_summary(node_df)
        output['node_classification'] = {
            'raw': node_df,
            'summary': node_summary,
        }
        print(f"Node Classification: {len(node_summary)} datasets loaded")
    else:
        missing.append('node_classification')
        print("Node Classification: no results found")

    # Ablation
    ablation_df = load_ablation_results(output_dir)
    if not ablation_df.empty:
        output['ablation'] = {'raw': ablation_df}
        print(f"Ablation Study: {len(ablation_df)} rows loaded")
    else:
        missing.append('ablation')
        print("Ablation Study: no results found")

    # CFI
    cfi_results = load_cfi_results(output_dir)
    if cfi_results:
        cfi_summary = generate_cfi_summary(cfi_results)
        output['cfi'] = {
            'raw': cfi_results,
            'summary': cfi_summary,
        }
        metric_count = sum(len(v) if isinstance(v, dict) else 1 for v in cfi_results.values())
        print(f"CFI Experiments: {metric_count} metrics loaded")
    else:
        missing.append('cfi')
        print("CFI Experiments: no results found")

    return output, missing


def print_summary(output: Dict[str, Any]) -> None:
    """Print formatted summary to console."""
    print("\n" + "=" * 70)
    print("CSG-TRANSFORMER EVALUATION SUMMARY")
    print("=" * 70)

    for category in ['graph_classification', 'node_classification']:
        if category not in output:
            continue
        print(f"\n--- {category.replace('_', ' ').title()} ---")
        summary = output[category].get('summary')
        if summary is not None and not summary.empty:
            display_cols = [c for c in ['dataset', 'CSG-Transformer', 'Best Baseline', 'Gap to Best']
                           if c in summary.columns]
            print(summary[display_cols].to_string(index=False))
        else:
            print("  (no results)")

    if 'ablation' in output:
        print("\n--- Ablation Study ---")
        df = output['ablation'].get('raw')
        if df is not None and not df.empty:
            print(df.to_string(index=False))
        else:
            print("  (no results)")

    if 'cfi' in output:
        print("\n--- CFI Experiments ---")
        cfi_summary = output['cfi'].get('summary', {})
        acc_df = cfi_summary.get('base_accuracy')
        if acc_df is not None and not acc_df.empty:
            print("\nBase Graph Accuracy:")
            print(acc_df.to_string(index=False))

        layer_df = cfi_summary.get('layer_ablation')
        if layer_df is not None and not layer_df.empty:
            print("\nLayer Ablation:")
            print(layer_df.to_string(index=False))

        comp_df = cfi_summary.get('component_ablation')
        if comp_df is not None and not comp_df.empty:
            print("\nComponent Ablation:")
            print(comp_df.to_string(index=False))


def save_aggregated_output(
    output: Dict[str, Any],
    save_dir: str,
) -> str:
    """Save aggregated results to CSV/JSON files."""
    os.makedirs(save_dir, exist_ok=True)

    # Graph classification summary
    if 'graph_classification' in output:
        summary = output['graph_classification'].get('summary')
        if summary is not None and not summary.empty:
            path = os.path.join(save_dir, 'graph_classification_summary.csv')
            summary.to_csv(path, index=False)
            print(f"  Saved: {path}")

        raw = output['graph_classification'].get('raw')
        if raw is not None and not raw.empty:
            path = os.path.join(save_dir, 'graph_classification_raw.csv')
            raw.to_csv(path, index=False)
            print(f"  Saved: {path}")

    # Node classification summary
    if 'node_classification' in output:
        summary = output['node_classification'].get('summary')
        if summary is not None and not summary.empty:
            path = os.path.join(save_dir, 'node_classification_summary.csv')
            summary.to_csv(path, index=False)
            print(f"  Saved: {path}")

        raw = output['node_classification'].get('raw')
        if raw is not None and not raw.empty:
            path = os.path.join(save_dir, 'node_classification_raw.csv')
            raw.to_csv(path, index=False)
            print(f"  Saved: {path}")

    # Ablation
    if 'ablation' in output:
        raw = output['ablation'].get('raw')
        if raw is not None and not raw.empty:
            path = os.path.join(save_dir, 'ablation_aggregated.csv')
            raw.to_csv(path, index=False)
            print(f"  Saved: {path}")

    # CFI summary
    if 'cfi' in output:
        cfi_summary = output['cfi'].get('summary', {})
        acc_df = cfi_summary.get('base_accuracy')
        if acc_df is not None and not acc_df.empty:
            path = os.path.join(save_dir, 'cfi_base_accuracy.csv')
            acc_df.to_csv(path, index=False)
            print(f"  Saved: {path}")

        layer_df = cfi_summary.get('layer_ablation')
        if layer_df is not None and not layer_df.empty:
            path = os.path.join(save_dir, 'cfi_layer_ablation.csv')
            layer_df.to_csv(path, index=False)
            print(f"  Saved: {path}")

        comp_df = cfi_summary.get('component_ablation')
        if comp_df is not None and not comp_df.empty:
            path = os.path.join(save_dir, 'cfi_component_ablation.csv')
            comp_df.to_csv(path, index=False)
            print(f"  Saved: {path}")

    # Full JSON dump
    json_safe = _make_json_safe(output)
    json_path = os.path.join(save_dir, 'full_aggregation.json')
    with open(json_path, 'w') as f:
        json.dump(json_safe, f, indent=2, default=str)
    print(f"  Saved: {json_path}")

    return save_dir


def _make_json_safe(output: Dict[str, Any]) -> Dict[str, Any]:
    """Convert DataFrames to JSON-serializable dicts."""
    safe = {}
    for key, value in output.items():
        if isinstance(value, dict):
            safe[key] = {}
            for subkey, subval in value.items():
                if isinstance(subval, pd.DataFrame):
                    safe[key][subkey] = json.loads(subval.to_json(orient='records'))
                else:
                    safe[key][subkey] = subval
        else:
            safe[key] = value
    return safe


# ============================================================================
# Main Entry Point
# ============================================================================

def main():
    parser = argparse.ArgumentParser(
        description='Aggregate and summarize CSG-Transformer experiment results'
    )
    parser.add_argument('--output-dir', type=str, default=_DEFAULT_OUTPUT_DIR,
                        help='Directory containing experiment outputs')
    parser.add_argument('--save-dir', type=str, default='',
                        help='Directory to save aggregated results (default: {output_dir}/aggregated)')
    parser.add_argument('--mode', type=str, default='all',
                        choices=['all', 'graph', 'node', 'ablation', 'cfi'],
                        help='Which results to aggregate')
    parser.add_argument('--print', action='store_true', default=True,
                        help='Print summary to console')
    parser.add_argument('--save', action='store_true', default=True,
                        help='Save aggregated results')
    args = parser.parse_args()

    output_dir = args.output_dir
    save_dir = args.save_dir or os.path.join(output_dir, 'aggregated')

    if args.mode != 'all':
        # TODO: support per-mode loading
        print(f"Mode {args.mode} selected (aggregating all for now)")

    output, missing = aggregate_all(output_dir)

    if args.print:
        print_summary(output)

    if args.save:
        save_aggregated_output(output, save_dir)

    # Report missing sections
    if missing:
        print(f"\nNote: no results found for: {', '.join(missing)}")
        print("Run experiments first with run_all_experiments.py or individual scripts.")


if __name__ == '__main__':
    main()
