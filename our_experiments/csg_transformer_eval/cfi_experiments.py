"""
CFI Graph Experiments for CSG-Transformer.

This script evaluates CSG-Transformer's ability to distinguish CFI graphs,
which are classic hard cases for graph neural networks. It implements the
full experimental protocol from the specification (Section 6.9).

Experiments:
1. CFI graph distinguishing accuracy on various base graphs
2. Layer number (L) ablation for CFI distinguishing
3. Ablation study of key components for CFI tasks

References:
- Cai, Fürer, Immerman (1992): An optimal lower bound on the number of
  variables for graph identification
- Specification Section 6.9: CFI Graph Experiments
"""

import argparse
import json
import os
import sys
from datetime import datetime
from typing import Dict, List, Tuple

import networkx as nx
import numpy as np
import pandas as pd
import torch

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))

# Default output directory: our_experiments/csg_transformer_eval/outputs/cfi/
_CFI_DIR = os.path.dirname(os.path.abspath(__file__))
_DEFAULT_CFI_DIR = os.path.join(_CFI_DIR, 'outputs', 'cfi')

from cyclic_schema.csg_transformer import CSGTransformer, build_model_from_config
from our_experiments.csg_transformer_eval.utils import set_seed, get_device


# ============================================================================
# CFI Graph Construction
# ============================================================================

def build_cfi_graph(
    base_edges: List[Tuple[int, int]],
    twist_edges: set,
) -> nx.Graph:
    """
    Build a CFI graph from a base graph with optional twist.
    
    Standard CFI construction: for each edge e=(u,v) in the base graph, 
    create 6 nodes (x_e1,x_e2,x_e3,y_e1,y_e2,y_e3). Node parity is
    encoded by the twist set.
    
    Args:
        base_edges: List of edges in the base graph [(u,v), ...]
        twist_edges: Set of edges to be twisted
        
    Returns:
        CFI graph
    """
    G = nx.Graph()
    vertices = sorted(set(v for e in base_edges for v in e))
    
    # Create vertex nodes (2 per base vertex)
    for v in vertices:
        G.add_node(f"v{v}_0")
        G.add_node(f"v{v}_1")
    
    # Create edge nodes (6 per base edge) and their connections
    for (u, v) in base_edges:
        is_twisted = (u, v) in twist_edges or (v, u) in twist_edges
        
        # x nodes (even parity)
        G.add_node(f"e{u}_{v}_x1")
        G.add_node(f"e{u}_{v}_x2")
        G.add_node(f"e{u}_{v}_x3")
        
        # y nodes (odd parity)
        G.add_node(f"e{u}_{v}_y1")
        G.add_node(f"e{u}_{v}_y2")
        G.add_node(f"e{u}_{v}_y3")
        
        # Connect x_i to y_i (vertex constraint)
        for i in range(1, 4):
            G.add_edge(f"e{u}_{v}_x{i}", f"e{u}_{v}_y{i}")
        
        if not is_twisted:
            # x nodes connect to v0 and y nodes connect to v1
            G.add_edge(f"e{u}_{v}_x1", f"v{u}_0")
            G.add_edge(f"e{u}_{v}_x1", f"v{v}_0")
            G.add_edge(f"e{u}_{v}_x2", f"v{u}_0")
            G.add_edge(f"e{u}_{v}_x2", f"v{v}_0")
            G.add_edge(f"e{u}_{v}_x3", f"v{u}_0")
            G.add_edge(f"e{u}_{v}_x3", f"v{v}_0")
            
            G.add_edge(f"e{u}_{v}_y1", f"v{u}_1")
            G.add_edge(f"e{u}_{v}_y1", f"v{v}_1")
            G.add_edge(f"e{u}_{v}_y2", f"v{u}_1")
            G.add_edge(f"e{u}_{v}_y2", f"v{v}_1")
            G.add_edge(f"e{u}_{v}_y3", f"v{u}_1")
            G.add_edge(f"e{u}_{v}_y3", f"v{v}_1")
        else:
            # Twisted: x nodes connect to mixed parity
            G.add_edge(f"e{u}_{v}_x1", f"v{u}_0")
            G.add_edge(f"e{u}_{v}_x1", f"v{v}_1")
            G.add_edge(f"e{u}_{v}_x2", f"v{u}_0")
            G.add_edge(f"e{u}_{v}_x2", f"v{v}_1")
            G.add_edge(f"e{u}_{v}_x3", f"v{u}_0")
            G.add_edge(f"e{u}_{v}_x3", f"v{v}_1")
            
            G.add_edge(f"e{u}_{v}_y1", f"v{u}_1")
            G.add_edge(f"e{u}_{v}_y1", f"v{v}_0")
            G.add_edge(f"e{u}_{v}_y2", f"v{u}_1")
            G.add_edge(f"e{u}_{v}_y2", f"v{v}_0")
            G.add_edge(f"e{u}_{v}_y3", f"v{u}_1")
            G.add_edge(f"e{u}_{v}_y3", f"v{v}_0")
    
    return G


def build_cfi_pair(
    base_graph_type: str,
    twist_size: int = -1,
    n: int = 6,
) -> Tuple[nx.Graph, nx.Graph]:
    """
    Build a pair of non-isomorphic CFI graphs.
    
    Args:
        base_graph_type: Type of base graph. Options:
            - 'k4': Complete graph K4 (6 edges)
            - 'petersen': Petersen graph (15 edges)
            - 'dodecahedron': Dodecahedral graph (30 edges)
            - 'gp8_3': Generalized Petersen GP(8,3) (16 edges)
            - 'cycle': Cycle graph C_n (n edges)
            - 'complete': Complete graph K_n
            - 'regular_3': Random 3-regular graph with 2n nodes
        twist_size: Number of edges to twist (-1 for default half)
        n: Size parameter for cycle/complete/regular graphs
        
    Returns:
        (G1, G2): A pair of non-isomorphic CFI graphs
    """
    base_graphs = {
        'k4': [(1, 2), (1, 3), (1, 4), (2, 3), (2, 4), (3, 4)],
        'petersen': list(nx.petersen_graph().edges()),
        'dodecahedron': list(nx.dodecahedral_graph().edges()),
    }
    
    if base_graph_type == 'gp8_3':
        # Manually construct Generalized Petersen GP(8,3)
        gp_edges = []
        for i in range(8):
            gp_edges.append((i, (i + 1) % 8))        # outer cycle
            gp_edges.append((i + 8, ((i + 3) % 8) + 8))  # inner star
            gp_edges.append((i, i + 8))                # spokes
        base_edges = gp_edges
    elif base_graph_type == 'cycle':
        G_base = nx.cycle_graph(n)
        base_edges = list(G_base.edges())
    elif base_graph_type == 'complete':
        G_base = nx.complete_graph(n)
        base_edges = list(G_base.edges())
    elif base_graph_type == 'regular_3':
        G_base = nx.random_regular_graph(3, max(4, 2 * n), seed=42)
        base_edges = list(G_base.edges())
    elif base_graph_type in base_graphs:
        base_edges = base_graphs[base_graph_type]
    else:
        raise ValueError(
            f"Unknown base graph type: {base_graph_type}. "
            f"Options: k4, petersen, dodecahedron, gp8_3, cycle, complete, regular_3"
        )
    
    if twist_size < 0:
        twist_size = len(base_edges) // 2
    
    twist_edges = set(base_edges[:twist_size])
    
    G1 = build_cfi_graph(base_edges, set())
    G2 = build_cfi_graph(base_edges, twist_edges)
    
    return G1, G2


# ============================================================================
# Training and Evaluation
# ============================================================================

def _cfi_features(G: nx.Graph, device: torch.device) -> torch.Tensor:
    """Compute degree-based features for CFI graphs."""
    deg = torch.tensor([G.degree(v) for v in G.nodes()], dtype=torch.float32, device=device)
    max_deg = deg.max().item() if deg.numel() > 0 and deg.max().item() > 0 else 1.0
    norm_deg = deg / max_deg
    return torch.stack([deg, norm_deg], dim=1)


def train_cfi_distinguisher(
    model: CSGTransformer,
    cfi_pairs: List[Tuple[nx.Graph, nx.Graph]],
    device: torch.device,
    config: Dict,
    verbose: bool = True,
) -> CSGTransformer:
    """
    Train model to distinguish CFI graph pairs.
    
    Args:
        model: CSG-Transformer model
        cfi_pairs: List of (G1, G2) CFI graph pairs
        device: Device
        config: Training configuration
        verbose: Print progress
        
    Returns:
        Trained model
    """
    optimizer = torch.optim.AdamW(
        model.parameters(),
        lr=config.get('lr', 1e-3),
        weight_decay=config.get('weight_decay', 1e-4),
    )
    
    num_epochs = config.get('epochs', 100)
    
    for epoch in range(num_epochs):
        model.train()
        total_loss = 0.0
        
        for G1, G2 in cfi_pairs:
            features1 = _cfi_features(G1, device)
            features2 = _cfi_features(G2, device)
            
            emb1, logits1 = model(G1, features1)
            emb2, logits2 = model(G2, features2)
            
            if logits1.dim() == 1:
                logits1 = logits1.unsqueeze(0)
            if logits2.dim() == 1:
                logits2 = logits2.unsqueeze(0)
            
            logits = torch.cat([logits1, logits2], dim=0)
            labels = torch.tensor([0, 1], device=device)
            
            loss = torch.nn.functional.cross_entropy(logits, labels)
            
            optimizer.zero_grad()
            loss.backward()
            optimizer.step()
            
            total_loss += loss.item()
        
        if verbose and (epoch + 1) % 10 == 0:
            print(f"Epoch {epoch+1}/{num_epochs}: Loss={total_loss/len(cfi_pairs):.4f}")
    
    return model


def evaluate_cfi_distinguishing(
    model: CSGTransformer,
    cfi_pairs: List[Tuple[nx.Graph, nx.Graph]],
    device: torch.device,
) -> Dict:
    """
    Evaluate model's ability to distinguish CFI graphs.
    
    Returns accuracy, where a correct prediction means the model assigns
    different labels to the two graphs in a pair. Also tracks label
    consistency across pairs (model should always assign 0 to original and
    1 to twisted, not just different labels).
    
    Args:
        model: Trained model
        cfi_pairs: List of (G1, G2) CFI graph pairs
        device: Device
        
    Returns:
        Dict with accuracy, consistency, and detailed metrics
    """
    model.eval()
    correct = 0
    total = len(cfi_pairs)
    label_assignments = []
    
    with torch.no_grad():
        for G1, G2 in cfi_pairs:
            features1 = _cfi_features(G1, device)
            features2 = _cfi_features(G2, device)
            
            _, logits1 = model(G1, features1)
            _, logits2 = model(G2, features2)
            
            pred1 = logits1.argmax().item()
            pred2 = logits2.argmax().item()
            label_assignments.append((pred1, pred2))
            
            if pred1 != pred2:
                correct += 1
    
    accuracy = correct / total if total > 0 else 0.0
    
    # Consistency: model should assign G1→0, G2→1 (or opposite) consistently
    assignments = np.array(label_assignments)
    consistent_pattern = np.sum((assignments[:, 0] == 0) & (assignments[:, 1] == 1))
    opposite_pattern = np.sum((assignments[:, 0] == 1) & (assignments[:, 1] == 0))
    consistency = max(consistent_pattern, opposite_pattern) / total if total > 0 else 0.0
    
    return {
        'accuracy': accuracy,
        'consistency': consistency,
        'correct': correct,
        'total': total,
        'consistent_pattern': int(consistent_pattern),
        'opposite_pattern': int(opposite_pattern),
    }


# ============================================================================
# Experiment Runners
# ============================================================================

def run_cfi_experiment(
    base_graph_type: str = 'k4',
    num_pairs: int = 100,
    config: Dict = None,
    device: torch.device = None,
    save_dir: str = '',
) -> Dict:
    """
    Run CFI graph distinguishing experiment (Experiment 1 from spec §6.9.3).
    
    Args:
        base_graph_type: Type of base graph
        num_pairs: Number of CFI pairs
        config: Model configuration
        device: Device
        save_dir: Output directory
        
    Returns:
        Experiment results dict
    """
    if not save_dir:
        save_dir = _DEFAULT_CFI_DIR
    
    if config is None:
        config = {
            'hidden_dim': 128,
            'num_heads': 4,
            'dropout': 0.1,
            'L': 3,
            'T': 3,
            'I': 5,
            'lr': 1e-3,
            'weight_decay': 1e-4,
            'epochs': 100,
            'in_dim': 2,
        }
    
    if device is None:
        device = get_device()
    
    set_seed(42)
    
    print(f"\n{'='*60}")
    print(f"CFI Experiment: Base={base_graph_type}, Pairs={num_pairs}")
    print(f"{'='*60}")
    
    cfi_pairs = []
    for _ in range(num_pairs):
        G1, G2 = build_cfi_pair(base_graph_type)
        cfi_pairs.append((G1, G2))
    
    G_sample, _ = cfi_pairs[0]
    print(f"CFI graph size: {G_sample.number_of_nodes()} nodes, "
          f"{G_sample.number_of_edges()} edges")
    
    model = build_model_from_config(config)
    model = model.to(device)
    
    model = train_cfi_distinguisher(model, cfi_pairs, device, config)
    
    results = evaluate_cfi_distinguishing(model, cfi_pairs, device)
    
    results['base_graph'] = base_graph_type
    results['num_pairs'] = num_pairs
    results['config'] = str(config)
    results['timestamp'] = datetime.now().isoformat()
    
    os.makedirs(save_dir, exist_ok=True)
    filepath = os.path.join(save_dir, f"cfi_{base_graph_type}_results.json")
    with open(filepath, 'w') as f:
        json.dump(results, f, indent=2)
    
    print(f"\nAccuracy: {results['accuracy']:.4f} "
          f"({results['correct']}/{results['total']})")
    print(f"Results: {filepath}")
    
    return results


def run_cfi_layer_ablation(
    base_graph_type: str = 'k4',
    num_pairs: int = 100,
    config: Dict = None,
    device: torch.device = None,
    save_dir: str = '',
) -> pd.DataFrame:
    """
    Run CFI layer ablation experiment (Experiment 3 from spec §6.9.3).
    
    Tests L = 1, 2, 3, 4 to measure impact of CSG depth on CFI distinguishing.
    
    Args:
        base_graph_type: Type of base graph
        num_pairs: Number of CFI pairs
        config: Base configuration
        device: Device
        save_dir: Output directory
        
    Returns:
        DataFrame with results for each L value
    """
    if not save_dir:
        save_dir = _DEFAULT_CFI_DIR
    
    if config is None:
        config = {
            'hidden_dim': 128,
            'num_heads': 4,
            'dropout': 0.1,
            'T': 3,
            'I': 5,
            'lr': 1e-3,
            'weight_decay': 1e-4,
            'epochs': 100,
            'in_dim': 2,
        }
    
    if device is None:
        device = get_device()
    
    set_seed(42)
    
    print(f"\n{'='*60}")
    print(f"CFI Layer Ablation: Base={base_graph_type}")
    print(f"{'='*60}")
    
    cfi_pairs = []
    for _ in range(num_pairs):
        G1, G2 = build_cfi_pair(base_graph_type)
        cfi_pairs.append((G1, G2))
    
    results = []
    L_values = [1, 2, 3, 4]
    
    for L in L_values:
        cfg = config.copy()
        cfg['L'] = L
        
        model = build_model_from_config(cfg)
        model = model.to(device)
        
        import time
        t0 = time.time()
        model = train_cfi_distinguisher(model, cfi_pairs, device, cfg, verbose=False)
        train_time = time.time() - t0
        
        eval_results = evaluate_cfi_distinguishing(model, cfi_pairs, device)
        
        results.append({
            'L': L,
            'accuracy': eval_results['accuracy'],
            'consistency': eval_results['consistency'],
            'correct': eval_results['correct'],
            'total': eval_results['total'],
            'train_time': train_time,
        })
        
        print(f"L={L}: Accuracy={eval_results['accuracy']:.4f}, "
              f"Time={train_time:.1f}s")
    
    df = pd.DataFrame(results)
    
    os.makedirs(save_dir, exist_ok=True)
    filepath = os.path.join(save_dir, "cfi_layer_ablation.csv")
    df.to_csv(filepath, index=False)
    print(f"\nLayer ablation results: {filepath}")
    
    return df


def run_cfi_ablation_study(
    base_graph_type: str = 'k4',
    num_pairs: int = 100,
    config: Dict = None,
    device: torch.device = None,
    save_dir: str = '',
) -> pd.DataFrame:
    """
    Run comprehensive CFI ablation study (per spec §6.9.4).
    
    Tests various component ablations on CFI distinguishing accuracy.
    
    Args:
        base_graph_type: Type of base graph
        num_pairs: Number of CFI pairs
        config: Base configuration
        device: Device
        save_dir: Output directory
        
    Returns:
        DataFrame with ablation results
    """
    if not save_dir:
        save_dir = _DEFAULT_CFI_DIR
    if config is None:
        config = {
            'hidden_dim': 128,
            'num_heads': 4,
            'dropout': 0.1,
            'L': 3,
            'T': 3,
            'I': 5,
            'lr': 1e-3,
            'weight_decay': 1e-4,
            'epochs': 100,
            'in_dim': 2,
        }
    
    if device is None:
        device = get_device()
    
    set_seed(42)
    
    print(f"\n{'='*60}")
    print(f"CFI Ablation Study: Base={base_graph_type}")
    print(f"{'='*60}")
    
    cfi_pairs = []
    for _ in range(num_pairs):
        G1, G2 = build_cfi_pair(base_graph_type)
        cfi_pairs.append((G1, G2))
    
    ablation_configs = {
        'Full Model': config,
        'w/o Backward Attention': {**config, 'use_backward': False},
        'w/o TNA-Attention': {**config, 'use_tna': False},
        'w/o PE': {**config, 'pe_type': 'none', 'use_pe': False},
        'Single Layer': {**config, 'L': 1},
        'w/o Multi-layer CSG': {**config, 'L': 1},
        'Random PE': {**config, 'pe_type': 'lpe'},
    }
    
    results = []
    
    for variant_name, cfg in ablation_configs.items():
        model = build_model_from_config(cfg)
        model = model.to(device)
        
        model = train_cfi_distinguisher(
            model, cfi_pairs, device, cfg, verbose=False
        )
        
        eval_results = evaluate_cfi_distinguishing(model, cfi_pairs, device)
        
        results.append({
            'variant': variant_name,
            'accuracy': eval_results['accuracy'],
            'consistency': eval_results['consistency'],
            'correct': eval_results['correct'],
            'total': eval_results['total'],
        })
        
        print(f"{variant_name}: Accuracy={eval_results['accuracy']:.4f}")
    
    df = pd.DataFrame(results)
    
    os.makedirs(save_dir, exist_ok=True)
    filepath = os.path.join(save_dir, "cfi_ablation_study.csv")
    df.to_csv(filepath, index=False)
    print(f"\nAblation results: {filepath}")
    
    return df


def run_all_cfi_experiments(
    config: Dict = None,
    device: torch.device = None,
    save_dir: str = '',
) -> Dict:
    """
    Run all CFI experiments from the specification.
    
    Includes:
    - Experiment 1: CFI distinguishing accuracy on multiple base graphs
    - Experiment 2: Layer ablation
    - Experiment 3: Component ablation
    
    Args:
        config: Base configuration
        device: Device
        save_dir: Output directory
        
    Returns:
        Dict with all results
    """
    if not save_dir:
        save_dir = _DEFAULT_CFI_DIR
    print("\n" + "=" * 60)
    print("RUNNING ALL CFI EXPERIMENTS")
    print("=" * 60)
    
    all_results = {}
    
    # Experiment 1: Multiple base graphs
    print("\n--- Experiment 1: CFI Distinguishing by Base Graph ---")
    base_graphs = ['k4', 'petersen', 'dodecahedron', 'gp8_3',
                   'cycle', 'complete', 'regular_3']
    for base in base_graphs:
        result = run_cfi_experiment(base, 100, config, device, save_dir)
        all_results[f'cfi_{base}'] = result['accuracy']
    
    # Experiment 2: Layer ablation
    print("\n--- Experiment 2: Layer Ablation ---")
    layer_df = run_cfi_layer_ablation('petersen', 100, config, device, save_dir)
    all_results['layer_ablation'] = layer_df.to_dict()
    
    # Experiment 3: Component ablation
    print("\n--- Experiment 3: Component Ablation ---")
    abl_df = run_cfi_ablation_study('petersen', 100, config, device, save_dir)
    all_results['component_ablation'] = abl_df.to_dict()
    
    summary_file = os.path.join(save_dir, "all_cfi_results.json")
    with open(summary_file, 'w') as f:
        json.dump(all_results, f, indent=2, default=str)
    
    print(f"\nAll CFI results saved to {summary_file}")
    
    return all_results


# ============================================================================
# Main Entry Point
# ============================================================================

def main():
    """Main entry point for CFI experiments."""
    parser = argparse.ArgumentParser(description='CFI Graph Experiments')
    
    parser.add_argument('--experiment', type=str, default='distinguish',
                        choices=['distinguish', 'layer_ablation',
                                 'ablation_study', 'all'],
                        help='CFI experiment type')
    parser.add_argument('--base_graph', type=str, default='petersen',
                        choices=['k4', 'petersen', 'dodecahedron', 'gp8_3',
                                 'cycle', 'complete', 'regular_3'],
                        help='Base graph for CFI construction')
    parser.add_argument('--num_pairs', type=int, default=100,
                        help='Number of CFI graph pairs')
    parser.add_argument('--hidden_dim', type=int, default=128)
    parser.add_argument('--L', type=int, default=3,
                        help='CSG layers')
    parser.add_argument('--T', type=int, default=3,
                        help='TNA rounds')
    parser.add_argument('--I', type=int, default=5,
                        help='Global iterations')
    parser.add_argument('--epochs', type=int, default=100,
                        help='Training epochs')
    parser.add_argument('--seed', type=int, default=42)
    parser.add_argument('--save_dir', type=str, default='')
    
    args = parser.parse_args()
    
    config = {
        'hidden_dim': args.hidden_dim,
        'num_heads': 4,
        'dropout': 0.1,
        'L': args.L,
        'T': args.T,
        'I': args.I,
        'epochs': args.epochs,
        'lr': 1e-3,
        'weight_decay': 1e-4,
        'in_dim': 2,
    }
    
    device = get_device()
    set_seed(args.seed)
    if not args.save_dir:
        args.save_dir = _DEFAULT_CFI_DIR
    
    if args.experiment == 'distinguish':
        run_cfi_experiment(args.base_graph, args.num_pairs, config,
                           device, args.save_dir)
    elif args.experiment == 'layer_ablation':
        run_cfi_layer_ablation(args.base_graph, args.num_pairs, config,
                               device, args.save_dir)
    elif args.experiment == 'ablation_study':
        run_cfi_ablation_study(args.base_graph, args.num_pairs, config,
                               device, args.save_dir)
    elif args.experiment == 'all':
        run_all_cfi_experiments(config, device, args.save_dir)


if __name__ == '__main__':
    main()
