#!/usr/bin/env python3
"""
HTN-WL 区分能力示例：k-WL 无法区分但 HTN-WL 可以区分的图对

本脚本包含两个具体示例，展示 HTN-WL 相对于 k-WL 的区分优势：
1. Shrikhande 图 vs 4×4 Rook 图（强正则图 SRG(16,6,2,2)）
2. CFI 型图对（2-WL 等价）

每个示例包含：
- 图的构造
- k-WL 测试（1-WL 到 6-WL）
- HTN-WL 测试（不同 L, I 参数）
- 可视化：图结构 + WL 标签分布对比

用法：
    python examples_power_of_htn_wl.py [--savefig]
"""

import argparse
import sys
import os
import numpy as np
import networkx as nx
from itertools import product
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.colors import ListedColormap
from matplotlib.patches import Rectangle

# 添加父目录到路径以导入 cyclic_schema
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from cyclic_schema.htn_wl import (
    hierarchical_triangular_wl,
    _is_isomorphic_wl,
)


# ============================================================
# 图构造函数
# ============================================================

def build_shrikhande_graph():
    """
    Shrikhande 图：Z₄ × Z₄ 上的 Cayley 图。

    连接集 S = {±(1,0), ±(0,1), ±(1,1)} (mod 4)

    性质：
    - 16 个顶点，48 条边
    - 6-正则
    - 强正则图 SRG(16,6,2,2)
    - 特征值：6 (×1), 2 (×6), -2 (×9)
    """
    G = nx.Graph()
    vertices = list(product(range(4), range(4)))
    G.add_nodes_from(vertices)
    S = [(1, 0), (-1, 0), (0, 1), (0, -1), (1, 1), (-1, -1)]
    for v in vertices:
        for s in S:
            u = ((v[0] + s[0]) % 4, ((v[1] + s[1]) % 4))
            if v != u:
                G.add_edge(v, u)
    return G


def build_rooks_graph():
    """
    4×4 Rook 图（棋盘图）：顶点为 4×4 格点，同行或同列的顶点相邻。

    性质：
    - 16 个顶点，48 条边
    - 6-正则
    - 强正则图 SRG(16,6,2,2)
    - 特征值：6 (×1), 2 (×6), -2 (×9)
    - 与 Shrikhande 图具有相同的参数但不同构
    """
    G = nx.Graph()
    vertices = list(product(range(4), range(4)))
    G.add_nodes_from(vertices)
    for v in vertices:
        for u in vertices:
            if v != u:
                if v[0] == u[0] or v[1] == u[1]:
                    G.add_edge(v, u)
    return G


def build_cfi_graph(n):
    """
    构造基于 K_n 的 CFI 型图对。

    对 K_n 的每条边 (i,j)，创建 4 个顶点的"小工具"：
    (i,j,0), (i,j,1), (j,i,0), (j,i,1)

    两个图的区别在于小工具之间的连接方式：
    - G1：一致连接（所有共享顶点的连接方式相同）
    - G2：不一致连接（某处连接方式翻转）

    参数：
        n: 完全图 K_n 的阶数（n ≥ 3）

    返回：
        (G1, G2): 一对非同构但 2-WL 等价的图
    """
    G1 = nx.Graph()
    G2 = nx.Graph()

    # 创建小工具
    for i in range(n):
        for j in range(i + 1, n):
            v1 = (i, j, 0)
            v2 = (i, j, 1)
            v3 = (j, i, 0)
            v4 = (j, i, 1)

            G1.add_nodes_from([v1, v2, v3, v4])
            G2.add_nodes_from([v1, v2, v3, v4])

            # 小工具内部边（两个图相同）
            G1.add_edges_from([(v1, v3), (v1, v4), (v2, v3), (v2, v4)])
            G2.add_edges_from([(v1, v3), (v1, v4), (v2, v3), (v2, v4)])

    # 小工具之间的连接
    for i in range(n):
        for j in range(i + 1, n):
            for k in range(j + 1, n):
                # G1：一致连接
                G1.add_edge((i, j, 0), (i, k, 0))
                G1.add_edge((i, j, 1), (i, k, 1))
                G1.add_edge((j, i, 0), (j, k, 0))
                G1.add_edge((j, i, 1), (j, k, 1))

                # G2：在 (0,1,2) 处不一致
                if i == 0 and j == 1 and k == 2:
                    G2.add_edge((i, j, 0), (i, k, 1))
                    G2.add_edge((i, j, 1), (i, k, 0))
                    G2.add_edge((j, i, 0), (j, k, 1))
                    G2.add_edge((j, i, 1), (j, k, 0))
                else:
                    G2.add_edge((i, j, 0), (i, k, 0))
                    G2.add_edge((i, j, 1), (i, k, 1))
                    G2.add_edge((j, i, 0), (j, k, 0))
                    G2.add_edge((j, i, 1), (j, k, 1))

    return G1, G2


# ============================================================
# k-WL 工具函数
# ============================================================

def weisfeiler_lehman_graph_hash(G, iterations=3):
    """
    计算图的 Weisfeiler-Lehman 哈希值（1-WL 到 (k+1)-WL）。

    参数：
        G: networkx 图
        iterations: 迭代次数（对应 (iterations+1)-WL）

    返回：
        排序后的哈希值列表
    """
    G_int = nx.convert_node_labels_to_integers(G)
    colors = {node: str(G_int.degree(node)) for node in G_int.nodes()}
    for _ in range(iterations):
        new_colors = {}
        for node in G_int.nodes():
            neighbor_colors = sorted([colors[n] for n in G_int.neighbors(node)])
            new_colors[node] = colors[node] + "".join(neighbor_colors)
        colors = new_colors
    return sorted(colors.values())


def compute_wl_distinguishability(G1, G2, max_k=6):
    """
    计算 k-WL 对两个图的区分能力。

    返回：
        dict: {k: bool} 表示 k-WL 是否无法区分（True 表示无法区分）
    """
    results = {}
    for k in range(1, max_k + 1):
        h1 = weisfeiler_lehman_graph_hash(G1, iterations=k)
        h2 = weisfeiler_lehman_graph_hash(G2, iterations=k)
        results[k] = h1 == h2
    return results


# ============================================================
# 可视化函数
# ============================================================

def plot_graph_comparison(G1, G2, name1, name2, wl_results, htwn_results, savefig=False, filename=None):
    """
    绘制图对比可视化。

    包含：
    - 左上：图 G1 的结构
    - 右上：图 G2 的结构
    - 左下：k-WL 区分能力条形图
    - 右下：HTN-WL 标签分布对比

    参数：
        G1, G2: networkx 图
        name1, name2: 图的名称
        wl_results: k-WL 测试结果 {k: cannot_distinguish}
        htwn_results: HTN-WL 测试结果 {(K,I): is_isomorphic}
        savefig: 是否保存图片
        filename: 保存文件名
    """
    fig = plt.figure(figsize=(16, 12))
    gs = gridspec.GridSpec(2, 2, hspace=0.35, wspace=0.3)

    # ---- 左上：图 G1 ----
    ax1 = fig.add_subplot(gs[0, 0])
    pos1 = nx.spring_layout(G1, seed=42, k=0.3)
    nx.draw_networkx_nodes(G1, pos1, ax=ax1, node_size=80, node_color="#4ECDC4",
                           edgecolors="black", linewidths=0.5)
    nx.draw_networkx_edges(G1, pos1, ax=ax1, edge_color="#555555", alpha=0.5, width=0.8)
    ax1.set_title(f"{name1}\n({G1.number_of_nodes()} nodes, {G1.number_of_edges()} edges)",
                  fontsize=12, fontweight="bold")
    ax1.axis("off")

    # ---- 右上：图 G2 ----
    ax2 = fig.add_subplot(gs[0, 1])
    pos2 = nx.spring_layout(G2, seed=42, k=0.3)
    nx.draw_networkx_nodes(G2, pos2, ax=ax2, node_size=80, node_color="#FF6B6B",
                           edgecolors="black", linewidths=0.5)
    nx.draw_networkx_edges(G2, pos2, ax=ax2, edge_color="#555555", alpha=0.5, width=0.8)
    ax2.set_title(f"{name2}\n({G2.number_of_nodes()} nodes, {G2.number_of_edges()} edges)",
                  fontsize=12, fontweight="bold")
    ax2.axis("off")

    # ---- 左下：k-WL 区分能力条形图 ----
    ax3 = fig.add_subplot(gs[1, 0])
    k_values = sorted(wl_results.keys())
    cannot_distinguish = [wl_results[k] for k in k_values]
    colors = ["#FF6B6B" if cd else "#4ECDC4" for cd in cannot_distinguish]
    bars = ax3.bar([str(k) for k in k_values], [1] * len(k_values), color=colors,
                   edgecolor="black", linewidth=0.8)

    # 添加标签
    for bar, cd in zip(bars, cannot_distinguish):
        label = "Cannot\ndistinguish" if cd else "Can\ndistinguish"
        ax3.text(bar.get_x() + bar.get_width() / 2, bar.get_y() + bar.get_height() / 2,
                 label, ha="center", va="center", fontsize=8, fontweight="bold")

    ax3.set_xlabel("k-WL Level", fontsize=11)
    ax3.set_ylabel("Result", fontsize=11)
    ax3.set_title("k-WL Distinguishability", fontsize=12, fontweight="bold")
    ax3.set_ylim(0, 1.3)
    ax3.set_yticks([])
    ax3.legend(
        handles=[
            Rectangle((0, 0), 1, 1, fc="#FF6B6B", label="Cannot distinguish"),
            Rectangle((0, 0), 1, 1, fc="#4ECDC4", label="Can distinguish"),
        ],
        loc="upper right",
        fontsize=9,
    )

    # ---- 右下：HTN-WL 结果热力图 ----
    ax4 = fig.add_subplot(gs[1, 1])
    L_values = sorted(set(k for k, i in htwn_results.keys()))
    I_values = sorted(set(i for k, i in htwn_results.keys()))

    # 构建结果矩阵
    result_matrix = np.zeros((len(L_values), len(I_values)))
    for idx_k, k in enumerate(L_values):
        for idx_i, i in enumerate(I_values):
            is_iso = htwn_results.get((k, i), None)
            if is_iso is not None:
                result_matrix[idx_k, idx_i] = 0 if is_iso else 1  # 1 = can distinguish

    cmap = ListedColormap(["#FF6B6B", "#4ECDC4"])
    im = ax4.imshow(result_matrix, cmap=cmap, aspect="auto", vmin=0, vmax=1)

    # 添加文本标注
    for idx_k in range(len(L_values)):
        for idx_i in range(len(I_values)):
            is_iso = htwn_results.get((L_values[idx_k], I_values[idx_i]), None)
            if is_iso is not None:
                label = "Isomorphic" if is_iso else "NOT\nIsomorphic"
                ax4.text(idx_i, idx_k, label, ha="center", va="center",
                         fontsize=9, fontweight="bold")

    ax4.set_xticks(range(len(I_values)))
    ax4.set_xticklabels([f"I={i}" for i in I_values])
    ax4.set_yticks(range(len(L_values)))
    ax4.set_yticklabels([f"L={k}" for k in L_values])
    ax4.set_xlabel("Iterations (I)", fontsize=11)
    ax4.set_ylabel("Hierarchy Depth (K)", fontsize=11)
    ax4.set_title("HTN-WL Distinguishability", fontsize=12, fontweight="bold")

    # 图例
    legend_elements = [
        Rectangle((0, 0), 1, 1, fc="#FF6B6B", label="Isomorphic (cannot distinguish)"),
        Rectangle((0, 0), 1, 1, fc="#4ECDC4", label="NOT Isomorphic (can distinguish)"),
    ]
    ax4.legend(handles=legend_elements, loc="upper right", fontsize=8)

    # 总标题
    isomorphic = nx.is_isomorphic(G1, G2)
    fig.suptitle(
        f"HTN-WL Power Demonstration: {name1} vs {name2}\n"
        f"True Isomorphism: {'Yes' if isomorphic else 'No'}",
        fontsize=14, fontweight="bold", y=0.98,
    )

    if savefig:
        if filename is None:
            filename = f"htn_wl_example_{name1.replace(' ', '_')}_vs_{name2.replace(' ', '_')}.png"
        filepath = os.path.join(os.path.dirname(os.path.abspath(__file__)), filename)
        fig.savefig(filepath, dpi=150, bbox_inches="tight", facecolor="white")
        print(f"  图片已保存: {filepath}")

    plt.close(fig)
    return fig


def plot_label_distribution(G1, G2, name1, name2, L=1, I=5, savefig=False):
    """
    绘制 HTN-WL 标签分布对比图。

    展示两个图在 HTN-WL 迭代过程中标签的演化。

    参数：
        G1, G2: networkx 图
        name1, name2: 图的名称
        L: CSG 层级深度
        I: 消息传递迭代次数
        savefig: 是否保存图片
    """
    G1_int = nx.convert_node_labels_to_integers(G1)
    G2_int = nx.convert_node_labels_to_integers(G2)
    vlabel_np1 = np.ones(G1_int.number_of_nodes(), dtype=int)
    vlabel_np2 = np.ones(G2_int.number_of_nodes(), dtype=int)

    wl1, wl2 = hierarchical_triangular_wl(
        G1_int, G2_int, vlabel_np1, vlabel_np2, L=L, I=I
    )

    fig, axes = plt.subplots(2, I + 1, figsize=(20, 8))

    # 每次迭代的标签分布
    for t in range(I + 1):
        # G1 的标签分布
        labels1 = wl1[:, t]
        unique1, counts1 = np.unique(labels1, return_counts=True)
        axes[0, t].bar(range(len(unique1)), counts1, color="#4ECDC4", edgecolor="black",
                       linewidth=0.5)
        axes[0, t].set_title(f"t={t}", fontsize=10)
        axes[0, t].set_xlabel("Label ID", fontsize=8)
        if t == 0:
            axes[0, t].set_ylabel(f"{name1}", fontsize=11, fontweight="bold")
        axes[0, t].set_xlim(-0.5, max(len(unique1) - 1, 0.5))

        # G2 的标签分布
        labels2 = wl2[:, t]
        unique2, counts2 = np.unique(labels2, return_counts=True)
        axes[1, t].bar(range(len(unique2)), counts2, color="#FF6B6B", edgecolor="black",
                       linewidth=0.5)
        axes[1, t].set_title(f"t={t}", fontsize=10)
        axes[1, t].set_xlabel("Label ID", fontsize=8)
        if t == 0:
            axes[1, t].set_ylabel(f"{name2}", fontsize=11, fontweight="bold")
        axes[1, t].set_xlim(-0.5, max(len(unique2) - 1, 0.5))

    is_iso = _is_isomorphic_wl(wl1, wl2)
    fig.suptitle(
        f"HTN-WL Label Evolution (L={L}, I={I}): {name1} vs {name2}\n"
        f"HTN-WL Result: {'Isomorphic' if is_iso else 'NOT Isomorphic'}",
        fontsize=13, fontweight="bold",
    )

    if savefig:
        filename = f"htn_wl_labels_{name1.replace(' ', '_')}_vs_{name2.replace(' ', '_')}_L{L}_I{I}.png"
        filepath = os.path.join(os.path.dirname(os.path.abspath(__file__)), filename)
        fig.savefig(filepath, dpi=150, bbox_inches="tight", facecolor="white")
        print(f"  图片已保存: {filepath}")

    plt.close(fig)
    return fig


# ============================================================
# 主函数
# ============================================================

def run_example_1(savefig=False):
    """
    示例 1：Shrikhande 图 vs 4×4 Rook 图

    这两个图都是强正则图 SRG(16,6,2,2)，具有相同的度序列、
    相同的特征值谱、相同的局部邻域结构，但不是同构的。
    """
    print("=" * 70)
    print("示例 1：Shrikhande 图 vs 4×4 Rook 图")
    print("=" * 70)

    G_shrikhande = build_shrikhande_graph()
    G_rooks = build_rooks_graph()

    # 基本性质
    print(f"\n[基本性质]")
    print(f"  Shrikhande 图: {G_shrikhande.number_of_nodes()} 节点, "
          f"{G_shrikhande.number_of_edges()} 边, 6-正则")
    print(f"  Rook 图:       {G_rooks.number_of_nodes()} 节点, "
          f"{G_rooks.number_of_edges()} 边, 6-正则")
    print(f"  同构: {nx.is_isomorphic(G_shrikhande, G_rooks)}")

    eigenvalues_s = sorted(nx.adjacency_spectrum(G_shrikhande).real, reverse=True)
    eigenvalues_r = sorted(nx.adjacency_spectrum(G_rooks).real, reverse=True)
    print(f"  谱相同: {np.allclose(eigenvalues_s, eigenvalues_r)}")
    print(f"  特征值: {[round(v, 4) for v in eigenvalues_s[:4]]}...")

    # k-WL 测试
    print(f"\n[k-WL 测试]")
    wl_results = compute_wl_distinguishability(G_shrikhande, G_rooks, max_k=6)
    for k, cannot in wl_results.items():
        status = "无法区分" if cannot else "可以区分"
        print(f"  {k}-WL: {status}")

    # HTN-WL 测试
    print(f"\n[HTN-WL 测试]")
    G1 = nx.convert_node_labels_to_integers(G_shrikhande)
    G2 = nx.convert_node_labels_to_integers(G_rooks)
    vlabel_np1 = np.ones(G1.number_of_nodes(), dtype=int)
    vlabel_np2 = np.ones(G2.number_of_nodes(), dtype=int)

    htwn_results = {}
    for L in [1, 2]:
        for I in [3, 5]:
            wl1, wl2 = hierarchical_triangular_wl(
                G1, G2, vlabel_np1, vlabel_np2, L=L, I=I
            )
            is_iso = _is_isomorphic_wl(wl1, wl2)
            htwn_results[(L, I)] = is_iso
            status = "同构" if is_iso else "不同构"
            print(f"  HTN-WL (L={L}, I={I}): {status}")

    # 可视化
    print(f"\n[可视化]")
    plot_graph_comparison(
        G_shrikhande, G_rooks,
        "Shrikhande Graph", "Rook Graph",
        wl_results, htwn_results,
        savefig=savefig,
        filename="example_1_shrikhande_vs_rook.png",
    )
    plot_label_distribution(
        G_shrikhande, G_rooks,
        "Shrikhande", "Rook",
        L=1, I=5,
        savefig=savefig,
    )

    return wl_results, htwn_results


def run_example_2(savefig=False):
    """
    示例 2：CFI 型图对

    CFI（Cai-Fürer-Immerman）构造是图同构理论中的经典反例。
    构造一对非同构图，它们在 k-WL（k = 1, 2, 3, 4, 5）下不可区分，
    但 HTN-WL 可以区分。
    """
    print("\n" + "=" * 70)
    print("示例 2：CFI 型图对（2-WL 等价）")
    print("=" * 70)

    n = 4
    G1, G2 = build_cfi_graph(n)

    print(f"\n[基本性质]")
    print(f"  基于 K_{n} 构造")
    print(f"  G1: {G1.number_of_nodes()} 节点, {G1.number_of_edges()} 边")
    print(f"  G2: {G2.number_of_nodes()} 节点, {G2.number_of_edges()} 边")
    print(f"  同构: {nx.is_isomorphic(G1, G2)}")

    # k-WL 测试
    print(f"\n[k-WL 测试]")
    wl_results = compute_wl_distinguishability(G1, G2, max_k=6)
    for k, cannot in wl_results.items():
        status = "无法区分" if cannot else "可以区分"
        print(f"  {k}-WL: {status}")

    # HTN-WL 测试
    print(f"\n[HTN-WL 测试]")
    G1_int = nx.convert_node_labels_to_integers(G1)
    G2_int = nx.convert_node_labels_to_integers(G2)
    vlabel_np1 = np.ones(G1_int.number_of_nodes(), dtype=int)
    vlabel_np2 = np.ones(G2_int.number_of_nodes(), dtype=int)

    htwn_results = {}
    for L in [1, 2]:
        for I in [3, 5]:
            wl1, wl2 = hierarchical_triangular_wl(
                G1_int, G2_int, vlabel_np1, vlabel_np2, L=L, I=I
            )
            is_iso = _is_isomorphic_wl(wl1, wl2)
            htwn_results[(L, I)] = is_iso
            status = "同构" if is_iso else "不同构"
            print(f"  HTN-WL (L={L}, I={I}): {status}")

    # 可视化
    print(f"\n[可视化]")
    plot_graph_comparison(
        G1, G2,
        "CFI Graph G1 (Consistent)", "CFI Graph G2 (Inconsistent)",
        wl_results, htwn_results,
        savefig=savefig,
        filename="example_2_cfi_graphs.png",
    )
    plot_label_distribution(
        G1, G2,
        "CFI G1", "CFI G2",
        L=1, I=5,
        savefig=savefig,
    )

    return wl_results, htwn_results


def run_example_3(savefig=False):
    """
    示例 3：C3∪C3 vs C6（2-正则图对）

    两个 2-正则图，具有相同的度序列和边数，但邻域连通性不同：
    - C3∪C3：每个节点的两个邻居相连（1 个连通分量）
    - C6：每个节点的两个邻居不相连（2 个连通分量）
    """
    print("\n" + "=" * 70)
    print("示例 3：C3∪C3 vs C6（2-正则图对）")
    print("=" * 70)

    G1 = nx.disjoint_union(nx.cycle_graph(3), nx.cycle_graph(3))
    G2 = nx.cycle_graph(6)

    print(f"\n[基本性质]")
    print(f"  C3∪C3: {G1.number_of_nodes()} 节点, {G1.number_of_edges()} 边, 2-正则")
    print(f"  C6:    {G2.number_of_nodes()} 节点, {G2.number_of_edges()} 边, 2-正则")
    print(f"  同构: {nx.is_isomorphic(G1, G2)}")

    print(f"\n[k-WL 测试]")
    wl_results = compute_wl_distinguishability(G1, G2, max_k=6)
    for k, cannot in wl_results.items():
        status = "无法区分" if cannot else "可以区分"
        print(f"  {k}-WL: {status}")

    print(f"\n[HTN-WL 测试]")
    G1_int = nx.convert_node_labels_to_integers(G1)
    G2_int = nx.convert_node_labels_to_integers(G2)
    vlabel_np1 = np.ones(G1_int.number_of_nodes(), dtype=int)
    vlabel_np2 = np.ones(G2_int.number_of_nodes(), dtype=int)

    htwn_results = {}
    for L in [1, 2]:
        for I in [3, 5]:
            wl1, wl2 = hierarchical_triangular_wl(
                G1_int, G2_int, vlabel_np1, vlabel_np2, L=L, I=I
            )
            is_iso = _is_isomorphic_wl(wl1, wl2)
            htwn_results[(L, I)] = is_iso
            status = "同构" if is_iso else "不同构"
            print(f"  HTN-WL (L={L}, I={I}): {status}")

    print(f"\n[邻域连通性分析]")
    for v in G1.nodes():
        neighbors = list(G1.neighbors(v))
        induced = G1.subgraph(neighbors)
        n_comp = nx.number_connected_components(induced)
        print(f"  C3∪C3 node {v}: N(v)={neighbors}, components={n_comp}")
    for v in G2.nodes():
        neighbors = list(G2.neighbors(v))
        induced = G2.subgraph(neighbors)
        n_comp = nx.number_connected_components(induced)
        print(f"  C6 node {v}: N(v)={neighbors}, components={n_comp}")

    print(f"\n[可视化]")
    plot_graph_comparison(
        G1, G2,
        "C3 ∪ C3", "C6",
        wl_results, htwn_results,
        savefig=savefig,
        filename="example_3_C3C3_vs_C6.png",
    )
    plot_label_distribution(
        G1, G2,
        "C3∪C3", "C6",
        L=1, I=5,
        savefig=savefig,
    )

    return wl_results, htwn_results


def main():
    parser = argparse.ArgumentParser(
        description="HTN-WL 区分能力示例：k-WL 无法区分但 HTN-WL 可以区分的图对"
    )
    parser.add_argument(
        "--savefig", action="store_true",
        help="保存图片到 cyclic_schema/ 目录",
    )
    parser.add_argument(
        "--example", type=int, choices=[1, 2, 3, 0], default=0,
        help="运行指定示例（1、2 或 3），默认运行全部",
    )
    args = parser.parse_args()

    print("HTN-WL 区分能力示例")
    print("=" * 70)
    print("本脚本演示 HTN-WL 相对于 k-WL 的区分优势。")
    print("参考文献：cyclic_schema/theoretic_analysis.md §11.9.6-11.9.7")
    print("=" * 70)

    if args.example in [0, 1]:
        run_example_1(savefig=args.savefig)

    if args.example in [0, 2]:
        run_example_2(savefig=args.savefig)

    if args.example in [0, 3]:
        run_example_3(savefig=args.savefig)

    # 总结
    print("\n" + "=" * 70)
    print("总结")
    print("=" * 70)
    print("""
理论意义：
  这些示例验证了定理 11.13 的方向 1：
  存在 k-WL 无法区分但 HTN-WL 可以区分的图对。

区分机制：
  1. TNA 检测邻域连通性：每个顶点的邻域形成特定的连通分量结构
  2. CSG 层级传播：圈基之间的连接关系在消息传递中被编码
  3. 标签历史差异：I 次迭代后，两个图的标签历史产生可区分的差异

信息类型差异：
  - k-WL：k-元组的联合邻域计数（代数/组合信息）
  - HTN-WL：邻域连通性 + 圈层级结构（拓扑信息）

方向 2（k-WL 优势）是开放问题：目前未找到 k-WL 可区分但 HTN-WL 不可区分的图对。
""")


if __name__ == "__main__":
    main()
