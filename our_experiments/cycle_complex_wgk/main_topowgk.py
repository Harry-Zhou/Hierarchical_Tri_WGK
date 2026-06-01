from itertools import product
from multiprocessing import Pool, cpu_count
import numpy as np
import pandas as pd
import os
import pickle
from sklearn.preprocessing import RobustScaler, StandardScaler, MinMaxScaler
from sklearn.svm import NuSVC
from sklearn.metrics import accuracy_score, f1_score, roc_auc_score
import csv
from sklearn.model_selection import train_test_split, GridSearchCV, PredefinedSplit

import sys
cur_dir = os.path.dirname(__file__)
parent_dir = os.path.dirname(cur_dir)
parent_dir = os.path.dirname(parent_dir)
sys.path.append(parent_dir)
from utils import get_start_end_indices, get_recommended_nproc
from topo_wasserstein_graph_kernel import TopoWassersteinGraphKernel
from process_dataset import download_dataset, construct_dataset

def gauss_kernel(dist_np, gk_gamma):
    sim_np = np.exp(-1.0 * gk_gamma * dist_np)
    sim_np = sim_np * (sim_np >= 1e-8).astype(np.int32)
    # Removed scaler to keep kernel values in [0,1] range for SVM compatibility
    return sim_np

def train_and_eval_worker(
    dname, nu, tol, max_iter, alpha_list, gk_gamma, 
    ot_sim_np, wl_sim_np, y, run_time, random_state = 42
):
    result_list = []
    
    indices = np.arange(len(y))
    # 第一次 split：train_val (85%) 和 test (15%)
    train_val_idx, test_idx = train_test_split(
        indices, test_size=0.15, stratify=y, random_state=random_state
    )
    # 第二次 split：train (70%) 和 val (15%)，基于 train_val
    train_idx, val_idx = train_test_split(
        train_val_idx, test_size=0.1765, stratify=y[train_val_idx], random_state=random_state
    )  # 15/85 ≈ 0.1765
    
    y_test = y[test_idx]
    y_train_val = y[train_val_idx]
    
    for alpha in alpha_list:
        X = alpha * ot_sim_np + (1 - alpha) * wl_sim_np
        # 构建 kernel 矩阵（恢复为对切片做 sanitize + max-abs 标准化）
        X_train_val = X[np.ix_(train_val_idx, train_val_idx)]
        X_test = X[np.ix_(test_idx, train_val_idx)]
        # sanitize numeric issues: replace nan/inf, ensure float64
        X_train_val = np.array(np.nan_to_num(X_train_val, nan=0.0, posinf=1e12, neginf=-1e12), dtype=np.float64)
        X_test = np.array(np.nan_to_num(X_test, nan=0.0, posinf=1e12, neginf=-1e12), dtype=np.float64)
        # normalize to avoid extremely large values causing solver instability
        max_abs = np.max(np.abs(X_train_val))
        if np.isfinite(max_abs) and max_abs > 0:
            X_train_val = X_train_val / max_abs
            X_test = X_test / max_abs
        # add tiny jitter on diagonal for numerical stability
        X_train_val = X_train_val + np.eye(X_train_val.shape[0]) * 1e-8
        
        # 超参数选择：使用 GridSearchCV 在 train_val 上搜索 svm_gamma
        # 定义 PredefinedSplit：train 部分为 fold 0，val 部分为 -1（外部验证）
        test_fold = np.full(len(train_val_idx), 0)
        test_fold[len(train_idx):] = -1  # val 部分标记为 -1
        ps = PredefinedSplit(test_fold)
        
        num_cls = y_train_val.max() + 1
        clf = NuSVC(
            kernel='precomputed', nu=nu, tol=tol, probability=True, 
            max_iter=max_iter, random_state=random_state, 
            decision_function_shape='ovo' if num_cls > 2 else 'ovr'
        )
        
        # 1. 计算每个特征的方差及其平均值
        var_per_feature = np.var(X, axis=0, ddof=0)
        mean_var = np.mean(var_per_feature)
        if not np.isfinite(mean_var) or mean_var <= 0:
            mean_var = 1.0
        svm_gamma = 1 / (X.shape[1] * mean_var)
        min_gamma = svm_gamma * 0.1  # 更小
        max_gamma = svm_gamma * 10   # 更大
        param_grid = {
            'gamma': ['scale', 'auto'] + list(
                np.logspace(
                    np.log10(min_gamma), np.log10(max_gamma), num = 18
                )
            )
        }
        # 超参数选择：使用 GridSearchCV 在 train_val 上搜索 svm_gamma
        # 定义 PredefinedSplit：train 部分为 fold 0，val 部分为 -1（外部验证）
        test_fold = np.full(len(train_val_idx), 0)
        test_fold[len(train_idx):] = -1  # val 部分标记为 -1
        ps = PredefinedSplit(test_fold)
        grid = GridSearchCV(
            clf, param_grid,
            cv=ps, scoring='accuracy'  # 基于准确率选择最佳
        )
        grid.fit(X_train_val, y_train_val)
        best_svm_gamma = grid.best_params_['gamma']
        
        # 最终评估：用最佳模型在 test 上评估（GridSearchCV 默认 refit 最佳模型）
        clf = grid.best_estimator_
        logits = clf.predict_proba(X_test)
        y_pred = logits.argmax(axis=1)
        acc = float(accuracy_score(y_test, y_pred))
        f1 = float(f1_score(y_test, y_pred, average='weighted'))
        roc_auc = 0.0
        if num_cls > 2:
            roc_auc = float(roc_auc_score(y_test, logits, multi_class='ovo'))
        else:
            roc_auc = float(roc_auc_score(y_test, logits[:, 1]))
        
        # 记录结果
        print(f'dname = {dname}, alpha = {round(alpha, 6)}, gk_gamma = {round(gk_gamma, 8)}, best_svm_gamma = {best_svm_gamma}, acc = {round(acc, 6)}, f1 = {round(f1, 6)}, roc_auc = {round(roc_auc, 6)}')
        result_list.append(
            (
                round(alpha, 6), round(gk_gamma, 6), 
                round(acc, 4), round(f1, 4), round(roc_auc, 4), round(run_time, 4)
            )
        )
    return result_list

# Globals for multiprocessing to avoid passing large arrays per-task
_GLOBAL_OT_SIM = None
_GLOBAL_WL_SIM = None
_GLOBAL_Y = None
_GLOBAL_RUN_TIME = None

def _mp_worker_init(ot_sim_np, wl_sim_np, y, run_time):
    global _GLOBAL_OT_SIM, _GLOBAL_WL_SIM, _GLOBAL_Y, _GLOBAL_RUN_TIME
    _GLOBAL_OT_SIM = ot_sim_np
    _GLOBAL_WL_SIM = wl_sim_np
    _GLOBAL_Y = y
    _GLOBAL_RUN_TIME = run_time

def _train_worker_wrapper(args):
    # args: (dname, nu, tol, max_iter, alpha_list, gk_gamma, random_state)
    dname, nu, tol, max_iter, alpha_list, gk_gamma, random_state = args
    return train_and_eval_worker(
        dname, nu, tol, max_iter, alpha_list, gk_gamma,
        _GLOBAL_OT_SIM, _GLOBAL_WL_SIM, _GLOBAL_Y, _GLOBAL_RUN_TIME, random_state
    )

def multiprocess_unit(train_and_eval_worker, args_list):
    # Detect if args_list contains repeated large objects (ot_sim, wl_sim, y, run_time)
    use_shared = False
    if len(args_list) > 0 and len(args_list[0]) >= 10:
        # heuristic: positions 6,7,8,9 are ot_sim, wl_sim, y, run_time in current caller
        maybe_ot = args_list[0][6]
        maybe_wl = args_list[0][7]
        maybe_y = args_list[0][8]
        maybe_rt = args_list[0][9]
        if hasattr(maybe_ot, 'shape') and hasattr(maybe_wl, 'shape'):
            use_shared = True

    if use_shared:
        # extract shared large objects from first tuple and build trimmed args
        ot_sim_np = args_list[0][6]
        wl_sim_np = args_list[0][7]
        y = args_list[0][8]
        run_time = args_list[0][9]
        trimmed_args = [
            (a[0], a[1], a[2], a[3], a[4], a[5], a[10]) for a in args_list
        ]
        # spawn pool with initializer to set globals in workers
        nproc = get_recommended_nproc()
        pool = Pool(
            processes = nproc,
            initializer = _mp_worker_init,
            initargs = (ot_sim_np, wl_sim_np, y, run_time)
        )
        try:
            # use imap_unordered to stream results and reduce synchronization overhead
            chunk = max(1, len(trimmed_args) // (cpu_count() * 4))
            async_iter = pool.imap_unordered(_train_worker_wrapper, trimmed_args, chunksize=chunk)
            result_list = []
            for res in async_iter:
                result_list.extend(res)
        finally:
            pool.close()
            pool.join()
        return result_list
    else:
        pool = Pool(cpu_count())
        ms = pool.starmap_async(train_and_eval_worker, args_list)
        pool.close()
        pool.join()
        async_result = ms.get()
        result_list = []
        for ar in async_result:
            result_list.extend(ar)
        return result_list

def compute_svm_statistics(all_results):
    # 创建 DataFrame
    columns = ['alpha', 'gk_gamma', 'acc', 'f1_score', 'roc_auc_score', 'run_time']
    df = pd.DataFrame(all_results, columns=pd.Index(columns))
    
    # 分组计算 mean 和 std
    grouped = df.groupby(['alpha', 'gk_gamma'])
    agg_funcs = {
        'acc': ['mean', 'std'],
        'f1_score': ['mean', 'std'],
        'roc_auc_score': ['mean', 'std'],
        'run_time': 'mean'  # run_time 只取 mean
    }
    result = grouped.agg(agg_funcs)
    result.columns = ['acc', 'acc_std', 'f1_score', 'f1_score_std', 'roc_auc_score', 'roc_auc_score_std', 'run_time']
    result = result.reset_index()
    return result

# 多进程
def main(
    dataset_dict, niter_tn, niter_hcc, nu, alpha_list, tol, 
    gk_gamma_list, max_iter, n_repeats=10, random_state_base=42
):
    dname = dataset_dict['dname']
    result_base_dir = os.path.join(
        '.', 'our_experiments', 'cycle_complex_wgk', 'topowgk_outputs'
    )
    fout = open(
        os.path.join(result_base_dir, f'{dname}_svm_result_niter_tn_hcc_({niter_tn}_{niter_hcc}).csv'), 
        'w', encoding = 'utf-8', newline = ''
    )
    csv_writer = csv.writer(fout)
    csv_writer.writerow(
        [
            'alpha', 'gk_gamma', 'acc', 'acc_std', 'f1_score', 'f1_score_std', 
            'roc_auc_score', 'roc_auc_score_std', 'run_time'
        ]
    )
    
    # Cache logic for ot_dist_np, wl_sim_np, run_time
    cache_dir = os.path.join(os.path.join('.', 'our_experiments'), 'cache')
    os.makedirs(cache_dir, exist_ok=True)
    cache_file = os.path.join(cache_dir, f'{dname}_niter_tn_{niter_tn}_niter_hcc_{niter_hcc}.pkl')
    
    if os.path.exists(cache_file):
        print(f"Loading cached results from {cache_file}")
        with open(cache_file, 'rb') as f:
            ot_dist_np, wl_sim_np, y, run_time = pickle.load(f)
    else:
        graph_list, vlabel_list, edges_list, elabel_list, \
            vtx_hierarchical_cycle_contexts_list, vtx_triangulated_neighbors_list, deg_distr_list, y \
            = construct_dataset(dataset_dict)

        print(f"Computing results for {dname} with niter_tn={niter_tn}, niter_hcc={niter_hcc}")
        cyc_wl_gk = TopoWassersteinGraphKernel(niter_tn, niter_hcc, wl_normalized = True)
        ot_dist_np, wl_sim_np, run_time = cyc_wl_gk.fit_transform(
            dataset_dict['dataset_info'], graph_list, vlabel_list, edges_list, elabel_list, 
            vtx_hierarchical_cycle_contexts_list, vtx_triangulated_neighbors_list,
            deg_distr_list
        )
        # Save to cache
        with open(cache_file, 'wb') as f:
            pickle.dump((ot_dist_np, wl_sim_np, y, run_time), f)
        print(f"Cached results saved to {cache_file}")
    
    all_results = []
    for gk_gamma in gk_gamma_list:
        ot_sim_np = gauss_kernel(ot_dist_np, gk_gamma)
        for repeat in range(n_repeats):
            random_state = random_state_base + repeat
            alpha_sublist = get_start_end_indices(len(alpha_list))
            args_list = []
            for (alpha_start, alpha_end) in alpha_sublist:
                args_list.append(
                    (
                        dname, nu, tol, max_iter, alpha_list[alpha_start: alpha_end], 
                        gk_gamma, ot_sim_np, wl_sim_np, y, run_time, random_state
                    )
                )
            result_list = multiprocess_unit(train_and_eval_worker, args_list)
            all_results.extend(result_list)
    
    # 计算统计量
    result = compute_svm_statistics(all_results)
    
    # 保存结果
    for _, row in result.iterrows():
        csv_writer.writerow([
            round(float(row['alpha']), 6), round(float(row['gk_gamma']), 6),
            round(float(row['acc']), 4), round(float(row['acc_std']), 4),
            round(float(row['f1_score']), 4), round(float(row['f1_score_std']), 4),
            round(float(row['roc_auc_score']), 4), round(float(row['roc_auc_score_std']), 4),
            round(float(row['run_time']), 4)
        ])
    fout.flush()
    fout.close()

if __name__ == '__main__':
    dname_params = {
        'IMDB-MULTI': \
            {
                'niter_tn': 3, 'niter_hcc': 3, 'nu': 0.1, 'random_state': 42, 
                'alpha_list': np.logspace(-2, 0, num = 15, base = 10.0), 
                'gk_gamma_list': np.logspace(-4, -3, num = 20, base = 10.0)
            }, 
        'DHFR_MD': \
            {
                'niter_tn': 3, 'niter_hcc': 3, 'nu': 0.05, 'random_state': 42, 
                'alpha_list': np.logspace(-2, 0, num = 15, base = 10.0), 
                'gk_gamma_list': np.logspace(-4, -3, num = 20, base = 10.0)
            }, 
        'DHFR': \
            {
                'niter_tn': 3, 'niter_hcc': 3, 'nu': 0.1, 'random_state': 42, 
                'alpha_list': np.logspace(-2, 0, num = 15, base = 10.0), 
                'gk_gamma_list': np.logspace(-4, -2, num = 20, base = 10.0)
            }, 
        'ER_MD': \
            {
                'niter_tn': 3, 'niter_hcc': 3, 'nu': 0.05, 'random_state': 42, 
                'alpha_list': np.logspace(-2, 0, num = 15, base = 10.0), 
                'gk_gamma_list': np.logspace(-4, -2, num = 20, base = 10.0)
            }, 
        'BZR_MD': \
            {
                'niter_tn': 3, 'niter_hcc': 3, 'nu': 0.05, 'random_state': 42, 
                'alpha_list': np.logspace(-2, 0, num = 15, base = 10.0), 
                'gk_gamma_list': np.logspace(-3, -2, num = 20, base = 10.0), 
            }, 
        'IMDB-BINARY': \
            {
                'niter_tn': 3, 'niter_hcc': 3, 'nu': 0.05, 'random_state': 42, 
                'alpha_list': np.logspace(-2, 0, num = 15, base = 10.0), 
                'gk_gamma_list': np.logspace(-5, -3, num = 20, base = 10.0), 
            }, 
        'ENZYMES': \
            {
                'niter_tn': 3, 'niter_hcc': 3, 'nu': 0.01, 'random_state': 42, 
                'alpha_list': np.logspace(-2, 0, num = 15, base = 10.0), 
                'gk_gamma_list': np.logspace(-4, -3, num = 20, base = 10.0), 
            }, 
        'COX2': 
            {
                'niter_tn': 3, 'niter_hcc': 3, 'nu': 0.05, 'random_state': 256, 
                'alpha_list': np.logspace(-2, 0, num = 15, base = 10.0), 
                'gk_gamma_list': np.logspace(-5, -3, num = 20, base = 10.0), 
            }, 
        'PROTEINS': 
            {
                'niter_tn': 3, 'niter_hcc': 5, 'nu': 0.08, 'random_state': 3407,
                'alpha_list': np.logspace(-2, 0, num = 15, base = 10.0), 
                'gk_gamma_list': np.logspace(-7, -3, num = 100, base = 10.0), 
            }, 
        'Mutagenicity': 
            {
                'niter_tn': 3, 'niter_hcc': 3, 'nu': 0.05, 'random_state': 42, 
                'alpha_list': np.logspace(-2, 0, num = 15, base = 10.0), 
                'gk_gamma_list': np.logspace(-4, -3, num = 20, base = 10.0)
            },
        'NCI1': 
            {
                'niter_tn': 3, 'niter_hcc': 3, 'nu': 0.03, 'random_state': 3407, 
                'alpha_list': np.logspace(-2, 0, num = 15, base = 10.0), 
                'gk_gamma_list': np.logspace(-6, -3, num = 20, base = 10.0), 
            }, 
        'NCI109': 
            {
                'niter_tn': 3, 'niter_hcc': 3, 'nu': 0.01, 'random_state': 42, 
                'alpha_list': np.logspace(-2, 0, num = 15, base = 10.0), 
                'gk_gamma_list': np.logspace(-6, -3, num = 20, base = 10.0), 
            },  
        'BZR': \
            {
                'niter_tn': 3, 'niter_hcc': 3, 'nu': 0.1, 'random_state': 42, 
                'alpha_list': np.logspace(-2, 0, num = 15, base = 10.0), 
                'gk_gamma_list': np.logspace(-4, -2, num = 20, base = 10.0), 
            }, 
        'MSRC_21': \
            {
                'niter_tn': 3, 'niter_hcc': 3, 'nu': 0.03, 'random_state': 42, 
                'alpha_list': np.logspace(-2, 0, num = 15, base = 10.0), 
                'gk_gamma_list': np.logspace(-4, -3, num = 20, base = 10.0)
            }, 
        'MSRC_9': \
            {
                'niter_tn': 3, 'niter_hcc': 5, 'nu': 0.08, 'random_state': 3407, 
                'alpha_list': np.logspace(-2, 0, num = 15, base = 10.0), 
                'gk_gamma_list': np.logspace(-5, -2, num = 100, base = 10.0)
            }, 
        'COX2_MD': \
            {
                'niter_tn': 3, 'niter_hcc': 3, 'nu': 0.05, 'random_state': 42, 
                'alpha_list': np.logspace(-2, 0, num = 15, base = 10.0),
                'gk_gamma_list': np.logspace(-4, -3, num = 20, base = 10.0)
            }, 
    }    
    
    dname_baseline = {
        'BZR_MD': [0.619355, 0.109972, 0.618951, 0.167287, 0.670000, 0.108232],
        'MSRC_21': [0.865487, 0.043243, 0.862474, 0.041771, 0.987434, 0.032442],
        'MSRC_9': [0.904444, 0.067077, 0.902731, 0.067633, 0.985556, 0.028982],
        'ENZYMES': [0.563333, 0.065505, 0.557807, 0.067417, 0.834833, 0.046983],
        'COX2_MD': [0.575410, 0.091193, 0.573452, 0.119436, 0.600591, 0.100467],
        'DHFR_MD': [0.640506, 0.145119, 0.633295, 0.179665, 0.626222, 0.122307],
        'ER_MD': [0.701111, 0.162795, 0.701813, 0.217219, 0.758720, 0.164358],
        'IMDB-BINARY': [0.717000, 0.121821, 0.716502, 0.121155, 0.784195, 0.194736],
        'IMDB-MULTI': [0.402000, 0.054564, 0.396991, 0.070308, 0.593108, 0.051730],
        'COX2': [0.792553, 0.065194, 0.785728, 0.085700, 0.794260, 0.165638],
        'BZR': [0.841975, 0.164383, 0.841396, 0.155878, 0.870588, 0.118176],
        'DHFR': [0.806579, 0.249315, 0.804528, 0.248347, 0.889338, 0.253118],
        'PROTEINS': [0.772646, 0.032888, 0.770225, 0.034706, 0.835079, 0.032712],
        'Mutagenicity': [0.823733, 0.007661, 0.823752, 0.007742, 0.896813, 0.008614],
        'FIRSTMM_DB': [0.761538, 0.076495, 0.721538, 0.090519, 0.810227, 0.031787],
        'NCI1': [0.850608, 0.010618, 0.850312, 0.010691, 0.910675, 0.011267],
        'NCI109': [0.837772, 0.010161, 0.837347, 0.010254, 0.900845, 0.008164]
    }
    
    n_repeats = 10
    tol = 1e-6
    max_iter = -1 #1000000
    
    for dname, params in dname_params.items():
        dataset_dict = download_dataset(dname)
        for niter_hcc in range(0, 8):
            niter_tn = params['niter_tn']
            nu = params['nu']
            random_state = params['random_state']
            alpha_list = params['alpha_list']
            gk_gamma_list = params['gk_gamma_list']
            main(
                dataset_dict, niter_tn, niter_hcc, nu, alpha_list, tol, 
                gk_gamma_list, max_iter, n_repeats, random_state
            )
    # from postprocess_svm_results import filter_superior_results, select_best_results
    # for dname, params in dname_params.items():
    #     baseline_vals = dname_baseline[dname]
    #     niter_tn = params['niter_tn']
    #     niter_hcc = params['niter_hcc']
    #     filter_superior_results(dname, niter_tn, niter_hcc, baseline_vals)
    #     select_best_results(dname, niter_tn, niter_hcc)