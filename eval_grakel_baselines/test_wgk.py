from dgl.data import TUDataset
import dgl
import os.path
import numpy as np
from sklearn import svm
from sklearn.metrics import accuracy_score, f1_score, roc_auc_score
import csv
import argparse
import time
import sys
sys.path.append('..')

from wgk_utils import (
    comp_graph_similarity_mat,
    gen_train_test_dataset,
    get_param_setting,
)

def load_gc_test_datasets():
    raw_dir = os.path.join('..', 'datasets_gc')
    return {
        'MSRC_21': TUDataset(
            'MSRC_21', raw_dir = raw_dir
        ), 
        'MSRC_9': TUDataset(
            'MSRC_9', raw_dir = raw_dir
        ), 
        'ENZYMES': TUDataset(
            'ENZYMES', raw_dir = raw_dir
        ), 
        'AIDS': TUDataset(
            'AIDS', raw_dir = raw_dir, 
        ), 
        'BZR_MD': TUDataset(
            'BZR_MD', raw_dir = raw_dir
        ), 
        'COX2_MD': TUDataset(
            'COX2_MD', raw_dir = raw_dir
        ), 
        'DHFR_MD': TUDataset(
            'DHFR_MD', raw_dir = raw_dir
        ), 
        'ER_MD': TUDataset(
            'ER_MD', raw_dir = raw_dir
        ), 
        'MUTAG': TUDataset(
            'MUTAG', raw_dir = raw_dir
        ), 
        'PROTEINS': TUDataset(
            'PROTEINS', raw_dir = raw_dir
        ), 
        'DD': TUDataset(
            'DD', raw_dir = raw_dir
        ), 
        'Mutagenicity': TUDataset(
            'Mutagenicity', raw_dir = raw_dir
        ), 
        'NCI1': TUDataset(
            'NCI1', raw_dir = raw_dir
        ), 
        'NCI109': TUDataset(
            'NCI109', raw_dir = raw_dir
        ), 
    }

def train_test_wgk(
    graph_list, label_list, len_rw_trace, num_rws, is_normalized, 
    p, eigval_ndigits, prob_ndigits
):
    run_time = 0.
    start_time = time.time()
    gsim_mat = comp_graph_similarity_mat(
        graph_list, len_rw_trace, num_rws, is_normalized, 
        p, eigval_ndigits, prob_ndigits
    ).detach().cpu().numpy()
    end_time = time.time()
    run_time = end_time - start_time
    
    X_train, X_test, y_train, y_test = gen_train_test_dataset(
        gsim_mat, label_list
    ) # 以ndarray的形式返回
    clf = svm.SVC(
        kernel = 'rbf', 
        decision_function_shape = 'ovo', 
        probability = True
    )
    clf.fit(X_train, y_train)
    y_test_proba = clf.predict_proba(X_test)
    y_test_pred = np.argmax(y_test_proba, axis = 1)
    print(f'X_train: {X_train}')
    print(f'X_test: {X_test}')
    print(f'y_test_pred: {y_test_pred}')
    
    num_cls = y_test_proba.shape[1]
    acc, f1, roc_auc = 0.0, 0.0, 0.0
    acc = accuracy_score(y_test, y_test_pred)
    if num_cls > 2: # 多分类, ovr: One-Versus-Rest
        roc_auc = roc_auc_score(
            y_test, 
            y_test_proba, 
            multi_class = 'ovo'
        )
        # average = macro, micro, weighted
        f1 = f1_score(y_test, y_test_pred, average = 'macro')
    else: # 二分类
        roc_auc = roc_auc_score(y_test, y_test_proba[:, 1])
        f1 = f1_score(y_test, y_test_pred)
    print(f'acc = {acc}, f1 = {f1}, roc_auc = {roc_auc}')
    return acc, f1, roc_auc, run_time

def eval_wgk(len_rw_trace, num_rws):
    param_setting = get_param_setting()
    is_normalized = param_setting['is_normalized']
    p = param_setting['p']
    eigval_ndigits = param_setting['eigval_ndigits']
    prob_ndigits = param_setting['prob_ndigits']
    
    fout_test = open(
        os.path.join(
            '.', 
            'outputs', 
            f'wgk_result_len_rw_trace_{len_rw_trace}_num_rws_{num_rws}.csv'
        ), 
        'w', 
        encoding = 'utf-8', 
        newline = ''
    )
    csv_writer = csv.writer(fout_test)
    csv_writer.writerow(
        [
            'dataset_name', 
            'acc', 
            'f1_score', 
            'roc_auc_score', 
            'run_time'
        ]
    )
    for dname, dataset in load_gc_test_datasets().items():
        if dname != 'MSRC_21':
            continue 
        graph_list, label_list = [], []
        for g, l in dataset:
            g = dgl.to_simple(g)
            graph_list.append(g)
            label_list.append(l.item())
        print(f'graph number: {[g.num_nodes() for g in graph_list]}')
        
        acc, f1, roc_auc, run_time = train_test_wgk(
            graph_list, label_list, len_rw_trace, num_rws, is_normalized, 
            p, eigval_ndigits, prob_ndigits
        )
        csv_writer.writerow(
            [
                dname, 
                round(acc, 6), 
                round(f1, 6), 
                round(roc_auc, 6), 
                round(run_time, 6)
            ]
        )
        fout_test.flush()
    fout_test.close()

if __name__=='__main__':
    # parser = argparse.ArgumentParser(description = 'Wasserstein 图核')
    # parser.add_argument('--lenRwTrace', 
    #     type = int, 
    #     help = '每次随机游走的轨迹长度', 
    #     required = True)
    # parser.add_argument('--numRws', 
    #     type = int, 
    #     help = '在每个节点进行几次随机游走', 
    #     required = True)
    
    # args = parser.parse_args()
    # eval_wgk(args.lenRwTrace, args.numRws)
    
    eval_wgk(6, 50)
    
    # for dname, dataset in load_gc_test_datasets().items():
    #     if dname != 'MSRC_9':
    #         continue
    #     print(f'加载数据集: {dname}')
    #     graph_list, label_list = [], []
    #     for g, l in dataset:
    #         g = dgl.to_simple(g)
    #         graph_list.append(g)
    #         label_list.append(l.item())
    #     print(f'{dname}: {label_list}')
