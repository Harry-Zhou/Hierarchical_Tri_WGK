import os.path
from sklearn import svm
from sklearn.model_selection import train_test_split, StratifiedShuffleSplit
from sklearn.metrics import accuracy_score, f1_score, roc_auc_score
from sklearn.impute import SimpleImputer
import csv
import time
import numpy as np
from sklearn.preprocessing import MinMaxScaler

from grakel.datasets import fetch_dataset, get_dataset_info
from grakel.kernels import \
    WeisfeilerLehman, \
    WeisfeilerLehmanOptimalAssignment, \
    LovaszTheta, \
    VertexHistogram, \
    CoreFramework, \
    NeighborhoodSubgraphPairwiseDistance, \
    NeighborhoodHash, \
    OddSth, \
    PyramidMatch, \
    Propagation, \
    GraphletSampling, \
    HadamardCode

def get_graph_kernel(gk_name, dataset_info, normalize = True):
    if gk_name == 'WeisfeilerLehman':
        return WeisfeilerLehman(
            n_iter = 5, 
            base_graph_kernel = VertexHistogram, 
            normalize = normalize
        )
    elif gk_name == 'WeisfeilerLehmanOptimalAssignment': 
        return WeisfeilerLehmanOptimalAssignment(
            n_iter = 5, 
            normalize = normalize
        )
    elif gk_name == 'LovaszTheta':
        return LovaszTheta(
            normalize = True, 
            random_state = 42
        )
    elif gk_name == 'CoreFramework':
        return CoreFramework(
            normalize = normalize
        )
    elif gk_name == 'NeighborhoodSubgraphPairwiseDistance':
        return NeighborhoodSubgraphPairwiseDistance(
            normalize = normalize
        )
    elif gk_name == 'NeighborhoodHash':
        return NeighborhoodHash(
            normalize = normalize, 
            nh_type = 'count_sensitive'
        )
    elif gk_name == 'OddSth':
        return OddSth(
            normalize = normalize, 
            h = 4
        )
    elif gk_name == 'PyramidMatch':
        return PyramidMatch(
            normalize = normalize, 
            with_labels = dataset_info['nl']
        )
    elif gk_name == 'Propagation':
        return Propagation(
            normalize = normalize, 
            random_state = 42
        )
    elif gk_name == 'GraphletSampling':
        return GraphletSampling(
            normalize = normalize, 
            random_state = 42, 
            k = 5
        )
    elif gk_name == 'HadamardCode':
        return HadamardCode()
    else:
        raise Exception("Please input a valid gk_name!")
  
def train_test_gk_baselines(
    nu, 
    gk_name, 
    normalize, 
    dataset_info, 
    G_train, 
    G_test, 
    y_train, 
    y_test
):
    num_cls = np.unique(np.concatenate([y_train, y_test])).shape[0]
    gk = get_graph_kernel(gk_name, dataset_info, normalize)
    start_time = time.time()
    X_train = gk.fit_transform(G_train)
    end_time = time.time()
    X_test = gk.transform(G_test)
    # 处理缺失值
    imp_mean = SimpleImputer(missing_values = np.nan, strategy = 'mean')
    X_train = imp_mean.fit_transform(X_train)
    X_test = imp_mean.fit_transform(X_test)
    
    clf = None
    if num_cls == 2:
        clf = svm.NuSVC(
            nu = nu, 
            kernel = 'rbf', 
            gamma = 'scale', 
            tol = 1e-8, 
            probability = True, 
            max_iter = -1, 
            random_state = None
        )
    else:
        clf = svm.NuSVC(
            nu = nu, 
            kernel = 'rbf', 
            gamma = 'scale', 
            tol = 1e-8, 
            probability = True, 
            decision_function_shape = 'ovo',
            max_iter = -1, 
            random_state = None
        )
    
    clf.fit(X_train, y_train)
    y_test_proba = clf.predict_proba(X_test)
    y_test_pred = clf.predict(X_test) # np.argmax(y_test_proba, axis = 1)
    acc = accuracy_score(y_test, y_test_pred)
    f1 = f1_score(y_test, y_test_pred, average = 'weighted')  
    if num_cls > 2: # 多分类, ovo: One-Versus-Rest
        roc_auc = roc_auc_score(y_test, y_test_proba, multi_class = 'ovo')
    else: # 二分类
        roc_auc = roc_auc_score(y_test, y_test_proba[:, 1])
    
    return acc, f1, roc_auc, end_time - start_time

def eval_gk_baselines(dname, gk_name, nu, normalize, n_splits, train_size = 0.8, test_size = 0.2):
    dataset = fetch_dataset(
        dname, 
        produce_labels_nodes = True, 
        verbose=False
    )
    G, y = dataset.data, dataset.target
    # Splits the dataset into a training and a test set
    results = []
    cv = StratifiedShuffleSplit(n_splits = n_splits, train_size = train_size, test_size = test_size, random_state = 42)
    for fidx, (gidx_train, gidx_test) in enumerate(cv.split(G, y)):
        print(f'当前数据集: {dname}, 当前图核: {gk_name}, fidx: {fidx}')
        G_train = [G[i] for i in gidx_train]
        G_test = [G[i] for i in gidx_test]
        y_train = np.array([y[i] for i in gidx_train])
        y_test = np.array([y[i] for i in gidx_test])
        dataset_info = get_dataset_info(dname) # nl, el, na, ea
        acc, f1, roc_auc, run_time = train_test_gk_baselines(
            nu, gk_name, normalize, 
            dataset_info, G_train, G_test, y_train, y_test
        )
        results.append((dname, gk_name, normalize, fidx, acc, f1, roc_auc, run_time))
    
    return results

def comp_final_metrics(dnames, gk_names):
    import pandas as pd
    read_chunks = pd.read_csv(
        os.path.join('.', 'eval_grakel_baselines', 'outputs (nu=0.01)', 'baseline_grakels_result_addition.csv'), 
        encoding = 'utf-8', iterator = True, 
    )
    chunk_list = []
    for rchunk in read_chunks:
        chunk_list.append(rchunk)
    base_df = pd.concat(chunk_list, axis = 0, ignore_index = True)
    base_df[['acc', 'f1_score', 'roc_auc_score', 'run_time']] = \
        base_df[['acc', 'f1_score', 'roc_auc_score', 'run_time']].apply(pd.to_numeric)
    final_result = {
        'dataset': [], 
        'graph_kernel': [], 
        'normalize': [], 
        'acc_mean': [], 
        'acc_std': [], 
        'f1_score_mean': [], 
        'f1_score_std': [], 
        'roc_auc_score_mean': [], 
        'roc_auc_score_std': [], 
        'run_time': []
    }
    for dname in dnames:
        for gk_name in gk_names:
            final_result['dataset'].append(dname)
            final_result['graph_kernel'].append(gk_name)
            selected_df = base_df[(base_df['dataset'] == dname) & (base_df['graph_kernel'] == gk_name)]
            final_result['normalize'].append(True)
            final_result['acc_mean'].append(selected_df['acc'].mean())
            final_result['acc_std'].append(selected_df['acc'].std())
            final_result['f1_score_mean'].append(selected_df['f1_score'].mean())
            final_result['f1_score_std'].append(selected_df['f1_score'].std())
            final_result['roc_auc_score_mean'].append(selected_df['roc_auc_score'].mean())
            final_result['roc_auc_score_std'].append(selected_df['roc_auc_score'].std())
            final_result['run_time'].append(selected_df['run_time'].mean())
    final_df = pd.DataFrame(data = final_result)
    final_df.to_csv(
        os.path.join('.', 'eval_grakel_baselines', 'outputs (nu=0.01)', f'final_baseline_grakels_result_addition.csv'), 
        index = False, sep = ','
    )

if __name__ == '__main__':
    n_splits = 10
    nu = 0.01
    base_dir = os.path.join('.', 'eval_grakel_baselines', f'outputs (nu={nu})')
    fout_test = open(
        os.path.join(base_dir, 'baseline_grakels_result_addition.csv'), 
        'w', encoding = 'utf-8', newline = ''
    )
    csv_writer = csv.writer(fout_test)
    csv_writer.writerow(
        [
            'dataset', 
            'graph_kernel', 
            'normalize', 
            'fold_idx', 
            'acc', 
            'f1_score', 
            'roc_auc_score', 
            'run_time'
        ]
    )
    
    dnames = [
        # 'BZR_MD', 
        # 'MSRC_21', 
        # 'MSRC_9', 
        # 'ENZYMES', 
        # 'AIDS', 
        # 'COX2_MD', 
        # 'DHFR_MD', 
        # 'ER_MD', 
        # 'IMDB-BINARY', 
        # 'IMDB-MULTI', 
        # 'COX2', 
        # 'BZR', 
        # 'DHFR', 
        # 'PROTEINS', 
        # 'Mutagenicity', 
        'NCI1', 
        'NCI109'
    ]
    
    gk_names = [
        # 'CoreFramework', 
        'GraphletSampling', 
        'HadamardCode', 
        'WeisfeilerLehman', 
        'WeisfeilerLehmanOptimalAssignment', 
        'Propagation', 
        'NeighborhoodHash', 
        'PyramidMatch', 
        'OddSth', 
    ]
    for normalize in [True,]:
        for dname in dnames:
            train_size = 0.8 if dname != 'FIRSTMM_DB' else 0.7
            test_size = 1.0 - train_size
            for gk_name in gk_names:
                results = eval_gk_baselines(dname, gk_name, nu, normalize, n_splits, train_size, test_size)
                csv_writer.writerows(results)
            fout_test.flush()
    fout_test.close()
    comp_final_metrics(dnames, gk_names)
