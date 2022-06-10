import argparse
import copy
import random
import warnings
from collections import defaultdict

import networkx as nx
import numpy as np
import scipy.sparse as sp
import torch
import torch.nn.functional as F
from scipy import sparse as sp
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import accuracy_score, f1_score, roc_auc_score
from sklearn.model_selection import StratifiedKFold
from sklearn.multiclass import OneVsRestClassifier
from sklearn.utils import shuffle as skshuffle
from tqdm import tqdm
import os
import csv

from gcc.datasets.data_util import create_node_classification_dataset
from gcc.tasks import build_model

warnings.filterwarnings("ignore")

class NodeClassification(object):
    """
    Node classification task.
    """

    def __init__(self, dataset, model, hidden_size, num_shuffle, seed, **model_args):
        self.data = create_node_classification_dataset(dataset).data
        self.label_matrix = self.data.y
        self.num_nodes, self.num_classes = self.data.y.shape

        self.model = build_model(model, hidden_size, **model_args)
        self.hidden_size = hidden_size
        self.num_shuffle = num_shuffle
        self.seed = seed

    def train(self):
        G = nx.Graph()
        G.add_edges_from(self.data.edge_index.t().tolist())
        embeddings = self.model.train(G)

        # Map node2id
        features_matrix = np.zeros((self.num_nodes, self.hidden_size))
        for vid, node in enumerate(G.nodes()):
            features_matrix[node] = embeddings[vid]

        label_matrix = torch.Tensor(self.label_matrix)

        return self._evaluate(features_matrix, label_matrix, self.num_shuffle)

    def _evaluate(self, features_matrix, label_matrix, num_shuffle):
        # shuffle, to create train/test groups
        skf = StratifiedKFold(n_splits=10, shuffle=True, random_state=self.seed)
        idx_list = []
        labels = label_matrix.argmax(axis=1).squeeze().tolist()
        num_cls = label_matrix.max() + 1
        for idx in skf.split(np.zeros(len(labels)), labels):
            idx_list.append(idx)

        # score each train/test group
        # all_results_acc = defaultdict(list)
        # all_results_f1 = defaultdict(list)
        # all_results_auc = defaultdict(list)
        accuracies, f1_scores, roc_auc_scores = [], [], []

        for train_idx, test_idx in idx_list:

            X_train = features_matrix[train_idx]
            y_train = label_matrix[train_idx]

            X_test = features_matrix[test_idx]
            y_test = label_matrix[test_idx]

            clf = TopKRanker(LogisticRegression(C=1000))
            clf.fit(X_train, y_train)

            # find out how many labels should be predicted
            top_k_list = y_test.sum(axis=1).long().tolist()
            y_pred = clf.predict(X_test, top_k_list)
            y_prob = clf.predict_proba(X_test, top_k_list)
            acc_result = accuracy_score(y_test, y_pred)
            if num_cls > 2:
                f1_result = f1_score(y_test, y_pred, average="micro")
                roc_auc_result = roc_auc_score(y_test, y_prob, multi_class = 'ovr')
            else:
                f1_result = f1_score(y_test, y_pred)
                roc_auc_result = roc_auc_score(y_test, y_prob[:, 1])
            
            accuracies.append(acc_result)
            f1_scores.append(f1_result)
            roc_auc_scores.append(roc_auc_result)
            # all_results_acc[""].append(acc_result)
            # all_results_f1[""].append(f1_result)
            # all_results_auc[""].append(roc_auc_result)

        return np.mean(accuracies), np.mean(f1_scores), np.mean(roc_auc_scores)
        # return dict((f"micro-f1{train_percent}", \
        #         sum(all_results_f1[train_percent]) / len(all_results_f1[train_percent])) \
        #     for train_percent in sorted(all_results_f1.keys()))

class TopKRanker(OneVsRestClassifier):
    def predict(self, X, top_k_list):
        assert X.shape[0] == len(top_k_list)
        probs = np.asarray(super(TopKRanker, self).predict_proba(X))
        all_labels = sp.lil_matrix(probs.shape)
        
        for i, k in enumerate(top_k_list):
            probs_ = probs[i, :]
            labels = self.classes_[probs_.argsort()[-k:]].tolist()
            for label in labels:
                all_labels[i, label] = 1
        return all_labels

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--transferring-mode", type=str, required=True)
    parser.add_argument("--dataset", type=str, required=True)
    parser.add_argument("--model", type=str, required=True)
    parser.add_argument("--hidden-size", type=int, required=True)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--num-shuffle", type=int, default=10)
    parser.add_argument("--emb-path", type=str, default="")
    args = parser.parse_args()
    task = NodeClassification(
        args.dataset,
        args.model,
        args.hidden_size,
        args.num_shuffle,
        args.seed,
        emb_path=args.emb_path,
    )
    
    with open(os.path.join('.', 'nc_results', args.dataset + '_results.csv'), 'w', newline = '', encoding = 'utf-8') as fout:
        csv_writer = csv.writer(fout)
        csv_writer.writerow(['Dataset', \
            'Transferring Mode', \
            'Accuracy', \
            'Micro_F1', \
            'Roc_Auc'])
        acc, f1, roc_auc = task.train()
        csv_writer.writerow([[args.dataset, args.transferring_mode, acc, f1, roc_auc]])

