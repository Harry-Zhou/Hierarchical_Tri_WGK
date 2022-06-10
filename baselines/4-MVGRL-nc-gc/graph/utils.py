''' Code adapted from https://github.com/fanyun-sun/InfoGraph '''
import torch as th
import torch.nn.functional as F

import math
import numpy as np

from sklearn.svm import LinearSVC
from sklearn.metrics import accuracy_score, f1_score, roc_auc_score
from sklearn.model_selection import GridSearchCV, StratifiedKFold

def linearsvc(embeds, labels):
    x = embeds.cpu().numpy()
    y = labels.cpu().numpy()
    num_cls = y.max() + 1
    params = {'C': [0.001, 0.01, 0.1, 1, 10, 100, 1000]}
    kf = StratifiedKFold(n_splits=10, shuffle=True, random_state=None)
    accuracies, f1_scores, roc_auc_scores = [], [], []
    for train_index, test_index in kf.split(x, y):
        x_train, x_test = x[train_index], x[test_index]
        y_train, y_test = y[train_index], y[test_index]
        classifier = GridSearchCV(LinearSVC(), params, cv=5, scoring='accuracy', verbose=0)
        classifier.fit(x_train, y_train)

        y_pred = classifier.predict(x_test)
        y_prob = classifier.predict_proba(x_test)

        accuracies.append(accuracy_score(y_test, y_prob))
        if num_cls > 2:
            f1_scores.append(f1_score(y_test, y_pred, average = 'micro'))
            roc_auc_scores.append(roc_auc_score(y_test, y_prob, multi_class = 'ovr'))
        else:
            f1_scores.append(f1_score(y_test, y_pred))
            roc_auc_scores.append(roc_auc_score(y_test, y_prob[:, 1]))
    return np.mean(accuracies), np.mean(f1_scores), np.mean(roc_auc_scores)

def get_positive_expectation(p_samples, average=True):
    """Computes the positive part of a JS Divergence.
    Args:
        p_samples: Positive samples.
        average: Average the result over samples.
    Returns:
        th.Tensor
    """
    log_2 = math.log(2.)
    Ep = log_2 - F.softplus(- p_samples)

    if average:
        return Ep.mean()
    else:
        return Ep


def get_negative_expectation(q_samples, average=True):
    """Computes the negative part of a JS Divergence.
    Args:
        q_samples: Negative samples.
        average: Average the result over samples.
    Returns:
        th.Tensor
    """
    log_2 = math.log(2.)
    Eq = F.softplus(-q_samples) + q_samples - log_2

    if average:
        return Eq.mean()
    else:
        return Eq


def local_global_loss_(l_enc, g_enc, graph_id):

    num_graphs = g_enc.shape[0]
    num_nodes = l_enc.shape[0]

    device = g_enc.device

    pos_mask = th.zeros((num_nodes, num_graphs)).to(device)
    neg_mask = th.ones((num_nodes, num_graphs)).to(device)

    for nodeidx, graphidx in enumerate(graph_id):

        pos_mask[nodeidx][graphidx] = 1.
        neg_mask[nodeidx][graphidx] = 0.

    res = th.mm(l_enc, g_enc.t())

    E_pos = get_positive_expectation(res * pos_mask, average=False).sum()
    E_pos = E_pos / num_nodes
    E_neg = get_negative_expectation(res * neg_mask, average=False).sum()
    E_neg = E_neg / (num_nodes * (num_graphs - 1))

    return E_neg - E_pos
