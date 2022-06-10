from sklearn import preprocessing
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.manifold import TSNE
from sklearn.metrics import accuracy_score, f1_score, roc_auc_score
from sklearn.model_selection import GridSearchCV, KFold, StratifiedKFold
from sklearn.model_selection import cross_val_score
from sklearn.neural_network import MLPClassifier
from sklearn.svm import SVC, LinearSVC
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
import torch
import torch.nn as nn

def draw_plot(datadir, DS, embeddings, fname, max_nodes=None):
    return
    import seaborn as sns
    graphs = read_graphfile(datadir, DS, max_nodes=max_nodes)
    labels = [graph.graph['label'] for graph in graphs]

    labels = preprocessing.LabelEncoder().fit_transform(labels)
    x, y = np.array(embeddings), np.array(labels)
    print('fitting TSNE ...')
    x = TSNE(n_components=2).fit_transform(x)

    plt.close()
    df = pd.DataFrame(columns=['x0', 'x1', 'Y'])

    df['x0'], df['x1'], df['Y'] = x[:,0], x[:,1], y
    sns.pairplot(x_vars=['x0'], y_vars=['x1'], data=df, hue="Y", size=5)
    plt.legend()
    plt.savefig(fname)

class LogReg(nn.Module):
    def __init__(self, ft_in, nb_classes):
        super(LogReg, self).__init__()
        self.fc = nn.Linear(ft_in, nb_classes)

        for m in self.modules():
            self.weights_init(m)

    def weights_init(self, m):
        if isinstance(m, nn.Linear):
            torch.nn.init.xavier_uniform_(m.weight.data)
            if m.bias is not None:
                m.bias.data.fill_(0.0)

    def forward(self, seq):
        ret = self.fc(seq)
        return ret

def logistic_classify(x, y):
    nb_classes = np.unique(y).shape[0]
    xent = nn.CrossEntropyLoss()
    hid_units = x.shape[1]
    num_cls = y.max() + 1
    accuracies, f1_scores, roc_auc_scores = [], [], []
    kf = StratifiedKFold(n_splits=10, shuffle=True, random_state=None)
    for train_index, test_index in kf.split(x, y):
        train_embs, test_embs = x[train_index], x[test_index]
        train_lbls, test_lbls= y[train_index], y[test_index]

        train_embs, train_lbls = torch.from_numpy(train_embs).cuda(), torch.from_numpy(train_lbls).cuda()
        test_embs, test_lbls= torch.from_numpy(test_embs).cuda(), torch.from_numpy(test_lbls).cuda()

        log = LogReg(hid_units, nb_classes)
        log.cuda()
        opt = torch.optim.Adam(log.parameters(), lr=0.01, weight_decay=0.0)

        best_val = 0
        test_acc = None
        for it in range(100):
            log.train()
            opt.zero_grad()

            logits = log(train_embs)
            loss = xent(logits, train_lbls)

            loss.backward()
            opt.step()

        test_lbls = test_lbls.detach().cpu().numpy()
        logits = log(test_embs).detach().cpu().numpy()
        preds = torch.argmax(logits, dim=1).detach().cpu().numpy()

        #acc = torch.sum(preds == test_lbls).float() / test_lbls.shape[0]
        accuracies.append(accuracy_score(test_lbls, preds))
        if num_cls > 2:
            f1_scores.append(f1_score(test_lbls, preds, average = 'micro'))
            roc_auc_scores.append(roc_auc_score(test_lbls, logits, multi_class = 'ovr'))
        else:
            f1_scores.append(f1_score(test_lbls, preds))
            roc_auc_scores.append(roc_auc_score(test_lbls, logits[:, 1]))
    return np.mean(accuracies), np.mean(f1_scores), np.mean(roc_auc_scores)

def svc_classify(x, y, search):
    kf = StratifiedKFold(n_splits=10, shuffle=True, random_state=None)
    num_cls = y.max() + 1
    accuracies, f1_scores, roc_auc_scores = [], [], []
    for train_index, test_index in kf.split(x, y):

        x_train, x_test = x[train_index], x[test_index]
        y_train, y_test = y[train_index], y[test_index]
        # x_train, x_val, y_train, y_val = train_test_split(x_train, y_train, test_size=0.1)
        if search:
            params = {'C':[0.001, 0.01,0.1,1,10,100,1000]}
            classifier = GridSearchCV(SVC(), params, cv=5, scoring='accuracy', verbose=0)
        else:
            classifier = SVC(C=10)
        classifier.fit(x_train, y_train)
        accuracies.append(accuracy_score(y_test, classifier.predict(x_test)))
        if num_cls > 2:
            f1_scores.append(f1_score(y_test, classifier.predict(x_test), average = 'micro'))
            roc_auc_scores.append(roc_auc_scores(y_test, classifier.predict_proba(x_test), multi_class = 'ovr'))
        else:
            f1_scores.append(f1_score(y_test, classifier.predict(x_test)))
            roc_auc_scores.append(roc_auc_scores(y_test, classifier.predict_proba(x_test)[:, 1]))
    return np.mean(accuracies), np.mean(f1_scores), np.mean(roc_auc_scores)

def randomforest_classify(x, y, search):
    kf = StratifiedKFold(n_splits=10, shuffle=True, random_state=None)
    num_cls = y.max() + 1
    accuracies, f1_scores, roc_auc_scores = [], [], []
    for train_index, test_index in kf.split(x, y):

        x_train, x_test = x[train_index], x[test_index]
        y_train, y_test = y[train_index], y[test_index]
        if search:
            params = {'n_estimators': [100, 200, 500, 1000]}
            classifier = GridSearchCV(RandomForestClassifier(), params, cv=5, scoring='accuracy', verbose=0)
        else:
            classifier = RandomForestClassifier()
        classifier.fit(x_train, y_train)
        accuracies.append(accuracy_score(y_test, classifier.predict(x_test)))
        if num_cls > 2:
            f1_scores.append(f1_score(y_test, classifier.predict(x_test), average = 'micro'))
            roc_auc_scores.append(roc_auc_scores(y_test, classifier.predict_proba(x_test), multi_class = 'ovr'))
        else:
            f1_scores.append(f1_score(y_test, classifier.predict(x_test)))
            roc_auc_scores.append(roc_auc_scores(y_test, classifier.predict_proba(x_test)[:, 1]))
    return np.mean(accuracies), np.mean(f1_scores), np.mean(roc_auc_scores)

def linearsvc_classify(x, y, search):
    kf = StratifiedKFold(n_splits=10, shuffle=True, random_state=None)
    num_cls = y.max() + 1
    accuracies, f1_scores, roc_auc_scores = [], [], []
    for train_index, test_index in kf.split(x, y):

        x_train, x_test = x[train_index], x[test_index]
        y_train, y_test = y[train_index], y[test_index]
        if search:
            params = {'C':[0.001, 0.01,0.1,1,10,100,1000]}
            classifier = GridSearchCV(LinearSVC(), params, cv=5, scoring='accuracy', verbose=0)
        else:
            classifier = LinearSVC(C=10)
        classifier.fit(x_train, y_train)
        accuracies.append(accuracy_score(y_test, classifier.predict(x_test)))
        if num_cls > 2:
            f1_scores.append(f1_score(y_test, classifier.predict(x_test), average = 'micro'))
            roc_auc_scores.append(roc_auc_scores(y_test, classifier.predict_proba(x_test), multi_class = 'ovr'))
        else:
            f1_scores.append(f1_score(y_test, classifier.predict(x_test)))
            roc_auc_scores.append(roc_auc_scores(y_test, classifier.predict_proba(x_test)[:, 1]))
    return np.mean(accuracies), np.mean(f1_scores), np.mean(roc_auc_scores)

def evaluate_embedding(embeddings, labels, method, search=True):
    """
    method in ['logreg', 'svc', 'linearsvc', 'randomforest']
    """
    labels = preprocessing.LabelEncoder().fit_transform(labels)
    x, y = np.array(embeddings), np.array(labels)
    # print(x.shape, y.shape)

    if method == 'logreg':
        accuracies = [logistic_classify(x, y)[0] for _ in range(1)]
        f1_scores = [logistic_classify(x, y)[1] for _ in range(1)]
        roc_auc_scores = [logistic_classify(x, y)[2] for _ in range(1)]
        # print(logreg_accuracies)
        print('LogReg Acc: ', np.mean(accuracies))
        print('LogReg F1: ', np.mean(f1_scores))
        print('LogReg Roc Auc: ', np.mean(roc_auc_scores))
    elif method == 'svc':
        accuracies = [svc_classify(x,y, search)[0] for _ in range(1)]
        f1_scores = [svc_classify(x, y)[1] for _ in range(1)]
        roc_auc_scores = [svc_classify(x, y)[2] for _ in range(1)]
        # print(svc_accuracies)
        print('SVC Acc', np.mean(accuracies))
        print('SVC F1: ', np.mean(f1_scores))
        print('SVC Roc Auc: ', np.mean(roc_auc_scores))
    elif method == 'linearsvc':
        accuracies = [linearsvc_classify(x, y, search) for _ in range(1)]
        f1_scores = [linearsvc_classify(x, y)[1] for _ in range(1)]
        roc_auc_scores = [linearsvc_classify(x, y)[2] for _ in range(1)]
        # print(linearsvc_accuracies)
        print('LinearSvc', np.mean(accuracies))
        print('LinearSVC F1: ', np.mean(f1_scores))
        print('LinearSVC Roc Auc: ', np.mean(roc_auc_scores))
    elif method == 'randomforest':
        accuracies = [randomforest_classify(x, y, search) for _ in range(1)]
        f1_scores = [randomforest_classify(x, y)[1] for _ in range(1)]
        roc_auc_scores = [randomforest_classify(x, y)[2] for _ in range(1)]
        # print(randomforest_accuracies)
        print('randomforest', np.mean(accuracies))
        print('randomforest F1: ', np.mean(f1_scores))
        print('randomforest Roc Auc: ', np.mean(roc_auc_scores))

    return np.mean(accuracies), np.mean(f1_scores), np.mean(roc_auc_scores)

if __name__ == '__main__':
    evaluate_embedding('./data', 'ENZYMES', np.load('tmp/emb.npy'))
