import numpy as np
from sklearn.model_selection import train_test_split
from sklearn.svm import SVC
from sklearn.metrics import accuracy_score, roc_auc_score
from grakel.datasets import fetch_dataset
from sklearn.impute import SimpleImputer
from grakel import Graph
from grakel.kernels import \
    ShortestPath, \
    SubgraphMatching, \
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
    RandomWalkLabeled, \
    ShortestPathAttr, \
    GraphHopper, \
    GraphletSampling, \
    MultiscaleLaplacian, \
    RandomWalk, \
    EdgeHistogram, \
    HadamardCode, \
    SvmTheta
from grakel.datasets import get_dataset_info

dname = 'MSRC_21'
dataset = fetch_dataset(dname, produce_labels_nodes = True, verbose=False)
G, y = dataset.data, dataset.target
edges = G[0][0]
for u,v in list(edges):
    if (v,u) not in edges:
        print(f'{v}-{u} in edges')

# dataset_info = get_dataset_info(dname)
# print(dataset_info)

# Splits the dataset into a training and a test set
G_train, G_test, y_train, y_test = train_test_split(G, y, test_size=0.2, random_state=42)
# Uses the Weisfeiler-Lehman subtree kernel to generate the kernel matrices
# gk = WeisfeilerLehman(n_iter=4, base_graph_kernel=VertexHistogram, normalize=True)
# gk = ShortestPath(normalize = True)
# gk = SubgraphMatching(normalize = True)
# gk = WeisfeilerLehmanOptimalAssignment(n_iter = 5, normalize = True)
# gk = LovaszTheta(normalize = True, random_state = 42)
# gk = NeighborhoodSubgraphPairwiseDistance(normalize = True)
# gk = CoreFramework(normalize = True)
# gk = OddSth(normalize = True)
# gk = NeighborhoodHash(normalize = True, nh_type = 'count_sensitive')
# gk = PyramidMatch(normalize = True)
# gk = Propagation(normalize = True, random_state = 42)
# gk = RandomWalkLabeled(normalize = True)
for gk in [
    WeisfeilerLehman(n_iter=5, base_graph_kernel=VertexHistogram, normalize=True), 
    # ShortestPath(normalize = True), 
    # WeisfeilerLehmanOptimalAssignment(n_iter = 5, normalize = True), 
    # LovaszTheta(normalize = True, random_state = 42), 
    # NeighborhoodSubgraphPairwiseDistance(normalize = True), 
    # CoreFramework(normalize = True), 
    # OddSth(normalize = True), 
    # NeighborhoodHash(normalize = True, nh_type = 'count_sensitive'), 
    # PyramidMatch(normalize = True), 
    # Propagation(normalize = True, random_state = 42), 
    # HadamardCode(), 
    # NeighborhoodSubgraphPairwiseDistance(
    #     normalize = True
    # )
]:
    K_train = gk.fit_transform(G_train)
    K_test = gk.transform(G_test)
    # 处理缺失值
    imp_mean = SimpleImputer(missing_values = np.nan, strategy = 'mean')
    K_train = imp_mean.fit_transform(K_train)
    K_test = imp_mean.fit_transform(K_test)

    # Uses the SVM classifier to perform classification
    clf = SVC(kernel="rbf", probability = True)
    clf.fit(K_train, y_train)
    y_pred = clf.predict(K_test)
    y_proba = clf.predict_proba(K_test)

    # Computes and prints the classification accuracy
    acc = accuracy_score(y_test, y_pred)
    # roc_auc = roc_auc_score(y_test, y_proba[:, 1])
    roc_auc = roc_auc_score(y_test, y_proba, multi_class = 'ovo')
    print(f'acc = {acc}, roc_auc = {roc_auc}')
