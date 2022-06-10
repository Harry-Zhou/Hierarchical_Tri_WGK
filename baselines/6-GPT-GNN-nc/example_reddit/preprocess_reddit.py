from torch_geometric.datasets import Reddit, CoraFull, Planetoid, Amazon, Entities, Coauthor, PPI
from GPT_GNN.data import *
import os

def preprocess_data(dname):
    root = '/datadrive/dataset'
    if dname == 'reddit':
        dataset = Reddit(root = root)
    elif dname == 'cora_full':
        dataset = CoraFull(root = root)
    elif dname == 'pubmed':
        dataset = Planetoid(root = root, name = 'Pubmed')
    elif dname == 'citeseer':
        dataset = Planetoid(root = root, name = 'Citeseer')
    elif dname == 'amazon_cobuy_computer':
        dataset = Amazon(root = root, name = 'Computers')
    elif dname == 'mutag':
        dataset = Entities(root = root, name = 'MUTAG')
    elif dname == 'coauthor_cs':
        dataset = Coauthor(root = root, name = 'CS')
    elif dname == 'ppi':
        dataset = PPI(root = root)
    else:
        raise Exception("数据的数据集名称不在预定集合中")
    graph = Graph()
    el = defaultdict(  #target_id
                    lambda: defaultdict( #source_id(
                    lambda: int # time
                    ))
    for i, j in tqdm(dataset.data.edge_index.t()):
        el[i.item()][j.item()] = 1

    target_type = 'def'
    graph.edge_list['def']['def']['def'] = el
    n = list(el.keys())
    degree = np.zeros(np.max(n)+1)
    for i in n:
        degree[i] = len(el[i])
    x = np.concatenate((dataset.data.x.numpy(), np.log(degree).reshape(-1, 1)), axis=-1)
    graph.node_feature['def'] = pd.DataFrame({'emb': list(x)})

    idx = np.arange(len(graph.node_feature[target_type]))
    np.random.seed(43)
    np.random.shuffle(idx)

    graph.pre_target_nodes   = idx[ : int(len(idx) * 0.7)]
    graph.train_target_nodes = idx[int(len(idx) * 0.7) : int(len(idx) * 0.8)]
    graph.valid_target_nodes = idx[int(len(idx) * 0.8) : int(len(idx) * 0.9)]
    graph.test_target_nodes  = idx[int(len(idx) * 0.9) : ]

    graph.y = dataset.data.y
    dill.dump(graph, open(os.path.join('.', 'preprocessed_dataset', 'graph_'+dname+'.pk'), 'wb'))

if __name__ == '__main__':
    for dname in ('reddit', 'cora_full', 'pubmed', 'citeseer', 'amazon_cobuy_computer', 'mutag', 'coauthor_cs', 'ppi'):
        preprocess_data(dname)