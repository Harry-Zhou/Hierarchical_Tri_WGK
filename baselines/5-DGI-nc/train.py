import argparse, time
import numpy as np
import networkx as nx
import torch
import torch.nn as nn
import torch.nn.functional as F
from dgl import DGLGraph
from dgl.data import register_data_args, load_data
from dgi import DGI, Classifier
import dgl
import os.path
from sklearn.metrics import accuracy_score, f1_score, roc_auc_score
import csv

def evaluate(model, features, labels, mask):
    model.eval()
    num_cls = labels.max() + 1
    with torch.no_grad():
        logits = model(features)
        logits = logits[mask].detach().cpu().numpy()
        labels = labels[mask].detach().cpu().numpy()
        preds = logits.argmax(axis = 1)
        # _, indices = torch.max(logits, dim=1)

        acc = accuracy_score(labels, preds)
        if num_cls > 2:
            f1 = f1_score(labels, preds, average = 'micro')
            roc_auc = roc_auc_score(labels, preds, multi_class = 'ovr')
        else:
            f1 = f1_score(labels, preds)
            roc_auc = roc_auc_score(labels, preds)

        # correct = torch.sum(indices == labels)
        return acc, f1, roc_auc

def main(args):
    raw_dir = os.path.join('..', '..', 'datasets_finetune_nc') # 与 DualGNN 的测试数据保持一致
    dgl_reddit = dgl.data.RedditDataset(raw_dir = raw_dir)
    dgl_cora_full = dgl.data.CoraFullDataset(raw_dir = raw_dir)
    dgl_cora = dgl.data.CoraGraphDataset(raw_dir = raw_dir)
    dgl_pubmed = dgl.data.PubmedGraphDataset(raw_dir = raw_dir)
    dgl_citeseer = dgl.data.CiteseerGraphDataset(raw_dir = raw_dir)
    dgl_amazon_cobuy_computer = dgl.data.AmazonCoBuyComputerDataset(raw_dir = raw_dir)
    dgl_mutag = dgl.data.MUTAGDataset(raw_dir = raw_dir)
    dgl_coauthor_cs = dgl.data.CoauthorCSDataset(raw_dir = raw_dir)
    dgl_ppi = dgl.data.PPIDataset(raw_dir = raw_dir)

    datasets = {'reddit': dgl_reddit, \
        'cora_full': dgl_cora_full, \
        'cora': dgl_cora, \
        'pubmed': dgl_pubmed, \
        'citeseer': dgl_citeseer, \
        'amazon_cobuy_computer': dgl_amazon_cobuy_computer, \
        'mutag': dgl_mutag, \
        'coauthor_cs': dgl_coauthor_cs, \
        'ppi': dgl_ppi}
    fout = open(os.path.join('.', 'results.csv'), 'w', newline = '', encoding = 'utf-8')
    csv_writer = csv.writer(fout)
    csv_writer.writerow(['Dataset', \
        'Accuracy', \
        'Micro_F1', \
        'Roc_Auc'])
    for dname, data in datasets.items():
        data = datasets[args.dataset]
        g = data[0]
        features = g.ndata['feat']
        in_feats = features.shape[1]
        labels = g.ndata['label']
        train_mask = g.ndata['train_mask']
        val_mask = g.ndata['val_mask']
        test_mask = g.ndata['test_mask']
        n_classes = data.num_classes
        n_edges = g.num_edges()
        if args.self_loop:
            g.remove_self_loop()
            g.add_self_loop()
    
        # # load and preprocess dataset
        # data = load_data(args)
        # features = torch.FloatTensor(data.features)
        # labels = torch.LongTensor(data.labels)
        # if hasattr(torch, 'BoolTensor'):
        #     train_mask = torch.BoolTensor(data.train_mask)
        #     val_mask = torch.BoolTensor(data.val_mask)
        #     test_mask = torch.BoolTensor(data.test_mask)
        # else:
        #     train_mask = torch.ByteTensor(data.train_mask)
        #     val_mask = torch.ByteTensor(data.val_mask)
        #     test_mask = torch.ByteTensor(data.test_mask)
        # in_feats = features.shape[1]
        # n_classes = data.num_labels
        # n_edges = data.graph.number_of_edges()

        if args.gpu < 0:
            cuda = False
        else:
            cuda = True
            torch.cuda.set_device(args.gpu)
            features = features.cuda()
            labels = labels.cuda()
            train_mask = train_mask.cuda()
            val_mask = val_mask.cuda()
            test_mask = test_mask.cuda()

        # # graph preprocess
        # g = data.graph
        # # add self loop
        # if args.self_loop:
        #     g.remove_edges_from(nx.selfloop_edges(g))
        #     g.add_edges_from(zip(g.nodes(), g.nodes()))
        # g = DGLGraph(g)
        # n_edges = g.number_of_edges()

        if args.gpu >= 0:
            g = g.to(args.gpu)
        # create DGI model
        dgi = DGI(g,
                  in_feats,
                  args.n_hidden,
                  args.n_layers,
                  nn.PReLU(args.n_hidden),
                  args.dropout)

        if cuda:
            dgi.cuda()

        dgi_optimizer = torch.optim.Adam(dgi.parameters(),
                                         lr=args.dgi_lr,
                                         weight_decay=args.weight_decay)

        # train deep graph infomax
        cnt_wait = 0
        best = 1e9
        best_t = 0
        dur = []
        for epoch in range(args.n_dgi_epochs):
            dgi.train()
            if epoch >= 3:
                t0 = time.time()

            dgi_optimizer.zero_grad()
            loss = dgi(features)
            loss.backward()
            dgi_optimizer.step()

            if loss < best:
                best = loss
                best_t = epoch
                cnt_wait = 0
                torch.save(dgi.state_dict(), 'best_dgi.pkl')
            else:
                cnt_wait += 1

            if cnt_wait == args.patience:
                print('Early stopping!')
                break

            if epoch >= 3:
                dur.append(time.time() - t0)

            print("Epoch {:05d} | Time(s) {:.4f} | Loss {:.4f} | "
                    "ETputs(KTEPS) {:.2f}".format(epoch, np.mean(dur), \
                    loss.item(), n_edges / np.mean(dur) / 1000))

        # create classifier model
        classifier = Classifier(args.n_hidden, n_classes)
        if cuda:
            classifier.cuda()

        classifier_optimizer = torch.optim.Adam(classifier.parameters(),
                                                lr=args.classifier_lr,
                                                weight_decay=args.weight_decay)

        # train classifier
        print('Loading {}th epoch'.format(best_t))
        dgi.load_state_dict(torch.load('best_dgi.pkl'))
        embeds = dgi.encoder(features, corrupt=False)
        embeds = embeds.detach()
        dur = []
        for epoch in range(args.n_classifier_epochs):
            classifier.train()
            if epoch >= 3:
                t0 = time.time()

            classifier_optimizer.zero_grad()
            preds = classifier(embeds)
            loss = F.nll_loss(preds[train_mask], labels[train_mask])
            loss.backward()
            classifier_optimizer.step()
        
            if epoch >= 3:
                dur.append(time.time() - t0)

            acc, f1, roc_auc = evaluate(classifier, embeds, labels, val_mask)
            print("Epoch {:05d} | Time(s) {:.4f} | Loss {:.4f} | Val Accuracy {:.4f} | Val Micro_F1 {:.4f} | Val Roc_Auc {:.4f}"
                "ETputs(KTEPS) {:.2f}".format(epoch, np.mean(dur), loss.item(), \
                acc, f1, roc_auc, n_edges / np.mean(dur) / 1000))

        print()
        acc, f1, roc = evaluate(classifier, embeds, labels, test_mask)
        csv_writer.writerow([dname, \
            acc, f1, roc_auc])
        # print("Test Accuracy {:.4f}".format(acc))
    fout.close()

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='DGI')
    register_data_args(parser)
    parser.add_argument("--dropout", type=float, default=0.,
                        help="dropout probability")
    parser.add_argument("--gpu", type=int, default=-1,
                        help="gpu")
    parser.add_argument("--dgi-lr", type=float, default=1e-3,
                        help="dgi learning rate")
    parser.add_argument("--classifier-lr", type=float, default=1e-2,
                        help="classifier learning rate")
    parser.add_argument("--n-dgi-epochs", type=int, default=300,
                        help="number of training epochs")
    parser.add_argument("--n-classifier-epochs", type=int, default=300,
                        help="number of training epochs")
    parser.add_argument("--n-hidden", type=int, default=512,
                        help="number of hidden gcn units")
    parser.add_argument("--n-layers", type=int, default=1,
                        help="number of hidden gcn layers")
    parser.add_argument("--weight-decay", type=float, default=0.,
                        help="Weight for L2 loss")
    parser.add_argument("--patience", type=int, default=20,
                        help="early stop patience condition")
    parser.add_argument("--self-loop", action='store_true',
                        help="graph self-loop (default=False)")
    parser.set_defaults(self_loop=False)
    args = parser.parse_args()
    print(args)
    
    main(args)
