import argparse

from loader import MoleculeDataset
from torch_geometric.data import DataLoader
from torch_geometric.datasets import TUDataset

import torch
import torch.nn as nn
import torch.nn.functional as F
import torch.optim as optim

from tqdm import tqdm
import numpy as np

from model import GNN, GNN_graphpred
from sklearn.metrics import accuracy_score, f1_score, roc_auc_score

from splitters import random_split
import pandas as pd

import os
import shutil
import csv

from tensorboardX import SummaryWriter

criterion = nn.BCEWithLogitsLoss(reduction = "none")

def train(args, model, device, loader, optimizer):
    model.train()

    for step, batch in enumerate(tqdm(loader, desc="Iteration")):
        batch = batch.to(device)
        pred = model(batch.x, batch.edge_index, batch.edge_attr, batch.batch)
        y = batch.y.view(pred.shape).to(torch.float64)

        #Whether y is non-null or not.
        is_valid = y**2 > 0
        #Loss matrix
        loss_mat = criterion(pred.double(), (y+1)/2)
        #loss matrix after removing null target
        loss_mat = torch.where(is_valid, loss_mat, torch.zeros(loss_mat.shape).to(loss_mat.device).to(loss_mat.dtype))
            
        optimizer.zero_grad()
        loss = torch.sum(loss_mat)/torch.sum(is_valid)
        loss.backward()

        optimizer.step()

def eval(args, model, device, loader):
    model.eval()
    y_true = []
    # y_scores = []
    y_prob = []

    for _, batch in enumerate(tqdm(loader, desc="Iteration")):
        batch = batch.to(device)

        with torch.no_grad():
            pred = model(batch.x, batch.edge_index, batch.edge_attr, batch.batch)

        # y_true.append(batch.y.view(pred.shape))
        # y_scores.append(pred)
        y_true.append(batch.go_target_downstream.detach().cpu())
        y_prob.append(pred.detach().cpu())

    # y_true = torch.cat(y_true, dim = 0).cpu().numpy()
    # y_scores = torch.cat(y_scores, dim = 0).cpu().numpy()
    y_true = torch.cat(y_true, dim = 0).numpy()
    y_prob = torch.cat(y_prob, dim = 0).numpy()
    y_pred = y_prob.argmax(axis = 1)
    
    acc = accuracy_score(y_true, y_pred)
    num_cls = y_prob.shape[1]
    if num_cls > 2:
        roc_auc = roc_auc_score(y_true, y_prob, multi_class = 'ovr')
        f1 = f1_score(y_true, y_pred, average = 'micro') # micro, weighted
    else:
        roc_auc = roc_auc_score(y_true, y_prob[:, 1])
        f1 = f1_score(y_true, y_pred)

    # roc_list = []
    # for i in range(y_true.shape[1]):
    #     #AUC is only defined when there is at least one positive data.
    #     if np.sum(y_true[:,i] == 1) > 0 and np.sum(y_true[:,i] == -1) > 0:
    #         is_valid = y_true[:,i]**2 > 0
    #         roc_list.append(roc_auc_score((y_true[is_valid,i] + 1)/2, y_scores[is_valid,i]))

    # if len(roc_list) < y_true.shape[1]:
    #     print("Some target is missing!")
    #     print("Missing ratio: %f" %(1 - float(len(roc_list))/y_true.shape[1]))

    # return sum(roc_list)/len(roc_list) #y_true.shape[1]
    return acc, f1, roc_auc

def main(csv_writer):
    # Training settings
    parser = argparse.ArgumentParser(description='PyTorch implementation of pre-training of graph neural networks')
    parser.add_argument('--device', type=int, default=0,
                        help='which gpu to use if any (default: 0)')
    parser.add_argument('--batch_size', type=int, default=32,
                        help='input batch size for training (default: 32)')
    parser.add_argument('--epochs', type=int, default=100,
                        help='number of epochs to train (default: 100)')
    parser.add_argument('--lr', type=float, default=0.001,
                        help='learning rate (default: 0.001)')
    parser.add_argument('--lr_scale', type=float, default=1,
                        help='relative learning rate for the feature extraction layer (default: 1)')
    parser.add_argument('--decay', type=float, default=0,
                        help='weight decay (default: 0)')
    parser.add_argument('--num_layer', type=int, default=5,
                        help='number of GNN message passing layers (default: 5).')
    parser.add_argument('--emb_dim', type=int, default=300,
                        help='embedding dimensions (default: 300)')
    parser.add_argument('--dropout_ratio', type=float, default=0.5,
                        help='dropout ratio (default: 0.5)')
    parser.add_argument('--graph_pooling', type=str, default="mean",
                        help='graph level pooling (sum, mean, max, set2set, attention)')
    parser.add_argument('--JK', type=str, default="last",
                        help='how the node features across layers are combined. last, sum, max or concat')
    parser.add_argument('--gnn_type', type=str, default="gin")
    parser.add_argument('--dataset', type=str, default = 'tox21', help='root directory of dataset. For now, only classification.')
    parser.add_argument('--input_model_file', type=str, default = '', help='filename to read the model (if there is any)')
    # parser.add_argument('--filename', type=str, default = '', help='output filename')
    parser.add_argument('--seed', type=int, default=42, help = "Seed for splitting the dataset.")
    parser.add_argument('--runseed', type=int, default=0, help = "Seed for minibatch selection, random initialization.")
    parser.add_argument('--split', type = str, default="scaffold", help = "random or scaffold or random_scaffold")
    parser.add_argument('--eval_train', type=int, default = 0, help='evaluating training or not')
    parser.add_argument('--num_workers', type=int, default = 4, help='number of workers for dataset loading')
    args = parser.parse_args()

    torch.manual_seed(args.runseed)
    np.random.seed(args.runseed)
    device = torch.device("cuda:" + str(args.device)) if torch.cuda.is_available() else torch.device("cpu")
    if torch.cuda.is_available():
        torch.cuda.manual_seed_all(args.runseed)

    raw_dir = os.path.join('..', 'DualGNN4gc_datasets_finetune_gc')
    tu_dataset = {'reddit_multi_12k': TUDataset('reddit_multi_12k', raw_dir = raw_dir), \
        'reddit_multi_5k': TUDataset('reddit_multi_5k', raw_dir = raw_dir), \
        'reddit_binary': TUDataset('reddit_binary', raw_dir = raw_dir), \
        'github_stargazers': TUDataset('github_stargazers', raw_dir = raw_dir), \
        'collab': TUDataset('COLLAB', raw_dir = raw_dir), \
        'dd': TUDataset('DD', raw_dir = raw_dir), \
        'proteins_full': TUDataset('PROTEINS_full', raw_dir = raw_dir), \
        'uacc257h': TUDataset('UACC257H', raw_dir = raw_dir), \
        'mutagenicity': TUDataset('Mutagenicity', raw_dir = raw_dir), \
        'nci_h23h': TUDataset('NCI-H23H', raw_dir = raw_dir), \
        'p388h': TUDataset('P388H', raw_dir = raw_dir), \
        'pc_3h': TUDataset('PC-3H', raw_dir = raw_dir), \
        'sn12ch': TUDataset('SN12CH', raw_dir = raw_dir)}
    #Bunch of classification tasks
    # if args.dataset == "tox21":
    #     num_tasks = 12
    # elif args.dataset == "hiv":
    #     num_tasks = 1
    # elif args.dataset == "pcba":
    #     num_tasks = 128
    # elif args.dataset == "muv":
    #     num_tasks = 17
    # elif args.dataset == "bace":
    #     num_tasks = 1
    # elif args.dataset == "bbbp":
    #     num_tasks = 1
    # elif args.dataset == "toxcast":
    #     num_tasks = 617
    # elif args.dataset == "sider":
    #     num_tasks = 27
    # elif args.dataset == "clintox":
    #     num_tasks = 2
    # else:
    #     raise ValueError("Invalid dataset name.")
    
    num_tasks = 1
    dataset = tu_dataset[args.dataset]
    #set up dataset
    #dataset = MoleculeDataset("dataset/" + args.dataset, dataset=args.dataset)

    print(dataset)
    train_dataset, valid_dataset, test_dataset = random_split(dataset, null_value=0, frac_train=0.8,frac_valid=0.1, frac_test=0.1, seed = args.seed)
    # if args.split == "scaffold":
    #     smiles_list = pd.read_csv('dataset/' + args.dataset + '/processed/smiles.csv', header=None)[0].tolist()
    #     train_dataset, valid_dataset, test_dataset = scaffold_split(dataset, smiles_list, null_value=0, frac_train=0.8,frac_valid=0.1, frac_test=0.1)
    #     print("scaffold")
    # elif args.split == "random":
    #     train_dataset, valid_dataset, test_dataset = random_split(dataset, null_value=0, frac_train=0.8,frac_valid=0.1, frac_test=0.1, seed = args.seed)
    #     print("random")
    # elif args.split == "random_scaffold":
    #     smiles_list = pd.read_csv('dataset/' + args.dataset + '/processed/smiles.csv', header=None)[0].tolist()
    #     train_dataset, valid_dataset, test_dataset = random_scaffold_split(dataset, smiles_list, null_value=0, frac_train=0.8,frac_valid=0.1, frac_test=0.1, seed = args.seed)
    #     print("random scaffold")
    # else:
    #     raise ValueError("Invalid split option.")

    print(train_dataset[0])

    train_loader = DataLoader(train_dataset, batch_size=args.batch_size, shuffle=True, num_workers = args.num_workers)
    val_loader = DataLoader(valid_dataset, batch_size=args.batch_size, shuffle=False, num_workers = args.num_workers)
    test_loader = DataLoader(test_dataset, batch_size=args.batch_size, shuffle=False, num_workers = args.num_workers)

    #set up model
    model = GNN_graphpred(args.num_layer, args.emb_dim, num_tasks, JK = args.JK, drop_ratio = args.dropout_ratio, graph_pooling = args.graph_pooling, gnn_type = args.gnn_type)
    if not args.input_model_file == "":
        model.from_pretrained(args.input_model_file)
    
    model.to(device)

    #set up optimizer
    #different learning rate for different part of GNN
    model_param_group = []
    model_param_group.append({"params": model.gnn.parameters()})
    if args.graph_pooling == "attention":
        model_param_group.append({"params": model.pool.parameters(), "lr":args.lr*args.lr_scale})
    model_param_group.append({"params": model.graph_pred_linear.parameters(), "lr":args.lr*args.lr_scale})
    optimizer = optim.Adam(model_param_group, lr=args.lr, weight_decay=args.decay)
    print(optimizer)

    train_acc_list, train_f1_list, train_roc_auc_list = [], [], []
    val_acc_list, val_f1_list, val_roc_auc_list = [], [], []
    test_acc_list, test_f1_list, test_roc_auc_list = [], [], []
    
    # if not args.filename == "":
    #     fname = 'runs/finetune_cls_runseed' + str(args.runseed) + '/' + args.filename
    #     #delete the directory if there exists one
    #     if os.path.exists(fname):
    #         shutil.rmtree(fname)
    #         print("removed the existing file.")
    #     writer = SummaryWriter(fname)

    for epoch in range(1, args.epochs+1):
        print("====epoch " + str(epoch))
        
        train(args, model, device, train_loader, optimizer)

        print("====Evaluation")
        if args.eval_train:
            train_acc, train_f1, train_roc_auc = eval(args, model, device, train_loader)
        else:
            print("omit the training accuracy computation")
            train_acc, train_f1, train_roc_auc = 0, 0, 0
        val_acc, val_f1, val_roc_auc = eval(args, model, device, val_loader)
        test_acc, test_f1, test_roc_auc = eval(args, model, device, test_loader)

        print("train: %f val: %f test: %f" %(train_acc, val_acc, test_acc))

        # val_acc_list.append(val_acc)
        # test_acc_list.append(test_acc)
        # train_acc_list.append(train_acc)

        train_acc_list.append(train_acc)
        train_f1_list.append(train_f1)
        train_roc_auc_list.append(train_roc_auc)

        val_acc_list.append(val_acc)
        val_f1_list.append(val_f1)
        val_roc_auc_list.append(val_roc_auc)

        test_acc_list.append(test_acc)
        test_f1_list.append(test_f1)
        test_roc_auc_list.append(test_roc_auc)

        # if not args.filename == "":
        #     writer.add_scalar('data/train auc', train_acc, epoch)
        #     writer.add_scalar('data/val auc', val_acc, epoch)
        #     writer.add_scalar('data/test auc', test_acc, epoch)

        print("")
    best_val_roc_auc = np.array(val_roc_auc_list)
    best_val_roc_auc_idx = np.array(val_roc_auc_list).argmax()
    csv_writer.writerow([args.dataset, \
        np.array(val_acc_list)[best_val_roc_auc_idx], \
        np.array(val_f1_list)[best_val_roc_auc_idx], \
        best_val_roc_auc, \
        np.array(test_acc_list)[best_val_roc_auc_idx], \
        np.array(test_f1_list)[best_val_roc_auc_idx], \
        np.array(test_roc_auc_list)[best_val_roc_auc_idx]])
    
    # if not args.filename == "":
    #     writer.close()

if __name__ == "__main__":
    with open(os.path.join('.', 'results.csv'), 'w', newline = '', encoding = 'utf-8') as fout:
        csv_writer = csv.writer(fout)
        csv_writer.writerow(['Dataset', \
            'Val Accuracy', \
            'Val Micro_F1', \
            'Val Roc_Auc', \
            'Test Accuracy', \
            'Test Micro_F1', \
            'Test Roc_Auc'])
        main(csv_writer)
