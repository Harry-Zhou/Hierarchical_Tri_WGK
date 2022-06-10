import argparse

from loader import BioDataset
from dataloader import DataLoaderFinetune
from splitters import random_split, species_split

import torch
import torch.nn as nn
import torch.nn.functional as F
import torch.optim as optim
from torch_geometric.datasets import TUDataset

from tqdm import tqdm
import numpy as np

from model import GNN, GNN_graphpred
from sklearn.metrics import accuracy_score, f1_score, roc_auc_score

import pandas as pd

import os
import csv
import pickle

criterion = nn.BCEWithLogitsLoss()

def train(args, model, device, loader, optimizer):
    model.train()

    for step, batch in enumerate(tqdm(loader, desc="Iteration")):
        batch = batch.to(device)
        pred = model(batch)
        y = batch.go_target_downstream.view(pred.shape).to(torch.float64)

        optimizer.zero_grad()
        loss = criterion(pred.double(), y)
        loss.backward()

        optimizer.step()

def eval(args, model, device, loader):
    model.eval()
    y_true = []
    y_prob = []

    for _, batch in enumerate(tqdm(loader, desc="Iteration")):
        batch = batch.to(device)
        with torch.no_grad():
            pred = model(batch)
    
        #y_true.append(batch.go_target_downstream.view(pred.shape).detach().cpu())
        y_true.append(batch.go_target_downstream.detach().cpu())
        y_prob.append(pred.detach().cpu())

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

    #roc_list = []
    #for i in range(y_true.shape[1]):
    #    #AUC is only defined when there is at least one positive data.
    #    if np.sum(y_true[:,i] == 1) > 0 and np.sum(y_true[:,i] == 0) > 0:
    #        roc_list.append(roc_auc_score(y_true[:,i], y_prob[:,i]))
    #    else:
    #        roc_list.append(np.nan)

    #return np.array(roc_list) #y_true.shape[1]
    return acc, f1, roc_auc

def main():
    # Training settings
    parser = argparse.ArgumentParser(description='PyTorch implementation of pre-training of graph neural networks')
    parser.add_argument('--device', type=int, default=0,
                        help='which gpu to use if any (default: 0)')
    parser.add_argument('--batch_size', type=int, default=32,
                        help='input batch size for training (default: 32)')
    parser.add_argument('--epochs', type=int, default=50,
                        help='number of epochs to train (default: 50)')
    parser.add_argument('--lr', type=float, default=0.001,
                        help='learning rate (default: 0.001)')
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
    parser.add_argument('--model_file', type=str, default = '', help='filename to read the model (if there is any)')
    parser.add_argument('--filename', type=str, default = '', help='output filename')
    parser.add_argument('--gnn_type', type=str, default="gin")
    parser.add_argument('--seed', type=int, default=42, help = "Seed for splitting dataset.")
    parser.add_argument('--runseed', type=int, default=0, help = "Seed for running experiments.")
    parser.add_argument('--num_workers', type=int, default = 0, help='number of workers for dataset loading')
    parser.add_argument('--eval_train', type=int, default = 0, help='evaluating training or not')
    #parser.add_argument('--split', type=str, default = "random", help='Random or species split')
    args = parser.parse_args()

    torch.manual_seed(args.runseed)
    np.random.seed(args.runseed)
    device = torch.device("cuda:" + str(args.device)) if torch.cuda.is_available() else torch.device("cpu")
    if torch.cuda.is_available():
        torch.cuda.manual_seed_all(args.runseed)

    fout = open('graphcl_eval_result.csv', 'w', newline = '', encoding = 'utf-8')
    csv_writer = csv.writer(fout)
    csv_writer.writerow(['dataset', 'Val Accuracy', 'Val Micro_F1', 'Val Roc_Auc', \
        'Test Accuracy', 'Test Micro_F1', 'Test Roc_Auc'])
    
    # root_supervised = 'dataset/supervised'
    # dataset = BioDataset(root_supervised, data_type='supervised')
    root_dir = os.path.join('.', 'datasets_finetune_gc')
    for dname in ['reddit_multi_12k', 'reddit_multi_5k', 'reddit_binary', \
        'github_stargazers', 'COLLAB', 'DD', 'PROTEINS_full', 'UACC257H', \
        'Mutagenicity', 'NCI-H23H', 'P388H', 'PC-3H', 'SN12CH']:
        dataset = TUDataset(root = root_dir, name = dname)
        print('当前数据集: {}'.format(dname))

        node_num = 0
        edge_num = 0
        for d in dataset:
            node_num += d.x.size()[0]
            edge_num += d.edge_index.size()[1]
        print(node_num / len(dataset))
        print(edge_num / len(dataset))

        train_dataset, valid_dataset, test_dataset = random_split(dataset, seed = args.seed) 
        #if args.split == "random":
        #    print("random splitting")
        #    train_dataset, valid_dataset, test_dataset = random_split(dataset, seed = args.seed)
        #elif args.split == "species":
        #    trainval_dataset, test_dataset = species_split(dataset)
        #    train_dataset, valid_dataset, _ = random_split(trainval_dataset, seed = args.seed, frac_train=0.85, frac_valid=0.15, frac_test=0)
        #    test_dataset_broad, test_dataset_none, _ = random_split(test_dataset, seed = args.seed, frac_train=0.5, frac_valid=0.5, frac_test=0)
        #    print("species splitting")
        #else:
        #    raise ValueError("Unknown split name.")

        # splitters.py
        train_dataset, valid_dataset, test_dataset = random_split(dataset, seed = args.seed) 
        # dataloader.py
        train_loader = DataLoaderFinetune(train_dataset, batch_size=args.batch_size, shuffle=True, num_workers = args.num_workers)
        test_loader = DataLoaderFinetune(test_dataset, batch_size=10*args.batch_size, shuffle=False, num_workers = args.num_workers)
        val_loader = DataLoaderFinetune(valid_dataset, batch_size=10*args.batch_size, shuffle=False, num_workers = args.num_workers)
    
        #if args.split == "random":
        #    test_loader = DataLoaderFinetune(test_dataset, batch_size=10*args.batch_size, shuffle=False, num_workers = args.num_workers)
        #else:
        #    ### for species splitting
        #    test_easy_loader = DataLoaderFinetune(test_dataset_broad, batch_size=10*args.batch_size, shuffle=False, num_workers = args.num_workers)
        #    test_hard_loader = DataLoaderFinetune(test_dataset_none, batch_size=10*args.batch_size, shuffle=False, num_workers = args.num_workers)

        num_tasks = len(dataset[0].go_target_downstream)

        print(train_dataset[0])

        #set up model
        model = GNN_graphpred(args.num_layer, args.emb_dim, num_tasks, JK = args.JK, drop_ratio = args.dropout_ratio, graph_pooling = args.graph_pooling, gnn_type = args.gnn_type)

        if not args.model_file == "":
            model.from_pretrained(args.model_file)
    
        model.to(device)

        #set up optimizer
        optimizer = optim.Adam(model.parameters(), lr=args.lr, weight_decay=args.decay)

        train_acc_list = []
        train_f1_list = []
        train_roc_auc_list = []

        val_acc_list = []
        val_f1_list = []
        val_roc_auc_list = []
    
        ### for random splitting
        test_acc_list = []
        test_f1_list = []
        test_roc_auc_list = []
    
        #### for species splitting
        #test_acc_easy_list = []
        #test_acc_hard_list = []

        if not args.filename == "":
            if os.path.exists(args.filename):
                print("removed existing file!!")
                os.remove(args.filename)

        for epoch in range(1, args.epochs+1):
            print("====epoch " + str(epoch))
        
            train(args, model, device, train_loader, optimizer)

            print("====Evaluation")
            if args.eval_train:
                train_acc, train_f1, train_roc_auc = eval(args, model, device, train_loader)
            else:
                train_acc, train_f1, train_roc_auc = 0, 0, 0
                print("ommitting training evaluation")
            val_acc, val_f1, val_roc_auc = eval(args, model, device, val_loader)

            train_acc_list.append(train_acc)
            train_f1_list.append(train_f1)
            train_roc_auc_list.append(train_roc_auc)

            val_acc_list.append(val_acc)
            val_f1_list.append(val_f1)
            val_roc_auc_list.append(val_roc_auc)

            test_acc, test_f1, test_roc_auc = eval(args, model, device, test_loader)
            test_acc_list.append(test_acc)
            test_f1_list.append(test_f1)
            test_roc_auc_list.append(test_roc_auc)

            #if args.split == "random":
            #    test_acc = eval(args, model, device, test_loader)
            #    test_acc_list.append(test_acc)
            #else:
            #    test_acc_easy = eval(args, model, device, test_easy_loader)
            #    test_acc_hard = eval(args, model, device, test_hard_loader)
            #    test_acc_easy_list.append(np.mean(test_acc_easy))
            #    test_acc_hard_list.append(np.mean(test_acc_hard))
            #    print(val_acc_list[-1])
            #    print(test_acc_easy_list[-1])
            #    print(test_acc_hard_list[-1])

            print("")
        best_val_roc_auc = np.array(val_roc_auc_list)
        best_val_roc_auc_idx = np.array(val_roc_auc_list).argmax()
        csv_writer.writerow([dname, \
            np.array(val_acc_list)[best_val_roc_auc_idx], \
            np.array(val_f1_list)[best_val_roc_auc_idx], \
            best_val_roc_auc, \
            np.array(test_acc_list)[best_val_roc_auc_idx], \
            np.array(test_f1_list)[best_val_roc_auc_idx], \
            np.array(test_roc_auc_list)[best_val_roc_auc_idx]])

        #with open('graphcl_eval_result.log', 'a+') as f:
        #    f.write(str(args.runseed) + ' ' + str(np.array(test_acc_list)[np.array(val_acc_list).argmax()]) + ' ' + str(np.array(test_acc_hard_list)[np.array(val_acc_list).argmax()]))
        #    f.write('\n')
    fout.close()

if __name__ == "__main__":
    main()
