#!/bin/bash

# 执行脚本的方法
# chmod +x ./test.sh
# ./test.sh

# 进入到 'node' 路径下，执行 node classification 任务
cd cd ../4-MVGRL-nc-gc/node

# Cora with full graph
python main.py --dataname cora --gpu 0
# Cora with sampled subgraphs
python main_sample.py --dataname cora --gpu 0
# Citeseer with full graph
python main.py --dataname citeseer --wd1 0.001 --wd2 0.01 --epochs 200 --gpu 0
# Citeseer with sampled subgraphs
python main_sample.py --dataname citeseer --wd2 0.01 --gpu 0

# Reddit with sampled subgraphs
python main_sample.py --dataname reddit --sample_size 4000 --epochs 400 --patience 999 --gpu 0
# CoraFull with sampled subgraphs
python main_sample.py --dataname cora_full --sample_size 4000 --epochs 400 --patience 999 --gpu 0
# Pubmed with sampled subgraphs
python main_sample.py --dataname pubmed --sample_size 4000 --epochs 400 --patience 999 --gpu 0
# AmazonCoBuyComputer with sampled subgraphs
python main_sample.py --dataname amazon_cobuy_computer --sample_size 4000 --epochs 400 --patience 999 --gpu 0
# MUTAG with sampled subgraphs
python main_sample.py --dataname mutag --sample_size 4000 --epochs 400 --patience 999 --gpu 0
# CoauthorCS with sampled subgraphs
python main_sample.py --dataname coauthor_cs --sample_size 4000 --epochs 400 --patience 999 --gpu 0
# PPI with sampled subgraphs
python main_sample.py --dataname ppi --epochs 400 --patience 999 --gpu 0