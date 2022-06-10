#!/bin/bash

# 执行脚本的方法
# chmod +x ./test.sh
# ./test.sh

cd ../6-GPT-GNN-nc/example_reddit
# 数据预处理，预训练集与测试下游任务的数据集
python preprocess_reddit.py
# 预训练
python pretrain_reddit.py --conv_name hgt --n_layers 3 --pretrain_model_dir /datadrive/models/gta_all_cs3
# 迁移到下游任务中

datasets = ('reddit', 'cora_full', 'pubmed', 'citeseer', 'amazon_cobuy_computer', 'mutag', 'coauthor_cs', 'ppi')
for dname in ${datasets[*]}
do
    python finetune_reddit.py --dataset dname --use_pretrain --pretrain_model_dir /datadrive/models/gta_all_cs3 --n_layer 3 --data_percentage 0.1
done