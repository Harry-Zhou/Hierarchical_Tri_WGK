#!/bin/bash

# 执行脚本的方法
# chmod +x ./test.sh
# ./test.sh

cd ../5-DGI-nc/
datasets = ('reddit', 'cora_full', 'pubmed', 'citeseer', 'amazon_cobuy_computer', 'mutag', 'coauthor_cs', 'ppi')

for data in ${datasets[*]}
do
    python train.py --dataset ${data} --gpu 0 --self-loop
done