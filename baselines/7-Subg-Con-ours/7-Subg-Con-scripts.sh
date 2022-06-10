#!/bin/bash

# 执行脚本的方法
# chmod +x ./test.sh
# ./test.sh

cd ../7-Subg-Con-nc/

datasets = ('reddit', 'cora_full', 'pubmed', 'citeseer', 'amazon_cobuy_computer', 'mutag', 'coauthor_cs', 'ppi')
for dname in ${datasets[*]}
do
    python train.py --dataset dname
done