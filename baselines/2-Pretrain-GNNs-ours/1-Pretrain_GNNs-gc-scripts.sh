#!/bin/bash

# 执行脚本的方法
# chmod +x ./test.sh
# ./test.sh

cd ../2-Pretrain-GNNs-gc/chem

datasets = ('reddit_multi_12k', 'reddit_multi_5k', 'reddit_binary', 'github_stargazers', 'COLLAB', 'DD', 'PROTEINS_full', 'UACC257H', 'Mutagenicity', 'NCI-H23H', 'P388H', 'PC-3H', 'SN12CH')
for dname in ${datasets[*]}
do
    python finetune.py --input_model_file ./model_gin/${model_file}.pth --split random --device 0 --runseed 1 --gnn_type gin --dataset ${dname}
done
