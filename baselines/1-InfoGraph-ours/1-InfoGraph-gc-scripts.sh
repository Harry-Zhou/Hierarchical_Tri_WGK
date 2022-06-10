#!/bin/bash

# 执行脚本的方法
# chmod +x ./test.sh
# ./test.sh

cd ../1-InfoGraph-gc/unsupervised

datasets = ('reddit_multi_12k', 'reddit_multi_5k', 'reddit_binary', 'github_stargazers', 'COLLAB', 'DD', 'PROTEINS_full', 'UACC257H', 'Mutagenicity', 'NCI-H23H', 'P388H', 'PC-3H', 'SN12CH')
for dname in ${datasets[*]}
do
    python main.py --DS ${dname} --lr 0.001 --num-gc-layers 3
done