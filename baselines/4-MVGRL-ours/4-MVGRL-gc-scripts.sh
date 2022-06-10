#!/bin/bash

# 执行脚本的方法
# chmod +x ./test.sh
# ./test.sh

# 进入到 'graph' 路径下，执行 graph classification 任务
cd ../4-MVGRL-nc-gc/graph

# reddit_multi_12k
python main.py --dataname reddit_multi_12k --gpu 0 --hid_dim 64
# reddit_multi_5k
python main.py --dataname reddit_multi_5k --gpu 0 --hid_dim 64
# reddit_binary
python main.py --dataname reddit_binary --gpu 0 --hid_dim 64
# github_stargazers
python main.py --dataname github_stargazers --gpu 0 --hid_dim 64
# COLLAB
python main.py --dataname COLLAB --gpu 0 --hid_dim 64
# DD
python main.py --dataname DD --gpu 0 --hid_dim 64
# PROTEINS_full
python main.py --dataname PROTEINS_full --gpu 0 --hid_dim 64
# UACC257H
python main.py --dataname UACC257H --gpu 0 --hid_dim 64
# Mutagenicity
python main.py --dataname Mutagenicity --gpu 0 --hid_dim 64
# NCI-H23H
python main.py --dataname NCI-H23H --gpu 0 --hid_dim 64
# P388H
python main.py --dataname P388H --gpu 0 --hid_dim 64
# PC-3H
python main.py --dataname PC-3H --gpu 0 --hid_dim 64
# SN12CH
python main.py --dataname SN12CH --gpu 0 --hid_dim 64