#!/bin/bash

# 执行脚本的方法
# chmod +x ./test.sh
# ./test.sh
cd ../3-GCC-nc-gc/Scripts/
# Transferring Mode: Freeze
# <gpu> <load-path> <dataset>
bash generate.sh
    0
    saved/Pretrain_moco_False_dgl_gin_layer_5_lr_0.005_decay_1e-05_bsz_256_hid_64_samples_2000_nce_t_0.07_nce_k_32_rw_hops_256_restart_prob_0.8_aug_1st_ft_False_deg_16_pos_32_momentum_0.999
    reddit_multi_12k, reddit_multi_5k, reddit_binary, github_stargazers, collab, dd, proteins_full, uacc257h, mutagenicity, nci-h23h, p388h, pc-3h, sn12ch
# Evaluate
# <load_path> <hidden_size> <dataset>
bash graph_classification/ours.sh
    saved/Pretrain_moco_False_dgl_gin_layer_5_lr_0.005_decay_1e-05_bsz_256_hid_64_samples_2000_nce_t_0.07_nce_k_32_rw_hops_256_restart_prob_0.8_aug_1st_ft_False_deg_16_pos_32_momentum_0.999
    64
    freeze
    reddit_multi_12k, reddit_multi_5k, reddit_binary, github_stargazers, collab, dd, proteins_full, uacc257h, mutagenicity, nci-h23h, p388h, pc-3h, sn12ch

# Transferring Mode: Fine-Tune
# <load-path> <gpu> <dataset>
bash finetune.sh
    saved/Pretrain_moco_True_dgl_gin_layer_5_lr_0.005_decay_1e-05_bsz_32_hid_64_samples_2000_nce_t_0.07_nce_k_16384_rw_hops_256_restart_prob_0.8_aug_1st_ft_False_deg_16_pos_32_momentum_0.999
    0
    reddit_multi_12k, reddit_multi_5k, reddit_binary, github_stargazers, collab, dd, proteins_full, uacc257h, mutagenicity, nci-h23h, p388h, pc-3h, sn12ch
# Evaluate
# <load_path> <hidden_size> <dataset>
bash graph_classification/ours.sh
    saved/Pretrain_moco_True_dgl_gin_layer_5_lr_0.005_decay_1e-05_bsz_32_hid_64_samples_2000_nce_t_0.07_nce_k_16384_rw_hops_256_restart_prob_0.8_aug_1st_ft_False_deg_16_pos_32_momentum_0.999
    64
    fine-tune
    reddit_multi_12k, reddit_multi_5k, reddit_binary, github_stargazers, collab, dd, proteins_full, uacc257h, mutagenicity, nci-h23h, p388h, pc-3h, sn12ch
