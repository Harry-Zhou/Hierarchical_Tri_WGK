#!/bin/bash

# 执行脚本的方法
# chmod +x ./test.sh
# ./test.sh

cd ../0-GraphCL-gc/transferLearning_MoleculeNet_PPI/bio

### Using Pretrained GIN
for runseed in 2 4 6 8 10
do
    python finetune.py --model_file ./models_graphcl/graphcl_80.pth --epochs 50 --device 0 --runseed $runseed --gnn_type gin --lr 1e-3
done
