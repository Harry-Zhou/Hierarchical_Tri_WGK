#!/bin/bash

# 执行脚本的方法
# chmod +x ./test.sh
# ./test.sh

lenRwTraceList=(8)
numRwsList=(10 20 30 40 50)

for lenRwTrace in ${lenRwTraceList[@]}
do
    for numRws in ${numRwsList[@]}
    do
        python test_wgk.py --lenRwTrace $lenRwTrace --numRws $numRws
    done
done
