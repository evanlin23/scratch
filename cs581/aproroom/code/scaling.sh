#!/bin/bash
# Q3 scaling: one job at a time on an otherwise idle machine. Restartable.
H=$(cd "$(dirname "$0")" && pwd); B="/opt/mm/root/envs/gdl/bin/python $H/bench.py"; O=$H/../results/scaling.jsonl
export BENCH_TMP=/opt/tmp
for k in 100 200 500 1000; do
  D=/opt/data/scaling/taxa_$k; K="{\"data\":\"scaling\",\"cond\":\"taxa_$k\",\"rep\":\"01\",\"ngen\":1000,\"level\":\"est100\",\"thr\":THR}"
  $B --genes $D/genes.trees --true $D/s_tree.trees --out $O --key "${K/THR/1}" --methods astrid-pro,astrid-disco --threads 1
  $B --genes $D/genes.trees --true $D/s_tree.trees --out $O --key "${K/THR/4}" --methods astrid-pro --threads 4
  (ulimit -v 12000000; BENCH_TIMEOUT=1200 $B --genes $D/genes.trees --true $D/s_tree.trees --out $O --key "${K/THR/1}" --methods asteroid --threads 1)
  [ $k -le 500 ] && BENCH_TIMEOUT=1500 $B --genes $D/genes.trees --true $D/s_tree.trees --out $O --key "${K/THR/4}" --methods astral-pro3 --threads 4
done
D=/opt/data/disco/trees/gtrees_10000_l1/01
for n in 100 1000 10000; do
  K="{\"data\":\"scaling\",\"cond\":\"genes\",\"rep\":\"01\",\"ngen\":$n,\"level\":\"est100\",\"thr\":THR}"
  $B --genes $D/g_100.trees --ngen $n --true $D/s_tree.trees --out $O --key "${K/THR/1}" --methods astrid-pro,astrid-disco --threads 1
  $B --genes $D/g_100.trees --ngen $n --true $D/s_tree.trees --out $O --key "${K/THR/4}" --methods astrid-pro --threads 4
  [ $n -le 1000 ] && (ulimit -v 12000000; BENCH_TIMEOUT=1200 $B --genes $D/g_100.trees --ngen $n --true $D/s_tree.trees --out $O --key "${K/THR/1}" --methods asteroid --threads 1)
  [ $n -le 1000 ] && BENCH_TIMEOUT=1500 $B --genes $D/g_100.trees --ngen $n --true $D/s_tree.trees --out $O --key "${K/THR/4}" --methods astral-pro3 --threads 4
  [ $n -le 100 ] && BENCH_TIMEOUT=1500 $B --genes $D/g_100.trees --ngen $n --true $D/s_tree.trees --out $O --key "${K/THR/4}" --methods wqfm-gdl --threads 4
done
