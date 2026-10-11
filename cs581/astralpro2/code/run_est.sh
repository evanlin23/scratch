#!/bin/bash
P=/opt/mm/root/envs/gdl/bin/python
for r in 0 1 2 3 4; do $P est.py Rown $r 1000 ../results/est.jsonl; done
for r in 0 1 2; do $P est.py S25 $r 1000 ../results/est.jsonl; $P est.py S25het $r 1000 ../results/est.jsonl; done
