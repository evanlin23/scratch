#!/bin/bash
# restartable; each pool in its own process
M=true,ovl,own,recon_true,recon_first,gtp_dup,gtp_dl,apro_bin,astral_multi,disco_astral,wqfm_gdl,duploss2
P=/opt/mm/root/envs/gdl/bin/python
for p in "$@"; do
  ( for K in 1000 5000; do NB=$([ $K = 1000 ] && echo 10 || echo 4)
      $P runmeth.py /opt/runs/ap2/$p.nwk $K $NB $M ../results/meth_$p.jsonl; done ) > /opt/runs/ap2/log_$p.txt 2>&1 &
done
wait
