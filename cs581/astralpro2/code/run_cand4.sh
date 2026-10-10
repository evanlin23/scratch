#!/bin/bash
P=/opt/mm/root/envs/gdl/bin/python
M=true,ovl,own,recon_true,recon_first,apro_bin,astral_multi,disco_astral,wqfm_gdl,duploss2
$P runmeth.py /opt/runs/ap2/cand4.nwk 1000 10 $M ../results/meth_cand4.jsonl
$P runmeth.py /opt/runs/ap2/cand4.nwk 5000 4 apro_bin,astral_multi,disco_astral,wqfm_gdl,duploss2 ../results/meth_cand4.jsonl
