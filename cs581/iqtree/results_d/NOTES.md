AI-assisted (Claude), exploration code for CS581 project

# iqtree results_d: magus reproduction check

`magus` rebuilt with gcmvote/code/vote.py vs MAGUS's recorded score in the bank (magus.json):

| rep | rebuilt avgErr | bank avgErr | rebuilt SPFN / SPFP | bank SPFN / SPFP | rebuilt / bank LenEst |
|---|---|---|---|---|---|
| SIMHIGH_R7 | 0.27277 | 0.27277 | 0.34404 / 0.20150 | 0.34404 / 0.20150 | 9902 / 9902 |
| SIMHIGH_R8 | 0.26186 | 0.26252 | 0.36195 / 0.16177 | 0.36501 / 0.16004 | 9854 / 9852 |

R7 reproduces exactly. R8 differs by 0.0007 avgErr (2 columns in length); not investigated.
IQ-TREE 3.1.4 (bioconda), `-m LG+G4 --fast -T 2 -seed 1`, 2 runs at a time on 4 cores.
