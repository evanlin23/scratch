AI-assisted (Claude), exploration code for CS581 project

SIMHIGH_R18 (helper t18). aln.jsonl: gg.py rows (linsi = MAGUS merge, wsoft0.03:linsi&fftns2#es4 = recipe,
linsi#es3, linsi&fftns2-op3 = hard, linsi#es4) followed by vote.py rows from run.py (B = 10).
vote.py `magus` reproduces MAGUS's merge: same SPFN/SPFP/TC/LenEst as the MAGUS draw and gg.py `linsi`
(out.fasta not byte-identical to gg.py's, same scores). vote.py `es4` scores 18.15 vs gg.py `linsi#es4` 18.09.
trees.jsonl: FastTree -lg -gamma (protbench trees.py), true / magus / es4 / recipe / es3 / hard, then vote_*;
vote_hard+mask uses vote/hard+mask/out.masked.fasta (5 of 8417 columns masked).
No tree-set alignment was byte-identical (cmp) to another, so every tree was computed, none reused.
