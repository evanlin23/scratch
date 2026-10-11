# DISCO-R pilot code

* `discor.py` – DL-parsimony rooting of gene-family trees against a rooted species tree (all rootings, O(n) rerooting DP, NOTUNG-style costs dup 1.5 / loss 1; ties broken by DISCO's MinDL score), overlap (DISCO) or LCA tagging, DISCO v1.4.1 decomposition.
* `run_rep.py` – one replicate: ASTRAL-Pro 2, ASTRID/ASTRAL-DISCO, ASTRID/ASTRAL-DISCO-R (S0 = ASTRAL-Pro 2 tree rooted by min total DL), DISCO-R-lca, DISCO-R with the true species tree (oracle), one extra iteration; RF vs true species tree; tagging accuracy vs SimPhy locus trees.
* `batch_qr.sh`, `batch_disco.sh` – restartable batches (4 parallel jobs).
* `aggregate.py` – tables, Wilcoxon signed-rank, W/T/L, runtime.

Tools (built in /opt/tools): DISCO v1.4.1 (JSdoubleL/DISCO), ASTRAL-Pro 2 = ASTER `astral-pro` v1.16.2.4 (commit 0cf694f, last before ASTRAL-Pro3), ASTRAL-IV `astral` v1.25.4.8 (ASTER HEAD ddc3dc6), ASTRID = wASTRID/internode v0.0.7 `--preset vanilla` (patched to build with bindgen 0.69).
