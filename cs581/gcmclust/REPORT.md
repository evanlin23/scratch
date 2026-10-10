# Replacing MCL in MAGUS's Graph Clustering Merger (GCM)

Exploration pilot, CS581 (UIUC), Fall 2026. Instructor's boldfaced project: "Modify the MAGUS design by using a
different clustering method instead of the Markov Clustering algorithm within GCM." Merge-only: each replicate's
25 L-INS-i subset alignments and 10 L-INS-i backbones (one MAGUS draw) are fixed; only GCM's clustering step
changes, and MAGUS's own violation purge + `minclusters` trace + writer run unchanged on the new clusters.
Code: `code/` (`gc.py`, `build_graph.py`, `summarize.py`); raw rows: `results/*.results.jsonl`; protocol fixed
before results: `PREREG.md`; all tables: `results/tables_train.md`, `results/tables_held.md`.

__TLDR__

## 1. Prior art

- **MAGUS** (Smirnov & Warnow, *Bioinformatics* 37:1666, 2021). GCM = alignment graph from backbones → MCL
  (inflation 4; the paper says the inflation sweep "follows a standard recommendation" and puts it in the
  supplement: 1000M1 only, values from van Dongen's 1.4/2/4/6) → trace by A* `minclusters`. MCL is motivated
  as extending the consistency principle to longer paths.
- **MWT-AM** (Zaharias, Smirnov & Warnow, *IEEE/ACM TCBB* 20:1700, 2023; AlCoB 2021). Formalises merging as Maximum
  Weight Trace (NP-hard). Tried other *clustering* steps — MLR-MCL and Region Growing (a constraint-aware greedy
  agglomeration close to our `agglo`) — and other *trace* steps (FM, MWT-greedy/search, RG-fast, plus an optimizer).
  Its main win is in the trace (GCM(fm+opt) beat default GCM on every HomFam set); clustering alternatives were
  not better than MCL there. That work used MAGUS's raw graph and HomFam/ROSE, not a support-filtered graph, and
  did not test modularity/CPM community detection.
- **Recursive MAGUS** (Smirnov, *PLOS Comput Biol* 17:e1008950, 2021): recursion, guide trees, scale; GCM's
  clustering is unchanged (MCL); recursion buys scalability, not accuracy.
- **Community detection.** Louvain (Blondel et al. 2008), Leiden (Traag, Waltman & van Eck 2019; fixes
  Louvain's badly-connected communities; supports CPM, which has no resolution limit), label propagation
  (Raghavan et al. 2007). Our literature search (repo notes `cs581/literature/sota_review.md` §P5 and a web search)
  found **no alignment-graph merger that uses Leiden/Louvain/CPM**; residue-clustering MSA (MSARC, Modzelewski &
  Dojer 2014) uses hierarchical clustering of residues, not of subset-alignment columns.
- **In this repo.** `cs581/gcmgen` (branch claude/cs581-gcmgen): cross-subset GCM edges supported by < 4 of 10
  backbones are 2–25% correct; deleting them (`linsi#es4`) gives proteins −2.27 (7/1/2), DNA/RNA −0.02.

__BODY__
