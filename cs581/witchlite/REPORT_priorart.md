## Prior art and novelty check (searched 2026-10-10)

Searches: web (Google-style) for "UPP2 hierarchical search HMM selection", "WITCH faster HMM selection
TIPP3 2025 2026 preprint", "WITCH adjusted bit-score hierarchical top-k", "SEPP/TIPP/UPP early stop
ensemble", "bioRxiv 2025-2026 Warnow fast alignment reads marker genes"; read UPP2 (PMC9846425), TIPP3
(PMC11970662), the J-bandit paper (PMC11211838), the WITCH README and the source of witch-msa 1.0.10
and tipp3 0.5 (PyPI).

| work | what it does | overlap with WITCH-lite |
|---|---|---|
| UPP2 (Park, Ivanovic, Chu, Shen, Warnow, *Bioinformatics* 39(1):btad007, 2023) | Replaces UPP's all-against-all query-vs-HMM search by a descent of the HMM decomposition tree: score both children, move to the one with the higher adjusted bit-score ("evaluates at most two HMMs per level"); *EarlyStop* stops when both children score lower. Picks **one** HMM. 16S.B.ALL: ~10 h vs ~110 h, "did not seem to impact accuracy". | **Idea (a) is UPP2's search**, for top-1. Not applied to WITCH's top-k weighting; no read-length queries. |
| J-bandit (Mazooji & Shomorony, *Bioinformatics* 40 Suppl 1:i328, 2024) | k-mer "J-score" + multi-armed bandit (sequential halving) to find the top-1 HMM in UPP; exact scores for 10 finalists. 16S.B.ALL ~1.9k s vs 9.1k s UPP2. Accuracy "similar, though often slightly degraded". Notes weak behaviour on short queries (no shared k-mers). | Top-1 for UPP; does not mention WITCH; no reads. |
| WITCH-NG (Liu & Warnow, *Bioinform Adv* 2023) | Faster *merging* (exact two-alignment MWT instead of GCM). Still all-vs-all HMM scoring. | Orthogonal (merge step). Now WITCH's default mode (`-m witch-ng`). |
| HMMerge (Park & Warnow, *Bioinform Adv* 2023) | Uses a selected subset of the ensemble to build one merged HMM per query. | Selection is still by scoring all HMMs. |
| TIPP3 (Shen, Wedell, Pop, Warnow, *PLOS CB* 2025) | Calls WITCH on each marker gene's full reference alignment + decomposition tree; WITCH re-decomposes and searches all-vs-all on every run. States WITCH step is "the biggest contribution to runtime" and that methods "substantially faster but not much less accurate than WITCH" are "a promising direction". | The open problem. TIPP3-fast's answer was BLAST. |
| TIPP-SD (ACM-BCB 2025 / *PLOS CB* 2026) | Species detection on TIPP3-fast (BLAST alignment). | Uses BLAST because WITCH was too slow. |

Code check: witch-msa 1.0.10 has no option to restrict which HMMs are scored (options: -k, -w,
--save-weight, -A, -Z; `SearchAlgorithm.search` always runs `hmmsearch --max` of every query chunk
against every HMM). It does read a precomputed `weights.txt`, which is how this pilot plugs in.

Novelty verdict: (a) hierarchical descent is **not novel as an idea** (UPP2), only its use for
WITCH's top-k weighting and for reads would be. (b) BLAST-guided selection (score only the HMMs on
the decomposition path of the query's BLAST top hit; in TIPP3 the BLAST hit already exists from
the read-binning step) and (c) adaptive k by cumulative weight were not found in the literature.
Together these are an incremental-but-open engineering/benchmark contribution, not a new algorithmic
idea.
