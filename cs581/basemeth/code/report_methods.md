## Methods

### Design: change only the subset aligner

For each dataset, one MAGUS run with the paper's flags (25 subsets, PASTA-style decomposition with a
300-sequence skeleton, 10 backbones of 200 sequences aligned by MAFFT L-INS-i, GCM graph from the backbones, MCL,
minclusters trace, `-f 4`) supplies the 25 subsets and the 10 backbone alignments. BAliBASE and nucleotide
runs reuse earlier MAGUS runs with the same flags (`cs581/experiments/runs/*/inputs.tar.xz`); HomFam and AliSim
runs are new (`code/basemeth.py prep`). Each subset is then re-aligned from its unaligned sequences by each base
method (one thread per subset, 4 subsets at a time), and the GCM merge alone
(`gcmx.run_magus -s SUBSETS -b BACKBONES`, same backbones, same MCL / minclusters flags, 4 threads) produces the
final alignment. So the only difference between variants is the subset aligner. `linsi` re-runs MAGUS's own
subset command (`mafft --localpair --maxiterate 1000 --ep 0.123`) and is the control; `orig` merges MAGUS's
original subset alignments untouched (checks that re-running L-INS-i reproduces MAGUS).

Base methods (single thread, default options unless stated):

| key | tool | version | command |
|---|---|---|---|
| linsi | MAFFT L-INS-i (MAGUS-bundled binary) | 7.450 | `--localpair --maxiterate 1000 --ep 0.123 --anysymbol` |
| ginsi | MAFFT G-INS-i | 7.450 | `--globalpair --maxiterate 1000 --anysymbol` |
| muscle5 | MUSCLE | 5.3 | `-align` (PPP, no ensemble) |
| famsa | FAMSA | 2.4.1 | default |
| probcons | ProbCons | 1.12 | default (2 consistency passes, 100 refinement passes) |
| clustalo | Clustal Omega | 1.2.4 | default |
| kalign | Kalign | 3.6.0 | default |
| prank | PRANK | v.250331 (bioconda 251117) | default (`-f=fasta`, no `+F`), nucleotides only |

Scoring: FastSP (SPFN, SPFP, TC); SP error = (SPFN+SPFP)/2 in percent. HomFam alignments are scored on the
Homstrad seed sequences only (the estimate restricted to them, all-gap columns removed). Paired difference =
variant minus `linsi`, in SP-error points; W/T/L with a tie band of |d| < 0.05 points; two-sided Wilcoxon
signed-rank (scipy); Holm over the four pre-registered protein candidates. Subset-level accuracy = each
method's 25 subset alignments scored against the reference restricted to each subset, pooled by homology
counts.

Runtime: per-subset wall and CPU seconds (children rusage) of the aligner; merge wall/CPU. All runs shared a
4-core / 15 GB VM with other jobs of this pilot (MAGUS runs, PASTA), so **wall times are inflated by contention;
CPU seconds are the fair cost comparison**.

Whole-dataset baselines (`code/baselines.py`, 4 threads, wall cap): FAMSA 2.4.1 default; MUSCLE 5.3 `-align`
(`-super5` above 1,000 sequences); regressive T-Coffee 12.00.7fb08c2 `-reg -reg_nseq 100 -reg_tree nj
-reg_method clustalo_msa` (the mBed guide tree hung in this bioconda build and `famsa_msa` is not available in it,
so this is not the published configuration). PASTA 1.8.3 (commit 738bec5, the MAGUS paper's version) with
`--aligner mafft|probcons|prank`, its bundled tools (MAFFT 7.149b, ProbCons, PRANK v.100311 with `+F`), OPAL
merger, FastTree, `--iter-limit 1` for the aligner comparison (3 = default for the reference run).
