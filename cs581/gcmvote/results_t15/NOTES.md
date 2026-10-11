AI-assisted (Claude), exploration code for CS581 project

SIMHIGH_R15, one MAGUS draw (gcmtrees run.sh, paper flags), helper t15.
- aln.jsonl: first rows are gg.py rows (linsi = magus, wsoft0.03:linsi&fftns2#es4 = recipe, linsi#es3 = es3,
  linsi&fftns2-op3 = hard, linsi#es4 = es4); then gcmvote run.py rows (B=10) for magus, es4, hard, soft, soft2,
  soft4, hard-bb, soft-bb, hard+mask. Error = avgErr = (SPFN+SPFP)/2.
- vote.py `magus` reproduces MAGUS's merge: same 1000 rows with identical gapped sequences (8356 columns);
  the file differs from MAGUS's only in row order / line wrapping (not byte-identical), same SPFN/SPFP.
- Vote variants were run on /opt/work/gcmvote/reps/SIMHIGH_R15, whose inputs/ symlink MAGUS's subalignments and
  backbones from the gcmtrees rep.
- trees.jsonl: FastTree -lg -gamma via protbench trees.py, one process per tree. The `true` tree was built from
  /opt/data/sim/SIMHIGH/R15/sim.fa, which holds the same names and aligned sequences as the rep's true.fasta
  (formatting differs only). Byte-identical vote alignments reuse an existing tree (listed below if any).
- No vote alignment was byte-identical to an already-treed one; all 13 trees were built.
