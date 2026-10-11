AI-assisted (Claude), exploration code for CS581 project

Rep bank: inputs to redo merge-only runs (vote.py / gg.py) without re-running MAGUS. One tarball per replicate,
listed in MANIFEST.tsv. Each holds the gcmgen-layout rep dir minus merged outputs: inputs/ (subalignments +
backbones), sets/ and aligned/ (gg.py's backbone files), true.fasta (reference), unaligned.fasta, magus.json,
subsets.json (MAGUS's subset order from vote.py, paths made relative to the rep dir), true.tree (if simulated;
ROSE: rose.tt). Extract and pass the extracted dir as REP.
