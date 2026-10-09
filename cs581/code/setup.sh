#!/bin/bash
# Reproduce the cs581 MSA environment on Ubuntu 24.04 (what the pilots used).
# Usage: bash cs581/code/setup.sh   (run from the repository root)
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"

# 1. Base tools from apt (MAFFT 7.505, MCL, FastTree, Clustal Omega, HMMER, MUSCLE 5.1)
sudo apt-get update -q
sudo DEBIAN_FRONTEND=noninteractive apt-get install -y -q mafft mcl fasttree clustalo hmmer muscle openjdk-17-jre-headless

# 2. Newer aligners from bioconda via micromamba (IQ-TREE 3 / AliSim, FAMSA, MUSCLE 5.3, TWILIGHT, PASTA)
if [ ! -x /opt/mm/micromamba ]; then
  mkdir -p /opt/mm && curl -sSL https://conda.anaconda.org/conda-forge/linux-64/micromamba-2.9.0-0.tar.bz2 \
    | tar xj -C /opt/mm bin/micromamba && mv /opt/mm/bin/micromamba /opt/mm/micromamba
fi
MAMBA_ROOT_PREFIX=/opt/mm/root /opt/mm/micromamba create -y -q -n bio -c conda-forge -c bioconda \
  iqtree=3.1.4 famsa muscle=5.3 twilight pasta

# 3. FastSP (alignment accuracy: SPFN/SPFP/TC)
[ -d /opt/tools/FastSP ] || git clone -q https://github.com/smirarab/FastSP.git /opt/tools/FastSP

# 4. MAGUS, pinned to the commit the pilots used (bundles its own MAFFT/MCL binaries)
if [ ! -d "$HERE/MAGUS" ]; then
  git clone -q https://github.com/vlasmirnov/MAGUS.git "$HERE/MAGUS"
  git -C "$HERE/MAGUS" checkout -q 39041fc8da5dcb44c95e90c212c667cf225ec129
fi
pip install -q -e "$HERE/MAGUS" numpy scipy pandas matplotlib

# 5. Benchmark data from the MAGUS paper (Illinois Data Bank IDB-2643961, Datasets.zip, ~390 MB):
#    ROSE 1000-taxon (L1-3, M1-4, S1-3; 20 reps), RNASim 1k/10k, 16S (CRW), BAliBASE
if [ ! -d /opt/data/Datasets ]; then
  mkdir -p /opt/data
  curl -sSL -A "Mozilla/5.0" -o /opt/data/Datasets.zip "https://databank.illinois.edu/datafiles/u373n/download"
  (cd /opt/data && unzip -q Datasets.zip)
fi
echo "setup done"
