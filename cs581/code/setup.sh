#!/bin/bash
# Reproduce the cs581 MSA environment on Ubuntu 24.04 (what the pilots used).
# Usage: bash cs581/code/setup.sh   (run from the repository root)
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"

SUDO=""; [ "$(id -u)" = 0 ] || SUDO=sudo

# 1. Base tools from apt (MAFFT 7.505, MCL, FastTree, Clustal Omega, HMMER, MUSCLE 5.1)
$SUDO apt-get update -q
$SUDO env DEBIAN_FRONTEND=noninteractive apt-get install -y -q mafft mcl fasttree clustalo hmmer muscle openjdk-17-jre-headless unzip

# 2. Newer aligners from bioconda via micromamba (IQ-TREE 3 / AliSim, FAMSA, MUSCLE 5.3, TWILIGHT, PASTA)
if [ ! -x /opt/mm/micromamba ]; then
  $SUDO mkdir -p /opt/mm && $SUDO chown "$(id -u)" /opt/mm && curl -sSL https://conda.anaconda.org/conda-forge/linux-64/micromamba-2.9.0-0.tar.bz2 \
    | tar xj -C /opt/mm bin/micromamba && mv /opt/mm/bin/micromamba /opt/mm/micromamba
fi
if [ ! -x /opt/mm/root/envs/bio/bin/twilight ]; then
  MAMBA_ROOT_PREFIX=/opt/mm/root /opt/mm/micromamba create -y -q -n bio -c conda-forge -c bioconda \
    iqtree=3.1.4 famsa muscle=5.3 twilight pasta
fi

# 3. FastSP (alignment accuracy: SPFN/SPFP/TC)
if [ ! -d /opt/tools/FastSP ]; then
  $SUDO mkdir -p /opt/tools && $SUDO chown "$(id -u)" /opt/tools
  git clone -q https://github.com/smirarab/FastSP.git /opt/tools/FastSP
fi

# 4. MAGUS, pinned to the commit the pilots used (bundles its own MAFFT/MCL binaries)
if [ ! -d "$HERE/MAGUS" ]; then
  git clone -q https://github.com/vlasmirnov/MAGUS.git "$HERE/MAGUS"
  git -C "$HERE/MAGUS" checkout -q 39041fc8da5dcb44c95e90c212c667cf225ec129
fi
pip install -q -e "$HERE/MAGUS" numpy scipy pandas matplotlib

# 5. Benchmark data from the MAGUS paper (Illinois Data Bank IDB-2643961, Datasets.zip, ~390 MB):
#    ROSE 1000-taxon (L1-3, M1-4, S1-3; 20 reps), RNASim 1k/10k, 16S (CRW), BAliBASE
if [ ! -d /opt/data/Datasets/ROSE ]; then
  $SUDO mkdir -p /opt/data && $SUDO chown "$(id -u)" /opt/data
  curl -sSL -A "Mozilla/5.0 (X11; Linux x86_64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/120.0 Safari/537.36" -o /opt/data/Datasets.zip "https://databank.illinois.edu/datafiles/u373n/download"
  (cd /opt/data && unzip -q -o Datasets.zip && rm Datasets.zip)
fi
# Length-filtered BAliBASE references used by the paper are committed in cs581/data/balibase_clean/
echo "setup done"
