#!/bin/bash
# Toolchain + data for the fragscale pilot (Ubuntu 24.04, 4 cores). Idempotent.
#   bash cs581/fragscale/code/install.sh
set -uo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
MM=/opt/mm/micromamba; export MAMBA_ROOT_PREFIX=/opt/mm/root
UA="Mozilla/5.0 (X11; Linux x86_64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/120.0 Safari/537.36"

apt-get update -q && DEBIAN_FRONTEND=noninteractive apt-get install -y -q fasttree mafft hmmer time unzip \
  cmake build-essential autoconf automake libtool flex bison zlib1g-dev git

if [ ! -x $MM ]; then
  mkdir -p /opt/mm && curl -sSL https://conda.anaconda.org/conda-forge/linux-64/micromamba-2.9.0-0.tar.bz2 \
    | tar xj -C /opt/mm bin/micromamba && mv /opt/mm/bin/micromamba $MM
fi
# IQ-TREE 3 (same pin as cs581/code/setup.sh)
[ -x $MAMBA_ROOT_PREFIX/envs/bio/bin/iqtree3 ] || $MM create -y -q -n bio -c conda-forge -c bioconda iqtree=3.1.4
# RAxML-NG, gappa, stock EPA-ng 0.3.8 (as in the fragml2 pilot)
[ -x $MAMBA_ROOT_PREFIX/envs/fml/bin/raxml-ng ] || $MM create -y -q -n fml -c conda-forge -c bioconda \
  raxml-ng=2.0.3 gappa epa-ng=0.3.8 python=3.11 dendropy scipy
# UPP (SEPP) and WITCH, BSCAMPP (python packages, with their bioconda binaries)
[ -x $MAMBA_ROOT_PREFIX/envs/aln/bin/witch.py ] || $MM create -y -q -n aln -c conda-forge -c bioconda \
  python=3.11 sepp hmmer mafft pip && $MAMBA_ROOT_PREFIX/envs/aln/bin/pip install -q witch-msa bscampp

# patched EPA-ng (cs581/epang/code/epa-ng-fix.patch on v0.3.8)
if [ ! -x /opt/src/epa-ng-fix/bin/epa-ng ]; then
  mkdir -p /opt/src && cd /opt/src
  [ -d epa-ng-fix ] || git clone -q --recursive https://github.com/pierrebarbera/epa-ng.git epa-ng-fix
  cd epa-ng-fix && git checkout -q v0.3.8 && git submodule update -q --init --recursive
  git apply "$HERE/epa-ng-fix.patch" && make -j4 > /opt/src/epa-ng-fix.build.log 2>&1
fi

# data: MAGUS Datasets.zip (ROSE 1000M1, RNASim 1K/10K) and Park et al. 2021 1000M1-HF (IDB-7008049)
mkdir -p /opt/data && cd /opt/data
if [ ! -d /opt/data/Datasets/ROSE ]; then
  curl -sSL -A "$UA" -o Datasets.zip "https://databank.illinois.edu/datafiles/u373n/download" && unzip -q -o Datasets.zip && rm Datasets.zip
fi
if [ ! -d /opt/data/park2021/1000M1_HF ]; then
  mkdir -p park2021 && curl -sSL -A "$UA" -o park2021/a.tgz "https://databank.illinois.edu/datafiles/p4e01/download" \
    && tar xzf park2021/a.tgz -C park2021 && rm park2021/a.tgz
fi
echo install done
