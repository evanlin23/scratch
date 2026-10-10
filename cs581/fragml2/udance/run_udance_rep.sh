#!/usr/bin/env bash
# usage: run_udance_rep.sh R0   (single-thread, nice'd uDance run on /opt/data/mlcache/M1HF/<rep>)
set -euo pipefail
REP=$1
D=/home/user/scratch/cs581/fragml2/udance
WD=/opt/udance_work/$REP
export PATH=/opt/mm/root/envs/udance/bin:$PATH
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
[ -f $WD/backbone.nwk ] || python3 $D/prep_udance.py /opt/data/mlcache/M1HF/$REP $WD
sed "s|WORKDIR|$WD|" $D/config_template.yaml > $WD/config.yaml
cd /opt/udance_work/uDance
/opt/mm/root/envs/udance/bin/time -v -o $WD/time_v.txt nice -n 10 snakemake --cores 1 --configfile $WD/config.yaml --snakefile udance.smk all > $WD/snakemake.log 2>&1
