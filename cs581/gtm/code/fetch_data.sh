#!/bin/bash
# Download the GTM-pipeline data (Park, Zaharias & Warnow 2021, doi:10.13012/B2IDB-7008049_V1),
# the SATé 1000M1 true trees (doi:10.13012/B2IDB-5139418_V1) and the GTM source into $GTMDATA.
set -euo pipefail
D=${GTMDATA:-/opt/gtmdata}; mkdir -p $D; cd $D
UA="Mozilla/5.0 (X11; Linux x86_64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/120.0 Safari/537.36"
get() { [ -s "$2" ] || curl -sSL -A "$UA" -o "$2" "https://databank.illinois.edu/datafiles/$1/download"; }
get 1xjr4 README_ds.txt
get eu0zs RNASim1000.tar.gz; get vyuj4 RNASim1000_Analysis.tar.gz
get 6446x Cox1-HET.zip;      get 1k83z Cox1-Het_Analysis.tar.gz
get p4e01 1000M1_HF_Analysis.tar.gz
mkdir -p d cox m1 sate1000M1 m1_true
tar xzf RNASim1000.tar.gz -C d; tar xzf RNASim1000_Analysis.tar.gz -C d
unzip -q -o Cox1-HET.zip -d cox; tar xzf Cox1-Het_Analysis.tar.gz -C cox
tar xzf 1000M1_HF_Analysis.tar.gz -C m1
if [ ! -d "sate1000M1/1000M1" ]; then
  get x5h4q sate.zip; unzip -q -o sate.zip sate_journal/1000M1.tar.bz2; tar xjf sate_journal/1000M1.tar.bz2 -C sate1000M1; rm sate.zip
fi
S="sate1000M1/1000,1000,.0000082,.005,medium_gap_pdf,GTR+second,35,2.0,1.0"
for i in 0 1 2 3 4; do ln -sfn "$D/$S/R$i" m1_true/R$i; done
[ -d GTM_src ] || git clone -q https://github.com/vlasmirnov/GTM.git GTM_src
git -C GTM_src checkout -q 18e3bc9
echo "data ready in $D"
