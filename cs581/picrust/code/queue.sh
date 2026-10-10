#!/bin/bash
for ds in "$@"; do for v in stock fix; do /home/user/scratch/cs581/picrust/code/run_pipeline.sh $ds $v; done; done
