#!/bin/bash
# waits for a running pipeline, then runs jobs in order. Job tokens: <ds>:<v> (pipeline) or loo:<set> (stock+fix)
while pgrep -f "run_pipeline2?.sh" >/dev/null; do sleep 15; done
for j in "$@"; do
  a=${j%%:*}; b=${j#*:}
  if [ "$a" = loo ]; then
    for v in stock fix; do /home/user/scratch/cs581/picrust/code/run_loo.sh /opt/work/loo/$b $v; done
  else
    /home/user/scratch/cs581/picrust/code/run_pipeline2.sh $a $b
  fi
done
