# AI-assisted (Claude), exploration code for CS581 project
#!/bin/bash
# slot-gated tree: METHOD SRC  (starts only when fewer than 4 FastTree processes are running)
m=$1; src=$2; T=/opt/work/gcmvote/trees/SIMHIGH_R14_d0
ln -sf "$src" $T/$m.fasta
exec 9>/opt/work/t14/slot.lock
flock 9
while true; do until [ $(pgrep -cx FastTree) -lt 4 ]; do sleep 15; done; sleep 30; [ $(pgrep -cx FastTree) -lt 4 ] && break; done
n=$(pgrep -cx FastTree)
nice -n 10 /opt/mm/root/envs/pasta183/bin/python /home/user/scratch/cs581/protbench/code/trees.py $T \
  /opt/data/sim/SIMHIGH/R14/tree.nwk $T/trees.jsonl $m 9>&- &
until [ $(pgrep -cx FastTree) -gt $n ] || ! kill -0 $! 2>/dev/null; do sleep 2; done
flock -u 9; exec 9>&-
wait
