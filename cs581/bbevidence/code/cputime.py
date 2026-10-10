"""CPU seconds (user+sys of the child) to align ONE backbone (MAGUS's backbone_1 set) with each tool,
one thread each.   python3 cputime.py REP_DIR [...]  -> TSV lines: rep tool cpu_s wall_s"""
import os
import resource
import subprocess
import sys
import time

M = "/home/user/scratch/cs581/code/MAGUS/magus/tools/mafft/mafft"
TOOLS = {
    "linsi": [M, "--localpair", "--maxiterate", "1000", "--ep", "0.123", "--quiet", "--thread", "1", "--anysymbol"],
    "clustalo": ["clustalo", "--threads", "1", "--force", "--outfmt", "fa", "-i"],
    "fftns2": [M, "--retree", "2", "--maxiterate", "0", "--quiet", "--thread", "1", "--anysymbol"],
    "fftns2-op3": [M, "--retree", "2", "--maxiterate", "0", "--op", "3.0", "--quiet", "--thread", "1", "--anysymbol"],
    "famsa": ["/opt/mm/root/envs/bio/bin/famsa", "-t", "1"],
    "ginsi": [M, "--globalpair", "--maxiterate", "1000", "--quiet", "--thread", "1", "--anysymbol"],
}
for rep in sys.argv[1:]:
    f = os.path.join(rep, "sets", "s0", "backbone_1.fa")
    for name, argv in TOOLS.items():
        before = resource.getrusage(resource.RUSAGE_CHILDREN)
        start = time.time()
        args = argv + [f] + (["/dev/null"] if name == "famsa" else [])
        subprocess.run(args, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True)
        after = resource.getrusage(resource.RUSAGE_CHILDREN)
        cpu = after.ru_utime - before.ru_utime + after.ru_stime - before.ru_stime
        print("{}\t{}\t{:.1f}\t{:.1f}".format(os.path.basename(rep.rstrip("/")), name, cpu, time.time() - start), flush=True)
