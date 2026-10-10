"""Keep exactly 4 tree-search processes busy (counting jobs already running from earlier runners);
launch queued jobs in priority order; skip jobs with a .done file. Restartable.
Usage: python scheduler.py"""
import os, glob, subprocess, time, shlex
IQ = '/opt/work/bin/iqtree3'; RX = '/opt/mm/root/envs/bio/bin/raxml-ng'; W = '/opt/work/runs'; P = 4
alns = sorted(glob.glob('/opt/work/sim/*.phy')) + sorted(glob.glob('/opt/work/emp/*.fa'))
def b(f): return os.path.basename(f).rsplit('.', 1)[0]
def models(f):
    return ('LG+G4', 'LG+G4') if b(f).startswith('empAA') else ('GTR+F+I+G4', 'GTR+FO+IO+G4')
def wall(ds, tag):
    p = f'{W}/{ds}/{tag}.done'
    return float(open(p).read().split()[3]) if os.path.exists(p) else None
Q = []
for f in alns:
    M, _ = models(f); Q.append((b(f), 'iqfast', f'{IQ} -s {f} -m {M} -T 1 -seed 1 -pre {W}/{b(f)}/iqfast -redo --quiet --fast'))
for f in alns:
    M, MR = models(f)
    Q.append((b(f), 'iqdef1', f'{IQ} -s {f} -m {M} -T 1 -seed 1 -pre {W}/{b(f)}/iqdef1 -redo --quiet'))
    Q.append((b(f), 'rxfast', f'{RX} --fast --msa {f} --model {MR} --threads 1 --seed 1 --prefix {W}/{b(f)}/rxfast --redo --log PROGRESS'))
for f in alns:
    M, _ = models(f)
    if b(f) in ('rg42163_n94', 'rg51878_n319', 'rg6754_n258', 'rg21213_n117', 'emp16S_16S3_1_n150', 'empAA_RV100_BBA0117_n200'):
        Q.append((b(f), 'iqnstop20', f'{IQ} -s {f} -m {M} -T 1 -seed 1 -nstop 20 -pre {W}/{b(f)}/iqnstop20 -redo --quiet'))
for f in alns:
    M, MR = models(f)
    if b(f).startswith('empAA'): continue
    Q.append((b(f), 'rxclassic1', f'{RX} --search --tree pars{{1}} --opt-topology classic --msa {f} --model {MR} --threads 1 --seed 1 --prefix {W}/{b(f)}/rxclassic1 --redo --log PROGRESS'))
for f in alns:
    M, _ = models(f)
    Q.append((b(f), 'iqdef2', f'{IQ} -s {f} -m {M} -T 1 -seed 2 -pre {W}/{b(f)}/iqdef2 -redo --quiet'))

def busy():
    out = subprocess.run(['pgrep', '-f', r'^(/opt/work/bin/iqtree3 -s|/opt/mm/root/envs/bio/bin/raxml-ng --)'], capture_output=True, text=True).stdout
    return len(out.split())
def running_tags():
    out = subprocess.run(['pgrep', '-af', 'iqtree3 -s|raxml-ng --'], capture_output=True, text=True).stdout
    return out
WRAP = 's=$(date +%s.%N); {cmd} > {d}/{tag}.stdout 2>&1; rc=$?; e=$(date +%s.%N); echo "{ds} {tag} $rc $(echo "$e - $s" | bc)" > {d}/{tag}.done'
for ds, tag, cmd in Q:
    d = f'{W}/{ds}'; os.makedirs(d, exist_ok=True)
    if os.path.exists(f'{d}/{tag}.done'): continue
    if tag == 'iqdef2':
        w1 = wall(ds, 'iqdef1')
        if w1 is None or w1 > 700: continue
    pre = cmd.split('-pre ')[1].split()[0] if '-pre ' in cmd else cmd.split('--prefix ')[1].split()[0]
    if pre in running_tags(): continue   # already running from an older runner
    while busy() >= P: time.sleep(5)
    subprocess.Popen(['bash', '-c', WRAP.format(cmd=cmd, d=d, tag=tag, ds=ds)], start_new_session=True)
    print(time.strftime('%H:%M:%S'), 'launched', ds, tag, flush=True)
    time.sleep(2)
print('queue exhausted', flush=True)
