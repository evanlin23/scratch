"""Write an alignment with its sequences in a random order (FastTree search-noise control).
    python3 shuffle.py IN OUT SEED"""
import random, sys
sys.path.insert(0, "/home/user/scratch/cs581/code")
from gcmx import fasta
a = fasta.read(sys.argv[1]); t = list(a); random.Random(int(sys.argv[3])).shuffle(t)
fasta.write({x: a[x] for x in t}, sys.argv[2])
