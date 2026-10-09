"""Run Spectral Cluster Supertree (McArthur et al. 2024) on rooted source trees: scs_run.py in out"""
import sys
from sc_supertree import load_trees, construct_supertree
t = construct_supertree(load_trees(sys.argv[1]), pcg_weighting="branch")
open(sys.argv[2], "w").write(t.get_newick(with_distances=False) + "\n")
