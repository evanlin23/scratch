"""Orthology-call accuracy of estimated rooting+tagging on TRUE gene trees.
Usage: python tagacc.py TT.trees MAX_GENES   (TT = truetag.py output; leaves species_locus_indiv)
For every pair of leaves from different species, 'ortholog' iff their LCA is tagged S. Compares the truth with
 - astral-pro3 -T (ASTRAL-Pro3's own rooting/tagging), and
 - DISCO v1.4.1 get_min_root + tag (MinDup rooting + species-overlap tags, as used by DISCO / ASTRID-Pro).
Prints one JSON line: pair counts, false-ortholog rate (true paralog pairs called orthologs), false-paralog rate."""
import json, os, subprocess, sys, tempfile
import numpy as np
import treeswift
sys.path.insert(0, "/opt/src/DISCO")
import disco  # noqa: E402

sp = lambda x: x.split("_")[0]


def classes(t, idx):
    n = len(idx)
    M = np.zeros((n, n), dtype=np.int8)  # 1 = ortholog, 2 = paralog
    for v in t.traverse_postorder():
        if v.is_leaf():
            v.L = [idx[v.label]]
            continue
        ch = v.children
        tagv = 2 if (getattr(v, "tag", None) or v.label) == "D" else 1
        for a in range(len(ch)):
            for b in range(a + 1, len(ch)):
                A, B = np.array(ch[a].L), np.array(ch[b].L)
                M[np.ix_(A, B)] = tagv
                M[np.ix_(B, A)] = tagv
        v.L = [x for c in ch for x in c.L]
    return M


def main():
    tt, mx = sys.argv[1], int(sys.argv[2])
    lines = [l.strip() for l in open(tt) if ";" in l][:mx]
    td = tempfile.mkdtemp(dir=os.environ.get("BENCH_TMP"))
    plain = os.path.join(td, "p.trees")
    with open(plain, "w") as f:
        for l in lines:
            t = treeswift.read_tree_newick(l)
            for v in t.traverse_internal():
                v.label = None
            f.write(t.newick().replace("[&R] ", "") + "\n")
    labs = set()
    for l in lines:
        labs |= {x.label for x in treeswift.read_tree_newick(l).traverse_leaves()}
    mp = os.path.join(td, "map.txt")
    open(mp, "w").write("".join("%s %s\n" % (x, sp(x)) for x in labs))
    ap = subprocess.run(["/opt/mm/root/envs/gdl/bin/astral-pro3", "-T", "-a", mp, plain], capture_output=True,
                        text=True).stdout.strip().split("\n")
    tot = {"astral-pro3": [0, 0, 0, 0], "disco": [0, 0, 0, 0]}  # true-para, para->orth, true-orth, orth->para
    for i, l in enumerate(lines):
        T = treeswift.read_tree_newick(l)
        leaves = [x.label for x in T.traverse_leaves()]
        if len(leaves) > 1500:
            continue
        idx = {x: j for j, x in enumerate(leaves)}
        spp = np.array([sp(x) for x in leaves])
        diff = spp[:, None] != spp[None, :]
        np.fill_diagonal(diff, False)
        Mt = classes(T, idx)
        E1 = treeswift.read_tree_newick(ap[i])
        D = treeswift.read_tree_newick(l)
        for v in D.traverse_internal():
            v.label = None
        r, _, _ = disco.get_min_root(D, sp)
        disco.reroot_on_edge(D, r)
        disco.tag(D, sp)
        for name, E in (("astral-pro3", E1), ("disco", D)):
            Me = classes(E, idx)
            para, orth = diff & (Mt == 2), diff & (Mt == 1)
            c = tot[name]
            c[0] += para.sum() / 2; c[1] += (para & (Me == 1)).sum() / 2
            c[2] += orth.sum() / 2; c[3] += (orth & (Me == 2)).sum() / 2
    out = {"file": tt, "genes": len(lines)}
    for k, c in tot.items():
        out[k] = {"true_para_pairs": int(c[0]), "false_orth_rate": round(c[1] / max(c[0], 1), 4),
                  "true_orth_pairs": int(c[2]), "false_para_rate": round(c[3] / max(c[2], 1), 4)}
    print(json.dumps(out))


if __name__ == "__main__":
    main()
