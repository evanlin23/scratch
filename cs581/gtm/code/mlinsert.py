"""GTM-Blend-ML: blending disjoint tree merger by constrained ML placement.

Start from the largest subset tree. In rounds: RAxML-NG optimizes branch lengths of
the current tree (fixed GTR+G model), EPA-ng places all remaining taxa, and each taxon
is inserted on its highest-likelihood edge among the edges that keep its own subset
tree induced (T*|S_i == T_i restricted to the inserted part of S_i). Several taxa are
inserted per round; a taxon whose chosen edge disappeared or became infeasible waits
for the next round.
"""
import json
import os
import re
import subprocess

from phylo import read_tree
import insert

BIN = os.environ.get("BIOBIN", "/opt/mm/root/envs/bio/bin")


def jplace_edges(tree_str):
    """Map EPA-ng edge number -> frozenset of leaf names below that edge."""
    toks = re.findall(r"\(|\)|,|;|[^(),;]+", tree_str)
    stack = [[]]
    edges = {}
    last = None
    for t in toks:
        if t == "(":
            stack.append([])
            last = None
        elif t == ",":
            last = None
        elif t == ")":
            ch = stack.pop()
            s = set()
            for c in ch:
                s |= c
            last = s
            stack[-1].append(s)
        elif t == ";":
            break
        else:
            m = re.match(r"([^:{]*)(?::[^{]*)?(?:\{(\d+)\})?", t.strip())
            name, num = m.group(1), m.group(2)
            if last is None:  # leaf
                s = {name}
                stack[-1].append(s)
                last = s
            if num is not None:
                edges[int(num)] = frozenset(last)
    return edges


def write_fasta(seqs, names, path):
    with open(path, "w") as f:
        for n in names:
            f.write(">%s\n%s\n" % (n, seqs[n]))


def run(names, seqs, subsets, guide_tree, workdir, threads=1, log=print, max_rounds=200):
    os.makedirs(workdir, exist_ok=True)
    g = insert.Grow(names, None, None)
    idx = g.idx
    subsets = sorted(subsets, key=lambda t: -len(t.leaves()))
    g.init_from(subsets[0])
    sub_of = {}
    for i, t in enumerate(subsets):
        for v in t.leaves():
            sub_of[t.label[v]] = i
    inserted = {i: [] for i in range(len(subsets))}
    done = {names[v] for v in g.leaves()}
    # model parameters once, on the guide tree
    full = f"{workdir}/full.fa"
    write_fasta(seqs, names, full)
    if not os.path.exists(f"{workdir}/model.raxml.bestModel"):
        subprocess.run([f"{BIN}/raxml-ng", "--evaluate", "--msa", full, "--tree", guide_tree, "--model", "GTR+G",
                        "--prefix", f"{workdir}/model", "--threads", str(threads), "--redo", "--log", "ERROR"],
                       check=True, capture_output=True)
    model = open(f"{workdir}/model.raxml.bestModel").read().split(",")[0].strip()
    stats = dict(rounds=0, fallback=0, placed=0)
    for rnd in range(max_rounds):
        rest = [n for n in names if n not in done]
        if not rest:
            break
        stats["rounds"] += 1
        cur = list(sorted(done))
        write_fasta(seqs, cur, f"{workdir}/ref.fa")
        write_fasta(seqs, rest, f"{workdir}/q.fa")
        g.recompute(with_fitch=False)
        open(f"{workdir}/cur.tre", "w").write(g.to_newick() + "\n")
        subprocess.run([f"{BIN}/raxml-ng", "--evaluate", "--msa", f"{workdir}/ref.fa", "--tree", f"{workdir}/cur.tre",
                        "--model", model, "--opt-model", "off", "--prefix", f"{workdir}/rnd", "--threads", str(threads),
                        "--redo", "--log", "ERROR"], check=True, capture_output=True)
        subprocess.run([f"{BIN}/epa-ng", "--ref-msa", f"{workdir}/ref.fa", "--tree", f"{workdir}/rnd.raxml.bestTree",
                        "--query", f"{workdir}/q.fa", "--model", model, "--filter-max", "100000",
                        "--filter-min-lwr", "0", "--no-heur", "--outdir", workdir, "--redo", "--threads", str(threads)],
                       check=True, capture_output=True)
        jp = json.load(open(f"{workdir}/epa_result.jplace"))
        emap = jplace_edges(jp["tree"])
        fi = jp["fields"]
        iL, iE, iW = fi.index("likelihood"), fi.index("edge_num"), fi.index("like_weight_ratio")
        plist = []
        for pl in jp["placements"]:
            nm = pl["n"][0] if "n" in pl else pl["nm"][0][0]
            ps = sorted(pl["p"], key=lambda r: -r[iL])
            plist.append((-(ps[0][iW]), nm, ps))
        plist.sort()  # most confident first
        # EPA edge number -> current edge (p,c), fixed for this round
        g.recompute(with_fitch=False)
        sp2edge = {}
        for p, c in g.edges():
            b = g.C[c]
            sp2edge[b if not (b & 1 << g.root) else g.ALL ^ b] = (p, c)
        num2edge = {}
        for num, s in emap.items():
            b = 0
            for x in s:
                b |= 1 << idx[x]
            num2edge[num] = sp2edge[b if not (b & 1 << g.root) else g.ALL ^ b]
        touched = set()
        n_ins = 0
        for _, ln, ps in plist:
            i = sub_of[ln]
            Qn = inserted[i]
            Q = 0
            for q in Qn:
                Q |= 1 << idx[q]
            tau = None
            if len(Qn) >= 3:
                tau = insert.norm(insert.sibling_sides(subsets[i], ln, Qn, idx)[0], Q)
            feas = set(insert.feasible_edges(g, Q, tau))
            choice = None
            defer = False
            for r in ps:  # by decreasing likelihood
                e = num2edge[r[iE]]
                if e in touched:
                    defer = True  # best edge was split this round: re-place next round
                    break
                if e in feas:
                    choice = e
                    break
            if defer:
                continue
            if choice is None:
                stats["fallback"] += 1
                sides = insert.sibling_sides(guide_tree_obj(guide_tree), ln, list(done), idx)
                choice = min(feas, key=lambda pc: min((g.C[pc[1]] ^ r).bit_count() for r in sides))
                if choice in touched:
                    continue
            p, c = choice
            g.insert(idx[ln], p, c)
            touched.add(choice)
            done.add(ln)
            inserted[i].append(ln)
            n_ins += 1
            stats["placed"] += 1
            g.recompute(with_fitch=False)
        log(f"round {rnd + 1}: inserted {n_ins}, remaining {len(names) - len(done)}")
    g.recompute(with_fitch=False)
    return g, stats


_gcache = {}


def guide_tree_obj(path):
    if path not in _gcache:
        _gcache[path] = read_tree(path)
    return _gcache[path]
