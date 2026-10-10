# gdlcons code

Everything runs in the micromamba env `gdl` (`/opt/mm/root/envs/gdl`). Bulk data and logs live under `/opt/runs/gdlcons`; only small result files are copied to `../results`.

## Tools

| tool | how it was installed |
|---|---|
| ASTRAL-Pro3, astral4 (ASTER v1.25.3.8), FastME | bioconda `aster`, `fastme` |
| **astral-pro3-fixed** | ASTER git `ddc3dc6` plus `astral-pro-fixed.patch` (29 lines), built with `g++ -std=gnu++11 -O2 -pthread src/astral-pro.cpp -o bin/astral-pro3-fixed` |
| SimPhy 1.0.2 | GitHub release binary `simphy_lnx64` |
| DISCO (`disco.py`, v1.4.1) | github.com/JSdoubleL/DISCO; needs `treeswift` |
| FastMulRFS | github.com/ekmolloy/fastmulrfs with phylokit, phylonaut and FastRFS built from source |
| STAG | github.com/davidemms/STAG; Python 2.7 env `py27` |
| DupTree | not installed: its GitHub repository needs authentication and there is no bioconda package |

To build FastMulRFS:
- replace `SIGSTKSZ` with a constant in phylokit's `catch.hpp`;
- copy `Astral/` (ASTRAL 5.7.8 jar) next to the FastRFS binary;
- use `fm.tree.single` as the output tree.

### The `astral-pro3-fixed` patch

When `APRO_FIXED=1` is set, `annotateTree` keeps the input root instead of choosing the min-duplication root, and overrides the tags. Tags are encoded in branch lengths:

- non-root-child internal node: length = 1 + tag;
- root children a, b: their summed length minus 2 = tag(a) + 2·tag(b) + 8·tag(root).

`core.encode` writes this encoding. Without `APRO_FIXED` the binary behaves exactly like ASTRAL-Pro3. Branch lengths are therefore meaningless in fixed mode; only the topology is used.

Caveat: the `t1/t2/t3` labels in `freqQuad.csv` are **not** stable across runs (they depend on internal child order). Always parse the split label; `core.split_of` does this.

## Scripts

| script | purpose |
|---|---|
| `gdlsim.py`, `phylo.py`, `methods.py` | copied from the sibling branch `claude/cs581-gdl`. They provide the GDL simulator with true tags, tree utilities, ASTRID-multi/FastME and the root/tag DP. `gdlsim.simulate` gained a `max_tries` argument |
| `core.py` | error models (`true`, `ovl`, `rovl`, `rtrue`, `flip`, `d2s`, `s2d`, `d2sv`, `flipv`, `rsp-X`), the encoder, and ASTRAL-Pro runners (species tree and `-C -u 3` quartet scores) |
| `gen_data.py` | nested 20,000-family datasets. GDL settings use our simulator; DLCOAL settings use SimPhy |
| `run_curve.py` | error vs number of families for the atlas methods and the ASTRAL-Pro error models (restartable jsonl) |
| `analyze_curves.py` | `results/curves_summary.md` and plots |
| `search4.py`, `qgrid.py`, `rootbias.py` | random and grid searches on 4–5 taxa for negative ASTRAL-Pro quartet margins |
| `exact_assoc.py` | exact generating-function survival-set probabilities; association inequalities I1, I2 |
| `predict_margin.py` | **exact limiting ASTRAL-Pro scores** (4-taxon caterpillar; true / overlap / `d2sv(q)` tags) |
| `prevalence.py` | the formula applied to 20,000 random rate configurations |
| `probe4.py`, `block_se.py`, `curve4.py`, `atlas4.py` | one adversarial 4-taxon configuration: big-sample scores, block standard errors, species-tree error vs K, the method atlas |
| `tagacc.py` | ASTRAL-Pro3's own rooting/tagging vs the truth on true gene trees |
| `loc_fail.py` | which branch an error model gets wrong |
