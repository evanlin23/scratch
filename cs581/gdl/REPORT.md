# Pilot: Is ASTRID-multi consistent under GDL? Can an "ASTRID-Pro" distance fix it?

*CS581 project pilot, 2026-10-09. Branch `claude/cs581-gdl`. Code in `code/`, result files in `results/`.*

## 1. Question

Source: lecture deck *Phylogenomics part 2*, slide 35, "(Some) Open Questions" (see `../literature/slide_open_problems.md`, item 7):

- "Can we find a distance correction for GDL so that ASTRID and NJst are statistically consistent under GDL models?"
- "Unknown if distance-based species tree estimation (e.g., ASTRID-multi) is statistically consistent under GDL models."
- Related questions on the same slide: which other methods are consistent under GDL/DLCoal, and whether ASTRAL-Pro is consistent under a random model of rooting and tagging error.

The pilot tests one candidate answer, **ASTRID-Pro**:
1. Root and tag each gene-family tree.
2. For each species pair, keep only orthologous copy pairs (LCA is a speciation node).
3. Count only speciation nodes on the path.
4. Average within each gene, then across genes.
5. Build the tree with FastME.

__RESULTS__
