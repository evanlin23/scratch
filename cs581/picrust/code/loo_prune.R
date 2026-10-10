# Rscript loo_prune.R full.tre heldout.txt out.tre : drop held-out tips (castor merges the resulting
# degree-2 nodes, summing branch lengths)
suppressMessages(library(castor))
a <- commandArgs(TRUE)
tr <- read_tree(file = a[1])
drop <- readLines(a[2])
keep <- setdiff(tr$tip.label, drop)
sub <- get_subtree_with_tips(tr, only_tips = keep, collapse_monofurcations = TRUE)$subtree
write_tree(sub, file = a[3])
cat(length(tr$tip.label), "->", length(sub$tip.label), "tips\n")
