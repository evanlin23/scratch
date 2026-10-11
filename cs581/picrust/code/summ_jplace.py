#!/usr/bin/env python3
"""Summarise cmp_jplace.py output (A=stock, B=fix). Usage: summ_jplace.py file.tsv [aligned_queries.fasta]"""
import csv, sys, json, statistics as st
r = list(csv.DictReader(open(sys.argv[1]), delimiter='\t'))
ch = [x for x in r if x['edgeA'] != x['edgeB']]
dl = [float(x['loglB']) - float(x['loglA']) for x in r]
o = dict(n=len(r), logl_changed=sum(abs(d) > 1e-3 for d in dl), edge_changed=len(ch),
         edge_dist_median=st.median([int(x['edge_dist']) for x in ch]) if ch else 0,
         edge_dist_max=max([int(x['edge_dist']) for x in ch]) if ch else 0,
         lwr_mean_stock=st.mean(float(x['lwrA']) for x in r), lwr_mean_fix=st.mean(float(x['lwrB']) for x in r),
         lwr_changed_stock=st.mean(float(x['lwrA']) for x in ch) if ch else None,
         lwr_changed_fix=st.mean(float(x['lwrB']) for x in ch) if ch else None,
         lwr99_stock=sum(float(x['lwrA']) >= 0.99 for x in r), lwr99_fix=sum(float(x['lwrB']) >= 0.99 for x in r),
         dlogl_min=min(dl), dlogl_max=max(dl))
if len(sys.argv) > 2:
    seqs, n = {}, None
    for l in open(sys.argv[2]):
        l = l.strip()
        if l.startswith('>'): n = l[1:]; seqs[n] = ''
        else: seqs[n] += l
    st0 = [len(s) - len(s.lstrip('-.')) for s in seqs.values()]
    o.update(aln_cols=len(next(iter(seqs.values()))), first_col_min=min(st0), first_col_median=st.median(st0), n_start_col0=sum(x == 0 for x in st0))
print(json.dumps(o))
