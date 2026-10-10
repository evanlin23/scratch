# Prior art: adding sequences (UPP family) and long / over-long queries

Gathered by a web literature pass (≈40 min, October 2026). Items marked
"unverified" were not checked against the paper itself.

## Bottom line

I found **no published benchmark of UPP, UPP2, WITCH, WITCH-NG, EMMA or HMMerge on
queries that are longer than the backbone or family**, for example sequences with
non-homologous flanks or an extra domain.

- **Every simulated length-heterogeneity condition shortens sequences.** The HF, LF and
  UHF fragment protocols (Smirnov & Warnow 2021; the WITCH, UPP2 and HMMerge papers)
  cut a prefix and a suffix off full-length sequences, so every query is a substring of
  a homologous sequence.
- **The long side of the backbone rule is never studied.** All of these methods put a
  sequence in the backbone only if its length is within 25% of the median, so a
  sequence longer than 1.25 × the median automatically becomes a query. Biological sets
  (CRW, HomFam) probably contain a few such sequences, but no paper reports results
  for them separately.
- **The papers frame heterogeneity as short sequences.** WITCH-NG describes it as large
  indels and unassembled reads. The EMMA paper's stated limitation of HMM-based adding,
  that query letters with no homolog in the backbone cannot be aligned to each other,
  is exactly the "second-domain" case tested here. Linking the two is our inference.
- **MAFFT `--addlong`** exists for sequences "(much) longer than the sequences already
  aligned". It was listed as alpha on the MAFFT site, and none of the papers above
  benchmark it.
- **Error-prone reads** are covered by other work: Rohmer et al. 2024 (MSA of long
  reads) and TIPP3 2025 (placing ONT/PacBio reads for abundance profiling). Neither is
  a UPP-family MSA benchmark.

## Per-paper notes

**UPP** (Nguyen et al. 2015), doi:10.1186/s13059-015-0688-z
- Data: ROSE NT/AA, RNASim up to 1M sequences, INDELible 10000M2–M4, CRW, 10AA and
  HomFam. Fragmentary versions had 12.5–50% of sequences fragmented.
- Backbone: full-length sequences, within 25% of the median.
- Long queries: not tested.

**UPP2** (Park et al. 2023), doi:10.1093/bioinformatics/btad007
- Data: ROSE 1000S/M/L full-length and HF, RNASim1000(-HF), CRW and HomFam.
- Long queries: not tested.

**WITCH** (Shen, Park & Warnow 2022), doi:10.1089/cmb.2021.0585
- Data: ROSE-HF/LF, among others. Its weighted ensemble of HMMs beats UPP, most
  clearly at high rates of evolution.
- Long queries: not tested.

**WITCH-NG** (Liu & Warnow 2023), doi:10.1093/bioadv/vbad024
- Data: 1000M1–4-HF, RNASim, CRW, 10AA and HomFam. Same accuracy as WITCH, up to
  5× faster.
- Table 2 is used for validation in this pilot.
- Long queries: not tested.

**EMMA** (Shen, Liu, Williams & Warnow 2023), doi:10.1186/s13015-023-00247-x
- Data: ROSE 1000M, INDELible 5000M-het, RNASim, CRW and Rec/Res.
- Results are not broken down by query length.
- Long queries: not tested.

**HMMerge** (Park & Warnow 2023), doi:10.1093/bioadv/vbad052
- Data: HF versions of ROSE, RNASim and INDELible (the INDELible data are
  IDB-0900513), plus an RNASim UHF condition.
- It uses glocal HMM alignment so that short queries can align to any region.
- Long queries: not tested.

## Data

- **IDB-6128941**, doi:10.13012/B2IDB-6128941_V1: ROSE-HF and ROSE-LF.
  - Conditions: 1000L1/L3/L4/M3/S1/S2/S4, 20 replicates each.
  - HF: 50% of sequences become fragments with a mean length of 25% of the original.
  - LF: 25% of sequences become fragments with a mean length of 50% of the original.
- **IDB-0900513**, doi:10.13012/B2IDB-0900513_V1: INDELible HF, from the HMMerge paper.
  - Indel rates 0.001 and 0.005, 10 replicates each.
  - 501 full-length sequences of about 1025 bp and 500 fragments of about 260 bp.
- **Neither contains over-long sequences**, so a long-query condition has to be
  constructed (`code/make_long.py`).

## Other

- Smirnov & Warnow 2021, Syst Biol, doi:10.1093/sysbio/syaa058 (unverified).
- Shen, Zaharias & Warnow 2022, MAGUS+eHMMs, doi:10.1093/bioinformatics/btab788
  (unverified).
- Shen 2025, PhD thesis, "Multiple sequence alignment with sequence length
  heterogeneity and its applications", UIUC, hdl:2142/129918. Chapters not seen; it
  is worth checking for an over-long chapter before committing to a project.
- TIPP3: Shen et al. 2025, PLoS Comput Biol, doi:10.1371/journal.pcbi.1012593.
- Rohmer, Touzet & Limasset 2024, PeerJ, doi:10.7717/peerj.17731.
- TWILIGHT: Tseng et al. 2025, doi:10.1093/bioinformatics/btaf212.
- Katoh & Frith 2012, MAFFT `--add`, doi:10.1093/bioinformatics/bts578 (unverified).
