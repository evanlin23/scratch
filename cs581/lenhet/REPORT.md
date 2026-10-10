# Pilot: adding very long sequences into an alignment (UPP / WITCH / EMMA)

Branch `claude/cs581-lenhet`. Code is in `code/`; small result files are in `results/`;
the bulky runs live under `/opt/runs/lenhet` and are not committed.
Pilot wall-clock was about 4.5 h on a 4-core machine.

## 1. Question

The question comes from the CS581 project list (`cs581/notes/CS581-project-suggestions.txt`,
lines 113–119), under "Projects related to adding sequences into an alignment (e.g. UPP,
WITCH, EMMA, etc.)":

> Most of these methods have been designed for handling sequences that are short but not
> necessarily sequences that are long nor sequences that are reads (and so have errors).
> Explore one of these differences.
> Evaluate UPP and other existing MSA methods on datasets with sequence length
> heterogeneity where some sequences are very long (and possibly design a new method that
> addresses this problem better).

This pilot takes up "very long": query sequences that are longer than the backbone family.

## 2. Prior art

The full notes are in `results/prior_art.md`.

**Methods.**
- UPP: doi:10.1186/s13059-015-0688-z
- UPP2: doi:10.1093/bioinformatics/btad007
- WITCH: doi:10.1089/cmb.2021.0585
- WITCH-NG: doi:10.1093/bioadv/vbad024
- EMMA: doi:10.1186/s13015-023-00247-x
- HMMerge: doi:10.1093/bioadv/vbad052

**Has "very long" been benchmarked? Not that we could find.**
- Every length-heterogeneity benchmark in these papers (HF, LF, UHF, "-het") shortens
  sequences: the queries are fragments of full-length homologous sequences.
- All these methods put a sequence in the backbone only if it is within 25% of the
  median length. A sequence more than 1.25× the median is therefore always a query, but
  no paper evaluates that case.
- MAFFT has an `--addlong` mode for exactly this case, but it has not been benchmarked
  in these papers.
- The EMMA paper states one limitation of HMM-based adding: query letters that have no
  homolog in the backbone cannot be aligned to each other. That is precisely the
  "second domain" failure measured below. Linking the two is our observation.
- One gap in our search: Chengze Shen's 2025 UIUC thesis (hdl:2142/129918) is on length
  heterogeneity, and we could not read its chapters. Check it before committing to a
  project.

## 3. Data

### Published datasets (downloaded to `/opt/data`)
- **ROSE-HF/LF**, doi:10.13012/B2IDB-6128941_V1. Used here for validation (1000M3-HF).
- **INDELible-HF**, doi:10.13012/B2IDB-0900513_V1 (HMMerge data).
- Neither contains over-long sequences, because both were made by fragmenting
  full-length sequences.

### Constructed datasets (`code/make_long.py`)
The source is the MAGUS paper's ROSE 1000M4 true alignments (doi:10.13012/B2IDB-2643961_V1).
1000M4 has a median pairwise p-distance of 0.51. We first tried 1000M2, but it is close
to saturation (p = 0.69), as described in §5.4.

Each dataset is built as follows:
- 500 sequences are subsampled from one replicate.
- 10% of them (50 sequences) become **long queries**, in one of two ways:
  - **Random flank** (`rand`, +0.5×, +1× or +3× of the sequence's own length): i.i.d.
    letters with the dataset's base composition, split at random between the 5′ and
    3′ ends. In the reference every flank letter is in a column of its own, so it is
    homologous to nothing.
  - **Second domain** (`dom`): a whole sequence from a *different* ROSE replicate (a
    second "gene", about 1000 bp) is appended at the 3′ end. The second-domain pieces of
    different long queries are homologous to each other, and the reference holds their
    true alignment in a separate block of columns.
- **Control** (`m0`): the same 50 sequences, with no flank.

The remaining 450 sequences form the backbone, given as their true alignment, with a
FastTree tree. The add-methods therefore differ only in how they add queries. The long
queries are always queries, so every condition is paired with the control: same
replicate, same query set, same backbone.

Each condition has 3 replicates (R0–R2). De novo MAGUS and MAFFT L-INS-i were run on R0
only, for time.

### Scoring
- FastSP SPFN/SPFP is computed on all sequences and on the long queries only (the
  alignment induced on those 50 sequences).
- UPP and WITCH write insertion letters in lower case, meaning "not aligned". Our primary
  score masks them (FastSP `-ml`). We also report the unmasked SPFP, which is what a user
  gets by upper-casing the output.

### Validation of the harness
Target: WITCH-NG Table 2 on 1000M3-HF (all sequences, same backbone for every method):
UPP 0.054 / 0.038 and WITCH 0.048 / 0.039 (SPFN / SPFP, 20 replicates).
We reran this with a MAGUS backbone on the 500 full-length sequences
(`code/validate.sh`); the results are in §5.1.

## 4. Methods run

| name | what |
|---|---|
| upp | UPP from SEPP 4.4.0 (bioconda), `-a backbone -t tree`, defaults (HMMER 3.1b2) |
| witch | WITCH (witch-msa 1.0.10, pip), `-b -e -q`, defaults |
| emma | EMMA (GitHub c5shen/EMMA @68be952), `-b -e -q`, defaults |
| mafft-add | `mafft --add queries backbone` (MAFFT 7.505) |
| mafft-addlong | `mafft --addlong` (R0 only) |
| mafft | MAFFT L-INS-i de novo on all 500 sequences (R0 only); `--auto` picks FFT-NS-2, which fails on ROSE (§5.4) |
| magus | MAGUS de novo, defaults (R0 only) |
| **upp-tfa / emma-tfa** | **pilot fix** (`code/trim.py`), described below |

UPP2 and HMMerge were not run. UPP2 is not in the bioconda SEPP package, and the HMMerge
repository did not clone anonymously.

**Pilot fix, "trim + flank realignment" (tfa):**
1. Build one profile HMM on the whole backbone and `hmmsearch` the queries against it.
2. Cut each query to the span from its smallest envelope start to its largest envelope
   end, plus 10 letters. An overhang shorter than 50 letters is not cut.
3. Add the trimmed queries with the base method (UPP or EMMA).
4. Put the cut flanks back on their side:
   - Flanks that are homologous to other flanks (all-vs-all `nhmmer`, E < 1e-5) are
     aligned to each other with MAFFT and placed as one block of columns.
   - All other flank letters stay unaligned, each in its own column.

## 5. Results

The tables are generated by `code/aggregate.py` and stored in `results/tables.md`; the
per-replicate rows are in `results/scores.csv`. Each cell is the mean SPFN / SPFP. For
UPP, WITCH and the -tfa variants, lower-case insertion letters are masked.

- Data: ROSE 1000M4, n = 500, 50 long queries, true backbone.
- Replicates: R0–R2 for every condition, plus R3–R5 for the second domain.

### 5.1 Validation
WITCH-NG Table 2, 1000M3-HF, all sequences, MAGUS backbone (`results/validation.md`):

| method | ours (2 replicates), SPFN / SPFP | published (20 replicates), SPFN / SPFP |
|---|---|---|
| UPP | 0.055 / 0.044 | 0.054 / 0.038 |
| WITCH | 0.052 / 0.044 | 0.048 / 0.039 |

SPFN reproduces within 0.004. SPFP is about 0.005 high, and the published ordering of the
methods holds.

### 5.2 Long queries only (the alignment induced on the 50 lengthened sequences)

| method | control (no flank) | random +0.5× | random +1× | random +3× | second domain |
|---|---|---|---|---|---|
| upp | 0.003 / 0.001 | 0.003 / 0.001 | 0.003 / 0.001 | 0.003 / 0.001 | **0.502** / 0.001 |
| witch | 0.003 / 0.001 | 0.003 / 0.001 | 0.003 / 0.001 | 0.003 / 0.001 | **0.502** / 0.001 |
| emma | 0.001 / 0.001 | 0.001 / 0.009 | 0.001 / 0.017 | 0.001 / **0.044** | **0.487** / 0.002 |
| mafft --add | 0.002 / 0.001 | 0.017 / 0.176 | 0.030 / 0.308 | 0.020 / **0.551** | 0.040 / 0.028 |
| mafft --addlong (R0) | 0.002 / 0.002 | 0.002 / 0.172 | 0.002 / 0.302 | 0.002 / 0.550 | 0.033 / 0.028 |
| **upp-tfa** | 0.003 / 0.001 | 0.003 / 0.001 | 0.003 / 0.001 | 0.003 / 0.001 | **0.030 / 0.008** |
| **emma-tfa** | 0.001 / 0.001 | 0.001 / 0.002 | 0.001 / 0.001 | 0.001 / 0.001 | **0.029 / 0.008** |
DENOVO_ROWS

**UPP and WITCH output used as-is.** If the output is upper-cased, insertion letters count
as aligned. Long-query SPFP then becomes:

| method | control | random +0.5× | random +1× | random +3× | second domain |
|---|---|---|---|---|---|
| upp | 0.001 | 0.246 | 0.393 | 0.651 | 0.434 |
| witch | 0.001 | 0.255 | 0.407 | 0.664 | 0.435 |
| upp-tfa | 0.001 | 0.020 | 0.020 | 0.020 | 0.010 |

**Over all 500 sequences** the effects are diluted about 50-fold, because only 10% of the
sequences are long and the backbone is true. For example, on the second domain
(all-sequence SPFN) UPP scores 0.010 and upp-tfa 0.001, and with random +3× flanks
`mafft --add` has an SPFP of 0.013. The full table is in `results/tables.md`.

### 5.3 Paired tests
Error is (SPFN + SPFP) / 2 on the long queries, compared by replicate. In W/T/L, W means
the first-named method is better.

| comparison | conditions | n | mean diff | W/T/L | Wilcoxon p (two-sided) |
|---|---|---|---|---|---|
| upp-tfa vs upp | second domain | 6 | −0.232 | 6/0/0 | 0.031 |
| upp-tfa vs upp | random flanks (+0.5/1/3×) | 9 | −0.0000 | 1/8/0 | 1.0 |
| upp-tfa vs upp | control | 3 | 0 | 0/3/0 | – |
| emma-tfa vs emma | second domain | 6 | −0.226 | 6/0/0 | 0.031 |
| emma-tfa vs emma | random flanks | 9 | −0.011 | 9/0/0 | 0.0039 |
| emma-tfa vs emma | control | 3 | 0 | 0/3/0 | – |
| upp-tfa vs mafft --add | second domain | 6 | −0.015 | 6/0/0 | 0.031 |
| emma-tfa vs mafft --add | random flanks | 9 | −0.182 | 9/0/0 | 0.0039 |
| witch vs upp | all 5 conditions | 18 | +0.0001 | 1/0/17 | 0.0002 (difference negligible) |

With 6 pairs, 0.031 is the smallest two-sided p a Wilcoxon test can give. Over all 18
datasets combined, the long-query p-values are 0.016 for upp-tfa vs upp (7/11/0) and
6×10⁻⁵ for emma-tfa vs emma (15/3/0) (`results/tables.md`).

Flank detection worked in both directions:
- On the random-flank datasets, the nhmmer homology test flagged no flank at all
  (0 realigned in 9 datasets).
- On the second domain it flagged 50/50 flanks in 5 replicates and 49/50 in the sixth.

### 5.4 Things that went wrong or surprised us (tuning runs, `results/tuning_runs.txt`)
- **ROSE 1000M2 is close to saturation** (median pairwise p-distance 0.69).
  - With that dataset as the second domain, nhmmer could detect only 10–12 of the 50
    pieces as homologous to each other, even with an iterated profile search. The fix
    then does nothing.
  - So the fix depends on the extra region being detectably similar across queries. We
    moved the grid to 1000M4 (p = 0.51).
- **MAFFT's fast modes collapse on ROSE.**
  - `mafft --auto` selects FFT-NS-2, which gave SPFN 0.99 on 500 sequences of 1000M2
    (0.75 on 60 sequences, where L-INS-i gave 0.12).
  - On 1000M2, plain `mafft --add` gave a long-query SPFN of 0.69 even without flanks.
  - De novo MAFFT here therefore means L-INS-i.
- **A first version of trimming hurt.** With a 10-letter pad and no minimum cut, cutting
  queries to the HMM envelope lost real homologous ends: EMMA's long-query SPFN went from
  0.007 to 0.036 on the 1000M2 control. Only cutting overhangs of at least 50 letters
  removed the damage; the final grid shows 0 losses on the controls.

## 6. Runtime

Mean seconds for the add step (backbone given). Replicates R1–R2 ran 4 single-thread jobs
at once on 4 cores. `mafft --addlong` and the de novo methods are R0 only, and their R0
jobs ran in parallel with others.

| method | control | random +0.5× | random +1× | random +3× | second domain |
|---|---|---|---|---|---|
| upp | 121 | 136 | 170 | 253 | 166 |
| witch | 190 | 238 | 340 | 543 | 317 |
| emma | 81 | 93 | 129 | 229 | 109 |
| mafft --add | 2 | 3 | 7 | 30 | 4 |
| mafft --addlong | 215 | 288 | 653 | 807 | 434 |
| upp-tfa | 119 | 124 | 150 | 170 | 192 |
| emma-tfa | 81 | 85 | 105 | 129 | 158 |
DENOVO_TIME

- **Long flanks cost the HMM methods time.** UPP takes 2.1× and WITCH 2.9× longer at
  +3× than in the control.
- **Trimming first recovers most of that.** upp-tfa is 1.5× faster than UPP at +3×.
- **Flank realignment costs about 25–50 s** on the second domain (nhmmer plus MAFFT on
  50 sequences of about 1 kb).

## 7. Verdict: **unclear, leaning promising as a narrow benchmarking project**

### What the pilot shows
1. **Random flanks**, i.e. extra sequence homologous to nothing.
   - UPP and WITCH are essentially immune. SPFN/SPFP on the long queries is unchanged
     from the control up to +3× flanks, because the flanks become insertion letters.
   - Pitfall: in UPP's and WITCH's own output, those insertion letters are stacked into
     shared columns. Anyone who upper-cases the output, or does not mask lower case,
     gets a long-query SPFP of 0.25–0.66.
   - EMMA (MAFFT-linsi `--add` inside subsets) and `mafft --add`/`--addlong` really do
     align flank letters to one another: long-query SPFP is 0.044 for EMMA and 0.55 for
     MAFFT at +3×.
   - Trimming queries to the backbone-HMM envelope removes that error for EMMA
     (9/0/0 wins, p = 0.004). For UPP it changes nothing, since UPP was already fine.
2. **Second domain**, an extra region homologous among the long queries but absent from
   the backbone. This is a clear, large failure that is the same for every
   backbone-anchored adder:
   - UPP, WITCH and EMMA all leave the extra region unaligned, so the long queries lose
     about half their homologies (SPFN 0.49–0.50).
   - `mafft --add`, which also aligns queries to each other, keeps most of it
     (SPFN 0.040).
   - The pilot fix (trim → add → detect homologous flanks with nhmmer → align them with
     MAFFT) brings UPP and EMMA down to SPFN 0.03 / SPFP 0.008. That is 6/6 wins over its
     base method (p = 0.031, the minimum for n = 6) and 6/6 wins over `mafft --add`.
DENOVO_VERDICT
3. **Literature gap.** No UPP-family paper benchmarks queries longer than the family.
   Every "length heterogeneity" dataset we found consists of fragments.

### Why only "unclear"
- **Partly true by construction.** The second-domain scenario was built so that the
  backbone lacks the region, and the fix was designed for that case. The benchmark
  shows the failure mode exists. It does not show how often it happens in real data.
- **Detection limit.** On a divergent family (ROSE 1000M2, p ≈ 0.69) the fix cannot find
  the shared region, so it does nothing (§5.4). Its benefit is limited to extra regions
  that are detectably similar across queries.
- **Easy setting.** The backbone was the true alignment and 1000M4 is easy (control
  SPFN 0.003), so absolute errors are tiny everywhere except where a method is
  structurally wrong. With estimated backbones and harder models, the flank effects may
  be swamped by ordinary error, or may interact with it; that is untested.
- **Diluted overall.** Over all sequences the effects are ≤ 0.01, because only 10% of
  sequences are long. A paper would need to argue that per-query accuracy is what
  matters, for example for placement or for downstream trees on those taxa.
- **Thin algorithmic novelty.** The fix is about 120 lines of glue around
  hmmsearch, nhmmer and MAFFT. The contribution would be mainly the benchmark and the
  characterisation.

### What a 4-week project would look like
- **Week 1: realistic over-long data.** This is the critical step. Candidates:
  - (a) rRNA-operon long reads (PacBio/ONT 16S–ITS–23S, ~4.5 kb) added to a 16S backbone,
    for example CRW 16S. The 23S part is a real "second domain".
  - (b) Pfam/HomFam multi-domain proteins added to single-domain seed alignments.
  - (c) Genomic contigs carrying a gene plus flanks.
  - Build the reference for the extra region from a curated 23S alignment or Pfam.
- **Week 2: benchmark.** Run UPP, UPP2, WITCH-NG, EMMA, `mafft --add`/`--addlong`, MAGUS
  and MAFFT, with *estimated* backbones. Vary the fraction of long queries (5–50%),
  their extra length, divergence (ROSE M1–M4) and sequencing error. Score long queries
  separately.
- **Week 3: method.** Harden tfa:
  - jackhmmer-style iterative homology detection across flanks;
  - several extra regions per side;
  - handling of internal (not only terminal) non-backbone insertions, by realigning
    long insertion runs among queries;
  - optionally, merging the realigned block into backbone columns when part of it is
    homologous to the backbone.
- **Week 4:** ablations, runtime scaling to 10k or more queries, write-up.

### Risks
- **No convincing real dataset.** If real over-long sequences rarely share an extra
  region, the second-domain result is an artefact of the simulation, and the random-flank
  result is mostly "UPP is already fine".
- **Prior unpublished work.** Shen's thesis (hdl:2142/129918) or Warnow-lab manuscripts
  may already cover over-long queries. Check with the instructor in week 1.
- **Tool fragility.**
  - The bioconda SEPP needs a hand-made `upp.config`, and its HMMER is 3.1b2.
  - UPP2 and HMMerge were not installed here.
  - MAFFT's fast modes fail on ROSE.
  - Expect a day of setup.
- **Compute.** MAGUS and MAFFT L-INS-i de novo at 500 × 1–4 kb take more than 30 min
  each on 1 core. Larger de novo baselines need the course cluster.

## 8. Reproduce
```
bash cs581/code/setup.sh                      # MAFFT, MAGUS, HMMER, FastSP, MAGUS data
# UPP (SEPP 4.4), WITCH, EMMA: micromamba env "add" (python 3.9, sepp, dendropy 4.5.2) + pip witch-msa
#   + git clone c5shen/EMMA to /opt/src/EMMA; SEPP needs share/sepp/sepp/upp.config (tool paths)
bash cs581/lenhet/code/grid.sh 4 1            # main grid (restartable)
bash cs581/lenhet/code/extra_dom.sh           # second-domain R3-R5
bash cs581/lenhet/code/validate.sh            # WITCH-NG Table 2 check
python3 cs581/lenhet/code/aggregate.py        # -> results/tables.md, scores.csv, paired.json
```
