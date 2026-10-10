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
