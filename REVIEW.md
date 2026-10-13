# DadaS: review guide

DadaS is a fork of [DADA2](https://github.com/benjjneb/dada2) 1.36.0 (base commit `b61fe2e`).
Its goals:

1. **Same results, faster.** With default options, every stage gives output identical to DADA2.
2. **More robust.** Reproducible taxonomy and no R API calls from worker threads.
3. **Optional, off-by-default extensions** for issues we found in our own data: indel errors and
   binned quality scores.

This document lists every change, why it was made, how it was verified, and the open questions
where a second opinion is most useful.

## How to check it

```r
# Install next to dada2 (the package is named dadaS, so both can coexist)
R CMD INSTALL .
# Unit tests (32 expectations, a few minutes). The comparison against DADA2 runs only if dada2 is installed.
Rscript -e 'testthat::test_local(".")'
```

`tests/testthat/test-identical-to-dada2.R` runs both packages on the bundled example data. It
checks that `filterAndTrim`, `derepFastq`, `learnErrors`, `dada` (serial and parallel),
`mergePairs` and `removeBimeraDenovo` give identical output.

## 1. Package rename

The package is renamed to `dadaS` (`DESCRIPTION`, `NAMESPACE`, `inst/include/dadaS*.h`,
`package="dadaS"` in examples, regenerated `RcppExports`) so it can be installed alongside
`dada2` and compared against it.

The original authors are kept in `DESCRIPTION`, as the LGPL-2 license requires.

## 2. Performance changes (results unchanged)

| Change | Files | Why it is safe |
|---|---|---|
| **Vectorized alignment.** R builds with `-O2`, and with GCC < 14 that does **not** vectorize the "vectorized" Needleman-Wunsch loops, so on Linux they ran scalar. The anti-diagonal loop now walks `s1` reversed so that the diagonal scores are computed inside the same loop, and vectorization is forced for that file only (GCC attribute; no effect on clang). | `src/nwalign_vectorized.cpp` | Integer arithmetic, same operations and same tie-breaking order (`up > left > diag`, or swapped). |
| **Scalar Poisson CDF.** Uses `R::ppois` instead of the vectorized `Rcpp::ppois`, which allocated two R vectors per call. | `src/pval.cpp` | Same Rmath routine, so the values are identical. |
| **Serial path error matrix.** The single-thread `b_compare` received the `Rcpp::NumericMatrix` by value on every comparison. It now copies the matrix once into a C array, as the parallel path already did. Before this, the serial path was slower than the parallel path with one thread. | `src/cluster.cpp` | Same product in the same order. |
| **NEON SIMD for k-mer distances.** The SSE2 versions existed only for x86-64; ARM (Apple Silicon, Graviton) used scalar code. The NEON code reproduces the SSE2 lanes, saturating/wrapping adds and overflow detection exactly. | `src/kmers_simd.h`, `src/kmers.cpp`, `src/Rmain.cpp` | Bit-identical; checked against the scalar reference on 80,910 uniques. |
| **Vectorized quality validation in `dada()`.** Replaces a row-by-row `sapply`, which took ~20 % of `dada()` time on large samples. | `R/dada.R` | Same checks, same error messages. |

### Measured speed-ups

Seven sediment samples (NextSeq 2000, V3–V4, 2.05 M read pairs, 2 cores):

| Stage | DADA2 1.36.0 | DadaS | Speed-up |
|---|---|---|---|
| learnErrors (F+R) | 4 h 18 min | 2 h 41 min | ×1.60 |
| dada (F+R) | 63 min | 38 min | ×1.65 |
| removeBimeraDenovo | 38 min | 18 min | ×2.09 |
| Whole pipeline | 6 h 16 min | 3 h 55 min | ×1.60 |

Stages 1–7 gave identical output (`identical()` on every object).

K-mer filter micro-benchmark, SIMD vs scalar: ×1.54 on an Apple M-series Mac (R 4.5.1, clang) and
×1.66 on Linux x86-64.

## 3. Robustness changes

- **Reproducible taxonomy** (`src/taxonomy.cpp`). DADA2 created a `std::random_device` +
  `mt19937` on every bootstrap call (~101 per sequence). That was slow, and it ignored
  `set.seed()`, so bootstrap values changed from run to run. DadaS draws one seed per sequence
  from R's RNG and uses one generator per sequence. Results now respect `set.seed()` and do not
  depend on the number of threads.
  - **This is the only default-output difference from DADA2.** For a given seed, bootstrap values
    differ from DADA2's, which are not reproducible anyway. In the sediment run, 165 of 52,686
    bootstrap values differed by at most 5/100, and 1 of 8,781 genus assignments changed.
    Repeating DADA2 three times with the same seed gave 49, 48 and 54 for that genus.
- **Thread-safe errors** (`src/dada.h`, `dada_stop`). `Rcpp::stop` calls the R API, which is
  unsafe inside an `RcppParallel::Worker` and could crash R. `dada_stop` throws a
  `std::runtime_error`; TBB propagates it to the calling thread, where Rcpp converts it into a
  normal R error. Taxonomy workers set an error flag instead of calling R.
- **Possible NULL dereference** (`src/cluster.cpp`). `comp.hamming = sub->nsubs` is now guarded
  for `sub == NULL`.
- **C++ standard** (`src/Makevars*`). `CXX_STD = CXX11` was removed, so the package uses R's default
  standard (C++17); R >= 4.3 flags an explicit C++11 requirement in `R CMD check`.

## 4. New, optional features (off by default)

### `INDEL_RATE` (`setDadaOpt(INDEL_RATE = x)`)

DADA2 ignores indels when computing λ, the probability that a sequence arose by errors from its
partition center. Variants that differ only by indels are therefore never split, however
abundant they are.

With `INDEL_RATE > 0`:

- `al2subs` counts internal indel events: runs of gap columns, with terminal gaps excluded because
  the alignment is ends-free (`src/nwalign_endsfree.cpp`, `Sub.nindel`).
- λ is multiplied by `INDEL_RATE ^ nindel`, and each event adds one to the hamming distance
  (`src/cluster.cpp`, both serial and parallel paths).
- `INDEL_RATE = 0` reproduces DADA2 exactly, which is verified by the tests.

Helpers:

- `indelOnlyVariants(dadaObj, derep)`: lists variants that differ from their ASV only by indels.
- `estimateIndelRate(dadaObj, derep)`: returns a high quantile (default 0.99) of the
  relative abundance of single-indel variants around abundant ASVs. Forward and reverse reads give
  different values, so the workflow accepts separate values for each.

Evaluation:

- **Simulation** (V4 amplicons, 240,000 read pairs, 1 % of molecules with a PCR homopolymer indel error):

  | | Indel variants recovered | Base taxa | False ASVs |
  |---|---|---|---|
  | DADA2 | 0/31 | 281/282 | 0 |
  | `INDEL_RATE` 1e-2 | 11/31 | 282/282 | 0 |
  | `INDEL_RATE` 1e-3 | 16/31 | 282/282 | 0 |
  | `INDEL_RATE` 1e-4 | 18/31 | 282/282 | 0 |
  | `INDEL_RATE` 1e-5 | 22/31 | 282/282 | 0 |

- **Sediment data with `INDEL_RATE = 2e-3`**, final table after merging and chimera removal:
  - 8,781 → 9,711 ASVs (1,209 new, 279 lost), but the new ASVs hold only 0.9 % of reads (median 3
    reads).
  - The new ASVs rarely occur in more than one sample: 2.2 % do, against 18.6 % of the original ASVs
    with the same abundance (9–40 reads).
  - Only 4 differ from an existing ASV by indels alone.
  - With per-read rates estimated by `estimateIndelRate` (99th percentile: F 7e-3, R 1.8e-3) the
    picture is the same: 9,593 ASVs (1,074 new, 262 lost), 0.8 % of reads in new ASVs, 2.3 % of the
    9–40-read new ASVs in more than one sample.
  - So in these data, the option mostly adds low-abundance noise rather than resolving indel
    variants. Most new ASVs (855 of 1,074) are more than 10 differences away from any original
    ASV, which suggests that the change in partitioning lets through sequences that DADA2
    discarded, not that indel variants are being split. This is the main open question.

Open questions for the reviewer:

1. One gap run that switches from one sequence to the other is counted as a single event.
2. The rate does not depend on homopolymer length, although slippage is much more frequent in long
   homopolymers.
3. Indel transitions are not learned by `learnErrors`; the rate is a fixed user parameter.
4. The extra events also enter the hamming distance used by `MIN_HAMMING` and the greedy shortcut.

### `loessErrfunMonotone(trans)`

This is the community variant of `loessErrfun` for binned quality scores (NovaSeq/NextSeq), from
[dada2#1307](https://github.com/benjjneb/dada2/issues/1307): loess weights `log10(total)`,
`span = 2`, `degree = 1`, then forced non-increasing in quality.

Our evaluation on one sediment sample (forward reads, final ASVs) compared three error functions
against the default `loessErrfun`:

| Error function | ASVs | Jaccard vs `loessErrfun` |
|---|---|---|
| `loessErrfun` (default) | 2,210 | — |
| `makeBinnedQualErrfun(c(12, 24, 40))` | 2,043 | 0.887 |
| `loessErrfunMonotone` | 2,842 | 0.737 |

Recurrence across samples suggests that 60–90 % of the extra ASVs from `loessErrfunMonotone` behave
like true sequences, and the rest is noise. In a binned-quality simulation, both loess variants
recovered the same true sequences with 0 false ASVs.

**Root cause, not addressed yet.** `derepFastq` averages qualities per unique sequence. With only
three quality values, that average is strongly coupled with abundance: rare uniques sit at Q40.
As a result, the fitted error rate at Q40 is inflated (observed 1.57 % at Q40 vs 0.012 % at Q39 in
our data), and `learnErrors` hit the 10-round limit without converging. Learning errors from
per-read qualities would fix this, but it changes `derepFastq`/`dada` internals and is left for
discussion.

## Follow-up after the code review (branch `mejoras`)

A line-by-line review found bugs and hotspots; everything that would change an existing output was
deliberately left out (listed below). On the 7 sediment samples, all 8 stages are identical to the previous
version and the run takes 2 h 36 min (was 3 h 55 min; DADA2 6 h 16 min). Tests: 55, all passing.

- **Bugs fixed** (each with a regression test in `tests/testthat/test-fixes.R`): `USE_QUALS = FALSE`
  failed; `assignTaxonomy` crashed on sequences shorter than 8 nt; `setDadaOpt` accepted values of any
  class; `mergePairs` left the session with its alignment scores after an error; `dada()` leaked
  ~185 MB per failed or interrupted call; `birth_from` was uninitialized for prior births.
- **AlignCache** (`src/dada.h`, `src/cluster.cpp`): during `learnErrors` rounds, a comparison is skipped
  when the cached k-mer screen or substitution count proves it would be discarded
  (`emax^nsubs * reads < E_minmax`, where `emax` is the largest substitution rate at the unique's
  qualities). Exact by construction; `learnErrors` 1.9x faster on the sediment data.
- **Fused traceback** (`nwalign_vectorized2_sub`): builds the `Sub` directly; checked identical to
  `al2subs(nwalign_vectorized2())` on 724,000 pairs (`C_check_fused_subs`). ~4% faster.
- **`availableThreads()`**: honors `SLURM_CPUS_PER_TASK` and CPU affinity instead of `detectCores()`.
- **Left out on purpose** (they change outputs): input-order dependence of `b_bud` (inputs from
  `derepFastq` are always sorted), `collapseNoMismatch` ignoring `orderBy`, a phiX word-count off-by-one,
  the -9999 band sentinel for reads over ~2,500 nt, `birth_qave = NaN`, and `KDIST_CUTOFF` (0.30 is 2.2x
  faster in denoising with +4 of 1,179 ASVs on one sample).

## Where to look first

1. `src/cluster.cpp`, `src/nwalign_endsfree.cpp`, `src/dada.h`: INDEL_RATE, the only algorithmic
   change.
2. `src/taxonomy.cpp`: RNG and thread-safety changes.
3. `src/nwalign_vectorized.cpp`, `src/kmers_simd.h`: performance kernels.
4. `R/indels.R`, `R/errorModels.R` (end of file), `R/dada.R`: new R functions and options.

`R CMD check` (R 4.3.3, Linux): 0 errors, 0 warnings, 8 notes, all inherited from DADA2 1.36.0
(e.g. `ShortRead:::.set_omp_threads`, `.travis.yml`, undefined ggplot2 globals).

The staged pipeline, the read simulator and the analysis scripts behind the numbers above live on
the `dadaS` branch under `workflow/` (comments in Spanish; not part of the package).
