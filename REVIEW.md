# DadaS: review guide

DadaS is a fork of [DADA2](https://github.com/benjjneb/dada2) 1.36.0 (base commit `b61fe2e`).
Its goals:

1. **Same results, faster.** With default options, every stage gives output identical to DADA2.
2. **More robust.** Fixes DADA2 failures, reproducible taxonomy, no R API calls from worker threads.
3. **Optional, off-by-default extensions** for issues we found in our own data: indel errors and
   binned quality scores.

This document lists every change, why it was made, how it was verified, the measured effect, and the
open questions where a second opinion is most useful.

## At a glance

Seven sediment samples (NextSeq 2000, V3–V4, 2.05 M read pairs, binned qualities), 2 cores, same
machine. Times in seconds; all DadaS outputs identical to DADA2 (taxonomy: see section 4).

| Stage | DADA2 1.36.0 | DadaS 1.36.0.9000 | DadaS 1.36.0.9001 | 9001 vs DADA2 |
|---|---|---|---|---|
| 1 `filterAndTrim` | 113 | 127 | 106–191 ¹ | ≈ 1× |
| 2 `derepFastq` | 54 | 46 | 42 | ≈ 1× |
| 3 `learnErrors` (F+R) | 15,462 | 9,686 | **5,031** | **3.1×** |
| 4 `dada` (F+R) | 3,797 | 2,304 | 2,191 | 1.7× |
| 5 `mergePairs` | 102 | 102 | 89 | 1.1× |
| 7 `removeBimeraDenovo` | 2,264 | 1,082 | 1,066 | 2.1× |
| 8 `assignTaxonomy` | 752 | 765 | 722 | 1.0× |
| **Whole pipeline** | **22,544 (6 h 16 min)** | **14,113 (3 h 55 min)** | **9,334 (2 h 36 min)** | **2.4×** |

¹ Filtering is I/O-bound and noisy on this machine: 191 s in the full run, 106–110 s for both
versions when repeated back to back.

Where the time goes (`dada()` on one sample, 45,000 uniques, 1 thread, gperftools profile):
Needleman–Wunsch alignment 67 %, building the substitution list 11 %, λ 6 %, k-mer screen 6 %,
reassigning uniques between partitions 4 %. `learnErrors` repeats this up to 10 times.

## How to check it

```r
# Install next to dada2 (the package is named dadaS, so both can coexist). Always build clean:
# R CMD INSTALL reuses stale object files in src/, which once produced a 7x slower build here.
R CMD INSTALL --preclean .
# Unit tests (55 expectations, a few minutes). The comparison against DADA2 runs only if dada2 is installed.
Rscript -e 'testthat::test_local(".")'
```

- `tests/testthat/test-identical-to-dada2.R` runs both packages on the bundled example data and
  checks that `filterAndTrim`, `derepFastq`, `learnErrors`, `dada` (serial and parallel),
  `mergePairs` and `removeBimeraDenovo` give identical output.
- `tests/testthat/test-fixes.R` has one or more regression tests per bug fix (section 3) and
  checks that the alignment cache and the chimera parent search do not change results.
- `R CMD check` (R 4.3.3, Linux): 0 errors, 0 warnings, 8 notes, all inherited from DADA2 1.36.0
  (e.g. `ShortRead:::.set_omp_threads`, `.travis.yml`, undefined ggplot2 globals).

## 1. Package rename

The package is renamed to `dadaS` (`DESCRIPTION`, `NAMESPACE`, `inst/include/dadaS*.h`,
`package="dadaS"` in examples, regenerated `RcppExports`) so it can be installed alongside
`dada2` and compared against it. The original authors are kept in `DESCRIPTION`, as the LGPL-2
license requires.

## 2. Performance changes (results unchanged)

| Change | Files | Measured effect | Why it is safe |
|---|---|---|---|
| **Alignment cache across `learnErrors` rounds.** The alignment of a center and a unique does not depend on the error model, only λ does. For each (center, unique) pair the cache keeps the k-mer screen outcome and the number of substitutions (1 byte). In later rounds a comparison is skipped when that proves it would be discarded: a failed k-mer screen gives λ = 0 exactly as recomputing it would, and otherwise λ ≤ `emax`^nsubs (`emax` = largest substitution rate at the unique's quality scores), so if `emax`^nsubs × reads < `E_minmax` the comparison could never be stored. | `src/dada.h` (`AlignCache`), `src/cluster.cpp`, `src/Rmain.cpp`, `R/dada.R` | `learnErrors` on one sample: 520 s → 257 s (2.0×). Stage 3 on 7 samples: 9,686 s → 5,031 s (1.9×). Memory: 1 byte per (center, unique) pair. | Only comparisons whose result is discarded are skipped. Tested with and without the cache, 1 and 2 threads, and with pseudo-pooling; identical on the sediment data. `options(dadaS.alignCache = FALSE)` turns it off. |
| **Vectorized alignment.** With GCC < 14, `-O2` did not vectorize the "vectorized" Needleman–Wunsch loops. The anti-diagonal loop now walks `s1` reversed so the diagonal scores are computed inside the same loop, and vectorization is forced for that file only. | `src/nwalign_vectorized.cpp` | Part of the 1.6× of 9000. | Integer arithmetic, same operations and tie-breaking order. |
| **Substitutions built in the traceback.** The aligner no longer builds two aligned strings that `al2subs` then scans twice; DP buffers are reused per thread chunk instead of allocated per alignment. | `src/nwalign_vectorized.cpp` (`nwalign_vectorized2_sub`), `src/nwalign_endsfree.cpp` | `dada()` on one sample: 85.5 s → 82.2 s (~4 %). | Identical `Sub` on 724,000 pairs (bands 16, 4, unbanded; indels; different lengths), checked with `C_check_fused_subs`. |
| **Scalar Poisson CDF.** `R::ppois` instead of `Rcpp::ppois`, which allocated two R vectors per call. | `src/pval.cpp` | Part of the 1.6× of 9000. | Same Rmath routine. |
| **Serial path error matrix.** The single-thread `b_compare` received the `Rcpp::NumericMatrix` by value on every comparison; it now copies it once into a C array, as the parallel path did. | `src/cluster.cpp` | Serial path faster than parallel with 1 thread before this change. | Same product in the same order. |
| **NEON SIMD for k-mer distances** on ARM (Apple Silicon, Graviton), which used scalar code. | `src/kmers_simd.h`, `src/kmers.cpp` | K-mer screen 1.54× (Mac M-series, clang), 1.66× (Linux x86-64, SSE2 vs scalar). | Bit-identical; checked against the scalar reference on 80,910 uniques. |
| **Vectorized quality validation in `dada()`** (was a row-by-row `sapply`). | `R/dada.R` | ~20 % of `dada()` time on large samples. | Same checks, same messages. |
| **Chimera parent search.** `isBimeraDenovoTable` walks each sample's ASVs by decreasing count and stops at the first one that is too small, instead of scanning all ASVs for every candidate. | `src/chimera.cpp` | No change with 7 samples (dominated by alignments); avoids O(samples × ASVs²) checks (~10¹² with 500 samples and 50,000 ASVs). | The flags are maxima over the candidate parents, so the visiting order does not matter. Identical to DADA2 on a 30-sample table with three `minFoldParentOverAbundance` values. |
| **`mergePairs` pair counting** without the dense forward × reverse `table()`. | `R/paired.R` | 20,000 ASVs per side: 1.6 GB → a few MB. | Same integer counts (identical to DADA2 in the tests). |
| **`assignTaxonomy` bootstrap uniforms** read in place instead of copied to a C array. | `src/taxonomy.cpp` | −365 MB peak memory for 8,781 ASVs (grows with the number of ASVs). | Same values, read-only from the worker threads. |
| **Threads on HPC nodes.** `availableThreads()` = min(`SLURM_CPUS_PER_TASK`, CPU affinity of the process, `detectCores()`), used by `dada`, `learnErrors`, `filterAndTrim`, chimeras and taxonomy when `multithread = TRUE`. DADA2 used every core of the node. | `R/misc.R` and callers | Avoids oversubscribing shared nodes. | The thread count never changes DadaS results. |

Measured but not adopted: compiling with AVX2 (`-mavx2`) gave no gain in the aligner kernel, and
reassigning uniques between partitions (`b_shuffle2`), suspected to be quadratic, takes only 4 % of
the time.

## 3. Bug fixes and robustness

Each bug was reproduced before fixing it; the regression tests fail on the previous version.

| Problem | Effect before | Fix | Files |
|---|---|---|---|
| `USE_QUALS = FALSE` | `dada()` always passed `use_quals = TRUE` to C++ and stopped with "Rounded quality exceeded range"; `noqualErrfun` failed on a one-column matrix. | Pass the option; `drop = FALSE` in `noqualErrfun`; the transition-matrix write is bounds-checked. | `R/dada.R`, `R/errorModels.R`, `src/error.cpp`, `src/Rmain.cpp` |
| `assignTaxonomy` with a sequence shorter than 8 nt | The whole call failed ("negative length vectors"): unsigned underflow. Sequences longer than ~10,000 nt overflowed fixed stack arrays. | Short sequences skipped when sizing; k-mer buffers on the heap. | `src/taxonomy.cpp` |
| `setDadaOpt` | A misplaced parenthesis accepted any value of any class (e.g. `BAND_SIZE = "sixteen"`). | Validates the class; numeric options accept integer or double, as before. | `R/dada.R` |
| `mergePairs` | Changed the session's alignment scores and restored them afterwards; an error in between left `MATCH = 1, MISMATCH = -64, GAP_PENALTY = -64` for every later `dada()`. | Scores passed directly to `nwalign`. | `R/paired.R` |
| Memory when `dada()` fails or is interrupted | Everything allocated for the sample stayed allocated: +1,851 MB after 10 failed calls on one sediment sample; 30 failures exhausted memory. Now +89 MB (R's own, constant). | All allocations released on exit, also on errors and Ctrl-C. | `src/Rmain.cpp` |
| `birth_from` for prior-driven partitions | Uninitialized memory in the output. | Set. | `src/cluster.cpp` |
| Errors in worker threads | `Rcpp::stop` from an `RcppParallel::Worker` could crash R. | `dada_stop` throws `std::runtime_error`, converted to an R error on the main thread. | `src/dada.h`, `src/taxonomy.cpp` |
| Out-of-bounds reads | All-gap alignments in `C_eval_pair` and the chimera overlap scan. | Bounds checked first. | `src/evaluate.cpp`, `src/chimera.cpp` |
| NULL dereference | `comp.hamming = sub->nsubs` with `sub == NULL`. | Guarded. | `src/cluster.cpp` |

Also: DADA2 1.36 segfaults in single-thread mode when the error matrix has fewer columns than the
quality scores (the serial λ reads outside the matrix); DadaS stops with an error message. The
documentation of `minQ` and `minParentAbundance` was corrected to describe what the code does (no
behavior change).

## 4. Taxonomy reproducibility (the only default-output difference from DADA2)

DADA2 created a `std::random_device` + `mt19937` on every bootstrap call (~101 per sequence): slow,
and `set.seed()` was ignored, so bootstrap values changed from run to run. DadaS draws one seed per
sequence from R's RNG and uses one generator per sequence, so results respect `set.seed()` and do
not depend on the number of threads. For a given seed, bootstrap values differ from DADA2's: in the
sediment run, 165 of 52,686 bootstrap values differed by at most 5/100, and 1 of 8,781 genus
assignments changed (DADA2 itself gave 49, 48 and 54 for that genus in three runs with the same seed).

## 5. New, optional features (off by default)

### `INDEL_RATE` (`setDadaOpt(INDEL_RATE = x)`)

DADA2 ignores indels when computing λ, so variants that differ only by indels are never split.
With `INDEL_RATE > 0`, `al2subs` counts internal indel events (gap runs, terminal gaps excluded), λ
is multiplied by `INDEL_RATE ^ nindel`, and each event adds one to the hamming distance.
`INDEL_RATE = 0` reproduces DADA2 exactly. Helpers: `indelOnlyVariants()` and `estimateIndelRate()`.

- **Simulation** (V4, 240,000 read pairs, 1 % of molecules with a PCR homopolymer indel):

  | | Indel variants recovered | Base taxa | False ASVs |
  |---|---|---|---|
  | DADA2 | 0/31 | 281/282 | 0 |
  | `INDEL_RATE` 1e-2 | 11/31 | 282/282 | 0 |
  | `INDEL_RATE` 1e-3 | 16/31 | 282/282 | 0 |
  | `INDEL_RATE` 1e-4 | 18/31 | 282/282 | 0 |
  | `INDEL_RATE` 1e-5 | 22/31 | 282/282 | 0 |

- **Sediment data**, final table: with 2e-3, 8,781 → 9,711 ASVs (1,209 new, 279 lost); with per-read
  rates from `estimateIndelRate` (F 7e-3, R 1.8e-3), 9,593 ASVs (1,074 new, 262 lost). The new ASVs
  hold < 1 % of reads, rarely recur across samples (2.3 % vs 18.6 % of original ASVs of the same
  abundance), and only 4 differ from an existing ASV by indels alone; 855 of 1,074 are more than 10
  differences away from any original ASV. In these data the option mostly adds low-abundance noise.

Open questions: a gap run that switches sequence counts as one event; the rate does not depend on
homopolymer length; indel rates are not learned by `learnErrors`; indel events also enter
`MIN_HAMMING` and the greedy shortcut.

### `loessErrfunMonotone(trans)`

The community variant of `loessErrfun` for binned qualities
([dada2#1307](https://github.com/benjjneb/dada2/issues/1307)): weights `log10(total)`, `span = 2`,
`degree = 1`, forced non-increasing in quality. On one sediment sample (forward reads, final ASVs):

| Error function | ASVs | Jaccard vs `loessErrfun` |
|---|---|---|
| `loessErrfun` (default) | 2,210 | — |
| `makeBinnedQualErrfun(c(12, 24, 40))` | 2,043 | 0.887 |
| `loessErrfunMonotone` | 2,842 | 0.737 |

**Root cause.** `derepFastq` averages qualities per unique sequence. With three quality values that
average is coupled with abundance (rare uniques sit at Q40), so the fitted Q40 error rate is inflated
(1.57 % at Q40 vs 0.012 % at Q39 in our data) and `learnErrors` never converges. A per-read quality
error model is being added as an option (see `NEWS`).

## 6. Left out on purpose (they would change existing outputs)

| Item | What would change |
|---|---|
| Input-order dependence: `b_bud` assumes the first unique is the partition center | ASVs for derep objects built by hand in another order (`derepFastq` always sorts) |
| `collapseNoMismatch` reorders by abundance even with `orderBy = "nsamples"` | Column order |
| phiX filter skips one extra position when counting non-overlapping words | Which reads are removed |
| −9,999 band sentinel in the scalar aligners | Alignments of reads over ~2,500 nt |
| `birth_qave` is `NaN` for a partition born with 0 substitutions | One diagnostic value |
| `KDIST_CUTOFF` 0.30 instead of 0.42 | On one sample, denoising 82 s → 37 s (2.2×), 1,179 → 1,183 ASVs; all samples under evaluation |
| `learnErrors` convergence with a tolerance (binned data never converge exactly) | The error model |

## Where to look first

1. `src/cluster.cpp`, `src/dada.h`: `AlignCache` and `INDEL_RATE`, the two algorithmic changes.
2. `src/nwalign_vectorized.cpp`, `src/nwalign_endsfree.cpp`: alignment kernels and the fused traceback.
3. `src/Rmain.cpp`: memory ownership and cleanup in `dada_uniques`.
4. `src/taxonomy.cpp`: RNG, buffers, thread-safety.
5. `R/dada.R`, `R/paired.R`, `R/misc.R` (`availableThreads`), `R/indels.R`, `R/errorModels.R`.

The staged pipeline, the read simulator and the analysis scripts behind the numbers above live on
the `dadaS` branch under `workflow/` (comments in Spanish; not part of the package).
