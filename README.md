# DadaS

A faster, more robust fork of [DADA2](https://github.com/benjjneb/dada2) 1.36.0.

- **Same results by default.** With default options every step (filtering, error learning,
  denoising, merging, chimera removal) gives output identical to DADA2. The exception is taxonomy
  bootstrap values, which are reproducible with `set.seed()` and independent of the thread count.
- **Faster.** On 7 NextSeq sediment samples (2.05 M read pairs, 2 cores) the whole pipeline takes
  2 h 36 min, against 6 h 16 min with DADA2 (2.4×). Most of the gain is in `learnErrors` (3.1×),
  which reuses alignments between its rounds, and in chimera removal (2.1×). Includes NEON SIMD on
  Apple Silicon / ARM.
- **More robust.** Fixes several DADA2 failures: `USE_QUALS = FALSE`, `assignTaxonomy` with very short
  sequences, option validation in `setDadaOpt`, `mergePairs` leaving the session's alignment scores
  changed after an error, memory leaked when `dada()` fails or is interrupted, and errors raised from
  worker threads.
- **HPC-ready.** With `multithread = TRUE` the number of threads follows the Slurm allocation
  (`SLURM_CPUS_PER_TASK`) and the CPU affinity of the process, see `availableThreads()`, instead of
  every core of the node.
- **Optional extensions, off by default:**
  - `setDadaOpt(INDEL_RATE = ...)` lets abundant indel-only variants form their own ASVs.
  - `indelOnlyVariants()` and `estimateIndelRate()` diagnose indel variants and estimate a rate.
  - `loessErrfunMonotone()` is an error function for binned quality scores.

The package is named `dadaS`, so it installs alongside `dada2` and the two can be compared directly.
See [REVIEW.md](REVIEW.md) for every change, its rationale and how it was verified, and [NEWS](NEWS)
for the version history.

```r
# from a clone of this repository (a clean build: stale object files in src/ are rebuilt)
install.packages(".", repos = NULL, type = "source", INSTALL_opts = "--preclean")
library(dadaS)   # same API as dada2
```

Tests: `Rscript -e 'testthat::test_local(".")'` (55 expectations; the comparison against DADA2 runs
when `dada2` is installed).

All credit for DADA2 goes to its authors (Callahan et al., 2016). The original README follows.

---


[![Build Status](https://github.com/benjjneb/dada2/actions/workflows/pages/pages-build-deployment/badge.svg)](https://github.com/benjjneb/dada2/actions/workflows/pages/pages-build-deployment)

## dada2

Exact sample inference from high-throughput amplicon data. Resolves real variants differing by as little as one nucleotide. Visit [the DADA2 website](https://benjjneb.github.io/dada2/index.html) for the most detailed and up-to-date documentation.

### Installation

The dada2 package binaries are available through Bioconductor:

```S
## try http:// if https:// URLs are not supported
if (!requireNamespace("BiocManager", quietly=TRUE))
    install.packages("BiocManager")
BiocManager::install("dada2")
```

In order to install dada2 from source (and get the latest and greatest new features) see our [installation from source instructions](https://benjjneb.github.io/dada2/dada-installation.html).

### Documentation

The [tutorial walkthrough of the DADA2 pipeline on paired end Illumina Miseq data](https://benjjneb.github.io/dada2/tutorial.html). 

The [dada2 R package manual](https://www.bioconductor.org/packages/3.6/bioc/manuals/dada2/man/dada2.pdf).

Further documentation is available on [the DADA2 front page](http://benjjneb.github.io/dada2/). 

### DADA2 Articles

[DADA2: High resolution sample inference from Illumina amplicon data. Nature Methods, 2016.](http://dx.doi.org/10.1038/nmeth.3869) [(Open Access link.)](http://rdcu.be/ipGh)

[Bioconductor workflow for microbiome data analysis: from raw reads to community analyses. F1000 Research, 2016.](https://f1000research.com/articles/5-1492)

[Exact sequence variants should replace operational taxonomic units in marker-gene data analysis. ISMEJ, 2017.](http://dx.doi.org/10.1038/ismej.2017.119)

[High-throughput amplicon sequencing of the full-length 16S rRNA gene with single-nucleotide resolution. Nucleic Acids Research, 2019.](http://dx.doi.org/10.1093/nar/gkz569)

### Other Resources

Planned feature improvements are publicly catalogued at the main DADA2 development site on github, specifically on the "Issues" page for DADA2:

https://github.com/benjjneb/dada2/issues

If the feature you are hoping for is not listed, you are welcome to add it as a feature request "issue" on this page. This request will be publicly available and listed on the page.

Bugs and difficulties in using DADA2 are also welcome on [the issue tracker](https://github.com/benjjneb/dada2/issues).
