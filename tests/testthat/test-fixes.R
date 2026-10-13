# Regression tests for the fixes made after the code review (2026-10).
fnF <- system.file("extdata", "sam1F.fastq.gz", package = "dadaS")
fnR <- system.file("extdata", "sam1R.fastq.gz", package = "dadaS")
drpF <- derepFastq(fnF)
drpR <- derepFastq(fnR)

test_that("USE_QUALS = FALSE runs (it failed with 'Rounded quality exceeded range')", {
  e1 <- matrix(tperr1[, 41], ncol = 1, dimnames = list(rownames(tperr1), "0"))
  dd <- dada(drpF, err = e1, USE_QUALS = FALSE, verbose = 0)
  expect_s4_class(dd, "dada")
  expect_equal(ncol(dd$trans), 1L)
  expect_gt(nrow(dd$clustering), 0)
  sc <- dada(drpF, err = NULL, selfConsist = TRUE, USE_QUALS = FALSE, verbose = 0)
  expect_equal(ncol(getErrors(sc)), 1L)
})

test_that("setDadaOpt rejects values of the wrong class and still accepts numeric for integer options", {
  old <- getDadaOpt()
  on.exit(setDadaOpt(old))
  expect_warning(setDadaOpt(BAND_SIZE = "sixteen"), "not set")
  expect_equal(getDadaOpt("BAND_SIZE"), old$BAND_SIZE)
  expect_warning(setDadaOpt(USE_KMERS = 1), "not set")
  expect_silent(setDadaOpt(MATCH = 5))   # numeric for an integer option: accepted, as before
  expect_silent(setDadaOpt(HOMOPOLYMER_GAP_PENALTY = NULL))
})

test_that("an error inside mergePairs no longer leaves the session with the merge scores", {
  dF <- dada(drpF, err = tperr1, verbose = 0); dR <- dada(drpR, err = tperr1, verbose = 0)
  before <- getDadaOpt(c("MATCH", "MISMATCH", "GAP_PENALTY"))
  expect_error(mergePairs(dF, drpF, dR, drpR, band = "x"))
  expect_identical(getDadaOpt(c("MATCH", "MISMATCH", "GAP_PENALTY")), before)
})

test_that("assignTaxonomy handles sequences shorter than the k-mer size", {
  seqs <- getSequences(system.file("extdata", "example_seqs.fa", package = "dadaS"))[1:5]
  train <- system.file("extdata", "example_train_set.fa.gz", package = "dadaS")
  set.seed(1); ref <- assignTaxonomy(seqs, train)
  set.seed(1); tx <- suppressWarnings(assignTaxonomy(c(seqs, "ACGTA"), train))
  expect_true(all(is.na(tx["ACGTA", ])))
  expect_identical(unname(tx[seqs, ]), unname(ref)) # most probable assignment does not depend on the extra sequence
})

test_that("the alignment cache does not change learnErrors or pseudo-pooling", {
  on.exit(options(dadaS.alignCache = NULL))
  fls <- system.file("extdata", c("sam1F.fastq.gz", "sam2F.fastq.gz"), package = "dadaS")
  for (mt in c(FALSE, 2)) {
    options(dadaS.alignCache = TRUE);  a <- learnErrors(fls, multithread = mt, verbose = 0)
    options(dadaS.alignCache = FALSE); b <- learnErrors(fls, multithread = mt, verbose = 0)
    expect_identical(a, b)
  }
  drps <- derepFastq(fls)
  options(dadaS.alignCache = TRUE);  p1 <- dada(drps, err = tperr1, pool = "pseudo", verbose = 0)
  options(dadaS.alignCache = FALSE); p2 <- dada(drps, err = tperr1, pool = "pseudo", verbose = 0)
  for (i in seq_along(p1)) expect_identical(p1[[i]]$denoised, p2[[i]]$denoised)
})

test_that("isBimeraDenovoTable gives the same flags with any parent-search order (larger table)", {
  skip_if_not_installed("dada2")
  fls <- system.file("extdata", c("sam1F.fastq.gz", "sam2F.fastq.gz"), package = "dadaS")
  dd <- dada(derepFastq(fls), err = tperr1, verbose = 0)
  st <- makeSequenceTable(dd)
  set.seed(5)
  big <- st[sample(nrow(st), 30, replace = TRUE), ]          # 30 samples
  big[] <- as.integer(big * matrix(rpois(length(big), 3), nrow(big)))
  for (mf in c(1.5, 1, 0.5)) {
    expect_identical(isBimeraDenovoTable(big, minFoldParentOverAbundance = mf),
                     dada2::isBimeraDenovoTable(big, minFoldParentOverAbundance = mf))
  }
  expect_identical(isBimeraDenovoTable(big, allowOneOff = TRUE), dada2::isBimeraDenovoTable(big, allowOneOff = TRUE))
})

test_that("availableThreads respects the Slurm allocation", {
  old <- Sys.getenv("SLURM_CPUS_PER_TASK", NA)
  on.exit(if (is.na(old)) Sys.unsetenv("SLURM_CPUS_PER_TASK") else Sys.setenv(SLURM_CPUS_PER_TASK = old))
  Sys.setenv(SLURM_CPUS_PER_TASK = "1")
  expect_identical(availableThreads(), 1L)
  Sys.unsetenv("SLURM_CPUS_PER_TASK")
  expect_true(availableThreads() >= 1L && availableThreads() <= parallel::detectCores())
})
