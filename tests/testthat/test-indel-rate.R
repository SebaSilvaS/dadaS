# Synthetic data: an abundant sequence and a variant that differs only by a 1-bp deletion in a
# homopolymer. DADA2 (INDEL_RATE = 0) merges them; with INDEL_RATE > 0 the variant is its own ASV.
make_fastq <- function(seqs, counts) {
  f <- tempfile(fileext = ".fastq.gz"); con <- gzfile(f, "w")
  k <- 0
  for (i in seq_along(seqs)) for (j in seq_len(counts[i])) {
    k <- k + 1
    writeLines(c(paste0("@r", k), seqs[i], "+", strrep("I", nchar(seqs[i]))), con)
  }
  close(con); f
}
set.seed(3)
center <- paste(sample(c("A", "C", "G", "T"), 120, TRUE), collapse = "")
center <- paste0(substr(center, 1, 60), "AAAAA", substr(center, 61, 120))   # homopolymer in the middle
variant <- paste0(substr(center, 1, 60), "AAAA", substr(center, 66, 125))    # 1-bp deletion in it
fq <- make_fastq(c(center, variant), c(1000, 100))
drp <- derepFastq(fq)
err <- tperr1

test_that("INDEL_RATE = 0 reproduces DADA2: the indel-only variant is merged", {
  setDadaOpt(INDEL_RATE = 0)
  dd <- dada(drp, err = err, verbose = 0)
  expect_equal(nrow(dd$clustering), 1L)
  v <- indelOnlyVariants(dd, drp)
  expect_equal(nrow(v), 1L)
  expect_equal(v$n_indels, 1L); expect_equal(v$variant_reads, 100L); expect_equal(v$fraction, 100 / 1100)
})

test_that("INDEL_RATE > 0 resolves the indel-only variant", {
  on.exit(setDadaOpt(INDEL_RATE = 0))
  setDadaOpt(INDEL_RATE = 1e-3)
  dd <- dada(drp, err = err, verbose = 0)
  expect_setequal(dd$clustering$sequence, c(center, variant))
})

test_that("INDEL_RATE is validated", {
  on.exit(setDadaOpt(INDEL_RATE = 0))
  setDadaOpt(INDEL_RATE = 1.5)
  expect_error(dada(drp, err = err, verbose = 0), "INDEL_RATE")
})

test_that("estimateIndelRate returns a quantile of indel-only variant fractions", {
  setDadaOpt(INDEL_RATE = 0)
  dd <- dada(drp, err = err, verbose = 0)
  r <- suppressWarnings(estimateIndelRate(dd, drp, minAsvReads = 500))
  expect_equal(as.numeric(r), 100 / 1100)
  expect_equal(attr(r, "n_variants"), 1L)
})
