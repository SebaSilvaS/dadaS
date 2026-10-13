# Per-read quality model (PER_READ_QUALS). Off by default.
fnF <- system.file("extdata", "sam1F.fastq.gz", package = "dadaS")

test_that("the quality histogram matches the dereplicated reads", {
  d <- derepFastq(fnF, qualityHistogram = TRUE)
  h <- d$qhist
  expect_identical(h$index, which(d$uniques >= 2))
  expect_identical(dim(h$counts), c(length(h$index), ncol(d$quals), length(h$values)))
  # the average quality computed from the histogram equals derep$quals
  avg <- apply(h$counts, c(1, 2), function(cnt) sum(cnt * h$values)) / as.vector(d$uniques[h$index])
  lens <- nchar(names(d$uniques)[h$index])
  ok <- outer(lens, seq_len(ncol(d$quals)), ">=")
  expect_equal(avg[ok], d$quals[h$index, ][ok])
  # chunked reading gives the same histogram
  d2 <- derepFastq(fnF, n = 300, qualityHistogram = TRUE)
  expect_identical(d2$qhist$counts[order(names(d2$uniques)[d2$qhist$index]), , ],
                   h$counts[order(names(d$uniques)[h$index]), , ])
})

test_that("PER_READ_QUALS equals averaging when all reads of a unique share their qualities", {
  # Rewrite sam1F so that every read of a unique carries the qualities of its first read
  fq <- ShortRead::readFastq(fnF)
  sq <- as.character(ShortRead::sread(fq)); ql <- as.character(Biostrings::quality(Biostrings::quality(fq)))
  ql <- ql[match(sq, sq)]
  f2 <- tempfile(fileext = ".fastq.gz")
  ShortRead::writeFastq(ShortRead::ShortReadQ(ShortRead::sread(fq), Biostrings::BStringSet(ql), ShortRead::id(fq)), f2)
  d <- derepFastq(f2, qualityHistogram = TRUE)
  a <- dada(d, err = tperr1, verbose = 0)
  b <- dada(d, err = tperr1, verbose = 0, PER_READ_QUALS = TRUE)
  expect_identical(a$denoised, b$denoised)
  expect_identical(a$map, b$map)
  expect_identical(a$trans, b$trans)
  expect_equal(a$clustering$pval, b$clustering$pval)
})

test_that("PER_READ_QUALS runs end to end and counts transitions at the reads' own qualities", {
  d <- derepFastq(fnF, qualityHistogram = TRUE)
  b <- dada(d, err = tperr1, verbose = 0, PER_READ_QUALS = TRUE, OMEGA_C = 0)
  # with OMEGA_C = 0 every read is counted once at each aligned position (gap positions are not transitions)
  tot <- sum(d$uniques * nchar(names(d$uniques)))
  expect_lte(sum(b$trans), tot)
  expect_gt(sum(b$trans), 0.99 * tot)
  # only quality scores that occur in the reads have counts
  used <- as.integer(colnames(b$trans))[colSums(b$trans) > 0]
  expect_true(all(used %in% c(d$qhist$values, unique(as.vector(round(d$quals[d$uniques == 1, ]))))))
  expect_gt(nrow(b$clustering), 0)
  e <- learnErrors(fnF, verbose = 0, PER_READ_QUALS = TRUE)
  expect_equal(dim(getErrors(e)), dim(tperr1))
  expect_error(dada(derepFastq(fnF), err = tperr1, verbose = 0, PER_READ_QUALS = TRUE), "qualityHistogram")
})
