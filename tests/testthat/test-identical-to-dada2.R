# With default options, DadaS must give exactly the same results as DADA2 (the reference).
skip_if_not_installed("dada2")

fnF <- system.file("extdata", c("sam1F.fastq.gz", "sam2F.fastq.gz"), package = "dadaS")
fnR <- system.file("extdata", c("sam1R.fastq.gz", "sam2R.fastq.gz"), package = "dadaS")

test_that("filterAndTrim, derepFastq and learnErrors match DADA2", {
  fF <- tempfile(fileext = c(".1.fq.gz", ".2.fq.gz")); fR <- tempfile(fileext = c(".1.fq.gz", ".2.fq.gz"))
  fF2 <- tempfile(fileext = c(".1.fq.gz", ".2.fq.gz")); fR2 <- tempfile(fileext = c(".1.fq.gz", ".2.fq.gz"))
  a <- dadaS::filterAndTrim(fnF, fF, fnR, fR, truncLen = c(240, 160), maxEE = 2)
  b <- dada2::filterAndTrim(fnF, fF2, fnR, fR2, truncLen = c(240, 160), maxEE = 2)
  expect_identical(unname(a), unname(b))
  expect_identical(dadaS::derepFastq(fF[1])$uniques, dada2::derepFastq(fF2[1])$uniques)
  set.seed(1); eS <- dadaS::learnErrors(fF, multithread = 2, verbose = 0)
  set.seed(1); e2 <- dada2::learnErrors(fF2, multithread = 2, verbose = 0)
  expect_identical(dadaS::getErrors(eS), dada2::getErrors(e2))
})

test_that("dada, mergePairs and removeBimeraDenovo match DADA2 (serial and parallel)", {
  drpF <- dadaS::derepFastq(fnF); drpR <- dadaS::derepFastq(fnR)
  for (mt in c(FALSE, 2)) {
    dS <- dadaS::dada(drpF, err = dadaS::tperr1, multithread = mt, verbose = 0)
    d2 <- dada2::dada(drpF, err = dada2::tperr1, multithread = mt, verbose = 0)
    for (i in seq_along(dS)) {
      expect_identical(dS[[i]]$denoised, d2[[i]]$denoised)
      expect_identical(dS[[i]]$map, d2[[i]]$map)
      expect_equal(dS[[i]]$clustering, d2[[i]]$clustering)
    }
  }
  dSF <- dadaS::dada(drpF, err = dadaS::tperr1, verbose = 0); dSR <- dadaS::dada(drpR, err = dadaS::tperr1, verbose = 0)
  mS <- dadaS::mergePairs(dSF, drpF, dSR, drpR)
  m2 <- dada2::mergePairs(dSF, drpF, dSR, drpR)
  expect_equal(mS, m2)
  stS <- dadaS::makeSequenceTable(mS)
  expect_identical(dadaS::removeBimeraDenovo(stS, multithread = 2), dada2::removeBimeraDenovo(stS, multithread = 2))
})
