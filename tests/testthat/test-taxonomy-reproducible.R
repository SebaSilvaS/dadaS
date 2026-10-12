seqs <- getSequences(system.file("extdata", "example_seqs.fa", package = "dadaS"))
train <- system.file("extdata", "example_train_set.fa.gz", package = "dadaS")

test_that("assignTaxonomy bootstraps are reproducible with set.seed", {
  set.seed(42); a <- assignTaxonomy(seqs, train, outputBootstraps = TRUE)
  set.seed(42); b <- assignTaxonomy(seqs, train, outputBootstraps = TRUE)
  expect_identical(a, b)
})

test_that("assignTaxonomy results do not depend on the number of threads", {
  set.seed(7); a <- assignTaxonomy(seqs, train, outputBootstraps = TRUE, multithread = FALSE)
  set.seed(7); b <- assignTaxonomy(seqs, train, outputBootstraps = TRUE, multithread = 2)
  expect_identical(a, b)
})
