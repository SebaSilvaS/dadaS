#' Variants that differ from their ASV only by indels (DadaS)
#'
#' For each ASV in a \code{dada()} result, finds the unique sequences assigned to it that differ from
#' the ASV center only by internal indels (zero substitutions). With the default
#' \code{INDEL_RATE = 0}, DADA2 never splits such variants into their own ASV; this table shows how
#' many there are and how abundant they are.
#'
#' Alignments use the same scoring options as \code{dada()} (\code{MATCH}, \code{MISMATCH},
#' \code{GAP_PENALTY}, \code{BAND_SIZE}), ends-free, so terminal gaps are ignored. Consecutive gap
#' columns count as one indel event, matching how \code{INDEL_RATE} counts them.
#'
#' @param dadaObj A \code{dada-class} object, or a list of them (one per sample).
#' @param derep The matching \code{derep-class} object(s), same length and order as \code{dadaObj}.
#' @param minAbundance Minimum abundance of a variant to be reported (default 2).
#' @return A data.frame with columns \code{sample}, \code{asv} (index in \code{$clustering}),
#'   \code{asv_reads} (total reads of the ASV, including the variant), \code{variant_reads},
#'   \code{n_indels} and \code{fraction} (\code{variant_reads / asv_reads}). Empty data.frame if no variants are found.
#' @seealso \code{\link{estimateIndelRate}}, \code{\link{setDadaOpt}} (\code{INDEL_RATE})
#' @export
#' @examples
#' derep1 <- derepFastq(system.file("extdata", "sam1F.fastq.gz", package="dadaS"))
#' dada1 <- dada(derep1, err=tperr1, errorEstimationFunction=loessErrfun, selfConsist=FALSE)
#' indelOnlyVariants(dada1, derep1)
indelOnlyVariants <- function(dadaObj, derep, minAbundance = 2) {
  if (is(dadaObj, "dada")) dadaObj <- list(dadaObj)
  if (is(derep, "derep")) derep <- list(derep)
  if (length(dadaObj) != length(derep)) stop("dadaObj and derep must have the same length.")
  nm <- names(dadaObj); if (is.null(nm)) nm <- as.character(seq_along(dadaObj))
  # Number of substitutions and indel events between two aligned sequences, ignoring terminal gaps
  count_diffs <- function(a0, a1) {
    x <- strsplit(a0, "")[[1]]; y <- strsplit(a1, "")[[1]]
    inner <- which(x != "-" & y != "-"); if (!length(inner)) return(c(NA, NA))
    r <- range(inner); x <- x[r[1]:r[2]]; y <- y[r[1]:r[2]]
    gap <- x == "-" | y == "-"
    c(sum(!gap & x != y), sum(diff(c(FALSE, gap)) == 1))
  }
  o <- getDadaOpt(c("MATCH", "MISMATCH", "GAP_PENALTY", "BAND_SIZE"))
  rows <- lapply(seq_along(dadaObj), function(i) {
    res <- dadaObj[[i]]; der <- getDerep(derep[[i]])
    uniq <- names(der$uniques); ab <- unname(der$uniques); asv <- res$map
    centers <- res$clustering$sequence
    cand <- which(!is.na(asv) & ab >= minAbundance & uniq != centers[asv])
    if (!length(cand)) return(NULL)
    al <- matrix(nwalign(centers[asv[cand]], uniq[cand], match = o$MATCH, mismatch = o$MISMATCH,
                         gap = o$GAP_PENALTY, band = o$BAND_SIZE, endsfree = TRUE, vec = TRUE),
                 ncol = 2, byrow = TRUE)
    cnt <- matrix(unlist(mapply(count_diffs, al[, 1], al[, 2], SIMPLIFY = FALSE, USE.NAMES = FALSE)),
                  ncol = 2, byrow = TRUE)
    k <- which(cnt[, 1] == 0 & cnt[, 2] >= 1)
    if (!length(k)) return(NULL)
    data.frame(sample = nm[i], asv = asv[cand[k]],
               asv_reads = res$clustering$abundance[asv[cand[k]]],
               variant_reads = ab[cand[k]], n_indels = cnt[k, 2], stringsAsFactors = FALSE)
  })
  tab <- do.call(rbind, rows)
  if (is.null(tab)) return(data.frame())
  tab$fraction <- tab$variant_reads / tab$asv_reads
  tab
}

#' Estimate the per-event indel error rate from the data (DadaS)
#'
#' Uses the most abundant ASVs, where errors are well sampled: for each variant that differs from its
#' center by exactly one indel, computes its abundance relative to the ASV, and returns a high
#' quantile of that distribution as a conservative value for \code{setDadaOpt(INDEL_RATE = ...)}.
#' Most of these variants are errors (polymerase slippage in homopolymers); a high quantile avoids
#' turning errors into ASVs, at the cost of not resolving very rare true indel variants.
#'
#' Run it on a \code{dada()} result obtained with \code{INDEL_RATE = 0}, separately for forward and
#' reverse reads: their indel error profiles differ, so the same quantile usually gives different rates.
#'
#' @param dadaObj,derep As in \code{\link{indelOnlyVariants}}.
#' @param minAsvReads Only ASVs with at least this many reads are used (default 2000).
#' @param quantile Quantile of the relative abundances that is returned (default 0.99).
#' @return The estimated rate (numeric), with attributes \code{n_variants} and \code{quantiles}.
#' @seealso \code{\link{indelOnlyVariants}}, \code{\link{setDadaOpt}}
#' @export
estimateIndelRate <- function(dadaObj, derep, minAsvReads = 2000, quantile = 0.99) {
  tab <- indelOnlyVariants(dadaObj, derep, minAbundance = 2)
  if (!nrow(tab)) stop("No indel-only variants found; cannot estimate a rate.")
  x <- tab[tab$asv_reads >= minAsvReads & tab$n_indels == 1, ]
  if (nrow(x) < 20) warning("Only ", nrow(x), " variants available; the estimate is unreliable.")
  if (!nrow(x)) return(NA_real_)
  q <- stats::quantile(x$fraction, unique(c(0.5, 0.9, 0.99, quantile)), names = TRUE)
  r <- unname(stats::quantile(x$fraction, quantile))
  attr(r, "n_variants") <- nrow(x); attr(r, "quantiles") <- q
  r
}
