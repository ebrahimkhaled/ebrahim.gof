#' EDGE for Streaming Data: Calibration Monitoring Without Recomputation
#'
#' @description
#' Keeps the directed test of a frozen model up to date as new patients arrive, without
#' revisiting the records already seen. The external-mode statistic of
#' \code{\link{def.gof}} depends on the data only through four sums per risk group --
#' observed events, expected events, the binomial variance and the count -- so each new
#' record updates one group in constant time, and the test is recomputed from the
#' \eqn{G} group summaries alone.
#'
#' @details
#' \strong{The partition.} The groups are fixed in advance by cut points on the
#' predicted risk: either given as \code{breaks}, or the \eqn{G}-quantiles of a reference
#' set of predictions \code{p_ref} (for example the development data, or the first
#' batch). Because the cut points depend on predictions only and never on outcomes, the
#' statistic keeps its \eqn{\chi^2_{d+1}} reference under a calibrated model however the
#' risk distribution of later patients drifts: the groups need not stay of equal size.
#' Small groups weaken the protection against corrupted records that equal-frequency
#' groups give, so \code{summary()} reports the smallest group.
#'
#' \strong{Exactness.} After any sequence of updates the statistic equals the one-shot
#' external statistic computed on all records seen, with the same partition: the sums are
#' additive, so the order and the batching of the updates do not matter.
#'
#' \strong{Repeated looks.} One test at any time is valid. Testing after every batch and
#' acting on the first rejection is a sequential procedure, and the plain chi-squared
#' reference does not control its overall false-alarm rate; spend the level across looks
#' (for example, Bonferroni over a planned number of looks) or test at fixed times.
#'
#' @param p_ref Optional numeric vector of reference predicted probabilities; the cut
#'   points are its \code{G}-quantiles.
#' @param breaks Optional increasing numeric vector of interior cut points in (0, 1)
#'   (\code{G - 1} of them); used instead of \code{p_ref}.
#' @param G Number of risk groups (default \code{10}).
#' @param basis Calibration basis: \code{"poly3"} (default), \code{"poly2"},
#'   \code{"stukel"} or \code{"sym"}, as in \code{\link{def.gof}}.
#'
#' @return An object of class \code{edge_stream}. Add data with
#'   \code{update(object, y, p)}; read the test with \code{summary(object)}, a one-row
#'   \code{data.frame} in the format of \code{def.gof(..., external = TRUE)} with the
#'   number of records and the smallest group added.
#'
#' @seealso \code{\link{def.gof}}, \code{\link{run.all.external}}.
#'
#' @author Ebrahim Khaled Ebrahim \email{ebrahimkhaled@@alexu.edu.eg}
#'
#' @examples
#' set.seed(1)
#' p_dev <- plogis(rnorm(5000, -1.5, 1))          # predictions on the development data
#' s <- edge.stream(p_ref = p_dev, G = 10)
#' for (day in 1:30) {                            # a deployed model, 100 patients a day
#'   p <- plogis(rnorm(100, -1.5, 1))
#'   y <- rbinom(100, 1, p)
#'   s <- update(s, y, p)
#' }
#' summary(s)
#'
#' @export
edge.stream <- function(p_ref = NULL, breaks = NULL, G = 10, basis = "poly3") {
  basis <- match.arg(basis, c("poly3", "poly2", "stukel", "sym"))
  if (is.null(breaks)) {
    if (is.null(p_ref)) stop("Give either 'p_ref' (reference predictions) or 'breaks'.")
    p_ref <- as.numeric(p_ref)
    if (anyNA(p_ref) || any(p_ref <= 0 | p_ref >= 1)) stop("'p_ref' must lie strictly between 0 and 1.")
    breaks <- unique(as.numeric(stats::quantile(p_ref, seq_len(G - 1) / G, type = 7, names = FALSE)))
  }
  breaks <- sort(as.numeric(breaks))
  if (any(breaks <= 0 | breaks >= 1)) stop("'breaks' must lie strictly between 0 and 1.")
  G <- length(breaks) + 1L
  structure(list(breaks = breaks, G = G, basis = basis,
                 O = numeric(G), E = numeric(G), V = numeric(G), m = numeric(G)),
            class = "edge_stream")
}

#' @rdname edge.stream
#' @param object An \code{edge_stream} object.
#' @param y Binary (0/1) outcomes of the new records.
#' @param p Predicted probabilities of the new records, made without their outcomes.
#' @param ... Unused.
#' @exportS3Method stats::update
update.edge_stream <- function(object, y, p, ...) {
  y <- as.numeric(y); p <- as.numeric(p)
  if (length(y) != length(p)) stop("'y' and 'p' have different lengths.")
  if (!all(y %in% c(0, 1))) stop("'y' must be binary (0/1).")
  if (anyNA(p) || any(p <= 0 | p >= 1)) stop("'p' must lie strictly between 0 and 1.")
  if (!length(y)) return(object)
  p <- pmin(pmax(p, 1e-6), 1 - 1e-6)                          # the clamp of def.gof(), so the two agree exactly
  g <- findInterval(p, object$breaks, left.open = TRUE) + 1L   # group of each new record
  object$O <- object$O + tabulate_sum(g, y, object$G)
  object$E <- object$E + tabulate_sum(g, p, object$G)
  object$V <- object$V + tabulate_sum(g, p * (1 - p), object$G)
  object$m <- object$m + tabulate(g, object$G)
  object
}

# Internal: per-group sums of x, as a length-G vector (empty groups give 0).
tabulate_sum <- function(g, x, G) {
  out <- numeric(G)
  s <- rowsum(x, g)
  out[as.integer(rownames(s))] <- s[, 1]
  out
}

#' @rdname edge.stream
#' @exportS3Method base::summary
summary.edge_stream <- function(object, ...) {
  keep <- object$m > 0 & object$V > 0
  n <- sum(object$m)
  need <- c(poly3 = 4L, poly2 = 3L, stukel = 4L, sym = 2L)[[object$basis]]   # distinct group means the basis needs
  if (sum(keep) < need)
    return(data.frame(Test = "EDGE", Basis = object$basis, Test_Statistic = NA_real_,
                      df = NA_integer_, Method = "stream", p_value = NA_real_, n = n,
                      smallest_group = if (any(keep)) min(object$m[keep]) else 0, stringsAsFactors = FALSE))
  r    <- (object$O[keep] - object$E[keep]) / sqrt(object$V[keep])
  pbar <- object$E[keep] / object$m[keep]
  # the external-mode construction of def.gof(): a constant column joins the basis, Omega = I
  Z <- cbind(1, .def_basis(pbar, object$basis))
  Z <- Z[, colSums(abs(Z)) > 1e-8, drop = FALSE]
  Z <- Z / rep(sqrt(colSums(Z^2)), each = nrow(Z))
  Q <- qr(Z)
  Z <- Z[, Q$pivot[seq_len(Q$rank)], drop = FALSE]
  k <- ncol(Z)
  S <- sum(qr.fitted(qr(Z), r)^2)
  data.frame(Test = "EDGE", Basis = object$basis, Test_Statistic = S, df = k,
             Method = "stream", p_value = stats::pchisq(S, k, lower.tail = FALSE), n = n,
             smallest_group = min(object$m[keep]), stringsAsFactors = FALSE)
}

#' @rdname edge.stream
#' @param x An \code{edge_stream} object.
#' @exportS3Method base::print
print.edge_stream <- function(x, ...) {
  s <- summary(x)
  cat(sprintf("EDGE stream: %d records in %d groups (smallest %d), basis %s\n",
              as.integer(s$n), x$G, as.integer(s$smallest_group), x$basis))
  if (is.finite(s$Test_Statistic))
    cat(sprintf("  S = %.3f on %d df, p = %s\n", s$Test_Statistic, s$df, format.pval(s$p_value, digits = 3)))
  else cat("  not enough data yet\n")
  invisible(x)
}
