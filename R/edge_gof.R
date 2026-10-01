#' EDGE: Directed Goodness-of-Fit Test for Binary Logistic Regression
#'
#' @description
#' \code{edge.gof()} is the primary interface to the EDGE test (Efficient Directed
#' Grouped Examination): a grouped, \emph{directed} goodness-of-fit test
#' for binary logistic regression under sparse data. EDGE projects the grouped
#' standardized residuals onto a small pre-specified basis of calibration shapes
#' (cubic \code{"poly3"} by default) and refers the resulting quadratic form to
#' its closed-form weighted chi-squared null distribution -- no refit, no
#' resampling, no tuning.
#'
#' \code{edge.gof()} computes exactly the same statistic as the legacy name
#' \code{\link{def.gof}} (retained for backward compatibility); the returned
#' \code{Test} label is \code{"EDGE"}.
#'
#' @details
#' \strong{Two partitions, side by side.} By default the test is reported at two
#' partitions, as the EDGE paper recommends. The \emph{default partition},
#' \code{G = "auto"} (\eqn{\max(10, \lceil n/25 \rceil)} groups), refines with the sample and has
#' more power, above all against misfit at the extremes of risk. \emph{Ten groups},
#' \code{G = 10}, keep many records in the extreme groups and so tolerate more
#' corrupted records, at some cost in power. Choose which one decides in the analysis
#' plan, from what is known about how the data were collected, not after seeing the
#' results. The first row is marked \code{Role = "verdict"} and the second
#' \code{Role = "check"}: with the default \code{G} the default partition decides and
#' ten groups are the robustness check; to let ten groups decide, give
#' \code{G = c(10, "auto")}. If only the default partition rejects, the signal sits in the
#' extreme groups: check those records. Give a single \code{G} to get one row.
#'
#' @inheritParams def.gof
#' @param G Number of groups: \code{"auto"}, a number, or a vector of them; the
#'   default \code{c("auto", 10)} reports both partitions.
#' @param y Optional alias for \code{object} when frozen predictions are tested:
#'   \code{edge.gof(y = y, predicted_probs = p, external = TRUE)}.
#'
#' @return A \code{data.frame} with one row per partition and columns \code{Test}
#'   (\code{"EDGE"}), \code{Partition}, \code{Role} (\code{"verdict"} for the first
#'   row, \code{"check"} for the second), \code{Basis}, \code{Test_Statistic},
#'   \code{df}, \code{Method} and \code{p_value}, as documented in \code{\link{def.gof}}.
#'
#' @references
#' Ebrahim EK, El-Kotory A (2026). "EDGE: A Closed-Form Directed Goodness-of-Fit
#' Test for Sparse Logistic Regression." arXiv:2608.20511 [stat.ME].
#' \doi{10.48550/arXiv.2608.20511}
#'
#' Ebrahim EK, El-Kotory A (2026). "A Directional Hosmer-Lemeshow Goodness-of-Fit
#' Test for Sparse Logistic Regression." arXiv:2607.15454 [stat.ME].
#' \doi{10.48550/arXiv.2607.15454}
#'
#' @author Ebrahim Khaled Ebrahim \email{ebrahimkhaled@@alexu.edu.eg}
#'
#' @examples
#' set.seed(1)
#' x <- runif(500, -3, 3)
#' y <- rbinom(500, 1, plogis(0.6 * x))
#' fit <- glm(y ~ x, family = binomial())
#' edge.gof(fit)                      # cubic basis, at the default partition and at ten groups
#' edge.gof(fit, G = 10)              # one partition only
#' edge.gof(fit, basis = "stukel")    # Stukel-shape basis
#' edge.gof(fit, basis = "sym")       # Stukel's symmetric direction, one column
#' edge.gof(fit, basis = "sym", weights = "score")   # its score form
#' edge.gof(fit, G = "auto")          # the default partition alone: max(10, ceiling(n / 25)) = 20 here
#'
#' @seealso \code{\link{def.gof}} (legacy name), \code{\link{ef.gof}},
#'   \code{\link{def.ensemble.gof}}, \code{\link{run.all.gof}}.
#' @concept goodness-of-fit
#' @concept calibration
#' @concept logistic regression
#' @concept model diagnostics
#' @concept EDGE
#' @concept directed test
#' @concept sparse data
#' @export
edge.gof <- function(object, predicted_probs = NULL, X = NULL, G = c("auto", 10),
                     basis = "poly3", method = "satterthwaite", weights = "unit",
                     external = FALSE, y = NULL) {
  if (missing(object)) {
    if (is.null(y)) stop("supply a fitted glm or the outcome vector")
    object <- y
  }
  n <- if (inherits(object, "glm")) length(object$y) else length(object)
  one <- function(g) {
    out <- def.gof(object, predicted_probs = predicted_probs, X = X, G = g,
                   basis = basis, method = method, weights = weights, external = external)
    out$Test <- "EDGE"
    gn <- if (identical(g, "auto")) .def_auto_G(n) else as.numeric(g)
    out$Partition <- if (identical(g, "auto")) sprintf("default (G = %d)", as.integer(gn)) else
      if (gn == 10) "ten groups (G = 10)" else sprintf("G = %d", as.integer(gn))
    out[, c("Test", "Partition", setdiff(names(out), c("Test", "Partition")))]
  }
  Gs <- lapply(as.list(G), function(g) if (identical(g, "auto")) g else as.numeric(g))
  out <- do.call(rbind, lapply(Gs, one))
  ## when n <= 250 the default partition is itself ten groups: report it once
  if (nrow(out) == 2L && isTRUE(all.equal(out$Test_Statistic[1], out$Test_Statistic[2])))
    out <- out[1, , drop = FALSE]
  out$Role <- c("verdict", "check", rep("check", max(0L, nrow(out) - 2L)))[seq_len(nrow(out))]
  out <- out[, c("Test", "Partition", "Role", setdiff(names(out), c("Test", "Partition", "Role")))]
  rownames(out) <- NULL
  out
}
