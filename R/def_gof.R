#' Directed Ebrahim-Farrington (DEF) Goodness-of-Fit Test
#'
#' @description
#' Performs the Directed Ebrahim-Farrington (DEF) goodness-of-fit test for a
#' fitted binary logistic regression model. DEF concentrates its power on a small
#' set of calibration-curve "shape" directions by projecting the grouped
#' standardized residuals onto a low-dimensional basis and testing the squared
#' length of that projection.
#'
#' \strong{Naming note:} this test is published under the name \strong{EDGE}
#' (Ebrahim Directed Goodness-of-fit Evaluation), and \code{\link{edge.gof}} is
#' the primary interface going forward. \code{def.gof()} is retained, unchanged,
#' as a fully supported legacy name.
#'
#' @details
#' The observations are sorted by predicted probability and split into \code{G}
#' equal-frequency groups; the standardized grouped residual vector \eqn{r} is
#' projected onto a basis matrix \eqn{Z} of smooth shapes, giving
#' \eqn{S = (Z'r)'(Z'Z)^{-1}(Z'r)}. Its null distribution is a weighted sum of
#' \eqn{\chi^2_1} variables with weights equal to the eigenvalues of
#' \eqn{(Z'Z)^{-1}Z'\Omega Z}, where \eqn{\Omega = I - U(X'WX)^{-1}U'} is the
#' estimation-adjusted covariance of the grouped residuals. The p-value uses a
#' Satterthwaite scaled-\eqn{\chi^2} approximation (default) or Imhof's method
#' (if the \pkg{CompQuadForm} package is installed). Bases: \code{"poly2"},
#' \code{"poly3"} (default), \code{"stukel"}, \code{"sym"}; \code{"ensemble"} runs
#' \code{"poly2"}, \code{"poly3"} and \code{"stukel"} and combines them via
#' \code{\link{def.ensemble.gof}}.
#'
#' With \code{weights = "score"} each column of \eqn{Z} is multiplied by
#' \eqn{\sqrt{V_g}}, the square root of its group's variance, so that \eqn{Z'r}
#' becomes \eqn{\sum_g z_g (O_g - E_g)}: for a logit fit, the score for adding the
#' grouped shape to the model as a step covariate. Its information after adjusting
#' for the fitted coefficients is \eqn{Z'\Omega Z}, and the statistic
#' \eqn{u'I^{-1}u} is referred to a \eqn{\chi^2} law on the rank of that
#' information (the number of columns unless one is redundant), read from it after
#' scaling to a correlation matrix. For a logit fit this is the Rao score test for
#' adding the grouped columns, and it agrees with \code{anova(..., test = "Rao")} up
#' to glm's convergence tolerance; for other links it is a score-type test.
#' Before either form is computed, a basis column shorter than \eqn{10^{-6}} times
#' the longest one is dropped. The unit form is the statistic as published;
#' the score form keeps a shape on the logit scale from losing its signal when the
#' group variances differ strongly, as they do at high discrimination.
#'
#' With fewer events (or fewer non-events) than groups, the grouped reference
#' distribution is unreliable. The p-value is still returned, with a warning;
#' a smaller \code{G} avoids it.
#'
#' @param object A fitted binary logistic \code{\link[stats]{glm}}, or a binary
#'   (0/1) response vector \code{y} (then supply \code{predicted_probs}).
#' @param predicted_probs Numeric predicted probabilities; required when
#'   \code{object} is a \code{y} vector, ignored when it is a glm.
#' @param X Optional design matrix, used only with the \code{y}/\code{predicted_probs}
#'   form: it enables the exact estimation-adjusted (\eqn{\Omega}) calibration
#'   (logit working weights assumed). Without it the conservative \eqn{\chi^2_k}
#'   reference is used and a warning is issued. Ignored when \code{object} is a glm.
#' @param G Integer number of equal-frequency groups (default 10; must be >= 3),
#'   or \code{"auto"} for \code{max(10, round(n / 25))}, the partition rule of the
#'   EDGE paper.
#' @param basis One of \code{"poly3"} (default), \code{"poly2"}, \code{"stukel"},
#'   \code{"sym"}, or \code{"ensemble"}. \code{"sym"} is one column,
#'   \eqn{\eta|\eta|} at the logit \eqn{\eta} of each group's mean fitted risk:
#'   Stukel's (1988) symmetric direction, aimed at tails that are too heavy or too
#'   light on both sides, for example a probit or cauchit truth fitted by a logit.
#' @param method One of \code{"satterthwaite"} (default) or \code{"imhof"}.
#'   Ignored when \code{weights = "score"}.
#' @param weights \code{"unit"} (default) is the statistic as published,
#'   \eqn{S = r'P_Z r} referred to a weighted chi-squared law. \code{"score"}
#'   multiplies each column by the square root of its group's variance, which for a
#'   logit fit makes the statistic the score test for adding the grouped shape to
#'   the model (a score-type test for other links). It is referred to chi-squared on
#'   the rank of its information matrix, which is the number of columns unless one
#'   is redundant (see Details).
#'
#' @return A one-row \code{data.frame} with columns \code{Test}, \code{Basis},
#'   \code{Test_Statistic} (the statistic \eqn{S}), \code{df}, \code{Method}, and
#'   \code{p_value}. For \code{weights = "score"}, \code{Method} is \code{"score"}
#'   and \code{df} is the integer rank the statistic is referred to. When
#'   \code{basis = "ensemble"}, the return is that of \code{\link{def.ensemble.gof}}.
#'
#' @references
#' Ebrahim EK, El-Kotory A (2026). "A Directional Hosmer-Lemeshow Goodness-of-Fit
#' Test for Sparse Logistic Regression." arXiv:2607.15454 [stat.ME].
#' \doi{10.48550/arXiv.2607.15454}
#'
#' Ebrahim EK, El-Kotory A (2026). "EDGE: A Closed-Form Directed Goodness-of-Fit
#' Test for Sparse Logistic Regression." arXiv:2608.20511 [stat.ME].
#' \doi{10.48550/arXiv.2608.20511}
#'
#' @author Ebrahim Khaled Ebrahim \email{ebrahimkhaled@@alexu.edu.eg}
#'
#' @examples
#' ## gof_demo carries a documented smooth calibration misfit: the risk bends in age,
#' ## and a model linear in age misses it. The point of a directed test is to see that.
#' data("gof_demo", package = "ebrahim.gof")
#' wrong <- glm(outcome ~ age + bmi + sex + treatment,
#'              data = gof_demo, family = binomial())
#' def.gof(wrong)                       # default poly3 basis
#' def.gof(wrong, basis = "stukel")     # tail-shape basis
#' def.gof(wrong, basis = "sym")        # symmetric tail direction, one column
#' def.gof(wrong, weights = "score")    # score form of the poly3 basis
#' def.gof(wrong, basis = "ensemble")   # combine poly2, poly3 and stukel (CCT)
#'
#' ## give the model the term it was missing, and the same test stands down
#' right <- glm(outcome ~ poly(age, 2) + bmi + sex + treatment,
#'              data = gof_demo, family = binomial())
#' def.gof(right)
#'
#' @seealso \code{\link{ef.gof}}, \code{\link{def.ensemble.gof}}.
#' @importFrom stats fitted predict model.matrix qlogis poly pchisq
#' @concept goodness-of-fit
#' @concept calibration
#' @concept logistic regression
#' @concept model diagnostics
#' @concept DEF
#' @concept directed test
#' @export
def.gof <- function(object, predicted_probs = NULL, X = NULL, G = 10,
                    basis   = c("poly3", "poly2", "stukel", "sym", "ensemble"),
                    method  = c("satterthwaite", "imhof"),
                    weights = c("unit", "score")) {

  basis   <- match.arg(basis)
  method  <- match.arg(method)
  weights <- match.arg(weights)
  if (!identical(G, "auto") && (!is.numeric(G) || length(G) != 1 || G < 3)) {
    stop("'G' must be a single integer >= 3, or 'auto'.")
  }

  if (basis == "ensemble")
    return(def.ensemble.gof(object, predicted_probs = predicted_probs, X = X, G = G,
                            weights = weights))

  # --- accept either a fitted glm, OR (y, predicted_probs[, X]) ---
  if (inherits(object, "glm")) {
    if (object$family$family != "binomial")
      stop("'object' must be a binomial glm (or pass y, predicted_probs, X).")
    y   <- as.numeric(object$y)
    ph  <- pmin(pmax(as.numeric(stats::fitted(object)), 1e-6), 1 - 1e-6)
    eta <- as.numeric(stats::predict(object, type = "link"))
    dmu <- object$family$mu.eta(eta)
    X   <- stats::model.matrix(object)
    naive <- FALSE
  } else {
    if (!is.numeric(object))
      stop("'object' must be a fitted binomial glm or a numeric (0/1) y vector.")
    y <- as.numeric(object)
    if (is.null(predicted_probs))
      stop("Provide 'predicted_probs' when 'object' is not a glm.")
    ph  <- pmin(pmax(as.numeric(predicted_probs), 1e-6), 1 - 1e-6)
    dmu <- ph * (1 - ph)
    naive <- is.null(X)
    if (naive)
      warning("def.gof: no model/design matrix supplied; using the conservative ",
              "chi-square reference (Omega = I). Pass the fitted glm or X for the ",
              "exact estimation-adjusted test.")
  }
  n <- length(y)
  if (!all(y %in% c(0, 1))) stop("DEF needs a binary (0/1) response.")
  if (length(ph) != n) stop("'object' (y) and 'predicted_probs' lengths differ.")
  if (identical(G, "auto")) G <- .def_auto_G(n)
  if (G > n) stop("'G' cannot exceed the number of observations.")
  if (min(sum(y), n - sum(y)) < G) .def_warn_few_events(sum(y), n, G)

  V <- ph * (1 - ph)
  w <- dmu^2 / V

  # --- equal-frequency groups by fitted probability ---
  grp  <- pmin(ceiling(rank(ph, ties.method = "first") / (n / G)), G)
  idx  <- split(seq_len(n), grp)
  Gn   <- length(idx)
  og   <- vapply(idx, function(I) sum(y[I]),   numeric(1))
  eg   <- vapply(idx, function(I) sum(ph[I]),  numeric(1))
  Vg   <- vapply(idx, function(I) sum(V[I]),   numeric(1))
  pbar <- vapply(idx, function(I) mean(ph[I]), numeric(1))
  r    <- (og - eg) / sqrt(Vg)

  # --- estimation-adjusted covariance ---
  if (naive) {
    Omega <- diag(Gn)
  } else {
    U     <- t(vapply(idx, function(I) colSums(dmu[I] * X[I, , drop = FALSE]),
                      numeric(ncol(X)))) / sqrt(Vg)
    Omega <- diag(Gn) - U %*% solve(crossprod(X, w * X)) %*% t(U)
  }

  # --- shape basis Z ---
  Z <- .def_basis(pbar, basis)
  # A column is dropped when it is negligible next to the largest one (a Stukel half reaching
  # a single group whose mean risk is a hair above 0.5), and the rest are scaled to unit
  # length. Both statistics and the eigenvalues do not depend on the column scale; solve() does.
  nz <- sqrt(colSums(Z^2))
  ok <- nz > 1e-6 * max(nz)
  Z  <- Z[, ok, drop = FALSE] / rep(nz[ok], each = nrow(Z))
  if (ncol(Z) < 1)
    stop("The chosen basis is degenerate for this fit. Try basis = 'poly3' or a larger G.")

  if (weights == "score") {
    # Score form: each column times sqrt(V_g), so Z'r = sum_g z_g (O_g - E_g), the score for
    # adding the group-level step covariate to the model. Its information after adjusting for
    # the fitted coefficients is Z' Omega Z, and u'I^-1 u is chi-square on its rank. The rank is
    # read from I scaled to a correlation matrix, so the scale of a column does not decide it.
    Zs <- Z * sqrt(Vg)
    u  <- drop(crossprod(Zs, r))
    I  <- crossprod(Zs, Omega %*% Zs)
    d  <- sqrt(pmax(diag(I), 0))
    ok <- d > 0
    k  <- 0L
    if (any(ok)) {
      R   <- I[ok, ok, drop = FALSE] / outer(d[ok], d[ok])
      ev  <- eigen((R + t(R)) / 2, symmetric = TRUE)
      pos <- ev$values > 1e-8
      k   <- sum(pos)
      S   <- sum(drop(crossprod(ev$vectors[, pos, drop = FALSE], u[ok] / d[ok]))^2 / ev$values[pos])
    }
    return(data.frame(Test = "Directed Ebrahim-Farrington", Basis = basis,
                      Test_Statistic = if (k > 0L) S else NA_real_,
                      df = if (k > 0L) k else NA_real_, Method = "score",
                      p_value = if (k > 0L) stats::pchisq(S, k, lower.tail = FALSE) else NA_real_,
                      stringsAsFactors = FALSE))
  }

  ZtZ <- crossprod(Z)
  Zr  <- crossprod(Z, r)
  S   <- as.numeric(t(Zr) %*% solve(ZtZ) %*% Zr)

  lam <- Re(eigen(solve(ZtZ) %*% (t(Z) %*% Omega %*% Z), only.values = TRUE)$values)
  lam <- lam[lam > 1e-9]
  if (length(lam) == 0) {
    return(data.frame(Test = "Directed Ebrahim-Farrington", Basis = basis,
                      Test_Statistic = S, df = NA_real_, Method = method,
                      p_value = NA_real_, stringsAsFactors = FALSE))
  }

  pval <- .def_pvalue(S, lam, method)
  nu   <- sum(lam)^2 / sum(lam^2)

  data.frame(Test = "Directed Ebrahim-Farrington", Basis = basis,
             Test_Statistic = S, df = nu, Method = method, p_value = pval,
             stringsAsFactors = FALSE)
}

# Internal: build the shape-basis matrix Z from the group mean probabilities.
.def_basis <- function(pbar, basis) {
  if (basis %in% c("poly2", "poly3")) {
    deg <- if (basis == "poly2") 2L else 3L
    if (length(unique(round(pbar, 8))) < deg + 1)
      stop("Too few distinct group probabilities for basis = '", basis,
           "'. Use a smaller polynomial basis or a larger sample/G.")
    as.matrix(stats::poly(pbar, deg))
  } else {
    e <- stats::qlogis(pbar)
    if (basis == "sym") return(cbind(e * abs(e)))   # Stukel's symmetric direction (alpha1 = alpha2)
    cbind(e, e^2 * (e >= 0), -e^2 * (e < 0))
  }
}

# Internal: the number of groups for G = "auto", the partition rule of the EDGE paper.
.def_auto_G <- function(n) max(10, round(n / 25))

# Internal: warn that the grouped reference is unreliable with fewer events (or non-events)
# than groups. The warning has its own class, so the battery can put it in Note and the
# ensemble can raise it once rather than once per basis.
.def_warn_few_events <- function(ne, n, G) {
  what <- if (ne <= n - ne) "events" else "non-events"
  msg  <- sprintf(paste("def.gof: %d %s for G = %s groups; the grouped reference distribution",
                        "is unreliable with fewer %s than groups."),
                  as.integer(min(ne, n - ne)), what, format(G), what)
  warning(structure(class = c("def_few_events", "warning", "condition"),
                    list(message = msg, call = NULL)))
}

# Internal: p-value of S under sum_j lambda_j chi^2_1.
.def_pvalue <- function(S, lam, method) {
  if (method == "imhof") {
    if (!requireNamespace("CompQuadForm", quietly = TRUE))
      stop("method = 'imhof' requires the 'CompQuadForm' package. ",
           "Install it, or use method = 'satterthwaite'.")
    p <- CompQuadForm::imhof(S, lam)$Qq
    return(min(max(p, 0), 1))
  }
  cc <- sum(lam^2) / sum(lam)
  nu <- sum(lam)^2 / sum(lam^2)
  stats::pchisq(S / cc, df = nu, lower.tail = FALSE)
}
