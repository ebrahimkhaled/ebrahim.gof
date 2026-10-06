#' Ebrahim-Farrington Goodness-of-Fit Test for Logistic Regression
#'
#' @description
#' Performs the Ebrahim-Farrington (EF) goodness-of-fit test for logistic regression
#' models: Farrington's (1996) correction applied to the Hosmer-Lemeshow risk groups.
#' It keeps the protection of the risk-ordered table against gross recording errors and
#' is the score test for records of similar fitted risk sharing a miscalibration.
#'
#' @details
#' The records are sorted by fitted risk and cut into \code{G} groups of sizes that differ
#' by at most one. For group \eqn{g} with \eqn{n_g} records, \eqn{o_g} observed and
#' \eqn{e_g = n_g \bar\pi_g} expected events and \eqn{V_g = n_g\bar\pi_g(1-\bar\pi_g)},
#' the standardized residual is \eqn{r_g = (o_g - e_g)/\sqrt{V_g}} and
#' \eqn{HL = \sum_g r_g^2}. With \eqn{b_g = (1 - 2\bar\pi_g)/\sqrt{V_g}},
#' \deqn{EF = HL - C, \qquad C = \sum_g b_g r_g .}
#' \eqn{EF} keeps only the within-group products of the residuals: \eqn{HL = G + C + Q}
#' and \eqn{EF = G + Q} (Ebrahim, Khattab and El-Kotory, 2026).
#'
#' \strong{Reference distributions.} \code{reference = "chisq"} (the default) refers
#' \eqn{EF} to \eqn{\chi^2_{G-2}}, the reference of the Hosmer-Lemeshow test, which the
#' fixed-\eqn{G} theory justifies. \code{reference = "normal"} is the plug-in normal
#' reference for many groups (Theorem 3 of the paper): \eqn{z = (EF - \hat\mu)/\hat\sigma}
#' with model-based moments evaluated at the fit, referred to the upper tail of the
#' standard normal. It needs the model's design matrix: pass the fitted \code{glm}, or
#' \code{X}.
#'
#' \strong{Which reference, how many groups.}
#' \itemize{
#'   \item With \code{G = 10} and \eqn{\chi^2_{G-2}}, each group needs about one expected
#'     event and one expected non-event. When some group expects fewer than one half,
#'     the test is conservative and the Hosmer-Lemeshow reference is unreliable too: use
#'     fewer groups. A warning of class \code{ef_sparse_groups} says so.
#'   \item With many groups, as under the rule of Paul, Pennell and Lemeshow (2013)
#'     (\code{G = "paul"}), use \code{reference = "normal"}, not \eqn{\chi^2_{G-2}}, and
#'     prefer groups of about 25 records (\code{G = "auto"}, that is
#'     \code{max(10, ceiling(n / 25))}): they keep the protection against gross errors and
#'     do not dilute a smooth misfit over very many small groups.
#' }
#'
#' \strong{Ties.} Records with tied fitted risks must not be ordered by the response: a
#' group boundary inside a covariate pattern then manufactures misfit. Ties are broken at
#' random, so a result with tied fitted risks depends on the random number stream; call
#' \code{set.seed()} first to reproduce it. Nothing is drawn when no fitted risks are tied.
#' Alternatively, \code{groups = "pattern"} keeps each pattern (records with the same
#' fitted risk) in one group of its own, as for binomial data expanded to binary records,
#' or \code{groups} takes any vector of group labels.
#'
#' \strong{Directional read-out.} The returned \code{C} and the per-group terms
#' \eqn{c_g = b_g r_g} in \code{attr(, "groups")} describe the direction of the misfit.
#' Since \eqn{b_g} changes sign at \eqn{\bar\pi_g = 1/2}, \eqn{C < 0} says that, on
#' balance, the model underestimates how fast the risk approaches one in the high-risk
#' groups or approaches zero in the low-risk groups. The read-out is descriptive.
#'
#' For grouped binomial data (\code{G = NULL} with \code{m} and \code{model}), the
#' original Farrington test is applied with its full variance calculations.
#'
#' @param y A fitted binary logistic \code{glm} (then \code{predicted_probs} and \code{X}
#'   are taken from it), or a numeric vector of binary responses (0/1) for binary data /
#'   counts of successes for grouped data.
#' @param predicted_probs Numeric vector of predicted probabilities from the
#'   logistic regression model. Must be same length as \code{y}.
#' @param model Optional \code{glm} object. Required only for the original
#'   Farrington test with grouped data (when \code{m} is provided and \code{G} is NULL).
#' @param m Optional numeric vector of trial counts for each observation
#'   (for grouped data). If NULL, data is assumed to be binary.
#' @param G Number of risk groups: an integer (default 10); \code{"auto"} for groups of
#'   about 25 records, \code{max(10, ceiling(n / 25))}; or \code{"paul"} for the rule of
#'   Paul, Pennell and Lemeshow (2013), \eqn{\max\{10, \min(n_1/2, (n-n_1)/2,
#'   2 + 8(n/1000)^2)\}} rounded down, \eqn{n_1} the number of events (designed for
#'   \eqn{1000 < n \le 25000}). Ignored when \code{groups} is given. If NULL (and no
#'   \code{groups}), no grouping is performed and \code{m} must be provided.
#' @param method Deprecated. \code{"chisq"} is \code{reference = "chisq"};
#'   \code{"normal"} is the standardization \eqn{(EF - (G-2))/\sqrt{2(G-2)}} of package
#'   versions up to 2.9.0, kept, with a warning, to reproduce earlier results. It is not
#'   the normal reference of the EF paper, which is \code{reference = "normal"}.
#' @param reference \code{"chisq"} (default) or \code{"normal"}; see Details.
#' @param X Design matrix of the fitted model (with the intercept column), needed for
#'   \code{reference = "normal"} when \code{y} is not a fitted \code{glm}.
#' @param groups \code{NULL} (default) for equal-size groups by rank of the fitted risk,
#'   ties at random; \code{"pattern"} for one group per distinct fitted risk; or a vector
#'   of group labels, one per record.
#'
#' @return
#' For the grouped test, a one-row data frame with columns
#' \item{Test}{\code{"Ebrahim-Farrington"}}
#' \item{Test_Statistic}{\eqn{EF} for \code{reference = "chisq"}; the standardized
#'   \eqn{z} for \code{reference = "normal"} (and for the deprecated \code{method = "normal"})}
#' \item{df}{\eqn{G - 2} for the chi-squared reference, \code{NA} otherwise}
#' \item{Reference}{\code{"chisq"}, \code{"normal"} or \code{"thesis"} (the deprecated
#'   \code{method = "normal"})}
#' \item{p_value}{the p-value}
#' \item{G}{the number of groups used}
#' \item{HL}{the Hosmer-Lemeshow statistic on the same groups}
#' \item{C}{the removed directional part, \eqn{C = HL - EF}}
#' and an attribute \code{"groups"}: one row per group with \code{group}, \code{n},
#' \code{observed}, \code{expected}, \code{pbar}, \code{r} (\eqn{r_g}), \code{b}
#' (\eqn{b_g}) and \code{c} (\eqn{c_g = b_g r_g}, summing to \code{C}).
#' For the original Farrington test (\code{G = NULL}), a data frame with \code{Test},
#' \code{Test_Statistic} (standardized) and \code{p_value}.
#'
#' @note
#' \itemize{
#'   \item For grouped data (\code{m} provided \emph{and} \code{G} set to \code{NULL}):
#'         Use the original Farrington test, which requires the fitted model object.
#'         Leaving \code{G} at its default keeps the automatic grouping, and \code{m}
#'         and \code{model} are then ignored.
#'   \item For binary data with \code{m=1} for all observations and no grouping,
#'         the test is not applicable and will return a p-value of 1.
#' }
#'
#' @references
#' Ebrahim EK, Khattab IG, El-Kotory A (2026). "A Modified Hosmer-Lemeshow
#' Goodness-of-Fit Test for Asymmetric Links: Second-Order Power and Robustness."
#' arXiv:2607.15454 [stat.ME]. \doi{10.48550/arXiv.2607.15454}
#'
#' Farrington CP (1996). "On Assessing Goodness of Fit of Generalized Linear Models
#' to Sparse Data." \emph{Journal of the Royal Statistical Society, Series B},
#' \bold{58}(2), 349-360. \doi{10.1111/j.2517-6161.1996.tb02086.x}
#'
#' Hosmer DW, Lemeshow S (1980). "A goodness-of-fit test for the multiple logistic
#' regression model." \emph{Communications in Statistics - Theory and Methods},
#' \bold{9}(10), 1043-1069. \doi{10.1080/03610928008827941}
#'
#' Paul P, Pennell ML, Lemeshow S (2013). "Standardizing the power of the
#' Hosmer-Lemeshow goodness of fit test in large data sets." \emph{Statistics in
#' Medicine}, \bold{32}(1), 67-80. \doi{10.1002/sim.5525}
#'
#' Ebrahim EK (2026). "Goodness-of-Fit Tests and Calibration Machine-Learning
#' Algorithms for Logistic Regression with Sparse Data." M.Sc. thesis, Alexandria
#' University. arXiv:2608.11140 [stat.ME]. \doi{10.48550/arXiv.2608.11140}
#'
#' @author Ebrahim Khaled Ebrahim \email{ebrahimkhaled@@alexu.edu.eg}
#'
#' @examples
#' # Example 1: Binary data with automatic grouping (Ebrahim-Farrington test)
#' set.seed(123)
#' n <- 500
#' x <- rnorm(n)
#' linpred <- 0.5 + 1.2 * x
#' prob <- 1 / (1 + exp(-linpred))
#' y <- rbinom(n, 1, prob)
#'
#' # Fit logistic regression
#' model <- glm(y ~ x, family = binomial())
#'
#' # Ten groups, chi-square(G - 2) reference
#' result <- ef.gof(model)
#' print(result)
#' attr(result, "groups")              # r_g, b_g and the read-out c_g
#'
#' # Groups of about 25 records with the normal reference
#' ef.gof(model, G = "auto", reference = "normal")
#'
#' # Frozen predictions: y and the predicted probabilities
#' ef.gof(y, fitted(model), G = 10)
#'
#' # Example 2: Grouped data (original Farrington test)
#' set.seed(456)
#' n_groups <- 50
#' m_trials <- sample(5:20, n_groups, replace = TRUE)
#' x_grouped <- rnorm(n_groups)
#' linpred_grouped <- -0.5 + 1.0 * x_grouped
#' prob_grouped <- 1 / (1 + exp(-linpred_grouped))
#' y_grouped <- rbinom(n_groups, m_trials, prob_grouped)
#'
#' # Fit model for grouped data
#' data_grouped <- data.frame(successes = y_grouped, trials = m_trials, x = x_grouped)
#' model_grouped <- glm(cbind(successes, trials - successes) ~ x,
#'                      data = data_grouped, family = binomial())
#' predicted_probs_grouped <- fitted(model_grouped)
#'
#' # Original Farrington test. G = NULL is required: left at its default of 10 the
#' # call takes the automatic-grouping branch instead, which ignores 'model' and 'm'
#' # and refers binomial counts to the binary statistic.
#' result_grouped <- ef.gof(y_grouped, predicted_probs_grouped,
#'                          model = model_grouped, m = m_trials,
#'                          G = NULL)
#' print(result_grouped)
#'
#' @seealso
#' \code{\link[ResourceSelection]{hoslem.test}} for the Hosmer-Lemeshow test
#'
#' @concept goodness-of-fit
#' @concept calibration
#' @concept logistic regression
#' @concept model diagnostics
#' @concept Ebrahim-Farrington
#' @concept sparse data
#' @concept Hosmer-Lemeshow
#' @export
ef.gof <- function(y, predicted_probs = NULL, model = NULL, m = NULL, G = 10,
                   method = NULL, reference = c("chisq", "normal"), X = NULL,
                   groups = NULL) {

  reference <- match.arg(reference)
  thesis <- FALSE
  if (!is.null(method)) {
    method <- match.arg(method, c("chisq", "normal"))
    if (method == "normal") {
      thesis <- TRUE
      warning("ef.gof: method = \"normal\" is the standardization (EF - (G - 2)) / sqrt(2 (G - 2)) ",
              "of versions up to 2.9.0, kept to reproduce earlier results. The normal reference of ",
              "the EF paper is reference = \"normal\".", call. = FALSE)
    }
  }

  # The first argument may be a fitted binomial glm; extract y, predictions and design.
  if (inherits(y, "glm")) {
    if (y$family$family != "binomial") stop("ef.gof: the model must be a binomial glm.")
    if (is.null(X)) X <- stats::model.matrix(y)
    predicted_probs <- as.numeric(y$fitted.values)   # fitted() pads NAs under na.exclude
    y <- as.numeric(y$y)
  }
  if (is.null(predicted_probs)) {
    stop("Provide 'predicted_probs' (or pass a fitted glm as the first argument).")
  }

  # Input validation
  n <- length(y)
  if (n != length(predicted_probs)) {
    stop("Length mismatch: 'y' and 'predicted_probs' must have the same length")
  }

  if (!is.numeric(y) || !is.numeric(predicted_probs)) {
    stop("'y' and 'predicted_probs' must be numeric vectors")
  }

  if (any(predicted_probs < 0) || any(predicted_probs > 1)) {
    stop("'predicted_probs' must be between 0 and 1")
  }

  if (!all(y %in% c(0, 1)) && is.null(m)) {
    stop("For binary data, 'y' must contain only 0s and 1s")
  }

  # Determine which test to perform
  if (!is.null(groups) || !is.null(G)) {
    # Ebrahim-Farrington test on risk groups
    if (!is.null(model) || !is.null(m)) {
      warning("Ebrahim-Farrington grouping is in force (G is not NULL), so 'model' and ",
              "'m' are ignored. For the original Farrington test on grouped data, set ",
              "G = NULL.")
    }

    if (is.null(groups)) {
      if (identical(G, "auto")) G <- max(10, ceiling(n / 25))
      else if (identical(G, "paul")) G <- .ef_paul_G(n, y)
      if (!is.numeric(G) || length(G) != 1 || is.na(G) || G < 2 || G != round(G)) {
        stop("'G' must be a single integer >= 2, \"auto\" or \"paul\"")
      }
      if (G > n) {
        stop("Number of groups 'G' cannot exceed sample size 'n'")
      }
      g <- .ef_groups(predicted_probs, G)
    } else {
      g <- .ef_labels(groups, predicted_probs)
    }

    if (reference == "normal" && !thesis) {
      if (is.null(X))
        stop("reference = \"normal\" needs the model's design matrix: pass the fitted glm as ",
             "the first argument, or X.")
      X <- as.matrix(X)
      if (nrow(X) != n) stop("'X' must have one row per observation.")
    }

    return(.ebrahim_farrington_grouped(y, predicted_probs, g,
                                       if (thesis) "thesis" else reference, X))

  } else {
    # Original Farrington test for grouped data
    if (is.null(m)) {
      stop("For the original Farrington test without grouping, ",
           "you must provide the 'm' vector (trial counts). ",
           "For binary data, use G parameter for automatic grouping.")
    }

    if (is.null(model)) {
      stop("For the original Farrington test with grouped data, ",
           "you must provide the fitted 'model' object.")
    }

    return(.farrington_original(y, predicted_probs, model, m))
  }
}

# Equal-size groups by rank of the fitted risk: sizes differ by at most one, and tied risks are
# put in random order, so that a boundary inside a covariate pattern does not follow the row order
# (or the response, for binomial data expanded events first). The uniforms are drawn only when
# some risks are tied, so untied data leave the random number stream alone.
.ef_groups <- function(p, G) {
  n <- length(p)
  o <- if (anyDuplicated(p)) order(p, stats::runif(n)) else order(p)
  g <- integer(n)
  g[o] <- ceiling(seq_len(n) * G / n)
  g
}

# Group labels given by the user, or one group per distinct fitted risk ("pattern"), numbered
# 1..G in order of the label (for "pattern", of the risk).
.ef_labels <- function(groups, p) {
  if (identical(groups, "pattern")) groups <- p
  if (length(groups) != length(p) || anyNA(groups))
    stop("'groups' must be \"pattern\" or one label per observation, without missing values.")
  g <- as.integer(factor(groups))
  if (max(g) < 3)
    stop("'groups' gives fewer than 3 groups; the test needs at least 3.")
  if (max(g) == length(p))
    stop("'groups' puts every record in a group of its own; the grouped statistic is then ",
         "the sample size whatever the data.")
  g
}

# The number of groups of Paul, Pennell and Lemeshow (2013).
.ef_paul_G <- function(n, y) {
  n1 <- sum(y)
  floor(max(10, min(n1 / 2, (n - n1) / 2, 2 + 8 * (n / 1000)^2)))
}

# The size rule of the chi-square reference: about one expected event and non-event per group.
.ef_sparse_warning <- function(emin) {
  w <- simpleWarning(sprintf(paste(
    "ef.gof: the smallest expected number of events or non-events in a group is %.2f;",
    "below about one half the chi-square reference is unreliable and EF is conservative.",
    "Use fewer groups."), emin))
  class(w) <- c("ef_sparse_groups", class(w))
  warning(w)
}

# The grouped statistic, its reference and the per-group read-out.
.ebrahim_farrington_grouped <- function(y, predicted_probs, g, reference = "chisq", X = NULL) {
  G  <- max(g)
  ng <- tabulate(g, G)
  o  <- as.numeric(rowsum(y, g, reorder = TRUE))
  e  <- as.numeric(rowsum(predicted_probs, g, reorder = TRUE))

  # Apply numerical stability bounds to the group means
  pb <- pmax(pmin(e / ng, 0.999999), 0.000001)
  V  <- ng * pb * (1 - pb)
  r  <- (o - ng * pb) / sqrt(V)
  b  <- (1 - 2 * pb) / sqrt(V)

  HL <- sum(r^2)
  C  <- sum(b * r)
  EF <- HL - C

  df <- NA_real_
  if (reference == "chisq") {
    emin <- min(pmin(e, ng - e))
    if (emin < 0.5) .ef_sparse_warning(emin)
    stat <- EF
    df   <- as.numeric(G - 2)
    p_value <- stats::pchisq(EF, df = G - 2, lower.tail = FALSE)
  } else if (reference == "normal") {
    mom  <- .ef_normal_moments(predicted_probs, g, X, ng, pb, V)
    stat <- (EF - mom$mu) / sqrt(mom$var)
    p_value <- stats::pnorm(stat, lower.tail = FALSE)
  } else {
    # the standardization of versions up to 2.9.0 (method = "normal")
    if (G - 2 <= 0) {
      stat <- 0
      p_value <- 1
      warning("Variance too small for reliable test. Consider using more groups.")
    } else {
      stat <- (EF - (G - 2)) / sqrt(2 * (G - 2))
      p_value <- 1 - stats::pnorm(stat)
    }
  }

  out <- data.frame(Test = "Ebrahim-Farrington", Test_Statistic = stat, df = df,
                    Reference = reference, p_value = p_value, G = G, HL = HL, C = C,
                    stringsAsFactors = FALSE)
  attr(out, "groups") <- data.frame(group = seq_len(G), n = ng, observed = o, expected = e,
                                    pbar = pb, r = r, b = b, c = b * r)
  out
}

# Null mean and variance of EF with many groups (Theorem 3 of the EF paper), at the fitted risks:
# tau_g = sum v_i, U and K the group sums of v_i x_i and v_i (1 - 2 p_i) x_i, M = X'VX. The mean
# subtracts the trace of the estimation term, reduced by 1 / mean(n_g); the variance carries the
# within-group third and fourth moments and the estimation terms.
.ef_normal_moments <- function(p, g, X, ng, pb, V) {
  v   <- p * (1 - p)
  cc  <- 1 - 2 * pb
  tau <- as.numeric(rowsum(v, g, reorder = TRUE))
  k3  <- as.numeric(rowsum(v * (1 - 2 * p), g, reorder = TRUE))
  k4  <- as.numeric(rowsum(v * (1 - 6 * v), g, reorder = TRUE))
  U   <- rowsum(v * X, g, reorder = TRUE)
  K   <- rowsum(v * (1 - 2 * p) * X, g, reorder = TRUE)
  Mi  <- solve(crossprod(X, v * X))
  dtr <- sum(diag(crossprod(U, U / V) %*% Mi))
  k1  <- colSums((tau / V) * cc * U / V)
  dd  <- k1 - colSums(cc * U / V)
  k3v <- colSums((K - cc * U) / V)
  vEF <- sum((2 * tau^2 + k4 - 2 * cc * k3 + cc^2 * tau) / V^2) -
    2 * drop(dd %*% Mi %*% k3v) + drop(dd %*% Mi %*% dd)
  list(mu = sum(tau / V) - dtr * (1 - 1 / mean(ng)), var = vEF)
}

# Internal function for original Farrington test
.farrington_original <- function(y, predicted_probs, model, m) {
  n <- length(y)

  # Validate inputs
  if (length(m) != n) {
    stop("Length mismatch: 'y' and 'm' must have the same length")
  }

  if (!inherits(model, "glm")) {
    stop("'model' must be a glm object")
  }

  # Apply numerical stability bounds
  prob_bounded <- pmax(pmin(predicted_probs, 0.999999), 0.000001)

  # Get design matrix
  design <- stats::model.matrix(model)
  p <- ncol(design)

  # Calculate variance terms
  variance_terms <- m * prob_bounded * (1 - prob_bounded)

  # Calculate test statistic components
  pearson_components <- ((y - m * prob_bounded)^2) / variance_terms
  farrington_correction <- ((1 - 2 * prob_bounded) / variance_terms) *
                          (y - m * prob_bounded)

  test_statistic_raw <- sum(pearson_components) - sum(farrington_correction)

  # Calculate expected value and variance for original Farrington test
  # This requires the information matrix
  tryCatch({
    w_matrix <- diag(as.vector(variance_terms))
    xpxi <- solve(t(design) %*% w_matrix %*% design)
    varlinpr <- diag(design %*% xpxi %*% t(design))

    expected_value <- (n - p) + sum(prob_bounded * (1 - prob_bounded) * varlinpr)
    variance_value <- 2 * (1 - (p / n)) * sum((m - 1) / m)

  }, error = function(e) {
    warning("Could not compute full Farrington variance. Using simplified version.")
    expected_value <- n - p
    variance_value <- 2 * (n - p) * mean((m - 1) / m)
  })

  # Handle cases where test is not applicable
  if (variance_value <= 1e-8 || is.na(variance_value)) {
    z_statistic <- 0
    p_value <- 1
    if (all(m == 1)) {
      warning("Farrington test is not applicable to purely binary data (all m=1). ",
              "Consider using G parameter for automatic grouping.")
    }
  } else {
    z_statistic <- (test_statistic_raw - expected_value) / sqrt(variance_value)
    p_value <- 1 - stats::pnorm(z_statistic)
  }

  # Return results
  data.frame(
    Test = "Farrington-Original",
    Test_Statistic = z_statistic,
    p_value = p_value,
    stringsAsFactors = FALSE
  )
}
