#' Cubic Calibration Likelihood-Ratio Test for Logistic Regression
#'
#' @description
#' \code{cubic.calib.gof()} asks whether the calibration curve of a binary risk
#' model bends on the logit scale. It fits two logistic regressions of the outcome
#' on the logit of the model's own predicted risk, \eqn{\hat\eta}: one with
#' \eqn{\hat\eta} alone, and one that adds \eqn{\hat\eta^2} and \eqn{\hat\eta^3}.
#' Twice the drop in deviance is referred to a chi-squared distribution on 2
#' degrees of freedom.
#'
#' @details
#' \strong{Where the test comes from.} This is not a new test, and the package
#' does not claim it. Testing a fitted model by adding powers of its own fitted
#' values goes back to Tukey's (1949) one degree of freedom for non-additivity,
#' which adds the square of the fitted values, and Ramsey (1969) made it a
#' general specification test for linear regression, RESET, by adding their
#' square and cube. Pregibon (1980) carried the idea to generalized linear
#' models as a goodness-of-link test: refit the model with constructed
#' variables computed from its fitted linear predictor, and read the change in
#' deviance. Morgan (1985) studied the logistic model whose linear predictor is a
#' cubic polynomial, the cubic logistic model, whose extra terms bend the two
#' tails of the logistic curve. \code{cubic.calib.gof()} is this family's test
#' for a logistic model: the constructed variables are \eqn{\hat\eta^2} and
#' \eqn{\hat\eta^3}, and the test is the likelihood-ratio test of their two
#' coefficients. The GiViTI calibration belt (Finazzi et al., 2011; Nattino et
#' al., 2014) fits the same kind of polynomial in the logit of the predicted
#' risk, choosing its degree from the data and testing all of its
#' coefficients; this function fixes the degree at three and tests only the
#' curvature.
#'
#' \strong{What it tests.} When the model was fitted by maximum likelihood to the
#' same data, the logistic regression of \code{y} on \eqn{\hat\eta} reproduces
#' the fit exactly (intercept 0, slope 1), so the test asks only whether the
#' calibration curve is curved on the logit scale: the quadratic term picks up
#' an asymmetric bend, as from a complementary log-log truth fitted by a logit,
#' and the cubic term a symmetric one, as from tails that are too heavy or too
#' light. With the predictions of an external model (\code{predicted_probs}) the
#' intercept and slope are free as well, so the test still asks only about
#' curvature; calibration-in-the-large and the calibration slope need their own
#' tests (Cox, 1958).
#'
#' \strong{What to know before using it.} The test weighs every record by its own
#' linear predictor. A few records with corrupted covariates, whose predictions
#' are pushed to the extremes, can therefore drive its false-alarm rate far above
#' the nominal level; the grouped directed test \code{\link{edge.gof}} is
#' protected against that. The polynomial is fitted on the standardized
#' \eqn{\hat\eta}, which spans the same polynomials and leaves the statistic
#' unchanged but keeps the fit well conditioned.
#'
#' @param object Either a fitted \code{glm} with \code{family = binomial} and a
#'   0/1 response, or a numeric 0/1 outcome vector, in which case
#'   \code{predicted_probs} is required.
#' @param predicted_probs Predicted probabilities for \code{object} when it is
#'   an outcome vector, for example from an external model. Ignored when
#'   \code{object} is a \code{glm}.
#'
#' @return A one-row \code{data.frame} with columns \code{Test}
#'   (\code{"Cubic-LR"}), \code{Test_Statistic} (the likelihood-ratio
#'   statistic), \code{df} (2, or fewer if a column cannot be estimated),
#'   \code{Method} (\code{"likelihood ratio"}) and \code{p_value}. A fit that
#'   fails or does not converge gives an \code{NA} p-value with a warning.
#'
#' @references
#' Cox, D. R. (1958). Two further applications of a model for binary regression.
#' \emph{Biometrika}, 45(3-4), 562-565. \doi{10.1093/biomet/45.3-4.562}
#'
#' Finazzi, S., Poole, D., Luciani, D., Cogo, P. E. and Bertolini, G. (2011).
#' Calibration belt for quality-of-care assessment based on dichotomous
#' outcomes. \emph{PLoS ONE}, 6(2), e16110. \doi{10.1371/journal.pone.0016110}
#'
#' Morgan, B. J. T. (1985). The cubic logistic model for quantal assay data.
#' \emph{Journal of the Royal Statistical Society, Series C (Applied
#' Statistics)}, 34(2), 105-113. \doi{10.2307/2347362}
#'
#' Nattino, G., Finazzi, S. and Bertolini, G. (2014). A new calibration test and
#' a reappraisal of the calibration belt for the assessment of prediction models
#' based on dichotomous outcomes. \emph{Statistics in Medicine}, 33(14),
#' 2390-2407. \doi{10.1002/sim.6100}
#'
#' Pregibon, D. (1980). Goodness of link tests for generalized linear models.
#' \emph{Journal of the Royal Statistical Society, Series C (Applied
#' Statistics)}, 29(1), 15-23. \doi{10.2307/2346405}
#'
#' Ramsey, J. B. (1969). Tests for specification errors in classical linear
#' least-squares regression analysis. \emph{Journal of the Royal Statistical
#' Society, Series B}, 31(2), 350-371.
#' \doi{10.1111/j.2517-6161.1969.tb00796.x}
#'
#' Tukey, J. W. (1949). One degree of freedom for non-additivity.
#' \emph{Biometrics}, 5(3), 232. \doi{10.2307/3001938}
#'
#' @author Ebrahim Khaled Ebrahim \email{ebrahimkhaled@@alexu.edu.eg}
#'
#' @examples
#' set.seed(1)
#' x <- runif(1000, -3, 3)
#' y <- rbinom(1000, 1, plogis(0.6 * x))
#' fit <- glm(y ~ x, family = binomial())
#' cubic.calib.gof(fit)                             # correct model
#'
#' y2 <- rbinom(1000, 1, 1 - exp(-exp(0.6 * x)))       # complementary log-log truth
#' cubic.calib.gof(glm(y2 ~ x, family = binomial()))   # fitted by a logit: the curve bends
#'
#' ## frozen predictions, e.g. an external validation
#' cubic.calib.gof(y, predicted_probs = fitted(fit))
#'
#' @seealso \code{\link{edge.gof}} for the grouped directed test,
#'   \code{\link{run.all.gof}}, where this test is the row \code{"Cubic-LR"}.
#' @concept goodness-of-fit
#' @concept calibration
#' @concept logistic regression
#' @concept link test
#' @export
cubic.calib.gof <- function(object, predicted_probs = NULL) {
  if (inherits(object, "glm")) {
    if (object$family$family != "binomial")
      stop("cubic.calib.gof: the model must be a binomial glm.")
    y <- as.numeric(object$y)
    eta <- if (identical(object$family$link, "logit")) as.numeric(object$linear.predictors)
           else stats::qlogis(pmin(pmax(as.numeric(object$fitted.values), 1e-12), 1 - 1e-12))
  } else {
    if (!is.numeric(object) || is.null(predicted_probs))
      stop("cubic.calib.gof: pass a binomial glm, or a 0/1 outcome vector with 'predicted_probs'.")
    y <- as.numeric(object)
    if (length(predicted_probs) != length(y))
      stop("cubic.calib.gof: 'object' and 'predicted_probs' have different lengths.")
    if (any(!is.finite(predicted_probs)) || any(predicted_probs < 0 | predicted_probs > 1))
      stop("cubic.calib.gof: 'predicted_probs' must be probabilities.")
    eta <- stats::qlogis(pmin(pmax(as.numeric(predicted_probs), 1e-12), 1 - 1e-12))
  }
  if (!all(y %in% c(0, 1)))
    stop("cubic.calib.gof: the outcome must be binary (0/1).")
  r <- .cubic_lr(y, eta)
  if (nzchar(r$note)) warning("cubic.calib.gof: ", r$note, call. = FALSE)
  data.frame(Test = "Cubic-LR", Test_Statistic = r$stat, df = r$df,
             Method = "likelihood ratio", p_value = r$p, stringsAsFactors = FALSE)
}

# The likelihood-ratio statistic of eta^2 and eta^3 added to a logistic regression of y on eta, fitted on the
# standardized eta (the same column space, so the same deviances, and a better conditioned fit).
.cubic_lr <- function(y, eta) {
  none <- function(note) list(stat = NA_real_, df = NA_real_, p = NA_real_, note = note)
  s <- stats::sd(eta)
  if (!is.finite(s) || s == 0) return(none("the predicted risks are all equal"))
  z <- (eta - mean(eta)) / s
  base <- tryCatch(suppressWarnings(stats::glm.fit(cbind(1, z), y, family = stats::binomial())),
                   error = function(e) NULL)
  full <- tryCatch(suppressWarnings(stats::glm.fit(cbind(1, z, z^2, z^3), y, family = stats::binomial())),
                   error = function(e) NULL)
  if (is.null(base) || is.null(full)) return(none("a calibration fit failed"))
  if (!isTRUE(base$converged) || !isTRUE(full$converged))
    return(none("a calibration fit did not converge (separation?)"))
  df <- full$rank - base$rank
  if (df < 1) return(none("the added columns cannot be estimated"))
  stat <- max(base$deviance - full$deviance, 0)
  list(stat = stat, df = df, p = stats::pchisq(stat, df, lower.tail = FALSE), note = "")
}

# run.all.gof() row: the logit of the fitted risk, from the model's linear predictor when it is a logit fit
gof_cubic_lr <- function(ctx, opts = list()) {
  eta <- if (isTRUE(ctx$has_model) && identical(ctx$model$family$link, "logit")) as.numeric(ctx$model$linear.predictors)
         else stats::qlogis(ctx$ph)
  r <- .cubic_lr(as.numeric(ctx$y), eta)
  list(Statistic = r$stat, df = r$df, p_value = r$p,
       Note = if (nzchar(r$note)) paste("Not run:", r$note) else "eta^2 and eta^3 added to a logit fit on eta")
}
