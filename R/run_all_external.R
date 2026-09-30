#' Run the External-Validation Tests at Once
#'
#' @description
#' Checks the calibration of predictions that were made \emph{without} the data at
#' hand -- a published model, or any model, applied to new patients -- and returns one
#' tidy \code{data.frame}, one row per test, in the format of \code{\link{run.all.gof}}.
#' Only the outcomes \code{y} and the predicted probabilities \code{p} are needed, so the
#' predictions may come from a logistic regression, a random forest, a neural network or
#' a clinical score: nothing is refitted.
#'
#' @details
#' \strong{Why a separate battery.} \code{run.all.gof(y, predicted_probs = p)} treats the
#' predictions as fitted to \code{y}, and each reference distribution then allows for the
#' parameters the fit spent. When the predictions are frozen nothing was spent: the
#' Hosmer--Lemeshow statistic is referred to \eqn{\chi^2_G}, not \eqn{\chi^2_{G-2}}; Stukel's
#' terms are added to the frozen linear predictor as an offset; and the directed test
#' takes \eqn{\Omega = I} with a constant column, because no score equation absorbs the
#' overall level. Using the internal references on frozen predictions makes every one of
#' these tests conservative.
#'
#' \strong{The tests} (\code{Family} in brackets).
#' \itemize{
#'   \item \code{EDGE (G=...)} [Directed] -- the directed grouped test in external mode,
#'     \code{def.gof(y, predicted_probs = p, external = TRUE)}: the cubic basis plus a
#'     constant, four degrees of freedom. Run at \code{G} (ten by default, the setting
#'     that keeps its level when a few records carry corrupted predictors) and at
#'     \code{G = "auto"}, \eqn{\max(10, \lceil n/25 \rceil)}, which has more power on
#'     clean data.
#'   \item \code{Cox recalibration} [Calibration] -- the likelihood-ratio test that the
#'     calibration intercept is 0 and the slope 1 in \code{glm(y ~ logit(p))}
#'     (Cox 1958; Miller et al. 1991). The Note gives both estimates.
#'   \item \code{Calibration-in-the-large} [Calibration] -- the score test that the
#'     intercept is 0 with the slope held at 1: \eqn{(O - E)/\sqrt{\sum p(1-p)}}.
#'   \item \code{Spiegelhalter-z} [Calibration] -- Spiegelhalter's (1986) \eqn{z} test,
#'     two-sided.
#'   \item \code{GiViTI} [Calibration] -- the GiViTI calibration test in its external
#'     mode (Nattino et al. 2014), run in an isolated \pkg{callr} process; needs
#'     \pkg{givitiR} and \pkg{callr}. When it selects a polynomial of degree one it
#'     returns the same \eqn{p}-value as the Cox test.
#'   \item \code{Hosmer-Lemeshow (external)} [Partition] -- \eqn{G} equal-frequency risk
#'     groups, referred to \eqn{\chi^2_G}.
#'   \item \code{Stukel (offset)} [Directed] -- Stukel's two terms added to the frozen
#'     linear predictor, likelihood-ratio \eqn{\chi^2_2}.
#'   \item \code{le Cessie (external)} [Smoothing] -- le Cessie and van Houwelingen's
#'     kernel statistic over covariate space with \eqn{\Omega = I}; only when \code{X}
#'     is given and \code{include_slow = TRUE}, because it builds an \eqn{n \times n}
#'     kernel.
#' }
#' Three descriptive rows carry no \eqn{p}-value: the ratio of observed to expected
#' events, the calibration slope, and the c-statistic (area under the ROC curve).
#'
#' \strong{Reading the panel.} In a simulation of external validation (EDGE paper,
#' Supporting Information) the directed test at ten groups had the power of the Cox test
#' and the GiViTI belt on average, led them on curved departures and trailed them on a
#' shift, and kept its level with ten reversed predictions in 1000 records, where the
#' Cox test and the belt did not. The Cox test says whether the predictions need
#' recalibrating; the directed test says whether a recalibration would be enough.
#'
#' @param y Binary (0/1) outcomes of the validation sample.
#' @param p Predicted probabilities for the same records, produced without them.
#' @param G Number of risk groups for the directed and Hosmer--Lemeshow tests
#'   (default \code{10}).
#' @param X Optional covariate matrix or data frame, for le Cessie's test only.
#' @param include_slow Logical; run le Cessie's test when \code{X} is given
#'   (default \code{FALSE}).
#'
#' @return A \code{data.frame} of class \code{gof_battery} with columns \code{Test},
#'   \code{Family}, \code{Statistic}, \code{df}, \code{p_value} and \code{Note}, printed
#'   by the same method as \code{\link{run.all.gof}}.
#'
#' @seealso \code{\link{def.gof}} (its \code{external} argument), \code{\link{run.all.gof}}.
#'
#' @author Ebrahim Khaled Ebrahim \email{ebrahimkhaled@@alexu.edu.eg}
#'
#' @examples
#' set.seed(1)
#' n <- 1000
#' x <- rnorm(n)
#' p <- plogis(-1 + 0.8 * x)                 # a published model, frozen
#' y <- rbinom(n, 1, plogis(-1 + 0.6 * x))    # new patients: the model is overfitted
#' run.all.external(y, p)
#'
#' @export
run.all.external <- function(y, p, G = 10, X = NULL, include_slow = FALSE) {
  y <- as.numeric(y); p <- as.numeric(p)
  n <- length(y)
  if (!all(y %in% c(0, 1))) stop("'y' must be binary (0/1).")
  if (length(p) != n) stop("'y' and 'p' have different lengths.")
  if (anyNA(y) || anyNA(p)) stop("'y' and 'p' must not contain missing values.")
  if (any(p <= 0 | p >= 1)) stop("'p' must lie strictly between 0 and 1.")
  if (!is.null(X) && NROW(X) != n) stop("'X' must have one row per record.")

  lp <- stats::qlogis(p); v <- p * (1 - p)
  dev0 <- -2 * sum(y * log(p) + (1 - y) * log(1 - p))
  row <- function(test, family, stat, df, pv, note = "")
    data.frame(Test = test, Family = family, Statistic = stat, df = df, p_value = pv,
               Note = note, stringsAsFactors = FALSE)
  safe <- function(test, family, expr)
    tryCatch(expr, error = function(e) row(test, family, NA_real_, NA_real_, NA_real_,
                                           paste("Not run:", conditionMessage(e))))
  rows <- list()

  # --- the directed test in external mode, at G and at the rule G ---
  for (g in unique(list(G, "auto"))) {
    lab <- if (identical(g, "auto")) sprintf("EDGE (G=auto, %d)", .def_auto_G(n)) else sprintf("EDGE (G=%d)", g)
    rows[[lab]] <- safe(lab, "Directed", {
      e <- suppressWarnings(def.gof(y, predicted_probs = p, G = g, external = TRUE))
      row(lab, "Directed", e$Test_Statistic, e$df, e$p_value, "external mode, cubic basis + constant")
    })
  }

  # --- calibration: Cox recalibration, calibration in the large, Spiegelhalter, GiViTI ---
  rows$cox <- safe("Cox recalibration", "Calibration", {
    f1 <- suppressWarnings(stats::glm(y ~ lp, family = stats::binomial()))
    S <- dev0 - f1$deviance; b <- stats::coef(f1)
    row("Cox recalibration", "Calibration", S, 2, stats::pchisq(S, 2, lower.tail = FALSE),
        sprintf("intercept %.3f, slope %.3f", b[[1]], b[[2]]))
  })
  zl <- (sum(y) - sum(p)) / sqrt(sum(v))
  rows$citl <- row("Calibration-in-the-large", "Calibration", zl^2, 1,
                   stats::pchisq(zl^2, 1, lower.tail = FALSE), sprintf("O/E = %.3f", sum(y) / sum(p)))
  zs <- sum((y - p) * (1 - 2 * p)) / sqrt(sum((1 - 2 * p)^2 * v))
  rows$spz <- row("Spiegelhalter-z", "Calibration", zs, NA_real_, 2 * stats::pnorm(-abs(zs)), "two-sided z")
  rows$giviti <- gof_giviti(list(y = y, ph = p), list(devel = "external"))
  rows$giviti <- row("GiViTI", "Calibration", rows$giviti$Statistic, rows$giviti$df,
                     rows$giviti$p_value, rows$giviti$Note)

  # --- the partition test in its external form ---
  rows$hl <- safe("Hosmer-Lemeshow (external)", "Partition", {
    grp <- pmin(ceiling(rank(p, ties.method = "first") / (n / G)), G)
    O <- as.numeric(rowsum(y, grp)); E <- as.numeric(rowsum(p, grp)); ng <- tabulate(grp, G)
    S <- sum((O - E)^2 / (E * (1 - E / ng)))              # Hosmer and Lemeshow's denominator
    row("Hosmer-Lemeshow (external)", "Partition", S, G, stats::pchisq(S, G, lower.tail = FALSE),
        sprintf("%d groups, chi-square(%d)", G, G))
  })

  # --- Stukel's two terms on the frozen linear predictor ---
  rows$stukel <- safe("Stukel (offset)", "Directed", {
    za <- 0.5 * lp^2 * (lp >= 0); zb <- -0.5 * lp^2 * (lp < 0)
    f2 <- suppressWarnings(stats::glm(y ~ 0 + za + zb, offset = lp, family = stats::binomial()))
    S <- dev0 - f2$deviance
    row("Stukel (offset)", "Directed", S, 2, stats::pchisq(S, 2, lower.tail = FALSE), "terms added to logit(p)")
  })

  # --- le Cessie's kernel statistic with Omega = I (optional, O(n^2) memory) ---
  if (!is.null(X) && isTRUE(include_slow)) {
    rows$lc <- safe("le Cessie (external)", "Smoothing", {
      Xs <- scale(as.matrix(as.data.frame(lapply(as.data.frame(X), as.numeric))))
      D <- as.matrix(stats::dist(Xs) * sqrt(0.5))
      R <- pmax(1 - D / mean(D), 0)
      r <- y - p
      Q <- sum(as.numeric(r %*% R) * r)
      # moments of r'Rr for independent Bernoulli residuals (le Cessie and van Houwelingen 1995, A.6)
      EQ <- sum(diag(R) * v)
      VarQ <- sum(diag(R)^2 * (v * (1 - 3 * v) - 3 * v^2)) + 2 * as.numeric(v %*% (R * R) %*% v)
      stat <- Q * 2 * EQ / VarQ; df <- 2 * EQ^2 / VarQ
      row("le Cessie (external)", "Smoothing", stat, df, stats::pchisq(stat, df, lower.tail = FALSE),
          "kernel over covariate space, Omega = I")
    })
  } else if (!is.null(X)) {
    rows$lc <- row("le Cessie (external)", "Smoothing", NA_real_, NA_real_, NA_real_,
                   "Not run: set include_slow = TRUE (n x n kernel)")
  }

  # --- descriptive rows ---
  slope <- tryCatch(unname(stats::coef(suppressWarnings(stats::glm(y ~ lp, family = stats::binomial())))[2]),
                    error = function(e) NA_real_)
  rk <- rank(p); n1 <- sum(y); n0 <- n - n1
  auc <- if (n1 > 0 && n0 > 0) (sum(rk[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0) else NA_real_
  rows$oe  <- row("O/E ratio", "Descriptive", sum(y) / sum(p), NA_real_, NA_real_, "1 = calibrated in the large")
  rows$sl  <- row("Calibration slope", "Descriptive", slope, NA_real_, NA_real_, "1 = calibrated; < 1 = overfitted")
  rows$auc <- row("c-statistic (AUC)", "Descriptive", auc, NA_real_, NA_real_, "discrimination, not calibration")

  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  class(out) <- c("gof_battery", "data.frame")
  out
}
