#' Run a Battery of Goodness-of-Fit Tests at Once
#'
#' @description
#' Runs several goodness-of-fit tests for a binary logistic regression in one
#' call and returns one tidy \code{data.frame}, one row per test. Pass a fitted
#' \code{glm} to run the whole battery; pass \code{(y, predicted_probs)} to run
#' the tests that need only predictions. Each test is wrapped so that a failure of
#' one test never aborts the whole run.
#'
#' @details
#' The currently bundled tests are: \code{Pearson} and \code{Deviance} (global),
#' \code{HL} (Hosmer-Lemeshow deciles) and \code{HL-equalwidth} (partition),
#' \code{EF} (the omnibus Ebrahim-Farrington test), \code{DEF.poly2/poly3/stukel}
#' and \code{Stukel} (directed), plus the two ensemble rows
#' (\code{Ensemble.Vote(3DEF)} and \code{Ensemble.Univ(3DEF+EF)}) from the Cauchy
#' combination test. Further thesis tests (Osius-Rojek, McCullagh, information
#' matrix, Copas, Tsiatis, Xie, Pulkstenis-Robinson, le Cessie) are planned for a
#' later build.
#'
#' @param object A fitted binary logistic \code{\link[stats]{glm}}, or a binary
#'   (0/1) response vector \code{y} (then supply \code{predicted_probs}).
#' @param predicted_probs Numeric predicted probabilities; required when
#'   \code{object} is a \code{y} vector.
#' @param X Optional design matrix; lets the directed (DEF) tests run from the
#'   \code{(y, predicted_probs)} form.
#' @param tests Either \code{"all"} (default) or a character vector of test names
#'   to run (e.g. \code{c("EF","DEF.poly3","HL")}).
#' @param G Integer number of groups passed to the grouping tests (default 10).
#' @param include_slow Logical; reserved for opt-in slow (bootstrap/GAM) tests.
#'   No slow tests are bundled yet, so this currently has no effect.
#' @param control Optional named list of per-test options (reserved).
#'
#' @return A \code{data.frame} with columns \code{Test}, \code{Family},
#'   \code{Statistic}, \code{df}, \code{p_value}, and \code{Note}.
#'
#' @author Ebrahim Khaled Ebrahim \email{ebrahimkhaled@@alexu.edu.eg}
#'
#' @examples
#' set.seed(1)
#' n <- 500
#' x <- runif(n, -3, 3)
#' y <- rbinom(n, 1, 1 / (1 + exp(-(0.6 * x))))
#' fit <- glm(y ~ x, family = binomial())
#' run.all.gof(fit)                       # the whole battery + ensemble rows
#' run.all.gof(fit, tests = c("EF", "DEF.poly3", "HL"))
#' run.all.gof(y, fitted(fit))            # prediction-only tests
#'
#' @seealso \code{\link{ef.gof}}, \code{\link{def.gof}}, \code{\link{def.ensemble.gof}}.
#' @importFrom stats fitted predict model.matrix model.frame coef deviance pchisq binomial glm.fit
#' @export
run.all.gof <- function(object, predicted_probs = NULL, X = NULL,
                        tests = "all", G = 10, include_slow = FALSE,
                        control = list()) {

  ctx <- .gof_context(object, predicted_probs, X, G = G)
  sel <- if (identical(tests, "all")) names(.GOF_REGISTRY) else intersect(tests, names(.GOF_REGISTRY))
  if (length(sel) == 0) stop("run.all.gof: no known tests selected.")

  rows <- list(); skipped_model <- FALSE
  for (nm in sel) {
    e <- .GOF_REGISTRY[[nm]]
    if (isTRUE(e$slow) && !include_slow) next
    if (isTRUE(e$needs_model) && !ctx$has_model && is.null(ctx$X)) {
      skipped_model <- TRUE
      rows[[nm]] <- data.frame(Test = nm, Family = e$family, Statistic = NA_real_,
                               df = NA_real_, p_value = NA_real_,
                               Note = "needs a glm model", stringsAsFactors = FALSE)
      next
    }
    res <- tryCatch(e$fn(ctx, control[[nm]]),
                    error = function(err) list(Statistic = NA, df = NA, p_value = NA,
                                               Note = paste("error:", conditionMessage(err))))
    rows[[nm]] <- data.frame(Test = nm, Family = e$family,
                             Statistic = .gof_num(res$Statistic), df = .gof_num(res$df),
                             p_value = .gof_num(res$p_value),
                             Note = if (is.null(res$Note)) "" else res$Note,
                             stringsAsFactors = FALSE)
  }
  out <- do.call(rbind, rows)

  # ensemble rows (only when a model is available and running the full set)
  if (ctx$has_model && identical(tests, "all")) {
    v3 <- tryCatch(def.ensemble.gof(ctx$model, G = G)$p_value, error = function(e) NA_real_)
    vu <- tryCatch(def.ensemble.gof(ctx$model, add_ef = TRUE, G = G)$p_value, error = function(e) NA_real_)
    out <- rbind(out, data.frame(
      Test = c("Ensemble.Vote(3DEF)", "Ensemble.Univ(3DEF+EF)"), Family = "Ensemble",
      Statistic = NA_real_, df = NA_real_, p_value = c(v3, vu), Note = "CCT",
      stringsAsFactors = FALSE))
  }

  if (skipped_model)
    message("Only prediction-based tests were run. Pass the fitted glm (or X) ",
            "to also run the directed and refit-based tests.")
  rownames(out) <- NULL
  out
}

# Coerce a possibly-NULL/empty scalar to a numeric (NA if missing).
.gof_num <- function(x) if (is.null(x) || length(x) == 0) NA_real_ else as.numeric(x)[1]

# Build the one context object every test reads from.
.gof_context <- function(object, predicted_probs = NULL, X = NULL, G = 10) {
  if (inherits(object, "glm")) {
    if (object$family$family != "binomial")
      stop("run.all.gof: the model must be a binomial glm.")
    y  <- as.numeric(object$y)
    ph <- pmin(pmax(as.numeric(stats::fitted(object)), 1e-6), 1 - 1e-6)
    list(y = y, ph = ph, X = stats::model.matrix(object),
         data = tryCatch(stats::model.frame(object), error = function(e) NULL),
         model = object, G = G, n = length(y), has_model = TRUE,
         p = length(stats::coef(object)))
  } else {
    if (!is.numeric(object))
      stop("run.all.gof: 'object' must be a fitted binomial glm or a numeric (0/1) y vector.")
    y <- as.numeric(object)
    if (is.null(predicted_probs))
      stop("run.all.gof: supply 'predicted_probs' when 'object' is not a glm.")
    ph <- pmin(pmax(as.numeric(predicted_probs), 1e-6), 1 - 1e-6)
    if (length(y) != length(ph))
      stop("run.all.gof: 'object' (y) and 'predicted_probs' lengths differ.")
    if (!all(y %in% c(0, 1)))
      stop("run.all.gof: 'object' must be a binary (0/1) vector or a glm.")
    list(y = y, ph = ph, X = X, data = NULL, model = NULL, G = G, n = length(y),
         has_model = FALSE, p = if (is.null(X)) NA_integer_ else ncol(X))
  }
}

# Equal-frequency grouping of the predicted probabilities into G groups.
.gof_groups_ef <- function(ph, G) {
  n <- length(ph)
  pmin(ceiling(rank(ph, ties.method = "first") / (n / G)), G)
}

# Hosmer-Lemeshow statistic for a given grouping; drops empty/degenerate groups.
.gof_hl_stat <- function(y, ph, grp) {
  idx <- split(seq_along(y), grp)
  O  <- vapply(idx, function(I) sum(y[I]),  numeric(1))
  E  <- vapply(idx, function(I) sum(ph[I]), numeric(1))
  ng <- vapply(idx, length, numeric(1))
  keep <- ng > 0 & E > 1e-8 & E < ng - 1e-8
  hl <- sum((O[keep] - E[keep])^2 / (E[keep] * (1 - E[keep] / ng[keep])))
  list(stat = hl, df = sum(keep) - 2)
}

# ---- test wrappers (internal): each returns list(Statistic, df, p_value, Note) ----

gof_pearson <- function(ctx, opts = list()) {
  V  <- ctx$ph * (1 - ctx$ph)
  X2 <- sum((ctx$y - ctx$ph)^2 / V)
  df <- if (ctx$has_model) ctx$model$df.residual else ctx$n - 1
  list(Statistic = X2, df = df, p_value = stats::pchisq(X2, df, lower.tail = FALSE),
       Note = if (ctx$has_model) "" else "df = n-1 (no model)")
}

gof_deviance <- function(ctx, opts = list()) {
  if (ctx$has_model) {
    D <- ctx$model$deviance; df <- ctx$model$df.residual
  } else {
    D  <- -2 * sum(ctx$y * log(ctx$ph) + (1 - ctx$y) * log(1 - ctx$ph)); df <- ctx$n - 1
  }
  list(Statistic = D, df = df, p_value = stats::pchisq(D, df, lower.tail = FALSE),
       Note = if (ctx$has_model) "" else "df = n-1 (no model)")
}

gof_hl <- function(ctx, opts = list()) {
  h <- .gof_hl_stat(ctx$y, ctx$ph, .gof_groups_ef(ctx$ph, ctx$G))
  if (h$df < 1) return(list(Statistic = h$stat, df = h$df, p_value = NA_real_,
                            Note = "too few non-empty groups"))
  list(Statistic = h$stat, df = h$df,
       p_value = stats::pchisq(h$stat, h$df, lower.tail = FALSE), Note = "")
}

gof_hlw <- function(ctx, opts = list()) {
  br <- seq(0, 1, length.out = ctx$G + 1); br[1] <- -Inf; br[length(br)] <- Inf
  grp <- cut(ctx$ph, breaks = br, labels = FALSE)
  h <- .gof_hl_stat(ctx$y, ctx$ph, grp)
  if (h$df < 1) return(list(Statistic = h$stat, df = h$df, p_value = NA_real_,
                            Note = "too few non-empty equal-width bins"))
  list(Statistic = h$stat, df = h$df,
       p_value = stats::pchisq(h$stat, h$df, lower.tail = FALSE), Note = "")
}

gof_ef <- function(ctx, opts = list()) {
  r <- ef.gof(ctx$y, ctx$ph, G = ctx$G)
  list(Statistic = r$Test_Statistic, df = ctx$G - 2, p_value = r$p_value, Note = "")
}

gof_def <- function(ctx, opts = list()) {
  if (!ctx$has_model && is.null(ctx$X))
    return(list(Statistic = NA, df = NA, p_value = NA, Note = "needs a glm model or X"))
  b <- if (is.null(opts$basis)) "poly3" else opts$basis
  r <- if (ctx$has_model) def.gof(ctx$model, G = ctx$G, basis = b)
       else suppressWarnings(def.gof(ctx$y, ctx$ph, X = ctx$X, G = ctx$G, basis = b))
  list(Statistic = r$Test_Statistic, df = r$df, p_value = r$p_value, Note = "")
}

# Stukel (1988) two-direction link test: add sign-split squared-logit terms and
# do a likelihood-ratio test for their joint significance.
gof_stukel <- function(ctx, opts = list()) {
  if (!ctx$has_model)
    return(list(Statistic = NA, df = NA, p_value = NA, Note = "needs a glm model"))
  eta <- as.numeric(stats::predict(ctx$model, type = "link"))
  aug <- cbind(za = 0.5 * eta^2 * (eta >= 0), zb = -0.5 * eta^2 * (eta < 0))
  aug <- aug[, colSums(abs(aug)) > 1e-8, drop = FALSE]
  if (ncol(aug) == 0)
    return(list(Statistic = NA, df = NA, p_value = NA, Note = "Stukel terms degenerate"))
  Xaug <- cbind(stats::model.matrix(ctx$model), aug)
  m1 <- tryCatch(stats::glm.fit(Xaug, ctx$y, family = stats::binomial()),
                 error = function(e) NULL)
  if (is.null(m1))
    return(list(Statistic = NA, df = NA, p_value = NA, Note = "augmented fit failed"))
  LR <- as.numeric(ctx$model$deviance - m1$deviance)
  df <- ncol(aug)
  list(Statistic = LR, df = df, p_value = stats::pchisq(LR, df, lower.tail = FALSE), Note = "")
}

# Registry of the bundled tests (test wrappers are internal, not exported).
.GOF_REGISTRY <- list(
  "Pearson"       = list(fn = gof_pearson,  family = "Global",       needs_model = FALSE, slow = FALSE),
  "Deviance"      = list(fn = gof_deviance, family = "Global",       needs_model = FALSE, slow = FALSE),
  "HL"            = list(fn = gof_hl,       family = "Partition",    needs_model = FALSE, slow = FALSE),
  "HL-equalwidth" = list(fn = gof_hlw,      family = "Partition",    needs_model = FALSE, slow = FALSE),
  "EF"            = list(fn = gof_ef,       family = "Standardized", needs_model = FALSE, slow = FALSE),
  "DEF.poly2"     = list(fn = function(ctx, opts) gof_def(ctx, list(basis = "poly2")),  family = "Directed", needs_model = TRUE, slow = FALSE),
  "DEF.poly3"     = list(fn = function(ctx, opts) gof_def(ctx, list(basis = "poly3")),  family = "Directed", needs_model = TRUE, slow = FALSE),
  "DEF.stukel"    = list(fn = function(ctx, opts) gof_def(ctx, list(basis = "stukel")), family = "Directed", needs_model = TRUE, slow = FALSE),
  "Stukel"        = list(fn = gof_stukel,   family = "Directed",     needs_model = TRUE,  slow = FALSE)
)
