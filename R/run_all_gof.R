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
#' The currently bundled tests are: \code{Pearson}, \code{Deviance},
#' \code{Osius-Rojek}, and \code{Copas-RSS} (global / standardized);
#' \code{HL} (Hosmer-Lemeshow deciles) and \code{HL-equalwidth} (partition);
#' \code{EF} (the omnibus Ebrahim-Farrington test); \code{DEF.poly2/poly3/stukel}
#' and \code{Stukel} (directed); \code{Tsiatis}, \code{Xie}, and
#' \code{Pulkstenis-Robinson} (covariate-space); plus the two ensemble rows
#' (\code{Ensemble.Vote(3DEF)} and \code{Ensemble.Univ(3DEF+EF)}) from the Cauchy
#' combination test. Further thesis tests (McCullagh, information matrix,
#' le Cessie, and opt-in bootstrap/GAM) are planned for a later build.
#'
#' Notes: \code{Tsiatis} and \code{Xie} cluster the covariate space with k-means
#' (a fixed internal seed, so results are reproducible and the caller's RNG is
#' left untouched). \code{Xie} uses the corrected degrees of freedom
#' \eqn{G - k/2 - 1} with \eqn{k} the number of predictors. \code{Pulkstenis-Robinson}
#' auto-detects the categorical covariate (any factor/character/logical, or a
#' numeric with at most \code{getOption("ebrahim.gof.pr.maxlev", 6)} distinct
#' values); it returns \code{NA} with a note when none is present.
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
#' @importFrom stats fitted predict model.matrix model.frame coef deviance pchisq binomial glm.fit kmeans median
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

# Osius-Rojek normal-approximation correction of the Pearson statistic.
# Ported from goflogit (7 tests.R); xpxi = (X'WX)^{-1}, W = diag(p(1-p)).
gof_osius <- function(ctx, opts = list()) {
  if (is.null(ctx$X))
    return(list(Statistic = NA, df = NA, p_value = NA, Note = "needs the design matrix X"))
  X <- ctx$X; ph <- ctx$ph; y <- ctx$y; n <- length(ph)
  V <- ph * (1 - ph)
  xpxi <- tryCatch(solve(crossprod(X, V * X)), error = function(e) NULL)
  if (is.null(xpxi))
    return(list(Statistic = NA, df = NA, p_value = NA, Note = "singular information matrix"))
  pearson <- sum((y - ph)^2 / V)
  vosius  <- sum(1 / V - 4)                       # binary (m = 1): 2 + 1/V - 6
  cvec    <- crossprod(X, 1 - 2 * ph)             # X'(1 - 2p)
  qosius  <- as.numeric(t(cvec) %*% xpxi %*% cvec)
  varosius <- vosius - qosius
  if (!is.finite(varosius) || varosius <= 0)
    return(list(Statistic = NA, df = NA, p_value = NA, Note = "non-positive variance"))
  z <- (pearson - n) / sqrt(varosius)
  # one-sided upper, matching the validated goflogit implementation
  list(Statistic = z, df = NA_real_, p_value = stats::pnorm(z, lower.tail = FALSE), Note = "")
}

# Copas (1989) unweighted residual-sum-of-squares test. Binary expansion is
# trivial (one trial per observation). Ported from goflogit (7 tests.R).
gof_copas <- function(ctx, opts = list()) {
  if (is.null(ctx$X))
    return(list(Statistic = NA, df = NA, p_value = NA, Note = "needs the design matrix X"))
  X <- ctx$X; ph <- ctx$ph; y <- ctx$y
  V <- ph * (1 - ph)
  copas    <- sum((y - ph)^2)
  meacopas <- sum(V)
  W   <- diag(V)
  c12 <- 1 - 2 * ph
  M   <- tryCatch(diag(length(ph)) - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X),
                  error = function(e) NULL)
  if (is.null(M))
    return(list(Statistic = NA, df = NA, p_value = NA, Note = "singular information matrix"))
  varcopas <- as.numeric(t(c12) %*% M %*% W %*% c12)
  if (!is.finite(varcopas) || varcopas <= 0)
    return(list(Statistic = NA, df = NA, p_value = NA, Note = "non-positive variance"))
  z <- (copas - meacopas) / sqrt(varcopas)
  list(Statistic = z, df = 1, p_value = stats::pchisq(z^2, 1, lower.tail = FALSE), Note = "")
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

# Moore-Penrose pseudo-inverse via SVD (matches MASS::ginv; avoids a dependency).
.gof_ginv <- function(M, tol = sqrt(.Machine$double.eps)) {
  s <- svd(M)
  pos <- s$d > max(tol * s$d[1], 0)
  if (!any(pos)) return(matrix(0, ncol(M), nrow(M)))
  s$v[, pos, drop = FALSE] %*% (t(s$u[, pos, drop = FALSE]) / s$d[pos])
}

# k-means with a fixed seed (123) that does not disturb the caller's RNG state.
.gof_kmeans <- function(mat, centers, nstart = 1L) {
  has <- exists(".Random.seed", envir = .GlobalEnv)
  old <- if (has) get(".Random.seed", envir = .GlobalEnv) else NULL
  set.seed(123)
  cl <- stats::kmeans(mat, centers = centers, nstart = nstart)$cluster
  if (has) assign(".Random.seed", old, envir = .GlobalEnv)
  cl
}

# Tsiatis (1980) clustering score test: cluster the covariate space, then score-
# test the cluster indicators added to the model. Ported from Tsiatis.R.
gof_tsiatis <- function(ctx, opts = list()) {
  if (is.null(ctx$X))
    return(list(Statistic = NA, df = NA, p_value = NA, Note = "needs the design matrix X"))
  X <- ctx$X; ph <- ctx$ph; y <- ctx$y
  cov_mat <- X[, -1, drop = FALSE]                       # drop intercept column
  if (ncol(cov_mat) < 1)
    return(list(Statistic = NA, df = NA, p_value = NA, Note = "no covariates to cluster"))
  cl <- tryCatch(.gof_kmeans(cov_mat, ctx$G), error = function(e) NULL)
  if (is.null(cl))
    return(list(Statistic = NA, df = NA, p_value = NA, Note = "clustering failed"))
  Xc <- stats::model.matrix(~ factor(cl))[, -1, drop = FALSE]
  if (ncol(Xc) < 1)
    return(list(Statistic = NA, df = NA, p_value = NA, Note = "only one non-empty cluster"))
  W   <- ph * (1 - ph)
  U   <- colSums(Xc * (y - ph))
  V11 <- crossprod(X, W * X)
  V12 <- crossprod(X, W * Xc)
  V22 <- crossprod(Xc, W * Xc)
  V   <- V22 - t(V12) %*% .gof_ginv(V11) %*% V12
  Tstat <- as.numeric(t(U) %*% .gof_ginv(V) %*% U)
  rankV <- sum(eigen(V, symmetric = TRUE, only.values = TRUE)$values > 1e-8)
  list(Statistic = Tstat, df = rankV,
       p_value = stats::pchisq(Tstat, rankV, lower.tail = FALSE), Note = "")
}

# Xie covariate-space grouped chi-square (own group rule, fractional df). From Xie.R.
gof_xie <- function(ctx, opts = list()) {
  if (is.null(ctx$X))
    return(list(Statistic = NA, df = NA, p_value = NA, Note = "needs the design matrix X"))
  X <- ctx$X; ph <- ctx$ph; y <- ctx$y
  k <- ncol(X) - 1
  cov_mat <- X[, -1, drop = FALSE]
  if (ncol(cov_mat) < 1)
    return(list(Statistic = NA, df = NA, p_value = NA, Note = "no covariates to cluster"))
  G <- if (k < 5) 10 else k + 5
  cl <- tryCatch(.gof_kmeans(cov_mat, G, nstart = 25L), error = function(e) NULL)
  if (is.null(cl))
    return(list(Statistic = NA, df = NA, p_value = NA, Note = "clustering failed"))
  stat <- 0
  for (I in split(seq_along(y), cl)) {
    ng <- length(I); pbar <- mean(ph[I])
    if (pbar > 0 && pbar < 1)
      stat <- stat + (sum(y[I]) - ng * pbar)^2 / (ng * pbar * (1 - pbar))
  }
  df <- G - k / 2 - 1
  if (df <= 0) return(list(Statistic = stat, df = df, p_value = NA, Note = "df <= 0"))
  list(Statistic = stat, df = df, p_value = stats::pchisq(stat, df, lower.tail = FALSE), Note = "")
}

# Pulkstenis-Robinson: covariate patterns from categorical vars, split by median
# fitted prob, chi-square on the 2*M subgroups. Base-R port of PR_test_only.R.
gof_pr <- function(ctx, opts = list()) {
  if (!ctx$has_model || is.null(ctx$data))
    return(list(Statistic = NA, df = NA, p_value = NA,
                Note = "needs a glm model (for categorical covariates)"))
  # model.frame stores the response in column 1; the rest are the model terms.
  cand <- names(ctx$data)[-1]
  maxlev <- getOption("ebrahim.gof.pr.maxlev", 6)
  cat_vars <- if (!is.null(opts$cat_var)) opts$cat_var else
    cand[vapply(cand, function(v) {
      col <- ctx$data[[v]]
      is.factor(col) || is.character(col) || is.logical(col) ||
        (is.numeric(col) && length(unique(col)) <= maxlev)
    }, logical(1))]
  if (length(cat_vars) == 0)
    return(list(Statistic = NA, df = NA, p_value = NA, Note = "no categorical covariate"))

  y <- ctx$y; ph <- ctx$ph
  patt <- do.call(paste, c(lapply(cat_vars, function(v) as.character(ctx$data[[v]])), sep = "_"))
  M <- length(unique(patt))
  lev <- character(length(y))
  for (pp in unique(patt)) {
    ix <- which(patt == pp)
    lev[ix] <- ifelse(ph[ix] <= stats::median(ph[ix]), "low", "high")
  }
  idx <- split(seq_along(y), paste(patt, lev, sep = "::"))
  os <- vapply(idx, function(I) sum(y[I] == 1), numeric(1))
  of <- vapply(idx, function(I) sum(y[I] == 0), numeric(1))
  es <- vapply(idx, function(I) sum(ph[I]), numeric(1))
  ef <- vapply(idx, function(I) sum(1 - ph[I]), numeric(1))
  keep <- es > 0 & ef > 0
  chisq <- sum((os[keep] - es[keep])^2 / es[keep] + (of[keep] - ef[keep])^2 / ef[keep])
  df <- 2 * M - length(cat_vars) - 2
  if (df <= 0)
    return(list(Statistic = chisq, df = df, p_value = NA, Note = "df <= 0 (too few patterns)"))
  list(Statistic = chisq, df = df, p_value = stats::pchisq(chisq, df, lower.tail = FALSE),
       Note = paste0("cat: ", paste(cat_vars, collapse = ",")))
}

# Registry of the bundled tests (test wrappers are internal, not exported).
.GOF_REGISTRY <- list(
  "Pearson"       = list(fn = gof_pearson,  family = "Global",       needs_model = FALSE, slow = FALSE),
  "Deviance"      = list(fn = gof_deviance, family = "Global",       needs_model = FALSE, slow = FALSE),
  "Osius-Rojek"   = list(fn = gof_osius,    family = "Standardized", needs_model = TRUE,  slow = FALSE),
  "Copas-RSS"     = list(fn = gof_copas,    family = "Standardized", needs_model = TRUE,  slow = FALSE),
  "HL"            = list(fn = gof_hl,       family = "Partition",    needs_model = FALSE, slow = FALSE),
  "HL-equalwidth" = list(fn = gof_hlw,      family = "Partition",    needs_model = FALSE, slow = FALSE),
  "EF"            = list(fn = gof_ef,       family = "Standardized", needs_model = FALSE, slow = FALSE),
  "DEF.poly2"     = list(fn = function(ctx, opts) gof_def(ctx, list(basis = "poly2")),  family = "Directed", needs_model = TRUE, slow = FALSE),
  "DEF.poly3"     = list(fn = function(ctx, opts) gof_def(ctx, list(basis = "poly3")),  family = "Directed", needs_model = TRUE, slow = FALSE),
  "DEF.stukel"    = list(fn = function(ctx, opts) gof_def(ctx, list(basis = "stukel")), family = "Directed", needs_model = TRUE, slow = FALSE),
  "Stukel"        = list(fn = gof_stukel,   family = "Directed",     needs_model = TRUE,  slow = FALSE),
  "Tsiatis"             = list(fn = gof_tsiatis, family = "Covariate-space", needs_model = TRUE, slow = FALSE),
  "Xie"                 = list(fn = gof_xie,     family = "Covariate-space", needs_model = TRUE, slow = FALSE),
  "Pulkstenis-Robinson" = list(fn = gof_pr,      family = "Covariate-space", needs_model = TRUE, slow = FALSE)
)
