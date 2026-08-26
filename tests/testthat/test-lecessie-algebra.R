## The 2.5.0 rewrite of gof_lecessie() replaced two O(N^3) matrix products by exact
## algebraic identities (a rank-p expansion of (I-H)'R(I-H), and elementwise evaluation
## of 2 tr(MVMV), which is how le Cessie and van Houwelingen state the variance in
## eq. A.6 before collapsing it to a trace). Nothing about the statistic, the reference
## distribution or the p-value was allowed to move. This file pins that: the reference
## below is the pre-rewrite body, verbatim, and the shipped function must agree with it
## to 1e-12 relative.

## Pre-rewrite body: forms H explicitly, then takes a full matrix product to extract a
## trace. Kept here only as the thing the fast path has to reproduce.
lecessie_pre_rewrite <- function(y, ph, X, covs) {
  N <- length(y)
  resids <- y - ph
  dlist <- lapply(covs, function(x) {
    if (is.numeric(x)) {
      as.numeric(0.5 * (stats::dist(scale(x)))^2)
    } else {
      xx <- as.numeric(as.factor(x)); nc <- length(unique(xx))
      if (nc <= 1) rep(0, N * (N - 1) / 2)
      else as.numeric((stats::dist(xx, method = "manhattan") != 0) * nc / (nc - 1))
    }
  })
  dist.mat <- matrix(0, N, N)
  dist.mat[lower.tri(dist.mat)] <- sqrt(rowSums(as.data.frame(dlist)))
  dist.mat <- dist.mat + t(dist.mat)
  R.raw <- pmax(1 - dist.mat / mean(dist.mat), 0)
  Q.raw <- sum(as.numeric(resids %*% R.raw) * resids)

  mu2   <- ph * (1 - ph)
  hat   <- (mu2 * X) %*% solve(crossprod(X, mu2 * X)) %*% t(X)
  IH    <- diag(N) - hat
  R.cor <- crossprod(IH, R.raw %*% IH)                 # (I-H)' R (I-H)
  E.Q   <- sum(diag(R.cor) * mu2)
  mu4   <- mu2 * (1 - 3 * mu2)
  VarQ1 <- sum(diag(R.cor)^2 * (mu4 - 3 * mu2^2))
  R.tmp <- R.cor * rep(mu2, each = N)
  VarQ2 <- 2 * sum(diag(R.tmp %*% R.tmp))              # the O(N^3) trace
  VarQ  <- VarQ1 + VarQ2
  if (!is.finite(VarQ) || VarQ <= 0) return(c(NA_real_, NA_real_, NA_real_))
  Test <- Q.raw * 2 * E.Q / VarQ
  df   <- 2 * E.Q^2 / VarQ
  c(Test, df, stats::pchisq(Test, df, lower.tail = FALSE))
}

## Six data-generating processes: the null, two alternatives (so agreement is not an
## artefact of a well-fitting model), a factor covariate under both hypotheses, and a
## design whose fitted probabilities sit near 0 and 1 -- the worst case for the rank-p
## expansion, because that is where X'VX is closest to singular.
lecessie_case <- function(scenario, n, seed) {
  set.seed(seed)
  x1 <- runif(n, -3, 3); x2 <- rnorm(n)
  g  <- factor(sample(c("a", "b", "c"), n, replace = TRUE))
  eta <- switch(scenario,
    null      = 0.3 + 0.8 * x1 - 0.5 * x2,
    quadratic = 0.3 + 0.8 * x1 - 0.5 * x2 + 1.2 * x1^2,
    interact  = 0.3 + 0.8 * x1 - 0.5 * x2 + 1.5 * x1 * x2,
    factor_ok = 0.3 + 0.8 * x1 + 1.6 * (g == "c"),
    factor_lo = 0.3 + 0.8 * x1 + 1.6 * (g == "c") + 1.0 * x1^2,
    extreme   = -1.0 + 3.5 * x1 - 2.5 * x2)
  dat <- data.frame(y = rbinom(n, 1, plogis(eta)), x1 = x1, x2 = x2, g = g)
  fml <- if (grepl("^factor", scenario)) y ~ x1 + g else y ~ x1 + x2
  suppressWarnings(glm(fml, data = dat, family = binomial()))
}

test_that("the rewritten le Cessie moments equal the pre-rewrite ones to 1e-12", {
  scenarios <- c("null", "quadratic", "interact", "factor_ok", "factor_lo", "extreme")
  for (scenario in scenarios) {
    for (seed in 1:3) {
      fit <- lecessie_case(scenario, n = 120, seed = seed)
      ## Take y, ph, X and the covariates from the context the shipped function is
      ## given, so the two paths differ only in the algebra. ph is clamped in
      ## .gof_context(); reproducing that by hand would be a second chance to be wrong.
      ctx <- ebrahim.gof:::.gof_context(fit)
      old <- lecessie_pre_rewrite(ctx$y, ctx$ph, ctx$X, ctx$data[, -1, drop = FALSE])
      got <- ebrahim.gof:::gof_lecessie(ctx)
      new <- c(got$Statistic, got$df, got$p_value)
      expect_equal(new, old, tolerance = 1e-12,
                   info = paste(scenario, "seed", seed))
    }
  }
})

test_that("the le Cessie moment reference keeps the transpose that smwrStats drops", {
  ## H = VX(X'VX)^{-1}X' is not symmetric in a weighted fit, so (I-H)'R(I-H) and the
  ## smwrStats form (I-H)R(I-H) are different matrices and give different p-values.
  ## The fast path symmetrises its result; that is only valid for the transposed form,
  ## so pin the asymmetry of H and the disagreement of the two sandwiches.
  fit <- lecessie_case("null", n = 60, seed = 11)
  ctx <- ebrahim.gof:::.gof_context(fit)
  X   <- ctx$X; mu2 <- ctx$ph * (1 - ctx$ph); N <- length(mu2)
  H   <- (mu2 * X) %*% solve(crossprod(X, mu2 * X)) %*% t(X)
  expect_false(isTRUE(all.equal(H, t(H), tolerance = 1e-8)))

  IH <- diag(N) - H
  set.seed(3)
  Z <- matrix(rnorm(N * N), N); R <- crossprod(Z)          # any symmetric R
  M_correct <- crossprod(IH, R %*% IH)                     # (I-H)' R (I-H)
  M_smwr    <- IH %*% R %*% IH                             # (I-H)  R (I-H)
  expect_equal(M_correct, t(M_correct), tolerance = 1e-8)  # symmetric, as the algebra says
  expect_false(isTRUE(all.equal(M_correct, M_smwr, tolerance = 1e-8)))
})
