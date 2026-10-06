test_that("the belt reports exactly the p-value def.gof() returns", {
  set.seed(11)
  n <- 600
  x <- runif(n, -3, 3); d <- rbinom(n, 1, 0.5)
  y <- rbinom(n, 1, plogis(0.6 * x + 0.5 * d))
  fit <- glm(y ~ x + d, family = binomial())

  for (bs in c("poly3", "poly2", "stukel", "sym")) {
    for (wt in c("unit", "score")) {
      b <- edge.belt(fit, G = 20, basis = bs, weights = wt, nsim = 2000, seed = 1)
      g <- def.gof(fit, G = 20, basis = bs, weights = wt)
      expect_identical(b$test$p_value, g$p_value)
      expect_identical(b$test$Test_Statistic, g$Test_Statistic)
    }
  }
})

test_that("the drawn direction is the component the unit statistic squares", {
  set.seed(12)
  n <- 700
  x <- runif(n, -3, 3)
  y <- rbinom(n, 1, 1 - exp(-exp(0.8 * x)))         # log-log truth, fitted as logit
  fit <- glm(y ~ x, family = binomial())

  for (bs in c("poly2", "poly3", "stukel")) {
    b <- edge.belt(fit, G = 25, basis = bs, weights = "unit", nsim = 2000, seed = 2)
    # S = r'P r = ||P r||^2, and P r is what the belt draws
    expect_equal(sum(b$groups$direction^2), b$test$Test_Statistic, tolerance = 1e-10)
  }
})

test_that("the grouped quantities reproduce the fit", {
  set.seed(13)
  n <- 500
  x <- runif(n, -2, 2)
  y <- rbinom(n, 1, plogis(0.7 * x))
  fit <- glm(y ~ x, family = binomial())
  b <- edge.belt(fit, G = 20, nsim = 2000, seed = 3)

  expect_equal(sum(b$groups$observed), sum(y))
  expect_equal(sum(b$groups$expected), sum(fitted(fit)), tolerance = 1e-8)
  expect_equal(sum(b$groups$size), n)
  # the estimation adjustment shrinks every null SD below one
  expect_true(all(b$groups$sd_null < 1))
  # residual and its standardization agree
  expect_equal(b$groups$z, b$groups$residual / b$groups$sd_null, tolerance = 1e-12)
})

test_that("the simultaneous band is wider than the pointwise one and is honest under the null", {
  set.seed(14)
  n <- 1000
  x <- runif(n, -3, 3)
  y <- rbinom(n, 1, plogis(0.5 * x))
  fit <- glm(y ~ x, family = binomial())
  b <- edge.belt(fit, G = 40, nsim = 20000, seed = 4)

  expect_gt(b$critical$simultaneous, b$critical$pointwise)
  # a 40-group pointwise band is expected to flag about two groups by chance; the
  # simultaneous band is built so that flagging any group at all is a 5% event
  expect_lte(sum(b$groups$outside_simultaneous), 2)
})

test_that("without a design matrix the bands fall back to Omega = I", {
  set.seed(15)
  n <- 400
  x <- runif(n, -2, 2)
  y <- rbinom(n, 1, plogis(0.7 * x))
  p <- as.numeric(fitted(glm(y ~ x, family = binomial())))
  b <- suppressWarnings(edge.belt(y, predicted_probs = p, G = 15, nsim = 2000, seed = 5))

  expect_false(b$adjusted)
  expect_equal(b$groups$sd_null, rep(1, 15), tolerance = 1e-12)
})

test_that("print and plot run without error", {
  set.seed(16)
  n <- 500
  x <- runif(n, -3, 3)
  y <- rbinom(n, 1, 1 - exp(-exp(0.8 * x)))
  fit <- glm(y ~ x, family = binomial())
  b <- edge.belt(fit, G = 20, nsim = 2000, seed = 6)

  expect_output(print(b), "EDGE-Belt")
  pf <- tempfile(fileext = ".png")
  png(pf, width = 1000, height = 500)
  expect_silent(plot(b))
  dev.off()
  expect_true(file.exists(pf))
  unlink(pf)
})

test_that("external = TRUE treats Omega = I as exact rather than conservative", {
  set.seed(17)
  n <- 1200
  x <- runif(n, -2, 2)
  dev <- seq_len(n / 2); val <- setdiff(seq_len(n), dev)
  y <- rbinom(n, 1, plogis(0.7 * x))
  f <- glm(y[dev] ~ x[dev], family = binomial())
  p <- as.numeric(plogis(cbind(1, x[val]) %*% coef(f)))      # frozen on the held-out half

  expect_warning(edge.belt(y[val], predicted_probs = p, G = 10, nsim = 2000, seed = 7),
                 "conservative")
  b <- edge.belt(y[val], predicted_probs = p, G = 10, external = TRUE, nsim = 2000, seed = 7)
  expect_true(b$external)
  expect_false(b$adjusted)
  expect_equal(b$groups$sd_null, rep(1, 10), tolerance = 1e-12)
  expect_output(print(b), "exact")
  # the belt carries the external test itself (constant column, d + 1 df), and the curve is
  # the component that statistic squares
  ext <- def.gof(y[val], predicted_probs = p, G = 10, external = TRUE)
  expect_identical(b$test$p_value, ext$p_value)
  expect_identical(b$test$df, ext$df)
  expect_equal(sum(b$groups$direction^2), ext$Test_Statistic, tolerance = 1e-10)
  # a fitted glm cannot be an external validation of itself
  expect_error(edge.belt(f, external = TRUE), "frozen predictions")
})
