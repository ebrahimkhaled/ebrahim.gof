# cavity.gof(): the leave-one-out Hosmer-Lemeshow test (Ebrahim and Arashi, 2026)

cav_data <- function(n = 300, p = 60, seed = 20261101) {
  set.seed(seed)
  X <- matrix(stats::rnorm(n * p), n)
  y <- stats::rbinom(n, 1, 1 / (1 + exp(-drop(X %*% stats::rnorm(p, sd = 0.3)))))
  list(X = X, y = y)
}

test_that("cavity.gof reproduces the research engine of the paper", {
  d <- cav_data()
  r <- cavity.gof(d$X, d$y, lambda = 50, smooth.degree = 3)
  # values from cav_gof_fast() in the paper's archive (doi:10.5281/zenodo.23861454) on the same data
  expect_equal(r$statistic, 9.8431754905, tolerance = 1e-6)
  expect_equal(r$p.value, 0.2762059886, tolerance = 1e-6)
  expect_equal(r$p.exact, 0.2782102140, tolerance = 1e-5)
  expect_equal(r$p.smooth, 0.5678074821, tolerance = 1e-5)
  expect_s3_class(r, "cavity.gof")
  expect_identical(r$df, 8)
  expect_equal(r$twin.share, 0)
})

test_that("group table is consistent with the logistic recalibration", {
  d <- cav_data()
  r <- cavity.gof(d$X, d$y, lambda = 50, check.twins = FALSE)
  expect_equal(sum(r$table$observed), sum(d$y))
  expect_equal(sum(r$table$expected), sum(d$y), tolerance = 1e-6)   # score equation of the intercept
  expect_equal(sum(r$table$residual^2), r$statistic, tolerance = 1e-10)
  expect_equal(nrow(r$table), 10L)
})

test_that("the cavity index is close to the exact leave-one-out linear predictor", {
  d <- cav_data(n = 200, p = 40)
  r <- cavity.gof(d$X, d$y, lambda = 40, check.twins = FALSE)
  X1 <- cbind(1, d$X); D <- diag(c(0, rep(1, ncol(d$X))))
  for (i in c(3, 77, 150)) {
    b <- ebrahim.gof:::.cav_ridge(X1[-i, ], d$y[-i], 40, D)
    expect_equal(r$index[i], sum(X1[i, ] * b), tolerance = 0.02)
  }
})

test_that("lambda on the glmnet scale gives the same test", {
  d <- cav_data()
  a <- cavity.gof(d$X, d$y, lambda = 50, check.twins = FALSE)
  b <- cavity.gof(d$X, d$y, lambda = 50 / 300, lambda_scale = "glmnet", check.twins = FALSE)
  expect_equal(a$statistic, b$statistic, tolerance = 1e-8)
})

test_that("repeated rows trigger the twin warning", {
  d <- cav_data()
  Xd <- d$X[c(1:150, 1:150), ]
  set.seed(2)
  yd <- stats::rbinom(300, 1, 0.4)
  expect_warning(r <- cavity.gof(Xd, yd, lambda = 50), "twin")
  expect_gt(r$twin.share, 0.9)
})

test_that("split and simulated references return valid p-values", {
  d <- cav_data(n = 200, p = 30)
  set.seed(3)
  r <- cavity.gof(d$X, d$y, lambda = 30, split = TRUE, reference = "simulate", B = 19, check.twins = FALSE)
  for (p in c(r$p.split, r$p.crossfit, r$p.simulated)) expect_true(p >= 0 && p <= 1)
  expect_identical(r$B, 19L)
})

test_that("the logistic-normal recalibration runs", {
  d <- cav_data(n = 400, p = 40)
  r <- cavity.gof(d$X, d$y, lambda = 40, recalibration = "logistic-normal", check.twins = FALSE)
  expect_true(r$p.value >= 0 && r$p.value <= 1)
  expect_identical(r$recal.family, "logistic-normal")
})

test_that("a singular recalibration information falls back to the pseudo-inverse", {
  A <- matrix(c(2, 1, 1, 1, 1, 0.5, 1, 0.5, 0.5), 3)   # rank 2
  Ai <- ebrahim.gof:::.cav_inverse(A)
  expect_equal(A %*% Ai %*% A, A, tolerance = 1e-8)
  expect_equal(ebrahim.gof:::.cav_inverse(diag(2) * 3), diag(2) / 3)
})

test_that("the logistic-normal recalibration does not stop under a strong misfit", {
  skip_on_cran()
  # seed 76 made 3.0.0 stop with "system is computationally singular"
  set.seed(76)
  n <- 400; p <- 100
  S <- 0.7^abs(outer(1:p, 1:p, "-"))
  X <- matrix(rnorm(n * p), n, p) %*% chol(S)
  eta <- drop(X %*% (rep(c(0.5, -0.5), length.out = p) / sqrt(p) * 3))
  eta <- eta + 0.6 * (eta^3 - 3 * eta)
  y <- rbinom(n, 1, plogis(eta))
  r <- cavity.gof(X, y, lambda = 80, recalibration = "logistic-normal", check.twins = FALSE)
  expect_true(is.finite(r$p.value) && is.finite(r$p.exact))
})

test_that("the lasso path runs on the active set", {
  skip_if_not_installed("glmnet")
  d <- cav_data()
  r <- cavity.gof(d$X, d$y, lambda = 0.02, alpha = 1, lambda_scale = "glmnet")
  expect_gt(r$nnz, 0)
  expect_true(r$p.value >= 0 && r$p.value <= 1)
})

test_that("bad input is refused", {
  d <- cav_data(n = 50, p = 5)
  expect_error(cavity.gof(d$X, d$y + 1, 10), "binary")
  expect_error(cavity.gof(d$X, d$y, -1), "positive")
  expect_error(cavity.gof(d$X[-1, ], d$y, 10), "one row per response")
  expect_error(cavity.gof(d$X, d$y, 10, alpha = 2), "alpha")
})
