## localize.update() and localize.size(): the expected values were computed with the paper's scripts
## localize_update.R and localize_size.R, on the same data and seed.

upd_data <- function() {
  set.seed(101)
  n <- 300
  X <- data.frame(x1 = rnorm(n), x2 = runif(n, -2, 2), x3 = rbinom(n, 1, 0.4))
  p <- plogis(-0.4 + 0.7 * X$x1 - 0.5 * X$x2 + 0.6 * X$x3)
  y <- rbinom(n, 1, plogis(qlogis(p) + 0.5 * X$x1 * X$x2 + 0.3))
  list(X = X, p = p, y = y)
}

test_that("localize.update() reproduces the paper's script", {
  skip_on_cran()   # exact Monte Carlo counts: one can differ on another CPU
  d <- upd_data()
  r <- localize.update(d$y, d$p, d$X, B = 49, seed = 7)
  expect_equal(r$single, c(INTERCEPT = 0.14, SLOPE = 0.14, LINK = 0.02, COV = 0.02), tolerance = 1e-8)
  expect_equal(r$adjusted, c(INTERCEPT = 0.28, SLOPE = 0.28, LINK = 0.08, COV = 0.08), tolerance = 1e-8)
  expect_identical(r$named, character(0))
})

test_that("localize.size() reproduces the paper's script", {
  skip_on_cran()
  d <- upd_data()
  r <- localize.size(d$y, d$p, d$X, R = 10, seed = 8)
  expect_equal(r$parts$size, c(4.088865591, 3.217268295, 5.07106397, 12.27784344), tolerance = 1e-7)
  expect_equal(r$parts$chance, c(1.753319354, 1.103811582, 3.739805029, 7.063109025), tolerance = 1e-7)
  expect_equal(r$ici, 7.457518011, tolerance = 1e-7)
})

test_that("localize.update() returns a verdict that print() and plot() can read", {
  d <- upd_data()
  r <- localize.update(d$y, d$p, d$X, B = 19, seed = 1)
  expect_s3_class(r, "gof_localize")
  expect_identical(names(r$single), c("INTERCEPT", "SLOPE", "LINK", "COV"))
  expect_true("COV.main" %in% names(r$members))
  expect_length(r$intersection, 15L)
  expect_identical(names(r$adjusted)[r$adjusted <= r$alpha], r$named)
  expect_output(print(r), "updated-model reference")
  pdf(NULL); on.exit(dev.off())
  expect_silent(plot(r))
})

test_that("localize.update() options behave", {
  d <- upd_data()
  r <- localize.update(d$y, d$p, d$X, B = 19, main = FALSE, seed = 1)
  expect_false("COV.main" %in% names(r$members))
  cl <- rep(1:6, length.out = length(d$y))
  r <- localize.update(d$y, d$p, d$X, B = 19, cluster = cl, tau = 0.3, seed = 1)
  expect_match(r$calibration, "centre effect")
  r <- localize.update(d$y, d$p, d$X[, 1, drop = FALSE], B = 19, seed = 1)
  expect_false("COV" %in% names(r$single))   # one covariate: no covariate structure to test
  expect_error(localize.update(d$y, d$p, d$X, B = 19, tau = 0.3), "cluster")
  expect_error(localize.update(d$y, d$p, d$X, B = 19, tau = -1), "tau")
  expect_error(localize.update(d$y, d$p, d$X, B = 5), "at least 19")
})

test_that("the seed leaves the session's random stream alone", {
  d <- upd_data()
  set.seed(3); a <- runif(1); set.seed(3)
  localize.update(d$y, d$p, d$X, B = 19, seed = 1)
  localize.size(d$y, d$p, d$X, R = 2, seed = 1)
  expect_identical(runif(1), a)
})

test_that("localize.size() returns the four parts and prints", {
  d <- upd_data()
  r <- localize.size(d$y, d$p, d$X, R = 3, seed = 1)
  expect_s3_class(r, "gof_localize_size")
  expect_identical(r$parts$part, c("INTERCEPT", "SLOPE", "LINK", "COV"))
  expect_identical(dimnames(r$calibration)[[2]], c("estimate", "lower", "upper"))
  expect_output(print(r), "Size of each part")
  expect_error(localize.size(d$y, d$p, d$X, R = 1), "at least 2")
})
