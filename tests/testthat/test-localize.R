## localize.external() and localize.gof(): the closed-testing localization of misfit.
## The expected values below were computed with the paper's script localize_groups.R (Option A),
## on the same data and seed; the package must reproduce them to the last digit.

loc_data <- function() {
  set.seed(101)
  n <- 300
  X <- data.frame(x1 = rnorm(n), x2 = runif(n, -2, 2), x3 = rbinom(n, 1, 0.4))
  p <- plogis(-0.4 + 0.7 * X$x1 - 0.5 * X$x2 + 0.6 * X$x3)
  y <- rbinom(n, 1, plogis(qlogis(p) + 0.5 * X$x1 * X$x2 + 0.3))
  list(X = X, p = p, y = y)
}

test_that("localize.external() reproduces the paper's script: Monte Carlo reference", {
  d <- loc_data()
  r <- localize.external(d$y, d$p, d$X, M = 99, seed = 7)
  expect_s3_class(r, "gof_localize")
  expect_equal(r$intersection, c(
    INTERCEPT = 0.11, SLOPE = 0.09, LINK = 0.04, COV = 0.01, `INTERCEPT+SLOPE` = 0.1,
    `INTERCEPT+LINK` = 0.04, `INTERCEPT+COV` = 0.01, `SLOPE+LINK` = 0.05,
    `SLOPE+COV` = 0.01, `LINK+COV` = 0.01, `INTERCEPT+SLOPE+LINK` = 0.05,
    `INTERCEPT+SLOPE+COV` = 0.01, `INTERCEPT+LINK+COV` = 0.01, `SLOPE+LINK+COV` = 0.01,
    `INTERCEPT+SLOPE+LINK+COV` = 0.01), tolerance = 1e-12)
  expect_equal(r$members, c(
    INTERCEPT.cal = 0.122438318746109, SLOPE.slope = 0.105483399536703,
    LINK.poly = 0.0121229475675447, LINK.stukel = 0.0133652296584094,
    LINK.spline = 0.0332869149652293, COV.poly = 0.000111253743386382,
    COV.products = 1.54818167051849e-05, COV.spline = 0.000522833046577262), tolerance = 1e-12)
  expect_identical(r$named, c("LINK", "COV"))
  expect_identical(r$highest, "COV")
  expect_match(r$action, "revise the model")
  ## a group is named exactly when its adjusted p-value is at most alpha
  expect_identical(names(r$adjusted)[r$adjusted <= 0.05], r$named)
  expect_equal(r$draws, 99L)
})

test_that("localize.external() reproduces the paper's script: multiplier reference", {
  d <- loc_data()
  r <- localize.external(d$y, d$p, d$X, M = 99, calibration = "multiplier", seed = 7)
  expect_equal(unname(r$intersection[c("INTERCEPT", "SLOPE", "LINK", "COV", "LINK+COV")]),
               c(0.17, 0.16, 0.13, 0.01, 0.01), tolerance = 1e-12)
  expect_equal(r$members, c(
    INTERCEPT.cal = 0.136258452948913, SLOPE.slope = 0.15100511395836,
    LINK.poly = 0.0828404667401494, LINK.stukel = 0.0820980378268524,
    LINK.spline = 0.175734974852949, COV.poly = 0.00385180352570479,
    COV.products = 0.000150410313360771, COV.spline = 0.0106483683076818), tolerance = 1e-12)
  expect_identical(r$named, "COV")
  expect_identical(r$calibration, "multiplier")
})

test_that("localize.gof() reproduces the paper's script, with and without de-aliasing", {
  d <- loc_data()
  fit <- glm(d$y ~ x1 + x2 + x3, data = d$X, family = binomial())
  r0 <- localize.gof(fit, X = d$X, B = 19, seed = 8)
  expect_equal(r0$intersection, c(LINK = 0.1, COV = 0.05, `LINK+COV` = 0.05), tolerance = 1e-12)
  expect_equal(r0$members, c(
    LINK.poly = 0.0176989662681755, LINK.stukel = 0.0215221018414533,
    LINK.spline = 0.051660202108351, COV.poly = 0.000746762913450828,
    COV.products = 0.000150723637027739, COV.spline = 0.000638515389597895), tolerance = 1e-12)
  expect_identical(r0$named, "COV")
  expect_equal(r0$draws, 19L)

  r1 <- localize.gof(fit, X = d$X, B = 19, dealias = TRUE, seed = 8)
  expect_equal(r1$intersection, c(LINK = 0.5, COV = 0.05, `LINK+COV` = 0.05), tolerance = 1e-12)
  expect_equal(unname(r1$members[c("LINK.poly", "LINK.stukel", "LINK.spline")]),
               c(0.333624448256956, 0.520578114302126, 0.39011437346929), tolerance = 1e-12)
  expect_identical(r1$dealiased, c("x1", "x2"))

  ## X defaults to the numeric covariates of the model frame
  expect_equal(localize.gof(fit, B = 19, seed = 8)$intersection, r0$intersection)
})

test_that("on data with no misfit nothing is named", {
  set.seed(31)
  n <- 500
  X <- data.frame(x1 = rnorm(n), x2 = rnorm(n))
  p <- plogis(-0.5 + 0.8 * X$x1 + 0.6 * X$x2)
  y <- rbinom(n, 1, p)
  r <- localize.external(y, p, X, M = 199, seed = 1)
  expect_identical(r$named, character(0))
  expect_true(is.na(r$highest))
  expect_match(r$action, "keep the model")
  expect_true(all(r$intersection > 0 & r$intersection <= 1))
  expect_equal(r$intersection * 200, round(r$intersection * 200))   # on the 1/(M+1) grid
})

test_that("an omitted interaction names COV in-sample", {
  set.seed(2)
  n  <- 500
  x1 <- rnorm(n); x2 <- rnorm(n)
  y  <- rbinom(n, 1, plogis(-0.3 + 0.8 * x1 + 0.6 * x2 + 0.8 * x1 * x2))
  fit <- glm(y ~ x1 + x2, family = binomial())
  r <- localize.gof(fit, B = 49, seed = 1)
  expect_true("COV" %in% r$named)
  expect_identical(r$highest, "COV")
})

test_that("one covariate drops the COV group", {
  set.seed(41)
  n <- 300
  x <- rnorm(n)
  p <- plogis(0.2 + 0.9 * x)
  y <- rbinom(n, 1, p)
  r <- localize.external(y, p, data.frame(x = x), M = 49, seed = 1)
  expect_identical(names(r$single), c("INTERCEPT", "SLOPE", "LINK"))
  expect_false(any(grepl("COV", names(r$intersection))))
  fit <- glm(y ~ x, family = binomial())
  r2 <- localize.gof(fit, B = 19, seed = 1)
  expect_identical(names(r2$single), "LINK")
  expect_identical(names(r2$groups), "LINK")
})

test_that("predictions of exactly 0 or 1 are kept just inside (0, 1)", {
  set.seed(51)
  n <- 300
  X <- data.frame(x1 = rnorm(n), x2 = rnorm(n))
  p <- plogis(-0.2 + 0.7 * X$x1 + 0.5 * X$x2)
  y <- rbinom(n, 1, p)
  p[1:3] <- 0; y[1:3] <- 0
  p[4:6] <- 1; y[4:6] <- 1
  r <- localize.external(y, p, X, M = 49, seed = 1)
  expect_true(all(is.finite(r$members)))
  expect_true(all(r$intersection > 0 & r$intersection <= 1))
})

test_that("localize.gof() refuses a probit fit and points to external validation", {
  set.seed(61)
  x1 <- rnorm(200); x2 <- rnorm(200)
  y <- rbinom(200, 1, pnorm(0.2 + 0.6 * x1))
  fit <- glm(y ~ x1 + x2, family = binomial(link = "probit"))
  expect_error(localize.gof(fit, B = 19), "localize.external")
  expect_error(localize.gof(lm(y ~ x1)), "binomial")
})

test_that("the inputs are checked", {
  d <- loc_data()
  expect_error(localize.external(d$y + 1, d$p, d$X, M = 19), "0/1")
  expect_error(localize.external(d$y, d$p[-1], d$X, M = 19), "one probability")
  expect_error(localize.external(d$y, d$p + 1, d$X, M = 19), "one probability")
  expect_error(localize.external(d$y, d$p, d$X[-1, ], M = 19), "one row per outcome")
  expect_error(localize.external(d$y, d$p, cbind(d$X, k = 1), M = 19), "constant column")
  expect_error(localize.external(d$y, d$p, data.frame(f = factor(d$X$x3), x1 = d$X$x1), M = 19),
               "numeric")
  expect_error(localize.external(d$y, d$p, d$X, M = 5), "at least 19")
  expect_error(localize.external(d$y, d$p, d$X, M = 19, alpha = 1), "alpha")
})

test_that("print() shows the groups named and the recommended update", {
  d <- loc_data()
  r <- localize.external(d$y, d$p, d$X, M = 99, seed = 7)
  out <- capture.output(print(r))
  expect_true(any(grepl("groups named at FWER 0.05: LINK \\+ COV", out)))
  expect_true(any(grepl("recommended update: revise the model", out)))
  expect_true(any(grepl("^INTERCEPT +0\\.1100 +0\\.1100", out)))
  capture.output(expect_invisible(print(r)))
})
