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
  expect_true(any(grepl("groups named (closed testing) at FWER 0.05: LINK + COV", out, fixed = TRUE)))
  expect_true(any(grepl("recommended update: revise the model", out)))
  expect_true(any(grepl("^INTERCEPT +0\\.1100 +0\\.1100", out)))
  capture.output(expect_invisible(print(r)))
})

test_that("naming = 'holm' and 'bonferroni' reproduce the paper's script", {
  d <- loc_data()
  r <- localize.external(d$y, d$p, d$X, M = 99, naming = "holm", seed = 7)
  ## the intersections are unchanged; only the naming rule differs
  expect_equal(r$intersection, localize.external(d$y, d$p, d$X, M = 99, seed = 7)$intersection)
  expect_equal(r$adjusted, c(INTERCEPT = 0.18, SLOPE = 0.18, LINK = 0.12, COV = 0.04), tolerance = 1e-12)
  expect_identical(r$named, "COV")
  expect_identical(r$naming, "holm")
  rb <- localize.external(d$y, d$p, d$X, M = 99, naming = "bonferroni", seed = 7)
  expect_equal(rb$adjusted, c(INTERCEPT = 0.44, SLOPE = 0.36, LINK = 0.16, COV = 0.04), tolerance = 1e-12)
  expect_identical(rb$named, "COV")
  out <- capture.output(print(r))
  expect_true(any(grepl("groups named (Holm over the groups) at FWER 0.05: COV", out, fixed = TRUE)))
  ## in-sample: Holm over LINK (0.10) and COV (0.05) names nothing at B = 19
  fit <- glm(d$y ~ x1 + x2 + x3, data = d$X, family = binomial())
  ri <- localize.gof(fit, X = d$X, B = 19, naming = "holm", seed = 8)
  expect_identical(ri$named, character(0))
  expect_equal(ri$adjusted, c(LINK = 0.1, COV = 0.1), tolerance = 1e-12)
})

test_that("robust = TRUE reproduces the paper's script", {
  d <- loc_data()
  r <- localize.external(d$y, d$p, d$X, M = 99, robust = TRUE, seed = 7)
  expect_equal(unname(r$intersection[c("INTERCEPT", "SLOPE", "LINK", "COV", "LINK+COV")]),
               c(0.11, 0.07, 0.04, 0.01, 0.02), tolerance = 1e-12)
  expect_equal(r$members, c(
    INTERCEPT.cal = 0.122438318746109, SLOPE.slope = 0.0721015566648641,
    LINK.poly = 0.00985599609578549, LINK.stukel = 0.0134600239638279,
    LINK.spline = 0.035546393704911, COV.poly = 0.00025010706706531,
    COV.products = 0.000165031810536696, COV.spline = 0.000857523843869752), tolerance = 1e-12)
  expect_identical(r$named, c("LINK", "COV"))
  expect_true(r$robust)

  fit <- glm(d$y ~ x1 + x2 + x3, data = d$X, family = binomial())
  ri <- localize.gof(fit, X = d$X, B = 19, dealias = TRUE, robust = TRUE, seed = 8)
  expect_equal(ri$intersection, c(LINK = 0.75, COV = 0.05, `LINK+COV` = 0.05), tolerance = 1e-12)
  expect_equal(ri$members, c(
    LINK.poly = 0.483764091120531, LINK.stukel = 0.654156527825062,
    LINK.spline = 0.453869341724405, COV.poly = 0.00318839357151969,
    COV.products = 0.000732266376956144, COV.spline = 0.00223743917446776), tolerance = 1e-12)
  expect_identical(ri$dealiased, c("x1", "x2"))
  expect_true(any(grepl("normal scores", capture.output(print(ri)))))
})

test_that("robust = TRUE limits the influence of a single extreme record", {
  ## a correct model; one covariate value is then recorded 50 times too large
  set.seed(8)
  n <- 400
  X <- data.frame(x1 = rnorm(n), x2 = rnorm(n))
  p <- plogis(-0.3 + 0.8 * X$x1 + 0.5 * X$x2)
  y <- rbinom(n, 1, p)
  i <- which.max(abs(y - p))
  X2 <- X; X2$x1[i] <- X2$x1[i] * 50
  cov_p <- function(XX, rb) localize.external(y, p, XX, M = 199, robust = rb, seed = 1)$single[["COV"]]
  ## the default bases let the one record drive COV (0.40 -> 0.005 on this seed) ...
  expect_gt(cov_p(X, FALSE), 0.2)
  expect_lt(cov_p(X2, FALSE), 0.05)
  ## ... the normal-score bases barely move (0.41 -> 0.32)
  expect_gt(cov_p(X2, TRUE), 0.2)
  expect_lt(abs(cov_p(X2, TRUE) - cov_p(X, TRUE)), 0.15)
})

test_that("the new options are checked", {
  d <- loc_data()
  expect_error(localize.external(d$y, d$p, d$X, M = 19, naming = "hochberg"), "should be one of")
  expect_error(localize.external(d$y, d$p, d$X, M = 19, robust = NA), "TRUE or FALSE")
})

test_that("plot() draws a verdict in both settings and returns it invisibly", {
  d <- loc_data()
  re <- localize.external(d$y, d$p, d$X, M = 49, seed = 1)
  fit <- stats::glm(y ~ x1 + x2 + x3, family = binomial(), data = data.frame(y = d$y, d$X))
  ri <- localize.gof(fit, B = 19, seed = 1)
  f <- tempfile(fileext = ".pdf"); grDevices::pdf(f)
  on.exit({ grDevices::dev.off(); unlink(f) })
  expect_identical(withVisible(plot(re))$visible, FALSE)
  expect_identical(plot(re), re)
  expect_silent(plot(re, which = "compass", colour = FALSE))
  expect_silent(plot(ri, which = "lattice"))
  expect_silent(plot(ri))
  expect_error(plot(re, which = "radar"), "should be one of")
})

test_that("plot = TRUE draws the verdict as the result is computed", {
  d <- loc_data()
  f <- tempfile(fileext = ".pdf"); grDevices::pdf(f)
  on.exit({ grDevices::dev.off(); unlink(f) })
  before <- grDevices::recordPlot()
  r <- localize.external(d$y, d$p, d$X, M = 49, seed = 1, plot = TRUE)
  expect_s3_class(r, "gof_localize")
  expect_identical(r$named, localize.external(d$y, d$p, d$X, M = 49, seed = 1, plot = FALSE)$named)
  expect_error(localize.external(d$y, d$p, d$X, M = 19, plot = NA), "TRUE or FALSE")
  expect_false(interactive())   # so the default does not draw under R CMD check
})

