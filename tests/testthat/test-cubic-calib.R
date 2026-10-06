make_cubic_data <- function(n = 800, seed = 3, truth = function(eta) plogis(eta)) {
  set.seed(seed)
  x <- runif(n, -3, 3); d <- rbinom(n, 1, 0.5)
  y <- rbinom(n, 1, truth(0.6 * x + 0.5 * d))
  data.frame(y = y, x = x, d = d)
}

## the two logistic fits on the raw linear predictor, as the EDGE paper's simulations compute the test
raw_cubic <- function(y, eta) {
  d1 <- glm(y ~ eta, family = binomial())$deviance
  d3 <- glm(y ~ eta + I(eta^2) + I(eta^3), family = binomial())$deviance
  lr <- max(d1 - d3, 0)
  c(stat = lr, p = pchisq(lr, 2, lower.tail = FALSE))
}

test_that("cubic.calib.gof is the deviance drop of adding eta^2 and eta^3, on 2 df", {
  dat <- make_cubic_data()
  fit <- glm(y ~ x + d, family = binomial(), data = dat)
  res <- cubic.calib.gof(fit)
  ref <- raw_cubic(dat$y, fit$linear.predictors)
  expect_s3_class(res, "data.frame")
  expect_identical(res$Test, "Cubic-LR")
  expect_equal(res$df, 2)
  expect_equal(res$Test_Statistic, ref[["stat"]], tolerance = 1e-8)
  expect_equal(res$p_value, ref[["p"]], tolerance = 1e-8)
})

test_that("frozen predictions give the same test as the fitted model", {
  dat <- make_cubic_data(seed = 7)
  fit <- glm(y ~ x + d, family = binomial(), data = dat)
  a <- cubic.calib.gof(fit)
  b <- cubic.calib.gof(dat$y, predicted_probs = fitted(fit))
  expect_equal(b$p_value, a$p_value, tolerance = 1e-6)
})

test_that("the run.all.gof row is the same test", {
  dat <- make_cubic_data(seed = 11)
  fit <- glm(y ~ x + d, family = binomial(), data = dat)
  row <- run.all.gof(fit, tests = "Cubic-LR")
  expect_identical(row$Family, "Calibration")
  expect_equal(row$p_value, cubic.calib.gof(fit)$p_value, tolerance = 1e-12)
})

test_that("a bend in the calibration curve is detected", {
  cll <- function(eta) 1 - exp(-exp(eta))           # complementary log-log truth, fitted by a logit
  dat <- make_cubic_data(n = 4000, seed = 5, truth = cll)
  res <- cubic.calib.gof(glm(y ~ x + d, family = binomial(), data = dat))
  expect_lt(res$p_value, 0.05)
})

test_that("bad input is refused, and a constant prediction gives NA with a warning", {
  dat <- make_cubic_data()
  expect_error(cubic.calib.gof(glm(x ~ d, family = gaussian(), data = dat)), "binomial")
  expect_error(cubic.calib.gof(dat$y), "predicted_probs")
  expect_error(cubic.calib.gof(dat$y, predicted_probs = rep(0.5, 3)), "lengths")
  expect_error(cubic.calib.gof(dat$y + 1, predicted_probs = rep(0.5, nrow(dat))), "binary")
  expect_warning(r <- cubic.calib.gof(dat$y, predicted_probs = rep(0.4, nrow(dat))), "all equal")
  expect_true(is.na(r$p_value))
})
