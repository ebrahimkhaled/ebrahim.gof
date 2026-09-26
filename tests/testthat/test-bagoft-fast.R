# bagoft.fast() against the BAGofT package: same seed, identical output (identical(), not
# all.equal()). Reduced nsplits / nsim keep the tests quick; the code path is the same.

bagoft_pkg <- function(formula, data, nsplits, nsim) {
  r <- NULL
  invisible(utils::capture.output(suppressMessages(suppressWarnings(
    r <- BAGofT::BAGofT(BAGofT::testGlmBi(formula = formula, link = "logit"),
                        BAGofT::parRF(), data = data, nsplits = nsplits, nsim = nsim)))))
  r
}

expect_same_bagoft <- function(fast, pkg) {
  for (k in c("p.value", "p.value2", "p.value3", "pmean", "pmedian", "pmin"))
    expect_identical(unname(fast[[k]]), unname(pkg[[k]]), label = k)
  expect_identical(unname(fast$simRes), unname(pkg$simRes))
}

test_that("bagoft.fast reproduces BAGofT exactly: two covariates", {
  skip_if_not_installed("BAGofT")
  skip_if_not_installed("randomForest")
  set.seed(101)
  n <- 60
  d <- data.frame(x1 = rnorm(n), x2 = runif(n, -2, 2))
  d$y <- rbinom(n, 1, plogis(0.3 + 0.8 * d$x1 - 0.5 * d$x2^2))
  d <- d[c("y", "x1", "x2")]
  fit <- glm(y ~ x1 + x2, family = binomial(), data = d)
  set.seed(7); f <- bagoft.fast(fit, nsplits = 5, nsim = 5)
  set.seed(7); p <- bagoft_pkg(y ~ x1 + x2, d, nsplits = 5, nsim = 5)
  expect_s3_class(f, "bagoft_fast")
  expect_same_bagoft(f, p)
})

test_that("bagoft.fast reproduces BAGofT exactly: a factor and a column the model omits", {
  skip_if_not_installed("BAGofT")
  skip_if_not_installed("randomForest")
  set.seed(202)
  n <- 70
  d <- data.frame(y = 0, x1 = rnorm(n), g = factor(sample(c("a", "b", "c"), n, TRUE)),
                  z = rnorm(n))
  d$y <- rbinom(n, 1, plogis(0.5 * d$x1 + (d$g == "b") - 0.7 * d$z))
  fit <- glm(y ~ x1 + g, family = binomial(), data = d)   # z is left out of the model
  set.seed(8); f <- bagoft.fast(fit, data = d, nsplits = 5, nsim = 5)
  set.seed(8); p <- bagoft_pkg(y ~ x1 + g, d, nsplits = 5, nsim = 5)
  expect_same_bagoft(f, p)
})

test_that("bagoft.fast reproduces BAGofT exactly: seven covariates (dcPre pre-selection)", {
  skip_if_not_installed("BAGofT")
  skip_if_not_installed("randomForest")
  skip_if_not_installed("dcov")
  set.seed(303)
  n <- 80
  X <- matrix(rnorm(n * 7), n, 7, dimnames = list(NULL, paste0("x", 1:7)))
  d <- data.frame(y = rbinom(n, 1, plogis(drop(X %*% c(0.6, -0.4, 0.3, 0, 0, 0.2, 0)))), X)
  fml <- y ~ x1 + x2 + x3 + x4 + x5 + x6 + x7
  fit <- glm(fml, family = binomial(), data = d)
  set.seed(9); f <- suppressWarnings(bagoft.fast(fit, nsplits = 4, nsim = 4))
  set.seed(9); p <- bagoft_pkg(fml, d, nsplits = 4, nsim = 4)
  expect_true(f$settings$preselected)
  expect_same_bagoft(f, p)
})

test_that("the BAGofT row gives the same p-value with either engine", {
  skip_if_not_installed("BAGofT")
  skip_if_not_installed("randomForest")
  set.seed(404)
  n <- 60
  x <- rnorm(n)
  y <- rbinom(n, 1, plogis(0.2 + 0.9 * x))
  fit <- glm(y ~ x, family = binomial())            # one predictor: the constant column is added
  ctl <- function(engine) list(BAGofT = list(nsim = 4, nsplits = 4, engine = engine))
  set.seed(3); rf <- run.all.gof(fit, tests = "BAGofT", control = ctl("fast"))
  set.seed(3); rp <- run.all.gof(fit, tests = "BAGofT", control = ctl("package"))
  expect_identical(rf$p_value, rp$p_value)
  expect_match(rf$Note, "fast engine")
  expect_match(rp$Note, "BAGofT package")
  expect_false(is.na(rf$p_value))
})

test_that("bagoft.fast checks its input", {
  skip_if_not_installed("randomForest")
  set.seed(1)
  x <- rnorm(100); y <- rbinom(100, 1, 0.5)
  fit <- glm(y ~ x, family = binomial())
  expect_error(bagoft.fast(lm(y ~ x)), "binomial glm")
  expect_error(bagoft.fast(fit, data = data.frame(x = x)), "not a column")
  r <- suppressWarnings(bagoft.fast(fit, nsplits = 3, nsim = 0))   # statistics only; one covariate works
  expect_null(r$p.value)
  expect_true(r$pmean >= 0 && r$pmean <= 1)
})
