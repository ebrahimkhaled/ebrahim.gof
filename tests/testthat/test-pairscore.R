## pairscore.gof(): the pair score test of the squares and product of the DeepGOF-1 pair.

## The Liu et al. benchmark settings of the DeepGOF-1 paper, with their archived seeds, so the
## published p-values can be reproduced exactly.
liu_data <- function(setting, n, kind, ci, rep) {
  C <- c(0, 0.25, 0.5, 0.75, 1)[ci + (ci == 0)]
  RNGkind("Mersenne-Twister", "Inversion", "Rejection")
  set.seed(101000L + setting * 500000L + kind * 50000L + ci * 5000L + rep)
  if (setting == 2L) {
    x1 <- runif(n, -3, 3); x2 <- rnorm(n, 0, 1.5); x3 <- rt(n, 4)
    d <- data.frame(x1, x2, x3)
    eta <- -2 + 0.3 * x1 + 0.3 * x2 + 0.3 * x3 + C * x1^2
  } else if (setting == 3L) {
    x1 <- runif(n, -3, 3); x2 <- rbinom(n, 1, 0.5); d <- data.frame(x1, x2)
    eta <- -1.972 + 0.135 * x1 + 1.117 * x2 + C * x1 * x2
  } else {
    x1 <- runif(n, -3, 3); x2 <- runif(n, -3, 3); d <- data.frame(x1, x2)
    eta <- 0.3 * x1 + 0.3 * x2 + C * x1 * x2
  }
  d$y <- rbinom(n, 1, plogis(eta))
  fit <- suppressWarnings(glm(stats::reformulate(setdiff(names(d), "y"), "y"),
                              data = d, family = binomial()))
  ## the seed of the bootstrap
  set.seed(880000L + setting * 100000L + n * 10L + kind * 7L + ci + rep * 13L)
  list(fit = fit, events = sum(d$y))
}

test_that("pairscore.gof() reproduces the published score-test p-values", {
  skip_on_cran()
  ## setting, n, kind, ci, rep, events (fingerprint), DeepGOF-1 p, score bootstrap p, chi-squared p
  ref <- read.table(header = TRUE, text = "
    S   n kind ci  rep ev   p_dg p_score            p_asym df
    2 100    2  3    3 40  0.540   0.670 0.64085409348051103  3
    2 100    2  3    4 41  0.620   0.730 0.67460166432445334  3
    3 100    0  0 1001 18  0.360   0.720 0.71486957094914638  2
    3 100    0  0 1003 18  0.195   0.065 0.04803165488938043  2
    3  50    2  4    1 16  0.400   0.060 0.07412779286991578  2
    4  50    2  5    2 25  0.075   0.005 0.00304096416829261  3
  ")
  for (i in seq_len(nrow(ref))) {
    r <- ref[i, ]
    d <- liu_data(r$S, r$n, r$kind, r$ci, r$rep)
    expect_equal(d$events, r$ev)
    st <- d$fit
    seed <- .Random.seed
    ps <- pairscore.gof(st)
    expect_equal(ps$p.value, r$p_score, tolerance = 1e-12)
    expect_equal(ps$p.asymptotic, r$p_asym, tolerance = 1e-8)
    expect_equal(ps$df, r$df)
    ## the same refits give DeepGOF-1 its published p-value, and score.test changes nothing else
    assign(".Random.seed", seed, envir = globalenv())
    dg <- deepgof1(st, score.test = TRUE)
    expect_equal(dg$p.value, r$p_dg, tolerance = 1e-12)
    expect_equal(dg$score.test$p.value, ps$p.value)
    assign(".Random.seed", seed, envir = globalenv())
    plain <- deepgof1(st)
    expect_equal(plain$boot, dg$boot)
    expect_null(plain$score.test)
  }
})

test_that("pairscore.gof() holds its level on a correctly specified model", {
  skip_on_cran()
  set.seed(2026)
  p <- vapply(1:200, function(i) {
    n <- 100
    x1 <- runif(n, -3, 3); x2 <- rnorm(n); x3 <- rnorm(n)
    y <- rbinom(n, 1, plogis(-0.5 + 0.6 * x1 - 0.5 * x2 + 0.3 * x3))
    pairscore.gof(glm(y ~ x1 + x2 + x3, family = binomial()), B = 49)$p.value
  }, 0)
  ## 200 data sets: the rejection rate at .05 has a standard error of about 0.015
  expect_lt(abs(mean(p <= 0.05) - 0.05), 0.05)
  expect_true(all(p >= 1 / 50 & p <= 1))
})

test_that("pairscore.gof() leaves out the square of a binary covariate", {
  set.seed(7)
  n <- 120
  x1 <- runif(n, -3, 3); x2 <- rbinom(n, 1, 0.5)
  y <- rbinom(n, 1, plogis(-0.5 + 0.4 * x1 + 0.8 * x2))
  fit <- glm(y ~ x1 + x2, family = binomial())
  r <- pairscore.gof(fit, B = 19)
  expect_s3_class(r, "pairscore.gof")
  expect_equal(r$df, 2)
  expect_identical(sort(r$axes), c("x1", "x2"))
  expect_true(is.finite(r$statistic) && r$statistic >= 0)
  expect_equal(r$p.value * 20, round(r$p.value * 20))
  expect_output(print(r), "pair: x")
})

test_that("pairscore.gof() refuses models it cannot serve", {
  set.seed(8)
  x <- rnorm(60); y <- rbinom(60, 1, plogis(x))
  expect_error(pairscore.gof(glm(y ~ x, family = binomial()), B = 9), "two covariates")
  z <- rnorm(60)
  expect_error(pairscore.gof(glm(y ~ x + z, family = binomial()), B = 9, reading = "columns"),
               "columns")
})
