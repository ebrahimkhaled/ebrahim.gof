make_fit <- function(seed = 1, n = 500) {
  set.seed(seed)
  x <- runif(n, -3, 3)
  glm(rbinom(n, 1, 1 / (1 + exp(-(0.6 * x)))) ~ x, family = binomial())
}

test_that("run.all.gof returns a tidy battery for a model", {
  res <- run.all.gof(make_fit())
  expect_s3_class(res, "data.frame")
  expect_equal(names(res), c("Test", "Family", "Statistic", "df", "p_value", "Note"))
  expect_true(all(c("EF", "DEF.poly3", "HL", "Stukel") %in% res$Test))
  expect_true(any(grepl("Ensemble", res$Test)))
  pv <- res$p_value[is.finite(res$p_value)]
  expect_true(all(pv >= 0 & pv <= 1))
})

test_that("run.all.gof selects a subset", {
  res <- run.all.gof(make_fit(), tests = c("EF", "DEF.poly3", "HL"))
  expect_setequal(res$Test, c("EF", "DEF.poly3", "HL"))
})

test_that("prediction-only input runs the ph-only tests with a message", {
  fit <- make_fit(2)
  expect_message(res <- run.all.gof(as.numeric(fit$y), stats::fitted(fit)),
                 "prediction-based")
  expect_true(all(is.na(res$p_value[res$Test %in% c("DEF.poly2", "Stukel")])))
  expect_true(is.finite(res$p_value[res$Test == "EF"]))
})

test_that("Osius-Rojek and Copas match goflogit (golden values)", {
  # Verified 2026-06-09 against the goflogit() macro (7 tests.R) to ~1e-15.
  set.seed(42)
  n <- 400
  x1 <- rnorm(n); x2 <- runif(n, -2, 2)
  y  <- rbinom(n, 1, plogis(0.3 + 0.8 * x1 - 0.5 * x2))
  fit <- glm(y ~ x1 + x2, family = binomial())
  res <- run.all.gof(fit)
  o <- res[res$Test == "Osius-Rojek", ]
  cps <- res[res$Test == "Copas-RSS", ]
  expect_equal(o$Statistic, 0.985610, tolerance = 1e-4)   # Osius z
  expect_equal(o$p_value,  0.162162, tolerance = 1e-4)
  expect_equal(cps$p_value, 0.466596, tolerance = 1e-4)   # Copas p
})

test_that("covariate-space tests appear and return valid results", {
  set.seed(5)
  n <- 400
  x1 <- rnorm(n); d <- factor(sample(c("A", "B"), n, replace = TRUE))
  y  <- rbinom(n, 1, plogis(0.3 + 0.6 * x1 + ifelse(d == "B", 0.5, 0)))
  fit <- glm(y ~ x1 + d, family = binomial())
  res <- run.all.gof(fit)
  expect_true(all(c("Tsiatis", "Xie", "Pulkstenis-Robinson") %in% res$Test))
  for (tt in c("Tsiatis", "Xie", "Pulkstenis-Robinson")) {
    p <- res$p_value[res$Test == tt]
    expect_true(is.na(p) || (p >= 0 && p <= 1))
  }
})

test_that("Pulkstenis-Robinson matches its source (deterministic golden)", {
  # PR uses no clustering, so this is platform-stable. Verified vs PR_test_only.R.
  set.seed(99)
  n <- 500
  x1 <- rnorm(n); x2 <- runif(n, -2, 2); d <- factor(sample(c("A", "B", "C"), n, replace = TRUE))
  eta <- 0.2 + 0.7 * x1 - 0.5 * x2 + ifelse(d == "B", 0.6, ifelse(d == "C", -0.4, 0))
  y <- rbinom(n, 1, plogis(eta))
  fit <- glm(y ~ x1 + x2 + d, family = binomial())
  pr <- run.all.gof(fit)
  pr <- pr[pr$Test == "Pulkstenis-Robinson", ]
  expect_equal(pr$Statistic, 1.910442, tolerance = 1e-4)
  expect_equal(pr$p_value,  0.591201, tolerance = 1e-4)
})

test_that("Pulkstenis-Robinson returns NA when no categorical covariate", {
  set.seed(3)
  n <- 300; x1 <- rnorm(n); x2 <- rnorm(n)
  fit <- glm(rbinom(n, 1, plogis(0.5 * x1)) ~ x1 + x2, family = binomial())
  res <- run.all.gof(fit)
  expect_true(is.na(res$p_value[res$Test == "Pulkstenis-Robinson"]))
})

test_that("run.all.gof flags a cloglog link misfit", {
  set.seed(7)
  x <- runif(1500, -3, 3)
  fit <- glm(rbinom(1500, 1, 1 - exp(-exp(0.6 * x))) ~ x, family = binomial())
  res <- run.all.gof(fit)
  expect_lt(res$p_value[res$Test == "EF"], 0.05)
  expect_lt(res$p_value[res$Test == "Stukel"], 0.05)
})
