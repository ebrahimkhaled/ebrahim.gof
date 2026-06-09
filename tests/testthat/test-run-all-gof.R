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

test_that("run.all.gof flags a cloglog link misfit", {
  set.seed(7)
  x <- runif(1500, -3, 3)
  fit <- glm(rbinom(1500, 1, 1 - exp(-exp(0.6 * x))) ~ x, family = binomial())
  res <- run.all.gof(fit)
  expect_lt(res$p_value[res$Test == "EF"], 0.05)
  expect_lt(res$p_value[res$Test == "Stukel"], 0.05)
})
