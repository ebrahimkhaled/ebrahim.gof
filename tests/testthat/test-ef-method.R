test_that("ef.gof chisq (default) refers EF to chi-square(G - 2); method = 'normal' is the old standardization", {
  set.seed(1)
  n <- 400
  x <- rnorm(n)
  y  <- rbinom(n, 1, plogis(0.5 * x))
  ph <- plogis(0.5 * x)
  rc <- ef.gof(y, ph, reference = "chisq")
  expect_equal(rc$p_value, pchisq(rc$Test_Statistic, 8, lower.tail = FALSE))
  expect_identical(ef.gof(y, ph, method = "chisq"), rc)
  expect_warning(rn <- ef.gof(y, ph, method = "normal"), "versions up to 2.9.0")
  expect_equal(rn$Test_Statistic, (rc$Test_Statistic - 8) / 4)
  expect_equal(rn$p_value, pnorm(rn$Test_Statistic, lower.tail = FALSE))
  expect_true(rc$p_value >= 0 && rc$p_value <= 1)
})

test_that("ef.gof accepts a fitted glm as its first argument", {
  set.seed(2)
  n <- 400
  x <- rnorm(n)
  fit <- glm(rbinom(n, 1, plogis(0.5 * x)) ~ x, family = binomial())
  expect_equal(ef.gof(fit)$p_value,
               ef.gof(as.numeric(fit$y), stats::fitted(fit))$p_value)
  expect_equal(ef.gof(fit, reference = "normal")$p_value,
               ef.gof(as.numeric(fit$y), stats::fitted(fit), reference = "normal",
                      X = model.matrix(fit))$p_value)
})
