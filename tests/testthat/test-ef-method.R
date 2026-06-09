test_that("ef.gof chisq (default) and normal give the same statistic, different p", {
  set.seed(1)
  n <- 400
  x <- rnorm(n)
  y  <- rbinom(n, 1, plogis(0.5 * x))
  ph <- plogis(0.5 * x)
  rc <- ef.gof(y, ph, method = "chisq")
  rn <- ef.gof(y, ph, method = "normal")
  expect_equal(rc$Test_Statistic, rn$Test_Statistic)        # same standardized statistic
  expect_false(isTRUE(all.equal(rc$p_value, rn$p_value)))   # different reference -> different p
  expect_true(rc$p_value >= 0 && rc$p_value <= 1)
})

test_that("ef.gof accepts a fitted glm as its first argument", {
  set.seed(2)
  n <- 400
  x <- rnorm(n)
  fit <- glm(rbinom(n, 1, plogis(0.5 * x)) ~ x, family = binomial())
  expect_equal(ef.gof(fit)$p_value,
               ef.gof(as.numeric(fit$y), stats::fitted(fit))$p_value)
})
