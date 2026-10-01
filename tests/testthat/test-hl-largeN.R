test_that("HL-largeN reproduces Nattino, Pennell and Lemeshow's (2020) application", {
  ## their Section 5: n = 315,828, G = 10; C = 25.35 gives p = 0.010 and C = 16.66 gives p = 0.11
  a <- ebrahim.gof:::.hl_largeN_p(25.35, 315828, 8)
  b <- ebrahim.gof:::.hl_largeN_p(16.66, 315828, 8)
  expect_equal(round(a$p, 3), 0.010)
  expect_equal(round(b$p, 2), 0.11)
  expect_equal(round(a$eps0 * 1e3, 2), 2.74)          # their eps0 for G = 10: 2.74e-3
})

test_that("HL-largeN runs in the battery and is never below the ordinary test's p-value", {
  set.seed(31)
  n <- 3000
  x <- rnorm(n)
  y <- rbinom(n, 1, plogis(-1 + x + 0.15 * x^2))
  fit <- glm(y ~ x, family = binomial())
  r <- run.all.gof(fit, tests = c("HL", "HL-largeN"), install = "no")
  expect_true(all(c("HL", "HL-largeN") %in% r$Test))
  expect_equal(r$Statistic[r$Test == "HL"], r$Statistic[r$Test == "HL-largeN"])
  expect_gte(r$p_value[r$Test == "HL-largeN"], r$p_value[r$Test == "HL"])
  expect_match(r$Note[r$Test == "HL-largeN"], "eps_hat")
})
