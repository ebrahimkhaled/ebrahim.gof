# The parallel (PSOCK) path of run.all.gof: same seed + same ncores must give
# identical bootstrap p-values across runs (clusterSetRNGStream reproducibility).

test_that("parallel Stute-Zhu is reproducible across two same-seed runs", {
  skip_on_cran()

  set.seed(11)
  n <- 120
  x <- stats::rnorm(n)
  y <- stats::rbinom(n, 1, stats::plogis(0.5 * x))
  fit <- stats::glm(y ~ x, family = stats::binomial())
  ctl <- list("Stute-Zhu" = list(B = 25))

  set.seed(42)
  r1 <- suppressMessages(
    run.all.gof(fit, tests = "Stute-Zhu", include_slow = TRUE,
                parallel = TRUE, ncores = 2, install = "no", control = ctl))
  set.seed(42)
  r2 <- suppressMessages(
    run.all.gof(fit, tests = "Stute-Zhu", include_slow = TRUE,
                parallel = TRUE, ncores = 2, install = "no", control = ctl))

  p1 <- r1$p_value[r1$Test == "Stute-Zhu"]
  p2 <- r2$p_value[r2$Test == "Stute-Zhu"]
  expect_true(is.finite(p1))
  expect_identical(p1, p2)
  expect_match(r1$Note[r1$Test == "Stute-Zhu"], "parallel, 2 workers")
})

test_that("ncores = 1 (or parallel = FALSE) stays on the sequential path", {
  skip_on_cran()

  set.seed(3)
  n <- 100
  x <- stats::rnorm(n)
  y <- stats::rbinom(n, 1, stats::plogis(0.4 * x))
  fit <- stats::glm(y ~ x, family = stats::binomial())
  ctl <- list("Stute-Zhu" = list(B = 10))

  set.seed(5)
  rs <- suppressMessages(
    run.all.gof(fit, tests = "Stute-Zhu", include_slow = TRUE,
                parallel = TRUE, ncores = 1, install = "no", control = ctl))
  expect_false(grepl("parallel", rs$Note[rs$Test == "Stute-Zhu"]))

  # and the sequential path is unchanged by the new arguments
  set.seed(5)
  r0 <- suppressMessages(
    run.all.gof(fit, tests = "Stute-Zhu", include_slow = TRUE,
                parallel = FALSE, install = "no", control = ctl))
  expect_identical(rs$p_value[rs$Test == "Stute-Zhu"],
                   r0$p_value[r0$Test == "Stute-Zhu"])
})
