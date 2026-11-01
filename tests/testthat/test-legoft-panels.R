## legoft(): panels and rules of the panel paper.

make_fit <- function(n = 200, seed = 11) {
  set.seed(seed)
  d <- data.frame(x1 = runif(n, -3, 3), x2 = rnorm(n))
  d$y <- rbinom(n, 1, plogis(0.3 + 0.8 * d$x1 - 0.5 * d$x2))
  glm(y ~ x1 + x2, family = binomial(), data = d)
}

test_that("the default is still the frozen LEGofT rule", {
  skip_on_cran()
  fit <- make_fit()
  a <- legoft(fit, B = 25, seed = 3)
  b <- legoft(fit, B = 25, seed = 3, weights = ebrahim.gof:::.LEGOFT_WEIGHTS)
  expect_equal(a$p.value, b$p.value)
  expect_equal(a$statistic, b$statistic)
  expect_identical(a$panel, "legoft")
  expect_match(a$method, "frozen weights")
})

test_that("panels give equal weights over the right tests", {
  skip_on_cran()
  fit <- make_fit()
  r <- legoft(fit, B = 25, seed = 3, panel = "protected")
  expect_setequal(names(r$weights), c("HL", "EF", "Pigeon-Heyse", "le-Cessie"))
  expect_equal(unname(r$weights), rep(0.25, 4))
  expect_true(all(is.finite(r$members)))
  expect_true(r$p.value > 0 && r$p.value <= 1)
  expect_length(ebrahim.gof:::.LEGOFT_PANELS$full, 14)
})

test_that("Stukel and Osius-Rojek can be named in weights", {
  skip_on_cran()
  fit <- make_fit()
  w <- c(HL = 1, Stukel = 1, "Osius-Rojek" = 1) / 3
  r <- legoft(fit, B = 25, seed = 3, weights = w)
  expect_identical(r$panel, "custom")
  expect_true(all(is.finite(r$members)))
})

test_that("rules other than cauchy need equal weights", {
  fit <- make_fit()
  expect_error(legoft(fit, B = 25, rule = "fisher"), "equal weights")
  expect_error(legoft(fit, B = 25, weights = c(HL = 1, Bogus = 1)), "unknown: Bogus")
})

test_that("each rule returns a valid bootstrap p-value", {
  skip_on_cran()
  fit <- make_fit()
  for (rl in c("cauchy", "hmp", "fisher", "minp")) {
    r <- legoft(fit, B = 25, seed = 3, panel = "equal", rule = rl)
    expect_identical(r$rule, rl)
    expect_true(r$p.value >= 1 / 26 && r$p.value <= 1)
  }
})
