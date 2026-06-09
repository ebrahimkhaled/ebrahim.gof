make_fit <- function(seed = 1, n = 600) {
  set.seed(seed)
  x <- runif(n, -3, 3)
  glm(rbinom(n, 1, 1 / (1 + exp(-(0.6 * x)))) ~ x, family = binomial())
}

test_that("def.ensemble.gof returns the documented shape", {
  res <- def.ensemble.gof(make_fit())
  expect_equal(names(res), c("Test", "Combiner", "Components", "k", "p_value"))
  expect_identical(res$Combiner, "cct")
  expect_identical(res$Components, "poly2+poly3+stukel")
  expect_equal(res$k, 3L)
  expect_true(res$p_value >= 0 && res$p_value <= 1)
})

test_that("add_ef appends the EF component", {
  res <- def.ensemble.gof(make_fit(), add_ef = TRUE)
  expect_equal(res$k, 4L)
  expect_match(res$Components, "EF$")
})

test_that("basis='ensemble' delegates to def.ensemble.gof", {
  fit <- make_fit(3)
  expect_equal(def.gof(fit, basis = "ensemble")$p_value,
               def.ensemble.gof(fit)$p_value, tolerance = 1e-12)
})

test_that("combiner internals are correct", {
  f <- ebrahim.gof:::.combine_pvalues
  expect_equal(f(c(0.5, 0.5, 0.5), "cct"), 0.5, tolerance = 1e-8)
  expect_equal(f(c(0.04, 0.3, 0.5), "minp"), 1 - (1 - 0.04)^3, tolerance = 1e-10)
})

test_that("extra_pvalues are included", {
  res <- def.ensemble.gof(make_fit(), extra_pvalues = c(Tsiatis = 0.2))
  expect_equal(res$k, 4L)
  expect_match(res$Components, "Tsiatis")
})
