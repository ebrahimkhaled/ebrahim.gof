## A fit with na.action = na.exclude pads fitted(), residuals() and predict() with NA for the
## records it dropped; every test must read the unpadded fitted values, so that the result is
## the one of the same fit made with na.omit.

na_fits <- function() {
  set.seed(42)
  n <- 600
  d <- data.frame(x1 = rnorm(n), x2 = runif(n, -2, 2))
  d$y <- rbinom(n, 1, plogis(-0.3 + 0.7 * d$x1 - 0.5 * d$x2))
  d$x1[c(5, 77, 301)] <- NA
  list(ex = glm(y ~ x1 + x2, family = binomial(), data = d, na.action = na.exclude),
       om = glm(y ~ x1 + x2, family = binomial(), data = d, na.action = na.omit))
}

test_that("the fast battery is the same under na.exclude and na.omit", {
  f <- na_fits()
  expect_length(fitted(f$ex), 600)            # the padding the battery must not read
  set.seed(1); a <- run.all.gof(f$ex, include_slow = FALSE, install = "no")
  set.seed(1); b <- run.all.gof(f$om, include_slow = FALSE, install = "no")
  expect_identical(a$Test, b$Test)
  expect_equal(a$Statistic, b$Statistic)
  expect_equal(a$p_value, b$p_value)
  expect_false(any(grepl("error", a$Note)))
})

test_that("edge.gof, ef.gof, localize.gof and deepgof1 accept an na.exclude fit", {
  f <- na_fits()
  expect_equal(as.data.frame(edge.gof(f$ex))$p_value, as.data.frame(edge.gof(f$om))$p_value)
  expect_equal(ef.gof(f$ex)$p_value, ef.gof(f$om)$p_value)
  a <- localize.gof(f$ex, B = 19, seed = 1, plot = FALSE)
  b <- localize.gof(f$om, B = 19, seed = 1, plot = FALSE)
  expect_equal(a$adjusted, b$adjusted)
  skip_on_cran()
  set.seed(3); da <- deepgof1(f$ex, B = 19)
  set.seed(3); db <- deepgof1(f$om, B = 19)
  expect_true(is.finite(da$statistic))
  expect_equal(da$statistic, db$statistic)
  expect_equal(da$p.value, db$p.value)
})
