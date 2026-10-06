test_that("plot.deepgof1() draws the map of two covariates and of one", {
  set.seed(31)
  n  <- 300
  x1 <- rnorm(n); x2 <- rnorm(n)
  y  <- rbinom(n, 1, plogis(0.5 * x1 + 0.5 * x2))
  f <- tempfile(fileext = ".pdf")
  grDevices::pdf(f, width = 4.5, height = 4)
  on.exit({ grDevices::dev.off(); unlink(f) }, add = TRUE)

  d2 <- deepgof1(glm(y ~ x1 + x2, family = binomial()), B = 19)
  expect_identical(dim(d2$map), c(d2$K, d2$K))
  expect_invisible(plot(d2))
  expect_silent(plot(d2, colour = FALSE, digits = 2))

  d1 <- deepgof1(glm(y ~ x1, family = binomial()), B = 19)
  expect_length(d1$axes, 1L)
  expect_silent(plot(d1))

  bad <- d2
  bad$map <- NULL
  expect_error(plot(bad), "no residual map")
})
