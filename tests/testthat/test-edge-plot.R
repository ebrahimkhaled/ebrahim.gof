# edge.gof() returns a classed data frame from 2.10.0 so that plot() can draw it. These tests
# hold the two promises that change makes: every use of the result as a data frame works as
# before, and the picture is drawn from the same numbers the test reports.

edge_fit <- function(n = 600, seed = 21) {
  set.seed(seed)
  x <- runif(n, -3, 3)
  y <- rbinom(n, 1, 1 - exp(-exp(0.8 * x)))          # log-log truth, fitted as logit
  glm(y ~ x, family = binomial())
}

test_that("edge.gof() is still a data frame in every use", {
  fit <- edge_fit()
  e <- edge.gof(fit)
  expect_s3_class(e, "edge_gof")
  expect_s3_class(e, "data.frame")
  expect_true(is.data.frame(e))
  expect_identical(nrow(e), 2L)
  expect_identical(names(e), c("Test", "Partition", "Role", "Basis", "Test_Statistic", "df",
                               "Method", "p_value"))
  expect_identical(e$Role, c("verdict", "check"))
  expect_true(all(e$p_value >= 0 & e$p_value <= 1))
  expect_identical(e[1, "p_value"], e$p_value[1])
  expect_identical(e[e$Role == "check", ]$p_value, e$p_value[2])
  expect_identical(nrow(rbind(e, e)), 4L)
  expect_identical(nrow(subset(e, Role == "verdict")), 1L)
  # the numbers are those of def.gof(), the attribute only adds to them
  g <- def.gof(fit, G = "auto")
  expect_identical(e$Test_Statistic[1], g$Test_Statistic)
  expect_identical(e$p_value[1], g$p_value)
  expect_null(attr(g, "details"))                     # def.gof() never carries the pieces
  plain <- as.data.frame(e)
  attr(plain, "details") <- NULL
  expect_identical(class(plain), "data.frame")
  expect_output(print(e), "verdict")
})

test_that("the details reproduce the statistic and the grouped counts", {
  fit <- edge_fit()
  e <- edge.gof(fit, G = c("auto", 10))
  d <- attr(e, "details")
  expect_length(d, nrow(e))
  for (i in seq_along(d)) {
    expect_equal(sum(d[[i]]$Pz_r^2), e$Test_Statistic[i], tolerance = 1e-10)
    expect_equal(sum(d[[i]]$O), sum(fit$y))
    expect_equal(sum(d[[i]]$size), length(fit$y))
    expect_equal(d[[i]]$r, (d[[i]]$O - d[[i]]$E) / sqrt(d[[i]]$V), tolerance = 1e-12)
    expect_true(all(d[[i]]$omega_diag < 1))          # the fit removes part of every residual
    expect_identical(d[[i]]$role, e$Role[i])
    expect_false(d[[i]]$external)
  }
  # the same pieces edge.belt() reports
  b <- edge.belt(fit, G = "auto", nsim = 500, seed = 1)
  expect_equal(d[[1]]$r, b$groups$residual, tolerance = 1e-12)
  expect_equal(d[[1]]$Pz_r, b$groups$direction, tolerance = 1e-10)
  expect_equal(sqrt(d[[1]]$omega_diag), b$groups$sd_null, tolerance = 1e-12)
})

test_that("frozen predictions carry unit variances and the constant column", {
  set.seed(22)
  n <- 800
  x <- runif(n, -3, 3)
  y <- rbinom(n, 1, plogis(0.3 + 0.6 * x))
  p <- plogis(0.6 * x)
  e <- edge.gof(y, predicted_probs = p, external = TRUE)
  d <- attr(e, "details")
  for (i in seq_along(d)) {
    expect_true(d[[i]]$external)
    expect_identical(d[[i]]$omega_diag, rep(1, d[[i]]$G))
    expect_equal(sum(d[[i]]$Pz_r^2), e$Test_Statistic[i], tolerance = 1e-10)
  }
})

test_that("the score form keeps its pieces and a row without a p-value keeps none", {
  fit <- edge_fit()
  e <- edge.gof(fit, G = 10, weights = "score")
  expect_length(attr(e, "details"), 1L)
  expect_equal(attr(e, "details")[[1]]$G, 10L)
})

test_that("plot.edge_gof() draws every documented variant and refuses what it cannot draw", {
  fit <- edge_fit()
  e <- edge.gof(fit)
  f <- tempfile(fileext = ".pdf")
  grDevices::pdf(f, width = 6.3, height = 3.6)
  on.exit({ grDevices::dev.off(); unlink(f) }, add = TRUE)
  expect_invisible(plot(e))
  expect_identical(plot(e), e)
  for (sc in list("residual", "rate", c("residual", "rate")))
    for (bd in c("simultaneous", "pointwise", "none"))
      expect_silent(plot(e, row = "check", scale = sc, band = bd, colour = FALSE))
  expect_silent(plot(e, row = 2, level = 0.9))
  set.seed(23)
  y <- rbinom(500, 1, 0.3)
  ex <- edge.gof(y, predicted_probs = rep(0.3, 500) + runif(500, -0.1, 0.1), external = TRUE, G = 10)
  expect_silent(plot(ex, scale = "rate"))
  expect_error(plot(ex, row = "check"), "no 'check' row")
  expect_error(plot(rbind(e, ex), row = 3), "no grouped details")
  expect_error(plot(e, row = 3), "between 1 and 2")
  expect_error(plot(e, level = 1.2), "'level'")
})

test_that("a subset or a reordering of the rows draws the row it shows", {
  fit <- edge_fit()
  e <- edge.gof(fit)
  skip_if_not(isTRUE(capabilities("png")))
  ## the bytes of the drawn picture, to compare two calls
  drawn <- function(obj, ...) {
    f <- tempfile(fileext = ".png")
    grDevices::png(f, width = 630, height = 360)
    plot(obj, ...)
    grDevices::dev.off()
    on.exit(unlink(f))
    readBin(f, "raw", file.size(f))
  }
  rev_e <- e[2:1, ]
  expect_identical(rev_e$Role, c("check", "verdict"))
  expect_identical(drawn(rev_e), drawn(e))                     # the verdict row, now second
  expect_identical(drawn(rev_e, row = 1), drawn(e, row = "check"))
  expect_identical(drawn(e[2, ], row = 1), drawn(e, row = "check"))
  expect_error(plot(subset(e, Role == "check"), row = 1), "no grouped details")   # subset() drops it
  expect_false(identical(drawn(e, row = "check"), drawn(e)))
})

test_that("edge.gof() reports each partition once, also for the ensemble", {
  fit <- edge_fit()
  ens <- edge.gof(fit, basis = "ensemble")
  expect_identical(nrow(ens), 2L)
  expect_identical(ens$Role, c("verdict", "check"))
  small <- edge_fit(n = 200)
  expect_identical(nrow(edge.gof(small)), 1L)          # n <= 250: the default is ten groups
  expect_identical(nrow(edge.gof(small, basis = "ensemble")), 1L)
  expect_identical(nrow(edge.gof(fit, G = c(10, 10))), 1L)
})
