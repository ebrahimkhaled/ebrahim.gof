## cut points that reproduce def.gof()'s rank groups on a given set of predictions
rank_breaks <- function(p, G) {
  ps <- sort(p); n <- length(p)
  idx <- vapply(seq_len(G - 1), function(g) max(which(pmin(ceiling(seq_len(n) / (n / G)), G) == g)), 0L)
  (ps[idx] + ps[idx + 1]) / 2
}

test_that("the streamed statistic equals the one-shot external statistic, however it is batched", {
  set.seed(21)
  n <- 3000
  p <- plogis(rnorm(n, -1, 1))
  y <- rbinom(n, 1, plogis(qlogis(p) * 0.8))
  ref <- def.gof(y, predicted_probs = p, G = 10, external = TRUE)
  s1 <- update(edge.stream(breaks = rank_breaks(p, 10)), y, p)
  s2 <- edge.stream(breaks = rank_breaks(p, 10))
  for (b in split(seq_len(n), rep(1:30, length.out = n))) s2 <- update(s2, y[b], p[b])
  o1 <- summary(s1); o2 <- summary(s2)
  expect_equal(o1$Test_Statistic, ref$Test_Statistic, tolerance = 1e-10)
  expect_equal(o1$df, ref$df)
  expect_equal(o1$p_value, ref$p_value, tolerance = 1e-10)
  expect_equal(o2$Test_Statistic, o1$Test_Statistic, tolerance = 1e-10)
  expect_equal(o2$n, n)
})

test_that("the order of the records does not matter", {
  set.seed(22)
  p <- plogis(rnorm(2000, -1, 1)); y <- rbinom(2000, 1, p)
  s0 <- edge.stream(p_ref = p, G = 10)
  a <- summary(update(s0, y, p))
  o <- sample.int(2000)
  b <- summary(update(s0, y[o], p[o]))
  expect_equal(a$Test_Statistic, b$Test_Statistic, tolerance = 1e-10)
})

test_that("every basis runs and gives d + 1 degrees of freedom", {
  set.seed(23)
  p <- plogis(rnorm(4000, -1, 1)); y <- rbinom(4000, 1, p)
  for (b in c("poly3", "poly2", "stukel", "sym")) {
    o <- summary(update(edge.stream(p_ref = p, G = 10, basis = b), y, p))
    expect_equal(o$df, c(poly3 = 4, poly2 = 3, stukel = 4, sym = 2)[[b]])
    expect_true(o$p_value > 0 && o$p_value < 1)
  }
})

test_that("a miscalibrated stream is flagged and an empty one is not tested", {
  set.seed(24)
  s <- edge.stream(p_ref = plogis(rnorm(1000, -1, 1)), G = 10)
  expect_true(is.na(summary(s)$p_value))
  expect_output(print(s), "not enough data")
  for (i in 1:50) {
    p <- plogis(rnorm(100, -1, 1))
    s <- update(s, rbinom(100, 1, plogis(qlogis(p) * 0.6)), p)   # overfitted model
  }
  expect_lt(summary(s)$p_value, 1e-3)
  expect_equal(summary(s)$n, 5000)
})

test_that("bad input is refused", {
  expect_error(edge.stream(), "p_ref")
  expect_error(edge.stream(breaks = c(0.2, 1.2)), "between 0 and 1")
  s <- edge.stream(breaks = c(0.1, 0.2))
  expect_error(update(s, c(0, 1), c(0.5)), "lengths")
  expect_error(update(s, c(0, 2), c(0.5, 0.5)), "binary")
})
