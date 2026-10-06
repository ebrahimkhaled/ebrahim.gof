## ef.gof() against the reproduction code of the EF paper (helper-ef-reference.R): the same data and
## the same seed must give the same groups, statistics and p-values.

test_that("ten groups: EF, HL, C and the chi-square p-value match the paper's code", {
  for (n in c(1000, 1003, 487)) {
    d <- ef_ref_data(n, seed = n)
    set.seed(1); r <- ef.gof(d$fit)
    set.seed(1); s <- ref_ls_stats(d$y, d$p, d$X, 10)
    expect_equal(r$Test_Statistic, s[["EF"]], tolerance = 1e-10)
    expect_equal(r$HL, s[["HL"]], tolerance = 1e-10)
    expect_equal(r$C, s[["C"]], tolerance = 1e-10)
    expect_equal(r$p_value, s[["p_EF_chisq"]], tolerance = 1e-10)
    expect_identical(r$df, 8)
    expect_identical(r$Reference, "chisq")
  }
})

test_that("the read-out r_g, b_g and c_g matches ef_parts() on the same groups", {
  d <- ef_ref_data(1003, seed = 3)
  set.seed(2); r <- ef.gof(d$y, d$p, G = 10)
  set.seed(2); g <- ref_ls_groups(d$p, 10)
  z <- ref_ef_parts(d$y, d$p, d$X, groups = g)
  gr <- attr(r, "groups")
  expect_equal(gr$r, z$r, tolerance = 1e-10)
  expect_equal(gr$b, z$b, tolerance = 1e-10)
  expect_equal(gr$c, z$b * z$r, tolerance = 1e-10)
  expect_equal(sum(gr$c), r$C, tolerance = 1e-12)
  expect_equal(gr$n, tabulate(g, 10))
  expect_true(max(gr$n) - min(gr$n) <= 1)             # sizes differ by at most one
})

test_that("many groups: the normal reference matches the paper's code", {
  for (n in c(3000, 5000)) {
    d <- ef_ref_data(n, seed = n + 1)
    G <- max(10, ceiling(n / 25))
    set.seed(3); r <- ef.gof(d$fit, G = "auto", reference = "normal")
    set.seed(3); s <- ref_ls_stats(d$y, d$p, d$X, G)
    expect_equal(r$G, G)
    expect_equal(r$Test_Statistic, s[["zEF"]], tolerance = 1e-10)
    expect_equal(r$p_value, s[["p_EF_norm"]], tolerance = 1e-10)
    expect_true(is.na(r$df))
    ## the same with y, p and X
    set.seed(3); r2 <- ef.gof(d$y, d$p, G = G, reference = "normal", X = d$X)
    expect_equal(r2$p_value, r$p_value, tolerance = 1e-12)
  }
  ## Paul's rule
  d <- ef_ref_data(5000, seed = 9)
  set.seed(4); r <- ef.gof(d$fit, G = "paul", reference = "normal")
  set.seed(4); s <- ref_ls_stats(d$y, d$p, d$X, ref_paul_G(5000, d$y))
  expect_equal(r$G, s[["G"]])
  expect_equal(r$p_value, s[["p_EF_norm"]], tolerance = 1e-10)
})

test_that("tied fitted risks are broken at random, through set.seed()", {
  d <- ef_ref_data(1000, seed = 5, tied = TRUE)
  expect_gt(anyDuplicated(d$p), 0)
  for (k in 1:3) {
    set.seed(k); r <- ef.gof(d$fit)
    set.seed(k); s <- ref_ls_stats(d$y, d$p, d$X, 10)
    expect_equal(r$Test_Statistic, s[["EF"]], tolerance = 1e-10)
    expect_equal(r$p_value, s[["p_EF_chisq"]], tolerance = 1e-10)
    set.seed(k); rn <- ef.gof(d$fit, G = 40, reference = "normal")
    set.seed(k); sn <- ref_ls_stats(d$y, d$p, d$X, 40)
    expect_equal(rn$p_value, sn[["p_EF_norm"]], tolerance = 1e-10)
  }
  set.seed(7); a <- ef.gof(d$fit)
  set.seed(7); b <- ef.gof(d$fit)
  expect_identical(a, b)
  ## untied data draw nothing, so the stream is left as it was
  u <- ef_ref_data(500, seed = 6)
  set.seed(8); ef.gof(u$fit); x1 <- runif(1)
  set.seed(8); x2 <- runif(1)
  expect_identical(x1, x2)
})

test_that("one group per pattern matches the paper's code with the pattern labels", {
  d <- ef_ref_data(600, seed = 10, tied = TRUE)
  r <- ef.gof(d$fit, groups = "pattern")
  g <- as.integer(factor(d$p))
  s <- ref_ls_stats(d$y, d$p, d$X, g = g)
  expect_equal(r$G, max(g))
  expect_equal(r$Test_Statistic, s[["EF"]], tolerance = 1e-10)
  expect_equal(r$p_value, s[["p_EF_chisq"]], tolerance = 1e-10)
  z <- ref_ef_parts(d$y, d$p, d$X, groups = d$p)
  expect_equal(attr(r, "groups")$c, z$b * z$r, tolerance = 1e-10)
  rn <- ef.gof(d$fit, groups = "pattern", reference = "normal")
  expect_equal(rn$p_value, s[["p_EF_norm"]], tolerance = 1e-10)
  ## any labels
  expect_equal(ef.gof(d$fit, groups = g)$p_value, r$p_value, tolerance = 1e-12)
  expect_error(ef.gof(d$y, d$p, groups = rep(1:2, 300)), "fewer than 3")
  expect_error(ef.gof(d$y, d$p, groups = 1:5), "one label per observation")
})

test_that("the groups of 2.9.0 reproduce its result", {
  d <- ef_ref_data(1003, seed = 11)
  n <- length(d$y)
  old <- pmin(ceiling(rank(d$p, ties.method = "first") / floor(n / 10)), 10)
  r <- ef.gof(d$y, d$p, groups = old)
  v <- ef_v290(d$y, d$p, 10)
  expect_equal(r$Test_Statistic, v[["T_EF"]], tolerance = 1e-10)
  expect_equal(r$p_value, v[["p_value"]], tolerance = 1e-10)
  expect_warning(rt <- ef.gof(d$y, d$p, groups = old, method = "normal"), "versions up to 2.9.0")
  expect_equal(rt$Test_Statistic, ef_v290(d$y, d$p, 10, "normal")[["z"]], tolerance = 1e-10)
  expect_identical(rt$Reference, "thesis")
})

test_that("the size rule and the design-matrix requirement are enforced", {
  set.seed(12)
  x <- rnorm(150)
  y <- rbinom(150, 1, plogis(-3.5 + 0.5 * x))         # rare events: nearly empty groups
  fit <- glm(y ~ x, family = binomial())
  expect_warning(ef.gof(fit), class = "ef_sparse_groups")
  expect_error(ef.gof(y, fitted(fit), reference = "normal"), "design matrix")
  ## the battery puts the warning in the Note of the EF row
  res <- as.data.frame(run.all.gof(fit, tests = c("EF", "EF-normal"), include_slow = FALSE))
  expect_match(res$Note[res$Test == "EF"], "expected event")
  expect_true(is.finite(res$p_value[res$Test == "EF-normal"]))
})
