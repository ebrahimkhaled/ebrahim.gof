# The projection test (Escanciano 2006; Liu et al. 2024): the closed-form weight against
# its definition as an integral over the sphere, estimated by Monte Carlo over directions.

# A_ij by Monte Carlo: the average over M uniform directions w of
# sum_l I(X_i'w <= X_l'w) I(X_j'w <= X_l'w). Returns A and the per-direction T.
mc_weight <- function(X, e, M) {
  n <- nrow(X); p <- ncol(X)
  A <- matrix(0, n, n); Tw <- numeric(M)
  for (m in seq_len(M)) {
    w  <- stats::rnorm(p); w <- w / sqrt(sum(w * w))
    pr <- drop(X %*% w)
    I  <- outer(pr, pr, "<=") + 0          # I[i, l] = I(X_i'w <= X_l'w)
    Aw <- tcrossprod(I)
    A  <- A + Aw
    Tw[m] <- sum(e * (Aw %*% e)) / n^2
  }
  list(A = A / M, Tw = Tw)
}

test_that("the projection weight matches a Monte Carlo over random directions", {
  set.seed(11)
  n <- 12
  X <- cbind(rnorm(n), rnorm(n))
  X[2, ] <- X[1, ]                          # a duplicated row exercises the zero-difference rule
  e <- rnorm(n)
  A <- ebrahim.gof:::.proj_weight(X)
  expect_equal(A, t(A))
  expect_equal(diag(A)[3], n / 2 + 0.5)     # X_3 unique: 1/2 per l, plus 1/2 at l = 3
  set.seed(12)
  mc <- mc_weight(X, e, M = 20000)
  Tclosed <- sum(e * (A %*% e)) / n^2
  se <- stats::sd(mc$Tw) / sqrt(length(mc$Tw))
  expect_lt(abs(Tclosed - mean(mc$Tw)) / se, 4.5)
  expect_lt(max(abs(mc$A - A)), 0.2)        # entrywise, well inside Monte Carlo error
  # the orientation printed in Liu et al. (arccos without "pi -") is far outside it
  Aalt <- n / 2 - A                          # sum_l angle/(2 pi) instead of (pi - angle)/(2 pi), off the zero cases
  expect_gt(max(abs(mc$A - Aalt)), 1)
})

test_that("with one covariate the weight is exact over the two directions", {
  set.seed(3)
  x <- c(round(rnorm(15), 1), 0.3, 0.3)    # ties included
  n <- length(x)
  A <- ebrahim.gof:::.proj_weight(cbind(x))
  exact <- 0.5 * (tcrossprod(outer(x, x, "<=") + 0) + tcrossprod(outer(x, x, ">=") + 0))
  expect_equal(A, exact, tolerance = 1e-12)
  # the general (angle) path gives the same answer on the same points embedded in 2-D
  A2 <- ebrahim.gof:::.proj_weight(cbind(x, 0))
  expect_equal(A2, A, tolerance = 1e-10)
})

test_that("projection.gof returns a valid, reproducible htest", {
  set.seed(5)
  n <- 80
  x1 <- rnorm(n); x2 <- rnorm(n)
  y  <- rbinom(n, 1, plogis(0.4 * x1 - 0.3 * x2))
  fit <- glm(y ~ x1 + x2, family = binomial())
  set.seed(9); r1 <- projection.gof(fit, B = 49)
  set.seed(9); r2 <- projection.gof(fit, B = 49)
  expect_s3_class(r1, "htest")
  expect_identical(r1$p.value, r2$p.value)
  expect_true(r1$p.value > 0 && r1$p.value <= 1)
  expect_equal(unname(r1$parameter), 49L)
  expect_length(r1$boot, 49)
  # the statistic is e'Ae / n^2 on the raw covariates
  A <- ebrahim.gof:::.proj_weight(cbind(x1, x2))
  e <- y - fitted(fit)
  expect_equal(unname(r1$statistic), sum(e * (A %*% e)) / n^2)
  # p-value = (1 + #{T* >= T}) / (B + 1)
  expect_equal(r1$p.value, (1 + sum(r1$boot >= r1$statistic)) / 50)
  expect_error(projection.gof(glm(y ~ 1, family = binomial())), "no covariates")
})

test_that("projection.gof rejects an omitted interaction", {
  skip_on_cran()
  set.seed(21)
  n  <- 200
  x1 <- rnorm(n); x2 <- rnorm(n)
  y  <- rbinom(n, 1, plogis(0.3 * x1 - 0.3 * x2 + 2 * x1 * x2))
  fit <- glm(y ~ x1 + x2, family = binomial())
  expect_lt(projection.gof(fit, B = 99)$p.value, 0.05)
})

test_that("the Projection row of run.all.gof agrees with projection.gof", {
  set.seed(8)
  n <- 60
  x <- rnorm(n)
  y <- rbinom(n, 1, plogis(0.5 * x))
  fit <- glm(y ~ x, family = binomial())
  set.seed(1)
  row <- run.all.gof(fit, tests = "Projection", control = list(Projection = list(B = 19)))
  set.seed(1)
  r <- projection.gof(fit, B = 19)
  expect_equal(row$p_value, r$p.value)
  expect_equal(row$Family, "Bootstrap")
  expect_false("Projection" %in% run.all.gof(fit, include_slow = FALSE)$Test)
  big <- run.all.gof(fit, tests = "Projection", control = list(Projection = list(max_n = 10)))
  expect_true(is.na(big$p_value))
  expect_match(big$Note, "Not run")
})
