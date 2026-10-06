## McCullagh's and Copas' rows are computed from the p x p inverse; up to 2.10.0 (pre-release)
## they formed n x n matrices. The dense bodies are kept here verbatim, and the shipped rows must
## agree with them, sparse and grouped, to rounding.

dense_mccullagh <- function(X, ph, y, m = rep(1, length(y))) {
  n <- sum(m); p <- ncol(X); J <- length(y)
  V  <- m * ph * (1 - ph)
  X2 <- sum((y - m * ph)^2 / V)
  W  <- diag(V)
  inv <- solve(t(X) %*% W %*% X)
  H <- X %*% inv %*% t(X); h <- diag(H)
  u <- (1 - 2 * ph) / V; uhat <- as.numeric(H %*% W %*% u)
  k3 <- V * (1 - 2 * ph); k4 <- m * ph * (1 - ph) * (1 - 6 * ph * (1 - ph))
  E <- (J - p) - 0.5 * sum(k4 * h / V) + 0.5 * sum(uhat * k3 * h)
  RSSu <- as.numeric(t(u) %*% (W - W %*% H %*% W) %*% u)
  (X2 - E) / sqrt(2 * sum((m - 1) / m) + (1 - p / n) * RSSu)
}

dense_copas <- function(X, ph, y) {
  V <- ph * (1 - ph); W <- diag(V); c12 <- 1 - 2 * ph
  M <- diag(length(ph)) - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X)
  (sum((y - ph)^2) - sum(V)) / sqrt(as.numeric(t(c12) %*% M %*% W %*% c12))
}

test_that("McCullagh and Copas agree with their n x n forms", {
  for (s in 1:3) {
    set.seed(s)
    n <- 300
    d <- data.frame(x1 = rnorm(n), x2 = runif(n, -2, 2), g = factor(sample(letters[1:3], n, TRUE)))
    d$y <- rbinom(n, 1, plogis(-0.4 + 0.9 * d$x1 - 0.5 * d$x2 + 0.4 * d$x1^2 * (s == 2)))
    fit <- glm(y ~ x1 + x2 + g, family = binomial(), data = d)
    ph <- pmin(pmax(fit$fitted.values, 1e-6), 1 - 1e-6); X <- model.matrix(fit)
    res <- run.all.gof(fit, tests = c("McCullagh", "Copas-RSS"), install = "no")
    expect_equal(res$Statistic[res$Test == "McCullagh"], dense_mccullagh(X, ph, fit$y), tolerance = 1e-10)
    z <- dense_copas(X, ph, fit$y)
    expect_equal(res$Statistic[res$Test == "Copas-RSS"], z, tolerance = 1e-10)
    expect_equal(res$p_value[res$Test == "Copas-RSS"], pchisq(z^2, 1, lower.tail = FALSE), tolerance = 1e-10)
    expect_true(is.na(res$df[res$Test == "Copas-RSS"]))
  }
  ## grouped: a design with repeated covariate patterns
  set.seed(9)
  d <- data.frame(a = sample(1:6, 800, TRUE), b = rbinom(800, 1, 0.4))
  d$y <- rbinom(800, 1, plogis(-1 + 0.3 * d$a + 0.5 * d$b))
  fit <- glm(y ~ a + b, family = binomial(), data = d)
  ph <- pmin(pmax(fit$fitted.values, 1e-6), 1 - 1e-6); X <- model.matrix(fit)
  key <- paste(d$a, d$b); u <- !duplicated(key)
  m  <- as.numeric(table(key)[key[u]]); yg <- as.numeric(tapply(d$y, key, sum)[key[u]])
  res <- run.all.gof(fit, tests = "McCullagh (grouped)", install = "no")
  expect_equal(res$Statistic, dense_mccullagh(X[u, , drop = FALSE], ph[u], yg, m), tolerance = 1e-10)
})
