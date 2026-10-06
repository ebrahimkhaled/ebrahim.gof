## Reference code for the EF tests.
##
## ls_groups(), paul_G() and ls_stats() are copied from large_sample.R, and ef_parts() from
## ef_exact.R, of the reproduction code of Ebrahim, Khattab and El-Kotory (2026), "A modified
## Hosmer-Lemeshow goodness-of-fit test for asymmetric links: second-order power and robustness",
## arXiv:2607.15454. ls_stats() is trimmed to the EF and HL statistics and their references,
## ef_parts() always takes its groups, and paul_G() is written without the native pipe (the
## package supports R >= 3.5); otherwise the bodies are unchanged. ef_v290() is the grouped
## branch of ef.gof() in ebrahim.gof 2.9.0, kept to check that the old groups still reproduce
## an old result.

ref_ls_groups <- function(p, G) {
  n <- length(p); o <- order(p, runif(n))
  g <- integer(n); g[o] <- ceiling(seq_len(n) * G / n)
  g
}

ref_paul_G <- function(n, y) {
  n1 <- sum(y)
  floor(max(10, min(n1 / 2, (n - n1) / 2, 2 + 8 * (n / 1000)^2)))
}

ref_ls_stats <- function(y, p, X, G, g = NULL) {
  if (is.null(g)) g <- ref_ls_groups(p, G)
  G   <- max(g)
  ng  <- tabulate(g, G)
  v   <- p * (1 - p)
  o   <- rowsum(y, g, reorder = TRUE)[, 1]
  e   <- rowsum(p, g, reorder = TRUE)[, 1]
  pb  <- e / ng
  V   <- ng * pb * (1 - pb)
  cc  <- 1 - 2 * pb
  HL  <- sum((o - e)^2 / V)
  C   <- sum(cc * (o - e) / V)
  EF  <- HL - C
  tau <- rowsum(v, g, reorder = TRUE)[, 1]
  k3  <- rowsum(v * (1 - 2 * p), g, reorder = TRUE)[, 1]
  k4  <- rowsum(v * (1 - 6 * v), g, reorder = TRUE)[, 1]
  U   <- rowsum(v * X, g, reorder = TRUE)
  K   <- rowsum(v * (1 - 2 * p) * X, g, reorder = TRUE)
  Mi  <- solve(crossprod(X, v * X))
  dtr <- sum(diag(crossprod(U, U / V) %*% Mi))
  mu  <- sum(tau / V) - dtr
  k1  <- colSums((tau / V) * cc * U / V)
  k2  <- colSums(K / V)
  vHL <- sum((2 * tau^2 + k4) / V^2) - 2 * drop(k1 %*% Mi %*% k2) + drop(k1 %*% Mi %*% k1)
  dd  <- k1 - colSums(cc * U / V)
  k3v <- colSums((K - cc * U) / V)
  vEF <- sum((2 * tau^2 + k4 - 2 * cc * k3 + cc^2 * tau) / V^2) - 2 * drop(dd %*% Mi %*% k3v) + drop(dd %*% Mi %*% dd)
  muEF <- sum(tau / V) - dtr * (1 - 1 / mean(ng))
  zEF <- (EF - muEF) / sqrt(vEF); zHL <- (HL - mu) / sqrt(vHL)
  pEFn <- pnorm(zEF, lower.tail = FALSE)
  c(G = G, HL = HL, C = C, EF = EF,
    p_HL_chisq = pchisq(HL, G - 2, lower.tail = FALSE), p_EF_chisq = pchisq(EF, G - 2, lower.tail = FALSE),
    zHL = zHL, zEF = zEF, p_HL_norm = pnorm(zHL, lower.tail = FALSE), p_EF_norm = pEFn)
}

ref_ef_parts <- function(y, p, X, G = 10, groups = NULL) {
  g   <- if (is.null(groups)) stop("pass the groups") else as.integer(factor(groups))
  G   <- max(g)
  v   <- p * (1 - p)
  ng  <- tabulate(g, G)
  o   <- as.numeric(rowsum(y, g, reorder = TRUE))
  e   <- as.numeric(rowsum(p, g, reorder = TRUE))
  pb  <- e / ng
  V   <- ng * pb * (1 - pb)
  Lam <- as.numeric(rowsum(v, g, reorder = TRUE))
  U   <- rowsum(v * X, g, reorder = TRUE)
  M   <- crossprod(X, v * X)
  S   <- (diag(Lam, G) - U %*% solve(M, t(U))) / sqrt(outer(V, V))
  r   <- (o - e) / sqrt(V)
  list(r = r, b = (1 - 2 * pb) / sqrt(V), Sigma = (S + t(S)) / 2, pbar = pb, V = V, o = o, e = e)
}

ef_v290 <- function(y, predicted_probs, G, method = "chisq") {
  n <- length(y)
  sorted_indices <- order(predicted_probs)
  sorted_y <- y[sorted_indices]
  sorted_probs <- predicted_probs[sorted_indices]
  group_size <- floor(n / G)
  remainder <- n %% G
  group_assignments <- rep(1:G, each = group_size)
  if (remainder > 0) group_assignments <- c(group_assignments, rep(G, remainder))
  y_group <- as.numeric(tapply(sorted_y, group_assignments, sum))
  m_group <- as.numeric(tapply(sorted_y, group_assignments, length))
  prob_group <- as.numeric(tapply(sorted_probs, group_assignments, mean))
  prob_bounded <- pmax(pmin(prob_group, 0.999999), 0.000001)
  variance_terms <- m_group * prob_bounded * (1 - prob_bounded)
  pearson_components <- ((y_group - m_group * prob_bounded)^2) / variance_terms
  farrington_correction <- ((1 - 2 * prob_bounded) / variance_terms) *
    (y_group - m_group * prob_bounded)
  t_ef <- sum(pearson_components) - sum(farrington_correction)
  z <- (t_ef - (G - 2)) / sqrt(2 * (G - 2))
  p <- if (identical(method, "normal")) 1 - pnorm(z) else pchisq(t_ef, df = G - 2, lower.tail = FALSE)
  c(T_EF = t_ef, z = z, p_value = p)
}

## data of the paper's checks: x1 uniform, x2 binary, logistic truth
ef_ref_data <- function(n, seed, tied = FALSE) {
  set.seed(seed)
  x1 <- if (tied) sample(seq(-3, 3, by = 0.5), n, replace = TRUE) else runif(n, -3, 3)
  x2 <- rbinom(n, 1, 0.5)
  y  <- rbinom(n, 1, plogis(0.6 * x1 + 0.5 * x2))
  fit <- glm(y ~ x1 + x2, family = binomial())
  list(fit = fit, y = y, p = as.numeric(fitted(fit)), X = model.matrix(fit))
}
