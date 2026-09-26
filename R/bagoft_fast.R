# A fast re-expression of the BAGofT test (Zhang, Ding and Yang 2023).
#
# ATTRIBUTION. The procedure below is the one implemented by the CRAN package BAGofT 1.0.0
# (Jiawei Zhang, Jie Ding and Yuhong Yang; GPL-3): its functions BAGofT(), BAGofT_multi(),
# BAGofT_sin(), testGlmBi(), parRF() and dcPre(). The control flow, the defaults, the
# adaptive choice of the number of groups and the merge of undersized groups follow that
# code step by step, and parts of it (the Kmax loop, the grouping of the test set and the
# merge loop) are adapted from it nearly verbatim. That code is copyright its authors and is
# included here under the GPL-3, the licence of both packages.
#
# WHAT CHANGES. Only the entry points. Every random number is drawn by the same call, in the
# same order, as in the package: sample() for each split, randomForest() for each forest
# (importance = TRUE is kept because the permutation importance draws random numbers), and
# rbinom() for each simulated response. The deterministic steps are the same arithmetic on
# the same numbers, reached through cheaper routes:
#   glm(formula, data)            -> glm.fit(X, y) on the model matrix, built once
#   predict(glm, se.fit = TRUE)   -> linkinv(X %*% beta), pivoting as predict.lm does
#   resid(glm, "pearson")         -> the formula of residuals.glm
#   xtabs(v ~ f)                  -> tapply(v, f, sum, default = 0L), which xtabs calls
#   cut(x, breaks, TRUE)          -> .bincode(x, breaks, TRUE, TRUE), cut's own codes
#   randomForest(res ~ ., datRf)  -> randomForest(x, y), the formula method's target
#   mtry = floor(1/3) = 0         -> mtry = 1, the value randomForest resets 0 to
# The p-values are therefore identical to the package's, not merely close (the package
# tests check identical() against BAGofT on small datasets).

# One split: BAGofT_sin + testGlmBi + parRF(+ dcPre) for the response y.
.bagoft_split <- function(y, s) {
  nr <- s$nr; nt <- s$nt; ne <- s$ne
  trainIn <- sample(c(1:nr), nt)
  testIdx <- c(1:nr)[-trainIn]                 # datset[-trainIn, ] keeps the original order
  yT <- y[trainIn]
  yE <- y[testIdx]
  XT <- s$X[trainIn, , drop = FALSE]
  XE <- s$X[testIdx, , drop = FALSE]

  fit   <- stats::glm.fit(XT, yT, family = s$fam)
  predT <- .bagoft_pred(fit, XT, s$fam)
  predE <- .bagoft_pred(fit, XE, s$fam)
  mu    <- fit$fitted.values
  res   <- (fit$y - mu) * sqrt(fit$prior.weights) / sqrt(s$fam$variance(mu))

  ZT <- s$Z[trainIn, , drop = FALSE]
  if (s$presel) {                              # dcPre(npreSel = 5, type = "V")
    dc   <- t(dcov::mdcor(res, ZT))
    vsel <- s$xcol[order(-as.numeric(dc))[c(1:s$npreSel)]]
  } else {
    vsel <- s$xcol
  }
  mtry <- if (is.null(s$mtry)) floor((if (s$presel) length(vsel) else 1L) / 3) else s$mtry
  mtry <- max(1, min(length(vsel), round(mtry)))   # randomForest's own reset of mtry
  ZTs  <- ZT[, vsel, drop = FALSE]
  rf   <- randomForest::randomForest(x = ZTs, y = res, ntree = s$ntree, maxnodes = s$maxnodes,
                                     mtry = mtry, importance = TRUE)
  trainsetPred <- stats::predict(rf, newdata = ZTs)
  rt  <- round(trainsetPred, digits = 10)
  urt <- unique(rt)

  if (length(urt) == 1) {
    gup <- as.factor(rep(1, ne))
  } else {
    Kmax_adj    <- min(s$Kmax, length(urt))
    chitrainVec <- numeric(Kmax_adj)
    dT <- predT - yT
    vT <- predT * (1 - predT)
    for (gt in c(1:Kmax_adj)) {
      gupt    <- .bagoft_cut(rt, stats::quantile(urt, probs = seq(0, 1, 1/gt), names = FALSE))
      dift    <- abs(tapply(dT, gupt, sum, default = 0L))
      dent    <- tapply(vT, gupt, sum, default = 0L)
      contrit <- (dift)^2/dent
      chitrainVec[gt] <- sum(contrit)
      if (is.nan(sum(contrit))) chitrainVec[gt] <- chitrainVec[(gt - 1)]
    }
    chitrainVec4 <- chitrainVec[-1] - chitrainVec[-length(chitrainVec)]
    Ksel <- which.max(chitrainVec4) + 1
    testsetPred <- stats::predict(rf, newdata = s$Z[testIdx, vsel, drop = FALSE])
    gup <- .bagoft_cut(testsetPred, c(Inf, -Inf, stats::quantile(urt, probs = seq(0, 1, 1/Ksel),
                                                                  names = FALSE)))
    gup <- droplevels(gup)
    freqTab <- table(gup)
    while (min(freqTab) < s$nmin) {            # parRF's merge of undersized groups
      levels(gup)[which(levels(gup) == names(sort(freqTab))[1])] <- names(sort(freqTab))[2]
      freqTab <- table(gup)
    }
  }

  ngp    <- length(levels(gup))
  dif    <- abs(tapply(predE - yE, gup, sum, default = 0L))
  den    <- tapply(predE * (1 - predE), gup, sum, default = 0L)
  contri <- (dif)^2/den
  chisq  <- sum(contri)
  1 - stats::pchisq(chisq, ngp)
}

# predict.lm's product for a glm.fit result, pivoting as it does when rank-deficient
.bagoft_pred <- function(fit, Xn, fam) {
  p   <- fit$rank
  piv <- if (p) fit$qr$pivot[seq_len(p)]
  fam$linkinv(drop(Xn[, piv, drop = FALSE] %*% fit$coefficients[piv]))
}

# cut(x, breaks, include.lowest = TRUE) as a factor of bin codes (cut's labels are always
# distinct, so the grouping is carried by the codes alone)
.bagoft_cut <- function(x, breaks) {
  nb <- length(breaks <- sort.int(as.double(breaks)))
  if (anyDuplicated(breaks)) stop("'breaks' are not unique")
  factor(.bincode(x, breaks, TRUE, TRUE), levels = seq_len(nb - 1L))
}

# BAGofT_multi: nsplits splits, summarised by the mean, median and minimum split p-value
.bagoft_multi <- function(y, s) {
  pv <- numeric(s$nsplits)
  for (j in c(1:s$nsplits)) pv[j] <- .bagoft_split(y, s)
  c(mean(pv), stats::median(pv), min(pv))
}

# The whole test. X: the model matrix; y: the 0/1 response; Z: the data frame the forest
# partitions (every column of the data except the response, in data order); link: the link.
.bagoft_core <- function(X, y, Z, link = "logit", nsplits = 100L, nsim = 100L,
                         ne = NULL, ntree = 60, Kmax = NULL, nmin = NULL, mtry = NULL,
                         maxnodes = NULL) {
  nr <- length(y)
  if (is.null(ne)) ne <- floor(5 * nr^(1/2))
  nt <- nr - ne
  if (is.null(nmin)) nmin <- ceiling(sqrt(ne))
  if (is.null(Kmax)) Kmax <- floor(ne / nmin)
  if (is.null(maxnodes)) maxnodes <- min(nt, ceiling(5 * ncol(Z)))
  xcol <- names(Z)
  s <- list(nr = nr, nt = nt, ne = ne, X = X, Z = Z, xcol = xcol,
            fam = stats::binomial(link = link), presel = length(xcol) > 5L, npreSel = 5L,
            mtry = mtry, ntree = ntree, maxnodes = maxnodes, nmin = nmin, Kmax = Kmax,
            nsplits = as.integer(nsplits))
  obs <- .bagoft_multi(y, s)
  out <- list(pmean = obs[1], pmedian = obs[2], pmin = obs[3])
  if (nsim >= 1) {
    fit0  <- stats::glm.fit(X, y, family = s$fam)
    pdat2 <- .bagoft_pred(fit0, X, s$fam)
    sims  <- matrix(NA_real_, nsim, 3L)
    for (i in c(1:nsim)) {
      ydat2 <- sapply(pdat2, function(x) stats::rbinom(1, 1, x))   # the package's own draw
      sims[i, ] <- .bagoft_multi(ydat2, s)
    }
    out <- list(p.value  = mean(obs[1] > sims[, 1]),
                p.value2 = mean(obs[2] > sims[, 2]),
                p.value3 = mean(obs[3] > sims[, 3]),
                pmean = obs[1], pmedian = obs[2], pmin = obs[3],
                simRes = list(pmeanSim = sims[, 1], pmediansim = sims[, 2],
                              pminsim = sims[, 3]))
  }
  out$settings <- list(nsplits = s$nsplits, nsim = nsim, ne = ne, ntree = ntree, Kmax = Kmax,
                       nmin = nmin, maxnodes = maxnodes, preselected = s$presel)
  out
}

#' Fast BAGofT: the Binary Adaptive Goodness-of-Fit Test
#'
#' @description
#' \code{bagoft.fast()} computes the BAGofT test of Zhang, Ding and Yang (2023) for a
#' fitted binary \code{glm}, and returns the same p-values as the \pkg{BAGofT} package
#' (version 1.0.0) for the same random seed -- identical, not approximately equal. It is
#' the call
#' \preformatted{BAGofT::BAGofT(testGlmBi(formula, link), parRF(), data,
#'                nsplits = 100, ne = floor(5 * sqrt(n)), nsim = 100)}
#' with the package's own defaults, computed without its per-split overhead
#' (formula parsing, model frames, \code{xtabs}, \code{cut} labels and the
#' \code{predict} wrappers, repeated \eqn{(nsim + 1) \times nsplits} times).
#'
#' @details
#' \strong{The test.} Each split fits the model on \eqn{n - n_e} observations, grows a
#' random forest of the Pearson residuals on the covariates, uses the forest to choose an
#' adaptive partition of the covariate space, and computes a Hosmer-Lemeshow-type
#' chi-squared statistic on the \eqn{n_e} held-out observations. The split p-values are
#' averaged over \code{nsplits} splits, and the average is calibrated against \code{nsim}
#' responses simulated from the full-data fit. The calibrated p-value is \code{p.value};
#' \code{p.value2} and \code{p.value3} calibrate the median and the minimum instead.
#'
#' \strong{Read \code{p.value}, not \code{pmean} or \code{pmin}.} \code{pmean},
#' \code{pmedian} and \code{pmin} are the observed \emph{statistics} (summaries of the
#' split p-values), not p-values, and are not uniform under the null.
#'
#' \strong{What the forest sees.} As in the package's default \code{parRF(parVar = ".")},
#' the forest partitions on every column of \code{data} except the response -- not only
#' the model's terms. By default \code{data} is the model frame, so the forest sees the
#' model's variables; pass a wider \code{data} to let it look for structure in covariates
#' the model leaves out, exactly as the package does. With more than five such columns the
#' package's pre-selection runs first in each split (\code{dcPre()}: the five columns with
#' the largest distance correlation with the residuals, computed by \pkg{dcov}), and so it
#' does here.
#'
#' \strong{Differences from the package.} None in the result where the package runs. Where
#' it does not, this function still does: with a single covariate BAGofT 1.0.0 stops (its
#' \code{parRF()} drops the data frame to a vector), while here the forest simply uses the
#' one column; and a formula with transformed terms (e.g. \code{log(x)}) is evaluated
#' through the fitted model's own design matrix. Offsets and non-unit prior weights are
#' not supported, as in the package.
#'
#' \strong{Speed.} The work that remains is the forests themselves:
#' \eqn{(nsim + 1) \times nsplits = 10{,}100} forests with the defaults. At \eqn{n = 200}
#' with two covariates the default test took 99 s against 239 s for the package, and 188 s
#' at \eqn{n = 500} (package about 360 s); the forests themselves are about half of it. Reduce
#' \code{nsim} or \code{nsplits} for a quicker, noisier answer.
#'
#' @param object A fitted binary \code{\link[stats]{glm}} (\code{family = binomial}, any
#'   link) with unit prior weights and no offset.
#' @param data Optional data frame to run the test on, as in \code{BAGofT(data = ...)}.
#'   It must contain the response (numeric 0/1) and the model's variables; every other
#'   column is offered to the forest. Default: the model frame of \code{object}.
#' @param nsplits Number of random splits (default 100, the package default).
#' @param nsim Number of simulated responses used to calibrate the p-value (default 100,
#'   the package default). \code{nsim = 0} returns the statistics only.
#' @param ne Size of the held-out part of each split. Default
#'   \code{floor(5 * sqrt(n))}, the package default.
#' @param ntree,Kmax,nmin,mtry,maxnodes The random-forest partitioner's settings, as in
#'   \code{BAGofT::parRF()}: number of trees (default 60), maximum number of groups
#'   (default \code{floor(ne / nmin)}), minimum group size (default
#'   \code{ceiling(sqrt(ne))}), variables tried at each node (default 1, the value the
#'   package ends up with), and maximum number of terminal nodes (default
#'   \code{min(n - ne, 5 * ncol)}, \code{ncol} the number of covariate columns).
#'
#' @return An object of class \code{"bagoft_fast"}: a list with \code{p.value},
#'   \code{p.value2}, \code{p.value3} (the calibrated p-values for the mean, median and
#'   minimum split p-value), \code{pmean}, \code{pmedian}, \code{pmin} (the observed
#'   statistics), \code{simRes} (the simulated statistics), and \code{settings}. The
#'   elements shared with \code{BAGofT::BAGofT()} have the same names and values.
#'   \code{singleSplit.results} is not returned.
#'
#' @references
#' Zhang J, Ding J, Yang Y (2023). "Is a classification procedure good enough? A
#' goodness-of-fit assessment tool for classification learning." \emph{Journal of the
#' American Statistical Association}, 118(542), 1115--1125.
#' \doi{10.1080/01621459.2021.1979010}
#'
#' Zhang J, Ding J, Yang Y (2021). BAGofT: A Binary Regression Adaptive Goodness-of-Fit
#' Test. R package version 1.0.0. \url{https://CRAN.R-project.org/package=BAGofT}
#'
#' @author The procedure and the code it is adapted from are by Jiawei Zhang, Jie Ding and
#'   Yuhong Yang (package \pkg{BAGofT}, GPL-3). The fast re-expression is by Ebrahim Khaled
#'   Ebrahim \email{ebrahimkhaled@@alexu.edu.eg}.
#'
#' @examples
#' \donttest{
#' if (requireNamespace("randomForest", quietly = TRUE)) {
#'   set.seed(1)
#'   n  <- 100
#'   x1 <- rnorm(n); x2 <- rnorm(n)
#'   y  <- rbinom(n, 1, plogis(0.5 * x1 + 0.5 * x2))
#'   fit <- glm(y ~ x1 + x2, family = binomial())
#'   set.seed(2)
#'   r <- bagoft.fast(fit, nsplits = 20, nsim = 20)    # reduced for speed
#'   r
#'   ## the same numbers from the BAGofT package, same seed:
#'   ## set.seed(2)
#'   ## BAGofT::BAGofT(BAGofT::testGlmBi(y ~ x1 + x2, link = "logit"),
#'   ##                BAGofT::parRF(), data = data.frame(y, x1, x2),
#'   ##                nsplits = 20, nsim = 20)$p.value
#' }
#' }
#'
#' @seealso \code{\link{run.all.gof}}, where this is the engine of the row
#'   \code{"BAGofT"}.
#' @concept goodness-of-fit
#' @concept logistic regression
#' @concept random forest
#' @export
bagoft.fast <- function(object, data = NULL, nsplits = 100, nsim = 100, ne = NULL,
                        ntree = 60, Kmax = NULL, nmin = NULL, mtry = NULL, maxnodes = NULL) {
  if (!inherits(object, "glm") || object$family$family != "binomial")
    stop("bagoft.fast: 'object' must be a fitted binomial glm.")
  if (!requireNamespace("randomForest", quietly = TRUE))
    stop("bagoft.fast: needs the 'randomForest' package.")
  if (!is.null(object$offset) && any(object$offset != 0))
    stop("bagoft.fast: models with an offset are not supported.")
  if (any(object$prior.weights != 1))
    stop("bagoft.fast: needs unit prior weights (binary data).")
  if (is.null(data)) {
    data <- stats::model.frame(object)
    attr(data, "terms") <- NULL
    y <- as.numeric(object$y)
    X <- stats::model.matrix(object)
    Z <- data[-1L]
  } else {
    data <- as.data.frame(data)
    Rsp  <- as.character(stats::formula(object))[2]
    if (!Rsp %in% names(data))
      stop("bagoft.fast: the response '", Rsp, "' is not a column of 'data'.")
    y <- data[[Rsp]]
    if (!is.numeric(y) || !all(y %in% c(0, 1)))
      stop("bagoft.fast: the response must be numeric 0/1.")
    mt <- stats::terms(stats::formula(object), data = data)
    X  <- stats::model.matrix(mt, stats::model.frame(mt, data))
    Z  <- data[, setdiff(names(data), Rsp), drop = FALSE]
  }
  if (length(y) != nrow(X))
    stop("bagoft.fast: 'data' has missing values in the model's variables.")
  if (ncol(Z) == 0L)
    stop("bagoft.fast: no covariate columns for the forest.")
  if (ncol(Z) > 5L && !requireNamespace("dcov", quietly = TRUE))
    stop("bagoft.fast: with more than five covariate columns the pre-selection ",
         "needs the 'dcov' package.")
  r <- .bagoft_core(X, y, Z, link = object$family$link, nsplits = nsplits, nsim = nsim,
                    ne = ne, ntree = ntree, Kmax = Kmax, nmin = nmin, mtry = mtry,
                    maxnodes = maxnodes)
  r$data.name <- paste(deparse(stats::formula(object)), collapse = " ")
  class(r) <- "bagoft_fast"
  r
}

#' @export
print.bagoft_fast <- function(x, ...) {
  st <- x$settings
  cat("\n\tBAGofT (Zhang, Ding and Yang 2023), fast implementation\n\n")
  cat("model: ", x$data.name, "\n", sep = "")
  cat(sprintf("nsplits = %d, nsim = %d, ne = %d, ntree = %d%s\n", st$nsplits, st$nsim,
              as.integer(st$ne), as.integer(st$ntree),
              if (isTRUE(st$preselected)) ", five columns pre-selected per split" else ""))
  if (!is.null(x$p.value)) {
    cat(sprintf("p-value = %s   (mean split p-value; median: %s, minimum: %s)\n",
                format(x$p.value), format(x$p.value2), format(x$p.value3)))
  }
  cat(sprintf("statistics: pmean = %s, pmedian = %s, pmin = %s (not p-values)\n\n",
              format(x$pmean, digits = 4), format(x$pmedian, digits = 4),
              format(x$pmin, digits = 4)))
  invisible(x)
}
