# data-raw/make_gof_demo.R
# Reproducibly builds data/gof_demo.rda: a synthetic binary dataset whose TRUE
# data-generating process contains a smooth quadratic-in-age term. A model that
# omits it (outcome ~ age + bmi + sex + treatment) is therefore subtly
# misspecified through a low-dimensional calibration distortion -- exactly the
# kind of misfit the directed EF / EDGE test is built to detect, while classical
# omnibus tests are weaker. Re-run with:  Rscript data-raw/make_gof_demo.R
set.seed(2026)

n <- 800L
age       <- round(runif(n, 20, 70), 1)          # years
bmi       <- round(rnorm(n, 27, 4), 1)           # kg/m^2
sex       <- rbinom(n, 1, 0.5)                   # 0 = female, 1 = male
treatment <- rbinom(n, 1, 0.4)                   # 0 = control, 1 = treated

za  <- (age - 45) / 14                            # standardized age
zb  <- (bmi - 27) / 4                             # standardized bmi

# TRUE linear predictor -- note the -0.7 * za^2 curvature that a linear-in-age
# model cannot capture.
eta <- -0.6 + 0.8 * za - 0.7 * za^2 + 0.5 * zb + 0.4 * sex - 0.3 * treatment
p   <- plogis(eta)
outcome <- rbinom(n, 1, p)

gof_demo <- data.frame(outcome, age, bmi, sex, treatment)

cat(sprintf("n = %d, event rate = %.3f\n", n, mean(gof_demo$outcome)))

dir.create("data", showWarnings = FALSE)
save(gof_demo, file = "data/gof_demo.rda", version = 2, compress = "xz")
cat("wrote data/gof_demo.rda\n")

# ---------------------------------------------------------------------------
# gof_demo_grouped: the SAME data-generating process and seed discipline, but
# the covariates are coarsened BEFORE the linear predictor is computed -- age is
# rounded to 10-year bins and bmi to integers -- so the recorded covariates ARE
# the covariates the outcome was generated from, and many observations share a
# covariate pattern. This is the replicated-pattern regime in which the sparse
# (per-observation) and grouped (per-pattern) forms of the pattern-sensitive
# tests (Pearson, deviance, McCullagh) can disagree on the same fitted model.
# (5-year bins were tried first; at that setting the sparse and grouped p-values
# already differ visibly but no verdict flips at alpha = .05, so the shipped
# dataset uses the coarser 10-year bins, where the deviance test flips.)
# ---------------------------------------------------------------------------
set.seed(2026)

n <- 800L
age_raw   <- round(runif(n, 20, 70), 1)
bmi_raw   <- round(rnorm(n, 27, 4), 1)
sex       <- rbinom(n, 1, 0.5)                   # 0 = female, 1 = male
treatment <- rbinom(n, 1, 0.4)                   # 0 = control, 1 = treated

age <- 10 * round(age_raw / 10)                  # 10-year bins: 20, 30, ..., 70
bmi <- round(bmi_raw)                            # integer BMI

za  <- (age - 45) / 14                            # standardized BINNED age
zb  <- (bmi - 27) / 4                             # standardized BINNED bmi

# Same true linear predictor as gof_demo (with the -0.7 * za^2 curvature),
# computed from the binned covariates.
eta <- -0.6 + 0.8 * za - 0.7 * za^2 + 0.5 * zb + 0.4 * sex - 0.3 * treatment
p   <- plogis(eta)
outcome <- rbinom(n, 1, p)

gof_demo_grouped <- data.frame(outcome, age, bmi, sex, treatment)

n_pat <- nrow(unique(gof_demo_grouped[, c("age", "bmi", "sex", "treatment")]))
cat(sprintf("grouped: n = %d, event rate = %.3f, covariate patterns = %d\n",
            n, mean(gof_demo_grouped$outcome), n_pat))

save(gof_demo_grouped, file = "data/gof_demo_grouped.rda",
     version = 2, compress = "xz")
cat("wrote data/gof_demo_grouped.rda\n")
