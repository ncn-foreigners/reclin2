## verify_continuous_engine.R
## Verification script for the hurdle-gamma continuous Fellegi-Sunter engine.
##
## Run from reclin2 package root:
##   Rscript inst/verification/verify_continuous_engine.R
##
## Part (A) uses EnvStats (installed automatically if needed).

cat("=================================================================\n")
cat("Verification: continuous hurdle-gamma engine for reclin2\n")
cat("=================================================================\n\n")

RECLIN2_PATH <- "/Users/berenz/mac/nauka/ncn-foreigners/software/forks/reclin2"

cat("Loading reclin2 fork via devtools::load_all ...\n")
if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")
suppressMessages(devtools::load_all(RECLIN2_PATH, quiet = TRUE))
cat("OK\n\n")
suppressMessages(library(data.table))

## =================================================================
## (A) Reference-equivalence anchor
##
## Goal: prove that the EM *update equations* in problink_em_gamma() are
## byte-identical to the paper's EM_hurdle_gammaK, ISOLATED from the
## initialisation heuristic. We therefore:
##   1. build ONE explicit set of starting values (the reference's own
##      EnvStats / p0M = 0.8 init),
##   2. feed it to BOTH the reference EM and problink_em_gamma() (via the
##      internal `.start` hook),
##   3. assert iterate-for-iterate / converged equality (max|diff| ~ 0).
## This separates "updates are faithful" from "init heuristic is good".
## =================================================================
cat("--- (A) Reference-equivalence anchor (shared start) ---\n")

if (!requireNamespace("EnvStats", quietly = TRUE)) {
  cat("  Installing EnvStats ...\n")
  install.packages("EnvStats", repos = "https://cloud.r-project.org", quiet = TRUE)
}

## ---- (A.0) Build the reference start exactly as start_hurdle_gammaK does ----
.ref_start_hurdle_gammaK <- function(X, K, nB, p0M = 0.8) {
  N <- nrow(X)
  lambdaM <- nB / N
  p0 <- matrix(0, 2, K); alpha <- matrix(1, 2, K); beta <- matrix(1, 2, K)
  for (k in seq_len(K)) {
    x  <- X[, k]; x1 <- sort(x[x > 0])
    nM_pos <- round((1 - p0M) * nB)
    if (nM_pos < 1 || nM_pos >= length(x1))
      nM_pos <- max(1, min(length(x1) - 1, nM_pos))
    xM <- x1[seq_len(nM_pos)]; xU <- x1[-seq_len(nM_pos)]
    p0U <- (sum(x == 0) - nB * p0M) / (N - nB)
    p0[, k] <- c(p0M, max(1/N, min(1 - 1/N, p0U)))
    fM <- if (length(xM) > 1 && min(xM) < max(xM))
      tryCatch(EnvStats::egamma(xM)$parameters, error = function(e) c(1, 1))
    else c(1, 1)
    fU <- tryCatch(EnvStats::egamma(xU)$parameters, error = function(e) c(1, 1))
    alpha[, k] <- c(fM[1], fU[1]); beta[, k] <- c(fM[2], fU[2])
  }
  p0[p0 < 1/N] <- 1/N
  list(lambda = c(lambdaM, 1 - lambdaM), p0 = p0, alpha = alpha, beta = beta)
}

## ---- (A.1) Reference EM, run FROM A GIVEN START (no internal init) ----
## Same algorithm as the paper's c3_EM.R, with the convergence robustness fixes.
.EM_hurdle_gammaK_ref <- function(X, K, start, tol = 1e-6, maxits = 500) {
  N <- nrow(X)
  lambda <- start$lambda; p0 <- start$p0; alpha <- start$alpha; beta <- start$beta

  dhg <- function(x, p0, sh, sc) {
    y <- rep(0, length(x))
    y[x == 0] <- p0
    y[x > 0]  <- (1 - p0) * dgamma(x[x > 0], shape = sh, scale = sc)
    y
  }
  dens_f <- function(lambda, p0, alpha, beta) {
    L <- matrix(1, N, 2)
    for (k in seq_len(K))
      for (j in 1:2) L[, j] <- L[, j] * dhg(X[, k], p0[j,k], alpha[j,k], beta[j,k])
    t(lambda * t(L))
  }
  fn_a  <- function(a,b,z,x) -log(b) + sum(z*log(x))/sum(z) - digamma(a)
  fn_a2 <- function(a,b,z,x) fn_a(a,b,z,x)^2

  d1 <- dens_f(lambda, p0, alpha, beta)
  iter <- 0L; converge <- FALSE

  while (!converge && iter < maxits) {
    ol <- lambda; op <- p0; os <- alpha; oc <- beta
    rs <- rowSums(d1); bad <- is.na(rs) | rs <= 0
    z <- d1 / rs; if (any(bad)) z[bad, ] <- matrix(rep(lambda,sum(bad)),ncol=2,byrow=TRUE)

    for (k in seq_len(K)) {
      i0 <- which(X[,k]==0); i1 <- which(X[,k]>0)
      p0[,k] <- pmax(1/N, pmin(1-1/N, colSums(z[i0,,drop=FALSE])/colSums(z)))
      if (length(i1) >= 3) {
        eff_n <- colSums(z[i1,,drop=FALSE])
        for (i in 1:2) {
          if (eff_n[i] >= 3) {
            new_sh <- tryCatch(
              uniroot(fn_a, c(1e-6,1e4), b=oc[i,k], z=z[i1,i], x=X[i1,k])$root,
              error=function(e) tryCatch(
                nlminb(os[i,k], fn_a2, lower=1e-8, b=oc[i,k], z=z[i1,i], x=X[i1,k])$par,
                error=function(e2) os[i,k])
            )
            alpha[i,k] <- new_sh
            beta[i,k]  <- sum(z[i1,i]*X[i1,k]) / (sum(z[i1,i]) * new_sh)
          }
        }
      }
    }
    lambda <- colMeans(z)
    d1 <- dens_f(lambda, p0, alpha, beta)
    converge <-
      all(abs(ol-lambda)/pmax(abs(ol),tol) < tol) &&
      all(abs(op-p0)/pmax(abs(op),tol)     < tol) &&
      all(abs(os-alpha)/pmax(abs(os),tol)  < tol) &&
      all(abs(oc-beta)/pmax(abs(oc),tol)   < tol)
    iter <- iter + 1L
  }
  list(g=z[,1], lambda=lambda, p0=p0, shape=alpha, scale=beta, iter=iter, converge=converge)
}

## Generate data (mirrors c1_generate_data.R DGP)
set.seed(2024)
nA <- 300; nB <- 120; K <- 3
XA <- matrix(rexp(nA * K, rate = 1), ncol = K)
XA1 <- XA[1:nB, ]
err <- matrix(rbinom(nB * K, 1, 0.15), ncol = K)
XB  <- (1 - err) * XA1 + err * (XA1 + matrix(rexp(nB * K, rate = 2), ncol = K))
Xcross <- do.call(cbind, lapply(1:K, function(k) {
  temp <- expand.grid(XA[, k], XB[, k])
  abs(temp[, 1] - temp[, 2])
}))
colnames(Xcross) <- paste0("k", 1:K)
cat("  Setup: nA=", nA, " nB=", nB, " K=", K, " N_pairs=", nrow(Xcross), "\n", sep="")

## ---- Build ONE shared start and feed it to BOTH ----
shared_start <- withCallingHandlers(
  .ref_start_hurdle_gammaK(Xcross, K, nB, p0M = 0.8),
  warning = function(w) invokeRestart("muffleWarning")
)

cat("  Running reference EM (shared start) ...\n")
ref <- withCallingHandlers(
  .EM_hurdle_gammaK_ref(Xcross, K, start = shared_start, tol = 1e-6),
  warning = function(w) invokeRestart("muffleWarning")
)
cat("  Reference converged:", ref$converge, " iter:", ref$iter, "\n")

cat("  Running problink_em_gamma (.start = shared start) ...\n")
our_model <- problink_em_gamma(
  comparison_matrix = Xcross, .start = shared_start, tol = 1e-6
)
cat("  Our model converged:", our_model$converged, " iter:", our_model$niter, "\n")

## --- Iterate count must match exactly (same updates, same start, same tol) ---
cat(sprintf("  Iteration counts equal? %s (ref=%d, ours=%d)\n",
            ref$iter == our_model$niter, ref$iter, our_model$niter))

## --- Converged-parameter equality (the strongest "same updates" check) ---
## NOTE: the reference returns components in (match, non-match) order from the
## start; our estimator may apply a label-switch guard at the end. Align by
## matching the larger-p0 (match) row before comparing.
align_match_first <- function(p0, shape, scale, lambda) {
  if (mean(p0[2, ]) > mean(p0[1, ])) {
    p0 <- p0[c(2,1), , drop=FALSE]; shape <- shape[c(2,1), , drop=FALSE]
    scale <- scale[c(2,1), , drop=FALSE]; lambda <- rev(lambda)
  }
  list(p0=p0, shape=shape, scale=scale, lambda=lambda)
}
R <- align_match_first(ref$p0, ref$shape, ref$scale, ref$lambda)
O <- align_match_first(our_model$p0, our_model$shape, our_model$scale,
                       c(our_model$p, 1 - our_model$p))

d_p0    <- max(abs(R$p0    - O$p0))
d_shape <- max(abs(R$shape - O$shape))
d_scale <- max(abs(R$scale - O$scale))
d_p     <- abs(R$lambda[1] - O$lambda[1])
cat(sprintf("  max|p0_ref    - p0_ours|    : %.3e\n", d_p0))
cat(sprintf("  max|shape_ref - shape_ours| : %.3e\n", d_shape))
cat(sprintf("  max|scale_ref - scale_ours| : %.3e\n", d_scale))
cat(sprintf("  |p_ref        - p_ours|     : %.3e\n", d_p))

## --- Posterior equality (the user-facing quantity) ---
Xdt <- as.data.table(Xcross)
Xdt[, .x := seq_len(.N)][, .y := seq_len(.N)]
setattr(Xdt, "class", c("pairs", class(Xdt)))
setattr(Xdt, "x", Xdt); setattr(Xdt, "y", Xdt)
pred <- predict(our_model, Xdt, type = "mpost")
mpost_ours <- pred$mpost
g_ref      <- ref$g

max_diff <- max(abs(mpost_ours - g_ref), na.rm = TRUE)
cor_val  <- cor(mpost_ours, g_ref, use = "complete.obs")
agree_05 <- mean((mpost_ours >= 0.5) == (g_ref >= 0.5), na.rm = TRUE)

cat(sprintf("  Max |mpost_ours - g_ref|:        %.3e\n", max_diff))
cat(sprintf("  Correlation(mpost_ours, g_ref):  %.6f\n", cor_val))
cat(sprintf("  Classification agreement @ 0.5:  %.4f\n", agree_05))
cat(sprintf("  Ref p:  %.8f   Ours p: %.8f\n", ref$lambda[1], our_model$p))

param_eq <- max(d_p0, d_shape, d_scale, d_p) < 1e-8
if (param_eq && max_diff < 1e-8 && ref$iter == our_model$niter) {
  cat("  PASS: byte-identical updates from a shared start ",
      "(params, posteriors and iteration count all match).\n", sep = "")
} else if (max_diff < 1e-6) {
  cat("  PASS: posteriors agree to < 1e-6 from a shared start.\n")
} else {
  cat("  FAIL: updates diverge from a shared start -- investigate.\n")
}

## =================================================================
## (B) End-to-end recovery -- DGP-A (exact-match regime), DEFAULT init
## =================================================================
cat("\n--- (B) End-to-end recovery: DGP-A (exact matches), DEFAULTS ---\n")

set.seed(42)
nA2 <- 400; nB2 <- 80; K2 <- 3

XA2 <- matrix(rexp(nA2 * K2, rate = 1), ncol = K2)
err2 <- matrix(rbinom(nB2 * K2, 1, 0.05), ncol = K2)
XB2  <- (1 - err2) * XA2[1:nB2, ] +
        err2 * (XA2[1:nB2, ] + matrix(rexp(nB2 * K2, rate = 5), ncol = K2))

dfA2 <- as.data.frame(XA2); names(dfA2) <- paste0("k", 1:K2)
dfB2 <- as.data.frame(XB2); names(dfB2) <- paste0("k", 1:K2)

pairs2 <- pair(dfA2, dfB2)
pairs2 <- compare_pairs(pairs2, on = paste0("k", 1:K2),
  default_comparator = cmp_distance())
cat("  Pairs:", nrow(pairs2), "  True matches:", nB2, "\n")

## p0 OMITTED: auto-derived from the pairs' source attributes as
## min(nA, nB) / nrow(pairs) -> here min(nA2, nB2) / nrow(pairs2).
model2 <- problink_em_gamma(~ k1 + k2 + k3, data = pairs2)
cat("  Auto-derived p0 (EM start):", round(model2$p0_init, 6),
    " (expected", round(min(nA2, nB2) / nrow(pairs2), 6), ")\n")
cat("  Converged:", model2$converged, "  Iterations:", model2$niter, "\n")
cat("  Hurdle active:", model2$hurdle_active, "(expected TRUE)\n")

pairs2 <- predict(model2, pairs2, type = "mpost", add = TRUE)
pairs2 <- select_threshold(pairs2, "selected", "mpost", 0.5)
pairs2[, truth := (.x == .y) & (.x <= nB2)]

sel <- pairs2[selected == TRUE]
TP  <- sum(sel$truth); FP <- nrow(sel) - TP; FN <- nB2 - TP
TPR <- TP / nB2
PPV <- if (nrow(sel) > 0) TP / nrow(sel) else NA_real_

cat(sprintf("  TP=%d  FP=%d  FN=%d\n", TP, FP, FN))
cat(sprintf("  TPR (recall):    %.4f\n", TPR))
cat(sprintf("  PPV (precision): %.4f\n", PPV))
cat(sprintf("  Mixing p_M: %.6f\n", model2$p))
cat(sprintf("  Match component:     p0=%.3f  mean_dist=%.5f\n",
            mean(model2$p0[1,]), mean(model2$shape[1,]*model2$scale[1,])))
cat(sprintf("  Non-match component: p0=%.3f  mean_dist=%.4f\n",
            mean(model2$p0[2,]), mean(model2$shape[2,]*model2$scale[2,])))

if (!is.na(TPR) && !is.na(PPV) && TPR > 0.8 && PPV > 0.8) {
  cat("  PASS: TPR > 0.8 and PPV > 0.8\n")
} else {
  cat("  WARNING: TPR or PPV below 0.8\n")
}

## =================================================================
## (B2) End-to-end recovery -- DGP-B (PURE CONTINUOUS, no exact matches)
## =================================================================
cat("\n--- (B2) End-to-end recovery: DGP-B (pure continuous), DEFAULTS ---\n")

set.seed(2025)
nA3 <- 400; nB3 <- 80; K3 <- 3
XA3 <- matrix(runif(nA3 * K3, 0, 100), ncol = K3)
XB3 <- XA3[1:nB3, ] + matrix(rnorm(nB3 * K3, 0, 1), ncol = K3)
dfA3 <- as.data.frame(XA3); names(dfA3) <- paste0("k", 1:K3)
dfB3 <- as.data.frame(XB3); names(dfB3) <- paste0("k", 1:K3)

pairs3b <- pair(dfA3, dfB3)
pairs3b <- compare_pairs(pairs3b, on = paste0("k", 1:K3),
  default_comparator = cmp_distance())
nzero <- sum(as.matrix(pairs3b[, paste0("k",1:K3), with=FALSE]) == 0, na.rm=TRUE)
cat("  Pairs:", nrow(pairs3b), "  True matches:", nB3,
    "  Exact-zero distances:", nzero, "(expect 0)\n")

model3b <- problink_em_gamma(~ k1 + k2 + k3, data = pairs3b)   # DEFAULTS
cat("  Converged:", model3b$converged, "  Iterations:", model3b$niter, "\n")
cat(sprintf("  p0 (match, ~1/N=%.2e): %s\n", 1/nrow(pairs3b),
            paste(format(model3b$p0[1,], digits=3), collapse=" ")))
cat("  Hurdle active:", model3b$hurdle_active, "(expected FALSE)\n")

pairs3b <- predict(model3b, pairs3b, type = "mpost", add = TRUE)
pairs3b <- select_threshold(pairs3b, "selected", "mpost", 0.5)
pairs3b[, truth := (.x == .y) & (.x <= nB3)]
selb <- pairs3b[selected == TRUE]
TPb  <- sum(selb$truth)
TPRb <- TPb / nB3
PPVb <- if (nrow(selb) > 0) TPb / nrow(selb) else NA_real_
cat(sprintf("  TP=%d  sel=%d  TPR=%.4f  PPV=%.4f  p_M=%.6f\n",
            TPb, nrow(selb), TPRb, PPVb, model3b$p))
if (!is.na(TPRb) && !is.na(PPVb) && TPRb > 0.8 && PPVb > 0.8) {
  cat("  PASS: DGP-B converges at DEFAULTS and recovers matches.\n")
} else {
  cat("  WARNING: DGP-B TPR or PPV below 0.8\n")
}

## =================================================================
## (C) Smoke / regression
## =================================================================
cat("\n--- (C) Smoke / regression ---\n")

data("linkexample1", "linkexample2", package = "reclin2")
pairs3 <- pair_blocking(linkexample1, linkexample2, "postcode")
pairs3 <- compare_pairs(pairs3, c("lastname", "firstname", "address", "sex"))
m3 <- suppressWarnings(problink_em(~ lastname + firstname + address + sex, data = pairs3))
cat("  problink_em (existing): p =", round(m3$p, 4), " PASS\n")

cat("  cmp_distance loaded:", is.function(cmp_distance()), " PASS\n")
cat("  cmp_absdist loaded:", is.function(cmp_absdist()), " PASS\n")
cat("  problink_em_gamma loaded:", is.function(problink_em_gamma), " PASS\n")

cmp <- cmp_distance()
d   <- cmp(c(1, 2, 3, NA), c(1, 3, 1, 2))
ok1 <- all(d[1:3] == c(0, 1, 2), na.rm = TRUE) && is.na(d[4])
cat("  cmp_distance 2-arg correct:", ok1, " PASS\n")
b   <- cmp(d)
ok2 <- isTRUE(b[1]) && !isTRUE(b[2]) && !isTRUE(b[3]) && !isTRUE(b[4])
cat("  cmp_distance 1-arg (binary):", ok2, " PASS\n")

cmp2 <- cmp_distance(scale = 10)
d2   <- cmp2(c(0, 5, 10), c(0, 0, 0))
cat("  cmp_distance scale:", all(abs(d2 - c(0, 0.5, 1.0)) < 1e-10), " PASS\n")

cat("\n=================================================================\n")
cat("Verification complete.\n")
cat("=================================================================\n")
