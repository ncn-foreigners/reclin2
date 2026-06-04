
library(reclin2)
library(data.table)
source("helpers.R")

## ========================================================
## 1. problink_em_gamma: comparison_matrix interface
## ========================================================
set.seed(42)
nA <- 200; nB <- 40; K <- 3
XA <- matrix(rexp(nA * K, rate = 1), ncol = K)
XA1 <- XA[1:nB, ]
err <- matrix(rbinom(nB * K, 1, 0.05), ncol = K)
XB <- (1 - err) * XA1 + err * (XA1 + matrix(rexp(nB * K, rate = 4), ncol = K))

## Build full cross-product comparison matrix
Xcross <- do.call(cbind, lapply(1:K, function(k) {
  temp <- expand.grid(XA[, k], XB[, k])
  abs(temp[, 1] - temp[, 2])
}))
colnames(Xcross) <- paste0("k", 1:K)

model_mat <- problink_em_gamma(
  comparison_matrix = Xcross,
  p0 = nB / nrow(Xcross),
  p0M = 0.8,
  tol = 1e-6,
  maxits = 500
)

stopifnot(inherits(model_mat, "problink_em_gamma"))
stopifnot(model_mat$converged)
stopifnot(model_mat$p > 0 && model_mat$p < 1)
stopifnot(nrow(model_mat$p0) == 2 && ncol(model_mat$p0) == K)
stopifnot(nrow(model_mat$shape) == 2 && ncol(model_mat$shape) == K)
stopifnot(nrow(model_mat$scale) == 2 && ncol(model_mat$scale) == K)
stopifnot(all(model_mat$shape > 0))
stopifnot(all(model_mat$scale > 0))
cat("  [1] comparison_matrix interface: PASS\n")

## ========================================================
## 2. problink_em_gamma: formula + data interface (pairs)
## ========================================================
dfA <- as.data.frame(XA); names(dfA) <- paste0("k", 1:K)
dfB <- as.data.frame(XB); names(dfB) <- paste0("k", 1:K)
dfA$id <- seq_len(nA)
dfB$id <- seq_len(nB)

pairs <- pair(dfA, dfB)
pairs <- compare_pairs(pairs, on = paste0("k", 1:K),
  default_comparator = cmp_distance())

## p0 is OMITTED: it should be auto-derived from the pairs' source attributes
## as min(nA, nB) / nrow(pairs).
model_pairs <- problink_em_gamma(~ k1 + k2 + k3, data = pairs)

stopifnot(inherits(model_pairs, "problink_em_gamma"))
stopifnot(model_pairs$converged)
stopifnot(identical(model_pairs$by, paste0("k", 1:K)))
## Auto-derived p0 default: min(nA, nB) / N
expect_equal(model_pairs$p0_init, min(nA, nB) / nrow(pairs))
cat("  [2] formula+data interface (auto p0): PASS\n")

## ========================================================
## 2b. Explicit p0 override is still honoured
## ========================================================
model_p0 <- problink_em_gamma(~ k1 + k2 + k3, data = pairs, p0 = 0.123)
stopifnot(inherits(model_p0, "problink_em_gamma"))
expect_equal(model_p0$p0_init, 0.123)
## Invalid p0 must be rejected
expect_error(problink_em_gamma(~ k1 + k2 + k3, data = pairs, p0 = 0))
expect_error(problink_em_gamma(~ k1 + k2 + k3, data = pairs, p0 = 1))
cat("  [2b] explicit p0 override: PASS\n")

## ========================================================
## 3. predict.problink_em_gamma: all types
## ========================================================
res_w <- predict(model_pairs, pairs, type = "weights", add = FALSE)
stopifnot(is.data.table(res_w))
stopifnot(all(c(".x", ".y", "weights") %in% names(res_w)))

res_m <- predict(model_pairs, pairs, type = "mpost", add = FALSE)
stopifnot(all(c(".x", ".y", "mpost") %in% names(res_m)))
stopifnot(all(res_m$mpost >= 0 & res_m$mpost <= 1, na.rm = TRUE))

res_p <- predict(model_pairs, pairs, type = "probs", add = FALSE)
stopifnot(all(c("mprob", "uprob", "mpost", "upost") %in% names(res_p)))
expect_equal(res_p$mpost + res_p$upost, rep(1, nrow(res_p)), attributes = FALSE)

res_all <- predict(model_pairs, pairs, type = "all", add = FALSE)
stopifnot(all(c("weight", "mprob", "uprob", "mpost", "upost") %in% names(res_all)))

# add = TRUE
pairs_aug <- predict(model_pairs, pairs, type = "mpost", add = TRUE)
stopifnot("mpost" %in% names(pairs_aug))
stopifnot(nrow(pairs_aug) == nrow(pairs))
cat("  [3] predict types: PASS\n")

## ========================================================
## 4. End-to-end TPR/PPV
## ========================================================
pairs[, truth := (.x == .y) & (.x <= nB)]
pairs <- predict(model_pairs, pairs, type = "mpost", add = TRUE)
pairs <- select_threshold(pairs, "selected", "mpost", 0.5)

sel <- pairs[selected == TRUE]
TP <- sum(sel$truth)
TPR <- TP / nB
PPV <- if (nrow(sel) > 0) TP / nrow(sel) else 0

cat(sprintf("  [4] TPR=%.4f  PPV=%.4f", TPR, PPV))
## With the auto-derived p0 the model must still recover matches well.
stopifnot(TPR > 0.8 && PPV > 0.8)
cat("  PASS\n")

## ========================================================
## 5. print/summary don't error
## ========================================================
capture.output(print(model_pairs))
capture.output(summary(model_pairs))
cat("  [5] print/summary: PASS\n")

## ========================================================
## 6. NA robustness: rows with NA get neutral treatment
## ========================================================
pairs_na <- copy(pairs)
pairs_na[1, k1 := NA_real_]
res_na <- predict(model_pairs, pairs_na, type = "mpost", add = FALSE)
stopifnot(!is.na(res_na$mpost[1]))  # NA should be handled, not propagate
cat("  [6] NA robustness: PASS\n")

## ========================================================
## 7. Label-switching: match component has higher p0
## ========================================================
stopifnot(mean(model_pairs$p0[1, ]) >= mean(model_pairs$p0[2, ]) - 0.05)
cat("  [7] Label-switching guard: PASS\n")

## ========================================================
## 8. Initialisation unit test (TIER 1 / 5.3):
##    p0M_init and the non-match hurdle start (p0U_init) are valid
##    probabilities in [0, 1] in BOTH regimes -- never negative -- and
##    ~1/N when there are no exact-zero distances.
## ========================================================
init_hurdle <- reclin2:::.init_hurdle_gammaK

## (a) Hurdle regime: many exact zeros -> derived p0M well inside (0,1)
set.seed(7)
Nh <- 500
Xh <- cbind(
  k1 = c(rep(0, 120), abs(rexp(Nh - 120))),
  k2 = c(rep(0,  90), abs(rexp(Nh -  90)))
)
zk_h     <- colMeans(Xh == 0)
lam_h    <- 0.25
p0M_h    <- reclin2:::clamp01_(zk_h / lam_h, lo = 1/Nh)
init_h   <- init_hurdle(Xh, K = 2, p0_init = lam_h, p0M_init = p0M_h,
                        N = Nh, use_envstats = FALSE)
stopifnot(all(p0M_h  >= 0 & p0M_h  <= 1))      # derived p0M valid
stopifnot(all(init_h$p0 >= 0 & init_h$p0 <= 1)) # all hurdle starts valid
stopifnot(all(init_h$p0[2, ] >= 1/Nh))          # non-match start (p0U) >= 1/N
stopifnot(all(is.finite(init_h$alpha) & init_h$alpha > 0))
stopifnot(all(is.finite(init_h$beta)  & init_h$beta  > 0))

## (b) Pure-continuous regime: NO exact zeros -> p0M ~ 1/N, p0U = 1/N,
##     and never negative (the old (#zeros - nB*p0M)/(N-nB) formula went < 0).
set.seed(8)
Nc <- 500
Xc <- cbind(k1 = abs(rexp(Nc)) + 0.01,
            k2 = abs(rexp(Nc)) + 0.01)   # strictly positive: no zeros
stopifnot(sum(Xc == 0) == 0)
zk_c   <- colMeans(Xc == 0)              # all zero
lam_c  <- 0.05
p0M_c  <- reclin2:::clamp01_(zk_c / lam_c, lo = 1/Nc)
init_c <- init_hurdle(Xc, K = 2, p0_init = lam_c, p0M_init = p0M_c,
                      N = Nc, use_envstats = FALSE)
stopifnot(all(p0M_c >= 0 & p0M_c <= 1))
stopifnot(all(abs(p0M_c - 1/Nc) < 1e-12))        # ~1/N when no zeros
stopifnot(all(init_c$p0 >= 0 & init_c$p0 <= 1))   # never negative
stopifnot(all(abs(init_c$p0 - 1/Nc) < 1e-12))     # both rows pinned to 1/N
stopifnot(all(init_c$lambda >= 0 & init_c$lambda <= 1))
cat("  [8] Initialisation unit test (both regimes): PASS\n")

## ========================================================
## 9. Reference-equivalence (TIER 5.2): from a SHARED explicit start the
##    EM updates equal a self-contained reference EM (max|diff| ~ 0). This
##    isolates "updates are faithful" from the init heuristic via the
##    internal `.start` hook. (Method-of-moments init -> no EnvStats needed.)
## ========================================================
set.seed(99)
nA9 <- 150; nB9 <- 60; K9 <- 2
XA9 <- matrix(rexp(nA9 * K9, rate = 1), ncol = K9)
err9 <- matrix(rbinom(nB9 * K9, 1, 0.15), ncol = K9)
XB9  <- (1 - err9) * XA9[1:nB9, ] +
        err9 * (XA9[1:nB9, ] + matrix(rexp(nB9 * K9, rate = 2), ncol = K9))
X9 <- do.call(cbind, lapply(1:K9, function(k) {
  t9 <- expand.grid(XA9[, k], XB9[, k]); abs(t9[, 1] - t9[, 2])
}))
colnames(X9) <- paste0("k", 1:K9)
N9 <- nrow(X9)

## Build an explicit start (method-of-moments, p0M=0.8) and floor p0 at 1/N
## exactly as the reference does.
mk_start <- function(X, K, nB, p0M = 0.8) {
  N <- nrow(X); lambdaM <- nB / N
  p0 <- matrix(0,2,K); alpha <- matrix(1,2,K); beta <- matrix(1,2,K)
  mom <- function(x){x<-x[x>0]; m<-mean(x); v<-var(x)
    if(!is.finite(v)||v<=0) c(1,1) else c(m^2/v, v/m)}
  for (k in seq_len(K)) {
    x <- X[,k]; x1 <- sort(x[x>0])
    nM_pos <- max(1L, min(length(x1)-1L, round((1-p0M)*nB)))
    fM <- mom(x1[seq_len(nM_pos)]); fU <- mom(x1[-seq_len(nM_pos)])
    p0U <- (sum(x==0) - nB*p0M)/(N-nB)
    p0[,k] <- c(p0M, max(1/N, min(1-1/N, p0U)))
    alpha[,k] <- c(fM[1], fU[1]); beta[,k] <- c(fM[2], fU[2])
  }
  p0[p0 < 1/N] <- 1/N
  list(lambda=c(lambdaM, 1-lambdaM), p0=p0, alpha=alpha, beta=beta)
}

## Self-contained reference EM run FROM A GIVEN START.
ref_em <- function(X, K, start, tol = 1e-6, maxits = 500) {
  N <- nrow(X)
  lambda <- start$lambda; p0 <- start$p0; alpha <- start$alpha; beta <- start$beta
  dhg <- function(x,p0,sh,sc){y<-rep(0,length(x)); y[x==0]<-p0
    y[x>0]<-(1-p0)*dgamma(x[x>0],shape=sh,scale=sc); y}
  dens_f <- function(lambda,p0,alpha,beta){L<-matrix(1,N,2)
    for(k in seq_len(K)) for(j in 1:2)
      L[,j]<-L[,j]*dhg(X[,k],p0[j,k],alpha[j,k],beta[j,k]); t(lambda*t(L))}
  fn_a <- function(a,b,z,x) -log(b)+sum(z*log(x))/sum(z)-digamma(a)
  fn_a2<- function(a,b,z,x) fn_a(a,b,z,x)^2
  d1 <- dens_f(lambda,p0,alpha,beta); iter<-0L; conv<-FALSE
  while(!conv && iter<maxits){
    ol<-lambda; op<-p0; os<-alpha; oc<-beta
    rs<-rowSums(d1); bad<-is.na(rs)|rs<=0; z<-d1/rs
    if(any(bad)) z[bad,]<-matrix(rep(lambda,sum(bad)),ncol=2,byrow=TRUE)
    for(k in seq_len(K)){
      i0<-which(X[,k]==0); i1<-which(X[,k]>0)
      p0[,k]<-pmax(1/N,pmin(1-1/N,colSums(z[i0,,drop=FALSE])/colSums(z)))
      if(length(i1)>=3){en<-colSums(z[i1,,drop=FALSE])
        for(i in 1:2) if(en[i]>=3){
          ns<-tryCatch(uniroot(fn_a,c(1e-6,1e4),b=oc[i,k],z=z[i1,i],x=X[i1,k])$root,
            error=function(e) tryCatch(nlminb(os[i,k],fn_a2,lower=1e-8,
              b=oc[i,k],z=z[i1,i],x=X[i1,k])$par, error=function(e2) os[i,k]))
          alpha[i,k]<-ns; beta[i,k]<-sum(z[i1,i]*X[i1,k])/(sum(z[i1,i])*ns)}}
    }
    lambda<-colMeans(z); d1<-dens_f(lambda,p0,alpha,beta)
    conv<- all(abs(ol-lambda)/pmax(abs(ol),tol)<tol) &&
           all(abs(op-p0)/pmax(abs(op),tol)<tol) &&
           all(abs(os-alpha)/pmax(abs(os),tol)<tol) &&
           all(abs(oc-beta)/pmax(abs(oc),tol)<tol)
    iter<-iter+1L
  }
  list(g=z[,1], p0=p0, shape=alpha, scale=beta, lambda=lambda, iter=iter, conv=conv)
}

start9 <- mk_start(X9, K9, nB9, p0M = 0.8)
ref9 <- suppressWarnings(ref_em(X9, K9, start9, tol = 1e-6))
our9 <- problink_em_gamma(comparison_matrix = X9, .start = start9, tol = 1e-6)

## .start forces a single deterministic run.
stopifnot(our9$nstart == 1L)
## Iteration counts must match (same updates, same start, same tol).
stopifnot(ref9$iter == our9$niter)

## Align match-first (our estimator may apply a label-switch at the very end).
al <- function(p0, sh, sc, lam) {
  if (mean(p0[2, ]) > mean(p0[1, ])) {
    list(p0 = p0[2:1, , drop = FALSE], sh = sh[2:1, , drop = FALSE],
         sc = sc[2:1, , drop = FALSE], lam = rev(lam))
  } else list(p0 = p0, sh = sh, sc = sc, lam = lam)
}
R9 <- al(ref9$p0, ref9$shape, ref9$scale, ref9$lambda)
O9 <- al(our9$p0, our9$shape, our9$scale, c(our9$p, 1 - our9$p))
stopifnot(max(abs(R9$p0 - O9$p0))   < 1e-8)
stopifnot(max(abs(R9$sh - O9$sh))   < 1e-8)
stopifnot(max(abs(R9$sc - O9$sc))   < 1e-8)
stopifnot(abs(R9$lam[1] - O9$lam[1]) < 1e-8)

## Posterior equality through predict().
X9dt <- as.data.table(X9)
X9dt[, .x := seq_len(.N)][, .y := seq_len(.N)]
setattr(X9dt, "class", c("pairs", class(X9dt)))
setattr(X9dt, "x", X9dt); setattr(X9dt, "y", X9dt)
mpost9 <- predict(our9, X9dt, type = "mpost")$mpost
stopifnot(max(abs(mpost9 - ref9$g), na.rm = TRUE) < 1e-8)
cat("  [9] Reference-equivalence (shared start, max|diff|~0): PASS\n")

## ========================================================
## 10. Pure-continuous regression (TIER 5.1): DGP-B at DEFAULTS converges,
##     recovers matches (TPR/PPV thresholds), and p0 ~ 1/N.
## ========================================================
set.seed(2025)
nA10 <- 400; nB10 <- 80; K10 <- 3
XA10 <- matrix(runif(nA10 * K10, 0, 100), ncol = K10)
XB10 <- XA10[1:nB10, ] + matrix(rnorm(nB10 * K10, 0, 1), ncol = K10)
dfA10 <- as.data.frame(XA10); names(dfA10) <- paste0("k", 1:K10)
dfB10 <- as.data.frame(XB10); names(dfB10) <- paste0("k", 1:K10)
pairs10 <- pair(dfA10, dfB10)
pairs10 <- compare_pairs(pairs10, on = paste0("k", 1:K10),
                         default_comparator = cmp_distance())
## No exact-zero distances in this regime.
stopifnot(sum(as.matrix(pairs10[, paste0("k", 1:K10), with = FALSE]) == 0,
              na.rm = TRUE) == 0)

model10 <- problink_em_gamma(~ k1 + k2 + k3, data = pairs10)  # DEFAULTS
stopifnot(model10$converged)
## p0 must collapse to ~1/N (hurdle inactive) -- the correct MLE here.
N10 <- nrow(pairs10)
stopifnot(all(abs(model10$p0 - 1/N10) < 1e-6))
stopifnot(isFALSE(model10$hurdle_active))

pairs10 <- predict(model10, pairs10, type = "mpost", add = TRUE)
pairs10 <- select_threshold(pairs10, "selected", "mpost", 0.5)
pairs10[, truth := (.x == .y) & (.x <= nB10)]
sel10 <- pairs10[selected == TRUE]
TP10  <- sum(sel10$truth)
TPR10 <- TP10 / nB10
PPV10 <- if (nrow(sel10) > 0) TP10 / nrow(sel10) else 0
cat(sprintf("  [10] DGP-B defaults: TPR=%.4f PPV=%.4f p0~1/N PASS-check",
            TPR10, PPV10))
stopifnot(TPR10 > 0.8 && PPV10 > 0.8)
cat(" PASS\n")

## ========================================================
## 11. nstart determinism + diagnostics (TIER 4):
##     nstart=1 deterministic; nstart>1 reproducible under set.seed and
##     never worse in log-likelihood; print reports diagnostics.
## ========================================================
m_a <- problink_em_gamma(~ k1 + k2 + k3, data = pairs10, nstart = 1)
m_b <- problink_em_gamma(~ k1 + k2 + k3, data = pairs10, nstart = 1)
expect_equal(m_a$p, m_b$p)            # nstart=1 is deterministic

set.seed(123)
m_ms1 <- problink_em_gamma(~ k1 + k2 + k3, data = pairs10, nstart = 3)
set.seed(123)
m_ms2 <- problink_em_gamma(~ k1 + k2 + k3, data = pairs10, nstart = 3)
expect_equal(m_ms1$loglik, m_ms2$loglik)   # reproducible under set.seed
stopifnot(m_ms1$nstart == 3L)
## Multistart keeps the best: never worse than the single deterministic start.
stopifnot(m_ms1$loglik >= m_a$loglik - 1e-6)
## RNG state is restored (set_seed/on.exit): the user's stream is untouched.
set.seed(321); ref_stream <- runif(3)
set.seed(321); invisible(problink_em_gamma(~ k1 + k2 + k3, data = pairs10, nstart = 3))
post_stream <- runif(3)
expect_equal(ref_stream, post_stream)
## print() reports the diagnostics (hurdle flag, p0M init, iterations).
out <- capture.output(print(model10))
stopifnot(any(grepl("Hurdle point mass", out)))
stopifnot(any(grepl("Initial exact-match prob", out)))
expect_error(problink_em_gamma(~ k1 + k2 + k3, data = pairs10, nstart = 0))
cat("  [11] nstart determinism + diagnostics: PASS\n")

## ========================================================
## 12. PARTIAL MISSINGNESS (FIX 1): per-field NA is a NEUTRAL factor, so a
##     pair with one NA field + strong observed fields is NOT reduced to the
##     prior. loglik finite; recovery sane; nstart>1 works under missingness.
## ========================================================
set.seed(202)
nA12 <- 300; nB12 <- 60; K12 <- 3
XA12 <- matrix(rexp(nA12 * K12, rate = 1), ncol = K12)
err12 <- matrix(rbinom(nB12 * K12, 1, 0.05), ncol = K12)
XB12 <- (1 - err12) * XA12[1:nB12, ] +
        err12 * (XA12[1:nB12, ] + matrix(rexp(nB12 * K12, 4), ncol = K12))
X12 <- do.call(cbind, lapply(1:K12, function(k) {
  t <- expand.grid(XA12[, k], XB12[, k]); abs(t[, 1] - t[, 2])
}))
colnames(X12) <- paste0("k", 1:K12)
N12 <- nrow(X12)
## Punch NA into ~5% of the cells of each field (per-field, scattered).
set.seed(203)
X12na <- X12
for (k in 1:K12) {
  ii <- sample(seq_len(N12), round(0.05 * N12))
  X12na[ii, k] <- NA_real_
}
stopifnot(anyNA(X12na))

m12 <- problink_em_gamma(comparison_matrix = X12na, p0 = nB12 / N12, p0M = 0.8)
stopifnot(m12$converged)
stopifnot(is.finite(m12$loglik))                # loglik finite under NA
stopifnot(all(is.finite(m12$shape)), all(m12$shape > 0))
stopifnot(all(is.finite(m12$scale)), all(m12$scale > 0))

## Construct a pair that is an EXACT match on 2 observed fields but NA on the
## 3rd; its mpost must reflect the observed (strongly matching) fields, i.e. be
## clearly ABOVE the prior p, NOT collapsed to it.
probe <- matrix(c(0, 0, NA_real_), nrow = 1)   # 2 exact matches + 1 NA
colnames(probe) <- paste0("k", 1:K12)
pdt <- as.data.table(probe)
pdt[, .x := 1L][, .y := 1L]
setattr(pdt, "class", c("pairs", class(pdt)))
setattr(pdt, "x", pdt); setattr(pdt, "y", pdt)
probe_mpost <- predict(m12, pdt, type = "mpost")$mpost
stopifnot(is.finite(probe_mpost))
stopifnot(probe_mpost > m12$p + 0.1)            # NOT reduced to the prior

## nstart > 1 under missingness: loglik finite and best is selected (>= single).
m12a <- problink_em_gamma(comparison_matrix = X12na, p0 = nB12 / N12,
                          p0M = 0.8, nstart = 1)
set.seed(7)
m12b <- problink_em_gamma(comparison_matrix = X12na, p0 = nB12 / N12,
                          p0M = 0.8, nstart = 4)
stopifnot(is.finite(m12b$loglik))
stopifnot(m12b$nstart == 4L)
stopifnot(m12b$loglik >= m12a$loglik - 1e-6)    # multistart never worse
cat("  [12] partial missingness (neutral NA, nstart>1, loglik finite): PASS\n")

## ========================================================
## 13. PREDICT UNDERFLOW (FIX 6): with many fields the raw mprob/uprob products
##     underflow to 0, but mpost stays finite in (0,1) and equals
##     plogis(weights + qlogis(p)) (computed on the log scale).
## ========================================================
set.seed(204)
K13 <- 80                              # many fields -> product underflows
XA13 <- matrix(rexp(60 * K13, rate = 1), ncol = K13)
err13 <- matrix(rbinom(20 * K13, 1, 0.1), ncol = K13)
XB13 <- (1 - err13) * XA13[1:20, ] +
        err13 * (XA13[1:20, ] + matrix(rexp(20 * K13, 4), ncol = K13))
X13 <- do.call(cbind, lapply(1:K13, function(k) {
  t <- expand.grid(XA13[, k], XB13[, k]); abs(t[, 1] - t[, 2])
}))
colnames(X13) <- paste0("k", 1:K13)
m13 <- problink_em_gamma(comparison_matrix = X13, p0 = 20 / nrow(X13),
                         p0M = 0.8, maxits = 50)
X13dt <- as.data.table(X13)
X13dt[, .x := seq_len(.N)][, .y := seq_len(.N)]
setattr(X13dt, "class", c("pairs", class(X13dt)))
setattr(X13dt, "x", X13dt); setattr(X13dt, "y", X13dt)
pr13 <- predict(m13, X13dt, type = "probs")
## At least some rows must have underflowed the raw products to exactly 0.
stopifnot(any(pr13$mprob == 0) || any(pr13$uprob == 0))
## ... yet mpost is always finite and in [0, 1] (never NaN from 0/0).
stopifnot(all(is.finite(pr13$mpost)))
stopifnot(all(pr13$mpost >= 0 & pr13$mpost <= 1))
stopifnot(all(is.finite(pr13$upost)))
## ... and mpost == plogis(weights + qlogis(p)) exactly.
w13 <- predict(m13, X13dt, type = "weights")$weights
expect_equal(pr13$mpost, as.numeric(plogis(w13 + qlogis(m13$p))),
             attributes = FALSE)
cat("  [13] predict underflow -> mpost finite via plogis(weights): PASS\n")

## ========================================================
## 14. INPUT VALIDATION (FIX 4/5): negative distance -> error; N = 1 -> error.
## ========================================================
Xbad <- X12; Xbad[1, 1] <- -1
expect_error(problink_em_gamma(comparison_matrix = Xbad))
Xnan <- X12; Xnan[1, 1] <- Inf
expect_error(problink_em_gamma(comparison_matrix = Xnan))
expect_error(problink_em_gamma(comparison_matrix = X12[1, , drop = FALSE]))  # N=1
cat("  [14] input validation (negative dist, N=1): PASS\n")

cat("PASS: test_problink_em_gamma.R\n")
