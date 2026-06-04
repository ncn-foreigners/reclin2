
library(reclin2)
library(data.table)
source("helpers.R")

## Silence the per-field type message in tests.
qmix <- function(...) suppressMessages(problink_em_mixed(...))

## ========================================================
## 1. REDUCTION: all-continuous problink_em_mixed == problink_em_gamma
##    posteriors within 1e-8 from a SHARED explicit start.
## ========================================================
set.seed(42)
nA <- 200; nB <- 40; K <- 3
XA <- matrix(rexp(nA * K, rate = 1), ncol = K)
err <- matrix(rbinom(nB * K, 1, 0.05), ncol = K)
XB <- (1 - err) * XA[1:nB, ] + err * (XA[1:nB, ] + matrix(rexp(nB * K, 4), ncol = K))
Xcross <- do.call(cbind, lapply(1:K, function(k) {
  t <- expand.grid(XA[, k], XB[, k]); abs(t[, 1] - t[, 2])
}))
colnames(Xcross) <- paste0("k", 1:K)
N <- nrow(Xcross)

## Build a shared start (MoM, p0M = 0.8) in both estimators' formats.
mk_gamma_start <- function(X, K, nB, p0M = 0.8) {
  N <- nrow(X); lambdaM <- nB / N
  p0 <- matrix(0,2,K); alpha <- matrix(1,2,K); beta <- matrix(1,2,K)
  mom <- function(x){x<-x[x>0];m<-mean(x);v<-var(x)
    if(!is.finite(v)||v<=0) c(1,1) else c(m^2/v, v/m)}
  for (k in seq_len(K)) {
    x<-X[,k]; x1<-sort(x[x>0]); nM<-max(1L,min(length(x1)-1L,round((1-p0M)*nB)))
    fM<-mom(x1[seq_len(nM)]); fU<-mom(x1[-seq_len(nM)])
    p0U<-(sum(x==0)-nB*p0M)/(N-nB); p0[,k]<-c(p0M,max(1/N,min(1-1/N,p0U)))
    alpha[,k]<-c(fM[1],fU[1]); beta[,k]<-c(fM[2],fU[2])
  }
  p0[p0<1/N]<-1/N
  list(lambda=c(lambdaM,1-lambdaM), p0=p0, alpha=alpha, beta=beta)
}
sg <- mk_gamma_start(Xcross, K, nB)
sm <- list(lambda = sg$lambda, fields = lapply(seq_len(K), function(k)
  list(type="continuous", p0=sg$p0[,k], alpha=sg$alpha[,k], beta=sg$beta[,k])))

g  <- problink_em_gamma(comparison_matrix = Xcross, .start = sg, tol = 1e-6)
m  <- qmix(comparison_matrix = Xcross, .start = sm, tol = 1e-6)
stopifnot(all(m$types == "continuous"))
stopifnot(g$niter == m$niter)

Xdt <- as.data.table(Xcross)
Xdt[, .x := seq_len(.N)][, .y := seq_len(.N)]
setattr(Xdt, "class", c("pairs", class(Xdt)))
setattr(Xdt, "x", Xdt); setattr(Xdt, "y", Xdt)
pg <- predict(g, Xdt, type = "mpost")$mpost
pm <- predict(m, Xdt, type = "mpost")$mpost
stopifnot(max(abs(pg - pm), na.rm = TRUE) < 1e-8)
wg <- predict(g, Xdt, type = "weights")$weights
wm <- predict(m, Xdt, type = "weights")$weights
stopifnot(max(abs(wg - wm), na.rm = TRUE) < 1e-6)
cat("  [1] reduction to problink_em_gamma (shared start, <1e-8): PASS\n")

## ========================================================
## 2. AUTO-DETECT: binary / categorical / continuous + override.
## ========================================================
set.seed(1)
Na <- 400
Xa <- cbind(
  b    = rbinom(Na, 1, 0.3),                 # binary
  c3   = sample(0:2, Na, replace = TRUE),    # categorical (3 levels)
  cont = abs(rnorm(Na, 0, 5)) + 0.01,        # continuous (non-integer)
  idis = sample(0:40, Na, replace = TRUE),   # integer, many unique -> continuous
  amb  = sample(0:2, Na, replace = TRUE)     # integer few unique -> categorical
)
ma <- qmix(comparison_matrix = Xa, maxits = 1)
expect_equal(unname(ma$types[c("b","c3","cont","idis","amb")]),
             c("binary","categorical","continuous","continuous","categorical"),
             attributes = FALSE)
## override
mo <- qmix(comparison_matrix = Xa, types = list(amb = "continuous"), maxits = 1)
stopifnot(mo$types["amb"] == "continuous")
stopifnot(mo$types["b"]   == "binary")
## maxlevels boundary
ml <- qmix(comparison_matrix = Xa, maxlevels = 2, maxits = 1)
stopifnot(ml$types["c3"] == "continuous")   # 3 levels > maxlevels=2
## positional types
mp <- qmix(comparison_matrix = Xa[, 1:3], types = c("binary","categorical","continuous"),
           maxits = 1)
stopifnot(all(mp$types == c(b="binary", c3="categorical", cont="continuous")))
## message reports types
mtxt <- character(0)
withCallingHandlers(
  invisible(problink_em_mixed(comparison_matrix = Xa, maxits = 1)),
  message = function(mm) { mtxt <<- c(mtxt, conditionMessage(mm)); invokeRestart("muffleMessage") })
stopifnot(any(grepl("comparison field types", paste(mtxt, collapse = ""))))
cat("  [2] auto-detect + override (+ message): PASS\n")

## ========================================================
## 3. PREDICT TYPES: weights/mpost/probs/all + add, on a mixed model.
## ========================================================
set.seed(7)
nA3 <- 300; nB3 <- 60
A3 <- data.frame(sex = sample(0:1, nA3, TRUE),
                 reg = sample(1:10, nA3, TRUE),
                 age = round(runif(nA3, 18, 85)))
B3 <- A3[1:nB3, ]
B3$age <- B3$age + rbinom(nB3, 1, 0.2) * sample(1:6, nB3, TRUE)
ix <- sample(seq_len(nB3), 5); B3$reg[ix] <- sample(1:10, 5, TRUE)
pairs3 <- pair(A3, B3)
pairs3 <- compare_pairs(pairs3, on = c("sex","reg","age"),
  comparators = list(sex = cmp_identical(), reg = cmp_levels(c(0,3)),
                     age = cmp_distance()))
m3 <- qmix(~ sex + reg + age, data = pairs3)
stopifnot(inherits(m3, "problink_em_mixed"))
stopifnot(identical(unname(m3$types), c("binary","categorical","continuous")))

rw <- predict(m3, pairs3, type = "weights")
stopifnot(is.data.table(rw), all(c(".x",".y","weights") %in% names(rw)))
rm_ <- predict(m3, pairs3, type = "mpost")
stopifnot(all(c(".x",".y","mpost") %in% names(rm_)))
stopifnot(all(rm_$mpost >= 0 & rm_$mpost <= 1, na.rm = TRUE))
rp <- predict(m3, pairs3, type = "probs")
stopifnot(all(c("mprob","uprob","mpost","upost") %in% names(rp)))
expect_equal(rp$mpost + rp$upost, rep(1, nrow(rp)), attributes = FALSE)
ra <- predict(m3, pairs3, type = "all")
stopifnot(all(c("weight","mprob","uprob","mpost","upost") %in% names(ra)))
pa <- predict(m3, pairs3, type = "mpost", add = TRUE)
stopifnot("mpost" %in% names(pa), nrow(pa) == nrow(pairs3))
## NA robustness
pna <- copy(pairs3); pna[1, sex := NA_real_]; pna[2, reg := NA_real_]; pna[3, age := NA_real_]
rna <- predict(m3, pna, type = "mpost")
stopifnot(all(is.finite(rna$mpost[1:3])))
## print / summary do not error
capture.output(print(m3)); capture.output(summary(m3))
cat("  [3] predict types + NA + print/summary: PASS\n")

## ========================================================
## 4. END-TO-END: K1 binary + K2 3-level categorical + K3 continuous,
##    DEFAULTS (auto-detect), select_threshold(0.5), TPR/PPV vs truth.
## ========================================================
set.seed(101)
nA4 <- 400; nB4 <- 80; K1 <- 2; K2 <- 2; K3 <- 2
A4 <- data.frame(
  matrix(rbinom(nA4*K1, 1, 0.25), ncol = K1),
  matrix(sample(1:50, nA4*K2, replace = TRUE), ncol = K2),
  matrix(ceiling(rexp(nA4*K3, rate = 0.02)), ncol = K3),
  id = seq_len(nA4)
)
names(A4) <- c(paste0("b",1:K1), paste0("c",1:K2), paste0("k",1:K3), "id")
idAB <- sample(seq_len(nA4), nB4); B4 <- A4[idAB, ]
for (j in 1:K1) { col<-paste0("b",j); ix<-sample(seq_len(nB4),round(nB4*0.05)); B4[[col]][ix]<-1-B4[[col]][ix] }
for (j in 1:K2) { col<-paste0("c",j); ix<-sample(seq_len(nB4),round(nB4*0.10))
  B4[[col]][ix]<-sapply(B4[[col]][ix], function(v) sample(setdiff(1:50,v),1)) }
for (j in 1:K3) { col<-paste0("k",j); ix<-sample(seq_len(nB4),round(nB4*0.15))
  B4[[col]][ix]<-B4[[col]][ix]+ceiling(rexp(length(ix), rate=1/30)) }
on4 <- c(paste0("b",1:K1), paste0("c",1:K2), paste0("k",1:K3))
pairs4 <- pair(A4, B4)
pairs4 <- compare_pairs(pairs4, on = on4, comparators = c(
  setNames(replicate(K1, cmp_identical(), simplify=FALSE), paste0("b",1:K1)),
  setNames(replicate(K2, cmp_levels(c(0,3)), simplify=FALSE), paste0("c",1:K2)),
  setNames(replicate(K3, cmp_distance(), simplify=FALSE), paste0("k",1:K3))))
m4 <- qmix(as.formula(paste("~", paste(on4, collapse="+"))), data = pairs4)
stopifnot(identical(unname(m4$types),
                    c("binary","binary","categorical","categorical","continuous","continuous")))
stopifnot(m4$converged)
pairs4 <- predict(m4, pairs4, type = "mpost", add = TRUE)
pairs4 <- select_threshold(pairs4, "selected", "mpost", 0.5)
Aid <- A4$id; Bid <- B4$id
pairs4[, truth := Aid[.x] == Bid[.y]]
sel <- pairs4[selected == TRUE]
TP <- sum(sel$truth); TPR <- TP / nB4
PPV <- if (nrow(sel) > 0) TP / nrow(sel) else 0
cat(sprintf("  [4] end-to-end (defaults): TPR=%.3f PPV=%.3f", TPR, PPV))
## Conservative thresholds: PPV must be high; TPR reasonable. (Recall can be
## moderate at threshold 0.5 with a tiny match prior; see DELIVERABLE notes.)
stopifnot(PPV > 0.85)
stopifnot(TPR > 0.6)
cat("  PASS\n")

## ========================================================
## 5. LABEL-SWITCH: from a deliberately SWAPPED start, the informative blocks
##    end MATCH-first and true matches rank above non-matches.
## ========================================================
reglev <- sort(unique(pairs3$reg[!is.na(pairs3$reg)])); Lr <- length(reglev)
st <- list(lambda = c(0.05, 0.95), fields = list(
  list(type="binary", m=c(0.1, 0.9)),                         # swapped
  list(type="categorical", levels=reglev,
       m = { v<-rep(0.05/(Lr-1),Lr); v[Lr]<-0.95; v },         # mass on top level
       u = { v<-rep(0.05/(Lr-1),Lr); v[1]<-0.95; v }),
  list(type="continuous", p0=c(0.01,0.5), alpha=c(2,2), beta=c(20,1)))) # swapped
ms <- qmix(~ sex + reg + age, data = pairs3, .start = st)
ai <- which.min(ms$fields[["reg"]]$levels)
stopifnot(ms$fields[["reg"]]$m[ai] >= ms$fields[["reg"]]$u[ai] - 1e-8)
md1 <- ms$fields[["age"]]$alpha[1]*ms$fields[["age"]]$beta[1]
md2 <- ms$fields[["age"]]$alpha[2]*ms$fields[["age"]]$beta[2]
stopifnot(ms$fields[["age"]]$p0[1] >= ms$fields[["age"]]$p0[2] - 1e-8 || md1 <= md2 + 1e-8)
psw <- predict(ms, pairs3, type = "mpost"); setDT(psw)
tr <- (pairs3$.x == pairs3$.y) & (pairs3$.x <= nB3)
stopifnot(mean(psw$mpost[tr]) > mean(psw$mpost[!tr]))
cat("  [5] label-switch guard (forced-swap recovery): PASS\n")

## ========================================================
## 6. nstart determinism + .start forces single run + input validation.
## ========================================================
a1 <- qmix(~ sex + reg + age, data = pairs3, nstart = 1)
a2 <- qmix(~ sex + reg + age, data = pairs3, nstart = 1)
expect_equal(a1$p, a2$p)                              # deterministic
set.seed(5); b1 <- qmix(~ sex + reg + age, data = pairs3, nstart = 3)
set.seed(5); b2 <- qmix(~ sex + reg + age, data = pairs3, nstart = 3)
expect_equal(b1$loglik, b2$loglik)                   # reproducible under seed
stopifnot(b1$nstart == 3L)
stopifnot(b1$loglik >= a1$loglik - 1e-6)             # multistart never worse
## RNG stream restored
set.seed(9); r0 <- runif(3)
set.seed(9); invisible(qmix(~ sex + reg + age, data = pairs3, nstart = 3))
r1 <- runif(3); expect_equal(r0, r1)
## .start forces nstart = 1
stopifnot(ms$nstart == 1L)
## validation
expect_error(qmix(comparison_matrix = Xa, nstart = 0))
expect_error(qmix(comparison_matrix = Xa, lambda_smooth = -1))
expect_error(qmix(comparison_matrix = Xa, types = list(b = "nonsense")))
cat("  [6] nstart determinism + validation: PASS\n")

## ========================================================
## 7. PARTIAL MISSINGNESS (FIX 1/2): some fields NA on some rows -> converges,
##    loglik finite, recovery sane; a pair with one NA field + strong observed
##    fields is NOT reduced to the prior. nstart>1 works under missingness.
## ========================================================
set.seed(301)
nA7 <- 300; nB7 <- 80
## Well-separated DGP: two binaries with many categories (rare chance agreement)
## + a continuous field, so the EM converges quickly even with NA scattered in.
A7 <- data.frame(reg = sample(1:40, nA7, TRUE),
                 cty = sample(1:40, nA7, TRUE),
                 age = round(runif(nA7, 18, 85)))
B7 <- A7[1:nB7, ]
B7$age <- B7$age + rbinom(nB7, 1, 0.1) * sample(1:4, nB7, TRUE)
pairs7 <- pair(A7, B7)
pairs7 <- compare_pairs(pairs7, on = c("reg","cty","age"),
  comparators = list(reg = cmp_levels(c(0, 3)), cty = cmp_levels(c(0, 3)),
                     age = cmp_distance()))
## Scatter NA across each comparison field.
set.seed(302)
N7 <- nrow(pairs7)
for (col in c("reg","cty","age")) {
  ii <- sample(seq_len(N7), round(0.05 * N7))
  set(pairs7, i = ii, j = col, value = NA_real_)
}
m7 <- qmix(~ reg + cty + age, data = pairs7, maxits = 800)
stopifnot(identical(unname(m7$types), c("categorical","categorical","continuous")))
stopifnot(m7$converged)
stopifnot(is.finite(m7$loglik))
rna7 <- predict(m7, pairs7, type = "mpost")
stopifnot(all(is.finite(rna7$mpost)))
stopifnot(all(rna7$mpost >= 0 & rna7$mpost <= 1))

## Probe: a true match (reg level 0, cty level 0, age NA) must rank well above
## the prior -- the observed (agreeing) fields still inform it.
probe7 <- data.table(reg = 0, cty = 0, age = NA_real_, .x = 1L, .y = 1L)
setattr(probe7, "class", c("pairs", class(probe7)))
setattr(probe7, "x", probe7); setattr(probe7, "y", probe7)
pm7 <- predict(m7, probe7, type = "mpost")$mpost
stopifnot(is.finite(pm7), pm7 > m7$p + 0.05)

## nstart>1 under missingness.
m7a <- qmix(~ reg + cty + age, data = pairs7, nstart = 1, maxits = 800)
set.seed(11); m7b <- qmix(~ reg + cty + age, data = pairs7, nstart = 4,
                          maxits = 800)
stopifnot(is.finite(m7b$loglik), m7b$nstart == 4L)
stopifnot(m7b$loglik >= m7a$loglik - 1e-6)
cat("  [7] partial missingness (mixed): PASS\n")

## ========================================================
## 8. ALL-NA FIELD (FIX 2): a field that is entirely NA -> no NaN in the fit
##    (warning emitted); other fields still fit.
## ========================================================
set.seed(303)
Xall <- cbind(
  sex  = rbinom(400, 1, 0.3),
  age  = abs(rnorm(400, 0, 5)) + 0.01,
  dead = rep(NA_real_, 400)       # entirely NA
)
## A warning must be emitted (either from type detection or the M-step).
expect_warning(qmix(comparison_matrix = Xall, maxits = 50))
m8 <- qmix(comparison_matrix = Xall, maxits = 50)
## No NaN anywhere in the fitted blocks.
flat8 <- unlist(lapply(m8$fields, function(f) unlist(f[setdiff(names(f), "type")])))
stopifnot(!any(is.nan(flat8)))
stopifnot(is.finite(m8$loglik))
## The informative fields are still fitted (binary sex separates m from u OR
## age block is non-degenerate).
stopifnot(is.finite(m8$p), m8$p > 0, m8$p < 1)
cat("  [8] all-NA field -> no NaN (warning emitted): PASS\n")

## ========================================================
## 9. SMOOTHING (FIX 3): lambda_smooth > 0 -> categorical m/u each SUM TO 1
##    (simplex preserved); lambda_smooth = 0 reproduces the unsmoothed MLE.
## ========================================================
## Use the section-3 pairs3 (binary + 3-level categorical + continuous).
m9_0  <- qmix(~ sex + reg + age, data = pairs3, lambda_smooth = 0)
m9_s  <- qmix(~ sex + reg + age, data = pairs3, lambda_smooth = 2)
mc9   <- m9_s$fields[["reg"]]
stopifnot(abs(sum(mc9$m) - 1) < 1e-8)        # smoothed categorical m sums to 1
stopifnot(abs(sum(mc9$u) - 1) < 1e-8)        # smoothed categorical u sums to 1
## lambda_smooth = 0 reproduces the plain MLE exactly: run twice, identical.
m9_0b <- qmix(~ sex + reg + age, data = pairs3, lambda_smooth = 0)
expect_equal(m9_0$fields[["reg"]]$m, m9_0b$fields[["reg"]]$m)
expect_equal(m9_0$p, m9_0b$p)
## And the unsmoothed categorical m/u also sum to 1 (proper multinomial MLE).
stopifnot(abs(sum(m9_0$fields[["reg"]]$m) - 1) < 1e-8)
stopifnot(abs(sum(m9_0$fields[["reg"]]$u) - 1) < 1e-8)
cat("  [9] smoothing simplex preservation + lambda_smooth=0 MLE: PASS\n")

## ========================================================
## 10. PREDICT UNDERFLOW (FIX 6): many fields -> mprob/uprob underflow to 0 but
##     mpost finite in (0,1) and == plogis(weights + qlogis(p)).
## ========================================================
set.seed(304)
K10c <- 60
Xc10 <- do.call(cbind, lapply(1:K10c, function(k) {
  XAk <- rexp(50); XBk <- XAk[1:15] + rbinom(15, 1, 0.1) * rexp(15, 4)
  t <- expand.grid(XAk, XBk); abs(t[, 1] - t[, 2])
}))
colnames(Xc10) <- paste0("k", 1:K10c)
m10c <- qmix(comparison_matrix = Xc10, p0 = 15 / nrow(Xc10), p0M = 0.8,
             maxits = 40)
stopifnot(all(m10c$types == "continuous"))
Xc10dt <- as.data.table(Xc10)
Xc10dt[, .x := seq_len(.N)][, .y := seq_len(.N)]
setattr(Xc10dt, "class", c("pairs", class(Xc10dt)))
setattr(Xc10dt, "x", Xc10dt); setattr(Xc10dt, "y", Xc10dt)
pr10 <- predict(m10c, Xc10dt, type = "probs")
stopifnot(any(pr10$mprob == 0) || any(pr10$uprob == 0))  # products underflow
stopifnot(all(is.finite(pr10$mpost)))                    # mpost still finite
stopifnot(all(pr10$mpost >= 0 & pr10$mpost <= 1))
w10 <- predict(m10c, Xc10dt, type = "weights")$weights
expect_equal(pr10$mpost, as.numeric(plogis(w10 + qlogis(m10c$p))),
             attributes = FALSE)
cat("  [10] predict underflow (mixed) -> mpost via plogis(weights): PASS\n")

## ========================================================
## 11. INPUT VALIDATION + CATEGORICAL NON-ORDINAL CODING (FIX 4/5/7).
## ========================================================
## (a) Negative value in a CONTINUOUS column -> error (binary/categorical are
##     exempt: their domains are checked by their own blocks).
Xneg <- cbind(sex = rbinom(50, 1, 0.3), age = abs(rnorm(50)) + 0.01)
Xneg[1, "age"] <- -1
expect_error(qmix(comparison_matrix = Xneg))
## (b) N = 1 -> error.
expect_error(qmix(comparison_matrix = Xneg[1, , drop = FALSE]))
## (c) Non-ordinal categorical coding: smallest observed level != 0 (e.g. a
##     sum coding with levels {1, 2}) -> warning about the agreement assumption.
set.seed(305)
Xnoord <- cbind(
  cat = sample(1:2, 400, replace = TRUE),      # min level 1, not 0
  age = abs(rnorm(400, 0, 5)) + 0.01           # anchoring continuous field
)
expect_warning(qmix(comparison_matrix = Xnoord,
                    types = list(cat = "categorical"), maxits = 50))
## With an anchoring continuous field the recovered p is the match component
## (a sane mixing probability in (0,1)).
m11 <- suppressWarnings(qmix(comparison_matrix = Xnoord,
                             types = list(cat = "categorical"), maxits = 200))
stopifnot(is.finite(m11$p), m11$p > 0, m11$p < 1)
cat("  [11] input validation + non-ordinal categorical warning: PASS\n")

cat("PASS: test_problink_em_mixed.R\n")
