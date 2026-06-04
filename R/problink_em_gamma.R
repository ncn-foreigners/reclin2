
#' Continuous Fellegi-Sunter EM estimator using hurdle-gamma mixture
#'
#' @description
#' Estimates the parameters of a two-component hurdle-gamma mixture model for
#' continuous comparison distances (e.g., absolute differences produced by
#' \code{\link{cmp_distance}}).  Component 1 is the \strong{match} component
#' and component 2 is the \strong{non-match} component.  The model is:
#'
#' \deqn{f(g \mid C) = p_{0,C} \cdot \mathbf{1}[g=0] +
#'   (1 - p_{0,C}) \cdot \text{Gamma}(g; \alpha_C, \beta_C) \cdot \mathbf{1}[g>0]}
#'
#' The mixing proportion \eqn{p_M} is the overall probability that a pair is a
#' match.
#'
#' @param formula a formula of the form \code{~ var1 + var2 + ...} specifying
#'   the comparison distance columns in \code{data}.
#' @param data a \code{data.table} (or \code{data.frame}) with comparison
#'   distances, typically produced by \code{\link{compare_pairs}} with
#'   \code{\link{cmp_distance}} comparators.  Either \code{data} or
#'   \code{comparison_matrix} must be supplied.
#' @param comparison_matrix a numeric matrix whose columns are comparison
#'   distances for each field.  Use this as an alternative to \code{data}/
#'   \code{formula} when working with a plain matrix (e.g., for the
#'   reference-equivalence verification).
#' @param p0 optional initial estimate of the probability that a pair is a
#'   match (scalar in \eqn{(0,1)}).  When \code{NULL} (the default) a
#'   data-derived value is used: if \code{data} is a \code{pairs} object that
#'   carries its source tables as attributes, \code{p0} is set to
#'   \code{min(nA, nB) / nrow(pairs)} (the largest possible match prevalence,
#'   i.e. every record in the smaller table matches); for the
#'   \code{comparison_matrix} interface (or when the attributes are missing) it
#'   falls back to \code{0.05}.  The resolved value is clamped to
#'   \code{[1/N, 0.5]}.  \code{p0} is only the EM \emph{start} value -- the match
#'   prevalence is re-estimated by the algorithm -- but a sensible start
#'   improves the robustness of the mixture EM.  Supply an explicit value to
#'   override the default.
#' @param p0M initial estimate of \eqn{p_0} for the match component, i.e., the
#'   probability of an \emph{exact} match (distance zero) for a true match pair.
#'   When \code{NULL} (the default) the value is \emph{derived from the data per
#'   field}: for field \eqn{k} it is set to
#'   \code{clamp(z_k / lambdaM_init, 1/N, 1 - 1/N)}, where
#'   \code{z_k = mean(X[, k] == 0)} is the observed exact-zero fraction and
#'   \code{lambdaM_init} is the resolved match prevalence (\code{p0}).  This
#'   inflates the raw zero fraction by the (small) match prevalence because the
#'   zeros that do occur are predominantly contributed by true matches, giving a
#'   well-calibrated hurdle start.  In the \dQuote{pure continuous} regime (no
#'   exact matches at all) \eqn{z_k = 0} so \eqn{p_{0M} \approx 1/N} and the
#'   hurdle point-mass is effectively inactive (the model degrades gracefully to
#'   a plain two-component gamma mixture).  Supply an explicit numeric value in
#'   \eqn{(0, 1)} to override (e.g. \code{p0M = 0.8} reproduces the historical
#'   fixed start).
#' @param tol convergence tolerance.  The algorithm stops when the relative
#'   change in \emph{all} parameters is below \code{tol}.  Default \code{1e-6}.
#' @param maxits maximum number of EM iterations.  Default \code{500}.
#' @param use_envstats logical; if \code{TRUE} and \pkg{EnvStats} is available
#'   (via \code{requireNamespace}), use \code{EnvStats::egamma} for the initial
#'   gamma parameter estimates.  Otherwise (default \code{FALSE}) use the
#'   method-of-moments estimator.
#' @param nstart integer >= 1; number of EM starts (default \code{1}).  With
#'   \code{nstart = 1} the estimator is fully \strong{deterministic} (the single
#'   data-driven start described above) and behaves exactly as in previous
#'   versions.  With \code{nstart > 1} the data-driven start is used first and
#'   \code{nstart - 1} additional \emph{perturbed} starts are run; the fit with
#'   the highest observed log-likelihood is returned.  The perturbations use the
#'   package RNG and \emph{respect any seed set by the user} (the current RNG
#'   state is saved and restored on exit), so results are reproducible under
#'   \code{set.seed()} and the caller's random stream is left untouched.
#' @param lambda_smooth non-negative smoothing strength (default \code{0} = exact
#'   maximum-likelihood updates).  When \code{> 0}, the hurdle point-masses
#'   \eqn{p_0} are shrunk towards their (per-component, per-field) \emph{starting}
#'   values via a Dirichlet/Beta pseudo-count update
#'   \eqn{(\sum q\,\mathbf{1}[g=0] + \lambda p_0^{0}) / (\sum q + \lambda)}.  Only
#'   the hurdle probabilities are smoothed; the gamma shape and scale (which are
#'   not probabilities) are left untouched.  The default \code{lambda_smooth = 0}
#'   reproduces the unsmoothed fit exactly.  Mirrors the \code{lambda_smooth}
#'   argument of \code{\link{problink_em_mixed}} and is useful when a field has
#'   very few exact-zero (exact-match) distances.
#' @param .start optional internal/testing hook.  A list with elements
#'   \code{lambda} (length-2), \code{p0}, \code{alpha}, \code{beta} (each a
#'   \eqn{2 \times K} matrix) giving \emph{explicit} EM starting values that
#'   bypass the initialisation heuristic entirely.  Used by the
#'   reference-equivalence verification to run this estimator and the reference
#'   implementation from byte-identical starts so that the \emph{EM update
#'   equations} can be checked in isolation from the init heuristic.  Forces
#'   \code{nstart = 1}.  Not intended for general use.
#'
#' @details
#' The EM algorithm alternates between:
#' \describe{
#'   \item{E-step}{Compute pair-level posterior match probability
#'     \eqn{z_i = \Pr(\text{match} \mid g_i)}.}
#'   \item{M-step}{Update \eqn{p_{0,C}}, shape (\eqn{\alpha_C}), and scale
#'     (\eqn{\beta_C}) for each field and component via closed-form updates;
#'     shape is solved via \code{uniroot} with an \code{nlminb} fallback;
#'     the mixing proportion \eqn{p_M} is updated as the mean posterior.}
#' }
#' After convergence, a regime-aware label-switching guard ensures component 1
#' is the match component (it prefers the component with the larger hurdle mass
#' \eqn{p_0}; on a near-tie -- e.g. the pure-continuous regime where both
#' \eqn{p_0 \approx 1/N} -- it falls back to the component with the smaller mean
#' distance \eqn{\alpha\beta}).
#'
#' \strong{No exact matches.}  When the comparison distances contain no exact
#' zeros, the hurdle point-mass is not supported by the data and the maximum-
#' likelihood value of every \eqn{p_0} is \eqn{1/N} (its floor); the model then
#' reduces to an ordinary two-component gamma mixture.  This is the correct
#' behaviour, not a failure: \code{problink_em_gamma} detects it, reports a
#' \dQuote{degenerated to gamma mixture} flag from \code{print}, and recovers
#' matches from the gamma component separation alone.
#'
#' \strong{Numerical hardening.}  All M-step denominators are floored away from
#' zero; mixing weights are clamped to \code{[1/N, 1 - 1/N]}; a component whose
#' effective sample size collapses keeps its previous shape/scale; the shape
#' update degrades through \code{uniroot} -> \code{nlminb} -> weighted method-of-
#' moments -> previous value and can never error.  These guards are inactive in
#' the well-separated regime and therefore do \emph{not} change the converged
#' values when exact matches are present.
#'
#' \strong{Missing values.}  \code{NA} comparison distances are handled
#' \emph{per field} on the log scale: an \code{NA} in field \eqn{k} contributes a
#' neutral log-density of \code{0} (a factor of \code{1}) to that pair's E-step,
#' so the pair is still informed by its \emph{observed} fields rather than being
#' reduced to the prior.  The per-field hurdle (\eqn{p_0}) M-step uses only rows
#' observed on that field for both its numerator and denominator.  The
#' observed-data log-likelihood is accumulated on the log scale (log-sum-exp), so
#' it stays finite under missingness or underflow and \code{nstart} can select
#' the best start.
#'
#' @return
#' An object of class \code{problink_em_gamma} (a list) with components:
#' \describe{
#'   \item{\code{p0}}{2 x K matrix of hurdle probabilities (row 1 = match, row 2 = non-match).}
#'   \item{\code{shape}}{2 x K matrix of gamma shape parameters.}
#'   \item{\code{scale}}{2 x K matrix of gamma scale parameters.}
#'   \item{\code{p}}{scalar mixing probability (probability pair is a match).}
#'   \item{\code{p0_init}}{the initial match prevalence used as the EM start
#'     value (the resolved/clamped \code{p0}).}
#'   \item{\code{p0M_init}}{the resolved per-field initial hurdle (exact-match)
#'     probabilities (length-K, or a scalar when an explicit \code{p0M} was
#'     supplied).}
#'   \item{\code{loglik}}{the observed-data log-likelihood at the returned fit.}
#'   \item{\code{hurdle_active}}{logical flag; \code{TRUE} when the fitted hurdle
#'     point-mass is meaningfully above \eqn{1/N} (exact matches inform the fit),
#'     \code{FALSE} when the model has degenerated to a plain gamma mixture.}
#'   \item{\code{by}}{character vector of field names.}
#'   \item{\code{niter}}{number of EM iterations performed (for the returned start).}
#'   \item{\code{nstart}}{number of EM starts run.}
#'   \item{\code{converged}}{logical; \code{TRUE} if the algorithm converged.}
#'   \item{\code{use_envstats}}{the value of \code{use_envstats} used.}
#'   \item{\code{call}}{the matched call.}
#' }
#'
#' @references
#' Fellegi, I. and A. Sunter (1969). "A Theory for Record Linkage",
#' \emph{Journal of the American Statistical Association}. 64 (328): 1183-1210.
#' \doi{10.2307/2286061}.
#'
#' @examples
#' set.seed(1)
#' nA <- 200; nB <- 40; K <- 2
#' XA <- matrix(rexp(nA * K, rate = 1), ncol = K)
#' XB <- XA[1:nB, ] + matrix(rbinom(nB * K, 1, 0.2) *
#'   rexp(nB * K, rate = 2), ncol = K)
#' dfA <- as.data.frame(XA); names(dfA) <- paste0("v", 1:K)
#' dfB <- as.data.frame(XB); names(dfB) <- paste0("v", 1:K)
#' library(reclin2)
#' pairs <- pair(dfA, dfB)
#' pairs <- compare_pairs(pairs, on = paste0("v", 1:K),
#'   default_comparator = cmp_distance())
#' model <- problink_em_gamma(~ v1 + v2, data = pairs)
#' print(model)
#' res <- predict(model, pairs, type = "mpost", add = FALSE)
#'
#' \dontshow{gc()}
#'
#' @importFrom stats dgamma uniroot nlminb var rnorm plogis
#' @export
problink_em_gamma <- function(formula, data, comparison_matrix,
    p0 = NULL, p0M = NULL, tol = 1e-6, maxits = 500,
    use_envstats = FALSE, nstart = 1L, lambda_smooth = 0, .start = NULL) {

  cl <- match.call()

  # ----- Resolve comparison matrix and field names -----
  from_data <- FALSE
  if (!missing(comparison_matrix) && !is.null(comparison_matrix)) {
    X  <- as.matrix(comparison_matrix)
    by <- if (!is.null(colnames(X))) colnames(X) else paste0("V", seq_len(ncol(X)))
    colnames(X) <- by
  } else {
    if (missing(data) || is.null(data))
      stop("Either 'data' or 'comparison_matrix' must be supplied.")
    if (missing(formula) || is.null(formula))
      stop("'formula' must be supplied when 'data' is used.")
    by <- if (length(formula) == 3) all.vars(formula[[3]]) else all.vars(formula)
    if (!all(by %in% names(data)))
      stop("Not all variables in formula are present in data.")
    X <- as.matrix(as.data.frame(data)[, by, drop = FALSE])
    from_data <- TRUE
  }

  K <- ncol(X)
  N <- nrow(X)

  # A 2-component mixture is not identifiable from a single observation.
  if (N < 2L)
    stop("need at least 2 comparison pairs to fit a 2-component mixture.")

  # ----- Validate continuous distances -----
  # The hurdle-gamma model is defined on non-negative, finite distances. NA is
  # allowed (handled per field), but a negative or non-finite (Inf/NaN) value is
  # a sign the input is not a valid distance (e.g. a raw difference, not |x-y|).
  if (any(X < 0, na.rm = TRUE) || any(is.nan(X)) ||
      any(is.infinite(X)))
    stop("comparison distances must be non-negative and finite ",
         "(NA allowed); did you pass raw differences instead of |x - y|? ",
         "Use cmp_distance(), which returns |x - y| / scale.")

  nstart <- as.integer(nstart)
  if (length(nstart) != 1L || is.na(nstart) || nstart < 1L)
    stop("'nstart' must be a single integer >= 1.")
  if (!is.null(.start)) nstart <- 1L  # explicit start => single deterministic run
  if (!is.numeric(lambda_smooth) || length(lambda_smooth) != 1L ||
      !is.finite(lambda_smooth) || lambda_smooth < 0)
    stop("'lambda_smooth' must be a single non-negative number.")

  # ----- Resolve initial match prevalence p0 -----
  # p0 is only the EM *start* value for the match-mixing weight; the match
  # prevalence is re-estimated during the EM. A data-derived default makes the
  # mixture EM more robust than a fixed guess.
  if (!is.null(p0)) {
    # User supplied an explicit value: validate and use as-is.
    if (length(p0) != 1L || !is.finite(p0) || p0 <= 0 || p0 >= 1)
      stop("'p0' must be a single number in (0, 1).")
  } else {
    # Derive a smart default. When fitting from a 'pairs' object that carries
    # its source tables as attributes, the best case is that every record in
    # the smaller table matches, i.e. min(nA, nB) matches among N pairs.
    p0 <- 0.05  # fallback (e.g. comparison_matrix interface)
    if (from_data) {
      ax <- attr(data, "x")
      ay <- attr(data, "y")
      if (!is.null(ax) && !is.null(ay)) {
        nA <- nrow(ax)
        nB <- nrow(ay)
        if (length(nA) && length(nB) && nA > 0 && nB > 0)
          p0 <- min(nA, nB) / N
      }
    }
  }
  # Clamp to a sensible range for an EM start value.
  p0 <- min(max(p0, 1/N), 0.5)
  p0_init <- p0                 # resolved match prevalence (EM start)
  lambdaM_init <- p0_init       # alias: match mixing-weight start

  # ----- Resolve initial hurdle (exact-match) probability p0M (per field) -----
  # p0M[k] is the probability that a TRUE match has distance exactly zero in
  # field k. We derive it from the data: the observed exact-zero fraction z_k is
  # contributed mostly by true matches, so the *match-component* zero mass is
  # z_k inflated by 1/lambdaM_init, clamped to (0,1). In the pure-continuous
  # regime z_k = 0, so p0M[k] ~ 1/N and the hurdle is inactive (the model is a
  # plain 2-component gamma mixture). An explicit numeric p0M overrides this and
  # is applied to all fields (e.g. p0M = 0.8 reproduces the old fixed start).
  z_k <- colMeans(X == 0, na.rm = TRUE)
  z_k[!is.finite(z_k)] <- 0
  if (!is.null(p0M)) {
    if (length(p0M) != 1L || !is.finite(p0M) || p0M <= 0 || p0M >= 1)
      stop("'p0M' must be a single number in (0, 1).")
    p0M_vec <- rep(p0M, K)
  } else {
    p0M_vec <- clamp01_(z_k / lambdaM_init, lo = 1/N)
  }
  p0M_vec  <- clamp01_(p0M_vec, lo = 1/N)
  p0M_init <- if (!is.null(p0M)) p0M else p0M_vec
  names(p0M_vec) <- by

  # Minimum effective sample size before updating shape/scale
  MIN_EFF_N <- 3

  # ----- Resolve the deterministic start -----
  if (!is.null(.start)) {
    base_start <- .check_start(.start, K = K, N = N)
  } else {
    base_start <- .init_hurdle_gammaK(X, K = K, p0_init = p0,
                                      p0M_init = p0M_vec, N = N,
                                      use_envstats = use_envstats)
  }

  # ----- Run EM from one or more starts; keep best log-likelihood -----
  # The EM update equations live in .run_em_gammaK() and are NOT perturbed by
  # nstart: only the *starting values* differ. nstart = 1 (default) runs exactly
  # the single deterministic start, so behaviour is unchanged.
  run_one <- function(start)
    .run_em_gammaK(X, K = K, N = N, start = start, tol = tol,
                   maxits = maxits, MIN_EFF_N = MIN_EFF_N,
                   lambda_smooth = lambda_smooth)

  best <- run_one(base_start)

  if (nstart > 1L) {
    # Respect (and restore) the user's RNG state so perturbed starts are
    # reproducible under set.seed() without disturbing the caller's stream.
    # We capture the *current* seed (so perturbations draw deterministically
    # from the user's stream after set.seed()) and restore it on exit via
    # set_seed(<saved seed vector>), leaving the caller's stream untouched.
    if (!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
      set_seed(NULL)
    old_seed <- get(".Random.seed", envir = .GlobalEnv)
    on.exit(set_seed(old_seed), add = TRUE)
    for (s in seq_len(nstart - 1L)) {
      pert <- .perturb_start(base_start, X = X, K = K, N = N,
                             use_envstats = use_envstats)
      fit  <- tryCatch(run_one(pert), error = function(e) NULL)
      if (!is.null(fit) && is.finite(fit$loglik) &&
          fit$loglik > best$loglik)
        best <- fit
    }
  }

  p0_mle     <- best$p0
  shape_mle  <- best$shape
  scale_mle  <- best$scale
  lambda_mle <- best$lambda
  iter       <- best$iter
  converge   <- best$converged
  loglik     <- best$loglik

  # ----- Label-switching guard: ensure component 1 = MATCH -----
  # Regime-aware: the match component is the one with the larger hurdle mass p0
  # (true matches contribute the exact zeros). When p0 is ~tied -- in particular
  # the pure-continuous regime where both p0 ~ 1/N and the hurdle is inactive --
  # we fall back to the SMALLER mean distance (shape * scale), since the match
  # component's records are more similar. All quantities are finite (the EM
  # guards guarantee it), so we never index NA.
  mean_p0_1   <- mean(p0_mle[1, ])
  mean_p0_2   <- mean(p0_mle[2, ])
  mean_dist_1 <- mean(shape_mle[1, ] * scale_mle[1, ])
  mean_dist_2 <- mean(shape_mle[2, ] * scale_mle[2, ])

  p0_scale <- max(abs(c(mean_p0_1, mean_p0_2)), tol)
  p0_tied  <- abs(mean_p0_1 - mean_p0_2) <= tol * p0_scale
  # Swap if component 2 looks more like the match component: larger hurdle mass,
  # or -- on a near-tie in p0 -- smaller mean distance.
  needs_swap <- (!p0_tied && mean_p0_2 > mean_p0_1) ||
                (p0_tied && mean_dist_2 < mean_dist_1)
  if (!isTRUE(needs_swap)) needs_swap <- FALSE

  if (needs_swap) {
    p0_mle     <- p0_mle[c(2, 1), , drop = FALSE]
    shape_mle  <- shape_mle[c(2, 1), , drop = FALSE]
    scale_mle  <- scale_mle[c(2, 1), , drop = FALSE]
    lambda_mle <- rev(lambda_mle)
  }

  rownames(p0_mle) <- rownames(shape_mle) <- rownames(scale_mle) <-
    c("match", "non-match")
  colnames(p0_mle) <- colnames(shape_mle) <- colnames(scale_mle) <- by

  # Diagnostic: is the hurdle point-mass active, or has the model degenerated to
  # a plain gamma mixture? "Active" means the match-component p0 sits meaningfully
  # above its 1/N floor (i.e. exact matches genuinely inform the fit).
  hurdle_active <- any(p0_mle[1, ] > (1/N) * (1 + 1e-6) + 1e-12)

  structure(
    list(
      p0          = p0_mle,
      shape       = shape_mle,
      scale       = scale_mle,
      p           = lambda_mle[1],   # pM
      p0_init     = p0_init,         # initial match prevalence (EM start value)
      p0M_init    = p0M_init,        # initial / derived exact-match probability
      loglik      = loglik,          # observed-data log-likelihood at the fit
      hurdle_active = hurdle_active, # FALSE => degenerated to gamma mixture
      by          = by,
      niter       = iter,
      nstart      = nstart,
      converged   = converge,
      use_envstats = use_envstats,
      call        = cl
    ),
    class = "problink_em_gamma"
  )
}


# =========================================================================
# Internal helpers (file-local, not exported)
# =========================================================================

#' Hurdle-gamma density (vectorised, NA-safe)
#' @noRd
dhgamma_ <- function(x, p0, shape, scale) {
  y <- rep(NA_real_, length(x))
  ind0 <- which(x == 0)
  ind1 <- which(x > 0)
  if (length(ind0) > 0) y[ind0] <- p0
  if (length(ind1) > 0) y[ind1] <- (1 - p0) * dgamma(x[ind1], shape = shape, scale = scale)
  y
}

#' Log hurdle-gamma density for one component (vectorised, per-field NA-safe).
#'
#' Returns the natural log of \code{dhgamma_(g, ...)}. A genuine \code{NA}
#' comparison value gets log-density 0 (a NEUTRAL factor of 1, so an observed
#' field is never discarded just because some other field is NA); a finite but
#' underflowed/zero density (e.g. p0 = 0, or a gamma density that underflows) is
#' floored at \code{log(1e-300)} so the E-step can never see \code{-Inf}. This
#' mirrors \code{.logdens_field()} in the mixed estimator so the all-continuous
#' case reduces exactly.
#' @noRd
.log_dhgamma <- function(g, p0, shape, scale) {
  ld <- log(pmax(dhgamma_(g, p0, shape, scale), 1e-300))
  ld[!is.finite(ld)] <- log(1e-300)   # guard residual -Inf (e.g. log(0))
  ld[is.na(g)] <- 0                    # NA comparison -> neutral log-density 0
  ld
}

#' Method-of-moments gamma estimator for a positive vector
#' @noRd
init_gamma_mom <- function(x) {
  x <- x[!is.na(x) & x > 0]
  if (length(x) < 2) return(c(shape = 1, scale = 1))
  m <- mean(x)
  v <- var(x)
  if (v <= 0) return(c(shape = 1, scale = m))
  shape <- m^2 / v
  scale <- v / m
  c(shape = shape, scale = scale)
}

#' Fit gamma start values, optionally using EnvStats::egamma
#' @noRd
fit_gamma_start <- function(x, use_envstats = FALSE) {
  x <- x[!is.na(x) & x > 0]
  if (length(x) < 2) return(c(shape = 1, scale = 1))
  if (use_envstats && requireNamespace("EnvStats", quietly = TRUE)) {
    tryCatch({
      fit <- EnvStats::egamma(x)
      return(fit$parameters)   # named c(shape=, scale=)
    }, error = function(e) NULL)
  }
  init_gamma_mom(x)
}

#' Clamp values to [lo, 1-lo]
#'
#' \code{lo} is itself capped at \code{0.5 - 1e-8} so that the lower and upper
#' bounds never cross (which would happen for a tiny N where \code{lo = 1/N} can
#' exceed \code{1 - lo}); this keeps \code{lo < 1 - lo} for any input.
#' @noRd
clamp01_ <- function(p, lo = 1e-6) {
  lo <- min(lo, 0.5 - 1e-8)
  p <- pmax(p, lo)
  p <- pmin(p, 1 - lo)
  p
}

#' Validate an explicit start list (testing hook)
#' @noRd
.check_start <- function(start, K, N) {
  if (!is.list(start) ||
      !all(c("lambda", "p0", "alpha", "beta") %in% names(start)))
    stop("'.start' must be a list with elements lambda, p0, alpha, beta.")
  lambda <- as.numeric(start$lambda)
  p0     <- as.matrix(start$p0)
  alpha  <- as.matrix(start$alpha)
  beta   <- as.matrix(start$beta)
  if (length(lambda) != 2L)
    stop("'.start$lambda' must have length 2.")
  if (!all(dim(p0)    == c(2L, K)) || !all(dim(alpha) == c(2L, K)) ||
      !all(dim(beta)  == c(2L, K)))
    stop("'.start$p0', '$alpha', '$beta' must each be 2 x K matrices.")
  # Match the reference's p0 floor exactly so a shared start is byte-identical.
  p0[p0 < 1/N] <- 1/N
  list(lambda = lambda, p0 = p0, alpha = alpha, beta = beta)
}

#' Initialise all hurdle-gamma parameters for K fields
#'
#' nB-free seeding: split the positive distances of each field into a small
#' "match" tail and a "non-match" bulk according to the resolved match
#' prevalence (lambdaM_init = p0_init), then fit a gamma to each by method of
#' moments (or EnvStats). The non-match hurdle mass starts at 1/N (always valid;
#' the EM raises it if coincidental zeros exist). No nB-based zero-mass formula.
#' @noRd
.init_hurdle_gammaK <- function(X, K, p0_init, p0M_init, N, use_envstats) {
  lambdaM <- p0_init
  lambdaU <- 1 - lambdaM

  # p0M_init may be a per-field vector (data-driven) or recyclable scalar.
  p0M_vec <- rep_len(p0M_init, K)

  p0    <- matrix(0, nrow = 2, ncol = K)
  alpha <- matrix(1, nrow = 2, ncol = K)
  beta  <- matrix(1, nrow = 2, ncol = K)

  for (k in seq_len(K)) {
    xk  <- X[, k]
    xk  <- xk[!is.na(xk)]            # drop NA per field
    x1  <- sort(xk[xk > 0])          # positive distances, ascending
    n_pos <- length(x1)

    # ----- Hurdle masses -----
    # Match: data-driven p0M (per field). Non-match: 1/N (EM raises if needed).
    p0[1, k] <- clamp01_(p0M_vec[k], lo = 1/N)
    p0[2, k] <- 1 / N

    # ----- nB-free gamma seeding -----
    # n_match = expected #matches among the positive distances; take the
    # smallest n_match positive distances as the match component, the rest as
    # the non-match component.
    if (n_pos >= 4L) {
      n_match <- max(2L, round(lambdaM * n_pos))
      n_match <- min(n_match, n_pos - 2L)   # keep >=2 for the non-match set
      xM <- x1[seq_len(n_match)]
      xU <- x1[-seq_len(n_match)]
      paramM <- .fit_gamma_or_unit(xM, use_envstats)
      paramU <- .fit_gamma_or_unit(xU, use_envstats)
      alpha[1, k] <- paramM["shape"]; beta[1, k] <- paramM["scale"]
      alpha[2, k] <- paramU["shape"]; beta[2, k] <- paramU["scale"]
    } else {
      # Too few positives to split: one gamma, mild offset for the other.
      params <- .fit_gamma_or_unit(x1, use_envstats)
      alpha[1, k] <- params["shape"]; beta[1, k] <- params["scale"]
      alpha[2, k] <- max(1, params["shape"] * 0.5)
      beta[2,  k] <- params["scale"] * 2
    }
  }

  list(lambda = c(lambdaM, lambdaU), p0 = p0, alpha = alpha, beta = beta)
}

#' Gamma start fit that returns shape=scale=1 on a degenerate set
#' (var ~ 0 or fewer than 2 positive points).
#' @noRd
.fit_gamma_or_unit <- function(x, use_envstats) {
  x <- x[!is.na(x) & x > 0]
  if (length(x) < 2L || var(x) <= 0) return(c(shape = 1, scale = 1))
  p <- fit_gamma_start(x, use_envstats)
  if (!all(is.finite(p)) || any(p <= 0)) p <- c(shape = 1, scale = 1)
  p
}

#' Perturb a base start for a multistart EM run (uses the package RNG)
#' @noRd
.perturb_start <- function(base, X, K, N, use_envstats) {
  jit <- function(v, lo) {
    v <- v * exp(stats::rnorm(length(v), 0, 0.5))   # log-normal jitter
    pmax(v, lo)
  }
  lambda <- base$lambda
  lambda[1] <- clamp01_(lambda[1] * exp(stats::rnorm(1, 0, 0.5)), lo = 1/N)
  lambda[2] <- 1 - lambda[1]
  p0 <- clamp01_(base$p0 * exp(stats::rnorm(length(base$p0), 0, 0.5)), lo = 1/N)
  list(lambda = lambda,
       p0     = p0,
       alpha  = matrix(jit(base$alpha, 1e-6), nrow = 2),
       beta   = matrix(jit(base$beta,  1e-6), nrow = 2))
}

#' Core EM loop for the K-field hurdle-gamma mixture.
#'
#' Holds the EM *update equations*: E-step posterior, p0 = weighted zero
#' fraction, shape via digamma/uniroot with nlminb/MoM fallback, scale closed
#' form, lambda = colMeans(z). The E-step and the observed-data log-likelihood
#' are computed on the LOG scale (log-sum-exp over the two weighted component
#' log-densities) with per-field NA treated as a neutral log-density 0. This is
#' identical to the original product form to floating point when there are no NA
#' values and no underflow, but (i) keeps the contribution of a row's OBSERVED
#' fields when another field is NA, and (ii) makes the log-likelihood finite
#' under NA/underflow so nstart best-log-likelihood selection works. The p0
#' M-step denominator ranges over rows OBSERVED on each field (numerator and
#' denominator over the same variable). Only initialisation and numerical guards
#' live around these updates.
#' @noRd
.run_em_gammaK <- function(X, K, N, start, tol, maxits, MIN_EFF_N,
                           lambda_smooth = 0) {

  ncomp <- 2L
  lambda_mle <- start$lambda
  p0_mle     <- start$p0
  shape_mle  <- start$alpha
  scale_mle  <- start$beta

  # Clamp p0 from below (matches the reference start exactly).
  p0_mle[p0_mle < 1/N] <- 1/N

  # Per-component, per-field STARTING hurdle masses, used as the shrinkage target
  # when lambda_smooth > 0. Captured AFTER the 1/N floor so the smoothing target
  # is itself a valid probability. With lambda_smooth = 0 this is never read.
  p0_start <- p0_mle

  field_info <- lapply(seq_len(K), function(k) {
    xk <- X[, k]
    ind1 <- which(xk > 0)
    list(
      obs = which(!is.na(xk)),
      zero = which(xk == 0),
      pos = ind1,
      xpos = xk[ind1],
      logxpos = log(xk[ind1])
    )
  })

  # All denominators floored at TINY so a near-collapsed component cannot yield
  # NaN/Inf from a 0-weight division.
  TINY <- 1e-300
  LOG_DENS_FLOOR <- log(TINY)
  fn_alpha <- function(alpha, beta, z, x)
    -log(beta) + sum(z * log(x)) / max(sum(z), TINY) - digamma(alpha)
  fn_alpha2 <- function(alpha, beta, z, x)
    (-log(beta) + sum(z * log(x)) / max(sum(z), TINY) - digamma(alpha))^2
  fn_beta  <- function(z, x, alpha)
    sum(z * x) / (max(sum(z), TINY) * max(alpha, TINY))

  # Compute the N x 2 matrix of summed per-field LOG-densities:
  #   logL[, j] = sum_k log dhgamma(g_k | component j).
  # NA comparison values contribute a neutral log-density 0 (factor 1) PER FIELD,
  # so a row with an NA in one field still uses the information in its observed
  # fields (the previous product form turned the whole row into NA, which was
  # then replaced by the prior, discarding the observed fields).
  .logL_mat <- function(p0, alpha, beta) {
    LL <- matrix(0, nrow = N, ncol = ncomp)
    for (k in seq_len(K)) {
      info <- field_info[[k]]
      for (j in seq_len(ncomp)) {
        ld <- numeric(N)
        if (length(info$zero)) {
          ld[info$zero] <- log(pmax(p0[j, k], TINY))
        }
        if (length(info$pos)) {
          lp <- log1p(-p0[j, k]) +
            dgamma(info$xpos, shape = alpha[j, k], scale = beta[j, k],
                   log = TRUE)
          lp <- pmax(lp, LOG_DENS_FLOOR)
          lp[!is.finite(lp)] <- LOG_DENS_FLOOR
          ld[info$pos] <- lp
        }
        LL[, j] <- LL[, j] + ld
      }
    }
    LL
  }

  # E-step on the LOG scale: posterior q = softmax over the two weighted
  # component log-densities (log-sum-exp), and the observed-data log-likelihood
  # = sum_i logsumexp(a_i, b_i). This is finite even when every per-field density
  # underflows (so nstart log-likelihood selection works), and identical to the
  # product form to floating point when nothing underflows.
  .estep <- function(LL, loglam) {
    a <- LL[, 1L] + loglam[1L]   # log( lambda_M * f_M )
    b <- LL[, 2L] + loglam[2L]   # log( lambda_U * f_U )
    mx <- pmax(a, b)
    denom <- mx + log(exp(a - mx) + exp(b - mx))
    q <- exp(a - denom)
    # Rows where both a and b are -Inf (all densities underflowed) fall back to
    # the prior mixing weight; scrub any residual non-finite posterior too.
    bad <- !is.finite(denom)
    if (any(bad)) q[bad] <- lambda_mle[1]
    q[!is.finite(q)] <- lambda_mle[1]
    list(z = cbind(q, 1 - q),
         loglik = sum(denom[is.finite(denom)]))
  }

  loglam <- log(pmax(lambda_mle, TINY))
  LL <- .logL_mat(p0_mle, shape_mle, scale_mle)

  iter     <- 0L
  converge <- FALSE

  while (!converge && iter < maxits) {
    old_p0    <- p0_mle
    old_shape <- shape_mle
    old_scale <- scale_mle
    old_lam   <- lambda_mle

    # ----- E-step (log scale) -----
    z <- .estep(LL, loglam)$z

    # ----- M-step -----
    for (k in seq_len(K)) {
      info <- field_info[[k]]
      ind0  <- info$zero   # exact-zero, observed (NA excluded by ==)
      ind1  <- info$pos    # positive,    observed (NA excluded by >)
      obs_k <- info$obs

      # Update p0 for each component. The hurdle probability is P(g = 0 | g
      # observed), so BOTH numerator and denominator must range over rows
      # OBSERVED on field k (ind0 already excludes NA; restricting the
      # denominator to obs_k matches the supplement, where numerator and
      # denominator are over the same variable). Floor the denominator so a
      # collapsed component (effective N -> 0) cannot produce 0/0 = NaN; clamp
      # the result to [1/N, 1 - 1/N]. With no NA, obs_k is every row and this is
      # identical to colSums(z). When there are no exact-zero distances (ind0
      # empty) the numerator is 0 and p0 is pinned at its lower floor 1/N.
      p0_num <- colSums(z[ind0,  , drop = FALSE])
      p0_den <- pmax(colSums(z[obs_k, , drop = FALSE]), 1/N)
      # Optional smoothing: shrink each component's hurdle mass towards its
      # STARTING value via a Beta pseudo-count of strength lambda_smooth (the
      # numerator adds lambda_smooth * p0_start, the denominator adds
      # lambda_smooth). This is a purely additive operation, so lambda_smooth = 0
      # reproduces the unsmoothed update p0_num / p0_den exactly. Only the hurdle
      # probabilities are smoothed -- shape/scale are not probabilities and are
      # left untouched.
      p0_mle[, k] <- clamp01_(
        (p0_num + lambda_smooth * p0_start[, k]) / (p0_den + lambda_smooth),
        lo = 1/N)

      if (length(ind1) >= 2) {
        # Update shape via uniroot; fallback to nlminb
        # Guard: only update shape/scale when effective sample size is adequate
        eff_n <- colSums(z[ind1, , drop = FALSE])  # length-ncomp vector

        temp <- old_shape[, k]  # default: keep old values
        sc   <- old_scale[, k]

        for (i in seq_len(ncomp)) {
          if (eff_n[i] >= MIN_EFF_N) {
            zi <- z[ind1, i]
            xi <- info$xpos
            logxi <- info$logxpos
            sw <- max(eff_n[i], TINY)
            sx <- sum(zi * xi)
            slogx <- sum(zi * logxi)
            old_beta_i <- old_scale[i, k]
            fn_alpha_i <- function(alpha)
              -log(old_beta_i) + slogx / sw - digamma(alpha)
            fn_alpha2_i <- function(alpha) fn_alpha_i(alpha)^2
            # Solve the shape MLE equation: uniroot -> nlminb -> weighted
            # method-of-moments -> previous value. Each stage is validated to
            # return a finite, positive shape so the EM never errors here.
            new_shape <- tryCatch(
              uniroot(fn_alpha_i, interval = c(1e-6, 1e4))$root,
              error = function(e) tryCatch(
                nlminb(old_shape[i, k], fn_alpha2_i, lower = 1e-8)$par,
                error = function(e2) NA_real_
              )
            )
            if (!is.finite(new_shape) || new_shape <= 0) {
              # Weighted method-of-moments fallback for the shape.
              mu  <- sx / sw
              v   <- sum(zi * (xi - mu)^2) / sw
              new_shape <- if (is.finite(mu) && is.finite(v) && v > 0)
                mu^2 / v else old_shape[i, k]
            }
            if (!is.finite(new_shape) || new_shape <= 0)
              new_shape <- old_shape[i, k]
            new_scale <- sx / (sw * max(new_shape, TINY))
            if (!is.finite(new_scale) || new_scale <= 0)
              new_scale <- old_scale[i, k]
            temp[i] <- new_shape
            sc[i]   <- new_scale
          }
          # else: leave shape/scale unchanged for this component
        }
        shape_mle[, k] <- temp
        scale_mle[, k] <- sc
      }
      # else: not enough positive observations -- keep old parameters
    }

    # Update mixing proportion. Floor each weight away from 0 (and renormalise)
    # so neither component can collapse to exactly zero weight -- a collapse
    # would zero out that component's density everywhere, drive all row sums to
    # underflow, and ultimately inject NaN into the next M-step. The floor 1/N
    # is the smallest credible prevalence (one record in N pairs).
    lambda_mle <- colMeans(z)
    lambda_mle <- pmax(lambda_mle, 1/N)
    lambda_mle <- lambda_mle / sum(lambda_mle)
    loglam     <- log(pmax(lambda_mle, TINY))

    LL <- .logL_mat(p0_mle, shape_mle, scale_mle)

    # Convergence check: relative change of ALL params < tol.
    # Use max(|old|, tol) in the denominator to guard against near-zero values,
    # and treat any non-finite relative change as "not converged" (isTRUE)
    # rather than letting NA reach the while() condition.
    rel_chg <- function(new_v, old_v) {
      r <- max(abs(new_v - old_v) / pmax(abs(old_v), tol))
      if (!is.finite(r)) Inf else r
    }

    converge <- isTRUE(
      rel_chg(lambda_mle, old_lam)   < tol &&
      rel_chg(p0_mle,     old_p0)    < tol &&
      rel_chg(shape_mle,  old_shape) < tol &&
      rel_chg(scale_mle,  old_scale) < tol
    )

    iter <- iter + 1L
  }

  # Observed-data log-likelihood = sum_i logsumexp(log lambda_M + log f_M,
  # log lambda_U + log f_U). Computed on the log scale so it is finite even with
  # NA fields / underflow (the previous rowSums(dens) form gave -Inf, which broke
  # nstart best-log-likelihood selection).
  loglik <- .estep(LL, loglam)$loglik
  if (!is.finite(loglik)) loglik <- -Inf

  list(lambda = lambda_mle, p0 = p0_mle, shape = shape_mle,
       scale = scale_mle, iter = iter, converged = converge,
       loglik = loglik)
}


# =========================================================================
# S3 methods
# =========================================================================

#' @export
print.problink_em_gamma <- function(x, ...) {
  cat("Hurdle-gamma EM model (problink_em_gamma)\n")
  cat("Fields:", paste(x$by, collapse = ", "), "\n\n")
  K <- length(x$by)
  tab <- data.frame(
    Field     = rep(x$by, each = 2),
    Component = rep(c("match", "non-match"), times = K),
    p0        = as.vector(t(x$p0)),
    shape     = as.vector(t(x$shape)),
    scale     = as.vector(t(x$scale)),
    stringsAsFactors = FALSE
  )
  print(tab, row.names = FALSE, digits = 4)
  if (!is.null(x$p0_init))
    cat("\nInitial match rate (p0, EM start): ", round(x$p0_init, 6), "\n")
  if (!is.null(x$p0M_init)) {
    pm <- x$p0M_init
    cat("Initial exact-match prob (p0M",
        if (length(pm) > 1) " per field" else "", "): ",
        paste(round(pm, 6), collapse = ", "), "\n", sep = "")
  }
  cat("Mixing probability (match): ", round(x$p, 6), "\n")
  if (!is.null(x$nstart) && x$nstart > 1L)
    cat("EM starts (nstart): ", x$nstart, "\n", sep = "")
  cat("Iterations:", x$niter, " Converged:", x$converged, "\n")
  # Hurdle activity diagnostic (Tier 4.2): tell the user whether the point mass
  # at zero is doing any work, or the model has degenerated to a gamma mixture.
  if (!is.null(x$hurdle_active)) {
    if (isTRUE(x$hurdle_active)) {
      cat("Hurdle point mass: ACTIVE (exact matches inform the fit)\n")
    } else {
      cat("Hurdle point mass: inactive -- model degenerated to a plain ",
          "2-component gamma mixture\n", sep = "")
      cat("  (no exact matches in the data; p0 -> 1/N is the correct MLE here)\n")
    }
  }
  invisible(x)
}

#' @export
summary.problink_em_gamma <- function(object, ...) {
  cat("Summary of hurdle-gamma EM model\n")
  cat("===================================\n")
  print(object, ...)
  cat("\nMean distance (shape * scale):\n")
  K <- length(object$by)
  means <- matrix(object$shape * object$scale, nrow = 2,
                  dimnames = list(c("match", "non-match"), object$by))
  print(round(means, 4))
  if (!is.null(object$loglik))
    cat("\nLog-likelihood:", round(object$loglik, 4), "\n")
  invisible(object)
}
