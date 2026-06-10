
#' Mixed-type Fellegi-Sunter EM estimator (binary + categorical + continuous)
#'
#' @description
#' Estimates a two-component Fellegi-Sunter mixture model in which each field
#' may have a different \strong{comparison type}.  Under conditional
#' independence the joint comparison density factorises into a product of
#' per-field blocks, one of three kinds:
#'
#' \describe{
#'   \item{\strong{binary} (Bernoulli)}{the comparison value \eqn{\gamma} is in
#'     \eqn{\{0, 1\}} (agree / disagree).  Parameters: \eqn{m, u} per field.}
#'   \item{\strong{categorical} (multinomial, general \eqn{L} levels)}{the
#'     comparison value \eqn{\gamma} is in \eqn{\{0, \dots, L-1\}} (an ordinal
#'     agreement level, e.g. from \code{\link{cmp_levels}}).  Parameters:
#'     \eqn{m_s, u_s} with \eqn{\sum_s m_s = 1}.  This generalises the paper's
#'     fixed 3-level block to an arbitrary number of observed levels.}
#'   \item{\strong{continuous} (hurdle gamma)}{the comparison value is a
#'     non-negative distance (e.g. from \code{\link{cmp_distance}}); a distance
#'     of \code{0} is an exact match.  Parameters: \eqn{p_0} (point mass at
#'     zero), gamma \eqn{shape} and \eqn{scale}.  This block is identical to the
#'     one fitted by \code{\link{problink_em_gamma}}, so a model in which
#'     \emph{all} fields are continuous reduces \strong{exactly} to that
#'     estimator.}
#' }
#'
#' Component 1 is the \strong{match} component; component 2 is the
#' \strong{non-match} component.  The mixing proportion \eqn{p} is the overall
#' probability that a pair is a match.
#'
#' @param formula a formula of the form \code{~ var1 + var2 + ...} specifying the
#'   comparison columns in \code{data}.
#' @param data a \code{data.table}/\code{data.frame} of comparison values,
#'   typically produced by \code{\link{compare_pairs}}.  Either \code{data} or
#'   \code{comparison_matrix} must be supplied.
#' @param comparison_matrix a numeric matrix whose columns are comparison values
#'   for each field.  An alternative to \code{data}/\code{formula} (e.g. for the
#'   reference-equivalence verification).
#' @param types optional named (or positional) list overriding the per-field
#'   auto-detected comparison type.  Each element must be one of
#'   \code{"binary"}, \code{"categorical"} or \code{"continuous"}, e.g.
#'   \code{types = list(age = "continuous", sex = "binary", dob = "categorical")}.
#'   Fields not named keep their auto-detected type.
#' @param p0 optional scalar initial estimate of the match prevalence (EM start).
#'   When \code{NULL} (default) a data-derived value is used: for a \code{pairs}
#'   object carrying its source tables it is \code{min(nA, nB) / N}; otherwise it
#'   falls back to \code{0.05}.  Clamped to \code{[1/N, 0.5]}.
#' @param p0M optional scalar initial hurdle (exact-match) probability for the
#'   \strong{continuous} blocks (see \code{\link{problink_em_gamma}}).  When
#'   \code{NULL} it is derived per continuous field from the data.
#' @param tol convergence tolerance (relative change of all parameters).
#'   Default \code{1e-6}.
#' @param maxits maximum number of EM iterations.  Default \code{500}.
#' @param nstart integer \eqn{\ge 1}; number of EM starts (default \code{1}).
#'   With \code{nstart = 1} the estimator is deterministic.  With
#'   \code{nstart > 1} the data-driven start is run first and \code{nstart - 1}
#'   perturbed starts follow; the fit with the highest log-likelihood is kept.
#'   Perturbations use the package RNG and respect \code{set.seed()} (the RNG
#'   state is saved and restored on exit).
#' @param lambda_smooth non-negative smoothing strength (default \code{0} = exact
#'   maximum-likelihood updates).  When \code{> 0}, every probability update is
#'   shrunk towards its initial value via a symmetric-Dirichlet pseudo-count:
#'   \eqn{m_s = (\sum q\, \mathbf{1} + \lambda m^0_s) / (\sum q + \lambda)} (and
#'   likewise for \eqn{u_s} and the binary \eqn{m, u}).  Because the target
#'   \eqn{m^0} is itself a simplex (\eqn{\sum_s m^0_s = 1}), the numerator's
#'   pseudo-counts sum to \eqn{\lambda} and the denominator adds exactly
#'   \eqn{\lambda}, so the smoothed categorical \eqn{m_s} (and \eqn{u_s}) still
#'   \strong{sum to 1}.  Useful with sparse levels.
#' @param maxlevels integer; the largest number of distinct non-negative integer
#'   values for which a field auto-detects as \code{"categorical"} (default
#'   \code{10}).  Fields with more distinct integer values auto-detect as
#'   \code{"continuous"}.
#' @param use_envstats logical; passed to the continuous-block initialiser (uses
#'   \code{EnvStats::egamma} for the gamma start when \code{TRUE} and the package
#'   is available).  Default \code{FALSE}.
#' @param .start optional internal/testing hook supplying explicit EM starting
#'   parameters (a list, see Details) and forcing \code{nstart = 1}.  Used to run
#'   this estimator and a reference implementation from byte-identical starts.
#'   Not intended for general use.
#'
#' @details
#' \strong{Type auto-detection.}  For each field's comparison column (non-\code{NA}
#' values \eqn{v}):
#' \itemize{
#'   \item if \eqn{v \subseteq \{0, 1\}} -> \code{"binary"};
#'   \item else if all \eqn{v} are non-negative integers and the number of
#'     distinct values is \eqn{\le} \code{maxlevels} -> \code{"categorical"}
#'     (levels = the sorted unique observed values);
#'   \item else -> \code{"continuous"}.
#' }
#' The inferred type of every field is reported via \code{message()} so it can be
#' overridden.  \strong{The one ambiguous case}: an integer-valued distance with
#' few distinct values (e.g. integer date lags binned to \code{0, 1, 2, ...})
#' auto-detects as \code{categorical}; if it should be treated as a continuous
#' distance, override it explicitly, e.g.
#' \code{types = list(datelag = "continuous")}.
#'
#' \strong{EM.}  The E-step computes the per-pair posterior match probability
#' \eqn{q_i = \Pr(\text{match} \mid \gamma_i)}.  Per-field log-densities are
#' summed on the \strong{log scale} and combined with a log-sum-exp, for
#' underflow safety with many fields.  The M-step updates \eqn{p} as the mean
#' posterior, then each field by its block (closed form for binary and
#' categorical; the same \code{uniroot}/\code{nlminb} shape solver as
#' \code{\link{problink_em_gamma}} for continuous).  A regime-aware
#' label-switching guard relabels so that component 1 is the match component,
#' combining binary (\eqn{m > u}), continuous (smaller mean distance / larger
#' \eqn{p_0}) and categorical (more mass on the exact-agreement level) signals.
#'
#' \strong{Categorical coding assumption.}  The categorical block's
#' initialisation, clamping and label-switching all assume the \emph{smallest}
#' level is the \strong{exact-agreement} level, i.e. \strong{level 0 = exact
#' match}.  This is exactly the contract of \code{\link{cmp_levels}} (which always
#' codes an exact match as level 0), so the recommended workflow -- build
#' categorical comparison columns with \code{cmp_levels} -- is always correct.
#' For an \emph{arbitrary} categorical coding where the smallest code is not the
#' agreement level (e.g. a sum coding such as the paper's
#' \code{temp[, 1] + temp[, 2]} where the \emph{largest} level marks rare-value
#' agreement) this assumption need not hold; a fully coding-agnostic label rule is
#' not well-defined.  A \code{warning} is emitted when a categorical field's
#' smallest observed level is not 0.  In that case either recode the field (or set
#' \code{types =}) so level 0 is agreement, or include an informative
#' binary/continuous field, which anchors the discriminative-power-weighted
#' label-switch so the unambiguous blocks determine the match component.
#'
#' \strong{Explicit start (\code{.start}).}  A list with elements
#' \code{lambda} (length-2 mixing weights), and per-field starting parameters in
#' \code{fields}: a list parallel to the fields, each element a list with the
#' block's parameters -- binary: \code{m}, \code{u} (length-2 component vectors);
#' categorical: \code{levels}, \code{m}, \code{u} (each \eqn{2 \times L});
#' continuous: \code{p0}, \code{alpha}, \code{beta} (each length-2).  See the
#' verification script for an example.
#'
#' @return
#' An object of class \code{problink_em_mixed} (a list) with components:
#' \describe{
#'   \item{\code{fields}}{list of per-field fitted blocks (type + parameters).}
#'   \item{\code{types}}{named character vector of the comparison type per field.}
#'   \item{\code{p}}{scalar mixing probability (probability a pair is a match).}
#'   \item{\code{p_init}}{the initial match prevalence used as the EM start.}
#'   \item{\code{loglik}}{observed-data log-likelihood at the returned fit.}
#'   \item{\code{by}}{character vector of field names.}
#'   \item{\code{niter}, \code{nstart}, \code{converged}}{EM diagnostics.}
#'   \item{\code{blocks_present}}{which of binary/categorical/continuous occur.}
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
#' nA <- 150; nB <- 40
#' # one binary, one 5-level categorical, one continuous field
#' A <- data.frame(
#'   sex = sample(0:1, nA, replace = TRUE),
#'   reg = sample(1:5, nA, replace = TRUE),
#'   age = round(runif(nA, 18, 80))
#' )
#' B <- A[1:nB, ]
#' B$age <- B$age + rbinom(nB, 1, 0.2) * sample(1:5, nB, replace = TRUE)
#' library(reclin2)
#' pairs <- pair(A, B)
#' pairs <- compare_pairs(pairs, on = c("sex", "reg", "age"),
#'   comparators = list(sex = cmp_identical(), reg = cmp_levels(c(0)),
#'                      age = cmp_distance()))
#' model <- problink_em_mixed(~ sex + reg + age, data = pairs)
#' print(model)
#' pred <- predict(model, pairs, type = "mpost")
#'
#' \dontshow{gc()}
#'
#' @importFrom stats dgamma uniroot nlminb var rnorm plogis
#' @export
problink_em_mixed <- function(formula, data, comparison_matrix, types = NULL,
    p0 = NULL, p0M = NULL, tol = 1e-6, maxits = 500, nstart = 1L,
    lambda_smooth = 0, maxlevels = 10L, use_envstats = FALSE, .start = NULL) {

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
  # logical comparison columns (e.g. cmp_identical output) -> 0/1 numeric
  storage.mode(X) <- "double"

  K <- ncol(X)
  N <- nrow(X)

  # A 2-component mixture is not identifiable from a single observation.
  if (N < 2L)
    stop("need at least 2 comparison pairs to fit a 2-component mixture.")

  nstart <- as.integer(nstart)
  if (length(nstart) != 1L || is.na(nstart) || nstart < 1L)
    stop("'nstart' must be a single integer >= 1.")
  if (!is.null(.start)) nstart <- 1L
  if (!is.numeric(lambda_smooth) || length(lambda_smooth) != 1L ||
      !is.finite(lambda_smooth) || lambda_smooth < 0)
    stop("'lambda_smooth' must be a single non-negative number.")
  maxlevels <- as.integer(maxlevels)
  if (length(maxlevels) != 1L || is.na(maxlevels) || maxlevels < 2L)
    stop("'maxlevels' must be a single integer >= 2.")

  # ----- Type auto-detection + override -----
  field_types <- .detect_field_types(X, by = by, types = types,
                                     maxlevels = maxlevels)
  # Tell the user the inferred (post-override) type per field.
  message("problink_em_mixed: comparison field types:\n",
          paste0("  ", by, ": ", field_types, collapse = "\n"))

  # ----- Validate continuous distances -----
  # The hurdle-gamma block is defined on non-negative, finite distances. Apply
  # this check ONLY to columns resolved as "continuous" -- binary (0/1) and
  # categorical (level codes) have their own domains. NA stays allowed.
  cont_k <- which(field_types == "continuous")
  for (k in cont_k) {
    xk <- X[, k]
    if (any(xk < 0, na.rm = TRUE) || any(is.nan(xk)) || any(is.infinite(xk)))
      stop("continuous comparison field '", by[k], "' must have non-negative, ",
           "finite distances (NA allowed); did you pass raw differences ",
           "instead of |x - y|? Use cmp_distance().")
  }

  # ----- Categorical-coding sanity check (see Details/FIX 7) -----
  # The categorical init, clamp and label-switch all treat the SMALLEST observed
  # level as the exact-agreement level -- which is the cmp_levels() contract
  # (level 0 = exact match) but NOT guaranteed for an arbitrary categorical
  # coding. If the smallest observed level is not 0, warn: the "agreement =
  # smallest level" assumption may not hold.
  for (k in which(field_types == "categorical")) {
    lev <- .field_levels(X[, k])
    if (length(lev) && min(lev) != 0)
      warning("categorical field '", by[k], "': the smallest observed level is ",
              min(lev), ", not 0. This estimator assumes level 0 is exact ",
              "agreement (the cmp_levels() contract); for a non-ordinal coding ",
              "this may mislabel components. Use cmp_levels() (which guarantees ",
              "level 0 = agreement), recode the field / set types= explicitly, ",
              "or include an informative binary/continuous field to anchor the ",
              "labels.", call. = FALSE)
  }

  blocks_present <- sort(unique(field_types))

  # ----- Resolve initial match prevalence p0 (EM start) -----
  if (!is.null(p0)) {
    if (length(p0) != 1L || !is.finite(p0) || p0 <= 0 || p0 >= 1)
      stop("'p0' must be a single number in (0, 1).")
  } else {
    p0 <- 0.05
    if (from_data) {
      ax <- attr(data, "x"); ay <- attr(data, "y")
      if (!is.null(ax) && !is.null(ay)) {
        nA <- nrow(ax); nB <- nrow(ay)
        if (length(nA) && length(nB) && nA > 0 && nB > 0)
          p0 <- min(nA, nB) / N
      }
    }
  }
  p0 <- min(max(p0, 1/N), 0.5)
  p_init <- p0

  # ----- Resolve the start (per-field block parameters) -----
  if (!is.null(.start)) {
    base_start <- .check_start_mixed(.start, field_types = field_types,
                                     by = by, N = N)
  } else {
    base_start <- .init_mixed(X, by = by, field_types = field_types,
                              p_init = p_init, p0M = p0M, N = N,
                              use_envstats = use_envstats)
  }

  # ----- Run EM from one or more starts; keep best log-likelihood -----
  if (all(field_types == "continuous")) {
    run_one <- function(start) {
      fit <- .run_em_gammaK(
        X, K = K, N = N, start = .mixed_start_to_gamma(start, by = by),
        tol = tol, maxits = maxits, MIN_EFF_N = 3,
        lambda_smooth = lambda_smooth
      )
      .gamma_fit_to_mixed_result(fit, by = by)
    }
  } else {
    run_one <- function(start)
      .run_em_mixed(X, by = by, field_types = field_types, N = N, start = start,
                    tol = tol, maxits = maxits, lambda_smooth = lambda_smooth)
  }

  best <- run_one(base_start)

  if (nstart > 1L) {
    if (!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
      set_seed(NULL)
    old_seed <- get(".Random.seed", envir = .GlobalEnv)
    on.exit(set_seed(old_seed), add = TRUE)
    for (s in seq_len(nstart - 1L)) {
      pert <- .perturb_start_mixed(base_start, N = N)
      fit  <- tryCatch(run_one(pert), error = function(e) NULL)
      if (!is.null(fit) && is.finite(fit$loglik) && fit$loglik > best$loglik)
        best <- fit
    }
  }

  fields    <- best$fields
  lambda    <- best$lambda
  iter      <- best$iter
  converge  <- best$converged
  loglik    <- best$loglik

  # ----- Label-switching guard: ensure component 1 = MATCH -----
  needs_swap <- .mixed_needs_swap(fields, lambda, tol = tol)
  if (needs_swap) {
    fields <- lapply(fields, .swap_block)
    lambda <- rev(lambda)
  }

  names(field_types) <- by
  names(fields) <- by

  structure(
    list(
      fields        = fields,
      types         = field_types,
      p             = lambda[1],
      p_init        = p_init,
      loglik        = loglik,
      by            = by,
      niter         = iter,
      nstart        = nstart,
      converged     = converge,
      blocks_present = blocks_present,
      lambda_smooth = lambda_smooth,
      use_envstats  = use_envstats,
      call          = cl
    ),
    class = "problink_em_mixed"
  )
}


# =========================================================================
# Internal helpers (file-local, not exported)
#
# These build on the package-internal helpers defined in problink_em_gamma.R
# (dhgamma_, clamp01_, the gamma shape/scale M-step machinery) which are reused
# verbatim for the continuous block so the all-continuous case reduces exactly.
# =========================================================================

#' Auto-detect (and override) per-field comparison types
#' @noRd
.detect_field_types <- function(X, by, types, maxlevels) {
  K <- ncol(X)
  ft <- character(K)
  for (k in seq_len(K)) {
    v <- X[, k]
    v <- v[!is.na(v)]
    if (length(v) == 0L) {
      # Degenerate: the field is entirely NA, so its type cannot be inferred and
      # it carries no information. Classify as binary and WARN so the caller can
      # drop it (its M-step is skipped, leaving the initial block unchanged).
      warning("field '", by[k], "' is entirely NA; it carries no information ",
              "and will not be updated (classified as 'binary'). Consider ",
              "removing it.", call. = FALSE)
      ft[k] <- "binary"   # degenerate: no info, treat as binary
      next
    }
    uv <- unique(v)
    if (all(uv %in% c(0, 1))) {
      ft[k] <- "binary"
    } else if (all(v >= 0) && all(v == floor(v)) &&
               length(uv) <= maxlevels) {
      ft[k] <- "categorical"
    } else {
      ft[k] <- "continuous"
    }
  }
  names(ft) <- by

  # Apply user overrides.
  if (!is.null(types)) {
    if (!is.list(types) && !is.character(types))
      stop("'types' must be a (named) list or character vector.")
    tnames <- names(types)
    valid  <- c("binary", "categorical", "continuous")
    if (is.null(tnames) || any(tnames == "")) {
      # positional: must match number of fields
      if (length(types) != K)
        stop("Unnamed 'types' must have one element per field (", K, ").")
      tnames <- by
    }
    for (i in seq_along(types)) {
      nm <- tnames[i]
      tv <- as.character(types[[i]])
      if (!nm %in% by)
        stop("types: unknown field '", nm, "'.")
      if (length(tv) != 1L || !tv %in% valid)
        stop("types['", nm, "'] must be one of: ",
             paste(valid, collapse = ", "), ".")
      ft[nm] <- tv
    }
  }
  ft
}

#' Sorted unique observed levels of a (possibly NA-containing) field
#' @noRd
.field_levels <- function(x) {
  x <- x[!is.na(x)]
  sort(unique(x))
}

#' Initialise per-field block parameters for the mixed EM
#'
#' Reuses the gamma initialiser (.init_hurdle_gammaK) for continuous fields so
#' the all-continuous start is byte-identical to problink_em_gamma's.
#' @noRd
.init_mixed <- function(X, by, field_types, p_init, p0M, N, use_envstats) {
  K <- ncol(X)
  lambdaM <- p_init
  fields  <- vector("list", K)
  names(fields) <- by

  # Pre-resolve per-field hurdle start for the continuous block, mirroring
  # problink_em_gamma's data-driven p0M (or an explicit scalar override).
  z_k <- colMeans(X == 0, na.rm = TRUE)
  z_k[!is.finite(z_k)] <- 0
  if (!is.null(p0M)) {
    if (length(p0M) != 1L || !is.finite(p0M) || p0M <= 0 || p0M >= 1)
      stop("'p0M' must be a single number in (0, 1).")
    p0M_vec <- rep(p0M, K)
  } else {
    p0M_vec <- clamp01_(z_k / lambdaM, lo = 1/N)
  }
  p0M_vec <- clamp01_(p0M_vec, lo = 1/N)
  names(p0M_vec) <- by

  for (k in seq_len(K)) {
    type <- field_types[k]
    x    <- X[, k]
    if (type == "binary") {
      # m = agreement prob per component c(P(agree|match), P(agree|non-match));
      # u = the complementary disagreement probs. Start: match agrees w.p. 0.9,
      # non-match agrees w.p. 0.1 (the paper's m11=0.9, u11=0.1).
      m_bin <- clamp01_(c(0.9, 0.1), lo = 1/N)
      fields[[k]] <- list(type = "binary", m = m_bin, u = 1 - m_bin)
    } else if (type == "categorical") {
      lev <- .field_levels(x)
      L   <- length(lev)
      # u_s: empirical level frequencies among ALL pairs (non-match proxy).
      tab <- tabulate(match(x[!is.na(x)], lev), nbins = L)
      u_s <- tab / max(sum(tab), 1)
      u_s <- .normalise_simplex(clamp01_(u_s, lo = 1/N))
      # m_s: concentrated on the exact-agreement level (level == min == 0),
      # with light smoothing so every level has positive mass.
      agree_idx <- which.min(lev)         # smallest level == best agreement
      m_s <- rep(0.05 / max(L - 1, 1), L)
      m_s[agree_idx] <- 0.95
      m_s <- .normalise_simplex(clamp01_(m_s, lo = 1/N))
      fields[[k]] <- list(type = "categorical", levels = lev,
                          m = m_s, u = u_s)
    } else { # continuous: reuse the gamma initialiser for a single field
      g1 <- .init_hurdle_gammaK(X[, k, drop = FALSE], K = 1L,
                                p0_init = p_init, p0M_init = p0M_vec[k],
                                N = N, use_envstats = use_envstats)
      fields[[k]] <- list(type = "continuous",
                          p0    = g1$p0[, 1],
                          alpha = g1$alpha[, 1],
                          beta  = g1$beta[, 1])
    }
  }
  list(lambda = c(lambdaM, 1 - lambdaM), fields = fields)
}

#' Normalise a non-negative vector to sum to 1 (safe)
#' @noRd
.normalise_simplex <- function(p) {
  s <- sum(p)
  if (!is.finite(s) || s <= 0) return(rep(1 / length(p), length(p)))
  p / s
}

#' Clamp a multinomial probability vector keeping a chosen residual level.
#'
#' Faithful general-L generalisation of the paper's M2 floor convention: clamp
#' every NON-residual level to [1/N, 1-1/N] (the paper's check1), set the
#' residual level (index \code{ref}, the exact-agreement level) to
#' \eqn{1 - \sum(\text{others})}, and -- if that residual falls below 1/N --
#' shave the deficit equally from the non-residual levels so the residual sits at
#' 1/N. For \code{L = 3} this reproduces M2's
#' \code{m20 = 1 - m22 - m21} plus its \code{1/(2N)} shave exactly. When no level
#' hits a bound, the input (already summing to 1) is returned unchanged.
#' @noRd
.multinom_clamp <- function(p, ref, N) {
  L <- length(p)
  if (L == 1L) return(1)
  nonref <- setdiff(seq_len(L), ref)
  p[nonref] <- clamp01_(p[nonref], lo = 1/N)
  resid <- 1 - sum(p[nonref])
  if (resid < 1/N) {
    deficit <- (1/N) - resid
    p[nonref] <- p[nonref] - deficit / length(nonref)
    resid <- 1/N
  }
  p[ref] <- resid
  p
}

#' Validate an explicit mixed start (testing hook)
#' @noRd
.check_start_mixed <- function(start, field_types, by, N) {
  if (!is.list(start) || !all(c("lambda", "fields") %in% names(start)))
    stop("'.start' must be a list with elements 'lambda' and 'fields'.")
  lambda <- as.numeric(start$lambda)
  if (length(lambda) != 2L) stop("'.start$lambda' must have length 2.")
  fl <- start$fields
  K  <- length(by)
  if (length(fl) != K)
    stop("'.start$fields' must have one element per field (", K, ").")
  out <- vector("list", K); names(out) <- by
  for (k in seq_len(K)) {
    type <- field_types[k]
    f    <- fl[[k]]
    if (type == "binary") {
      # m = agreement prob per component c(m11, u11) (length 2). u is derived as
      # 1 - m, so only m is required (u, if supplied, is ignored).
      m <- as.numeric(f$m)
      if (length(m) != 2L)
        stop("binary start for '", by[k], "' needs length-2 m (c(m11, u11)).")
      mb <- clamp01_(m, lo = 1/N)
      out[[k]] <- list(type = "binary", m = mb, u = 1 - mb)
    } else if (type == "categorical") {
      lev <- as.numeric(f$levels)
      m   <- as.numeric(f$m); u <- as.numeric(f$u)
      L   <- length(lev)
      if (length(m) != L || length(u) != L)
        stop("categorical start for '", by[k], "' needs length-L m and u.")
      out[[k]] <- list(type = "categorical", levels = lev,
                       m = .normalise_simplex(clamp01_(m, lo = 1/N)),
                       u = .normalise_simplex(clamp01_(u, lo = 1/N)))
    } else {
      p0 <- as.numeric(f$p0); a <- as.numeric(f$alpha); b <- as.numeric(f$beta)
      if (length(p0) != 2L || length(a) != 2L || length(b) != 2L)
        stop("continuous start for '", by[k], "' needs length-2 p0/alpha/beta.")
      p0[p0 < 1/N] <- 1/N
      out[[k]] <- list(type = "continuous", p0 = p0, alpha = a, beta = b)
    }
  }
  list(lambda = lambda, fields = out)
}

#' Perturb a mixed start for a multistart run (uses the package RNG)
#' @noRd
.perturb_start_mixed <- function(base, N) {
  lambda <- base$lambda
  lambda[1] <- clamp01_(lambda[1] * exp(stats::rnorm(1, 0, 0.5)), lo = 1/N)
  lambda[2] <- 1 - lambda[1]
  jit_prob <- function(v) clamp01_(v * exp(stats::rnorm(length(v), 0, 0.5)),
                                   lo = 1/N)
  jit_pos  <- function(v) pmax(v * exp(stats::rnorm(length(v), 0, 0.5)), 1e-6)
  fields <- lapply(base$fields, function(f) {
    if (f$type == "binary") {
      f$m <- jit_prob(f$m); f$u <- 1 - f$m
    } else if (f$type == "categorical") {
      f$m <- .normalise_simplex(jit_prob(f$m))
      f$u <- .normalise_simplex(jit_prob(f$u))
    } else {
      f$p0    <- jit_prob(f$p0)
      f$alpha <- jit_pos(f$alpha)
      f$beta  <- jit_pos(f$beta)
    }
    f
  })
  list(lambda = lambda, fields = fields)
}

#' Convert an all-continuous mixed start to the gamma EM start format
#' @noRd
.mixed_start_to_gamma <- function(start, by) {
  K <- length(by)
  p0 <- alpha <- beta <- matrix(NA_real_, nrow = 2L, ncol = K)
  for (k in seq_len(K)) {
    f <- start$fields[[k]]
    p0[, k]    <- f$p0
    alpha[, k] <- f$alpha
    beta[, k]  <- f$beta
  }
  list(lambda = start$lambda, p0 = p0, alpha = alpha, beta = beta)
}

#' Convert a gamma EM fit to the mixed EM result format
#' @noRd
.gamma_fit_to_mixed_result <- function(fit, by) {
  K <- length(by)
  fields <- vector("list", K)
  names(fields) <- by
  for (k in seq_len(K)) {
    fields[[k]] <- list(
      type = "continuous",
      p0 = fit$p0[, k],
      alpha = fit$shape[, k],
      beta = fit$scale[, k]
    )
  }
  list(lambda = fit$lambda, fields = fields, iter = fit$iter,
       converged = fit$converged, loglik = fit$loglik)
}

#' Per-field log-density of one component (length-N), NA -> 0 (neutral)
#'
#' Returns log f_j(gamma_k) for component j of field k. NA comparison values get
#' log-density 0 (a neutral factor of 1 in the product / weight).
#' @noRd
.logdens_field <- function(x, field, comp) {
  type <- field$type
  if (type == "binary") {
    m <- field$m[comp]
    # gamma in {0,1}: m^gamma (1-m)^(1-gamma)
    ld <- ifelse(x == 1, log(m), log1p(-m))
  } else if (type == "categorical") {
    lev <- field$levels
    pr  <- if (comp == 1L) field$m else field$u
    idx <- match(x, lev)              # NA for unseen / NA values
    ld  <- log(pr[idx])
  } else { # continuous
    d  <- dhgamma_(x, field$p0[comp], field$alpha[comp], field$beta[comp])
    ld <- log(pmax(d, 1e-300))
  }
  # Guard -Inf from a zero probability (e.g. an unseen categorical level),
  # then set genuine NA comparison values to the neutral log-density 0.
  ld[!is.finite(ld)] <- log(1e-300)
  ld[is.na(x)] <- 0
  ld
}

#' Core EM loop for the mixed-type mixture.
#'
#' E-step combines per-field log-densities (log-sum-exp over the 2 components)
#' for the posterior. M-step updates p (mean posterior) then each field by block:
#'   - binary / categorical: closed form (optionally smoothed by lambda_smooth);
#'   - continuous: the SAME uniroot/nlminb shape solver and closed-form scale as
#'     problink_em_gamma (so the all-continuous case reduces exactly).
#' @noRd
.run_em_mixed <- function(X, by, field_types, N, start, tol, maxits,
                          lambda_smooth) {
  K      <- ncol(X)
  lambda <- start$lambda
  fields <- start$fields
  # initial values (for smoothing shrinkage targets)
  fields0 <- fields

  field_info <- lapply(seq_len(K), function(k) {
    x <- X[, k]
    obs <- which(!is.na(x))
    f <- fields[[k]]
    if (f$type == "binary") {
      list(obs = obs, agree = x[obs] == 1)
    } else if (f$type == "categorical") {
      list(obs = obs, code = match(x[obs], f$levels))
    } else {
      pos <- which(x > 0)
      list(obs = obs, zero = which(x == 0), pos = pos,
           xpos = x[pos], logxpos = log(x[pos]))
    }
  })

  TINY <- 1e-300
  LOG_DENS_FLOOR <- log(TINY)
  loglam <- log(pmax(lambda, TINY))

  # ---- log-density bookkeeping: logL[, j] = sum_k log f_j(gamma_k) ----
  .logL_mat <- function(fields) {
    LL <- matrix(0, nrow = N, ncol = 2L)
    for (k in seq_len(K)) {
      f <- fields[[k]]
      info <- field_info[[k]]
      if (f$type == "binary") {
        if (length(info$obs)) {
          for (j in 1:2) {
            ld <- ifelse(info$agree, log(f$m[j]), log1p(-f$m[j]))
            ld[!is.finite(ld)] <- LOG_DENS_FLOOR
            LL[info$obs, j] <- LL[info$obs, j] + ld
          }
        }
      } else if (f$type == "categorical") {
        if (length(info$obs)) {
          valid <- !is.na(info$code)
          for (j in 1:2) {
            pr <- if (j == 1L) f$m else f$u
            ld <- rep(LOG_DENS_FLOOR, length(info$obs))
            if (any(valid)) ld[valid] <- log(pr[info$code[valid]])
            ld[!is.finite(ld)] <- LOG_DENS_FLOOR
            LL[info$obs, j] <- LL[info$obs, j] + ld
          }
        }
      } else {
        for (j in 1:2) {
          ld <- numeric(N)
          if (length(info$zero)) {
            ld[info$zero] <- log(pmax(f$p0[j], TINY))
          }
          if (length(info$pos)) {
            # Inlined log dgamma(g; shape, scale) using the cached log(g)
            # (see .run_em_gammaK): avoids recomputing log(g) inside
            # dgamma(log = TRUE) every iteration, the EM's dominant cost.
            lp <- log1p(-f$p0[j]) +
              (f$alpha[j] - 1) * info$logxpos - info$xpos / f$beta[j] -
              (f$alpha[j] * log(f$beta[j]) + lgamma(f$alpha[j]))
            lp <- pmax(lp, LOG_DENS_FLOOR)
            lp[!is.finite(lp)] <- LOG_DENS_FLOOR
            ld[info$pos] <- lp
          }
          LL[, j] <- LL[, j] + ld
        }
      }
    }
    LL
  }

  # log-sum-exp of the two weighted component log-densities, per row
  .estep <- function(LL, loglam) {
    a <- LL[, 1L] + loglam[1L]   # log( lambdaM * fM )
    b <- LL[, 2L] + loglam[2L]   # log( lambdaU * fU )
    mx <- pmax(a, b)
    # posterior q = exp(a) / (exp(a) + exp(b)); numerically via mx
    denom <- mx + log(exp(a - mx) + exp(b - mx))
    q <- exp(a - denom)
    # rows where both a and b are -Inf (all densities underflowed): fall back
    bad <- !is.finite(denom)
    if (any(bad)) q[bad] <- lambda[1]
    q[!is.finite(q)] <- lambda[1]
    list(q = q, loglik = sum(denom[is.finite(denom)]))
  }

  # ---- gamma shape/scale M-step (same equations as problink_em_gamma; the
  # per-field closures fn_alpha_i/fn_alpha2_i are built inline in the loop) ----
  MIN_EFF_N <- 3

  LL <- .logL_mat(fields)
  iter <- 0L; converge <- FALSE; loglik <- -Inf
  # Track which all-NA fields we have already warned about (warn once per run).
  warned_empty <- logical(K)

  # helper to flatten all params for the convergence check
  .flatten <- function(lambda, fields) {
    v <- lambda
    for (f in fields) {
      if (f$type == "binary")        v <- c(v, f$m, f$u)
      else if (f$type == "categorical") v <- c(v, as.vector(f$m), as.vector(f$u))
      else                            v <- c(v, f$p0, f$alpha, f$beta)
    }
    v
  }

  while (!converge && iter < maxits) {
    old_flat <- .flatten(lambda, fields)

    # ----- E-step -----
    es <- .estep(LL, loglam)
    q  <- es$q
    sum_q  <- sum(q)
    sum_1q <- N - sum_q

    # ----- M-step: mixing proportion -----
    lambda <- c(sum_q / N, sum_1q / N)
    lambda <- pmax(lambda, 1/N)
    lambda <- lambda / sum(lambda)
    loglam <- log(pmax(lambda, TINY))

    # ----- M-step: per-field blocks -----
    for (k in seq_len(K)) {
      f <- fields[[k]]
      f0 <- fields0[[k]]
      info <- field_info[[k]]
      obs <- info$obs
      qk  <- q[obs]
      oneqk <- 1 - qk
      sqk  <- sum(qk); s1qk <- sum(oneqk)

      # A field with NO observed rows (entirely NA) carries no information: skip
      # its update so the initial block stays unchanged (no NaN from 0/0). Warn
      # at most once per field for this EM run.
      if (length(obs) == 0L) {
        if (!warned_empty[k]) {
          warning("field '", by[k], "' has no observed (non-NA) comparison ",
                  "values; its parameters are left at their initial values.",
                  call. = FALSE)
          warned_empty[k] <- TRUE
        }
        next
      }

      if (f$type == "binary") {
        # Agreement probabilities per component:
        #   m11 = P(agree | match)     = sum(q * 1[gamma=1]) / sum(q)
        #   u11 = P(agree | non-match) = sum((1-q) * 1[gamma=1]) / sum(1-q)
        # f$m holds c(m11, u11) (agreement prob per component); f$u holds the
        # complementary disagreement probs c(1-m11, 1-u11). The shrinkage target
        # f0$m is the length-2 initial agreement-prob vector. Floor each
        # denominator away from zero (a component can have zero weight on the
        # observed rows even when the field itself is observed) so no 0/0 occurs.
        agree <- as.numeric(info$agree)
        m11 <- (sum(qk * agree)   + lambda_smooth * f0$m[1]) /
               max(sqk            + lambda_smooth, TINY)
        u11 <- (sum(oneqk * agree) + lambda_smooth * f0$m[2]) /
               max(s1qk           + lambda_smooth, TINY)
        f$m <- clamp01_(c(m11, u11), lo = 1/N)
        f$u <- 1 - f$m

      } else if (f$type == "categorical") {
        lev <- f$levels; L <- length(lev)
        idx <- info$code
        valid <- !is.na(idx)
        # weighted counts per level for match (q) and non-match (1-q)
        wm <- numeric(L); wu <- numeric(L)
        if (any(valid)) {
          wm_tab <- rowsum(qk[valid], idx[valid], reorder = FALSE)
          wu_tab <- rowsum(oneqk[valid], idx[valid], reorder = FALSE)
          wm[as.integer(rownames(wm_tab))] <- wm_tab[, 1]
          wu[as.integer(rownames(wu_tab))] <- wu_tab[, 1]
        }
        # Multinomial MLE: weighted relative frequency per level (this already
        # sums to 1). Optional smoothing is a proper symmetric-Dirichlet
        # pseudo-count update: the numerator adds lambda_smooth * f0 (which sums
        # to lambda_smooth, since f0 is a simplex) and the denominator adds
        # lambda_smooth (NOT lambda_smooth * L), so the smoothed vector still
        # SUMS TO 1. Floor the denominator away from zero so a component with no
        # weight on the observed rows cannot produce 0/0.
        m_s <- (wm + lambda_smooth * f0$m) / max(sqk  + lambda_smooth, TINY)
        u_s <- (wu + lambda_smooth * f0$u) / max(s1qk + lambda_smooth, TINY)
        # Floor/cap each level, treating the exact-agreement (smallest) level as
        # the residual so the vector still sums to exactly 1 -- the faithful
        # general-L generalisation of the paper's M2 (where m20 = 1-m22-m21 and a
        # below-1/N residual is corrected by shaving the other levels). This
        # matches reference behaviour exactly when clamping is inactive AND when a
        # rare level hits the 1/N floor.
        f$m <- .multinom_clamp(m_s, ref = which.min(lev), N = N)
        f$u <- .multinom_clamp(u_s, ref = which.min(lev), N = N)

      } else { # continuous: identical update to problink_em_gamma
        ind0 <- info$zero
        ind1 <- info$pos
        # p0 per component
        p0_num <- c(sum(q[ind0]), sum(1 - q[ind0]))
        p0_den <- pmax(c(sum(q[obs]), sum(1 - q[obs])), 1/N)
        f$p0 <- clamp01_(p0_num / p0_den, lo = 1/N)
        if (length(ind1) >= 2) {
          eff_n <- c(sum(q[ind1]), sum(1 - q[ind1]))
          new_alpha <- f$alpha; new_beta <- f$beta
          for (i in 1:2) {
            if (eff_n[i] >= MIN_EFF_N) {
              zi <- if (i == 1L) q[ind1] else 1 - q[ind1]
              xi <- info$xpos
              logxi <- info$logxpos
              sw <- max(eff_n[i], TINY)
              sx <- sum(zi * xi)
              slogx <- sum(zi * logxi)
              old_beta_i <- f$beta[i]
              fn_alpha_i <- function(alpha)
                -log(old_beta_i) + slogx / sw - digamma(alpha)
              fn_alpha2_i <- function(alpha) fn_alpha_i(alpha)^2
              ns <- tryCatch(
                uniroot(fn_alpha_i, interval = c(1e-6, 1e4))$root,
                error = function(e) tryCatch(
                  nlminb(f$alpha[i], fn_alpha2_i, lower = 1e-8)$par,
                  error = function(e2) NA_real_))
              if (!is.finite(ns) || ns <= 0) {
                mu <- sx / sw
                v  <- sum(zi * (xi - mu)^2) / sw
                ns <- if (is.finite(mu) && is.finite(v) && v > 0) mu^2 / v
                      else f$alpha[i]
              }
              if (!is.finite(ns) || ns <= 0) ns <- f$alpha[i]
              nb <- sx / (sw * max(ns, TINY))
              if (!is.finite(nb) || nb <= 0) nb <- f$beta[i]
              new_alpha[i] <- ns; new_beta[i] <- nb
            }
          }
          f$alpha <- new_alpha; f$beta <- new_beta
        }
      }
      fields[[k]] <- f
    }

    LL <- .logL_mat(fields)

    # convergence: relative change of ALL params
    new_flat <- .flatten(lambda, fields)
    rel <- max(abs(new_flat - old_flat) / pmax(abs(old_flat), tol))
    if (!is.finite(rel)) rel <- Inf
    converge <- isTRUE(rel < tol)
    iter <- iter + 1L
  }

  es <- .estep(LL, loglam)
  loglik <- es$loglik
  if (!is.finite(loglik)) loglik <- -Inf

  list(lambda = lambda, fields = fields, iter = iter, converged = converge,
       loglik = loglik)
}

#' Decide whether components must be swapped so component 1 = MATCH
#'
#' Combines signals across the present block types into a single signed,
#' \strong{discriminative-power-weighted} score. Each field contributes
#' \code{weight * direction} where \code{direction > 0} means "component 1 already
#' looks like the match" (keep) and \code{< 0} means "swap"; \code{weight} is the
#' size of that field's match/non-match separation, so an uninformative field
#' (m ~ u, weight ~ 0) cannot outvote a discriminating one. Signals:
#' \itemize{
#'   \item binary: agreement prob, direction = sign(m1 - m2), weight = |m1 - m2|;
#'   \item categorical: mass on the exact-agreement (smallest) level,
#'     direction = sign(m_agree - u_agree), weight = |m_agree - u_agree|;
#'   \item continuous: hurdle mass p0 (weight |p01 - p02|) plus a mean-distance
#'     term (smaller distance = match), weight = |d1 - d2| / (|d1| + |d2|).
#' }
#' A total score >= 0 keeps the labelling; < 0 swaps. On an exact tie (no
#' informative field) the smaller mixing component is taken as the match.
#' @noRd
.mixed_needs_swap <- function(fields, lambda, tol) {
  score <- 0  # > 0 => keep (comp 1 = match); < 0 => swap
  for (f in fields) {
    if (f$type == "binary") {
      score <- score + (f$m[1] - f$m[2])
    } else if (f$type == "categorical") {
      ai <- which.min(f$levels)   # exact-agreement level
      score <- score + (f$m[ai] - f$u[ai])
    } else {
      p01 <- f$p0[1]; p02 <- f$p0[2]
      d1  <- f$alpha[1] * f$beta[1]; d2 <- f$alpha[2] * f$beta[2]
      # hurdle-mass signal (match has the larger point mass at distance 0)
      score <- score + (p01 - p02)
      # mean-distance signal (match is the SMALLER mean distance), normalised so
      # it is on a comparable [-1, 1] scale to the probability signals.
      dscale <- abs(d1) + abs(d2)
      if (is.finite(dscale) && dscale > 0)
        score <- score + (d2 - d1) / dscale
    }
  }
  if (abs(score) > tol) return(score < 0)
  # Exact tie (no informative field): the match component is the smaller one.
  isTRUE(lambda[2] < lambda[1])
}

#' Swap the two components of a fitted block (match <-> non-match)
#' @noRd
.swap_block <- function(f) {
  if (f$type == "binary") {
    f$m <- f$m[c(2, 1)]; f$u <- f$u[c(2, 1)]
  } else if (f$type == "categorical") {
    # m = match level-probs, u = non-match level-probs: swap the two vectors.
    tmp <- f$m; f$m <- f$u; f$u <- tmp
  } else {
    f$p0 <- f$p0[c(2, 1)]; f$alpha <- f$alpha[c(2, 1)]; f$beta <- f$beta[c(2, 1)]
  }
  f
}


# =========================================================================
# S3 methods
# =========================================================================

#' @export
print.problink_em_mixed <- function(x, ...) {
  cat("Mixed-type Fellegi-Sunter EM model (problink_em_mixed)\n")
  cat("Fields:", paste0(x$by, " [", x$types, "]", collapse = ", "), "\n")
  cat("Blocks present:", paste(x$blocks_present, collapse = ", "), "\n\n")
  for (nm in x$by) {
    f <- x$fields[[nm]]
    cat("Field '", nm, "' (", f$type, "):\n", sep = "")
    if (f$type == "binary") {
      tab <- data.frame(Component = c("match", "non-match"),
                        m = f$m, u = f$u)
      print(tab, row.names = FALSE, digits = 4)
    } else if (f$type == "categorical") {
      mm <- rbind(level = f$levels, m_match = f$m, u_nonmatch = f$u)
      colnames(mm) <- paste0("L", f$levels)
      print(round(mm, 4))
    } else {
      tab <- data.frame(Component = c("match", "non-match"),
                        p0 = f$p0, shape = f$alpha, scale = f$beta)
      print(tab, row.names = FALSE, digits = 4)
    }
    cat("\n")
  }
  if (!is.null(x$p_init))
    cat("Initial match rate (p, EM start): ", round(x$p_init, 6), "\n")
  cat("Mixing probability (match): ", round(x$p, 6), "\n")
  if (!is.null(x$lambda_smooth) && x$lambda_smooth > 0)
    cat("Smoothing (lambda_smooth): ", x$lambda_smooth, "\n", sep = "")
  if (!is.null(x$nstart) && x$nstart > 1L)
    cat("EM starts (nstart): ", x$nstart, "\n", sep = "")
  cat("Log-likelihood:", round(x$loglik, 4), "\n")
  cat("Iterations:", x$niter, " Converged:", x$converged, "\n")
  invisible(x)
}

#' @export
summary.problink_em_mixed <- function(object, ...) {
  cat("Summary of mixed-type Fellegi-Sunter EM model\n")
  cat("=============================================\n")
  print(object, ...)
  invisible(object)
}
