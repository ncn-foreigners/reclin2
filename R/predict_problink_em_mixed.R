
#' Calculate weights and posterior match probabilities for mixed-type pairs
#'
#' @description
#' \code{predict} method for objects of class \code{\link{problink_em_mixed}}.
#' Produces the same output columns as \code{\link{predict.problink_em}} and
#' \code{\link{predict.problink_em_gamma}} so that downstream functions such as
#' \code{\link{select_threshold}}, \code{\link{select_n_to_m}} and
#' \code{\link{link}} work unchanged.
#'
#' Each field contributes its per-type log-likelihood ratio \eqn{\log(d_M/d_U)}
#' to the weight, and its density ratio to the m-/u-probability products:
#' \itemize{
#'   \item \strong{binary}: \eqn{m^{\gamma}(1-m)^{1-\gamma}} for each component;
#'   \item \strong{categorical}: \eqn{m_{\gamma}} resp. \eqn{u_{\gamma}};
#'   \item \strong{continuous}: the hurdle-gamma density (the same as
#'     \code{\link{predict.problink_em_gamma}}).
#' }
#' Pairs with \code{NA} comparison values receive a neutral contribution
#' (weight \code{0}, density factor \code{1}).
#'
#' @param object an object of class \code{problink_em_mixed}.
#' @param pairs a \code{pairs} object (data.table) containing the comparison
#'   columns referenced by \code{object$by}.
#' @param newdata alternative name for \code{pairs}.
#' @param type character scalar: \code{"weights"}, \code{"mpost"} (default),
#'   \code{"probs"} (\code{mprob}, \code{uprob}, \code{mpost}, \code{upost}) or
#'   \code{"all"}.
#' @param add logical; if \code{TRUE} add new columns to \code{pairs}, otherwise
#'   return a data.table with only \code{.x}, \code{.y} and the requested
#'   columns.
#' @param inplace logical; if \code{TRUE} and \code{add = TRUE} modify
#'   \code{pairs} in place (no copy).
#' @param ... unused.
#'
#' @return
#' A \code{data.table} with columns \code{.x}, \code{.y} (and, if
#' \code{add = TRUE}, all columns of \code{pairs}) plus the columns determined by
#' \code{type}.
#'
#' @examples
#' set.seed(1)
#' nA <- 200; nB <- 40
#' A <- data.frame(sex = sample(0:1, nA, TRUE), age = round(runif(nA, 18, 80)))
#' B <- A[1:nB, ]; B$age <- B$age + rbinom(nB, 1, 0.2) * 3
#' library(reclin2)
#' pairs <- pair(A, B)
#' pairs <- compare_pairs(pairs, on = c("sex", "age"),
#'   comparators = list(sex = cmp_identical(), age = cmp_distance()))
#' model <- problink_em_mixed(~ sex + age, data = pairs)
#' res <- predict(model, pairs, type = "mpost", add = FALSE)
#' head(res)
#'
#' \dontshow{gc()}
#'
#' @import data.table
#' @export
predict.problink_em_mixed <- function(object, pairs = newdata, newdata = NULL,
    type = c("weights", "mpost", "probs", "all"),
    add = FALSE, inplace = FALSE, ...) {
  type <- match.arg(type)
  if (is.null(pairs)) pairs <- newdata
  if (is.null(pairs)) stop("Missing pairs or newdata.")
  predict_problinkemmixed(pairs, object, type, add, inplace, ...)
}


predict_problinkemmixed <- function(pairs, model, type, add, inplace, ...) {
  UseMethod("predict_problinkemmixed")
}

#' @import data.table
predict_problinkemmixed.pairs <- function(pairs, model, type, add,
    inplace = FALSE, ...) {

  on <- model$by
  N  <- nrow(pairs)

  weights <- rep(0, N)
  mprobs  <- rep(1, N)
  uprobs  <- rep(1, N)

  for (col in on) {
    g <- pairs[[col]]
    f <- model$fields[[col]]

    if (f$type == "binary") {
      gg <- as.numeric(g)
      dM <- ifelse(gg == 1, f$m[1], 1 - f$m[1])
      dU <- ifelse(gg == 1, f$m[2], 1 - f$m[2])
    } else if (f$type == "categorical") {
      # f$m = match level-probs, f$u = non-match level-probs (length L).
      idx <- match(g, f$levels)            # NA for NA / unseen levels
      dM  <- f$m[idx]
      dU  <- f$u[idx]
    } else { # continuous: identical to predict.problink_em_gamma
      dM <- dhgamma_(g, f$p0[1], f$alpha[1], f$beta[1])
      dU <- dhgamma_(g, f$p0[2], f$alpha[2], f$beta[2])
    }

    # Log-weight contribution; NA pairs get weight 0 (neutral)
    w <- log(pmax(dM, 1e-300)) - log(pmax(dU, 1e-300))
    w[is.na(w)] <- 0

    # Density products; NA pairs get factor 1 (neutral)
    dM[is.na(dM)] <- 1
    dU[is.na(dU)] <- 1

    weights <- weights + w
    mprobs  <- mprobs * dM
    uprobs  <- uprobs * dU
  }

  res <- if (add) {
    if (inplace) pairs else copy(pairs)
  } else {
    if (inplace)
      warning("inplace = TRUE is only relevant when add = TRUE. inplace is ignored.")
    pairs[, list(.x, .y)]
  }

  p <- model$p

  # Posterior on the LOG scale: mpost = plogis(weights + qlogis(p)). Identical to
  # mprobs*p / (mprobs*p + uprobs*(1-p)) (weights = log(mprobs/uprobs)) but cannot
  # underflow to NaN with many fields. mprob/uprob keep the raw likelihood
  # products (may underflow); mpost/upost are derived from the log-safe weights.
  mpost <- stats::plogis(weights + log(p) - log1p(-p))

  if (type == "weights") {
    res[, weights := weights]
  } else if (type == "mpost") {
    res[, mpost := mpost]
  } else {
    res[, mprob := mprobs]
    res[, uprob := uprobs]
    res[, mpost := mpost]
    res[, upost := 1 - mpost]
    if (type == "all") res[, weight := weights]
  }

  if (inplace) invisible(res) else res
}
