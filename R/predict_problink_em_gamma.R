
#' Calculate weights and posterior match probabilities for continuous pairs
#'
#' @description
#' \code{predict} method for objects of class \code{\link{problink_em_gamma}}.
#' Produces the same output columns as \code{\link{predict.problink_em}} so that
#' downstream functions such as \code{\link{select_threshold}},
#' \code{\link{select_n_to_m}} and \code{\link{link}} work unchanged.
#'
#' @param object an object of class \code{problink_em_gamma} as produced by
#'   \code{\link{problink_em_gamma}}.
#' @param pairs a \code{pairs} object (data.table) containing the comparison
#'   distance columns referenced by \code{object$by}.
#' @param newdata alternative name for \code{pairs}.
#' @param type character scalar specifying what to return:
#'   \describe{
#'     \item{\code{"weights"}}{log-likelihood ratio weight per pair.}
#'     \item{\code{"mpost"}}{posterior probability of match (default).}
#'     \item{\code{"probs"}}{columns \code{mprob}, \code{uprob}, \code{mpost},
#'       \code{upost}.}
#'     \item{\code{"all"}}{all of the above plus \code{weight}.}
#'   }
#' @param add logical; if \code{TRUE} add new columns to \code{pairs}, otherwise
#'   return a data.table with only \code{.x}, \code{.y} and the requested
#'   columns.
#' @param inplace logical; if \code{TRUE} and \code{add = TRUE} modify
#'   \code{pairs} in place (no copy).
#' @param ... unused.
#'
#' @return
#' A \code{data.table} with columns \code{.x}, \code{.y} (and, if
#' \code{add = TRUE}, all columns of \code{pairs}) plus the columns determined
#' by \code{type}.  Pairs with all-\code{NA} comparison distances receive a
#' contribution of 0 to log-weights and 1 to likelihood products (neutral).
#'
#' @examples
#' set.seed(42)
#' nA <- 200; nB <- 40
#' XA <- data.frame(v1 = rexp(nA), v2 = rexp(nA))
#' XB <- XA[1:nB, ] + data.frame(
#'   v1 = rbinom(nB, 1, 0.2) * rexp(nB, 2),
#'   v2 = rbinom(nB, 1, 0.2) * rexp(nB, 2)
#' )
#' library(reclin2)
#' pairs <- pair(XA, XB)
#' pairs <- compare_pairs(pairs, on = c("v1", "v2"),
#'   default_comparator = cmp_distance())
#' model <- problink_em_gamma(~ v1 + v2, data = pairs)
#' res <- predict(model, pairs, type = "mpost", add = FALSE)
#' head(res)
#'
#' \dontshow{gc()}
#'
#' @import data.table
#' @export
predict.problink_em_gamma <- function(object, pairs = newdata, newdata = NULL,
    type = c("weights", "mpost", "probs", "all"),
    add = FALSE, inplace = FALSE, ...) {
  type <- match.arg(type)
  if (is.null(pairs)) pairs <- newdata
  if (is.null(pairs)) stop("Missing pairs or newdata.")
  predict_problinkemgamma(pairs, object, type, add, inplace, ...)
}


predict_problinkemgamma <- function(pairs, model, type, add, inplace, ...) {
  UseMethod("predict_problinkemgamma")
}

#' @import data.table
predict_problinkemgamma.pairs <- function(pairs, model, type, add,
    inplace = FALSE, ...) {

  on <- model$by
  N  <- nrow(pairs)

  # Initialise accumulators
  weights <- rep(0,  N)
  mprobs  <- rep(1,  N)
  uprobs  <- rep(1,  N)

  for (col in on) {
    g  <- pairs[[col]]

    # Densities under match and non-match components
    dM <- dhgamma_(g, model$p0[1, col], model$shape[1, col], model$scale[1, col])
    dU <- dhgamma_(g, model$p0[2, col], model$shape[2, col], model$scale[2, col])

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

  # Build output table
  res <- if (add) {
    if (inplace) pairs else copy(pairs)
  } else {
    if (inplace)
      warning("inplace = TRUE is only relevant when add = TRUE. inplace is ignored.")
    pairs[, list(.x, .y)]
  }

  p <- model$p

  # Posterior on the LOG scale: mpost = plogis(weights + qlogis(p)). This is
  # algebraically identical to mprobs*p / (mprobs*p + uprobs*(1-p)) -- since
  # weights = sum_k (log dM - log dU) = log(mprobs/uprobs) -- but cannot underflow
  # to 0/0 = NaN when there are many fields (mprobs, uprobs -> 0). The mprob/uprob
  # columns keep the raw likelihood products (they may underflow for very many
  # fields); mpost/upost are derived from the log-safe weights.
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
