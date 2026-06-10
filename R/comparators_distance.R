
#' Distance-based comparison functions
#'
#' @description
#' \code{cmp_distance} returns a comparator function that computes the scaled
#' absolute distance \eqn{|x - y| / \text{scale}} between two numeric vectors.
#' This is intended for use with the continuous Fellegi-Sunter model
#' (\code{\link{problink_em_gamma}}) where a distance of \strong{0 indicates an
#' exact match} and larger values indicate greater dissimilarity.
#'
#' \code{cmp_absdist} is a convenience wrapper with \code{scale = 1} and no
#' transformation, equivalent to raw absolute difference.
#'
#' @param scale numeric scalar (default \code{1}). The computed distance
#'   \eqn{|x - y|} is divided by \code{scale} before being returned.  Use this
#'   to normalise variables onto a common scale, e.g., set \code{scale} to the
#'   expected range of the variable.
#' @param transform optional function applied to the (scaled) distance after
#'   computation, e.g., \code{log1p}.  When \code{NULL} (the default), no
#'   transformation is applied.
#'
#' @details
#' Like all comparators in \pkg{reclin2}, the returned function is polymorphic:
#' \itemize{
#'   \item \strong{Two-argument form} \code{cmp(x, y)}: computes the scaled
#'     (and optionally transformed) absolute distance element-wise.  Both
#'     vectors are coerced to \code{numeric} via \code{as.numeric}.
#'   \item \strong{One-argument form} \code{cmp(x)}: binary exact-match
#'     fallback -- returns \code{TRUE} where \code{x == 0} and \code{!is.na(x)}.
#'     This allows the comparator to be used with methods that require binary
#'     comparison vectors, such as the classic \code{\link{problink_em}}.
#' }
#'
#' Missing values (\code{NA}) in either \code{x} or \code{y} propagate as
#' \code{NA} in the output (standard R arithmetic behaviour).
#'
#' @return
#' A function of class \code{"function"} with the polymorphic interface
#' described in Details.
#'
#' @examples
#' cmp <- cmp_distance()
#' # Two-argument form: absolute distance
#' cmp(c(1.0, 2.5, 3.0, NA), c(1.0, 3.0, 1.0, 2.0))
#'
#' # One-argument form: exact-match binary
#' d <- cmp(c(1.0, 2.5, 3.0, NA), c(1.0, 3.0, 1.0, 2.0))
#' cmp(d)
#'
#' # With scaling
#' cmp2 <- cmp_distance(scale = 10)
#' cmp2(c(0, 5, 10), c(0, 0, 0))
#'
#' # cmp_absdist is shorthand for cmp_distance()
#' cmp3 <- cmp_absdist()
#' cmp3(c(3, 7), c(1, 7))
#'
#' \dontshow{gc()}
#'
#' @rdname comparators_distance
#' @export
cmp_distance <- function(scale = 1, transform = NULL) {
  # The hurdle-gamma model divides by 'scale' and expects non-negative distances,
  # so 'scale' must be a single positive finite number.
  if (length(scale) != 1L || !is.numeric(scale) || !is.finite(scale) ||
      scale <= 0)
    stop("'scale' must be a single positive finite number.")
  function(x, y) {
    if (!missing(y)) {
      # Two-argument form: return scaled absolute distance (0 = exact match)
      d <- abs(as.numeric(x) - as.numeric(y)) / scale
      if (!is.null(transform)) {
        d <- transform(d)
        # The continuous Fellegi-Sunter model (problink_em_gamma / the
        # continuous block of problink_em_mixed) requires non-negative distances;
        # a transform that produces negatives (e.g. log on values < 1) breaks it.
        if (any(d < 0, na.rm = TRUE))
          warning("cmp_distance: 'transform' produced negative values; the ",
                  "hurdle-gamma model requires non-negative distances.",
                  call. = FALSE)
      }
      d
    } else {
      # One-argument form: binary exact-match fallback
      !is.na(x) & x == 0
    }
  }
}

#' @rdname comparators_distance
#' @export
cmp_absdist <- function() cmp_distance(scale = 1)


#' Ordinal (binned-distance) comparison function
#'
#' @description
#' \code{cmp_levels} returns a comparator function that discretises the absolute
#' distance \eqn{|x - y|} into an \strong{ordinal level} in
#' \eqn{\{0, 1, \dots, L\}} using a vector of cut points \code{breaks}.  Level
#' \code{0} always denotes an exact agreement (distance \code{0}); higher levels
#' denote increasingly large disagreements.  This produces the multi-category
#' comparison vectors consumed by the categorical block of
#' \code{\link{problink_em_mixed}}.
#'
#' For \code{L} cut points \code{breaks}, there are \code{L + 1} ordinal levels:
#' a distance \eqn{d} maps to level \eqn{\sum_j \mathbf{1}[d > \text{breaks}_j]}.
#' For example \code{breaks = c(0, 3)} reproduces the paper's \code{compare3}:
#' level \code{0} when \eqn{d = 0}, level \code{1} when \eqn{0 < d \le 3}, and
#' level \code{2} when \eqn{d > 3}.
#'
#' @param breaks numeric vector of \strong{ascending} cut points (e.g.
#'   \code{c(0, 3)}).  The first cut point is typically \code{0} so that an exact
#'   match (distance \code{0}) is its own level (level \code{0}).
#'
#' @details
#' Like all comparators in \pkg{reclin2}, the returned function is polymorphic:
#' \itemize{
#'   \item \strong{Two-argument form} \code{cmp(x, y)}: returns the integer
#'     ordinal level of \eqn{|x - y|} as described above.  Both vectors are
#'     coerced to \code{numeric}.  Missing values propagate as \code{NA}.
#'   \item \strong{One-argument form} \code{cmp(x)}: binary agreement fallback --
#'     returns \code{TRUE} where the level is \code{0} (exact agreement) and is
#'     not \code{NA}.  This lets the comparator be used with methods that need a
#'     binary comparison vector, such as the classic \code{\link{problink_em}}.
#' }
#'
#' @return
#' A function with the polymorphic interface described in Details.
#'
#' @examples
#' # Reproduce the paper's 3-level compare3 (breaks = c(0, 3)):
#' cmp <- cmp_levels(c(0, 3))
#' cmp(c(0, 2, 3, 10, NA), c(0, 0, 0, 0, 0))   # -> 0 1 1 2 NA
#'
#' # One-argument form: TRUE only for exact agreement (level 0)
#' lv <- cmp(c(0, 2, 3, 10, NA), c(0, 0, 0, 0, 0))
#' cmp(lv)                                       # -> TRUE FALSE FALSE FALSE FALSE
#'
#' \dontshow{gc()}
#'
#' @rdname comparators_distance
#' @export
cmp_levels <- function(breaks) {
  if (missing(breaks) || length(breaks) < 1 || anyNA(breaks) ||
      !is.numeric(breaks))
    stop("'breaks' must be a non-empty numeric vector of cut points.")
  if (is.unsorted(breaks))
    stop("'breaks' must be in ascending order.")
  function(x, y) {
    if (!missing(y)) {
      # Two-argument form: ordinal level of |x - y| via cut points.
      d <- abs(as.numeric(x) - as.numeric(y))
      # level = number of breaks the distance strictly exceeds (0..length(breaks))
      lv <- rep(0L, length(d))
      na <- is.na(d)
      for (b in breaks) lv <- lv + (d > b)
      lv[na] <- NA_integer_
      as.integer(lv)
    } else {
      # One-argument form: binary agreement (level 0), NA-safe.
      !is.na(x) & x == 0
    }
  }
}
