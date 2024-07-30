#' @title Lp distance on the Lebesgue measure.
#'
#' @description Lp distance to be used for comparing distribution functions.
#'
#' @param pdist1 First distribution function.
#' @param pdist2 Second distribution function.
#' @param p Index of the Lp distance to be used, default is 1 (L1 distance).
#' @param lower Lower bound of the variable, default is -Inf
#' @param upper Upper bound of the variable, default is -Inf
#'
#' @return Lp distance between two functions.
#'
#' @examples
#' Lp_dist(pdist1 = pnorm, pdist2 = function(x) pt(x, df = 10), p = 2)
#'
#' @export
Lp_dist <- function(pdist1, pdist2, p = 1, lower = -Inf, upper = Inf) {
  (stats::integrate(function(x) (abs(pdist1(x) - pdist2(x)))^p,
    lower = lower, upper = upper
  )$value)^(1 / p)
}
