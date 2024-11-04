#' @title Two-sample test
#'
#' @description PROTEST version of the two-sample test. For two list of
#' distribution functions, one for each datasets, it tests the hypothesis that
#' the true distribution function is the same for both sets.
#'
#' @param p1_list List of distribution functions based on the first variable.
#' @param p2_list List of distribution functions based on the first variable.
#' @param dist Distance function to be applied, default is the L1 distance.
#' @param alpha Significance level, default is 0.05. When given two values, a
#' three-way test is performed.
#' @param epsilon Threshold value.
#' @param verbose Should the decision be returned as a message? Default if TRUE.
#' @param plot Should the test return a plot of the decision for each
#' (alpha, epsilon)? Default is FALSE
#'
#' @return Test results of the two-sample test.
#'
#' @examples
#' p1_list <- list(
#'   pexp,
#'   function(x) pgamma(x, shape = 2),
#'   function(x) pexp(x, rate = 2)
#' )
#'
#' p2_list <- list(
#'   function(x) pexp(x, rate = 1.1),
#'   function(p) pgamma(p, shape = 2.1),
#'   function(p) pexp(p, rate = 2.1)
#' )
#'
#' my_dist <- function(pdist1, pdist2) {
#'   protest::Lp_dist(pdist1, pdist2, p = 2, lower = 0)
#' }
#'
#' twosample_test(p1_list, p2_list,
#'   dist = my_dist, alpha = 0.05, epsilon = 0.1,
#'   verbose = TRUE, plot = TRUE
#' )
#'
#' @export
twosample_test <- function(
    p1_list, p2_list, epsilon,
    dist = protest::Lp_dist,
    alpha = 0.05, verbose = TRUE, plot = FALSE) {
  dissim <- sapply(
    p1_list,
    function(p1) {
      sapply(
        p2_list,
        function(p2) dist(p1, p2)
      )
    }
  )
  protest::test_results(as.numeric(dissim), epsilon, alpha, verbose, plot)
}
