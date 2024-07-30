#' @title Goodness-of-fit test
#'
#' @description PROTEST version of the goodness-of-fit test. For a list of
#' distribution functions and a distribution function/parametric family being
#' tested, it tests the hypothesis that the true distribution of the data is
#' described in such function/family.
#'
#' @param p_list List of distribution functions.
#' @param H0_pdist Distribution function of the null hypothesis.
#' @param dissim Dissimilarity function to be applied, default is the L1
#' distance.
#' @param alpha Significance level, default is 0.05. When given two values, a
#' three-way test is performed.
#' @param epsilon Threshold value.
#' @param verbose Should the decision be returned as a message? Default if TRUE.
#' @param plot Should the test return a plot of the decision for each
#' (alpha, epsilon)? Default is FALSE
#' @param dissim_min Smallest value the dissimilarity can assume, required only
#' if plot = TRUE. Default is 0.
#'
#' @return Test results of the goodness-of-fit test.
#'
#' @examples
#' p_list <- list(
#'   pnorm, function(x) pt(x, df = 5), function(x) pt(x, df = 10),
#'   function(x) pt(x, df = 20), function(x) pt(x, df = 30)
#' )
#'
#' my_dissim <- function(pdist1, pdist2) {
#'   protest::Lp_dist(pdist1, pdist2, p = 2)
#' }
#'
#' gof_test(p_list,
#'   H0_pdist = pnorm, dissim = my_dissim, alpha = 0.05,
#'   epsilon = 0.01, verbose = TRUE, plot = TRUE
#' )
#'
#' @export
gof_test <- function(
    p_list, H0_pdist,
    dissim = protest::Lp_dist, alpha = 0.05,
    epsilon, verbose = TRUE, plot = FALSE, dissim_min = 0) {
  dissim_sam <- sapply(p_list, function(pdist) dissim(pdist, H0_pdist))
  protest::test_results(dissim_sam, alpha, epsilon, verbose, plot, dissim_min)
}
