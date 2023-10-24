#' @title L1 distance for the quantile test
#' 
#' @description L1 distance between a distribution function F and the function G
#' that best appoximates and obeys the condition G(x0) = p0.
#' 
#' @param pdist Distribution function that receives x.
#' @param qdist Quantile distribution associated to pdist.
#' @param x0 Value from the null hypothesis, default is 0.
#' @param p0 Probability from the null hypothesis, default is 0.5.
#'
#' @return Distance between distribution function F and the class of
#' distribution functions such that G(x0) = p0.
#' 
#' @examples
#' quant.dist(pdist = pexp, qdist = qexp, x0 = 1, p0 = 0.5)
#' 
#' @export
quant.dist <- function(pdist, qdist, x0 = 0, p0 = 0.5){
  F.quant <- qdist(p = p0)
  a = min(F.quant, x0)
  b = max(F.quant, x0)
  dissim <- stats::integrate(function(x) abs(p0 - pdist(x)), lower = a, upper = b)$value
  return(dissim)
}

#' @title Quantile test
#' 
#' @description PROTEST version of the quantile test. For a given x0 and p0, it
#' checks if the data is such that F(x0) = p0, where F is the unknown true
#' distribution function.
#' 
#' @param p_list List of distribution functions.
#' @param q_list List of quantile distributions associated to p_list.
#' @param x0 Value from the null hypothesis, default is 0.
#' @param p0 Probability from the null hypothesis, default is 0.5.
#' @param alpha Significance level, default is 0.05. When given two values, a
#' three-way test is performed.
#' @param epsilon Threshold value.
#' @param verbose Should the decision be returned as a message? Default if TRUE.
#' @param plot Should the test return a plot of the decision for each
#' (alpha, epsilon)? Default is FALSE
#'
#' @return Test results of the quantile test.
#' 
#' @examples
#' p_list <- list(pexp, function(x) pgamma(x, shape = 2),
#' function(x) pexp(x, rate = 2))
#' 
#' q_list <- list(qexp, function(p) qgamma(p, shape = 2),
#' function(p) qexp(p, rate = 2))
#' 
#' quant.test(p_list, q_list, alpha = 0.05, epsilon = 0.1, x0 = 1,
#' verbose = TRUE, plot = TRUE)
#' 
#' @export
quant.test <- function(p_list, q_list, x0 = 0, p0 = 0.5, alpha = 0.05, epsilon,
                       verbose = T, plot = F){
  n <- length(p_list)
  dissim <- sapply(1:length(p_list),
                   function(i) protest::quant.dist(p_list[[i]], q_list[[i]],
                                                     x0 = 0, p0 = 0.5))
  protest::test.results(dissim, alpha, epsilon, verbose, plot)
}
