#' @title Linear model L2 distance
#' 
#' @description L2 distance between a function f and the covariance matrix X.
#' While the distance allows for the use of any probability measure to be
#' applied to X, we used the empirical distribution for simplicity.
#' 
#' @param f Regression function that receives x
#' @param X Covariate matrix
#'
#' @return Distance between function f and the linear model specified by X
#' 
#' @examples
#' f <- function(x) sum(x)
#' set.seed(42)
#' n <- 100
#' X <- cbind(1, runif(n, 1, 2), runif(n, -1, 1))
#' lm.dist(f, X) #f is linear, so distance is close to 0
#' 
#' @export
lm.dist <-function(f, X){
  Ax <- t(X)%*%X
  fx <- base::apply(X, 1, function(x) f(x))
  Xfx <- t(X)%*%fx
  beta_hat <- base::solve(Ax, Xfx)
  sqrt(sum((X%*%beta_hat - fx)^2))
}

#' @title Linear model test
#' 
#' @description PROTEST version of a test that checks the adherence of linear
#' models to data. The empirical distribution of X is used as the probability
#' measure for simplicity. Results are a function of the significance level
#' alpha and the threshold epsilon.
#' 
#' @param f_list List of regression functions that receives x
#' @param X Covariate matrix
#' @param alpha Significance level, default is 0.05.
#' @param epsilon Threshold value
#' @param verbose Should the decision be returned as a message? Default if TRUE
#' @param plot Should the test return a plot of the decision for each
#' (alpha, epsilon)? Default is FALSE
#'
#' @return Test results of the linear model test
#' 
#' @examples
#' f_list <- list(f1 = function(x) sum(x), f2 = function(x) 2*sum(x))
#' set.seed(42)
#' n <- 100
#' X <- cbind(1, runif(n, 1, 2), runif(n, -1, 1))
#' lm.test(f_list, X, alpha = 0.05, epsilon = 0.1, verbose = TRUE, plot = TRUE)
#' 
#' @export
lm.test <- function(f_list, X, alpha = 0.05, epsilon, verbose = T, plot = F){
  dissim <- sapply(f_list, function(x) protest::lm.dist(x, X = X))
  protest::test.results(dissim, alpha, epsilon, verbose, plot)
}