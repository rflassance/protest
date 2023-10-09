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
  fx <- (apply(X, 1, function(x) f(x)))
  Xfx <- t(X)%*%fx
  beta_hat <- solve(Ax, Xfx)
  sqrt(sum((X%*%beta_hat - fx)^2))
}

lm.test <- function(f_list, X, alpha = 0.05, epsilon, verbose = T, plot = F){
  f_len <- length(f_list)
  dissim <- lapply(f_list, lm.dist, X = X)
  test.results(dissim, alpha, epsilon, verbose, plot)
}