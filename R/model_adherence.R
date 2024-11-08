#' @title Linear model L2 distance
#'
#' @description L2 distance between a function f and the linear model is closest
#' to it on the domain of covariance matrix X. While the distance allows for the
#' use of any probability measure to be applied to X, we use the empirical
#' distribution for simplicity.
#'
#' @param f Regression function that receives x.
#' @param X Covariate matrix.
#' @param g Link function. Default is the identity function.
#' @param vars Either 'all' or a vector of the columns of X to be added to the
#' test. Default is 'all'.
#' @param add_intercept Should the intercept be added to X? Default is FALSE.
#'
#' @return Distance between function f and the linear model specified by X.
#'
#' @examples
#' f <- function(x) sum(x)
#' set.seed(42)
#' n <- 100
#' X <- cbind(1, runif(n, 1, 2), runif(n, -1, 1))
#' lm_dist(f, X) # f is linear, so distance is close to 0
#'
#' @export
lm_dist <- function(
    f,
    X, # nolint
    g = function(x) x, vars = "all", add_intercept = FALSE) {
  if (vars == "all") {
    XX <- X # nolint
  } else {
    XX <- X[, vars, drop = FALSE]
  }
  if (add_intercept) XX <- base::cbind(1, XX)
  Ax <- t(XX) %*% XX
  fx <- base::apply(X, 1, function(x) g(f(x)))
  Xfx <- t(XX) %*% fx
  if (dim(XX)[2] > 1) {
    beta_hat <- base::solve(Ax, Xfx)
  } else {
    beta_hat <- Xfx / Ax
  }
  sqrt(sum((XX %*% beta_hat - fx)^2) / dim(XX)[1])
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
#' @param g Link function. Default is the identity function.
#' @param vars Either 'all' or a vector of the columns of X to be added to the
#' test. Default is 'all'.
#' @param add_intercept Should the intercept be added to X? Default is FALSE.
#' @param alpha Significance level, default is 0.05. When given two values, a
#' three-way test is performed.
#' @param epsilon Threshold value.
#' @param verbose Should the decision be returned as a message? Default if TRUE.
#' @param plot Should the test return a plot of the decision for each
#' (alpha, epsilon)? Default is FALSE
#'
#' @return Test results of the linear model test
#'
#' @examples
#' f_list <- list(f1 = function(x) sum(x), f2 = function(x) 2 * sum(x))
#' set.seed(42)
#' n <- 100
#' X <- cbind(1, runif(n, 1, 2), runif(n, -1, 1))
#' lm_test(f_list, X, alpha = 0.05, epsilon = 0.1, verbose = TRUE, plot = TRUE)
#'
#' @export
lm_test <- function(f_list, X, epsilon, g = function(x) x, vars = "all",
                    add_intercept = FALSE, alpha = 0.05, verbose = TRUE,
                    plot = FALSE) {
  dissim <- sapply(
    f_list,
    function(x) {
      protest::lm_dist(x,
        X = X, g = g, vars = vars,
        add_intercept = add_intercept
      )
    }
  )
  protest::test_results(dissim, epsilon, alpha, verbose, plot)
}
