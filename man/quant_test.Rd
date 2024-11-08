% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/quantile.R
\name{quant_test}
\alias{quant_test}
\title{Quantile test}
\usage{
quant_test(
  p_list,
  q_list,
  x0 = 0,
  p0 = 0.5,
  alpha = 0.05,
  epsilon,
  verbose = TRUE,
  plot = FALSE
)
}
\arguments{
\item{p_list}{List of distribution functions.}

\item{q_list}{List of quantile distributions associated to p_list.}

\item{x0}{Value from the null hypothesis, default is 0.}

\item{p0}{Probability from the null hypothesis, default is 0.5.}

\item{alpha}{Significance level, default is 0.05. When given two values, a
three-way test is performed.}

\item{epsilon}{Threshold value.}

\item{verbose}{Should the decision be returned as a message? Default if TRUE.}

\item{plot}{Should the test return a plot of the decision for each
(alpha, epsilon)? Default is FALSE}
}
\value{
Test results of the quantile test.
}
\description{
PROTEST version of the quantile test. For a given x0 and p0, it
checks if the data is such that F(x0) = p0, where F is the unknown true
distribution function.
}
\examples{
p_list <- list(
  pexp, function(x) pgamma(x, shape = 2),
  function(x) pexp(x, rate = 2)
)

q_list <- list(
  qexp, function(p) qgamma(p, shape = 2),
  function(p) qexp(p, rate = 2)
)

quant_test(p_list, q_list,
  alpha = 0.05, epsilon = 0.1, x0 = 1,
  verbose = TRUE, plot = TRUE
)

}
