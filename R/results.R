#' @title PROTEST results
#' 
#' @description Test results based on the dissimilarity, significance alpha and
#' threshold epsilon. The procedure can also return a plot that shows the
#' decision based on each combination of (alpha, epsilon).
#' 
#' @param dissim Vector of dissimilarities
#' @param alpha Significance level, default is 0.05
#' @param epsilon Threshold
#' @param verbose Should the decision be returned as a message? Default if TRUE
#' @param plot Should the test return a plot of the decision for each
#' (alpha, epsilon)? Default is FALSE
#' @param dissim_min Smallest value possible for the dissimilarity function,
#' default is 0.
#' 
#' @return Test results (decision, dissimilarity vector, significance level,
#' threshold and plot).
#' 
#' @examples
#' dissim <- seq(.5, 1, length.out = 100)
#' test.results(dissim = dissim, epsilon = .5, verbose = TRUE, plot = TRUE, dissim_min = 0.5)
#' 
#' @export

test.results <- function(dissim, alpha = 0.05, epsilon, verbose = T, plot = F,
                         dissim_min = 0){
  prag.prob <- sapply(epsilon, function(x) mean(dissim <= x))
  decision <- ifelse(prag.prob <= alpha, 'Reject', 'Do not reject')
  if(verbose){
    cat("PROTEST conclusion: ", decision, " H0.\n")
    message("Probability of Pg(H0) = ", prag.prob, ".")
  }
  if(plot){
    x <- seq(from = 0, to = 1, length.out = 1000)
    get_eps <- function(x) stats::quantile(dissim, probs = x)
    y <- c(dissim_min, get_eps(x[-1]))
    g <- ggplot2::ggplot(data = data.frame(cbind(x,y)),
                         ggplot2::aes(x = y, y = x)) + ggplot2::theme_bw() +
      ggplot2::geom_line() +
      ggplot2::geom_ribbon(ggplot2::aes(ymin=-Inf, ymax=x, fill='nrej'), alpha=0.2) +
      ggplot2::geom_ribbon(ggplot2::aes(ymin=x, ymax=Inf, fill='rej'), alpha=0.2) +
      ggplot2::xlab(expression(epsilon)) + ggplot2::ylab(expression(alpha)) +
      ggplot2::scale_x_continuous(limits = c(dissim_min, max(y)), expand = c(0, 0)) +
      ggplot2::scale_y_continuous(limits = c(0,1), expand = c(0, 0)) +
      ggplot2::theme(legend.justification = c("right", "top"),
                     legend.box.just = "right",
                     legend.title = ggplot2::element_text(face = "bold")) +
      ggplot2::scale_fill_manual(name = 'Decision',
                                 labels = c("Reject", "Do not reject"),
                                 values = c("rej" = "red", "nrej" = "green"))
    
  } else g <- NULL
  return(list(decision = decision, dissim = dissim,
              alpha = alpha, epsilon = epsilon, plot = g))
}