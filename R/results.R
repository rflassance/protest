#' @title PROTEST results
#'
#' @description Test results based on the dissimilarity, significance alpha and
#' threshold epsilon. The procedure can also return a plot that shows the
#' decision based on each combination of (alpha, epsilon).
#'
#' @param dissim Vector of dissimilarities.
#' @param alpha Significance level, default is 0.05.
#' @param epsilon Threshold.
#' @param verbose Should the decision be returned as a message? Default is TRUE.
#' @param plot Should the test return a plot of the decision for each
#' (alpha, epsilon)? Default is FALSE.
#' @param dissim_min Smallest value possible for the dissimilarity function,
#' default is 0.
#' @param three_way Should the results provided be of a three-way test? Default
#' is FALSE. If two values for alpha are provided, this option is ignored.
#'
#' @return Test results (decision, dissimilarity vector, significance level,
#' threshold and plot). If plotting on the three-way case, the plot assumes that
#' alpha2 = 1 - alpha1.
#'
#' @examples
#' dissim <- seq(.5, 1, length.out = 100)
#' test_results(
#'   dissim = dissim,
#'   epsilon = .8,
#'   verbose = TRUE, plot = TRUE, dissim_min = 0.5
#' )
#'
#' @export

test_results <- function(
    dissim,
    epsilon,
    alpha = 0.05,
    verbose = TRUE,
    plot = FALSE,
    dissim_min = 0,
    three_way = FALSE) {
  prag_prob <- sapply(epsilon, function(x) mean(dissim <= x))
  if (three_way || length(alpha) == 2) {
    alpha <- sort(alpha)
    if (length(alpha) == 1) {
      if (alpha > .5) stop("Please specify a value for
      alpha below or equal to 0.5 or set alpha2 yourself.")
      alpha <- c(alpha, 1 - alpha)
    }
    decision <- ifelse(prag_prob <= alpha[1], "Reject",
      ifelse(prag_prob <= alpha[2], "Remain agnostic", "Accept")
    )
    if (verbose) {
      cat("PROTEST conclusion: ", ifelse(decision != "Remain agnostic",
        decision, "Remain agnostic about"
      ),
      " H0.\n",
      sep = ""
      )
      message("Prob(Pg(H0)) = ", prag_prob, ".")
    }
    if (plot) {
      x <- seq(from = 0, to = .5, length.out = 1000)
      get_eps1 <- function(x) stats::quantile(dissim, probs = x)
      get_eps2 <- function(x) stats::quantile(dissim, probs = 1 - x)
      y <- c(dissim_min, get_eps1(x[-1]))
      z <- get_eps2(x)
      g <- ggplot2::ggplot(
        data = data.frame(cbind(x, y)),
        ggplot2::aes(x = y, y = x)
      ) +
        ggplot2::theme_bw() +
        ggplot2::geom_line() +
        ggplot2::geom_line(ggplot2::aes(x = z, y = x)) +
        ggplot2::geom_ribbon(ggplot2::aes(xmin = z, xmax = Inf, fill = "acc"),
          alpha = 0.2
        ) +
        ggplot2::geom_ribbon(ggplot2::aes(xmin = y, xmax = z, fill = "agn"),
          alpha = 0.2
        ) +
        ggplot2::geom_ribbon(ggplot2::aes(xmin = -Inf, xmax = y, fill = "rej"),
          alpha = 0.2
        ) +
        ggplot2::xlab(expression(epsilon)) +
        ggplot2::ylab(expression(Prob(Pg(H[0]) * "|" * bold(X) == bold(x)))) +
        ggplot2::scale_x_continuous(
          limits = c(dissim_min, max(z)),
          expand = c(0, 0)
        ) +
        ggplot2::scale_y_continuous(limits = c(0, .5), expand = c(0, 0)) +
        ggplot2::theme(
          legend.justification = c("right", "top"),
          legend.box.just = "right",
          legend.title = ggplot2::element_text(face = "bold")
        ) +
        ggplot2::scale_fill_manual(
          name = "Decision",
          labels = c(
            "rej" = "Reject",
            "agn" = "Remain agnostic",
            "acc" = "Accept"
          ),
          values = c(
            "rej" = "red", "agn" = "yellow",
            "acc" = "green"
          )
        )
    } else {
      g <- NULL
    }
  } else {
    decision <- ifelse(prag_prob <= alpha, "Reject", "Do not reject")
    if (verbose) {
      cat("PROTEST conclusion: ", decision, " H0.\n", sep = "")
      message("Prob(Pg(H0)|x) = ", prag_prob, ".")
    }
    if (plot) {
      x <- seq(from = 0, to = 1, length.out = 1000)
      get_eps <- function(x) stats::quantile(dissim, probs = x)
      y <- c(dissim_min, get_eps(x[-1]))
      g <- ggplot2::ggplot(
        data = data.frame(cbind(x, y)),
        ggplot2::aes(x = y, y = x)
      ) +
        ggplot2::theme_bw() +
        ggplot2::geom_line() +
        ggplot2::geom_ribbon(ggplot2::aes(ymin = -Inf, ymax = x, fill = "nrej"),
          alpha = 0.2
        ) +
        ggplot2::geom_ribbon(ggplot2::aes(ymin = x, ymax = Inf, fill = "rej"),
          alpha = 0.2
        ) +
        ggplot2::xlab(expression(epsilon)) +
        ggplot2::ylab(expression(alpha)) +
        ggplot2::scale_x_continuous(
          limits = c(dissim_min, max(y)),
          expand = c(0, 0)
        ) +
        ggplot2::scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
        ggplot2::theme(
          legend.justification = c("right", "top"),
          legend.box.just = "right",
          legend.title = ggplot2::element_text(face = "bold")
        ) +
        ggplot2::scale_fill_manual(
          name = "Decision",
          labels = c(
            "rej" = "Reject",
            "nrej" = "Do not reject"
          ),
          values = c("rej" = "red", "nrej" = "green")
        )
    } else {
      g <- NULL
    }
  }
  list(
    decision = decision, dissim = dissim, alpha = alpha, epsilon = epsilon,
    plot = g
  )
}
