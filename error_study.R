################################################################################
################################################################################
##################_____________________________________________#################
##################                                             #################
##################          Decision error in PROTEST          #################
##################_____________________________________________#################
##################                                             #################
################################################################################
################################################################################

library(kernlab)
library(MASS)
library(protest)

################################################################################
#####################          Simulation setting          #####################
################################################################################

# Custom functions for drawing from the GP
my_inv <- function(mat){
  mat_inv <- try(chol2inv(chol(mat)), silent = T)
  if(any(class(mat_inv) == 'try-error')) mat_inv <- ginv(mat)
  mat_inv
}

GP_pred <- function(Y, X, X_new, sigma2, mu_fun, ker_fun){
  mu_X <- mu_fun(X)
  mu_Xn <- mu_fun(X_new)
  ker_XX <- ker_fun(X,X)
  mat_inv <- my_inv(ker_XX + diag(sigma2, length(Y)))
  ker_XnX <- ker_fun(X_new,X)
  ker_XnXn <- ker_fun(X_new,X_new)
  f_mean <- mu_Xn + ker_XnX%*%mat_inv%*%(Y - mu_X)
  f_covmat <- ker_XnXn - ker_XnX%*%mat_inv%*%t(ker_XnX)
  list(f_mean = f_mean, f_covmat = f_covmat)
}

mu_fun <- function(X) rep(0, dim(X)[1]) #Mean function

my_ker <- rbfdot(sigma = 2)
ker_fun <- function(X1, X2)  kernelMatrix(my_ker, X1, X2) #Covariance function

# Simulation study
sample_size <- seq(50, 1500, length.out = 30) #Sample size
m <- 500 #Data draws
N <- 2000 #GP samples

sigma2 <- .5 #Variance

pol <- c(1, 1.25, 1.5, 1.75, 2) #Polynomial terms
pol_prag <- 1.4 #Polynomial term negligibly different from H0

ground_truth <- ifelse(pol_prag >= pol, "Do not reject", "Reject") #True state

eps <- lm_dist(f = function(x) x^pol_prag, X = cbind(runif(20000000, 0, 3)),
               add_intercept = T) #Getting epsilon

set.seed(42)
error_prob <- matrix(NA, nrow = length(pol), ncol = length(sample_size))
for(n in sample_size){ #Sample size
  for(i in 1:length(pol)){ #Polynomial term
    dec <- character(m)
    for(j in 1:m){ #Data draws
      # Make the data (Y,X)
      X <- cbind(runif(n, 0, 3))
      Y <- rnorm(n, mean = X^pol[i], sd = sqrt(sigma2))
      X_new <- cbind(seq(from = 0, to = 3, length.out = 61))
      
      # Model the data
      parms <- GP_pred(Y, X, X_new, sigma2, mu_fun, ker_fun)
      
      # List of functions
      f_list <- list()
      n_sims <- mvrnorm(n = N, mu = parms$f_mean, Sigma = parms$f_covmat)
      colnames(n_sims) <- X_new
      my_fun <- function(k){
        force(k)
        function(x) n_sims[k,as.character(x)]
      }
      for(k in 1:N) f_list[[k]] <- my_fun(k) 
      
      # Testing
      dec[j] <- lm_test(f_list, X_new, epsilon = eps, add_intercept = T,
                        alpha = 0.05, verbose = FALSE)$decision
    }
    #Mean error rate
    error_prob[i, which(sample_size == n)] <- mean(dec != ground_truth[i])
  }
}

colnames(error_prob) <- sample_size
row.names(error_prob) <- pol

# saveRDS(error_prob, file = "error_prob.RDS") #Saving the data


################################################################################
#####################          Simulation results          #####################
################################################################################
library(ggplot2)
library(scales)
library(tidyr)

error_prob <- readRDS("error_prob.RDS") #Loading the data

pol <- as.numeric(row.names(error_prob)) #Polynomials terms
sample_size <- as.numeric(colnames(error_prob)) #Sample sizes

error_prob <- as.data.frame(error_prob) #Turn do data frame
error_prob$pol <- pol #Add polynomial term as new column

#Pivoting
error_prob <- error_prob %>%
  pivot_longer(cols = as.character(sample_size), names_to = "sample_size",
               values_to='mean_error')
error_prob$sample_size <- as.numeric(error_prob$sample_size)

# Plotting results in paper
my_theme <- theme(legend.justification = c("right", "top"),
                  legend.box.just = "right",
                  legend.margin = margin(6, 6, 6, 6),
                  legend.box.background = element_rect(colour = "black", size = .5),
                  legend.text.align = 0,
                  legend.title = element_text(size = 20, face = "bold"),
                  legend.text = element_text(size=15),
                  axis.text = element_text(size = 15),
                  axis.title = element_text(size = 25))
ggplot2::ggplot(error_prob,
                aes(x=sample_size, y=mean_error, color=as.factor(pol))) +
  ggplot2::geom_line(size = 1) +
  ggplot2::geom_hline(yintercept = 0.05, linetype = "dashed") +
  ggplot2::scale_y_continuous(breaks = c(0.4, 0.3, 0.2, 0.1, 0.05, 0)) +
  ggplot2::scale_x_continuous(breaks = seq(0, 1500, by = 250)) +
  ggplot2::xlab("Sample size") + ggplot2::ylab("Prob(error)") +
  ggplot2::theme_bw() +
  ggplot2::guides(color = guide_legend(title = "X raised to:")) +
  my_theme
