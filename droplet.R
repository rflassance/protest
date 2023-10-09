library(GauPro)
library(ggpubr)
devtools::install()

elapsed_time <- seq(from = 0, to = 7, by = .5)
s <- c(4.7, 8.1, 11.5, 14.5, 17.4, 19.9, 22.4, 24.6, 26.6, 28.5, 30.1, 31.5, 32.8, 33.8, 34.6)
delta_s <- c(NA, diff(s))
Ks <- c(8.446, 8.569, 8.694)
drop_rad <- sqrt(2*delta_s*Ks[1])

plot(elapsed_time, drop_rad)

x <- elapsed_time
y <- drop_rad
epsilon <- 0.621835

gp <- GauPro(x[-1], y[-1])
kern <- Matern52$new(0)
gpm <- GauPro_kernel_model$new(matrix(x[-1], ncol=1), y[-1], kernel=kern)
kern.exp <- Exponential$new(0)
gpe <- GauPro_kernel_model$new(matrix(x[-1], ncol=1), y[-1], kernel=kern.exp)

n <- 100

get_flist <- function(gp, X, n, seed = NULL){
  f_list <- list()
  set.seed(seed)
  gp_sam <- gp$sample(x, n)
  my_fun <- function(i){
    force(i)
    function(x) gp_sam[i, which(X == tail(x, 1))]
  }
  for(i in 1:n){
    f_list[[i]] <- my_fun(i)
  }
  f_list
}

f_list1 <- get_flist(gp, x, n = 2000, seed = 42)

lm_res <- protest::lm.test(f_list1, X = cbind(1, x),
                           alpha = 0.05, epsilon = epsilon,
                           verbose = T, plot = T)
g1 <- lm_res$plot
g1

lm_res <- protest::lm.test(f_list1, X = cbind(1, x), vars = 1,
                           alpha = 0.05, epsilon = epsilon,
                           verbose = T, plot = T)
r1 <- lm_res$plot
r1

f_list2 <- get_flist(gpm, x, n = 2000, seed = 42)

lm_res <- protest::lm.test(f_list2, X = cbind(1, x),
                           alpha = 0.05, epsilon = epsilon,
                           verbose = T, plot = T)
g2 <- lm_res$plot
g2

lm_res <- protest::lm.test(f_list2, X = cbind(1, x), vars = 1,
                           alpha = 0.05, epsilon = epsilon,
                           verbose = T, plot = T)
r2 <- lm_res$plot
r2


f_list3 <- get_flist(gpe, x, n = 2000, seed = 42)

lm_res <- protest::lm.test(f_list3, X = cbind(1, x),
                           alpha = 0.05, epsilon = epsilon,
                           verbose = T, plot = T)
g3 <- lm_res$plot
g3

lm_res <- protest::lm.test(f_list3, X = cbind(1, x), vars = 1,
                           alpha = 0.05, epsilon = epsilon,
                           verbose = T, plot = T)
r3 <- lm_res$plot
r3

ggarrange(g1, g2, g3)


ggarrange(r1, r2, r3)
