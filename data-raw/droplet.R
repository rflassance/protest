library(GauPro)
library(ggpubr)
devtools::install()

elapsed_time <- seq(from = 0, to = 7, by = .5)
s <- c(4.7, 8.1, 11.5, 14.5, 17.4, 19.9, 22.4, 24.6, 26.6, 28.5, 30.1, 31.5,
       32.8, 33.8, 34.6)
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
                           alpha = 0.05, epsilon = epsilon/sqrt(length(x)),
                           verbose = T, plot = T)
d11 <- lm_res$dissim
g1 <- lm_res$plot
g1

lm_res <- protest::lm.test(f_list1, X = cbind(1, x), vars = 1,
                           alpha = 0.05, epsilon = epsilon/sqrt(length(x)),
                           verbose = T, plot = T)
d12 <- lm_res$dissim
r1 <- lm_res$plot
r1

f_list2 <- get_flist(gpm, x, n = 2000, seed = 42)

lm_res <- protest::lm.test(f_list2, X = cbind(1, x),
                           alpha = 0.05, epsilon = epsilon/sqrt(length(x)),
                           verbose = T, plot = T)
d21 <- lm_res$dissim
g2 <- lm_res$plot
g2

lm_res <- protest::lm.test(f_list2, X = cbind(1, x), vars = 1,
                           alpha = 0.05, epsilon = epsilon/sqrt(length(x)),
                           verbose = T, plot = T)
d22 <- lm_res$dissim
r2 <- lm_res$plot
r2


f_list3 <- get_flist(gpe, x, n = 2000, seed = 42)

lm_res <- protest::lm.test(f_list3, X = cbind(1, x),
                           alpha = 0.05, epsilon = epsilon/sqrt(length(x)),
                           verbose = T, plot = T)
d31 <- lm_res$dissim
g3 <- lm_res$plot
g3

lm_res <- protest::lm.test(f_list3, X = cbind(1, x), vars = 1,
                           alpha = 0.05, epsilon = epsilon/sqrt(length(x)),
                           verbose = T, plot = T)
d32 <- lm_res$dissim
r3 <- lm_res$plot
r3

ggarrange(g1, g2, g3)


ggarrange(r1, r2, r3)



y <- seq(from = 0, to = 1, length.out = 1000)
dec_df <- data.frame(x = numeric(), y = numeric(), kernel = character(),
                     test = character())
dis_list <- list(list(d11, d12), list(d21, d22), list(d31, d32))
ker <- c("Gaussian", "Matérn", "Exponential")
test <- c("First test", "Second test")
get_eps <- function(x, dissim) stats::quantile(dissim, probs = x)
for(i in 1:length(dis_list)){
  dissims <- dis_list[[i]]
  for(j in 1:length(dissims)){
    my_dissim <- dissims[[j]]
    x <- c(0, get_eps(y, my_dissim)[-1])
    dec_df <- rbind(dec_df, data.frame(x, y, kernel = ker[i], test = test[j]))
  }
}

my_theme <- theme(legend.justification = c("right", "top"),
              legend.box.just = "right",
              legend.margin = margin(6, 6, 6, 6),
              legend.box.background = element_rect(colour = "black", size = .5),
              legend.text.align = 0,
              legend.title = element_text(size = 20, face = "bold"),
              legend.text = element_text(size=15),
              axis.text = element_text(size = 15),
              axis.title = element_text(size = 25),
              strip.text = element_text(size=25, face = "bold"),
              strip.background = element_blank())

dec_df %>% dplyr::filter(test == "First test") %>%
  ggplot(aes(x = x, y = y, color = as.factor(kernel))) +
  geom_line() + theme_bw() + xlab(expression(epsilon)) +
  ylab(expression(Prob(Pg(H[0])*"|"*bold(X)==bold(x)))) +
  scale_color_discrete(name = expression(bold(Kernel))) +
  geom_vline(aes(xintercept = epsilon/sqrt(length(elapsed_time))),
             linetype = "dashed") +
  my_theme

dec_df %>% dplyr::filter(test == "Second test") %>%
  ggplot(aes(x = x, y = y, color = as.factor(kernel))) +
  geom_line() + theme_bw() + xlab(expression(epsilon)) +
  ylab(expression(Prob(Pg(H[0])*"|"*bold(X)==bold(x)))) +
  scale_color_discrete(name = expression(bold(Kernel))) +
  geom_vline(aes(xintercept = epsilon/sqrt(length(elapsed_time))),
             linetype = "dashed") +
  my_theme
