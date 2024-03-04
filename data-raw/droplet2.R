library(GauPro)
library(ggpubr)
library(protest)
library(dplyr)

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

mod <- lm(y ~ x)

preds <- predict(mod, data.frame(x))

pred_lower <- preds - epsilon
pred_upper <- preds + epsilon

gp <- GauPro(x[-1], y[-1])

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

f_list <- get_flist(gp, x, n = 20, seed = 42)

sam_gp <- data.frame(pred = numeric(), x = numeric(), ind = integer())
for(i in 1:length(f_list)){
  sam_gp <- bind_rows(sam_gp, data.frame(pred = f_list[[i]](t(x)), x = x,
                                         ind = i))
}

lm_res1 <- data.frame(preds = preds, y = y, x = x, prag_lower = pred_lower,
                      prag_upper = pred_upper)
lm_res2 <- data.frame(preds = mean(y, na.rm = T), y = y, x = x,
                      prag_lower = mean(y, na.rm = T) - epsilon,
                      prag_upper = mean(y, na.rm = T) + epsilon)

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

ggplot(data = lm_res1, aes(y = y, x = x)) +
  geom_ribbon(data = lm_res2, aes(ymin = prag_lower, ymax = prag_upper),
              fill = 'azure3', color = "royalblue4", alpha=0.5) +
  geom_ribbon(aes(ymin = prag_lower, ymax = prag_upper), fill = 'azure3',
              color = "firebrick4", alpha=0.5) +
  geom_line(data = lm_res2, aes(y = preds), linetype = 'dashed', size = 1.5,
            col = "royalblue4") +
  geom_line(aes(y = preds), linetype = 'dashed', size = 1.5,
            col = "firebrick4") +
  geom_point(aes(y = y), size = 2.5) +
  theme_bw() + xlab("Elapsed time") + ylab("Droplet radius") + my_theme

ggplot(data = lm_res1, aes(y = y, x = x)) +
  geom_ribbon(data = lm_res2, aes(ymin = prag_lower, ymax = prag_upper),
              fill = 'azure3', color = "royalblue4", alpha=0.5) +
  geom_ribbon(aes(ymin = prag_lower, ymax = prag_upper), fill = 'azure3',
              color = "firebrick4", alpha=0.5) +
  geom_line(data = sam_gp, aes(y = pred, x = x, group = as.factor(ind)),
            col = "chartreuse4", size = .5, alpha = .8, show.legend = FALSE) +
  geom_point(data = sam_gp, aes(y = pred, x = x, group = as.factor(ind)),
             col = "chartreuse4", size = 1, alpha = .8, show.legend = FALSE) +
  geom_line(data = lm_res2, aes(y = preds), linetype = 'dashed', size = 1.5,
            col = "royalblue4") +
  geom_line(aes(y = preds), linetype = 'dashed', size = 1.5,
            col = "firebrick4") +
  theme_bw() + xlab("Elapsed time") + ylab("Droplet radius") + my_theme


################################################################################

library(GauPro)
library(ggpubr)
library(protest)
library(dplyr)
library(reshape2)

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

mod <- lm(y ~ x)

gp <- GauPro(x[-1], y[-1])

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

f_list <- get_flist(gp, x, n = 10000, seed = 42)

GP_mat <- matrix(NA, nrow = length(f_list), ncol = length(x))
for(i in 1:length(f_list)){
  GP_mat[i,] <- f_list[[i]](t(x))
}

min_ind <- which.min(GP_mat[,1])
max_ind <- which.max(GP_mat[,1])

lm.beta <-function(f, X, g = function(x) x, vars = 'all', add_intercept = F){
  if(vars == 'all'){
    XX <- X
  } else XX <- X[,vars, drop = F]
  if(add_intercept) XX <- base::cbind(1,XX)
  Ax <- t(XX)%*%XX
  fx <- base::apply(X, 1, function(x) g(f(x)))
  Xfx <- t(XX)%*%fx
  if(dim(XX)[2] > 1){
    beta_hat <- base::solve(Ax, Xfx)
  } else beta_hat <- Xfx/Ax
  XX%*%beta_hat
}

c(lm.dist(f_list[[min_ind]],t(t(x)),add_intercept = T),
lm.dist(f_list[[max_ind]],t(t(x)),add_intercept = T))
c(lm.dist(f_list[[min_ind]],t(t(x)),add_intercept = T,
        vars = 0),
lm.dist(f_list[[max_ind]],t(t(x)),add_intercept = T,
        vars = 0))

df <- data.frame(time = x,
                 vals = c(GP_mat[min_ind,],
                          lm.beta(f_list[[min_ind]],t(t(x)),add_intercept = T),
                          lm.beta(f_list[[min_ind]],t(t(x)),add_intercept = T,
                                  vars = 0),
                          GP_mat[max_ind,],
                          lm.beta(f_list[[max_ind]],t(t(x)),add_intercept = T),
                          lm.beta(f_list[[max_ind]],t(t(x)),add_intercept = T,
                                  vars = 0)),
                ref = c(rep(1, 3*length(x)), rep(2, 3*length(x))),
                fun = rep(c(rep('GP', length(x)), rep('H01', length(x)),
                            rep('H02', length(x))), 2)
                )

my_theme <- theme(legend.justification = c("right", "top"),
                  legend.box.just = "right",
                  legend.margin = margin(6, 6, 6, 6),
                  legend.box.background = element_rect(colour = "black", size = .5),
                  legend.text.align = 0,
                  legend.title = element_text(size = 20, face = "bold"),
                  legend.text = element_text(size=18),
                  axis.text = element_text(size = 15),
                  axis.title = element_text(size = 25),
                  strip.text = element_text(size=25, face = "bold"),
                  strip.background = element_blank(),
                  strip.text.x = element_blank())

ggplot(df, aes(y = vals, x = time, color = as.factor(ref), fill = fun)) +
  geom_line(aes(linetype = fun), size = 1.5) + guides(color = 'none') +
  scale_linetype_discrete("Functions",
                          labels = c("GP" = expression(italic(a)(italic(t))),
                                     "H01" = expression(argmin[italic(f) %in% H[0]^1](d(italic(f), italic(a)))),
                                     "H02" = expression(argmin[italic(f) %in% H[0]^2](d(italic(f), italic(a)))))) +
  labs(x = "Time", y = "Droplet radius") +
  guides(linetype = guide_legend(override.aes = list(size = .7))) +
  theme_bw() + my_theme

ggplot(df, aes(y = vals, x = time, color = as.factor(ref))) +
  geom_line(aes(linetype = fun), size = 1.5) + guides(color = 'none') +
  scale_linetype_discrete("Functions",
                          labels = c("GP" = expression("Candidate for"~italic(a)(italic(t))),
                                     "H01" = expression(argmin[italic(f) %in% H[0]^1](d(italic(f), italic(a)))),
                                     "H02" = expression(argmin[italic(f) %in% H[0]^2](d(italic(f), italic(a)))))) +
  labs(x = "Time", y = "Droplet radius") +
  guides(linetype = guide_legend(override.aes = list(size = .7))) +
  theme_bw() + my_theme + facet_wrap(.~as.factor(ref))



lm_res <- lm.test(f_list, X = cbind(1,x),
                  alpha = 0.05, epsilon = epsilon/sqrt(length(x)),
                  verbose = T, plot = T)
d1 <- lm_res$dissim
lm_res <- lm.test(f_list, X = cbind(1,x), vars = 1,
                  alpha = 0.05, epsilon = epsilon/sqrt(length(x)),
                  verbose = T, plot = T)
d2 <- lm_res$dissim

y <- seq(from = 0, to = 1, length.out = 1000)
dec_df <- data.frame(x = numeric(), y = numeric(),
                     test = character())
dis_list <- list(d1, d2)
test <- c("H01", "H02")
get_eps <- function(x, dissim) stats::quantile(dissim, probs = x)
for(i in 1:length(dis_list)){
  my_dissim <- dis_list[[i]]
  x <- c(0, get_eps(y, my_dissim)[-1])
  dec_df <- rbind(dec_df, data.frame(x, y, test = test[i]))
}


dec_df %>%
  ggplot(aes(x = x, y = y, color = as.factor(test))) +
  geom_line(size = 1.5) + theme_bw() + xlab(expression(epsilon)) +
  ylab(expression(Prob(Pg(H[0])*"|"*bold(X)==bold(x)))) +
  scale_color_discrete(name = expression(bold(Hypothesis)),
                       labels = c("H01" = expression(H[0]^1*":"~ a(t) *"="* beta[0] ~"+"~ beta[1]*t),
                                  "H02" = expression(H[0]^2*":"~ a(t) *"="* beta[0]))) +
  geom_vline(aes(xintercept = epsilon/sqrt(length(elapsed_time))),
             linetype = "dashed") +
  geom_hline(aes(yintercept = 0.05), linetype = "dashed") +
  geom_point(x = epsilon/sqrt(length(elapsed_time)), y = 0.05, color = 'black', size = 3) +
  labs(tag = expression(alpha)) +
  theme(plot.tag.position = c(0.07, 0.23),
        plot.tag = element_text(size = 25, face = "bold"))+
  my_theme



mH1 <- cbind(GP_mat[min_ind,],
             lm.beta(f_list[[min_ind]],t(t(x)),add_intercept = T))
mH1 <- cbind(mH1, t(apply(mH1, 1, function(x) c(min(x), max(x)))))

mH2 <- cbind(GP_mat[min_ind,],
             lm.beta(f_list[[min_ind]],t(t(x)),add_intercept = T, vars = 0))
mH2 <- cbind(mH2, t(apply(mH2, 1, function(x) c(min(x), max(x)))))

MH1 <- cbind(GP_mat[max_ind,],
             lm.beta(f_list[[max_ind]],t(t(x)),add_intercept = T))
MH1 <- cbind(MH1, t(apply(MH1, 1, function(x) c(min(x), max(x)))))

MH2 <- cbind(GP_mat[max_ind,],
             lm.beta(f_list[[max_ind]],t(t(x)),add_intercept = T, vars = 0))
MH2 <- cbind(MH2, t(apply(MH2, 1, function(x) c(min(x), max(x)))))

df2 <- data.frame(rbind(mH1, mH2, MH1, MH2))
names(df2) <- c("vals", "argmin", "lower", "upper")
df2 <- df2[,c(2,1,3:dim(df2)[2])]

df2$hip <- factor(c(rep('H01', length(x)), rep('H02', length(x)),
                    rep('H01', length(x)), rep('H02', length(x))))
df2$fun <- factor(c(rep('f1', 2*length(x)), rep('f2', 2*length(x))))

df2$time <- x

df2$diff <- df2$upper - df2$lower

df2$diff2 <- (df2$upper - df2$lower)^2/length(x)

my_theme2 <- theme(legend.justification = c("right", "top"),
                  legend.box.just = "right",
                  legend.margin = margin(6, 6, 6, 6),
                  legend.box.background = element_rect(colour = "NA", size = .5),
                  legend.text.align = 0,
                  legend.title = element_text(size = 20, face = "bold"),
                  legend.text = element_text(size=18),
                  axis.text = element_text(size = 18),
                  axis.title = element_text(size = 25),
                  strip.text = element_text(size=20, face = "bold"),
                  strip.background = element_blank())


library(ggridges)

df2 %>%
  mutate(hip = recode(hip,
                      "H01" = "H[0]^1: a(t) == beta[0]+beta[1]*t",
                      "H02" = "H[0]^2: a(t) == beta[0]"),
         fun = recode(fun,
                      "f1" = "Candidate~1", "f2" = "Candidate~2")) %>%
  melt(id = c("hip", "fun", "time", "lower", "upper", "diff", "diff2")) %>%
  ggplot(aes(x = time, y = value, color = variable)) +
  geom_ridgeline_gradient(aes(y = lower, height = diff, fill = diff2, color = NA)) +
  geom_line(size = 1.25) +
  facet_grid(vars(hip), vars(fun),
             labeller = label_parsed, scales = "free_y") +
  theme_bw() + labs(x = "Time", y = "Droplet radius") + my_theme2 +
  scale_color_discrete("Functions",
                       guide = guide_legend(reverse = TRUE),
                       labels = c("vals" = expression("Candidate for"~italic(a)(italic(t))),
                                  "argmin" = expression(argmin[italic(f) %in% H[0]](d(italic(f), italic(a)))))) +
  scale_fill_continuous("(Squared error)/|T|", low = "darkgray", high = "azure2")
