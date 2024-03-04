library(ggplot2)
library(protest)

# Dirichlet process functions
rDP <- function(n = 1000, Y, rdist, alpha = 1, seed = NULL){
  N <- length(Y)
  alpha.post <- alpha + N
  
  # With prob alpha/alpha.post, draw from P0, else from sum delta_x
  set.seed(seed)
  p.delta <- rbinom(n, size = 1, prob = alpha/alpha.post)
  X <- numeric(n)
  X[p.delta == 1] <- rdist(sum(p.delta))
  X[p.delta == 0] <- sample(Y, size = sum(1-p.delta), replace = T)
  
  #Stick-breaking process
  theta <- rbeta(n, 1, alpha.post)
  theta_prod <- theta*c(1, sapply(2:n, function(x) prod(1 - theta[1:(x-1)])))
  
  dist.fun = function(y) sapply(y, function(x) sum(theta_prod*(X <= x)))
  
  points = sort(unique(X))
  
  return(list(dist.fun = dist.fun, points = points, probs = dist.fun(points)))
}

qDPdist <- function(p, dist.sample){
  dist.sample$points[which(dist.sample$probs >= p)[1]]
}

quant_dissimDP <- function(dist.samples, x0, p0 = 0.5){
  len_dist <- length(dist.samples)
  DP.quant <- lapply(dist.samples, function(x) qDPdist(p = p0, x))
  a = unlist(pmin(DP.quant, x0))
  b = unlist(pmax(DP.quant, x0))
  dissim <- numeric()
  int_fun <- function(model, a, b){
    integrate(function(x) abs(p0 - model(x)), lower = a, upper = b)
  }
  for(i in 1:len_dist){
    dissim[i] <- int_fun(dist.samples[[i]]$dist.fun, a[i], b[i])$value
  }
  return(dissim)
}



# Loading and cleaning data

paygap <- read.csv("data-raw/UNECE_PayGap.csv", header = T, na.strings = "..",
                   check.names = F)

paygap <- paygap$`2021` - paygap$`2020`
paygap <- paygap[!is.na(paygap)]

# DP model
B = 1000
rdist = function(n) rnorm(n, mean = 0, sd = sqrt(10))
set.seed(42)
dp.draws <- lapply(1:B, function(x) rDP(n = 1000, Y = paygap, rdist = rdist))
l1.dissim <- quant_dissimDP(dist.samples = dp.draws, x0 = 0, p0 = 0.25)

################################################################################
alpha = 0.05

quantile(l1.dissim, probs = alpha)

get_eps <- function(x) quantile(l1.dissim, probs = x)

tema <- theme(legend.position = c(.954, 0.645),
              legend.justification = c("right", "bottom"),
              legend.box.just = "right",
              legend.margin = margin(6, 6, 6, 6),
              legend.box.background = element_rect(colour = "black", size = .5),
              legend.text.align = 0,
              legend.title = element_text(size = 30, face = "bold"),
              legend.text = element_text(size=25),
              axis.text = element_text(size = 15),
              axis.title = element_text(size = 25))

ggplot(data = data.frame(`Pay gap` = paygap)) + theme_bw() +
  geom_function(fun = get_eps, col = 'red') + xlim(0,1) +
  geom_segment(aes(y = 0, yend = quantile(l1.dissim, probs = alpha), x = alpha,
                   xend = alpha),
               linetype = 'dashed') +
  geom_segment(aes(y = quantile(l1.dissim, probs = alpha),
                   yend = quantile(l1.dissim, probs = alpha),
                   x = 0, xend = alpha),
               linetype = 'dashed') +
  geom_point(aes(x = alpha, y = quantile(l1.dissim, probs = alpha)), size = 2) +
  coord_flip() +
  xlab(expression(alpha)) + ylab(expression(epsilon)) +
  tema

x <- c(0, seq(from = 0, to = 1, length.out = 100))
y <- c(0, get_eps(x)[-1])

data.frame(cbind(x,y)) %>%
  ggplot(aes(x = y, y = x)) + theme_bw() +
  # geom_segment(aes(y = 0, yend = quantile(l1.dissim, probs = alpha), x = alpha,
  #                  xend = alpha),
  #              linetype = 'dashed') +
  # geom_segment(aes(y = quantile(l1.dissim, probs = alpha),
  #                  yend = quantile(l1.dissim, probs = alpha),
  #                  x = 0, xend = alpha),
  #              linetype = 'dashed') +
  # geom_point(aes(x = alpha, y = quantile(l1.dissim, probs = alpha)), size = 2) +
  # coord_flip() +
  geom_line(size = 1) +
  geom_ribbon(ymin=-Inf, aes(ymax=x), fill='green', alpha=0.2) +
  xlab(expression(epsilon)) + ylab(expression(alpha)) +
  geom_ribbon(aes(ymin=x), ymax=Inf, fill='red', alpha=0.2) +
  scale_x_continuous(limits = c(0,max(y)), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0,1), expand = c(0, 0)) +
  tema +
  scale_colour_manual(name='', values=c("Reject" = "red", "Do not reject" = "green")) +
  scale_fill_manual(name = '',  values=c("Reject" = "red", "Do not reject" = "green"))


tema <- theme(#legend.position = c(.954, 0.645),
  legend.justification = c("right", "top"),
  legend.box.just = "right",
  legend.margin = margin(6, 6, 6, 6),
  legend.box.background = element_rect(colour = "black", size = .5),
  legend.text.align = 0,
  legend.title = element_text(size = 20, face = "bold"),
  legend.text = element_text(size=15),
  axis.text = element_text(size = 15),
  axis.title = element_text(size = 25))

data.frame(cbind(x,y)) %>%
  ggplot(aes(x = y, y = x)) + theme_bw() +
  geom_line(size = 1) +
  geom_ribbon(ymin=-Inf, aes(ymax=x, fill='notrej'), alpha=0.2) +
  geom_ribbon(aes(ymin=x, ymax=Inf, fill='rej'), alpha=0.2) +
  geom_segment(aes(x = 0, xend = quantile(l1.dissim, probs = alpha), y = alpha,
                   yend = alpha),
               linetype = 'dashed') +
  geom_segment(aes(x = quantile(l1.dissim, probs = alpha),
                   xend = quantile(l1.dissim, probs = alpha),
                   y = 0, yend = alpha),
               linetype = 'dashed') +
  geom_point(aes(y = alpha, x = quantile(l1.dissim, probs = alpha)), size = 3) +
  xlab(expression(epsilon)) +
  ylab(expression(alpha)) +
  scale_x_continuous(limits = c(0,max(y)), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0,1), expand = c(0, 0)) +
  tema +
  scale_fill_manual(name = 'Decision', labels = c("Reject", "Do not reject"),
                    values=c("rej" = "red", "notrej" = "green"))
