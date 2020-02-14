library(extraDistr)

obs <- c(17, 24, 22, 17, 24, 16)
obs2 <- c(22, 20, 30, 16, 16, 16)
value <- 1:6
names(obs) <- c("1", "2", "3", "4", "5", "6")
names(obs2) <- names(obs)
n1 <- sum(obs)
n2 <- sum(obs2)
ysum <- sum(obs * value)
#categorical
#jeffreys prior is dirichlet(a = 1/2)

#theta_1, ..., theta_6 ~ Dirichlet(a_1, ..., a_6)
#c_1, ..., c_6 ~ Multinomial(theta_1, ..., theta_6)
#a'_k = a_k + ysum

j_prior_params <- rep(1 / 2, 6)

post_params <- j_prior_params + obs
post_params_new <- j_prior_params + obs + obs2

post_pred <- post_params / sum(post_params)
post_pred_new <- post_params_new / sum(post_params_new)

alpha = 0.05
for (i in post_params){
  print(qbeta(c(0.025, 0.975), i, sum(post_params) - i))
}
for (i in post_params_new){
  print(qbeta(c(0.025, 0.975), i, sum(post_params_new) - i))
}

# n <- 1e6
# rope_r <- 0.05#1/30
# A <- rdirichlet(n, post_params)
# B <- as.matrix(A)
# f <- function(x) all((x > 1/6-rope_r)*(x < 1/6+rope_r))
# mean(apply(B, MARGIN = 1, FUN = f))

# nh_p <- rep(1/6, 6)
# n <- 1e2
# size <- sum(obs)
# nh_lik <- t(rmultinom(n, size, nh_p))

#pearson chi-squared

expected <- n1 * rep(1/6, 6)
observed <- obs
sq_diff <- (observed - expected)^2
pear_stat <- sq_diff / expected
chi_stat <- sum(pear_stat)
chi_dof <- length(obs) - 1
chi_p <- pchisq(chi_stat, chi_dof)
chi_p

nh <- rep(1/6, 6)
mle <- obs/n1
lr_stat <- -2 * sum(obs * log(nh/mle))
lr_chi_stat <- pchisq(lr_stat, chi_dof)
lr_chi_stat

n <- 1e5

C <- rdirmnom(n, 60, post_params)
D <- as.matrix(C)
ts <- D[,5] > D[,6]
sum(ts)/n
