library(coda)
library(rstan)
library(rstantools)
library(bayesplot)
library(ggplot2)

library(dplyr)
library(data.table)
library(modeest)

rstan_options(auto_write = TRUE)
options(datatable.fread.datatable = FALSE)
options(mc.cores = parallel::detectCores() - 1)
devAskNewPage(ask = FALSE)
par(ask = F)
#options(mc.cores = 1)

abalone <- fread('data/abalone.csv')
names(abalone) <- make.names(names(abalone))
abalone$Age <- abalone$Rings + 1.5
abalone <- subset(abalone, select = -c(Rings))

cmeans <-
  colMeans(Filter(is.numeric, subset(abalone, select = -c(Age))))
amean <- mean(abalone$Age)

scale2 <- function(x, na.rm = FALSE)
  (x - mean(x, na.rm = na.rm))
abalone <- mutate_if(abalone, is.numeric, scale2, na.rm = TRUE)
abalone$Sex <- as.factor(abalone$Sex)

X <-
  model.matrix(
    Age ~ Sex + Length + Diameter + Height + Whole.weight + Shucked.weight + Viscera.weight + Shell.weight,
    data = abalone
  )

xlevs <-
  lapply(abalone[, sapply(abalone, is.factor), drop = F], function(j) {
    levels(j)
  })

new_obs <-
  data.frame(
    Sex = 'M',
    Length = 0.515,
    Diameter = 0.400,
    Height = 0.133,
    Whole.weight = 0.531,
    Shucked.weight = 0.231,
    Viscera.weight = 0.122,
    Shell.weight = 0.168
  )

dfcmean <- as.data.frame(t(cmeans))

new_obs[-1] <- new_obs[-1] - dfcmean
mm_new <- model.matrix( ~ ., data = new_obs, xlev = xlevs)

f_lm <- lm(data = abalone, Age ~ .)
summary(f_lm)
pred_f_lm <- predict(f_lm, new_obs)
pred_f_lm
par(mfrow = c(2, 2))
plot(f_lm)
par(mfrow = c(1, 1))
plot(fitted(f_lm), abalone$Age, asp = 1)


b_lm_data <-
  list(
    N = nrow(X),
    K = ncol(X),
    X = X,
    y = abalone$Age,
    N_new = 1,
    x_new = as.array(mm_new)
  )

b_lm <- stan_model(file = 'a1q3.stan')

b_lm_fit <- sampling(
  b_lm,
  data = b_lm_data,
  chains = 7,
  control = list(adapt_delta = 0.8),
  iter = 8000
)
plot(b_lm_fit)
b_lm_fit

smps <- extract(b_lm_fit)
tp <- traceplot(b_lm_fit, pars = c("beta", "sigma"))
tp

b_lm_coda <- As.mcmc.list(b_lm_fit)
gelman.plot(b_lm_coda, ask=FALSE)
gelman.diag(b_lm_coda)
plot(b_lm_coda, ask = FALSE)

eap_beta <- matrix(colMeans(smps$beta), ncol = 1)

b_lm_age <- X %*% eap_beta

n_dens <-
  ggplot(data = data.frame(samp = smps$y_new), aes(x = samp)) + geom_density()
n_dens

pred_eap <- mean(smps$y_new)
pred_map <- mlv(smps$y_new, method = "venter")
pred_f_lm
pred_map
pred_eap

available_mcmc(pattern = "_nuts_")

posterior_lm <- as.array(b_lm_fit)
np_lm <- nuts_params(b_lm_fit)
lp_lm <- log_posterior(b_lm_fit)

mcmc_nuts_divergence(np_lm, lp_lm)
mcmc_nuts_acceptance(np_lm, lp_lm)

pairs(b_lm_fit, pars = c("beta[1]", "beta[2]", "beta[3]"))
pairs(b_lm_fit, pars = c("beta[4]", "beta[5]", "beta[6]"))
pairs(b_lm_fit, pars = c("beta[7]", "beta[8]", "beta[9]", "beta[10]"))
