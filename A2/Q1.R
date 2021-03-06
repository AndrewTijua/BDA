library(data.table)
library(ggplot2)

library(rstan)
rstan_options(auto_write = TRUE)
#options(mc.cores = parallel::detectCores())
Sys.setenv(LOCAL_CPPFLAGS = '-march=corei7 -mtune=corei7')
options(mc.cores = 1)
library(rstanarm)
library(coda)
library(bayesplot)


#####
#a
avalanches <- fread(file = "data/Avalanches.csv")
avalanches[, ':=' (EADS1 = (Season >= 1994 &
                              Season <= 2003),
                   EADS2 = (Season >= 2004))]

avalanches[Season %in% c(1986, 1994, 2004)]

avalanches[, EWS := 1 + EADS1 + 2 * EADS2]
avalanches[, EWS := as.factor(EWS)]

base_plot <-
  ggplot(data = as.data.frame(avalanches), aes(colour = EWS)) + theme_minimal()
base_plot + geom_line(aes(x = Season, y = Rep.events, group = F))
base_plot + geom_line(aes(x = Season, y = Deaths, group = F))
base_plot + geom_boxplot(aes(x = EWS, y = Deaths), colour = "black")

#avalanches <- avalanches[Rep.events > 0]
cor_boot <- function(data, index) {
  dt_s <- data[index, ]
  return(cor(dt_s))
}

cor(avalanches[(EADS1 == FALSE &
                  EADS2 == FALSE), .(Rep.events, Deaths)])
cor(avalanches[EADS1 == TRUE, .(Rep.events, Deaths)])
cor(avalanches[EADS2 == TRUE, .(Rep.events, Deaths)])

bs1 <- boot::boot(avalanches[(EADS1 == FALSE &
                                EADS2 == FALSE),
                             .(Rep.events, Deaths)]
                  , cor_boot, R = 1e3)
bs2 <- boot::boot(avalanches[(EADS1 == TRUE),
                             .(Rep.events, Deaths)]
                  , cor_boot, R = 1e3)
bs3 <- boot::boot(avalanches[(EADS2 == TRUE),
                             .(Rep.events, Deaths)]
                  , cor_boot, R = 1e3)
boot::boot.ci(bs1,
              index = 2,
              type = "perc",
              conf = 0.9)
boot::boot.ci(bs2,
              index = 2,
              type = "perc",
              conf = 0.9)
boot::boot.ci(bs3,
              index = 2,
              type = "perc",
              conf = 0.9)
#####
#b
to_model <- avalanches[, .(Rep.events, Deaths, EADS1, EADS2)]
model_mat <-
  model.matrix(Deaths ~ ., data = to_model)#no intercept as cannot have deaths without avalanche
#d_offset <- log(avalanches$Rep.events)
d_offset <- rep(0, nrow(avalanches))
model_mat <- model_mat[,]
out_names = colnames(model_mat)
#no need to centre as discrete

#new data

# X_new = matrix(c(1, 0, 1, 1, 0, 0, 1, 1, 0, 1, 0, 1),
#                nrow = 4,
#                byrow = T)

X_new = matrix(c(1, 20, 0, 1,
                 1, 1, 0, 0,
                 1, 1, 1, 0,
                 1, 1, 0, 1),
               nrow = 4,
               byrow = T)
#n_offset <- log(c(20, 1, 1, 1))
n_offset <- rep(0, nrow(X_new))

N_new = nrow(X_new)
#check, should be similar
f_glm <-
  glm(Deaths ~ ., data = to_model, family = poisson(link = "log"))


stan_poisson_glm <- stan_model(file = "stan/poisson_glm.stan")
stan_poisson_glm_data <-
  list(
    N = nrow(model_mat),
    P = ncol(model_mat),
    y = avalanches$Deaths,
    X = model_mat,
    n_params = c(0, 1e2),
    N_new = N_new,
    X_new = X_new,
    offset = d_offset,
    offset_new = n_offset
  )


stan_poisson_glm_s <-
  sampling(
    stan_poisson_glm,
    data = stan_poisson_glm_data,
    chains = 7,
    control = list(adapt_delta = 0.8),
    iter = 1e4,
    init_r = 0.1
  )

post_params <- extract(stan_poisson_glm_s, "lambda")[[1]]
colnames(post_params) <- out_names
exp_post_params <- exp(post_params)
apply(exp_post_params, 2, summary)
apply(post_params, 2, summary)

news_1 <- mean(exp(post_params[, 1]) > 1)
news_2 <- mean(exp(post_params[, 1] + post_params[, 2]) > 1)
news_3 <- mean(exp(post_params[, 1] + post_params[, 3]) > 1)


p_pred <- extract(stan_poisson_glm_s, "y_new")[[1]]
mean(p_pred[, 1] < 15)
mean(p_pred[, 2] > 1)
mean(p_pred[, 3] > 1)
mean(p_pred[, 4] > 1)

pp1 <- p_pred[,1] < 15

mean_boot <- function(data, index) {
  dt_s <- data[index]
  return(mean(dt_s))
}

bs4 <- boot::boot(pp1, mean_boot, R = 1e3)
boot::boot.ci(bs4, type = "perc", conf = 0.95)

data_pred <- extract(stan_poisson_glm_s, "data_ppred")[[1]]
apply(data_pred, 2, summary)

dpp_m1_plotdf <-
  data.frame(
    mean = apply(data_pred, 2, mean),
    lq = apply(data_pred, 2, quantile, 0.05),
    uq = apply(data_pred, 2, quantile, 0.95),
    Season = avalanches$Season
  )

lr_data <- extract(stan_poisson_glm_s, "log_rate")[[1]]
r_data <- exp(lr_data)
r_data_pe <- apply(r_data, 1, '/', avalanches$Rep.events)

r_data_b <- r_data_pe[avalanches$Season < 1994]
r_data_do <- r_data_pe[avalanches$EADS1 == TRUE]
r_data_o <- r_data_pe[avalanches$EADS2 == TRUE]

r_data_b <- unlist(r_data_b)
r_data_b <- r_data_b[!is.infinite(r_data_b)]
mean(r_data_b >= 1)
mean(r_data_b)

r_data_do <- unlist(r_data_do)
r_data_do <- r_data_do[!is.infinite(r_data_do)]
mean(r_data_do >= 1)
mean(r_data_do)

r_data_o <- unlist(r_data_o)
r_data_o <- r_data_o[!is.infinite(r_data_o)]
mean(r_data_o >= 1)
mean(r_data_o)
#####
#dic is bad
#formulae taken from https://en.wikipedia.org/wiki/Deviance_information_criterion
plikrar <- function(x, data) {
  sum(dpois(data, x, log = T))
}
sampling_rates <- extract(stan_poisson_glm_s, "rate")[[1]]
sr_like <-
  apply(sampling_rates, 1, plikrar, avalanches$Deaths)#calculate log likelihoods of each sampling
sr_like_mean <-
  mean(sr_like)#calculate mean log likelihood of samples
eap <-
  colMeans(sampling_rates)#calculate posterior means of rates (not parameters)
p_mean_like <-
  sum(dpois(avalanches$Deaths, eap, log = T))#calculate log likelihood of EAP
dbar <- -2 * sr_like_mean#expected deviance
pd <- dbar + 2 * p_mean_like#calculate penalty
dic <- pd + dbar#give dic
#####
#prior checking
# dp_av <- avalanches$Deaths/avalanches$Rep.events
# dp_av <- dp_av[!is.nan(dp_av)]
# m_deaths <- mean(dp_av)
# xm <- dp_av - m_deaths
# lnfactor <- 2/(xm)^2
# inffactor <- dp_av / m_deaths
# beta_p <-
# mfc <- exp(xm * inffactor)
# mfc_p <- plnorm(mfc, 0, 2)
avno <- avalanches$Rep.events
avde <- avalanches$Deaths
mede <- mean(avde)
psi <- avde / mede
beta <- log(psi) / (avno - mean(avno))
psi_p <- dlnorm(psi, 0, 2)
beta_p <- dnorm(beta, 0, (avno - mean(avno)) ^ (-2))
#####
stan_poisson_glm_exvar <-
  stan_model(file = "stan/poisson_glm_exvar.stan")

model_mat <- model_mat[,-1]#messes with exvar
out_names = colnames(model_mat)

# X_new = matrix(c(0, 1, 0, 0, 1, 0, 0, 1),
#                nrow = 4,
#                byrow = T)

X_new = matrix(c(20, 0, 1,
                 1, 0, 0,
                 1, 1, 0,
                 1, 0, 1),
               nrow = 4,
               byrow = T)

#n_offset <- log(c(20, 1, 1, 1))

ym <- data.frame(ym = as.factor(avalanches$Season))
yim <- model.matrix( ~ . - 1, ym)

stan_poisson_glm_exvar_data <-
  list(
    N = nrow(model_mat),
    P = ncol(model_mat),
    y = avalanches$Deaths,
    X = model_mat,
    n_params = c(0, sqrt(10)),
    N_new = N_new,
    X_new = X_new,
    yearindmat = yim,
    N_years = ncol(yim),
    offset = d_offset,
    offset_new = n_offset
  )


stan_poisson_glm_exvar_s <-
  sampling(
    stan_poisson_glm_exvar,
    data = stan_poisson_glm_exvar_data,
    chains = 4,
    control = list(adapt_delta = 0.99, max_treedepth = 15),
    iter = 8000,
    init_r = 0.05,
    pars = c("lambda", "theta", "data_ppred", "rate")
  )

post_params_exvar <-
  extract(stan_poisson_glm_exvar_s, c("lambda"))[[1]]
post_params_theta <- extract(stan_poisson_glm_exvar_s, "theta")[[1]]
colnames(post_params_exvar) <- out_names
names(post_params_theta) <- "theta"

bound <- cbind(post_params_exvar, post_params_theta)
colnames(bound) <- c(out_names, "theta")
apply(exp(bound), 2, summary)

dpp <- extract(stan_poisson_glm_exvar_s, "data_ppred")[[1]]
apply(dpp, 2, summary)

dpp_m2_plotdf <-
  data.frame(
    mean = apply(dpp, 2, mean),
    lq = apply(dpp, 2, quantile, 0.05),
    uq = apply(dpp, 2, quantile, 0.95),
    Season = avalanches$Season
  )
#####
plikrar <- function(x, data) {
  sum(dpois(data, x, log = T))
}
sampling_rates_exv <- extract(stan_poisson_glm_exvar_s, "rate")[[1]]
sr_like_exv <-
  apply(sampling_rates_exv, 1, plikrar, avalanches$Deaths)#calculate log likelihoods of each sampling
sr_like_mean_exv <-
  mean(sr_like_exv)#calculate mean log likelihood of samples
eap_exv <-
  colMeans(sampling_rates_exv)#calculate posterior means of rates (not parameters)
p_mean_like_exv <-
  sum(dpois(avalanches$Deaths, eap_exv, log = T))#calculate log likelihood of EAP
dbar_exv <- -2 * sr_like_mean_exv#expected deviance
pd_exv <- dbar_exv + 2 * p_mean_like_exv#calculate penalty
dic_exv <- pd_exv + dbar_exv#give dic
#####
ggplot(data = dpp_m1_plotdf, aes(x = Season)) + theme_minimal() +
  geom_ribbon(aes(ymin = lq, ymax = uq), alpha = 0.5) + labs(title = "Posterior Predictive for Model 1", y = "Number of Deaths") +
  geom_line(aes(y = mean), size = 2, colour = "red")

ggplot(data = dpp_m2_plotdf, aes(x = Season)) + theme_minimal() +
  geom_ribbon(aes(ymin = lq, ymax = uq), alpha = 0.5) + labs(title = "Posterior Predictive for Model 2 (extra variance)", y = "Number of Deaths") +
  geom_line(aes(y = mean), size = 2, colour = "red")

pp_mod_1 <- as.data.frame(exp_post_params)
pp_mod_1_long <- reshape2::melt(pp_mod_1)
pp_mod_2 <- as.data.frame(exp(bound))
pp_mod_2_long <- reshape2::melt(pp_mod_2)

ggplot(data = pp_mod_1_long, aes(x = variable, y = value)) + theme_minimal() +
  geom_boxplot() + labs(title = "Posterior summaries for model 1", y = "Parameter value", x = "Parameter") + coord_cartesian(ylim = c(0, 3))
ggplot(data = pp_mod_2_long, aes(x = variable, y = value)) + theme_minimal() +
  geom_boxplot() + labs(title = "Posterior summaries for model 2 (extra variance)", y = "Parameter value", x = "Parameter") + coord_cartesian(ylim = c(0, 3))
