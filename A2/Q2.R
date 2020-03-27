library(data.table)
library(ggplot2)
library(dplyr)

library(rstan)
rstan_options(auto_write = TRUE)
#options(mc.cores = parallel::detectCores())
Sys.setenv(LOCAL_CPPFLAGS = '-march=corei7 -mtune=corei7')
options(mc.cores = 1)
library(rstanarm)
library(coda)
library(bayesplot)

#####
#loading and eda
avalanches_prop <- fread(file = "data/Avalanches_part2.csv")
#avalanches_prop[, Event_ID := NULL]
avalanches_prop[, Snow_meters := Snow_total / 100]
avalanches_prop[, Snow_fnights := Snow_days / 14]
avalanches_prop[, Year := Season]
avalanches_prop[, death_prop := Deaths / Hit]
avalanches_prop[, Geo_space := as.factor(Geo_space)]
avalanches_prop[, Rec.station := as.factor(Rec.station)]
cor(avalanches_prop[, .(Season, Snow_meters, Snow_fnights)])
#####
stan_binomial_glm_reff <-
  stan_model(file = "stan/binomial_glm_randomeffects.stan")

submin <- function(x) {
  m <- min(x)
  x <- x - m
  attributes(x) <- list("scaled:submin" = m)
  return(x)
}

probcomp_geq <- function(x, value){
  mean(x >= value)
}

probcomp_leq <- function(x, value){
  mean(x <= value)
}

cont_vars <- c("Snow_meters", "Snow_fnights")#variables to centre
avalanches_prop[, (cont_vars) := lapply(.SD, scale, scale = FALSE), .SDcols = cont_vars]#centre variables
tm_vars <- c("Season")
avalanches_prop[, (tm_vars) := lapply(.SD, submin), .SDcols = tm_vars]


X_fixedeff <-
  model.matrix(death_prop ~ Season + Snow_meters + Snow_fnights - 1, data = avalanches_prop)
X_randomeff <-
  model.matrix(death_prop ~ Geo_space - 1, data = avalanches_prop)
success <- avalanches_prop[, Deaths]
trials <- avalanches_prop[, Hit]


stan_binomial_glm_reff_data <-
  list(
    success = success,
    trials = trials,
    X_f = X_fixedeff,
    X_r = X_randomeff,
    N = length(success),
    P_f = ncol(X_fixedeff),
    P_r = ncol(X_randomeff),
    n_params = c(0, sqrt(10))
  )

stan_binomial_glm_reff_s <-
  sampling(
    stan_binomial_glm_reff,
    data = stan_binomial_glm_reff_data,
    chains = 4,
    control = list(adapt_delta = 0.95),
    iter = 1e4,
    init_r = 0.1
  )

post_params_rand <-
  extract(stan_binomial_glm_reff_s, c("beta_r"))[[1]]
post_params_fixed <-
  extract(stan_binomial_glm_reff_s, c("beta_f"))[[1]]
post_params <- cbind(post_params_fixed, post_params_rand)
colnames(post_params) <-
  c(colnames(X_fixedeff), colnames(X_randomeff))
ilogit_post_params <- plogis(post_params)
apply(ilogit_post_params, 2, summary)
apply(post_params, 2, summary)

dpp_rand <- extract(stan_binomial_glm_reff_s, "data_ppred")[[1]]
dpp_prop <- apply(dpp_rand, 1, "/", avalanches_prop$Hit)
apply(dpp_prop, 1, summary)

reff_coda <-
  As.mcmc.list(stan_binomial_glm_reff_s, pars = c("beta_r", "beta_f"))
gelman.plot(reff_coda, ask = FALSE)

plot_diag_objects <- function(stanfit) {
  list(
    post = as.array(stanfit),
    lp = log_posterior(stanfit),
    np = nuts_params(stanfit)
  )
}

plot_diag <- function(stanfit, pars) {
  ps <- vars(starts_with(pars))
  post <- as.array(stanfit)
  lp <- log_posterior(stanfit)
  np <- nuts_params(stanfit)
  p1 <- mcmc_parcoord(post, np = np, pars = ps)
  p2 <- mcmc_pairs(post, np = np, pars = ps)
  p3 <- mcmc_trace(post, pars = ps, np = np)
  p4 <- mcmc_nuts_divergence(np, lp)
  p5 <- mcmc_nuts_energy(np)
  list(p1, p2, p3, p4, p5)
}

#mcmc_trace(stan_binomial_glm_reff_s, pars = vars(starts_with("beta")))

#####
#sans snow fortnights
varofint <- avalanches_prop[(Rec.station %in% c(1, 8, 10)) & (Year %in% c(2015, 2018))]
ids <- unique(varofint, by = c("Rec.station", "Year"))$Event_ID
index <- which(avalanches_prop$Event_ID %in% ids)

X_f_nsf <-
  model.matrix(death_prop ~ Season + Snow_meters - 1, data = avalanches_prop)

stan_binomial_glm_reff_nsf_data <-
  list(
    success = success,
    trials = trials,
    X_f = X_f_nsf,
    X_r = X_randomeff,
    N = length(success),
    P_f = ncol(X_f_nsf),
    P_r = ncol(X_randomeff),
    n_params = c(0, sqrt(10))
  )

stan_binomial_glm_reff_nsf_s <-
  sampling(
    stan_binomial_glm_reff,
    data = stan_binomial_glm_reff_nsf_data,
    chains = 4,
    control = list(adapt_delta = 0.95),
    iter = 10000,
    init_r = 0.1
  )



post_params_rand_ns <-
  extract(stan_binomial_glm_reff_nsf_s, c("beta_r"))[[1]]
post_params_fixed_ns <-
  extract(stan_binomial_glm_reff_nsf_s, c("beta_f"))[[1]]
post_params_ns <- cbind(post_params_fixed_ns, post_params_rand_ns)
colnames(post_params_ns) <-
  c(colnames(X_f_nsf), colnames(X_randomeff))
ilogit_post_params_ns <- plogis(post_params_ns)
apply(ilogit_post_params_ns, 2, summary)
apply(post_params_ns, 2, summary)

dpp_rand_nf <- extract(stan_binomial_glm_reff_nsf_s, "data_prop")[[1]]
apply(dpp_rand_nf, 2, summary)
dpp_ofint <- dpp_rand_nf[,index]
apply(dpp_ofint, 2, mean)
apply(dpp_ofint, 2, quantile, c(0.025, 0.975))
apply(dpp_ofint > 0.6, 2, mean)

#####
#hierarchical on station, sans snow fortnights
X_r_station <-
  model.matrix(death_prop ~ Rec.station - 1, data = avalanches_prop)

stan_binomial_glm_reff_station_data <-
  list(
    success = success,
    trials = trials,
    X_f = X_f_nsf,
    X_r = X_r_station,
    N = length(success),
    P_f = ncol(X_f_nsf),
    P_r = ncol(X_r_station),
    n_params = c(0, sqrt(10))
  )

stan_binomial_glm_reff_station_s <-
  sampling(
    stan_binomial_glm_reff,
    data = stan_binomial_glm_reff_station_data,
    chains = 4,
    control = list(adapt_delta = 0.9),
    iter = 10000#,
    #init_r = 0.1
  )

post_params_rand_ns_stat <-
  extract(stan_binomial_glm_reff_station_s, c("beta_r"))[[1]]
post_params_fixed_ns_stat <-
  extract(stan_binomial_glm_reff_station_s, c("beta_f"))[[1]]
post_params_ns_stat <- cbind(post_params_fixed_ns_stat, post_params_rand_ns_stat)
colnames(post_params_ns_stat) <-
  c(colnames(X_f_nsf), colnames(X_r_station))
ilogit_post_params_ns_stat <- plogis(post_params_ns_stat)
apply(ilogit_post_params_ns_stat, 2, summary)
apply(post_params_ns_stat, 2, summary)

dpp_rand_ns_stat <- extract(stan_binomial_glm_reff_station_s, "data_prop")[[1]]
apply(dpp_rand_ns_stat, 2, summary)
dpp_ofintns_stat <- dpp_rand_ns_stat[,index]
apply(dpp_ofintns_stat, 2, mean)
apply(dpp_ofintns_stat, 2, quantile, c(0.025, 0.975))
apply(dpp_ofintns_stat > 0.6, 2, mean)
