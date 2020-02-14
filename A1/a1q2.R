library(tidyverse)
library(HDInterval)
library(extraDistr)

par(mfrow = c(1, 1))

waiting_times  <-
  c(
    0.8,
    0.8,
    1.3,
    1.5,
    1.8,
    1.9,
    1.9,
    2.1,
    2.6,
    2.7,
    2.9,
    3.1,
    3.2,
    3.3,
    3.5,
    3.6,
    4.0,
    4.1,
    4.2,
    4.2,
    4.3,
    4.3,
    4.4,
    4.4,
    4.6,
    4.7,
    4.7,
    4.8,
    4.9,
    4.9,
    5,
    5.3,
    5.5,
    5.7,
    5.7,
    6.1,
    6.2,
    6.2,
    6.2,
    6.3,
    6.7,
    6.9,
    7.1,
    7.1,
    7.1,
    7.1,
    7.4,
    7.6,
    7.7,
    8,
    8.2,
    8.6,
    8.6,
    8.6,
    8.8,
    8.8,
    8.9,
    8.9,
    9.5,
    9.6,
    9.7,
    9.8,
    10.7,
    10.9,
    11,
    11,
    11.1,
    11.2,
    11.2,
    11.5,
    11.9,
    12.4,
    12.5,
    12.9,
    13,
    13.1,
    13.3,
    13.6,
    13.7,
    13.9,
    14.1,
    15.4,
    15.4,
    17.3,
    17.3,
    18.1,
    18.2,
    18.4,
    18.9,
    19,
    19.9,
    20.6,
    21.3,
    21.4,
    21.9,
    23.0,
    27,
    31.6,
    33.1,
    38.5
  )

gammaShRaFromMeanSD = function(mean , sd) {
  if (mean <= 0)
    stop("mean must be > 0")
  if (sd <= 0)
    stop("sd must be > 0")
  shape = mean ^ 2 / sd ^ 2
  rate = mean / sd ^ 2
  return(list(shape = shape , rate = rate))
}

gammaShRaFromModeSD = function(mode , sd) {
  if (mode <= 0)
    stop("mode must be > 0")
  if (sd <= 0)
    stop("sd must be > 0")
  rate = (mode + sqrt(mode ^ 2 + 4 * sd ^ 2)) / (2 * sd ^ 2)
  shape = 1 + mode * rate
  return(list(shape = shape , rate = rate))
}

iGammaShScFromMeanVar = function(mean, var) {
  shape <- (mean ^ 2 / var) + 2
  scale <- (mean ^ 3 / var) + mean
  return(list(shape = shape, scale = scale))
}

ex_1_prior <- iGammaShScFromMeanVar(7.5, (0.25 * (10 - 5)) ^ 2)
ex_2_prior <- iGammaShScFromMeanVar(12.5, (0.25 * (25 - 0)) ^ 2)
#prior is 0.5G(9, 0.6) + 0.5G(1,0.08)

n <- length(waiting_times)
mean_wait <- mean(waiting_times)

#posterior is Q*G(9+n, 0.6+n*xbar) + (1-Q)*G(12.5+n, 12.5+n*xbar)

p <- 0.5

find_mix_coefs_2gamma <- function(a1, b1, a2, b2, q1, n, xbar) {
  q2 <- 1 - q1
  a1post <- a1 + n
  a2post <- a2 + n
  b1post <- b1 + (n * xbar)
  b2post <- b2 + (n * xbar)
  log_c1 <-
    (a1 * log(b1) - lgamma(a1)) + (lgamma(a1post) - (a1post * log(b1post)))
  log_c2 <-
    (a2 * log(b2) - lgamma(a2)) + (lgamma(a2post) - (a2post * log(b2post)))
  log_Q <-
    (log(q1) + log_c1) - (log(q1) + log_c1 + log(1 + q2 / q1 * exp(log_c2 -
                                                                     log_c1)))
  #log_Q <- log(q1) + log_c1 - log(q1 * exp(log_c1) + q2 * exp(log_c2))
  Q <- exp(log_Q)
  # c1 <- ((b1^a1)/(gamma(a1))) * ((gamma(a1post))/ (b1post^a1post))
  # c2 <- ((b2^a2)/(gamma(a2))) * ((gamma(a2post))/ (b2post^a2post))
  # Q <- (q1 * c1)/(q1*c1 + q2*c2)
  return(c(Q, 1 - Q))
}

mcoefs <-
  find_mix_coefs_2gamma(ex_1_prior[[1]],
                        ex_1_prior[[2]],
                        ex_2_prior[[1]],
                        ex_2_prior[[2]],
                        p,
                        n,
                        mean_wait)

lambda <- seq(0, 1, len = 5000)

likelihood <-
  as.matrix(apply(as.array(waiting_times), 1, dexp, lambda))
likelihood <- apply(likelihood, 1, prod)

area = sfsmisc::integrate.xy(lambda, likelihood)
const = 1 / area
likelihood <- const * likelihood

prior <-
  p * dgamma(lambda, ex_1_prior[[1]], ex_1_prior[[2]]) + (1 - p) * dgamma(lambda, ex_2_prior[[1]], ex_2_prior[[2]])
posterior <-
  mcoefs[[1]] * dgamma(lambda, ex_1_prior[[1]] + n, ex_1_prior[[2]] + n *
                         mean_wait) + mcoefs[[2]] * dgamma(lambda, ex_2_prior[[1]] + n, ex_2_prior[[2]] + n *
                                                             mean_wait)

plot(
  lambda,
  prior,
  col = "darkgreen",
  ylab = "Density",
  xlab = expression(lambda),
  type = "l",
  ylim = c(0, 45),
  lwd = 2
)
lines(lambda, likelihood, col = "blue2", lwd = 2)
lines(lambda, posterior, col = "red", lwd = 2)
legend(
  "topright",
  legend = c("prior", "scaled likelihood", "posterior"),
  lty = c(1, 1, 1),
  lwd = c(2, 2, 2),
  col = c("darkgreen", "blue2", "red"),
  bty = "n"
)

#prior is relatively flat, posterior is very compatible with data and not very influenced by prior

post_mean <-
  mcoefs[[1]] * ((ex_1_prior[[1]] + n) / (ex_1_prior[[2]] + n * mean_wait)) + mcoefs[[2]] *
  ((ex_2_prior[[1]] + n) / (ex_2_prior[[2]] + n * mean_wait))
post_mean_wait <- 1 / post_mean

posterior_ci <-
  mcoefs[[1]] * qgamma(c(0.025, 0.975), ex_1_prior[[1]] + n, ex_1_prior[[2]] + n *
                         mean_wait) + mcoefs[[2]] * qgamma(c(0.025, 0.975), ex_2_prior[[1]] + n, ex_2_prior[[2]] + n *
                                                             mean_wait)

pci_area <-
  sfsmisc::integrate.xy(lambda, posterior, posterior_ci[[1]], posterior_ci[[2]])


time <- seq(0, 30, len = 5000)

post_wait <-
  mcoefs[[1]] * dinvgamma(time, ex_1_prior[[1]] + n, ex_1_prior[[2]] + n *
                            mean_wait) + mcoefs[[2]] * dinvgamma(time, ex_2_prior[[1]] + n, ex_2_prior[[2]] + n *
                                                                   mean_wait)

post_wait_ci <-
  mcoefs[[1]] * qinvgamma(c(0.025, 0.975), ex_1_prior[[1]] + n, ex_1_prior[[2]] + n *
                            mean_wait) + mcoefs[[2]] * qinvgamma(c(0.025, 0.975), ex_2_prior[[1]] + n, ex_2_prior[[2]] + n *
                                                                   mean_wait)
post_wait_20 <-
  mcoefs[[1]] * pinvgamma(20, ex_1_prior[[1]] + n, ex_1_prior[[2]] + n *
                            mean_wait) + mcoefs[[2]] * pinvgamma(20, ex_2_prior[[1]] + n, ex_2_prior[[2]] + n *
                                                                   mean_wait)

post_wait_sim <-
  mcoefs[[1]] * rinvgamma(1e5, ex_1_prior[[1]] + n, ex_1_prior[[2]] + n *
                            mean_wait) + mcoefs[[2]] * rinvgamma(1e5, ex_2_prior[[1]] + n, ex_2_prior[[2]] + n *
                                                                   mean_wait)

N <- 1e6

components <-
  sample(
    1:2,
    prob = c(mcoefs[[1]], mcoefs[[2]]),
    size = N,
    replace = TRUE
  )
alphas <- c(ex_1_prior[[1]] + n, ex_2_prior[[1]] + n)
betas <-
  c(ex_1_prior[[2]] + n * mean_wait, ex_2_prior[[2]] + n * mean_wait)
samples <- rinvgamma(N, alphas[components], betas[components])

plot(density(samples))
mean(samples)
quantile(samples, c(0.025, 0.975))

post_pred <-
  Renext::rlomax(N, betas[components], alphas[components])

plot(density(post_pred))
mean(post_pred)
quantile(post_pred, c(0.025, 0.975))
sum(post_pred > 20) / N

library(coda)
library(rstan)

rstan_options(auto_write = TRUE)

b_mm_data <-
  list(
    K = 2,
    N = n,
    y = waiting_times,
    theta = c(0.5, 0.5),
    alpha = c(6, 38),
    beta = c(62.5, 277.5)
  )

# b_data <- list(N=n, y = waiting_times, gprior=c(38,278))
#
# b_m <- stan_model(file = 'a1q2simp.stan')
#
# b_m_fit <- sampling(
#   b_m,
#   data = b_data,
#   chains = 7,
#   control = list(adapt_delta = 0.8),
#   iter = 1e5
# )
# plot(b_m_fit)
# b_m_fit
#
# smps <- extract(b_m_fit)
# tp <- traceplot(b_m_fit, pars = c("lambda"))
# tp

# b_m_coda <- As.mcmc.list(b_m_fit)
# gelman.plot(b_m_coda, ask=FALSE)
# gelman.diag(b_m_coda)
# plot(density(smps$postdraw))
# sum(smps$postdraw > 20) / (length(smps$postdraw))
# mean(smps$postdraw)
# quantile(smps$postdraw, c(0.025, 0.975))
# plot(density(smps$ewt))
# sum(smps$ewt > 20) / (length(smps$ewt))
# mean(smps$ewt)
# quantile(smps$ewt, c(0.025, 0.975))

b_mm <- stan_model(file = 'a1q2.stan')

b_mm_fit <- sampling(
  b_mm,
  data = b_mm_data,
  chains = 7,
  control = list(adapt_delta = 0.8),
  iter = 40000
)
plot(b_mm_fit)
b_mm_fit

smps <- extract(b_mm_fit)
tp <- traceplot(b_mm_fit, pars = c("lambda", "mwt"))
tp

b_mm_coda <- As.mcmc.list(b_mm_fit)
gelman.plot(b_mm_coda, ask = FALSE)
gelman.diag(b_mm_coda)

plot(density(smps$postdraw))
sum(smps$postdraw > 20) / (length(smps$postdraw))
mean(smps$postdraw)
quantile(smps$postdraw, c(0.025, 0.975))

plot(density(smps$mwt))
sum(smps$mwt > 20) / (length(smps$mwt))
mean(smps$mwt)
quantile(smps$mwt, c(0.025, 0.975))


require(rjags)

model = jags.model(
  file = "a1q2jags.jags",
  data =
    list(
      y = waiting_times,
      n = n,
      a = c(38, 6),
      b = c(277.5, 62.5),
      p = c(0.5, 0.5)
    ),
  n.chains = 10
)

# Burnin for 1000 samples
update(model, 100000, progress.bar = "none")

# Running the model
res = coda.samples(
  model,
  variable.names = c("lambda", "ypred", "texp"),
  n.iter = 200000,
  progress.bar = "none"
)
summary(res)

qqplot(smps$postdraw, waiting_times, xlim = c(0, 40))
qqplot(post_pred, waiting_times, xlim = c(0, 40))
qqplot(res[[1]][, "ypred"], waiting_times, xlim = c(0, 40))


b_mm_lind <- stan_model(file = 'a1q2lindley.stan')

b_mm_lind_fit <- sampling(
  b_mm_lind,
  data = b_mm_data,
  chains = 7,
  control = list(adapt_delta = 0.8),
  iter = 40000
)
plot(b_mm_lind_fit)
b_mm_lind_fit

smps_lind <- extract(b_mm_lind_fit)
tp_lind <- traceplot(b_mm_lind_fit, pars = c("lambda", "mwt"))
tp_lind

b_mm_lind_coda <- As.mcmc.list(b_mm_lind_fit)
gelman.plot(b_mm_lind_coda, ask = FALSE)
gelman.diag(b_mm_lind_coda)

plot(density(smps_lind$postdraw))
sum(smps_lind$postdraw > 20) / (length(smps_lind$postdraw))
mean(smps_lind$postdraw)
quantile(smps_lind$postdraw, c(0.025, 0.975))

plot(density(smps_lind$mwt))
sum(smps_lind$mwt > 20) / (length(smps_lind$mwt))
mean(smps_lind$mwt)
quantile(smps_lind$mwt, c(0.025, 0.975))

qqplot(smps_lind$postdraw, waiting_times, xlim = c(0, 40))
qqplot(smps$postdraw, waiting_times, xlim = c(0, 40))
ks.test(smps_lind$postdraw, waiting_times)
ks.test(smps$postdraw, waiting_times)


dens <- data.frame(postdraw = post_pred)
ggdens <-
  ggplot(data = dens, aes(x = postdraw)) +
  geom_density(size = 1, color = "darkblue", fill = "lightblue") +
  theme_minimal() +
  labs(x = "Posterior Waiting Time", y = "Density of Posterior Predictive", title = "Posterior Predictive Density (Exponential Likelihood)") +
  geom_vline(
    xintercept = mean(dens$postdraw),
    size = 2,
    color = "red"
  ) +
  geom_vline(
    xintercept = quantile(dens$postdraw, c(0.025, 0.975)),
    color = "darkred",
    size = 1
  ) + coord_cartesian(xlim = c(0, 50))
ggdens

lambda <- seq(0, 1, len = 5000)

likelihood <-
  as.matrix(apply(as.array(waiting_times), 1, dexp, lambda))
likelihood <- apply(likelihood, 1, prod)

area = sfsmisc::integrate.xy(lambda, likelihood)
const = 1 / area
likelihood <- const * likelihood

prior <-
  p * dgamma(lambda, ex_1_prior[[1]], ex_1_prior[[2]]) + (1 - p) * dgamma(lambda, ex_2_prior[[1]], ex_2_prior[[2]])
posterior <-
  mcoefs[[1]] * dgamma(lambda, ex_1_prior[[1]] + n, ex_1_prior[[2]] + n *
                         mean_wait) + mcoefs[[2]] * dgamma(lambda, ex_2_prior[[1]] + n, ex_2_prior[[2]] + n *
                                                             mean_wait)

prlikpost <-
  reshape2::melt(list(
    Prior = prior,
    Likelihood = likelihood,
    Posterior = posterior
  ))
ggcompplot <-
  ggplot(data = prlikpost) + geom_area(
    aes(
      x = rep(lambda, 3),
      y = value,
      colour = L1,
      fill = L1
    ),
    alpha = 0.3,
    size = 0.5,
    position = "identity"
  ) +
  theme_minimal() +
  labs(x = "Lambda", y = "Density", title = "Prior-Posterior Compatibility Plot")
ggcompplot

mcoefs[[1]] *
  Renext::qlomax(c(0.025, 0.975), betas[1], alphas[1]) +
  mcoefs[[2]] *
  Renext::qlomax(c(0.025, 0.975), betas[2], alphas[2])

dens_lind <- data.frame(postdraw = smps_lind$postdraw)
ggdens_lind <-
  ggplot(data = dens, aes(x = postdraw)) +
  geom_density(size = 1, color = "darkblue", fill = "lightblue") +
  theme_minimal() +
  labs(x = "Posterior Waiting Time", y = "Density of Posterior Predictive", title = "Posterior Predictive Density (Lindley Likelihood)") +
  geom_vline(
    xintercept = mean(dens_lind$postdraw),
    size = 2,
    color = "red"
  ) +
  geom_vline(
    xintercept = quantile(dens_lind$postdraw, c(0.025, 0.975)),
    color = "darkred",
    size = 1
  ) + coord_cartesian(xlim = c(0, 50))
ggdens_lind


quantile_set <- seq(0.01, 0.99, by = 0.01)
data_q <- quantile(waiting_times, quantile_set)
exp_q <- quantile(smps$postdraw, quantile_set)
lind_q <- quantile(smps_lind$postdraw, quantile_set)

quantile_df <- data.frame(obs = data_q, exp = exp_q, lind = lind_q)
ggqq <- ggplot(data = quantile_df, aes(y = obs)) + coord_cartesian(xlim = c(0, 45), ylim = c(0, 45)) +
  geom_point(aes(x = exp, colour = "Exponential")) + 
  geom_point(aes(x = lind, colour = "Lindley")) + 
  theme_minimal() +
  labs(x = "Fitted Quantiles", y = "Sample Quantiles") + geom_abline(slope = 1, intercept = 0) + scale_colour_discrete(name = "Distribution")
ggqq

