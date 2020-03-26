library(data.table)
library(ggplot2)

library(rjags)
library(coda)
library(bayesplot)


#####
#a
avalanches <- fread(file = "data/Avalanches.csv")
#avalanches <- avalanches[Rep.events > 0]
avalanches[, ':=' (EADS1 = (Season >= 1994 &
                              Season <= 2003),
                   EADS2 = (Season >= 2004))]

avalanches[Season %in% c(1986, 1994, 2004)]

avalanches[, EWS := 1 + EADS1 + 2 * EADS2]
avalanches[, EWS := as.factor(EWS)]

d_offset <- rep(0, nrow(avalanches))

pglm_data <-
  list(
    n = nrow(avalanches),
    w1 = avalanches$EADS1,
    w2 = avalanches$EADS2,
    rep = avalanches$Rep.events,
    death = avalanches$Deaths,
    offset = d_offset
  )

res.a <-
  jags.model(
    file = "jags/poisson.jags",
    data = pglm_data,
    n.chains = 4,
    quiet = T
  )
update(res.a, n.iter = 1e4)
res.b <-
  coda.samples(
    res.a,
    variable.names = c("intercept", "beta_w1", "beta_w2", "beta_rep"),
    n.iter = 1e4
  )
summary(res.b)
dic.samples(model = res.a,
            n.iter = 1e4,
            type = 'pD')

sm <- rbindlist(lapply(res.b, as.data.frame))

news_1_j <- mean(exp(sm$intercept) > 1)
news_2_j <- mean(exp(sm$beta_w1 + sm$intercept) > 1)
news_3_j <- mean(exp(sm$beta_w2 + sm$intercept) > 1)

res.a.ev <-
  jags.model(
    file = "jags/poisson_exvar.jags",
    data = pglm_data,
    n.chains = 4,
    quiet = T
  )
update(res.a, n.iter = 1e4)
res.b.ev <-
  coda.samples(
    res.a.ev,
    variable.names = c("beta_w1", "beta_w2", "beta_rep", "theta"),
    n.iter = 1e4
  )
summary(res.b.ev)
dic.samples(model = res.a.ev,
            n.iter = 1e4,
            type = 'pD')
