library(data.table)
library(ggplot2)

library(rjags)
library(coda)
library(bayesplot)

#####
#loading and eda
avalanches_prop <- fread(file = "data/Avalanches_part2.csv")
avalanches_prop[, Event_ID := NULL]
avalanches_prop[, Snow_meters := Snow_total / 100]
avalanches_prop[, Snow_fnights := Snow_days / 14]
avalanches_prop[, death_prop := Deaths / Hit]
avalanches_prop[, Geo_space := as.factor(Geo_space)]
avalanches_prop[, Rec.station := as.factor(Rec.station)]
cor(avalanches_prop[, .(Season, Snow_meters, Snow_fnights)])
#####

submin <- function(x) {
  m <- min(x)
  x <- x - m
  attributes(x) <- list("scaled:submin" = m)
  return(x)
}

cont_vars <- c("Snow_meters", "Snow_fnights")#variables to centre
avalanches_prop[, (cont_vars) := lapply(.SD, scale, scale = FALSE), .SDcols = cont_vars]#centre variables
tm_vars <- c("Season")
avalanches_prop[, (tm_vars) := lapply(.SD, submin), .SDcols = tm_vars]

snow <- avalanches_prop$Snow_meters
fnight <- avalanches_prop$Snow_fnights
season <- avalanches_prop$Season
n_eff <- length(unique(avalanches_prop$Geo_space))
eff <- as.integer(avalanches_prop$Geo_space)
n <- nrow(avalanches_prop)
deaths <- as.integer(avalanches_prop$Deaths)
hit <- as.integer(avalanches_prop$Hit)

bglm_data <-
  list(
    n = n,
    snow = snow,
    fnight = fnight,
    season = season,
    n_eff = n_eff,
    eff = eff,
    deaths = deaths,
    hit = hit
  )

res.a <-
  jags.model(
    file = "jags/binom_reff.jags",
    data = bglm_data,
    n.chains = 4,
    quiet = T
  )
update(res.a, n.iter = 1e4)
res.b <-
  coda.samples(
    res.a,
    variable.names = c("beta_snow", "beta_season", "beta_fnight", "reff"),
    n.iter = 1e4
  )

summary(res.b)
#####
snow <- avalanches_prop$Snow_meters
season <- avalanches_prop$Season
n_eff <- length(unique(avalanches_prop$Geo_space))
eff <- as.integer(avalanches_prop$Geo_space)
n <- nrow(avalanches_prop)
deaths <- as.integer(avalanches_prop$Deaths)
hit <- as.integer(avalanches_prop$Hit)

bglm_data_nf <-
  list(
    n = n,
    snow = snow,
    season = season,
    n_eff = n_eff,
    eff = eff,
    deaths = deaths,
    hit = hit
  )

res.a_nf <-
  jags.model(
    file = "jags/binom_reff_nofn.jags",
    data = bglm_data_nf,
    n.chains = 4,
    quiet = T
  )
update(res.a_nf, n.iter = 1e4)
res.b_nf <-
  coda.samples(
    res.a_nf,
    variable.names = c("beta_snow", "beta_season", "reff"),
    n.iter = 1e4
  )

summary(res.b_nf)
#####