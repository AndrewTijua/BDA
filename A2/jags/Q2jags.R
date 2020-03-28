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

dic.samples(model = res.a_nf,
            n.iter = 1e4,
            type = 'pD')
#####
snow <- avalanches_prop$Snow_meters
season <- avalanches_prop$Season
#n_eff <- length(unique(avalanches_prop$Geo_space))
#eff <- as.integer(avalanches_prop$Geo_space)
eff_stat <- as.integer(avalanches_prop$Rec.station)
n_eff_stat <- length(unique(eff_stat))
n <- nrow(avalanches_prop)
deaths <- as.integer(avalanches_prop$Deaths)
hit <- as.integer(avalanches_prop$Hit)

bglm_data_nf_stat <-
  list(
    n = n,
    snow = snow,
    season = season,
    n_eff = n_eff_stat,
    eff = eff_stat,
    deaths = deaths,
    hit = hit
  )

res.a_nf_stat <-
  jags.model(
    file = "jags/binom_reff_nofn.jags",
    data = bglm_data_nf_stat,
    n.chains = 4,
    quiet = T
  )
update(res.a_nf_stat, n.iter = 1e4)
res.b_nf_stat <-
  coda.samples(
    res.a_nf_stat,
    variable.names = c("beta_snow", "beta_season", "reff"),
    n.iter = 1e4
  )

summary(res.b_nf_stat)

dic.samples(model = res.a_nf_stat,
            n.iter = 1e4,
            type = 'pD')
#####
snow <- avalanches_prop$Snow_meters
season <- avalanches_prop$Season
#n_eff <- length(unique(avalanches_prop$Geo_space))
#eff <- as.integer(avalanches_prop$Geo_space)
stations <- as.integer(avalanches_prop$Rec.station)
geos <- as.integer(avalanches_prop$Geo_space)
n_station <- length(unique(stations))
n_geo <- length(unique(geos))
n <- nrow(avalanches_prop)
deaths <- as.integer(avalanches_prop$Deaths)
hit <- as.integer(avalanches_prop$Hit)

bglm_data_nf_statgeo <-
  list(
    n = n,
    snow = snow,
    season = season,
    geos = geos,
    stations = stations,
    n_station = n_station,
    n_geo = n_geo,
    deaths = deaths,
    hit = hit
  )

res.a_nf_statgeo <-
  jags.model(
    file = "jags/binom_doublereff.jags",
    data = bglm_data_nf_statgeo,
    n.chains = 4,
    quiet = T
  )
update(res.a_nf_statgeo, n.iter = 1e4)
res.b_nf_statgeo <-
  coda.samples(
    res.a_nf_statgeo,
    variable.names = c("beta_snow", "beta_season", "r_eff_geo", "r_eff_statgeo"),
    n.iter = 1e4
  )

summary(res.b_nf_statgeo)

dic.samples(model = res.a_nf_statgeo,
            n.iter = 1e4,
            type = 'pD')
