## ----setup, include = FALSE------------------------------------------------------------------------------------------------
knitr::opts_chunk$set(echo = FALSE, cache = TRUE, warning = FALSE, message = FALSE)


## ----precode, include = FALSE----------------------------------------------------------------------------------------------
# set wd
setwd("/Users/austin/Desktop/india_mds_csmr/compare_multinomial_poisson")

library(SUMMER);
library(milliondeaths); library(tidyverse);
library(INLA);library(mapmisc);library(spdep);
library(ggpubr); library(broom); library(geofacet);
library(knitr); library(kableExtra);

#### functions used in this code ####


#### set up MDS info and directories ####

# MDS username/password
username <- "schumachera"
password <- "schumachera2402"

# dropbox directory to save data
datadir <- "../../../Dropbox/csmr_india/data/mds"
savedir <- "../../../Dropbox/csmr_india/sandbox"

#### parameters ####

# set years and ages for analysis
years <- 2004:2013
ages <- c(0,1,5)

# dimensions for data aggregation
aggregationDims <- c("none", "state", "year", "year-state")
# aggregationDims <- c("year")


##### load aggregated data (created from 'aggregate_format_data_for_exploration.R') #####
# this loads:
# - dat.df.agg.list: named list where each element is an aggregated data frame of MDS deaths by cause
# - allChildCodes: data frame with all child COD codes
# - allChildCodesBroad: data frame with broad child COD codes
# - india: shapefile for Indian states/regions
# - india.adj: neighborhood file for india states/regions for use with INLA spatial models
# - india.toplot: tidy() version of India map shapefile for plotting with ggplot
# - india_grid_MDS: grid file for use with geo_facet plots in ggplot with the geofacet package

load(paste0(datadir, "/aggregated_formatted_data_lists_and_codes_and_mapping_data.RData"))

##### single cause binomial models (with and without Poisson Trick) #####

# vector of cause names
causenames_v <- sort(allChildCodesBroad$cause)
# causenames_v <- "cause_1C01"

# create list of models to run for each aggregation dimension
aggDimModels <- vector(mode = "list", length = length(aggregationDims))
names(aggDimModels) <- aggregationDims
for (i in 1:length(aggDimModels)) {
    ad <- aggregationDims[i]
    if (ad == "none") aggDimModels[[i]] <- c("int-only")
    if (ad == "state") aggDimModels[[i]] <- c("bym2")
    if (ad == "year") aggDimModels[[i]] <- c("linear", "rw2")
    if (ad == "year-state") aggDimModels[[i]] <- c("bym2_rw2_no-int", "bym2_rw2_typeI", "bym2_rw2_typeIV")
}

##### format data for year only aggregation #####

# set aggregation dimension
aggdim <- "year"

# set data to use
mod_data <- dat.df.agg.list[[aggdim]]

# format data
mod_data$year_num <- mod_data$year - min(mod_data$year) + 1
mod_data$year_cent <- mod_data$year - mean(mod_data$year)
mod_data$year2_num <- mod_data$year2 - min(mod_data$year2) + 1
mod_data$year2_cent <- mod_data$year2 - mean(mod_data$year2)

# set cause
causename <- "cause_1C01"

# set CI level and quantiles for analysis
ci.level <- 0.95
res.quantiles <- c((1-ci.level) - (1-ci.level)/2, 0.5, ci.level + (1-ci.level)/2)

# set reference cause
reference.cause <- "cause_other"


## ----intonly, include=FALSE, cache = TRUE----------------------------------------------------------------------------------
# set model name
model_name <- "none"

# binomial model formula
f.bin <- paste(causename, "~ 1")
    
# hyperpars
fe.prec.bin <- list(prec.intercept = 0,
                    prec = 0)

# run binomial model intercept-only
mod.bin <- inla(as.formula(f.bin),
                control.fixed = fe.prec.bin,
                data = mod_data,
                family = "binomial",
                Ntrials = all.cause,
                quantiles = res.quantiles,
                control.predictor=list(compute=TRUE),
                control.compute=list(config = TRUE))

# create long data for poisson trick model
dat.df.agg <- mod_data %>%
    mutate(cause_other = all.cause - get(causename))
dat.df.agg.long <- dat.df.agg %>%
    pivot_longer(cols = c(all_of(causename), 
                          cause_other,
                          all.cause),
                 names_to = "cause",
                 values_to = "deaths") %>%
    filter(cause %in% c(causename, "cause_other")) %>%
    dplyr::select(-contains("cause_"))

# make cause a numeric vector, with cause 1 equal to the other causes 
#   - (so the coefficients from the binomial and Poisson models match)
dat.df.agg.long$cause_factor <- as.factor(ifelse(dat.df.agg.long$cause == causename, 2, 1))

f.ptrick.fe <- "deaths ~ 1 + cause_factor"
f.ptrick.fe <- paste(f.ptrick.fe, 
                     "+ f(year_num, model = 'iid',
                          hyper = list(prec = list(initial = log(0),
                                                   fixed = TRUE)))")

fe.prec.pois <- list(prec.intercept = 100000000,
                     prec = list(default = 0, 
                                 year_cent = 100000000))

mod.ptrick <- inla(as.formula(f.ptrick.fe),
                   control.fixed = fe.prec.pois,
                   family = "poisson",
                   data = dat.df.agg.long,
                   quantiles = res.quantiles,
                   control.predictor=list(compute=TRUE),
                   control.compute=list(config = TRUE))

allres <- data.frame(model = c("Multinomial", 
                               "Poisson"),
                     csmf_lower = c(expit(mod.bin$summary.fixed$`0.025quant`),
                                    expit(mod.ptrick$summary.fixed["cause_factor2", "0.025quant"])),
                     csmf_med = c(expit(mod.bin$summary.fixed$`0.5quant`),
                                    expit(mod.ptrick$summary.fixed["cause_factor2", "0.5quant"])),
                     csmf_upper = c(expit(mod.bin$summary.fixed$`0.975quant`),
                                    expit(mod.ptrick$summary.fixed["cause_factor2", "0.975quant"])))


## ----tabintonly------------------------------------------------------------------------------------------------------------
## compare results
summ_intonly <- rbind(mod.bin$summary.fixed,
                      mod.ptrick$summary.fixed)
summ_intonly <- cbind(model = c("Multinom", 
                                rep("", nrow(mod.bin$summary.fixed) - 1), 
                                "Pois", 
                                rep("", nrow(mod.ptrick$summary.fixed) - 1)), 
                      coefficient = c("$\\beta_0$","$\\tilde\\beta_0$", "$\\beta_0$"), 
                      summ_intonly)
rownames(summ_intonly) <- NULL
kable(summ_intonly, 
      caption = "Posterior distributions of parameters from Multinomial and Poisson models: fixed intercepts only",
      booktabs = TRUE, format = 'markdown', digits = 5) %>% 
    kable_styling(latex_options = "HOLD_position")


## ----alphajintonly---------------------------------------------------------------------------------------------------------
kable(mod.ptrick$summary.random, 
      caption = "Posterior distributions of nuisance parameters from Poisson model: fixed intercepts only",
      booktabs = TRUE, format = 'markdown', digits = 5) %>% 
    kable_styling(latex_options = "HOLD_position")


## ----get_csmfs_and_plot_intonly--------------------------------------------------------------------------------------------
# ggplot CSMFs
ggplot(data = allres, 
       aes(x = model, y = csmf_med, color = model)) +
    geom_point() +
    geom_errorbar(aes(ymin = csmf_lower, ymax = csmf_upper)) +
    labs(x="model", 
         y="Cause fraction",
         subtitle = paste0("Estimated CSMF (",ci.level*100,"% posterior intervals) from yearly national model with cause-specific intercepts only"),
         title=unique(allChildCodesBroad$broadCodeName[allChildCodesBroad$cause == causename])) +
    theme_light()


## ----lineartime, include=FALSE, cache = TRUE-------------------------------------------------------------------------------
##### format data, run models #####

# set aggregation dimension
aggdim <- "year"

# set data to use
mod_data <- dat.df.agg.list[[aggdim]]

# format data
mod_data$year_num <- mod_data$year - min(mod_data$year) + 1
mod_data$year_cent <- mod_data$year - mean(mod_data$year)
mod_data$year2_num <- mod_data$year2 - min(mod_data$year2) + 1
mod_data$year2_cent <- mod_data$year2 - mean(mod_data$year2)

# set cause
causename <- "cause_1C01"

# binomial model formula
f.bin <- paste(causename, "~ 1")
f.bin <- paste(f.bin, "+ year_cent")
    
# run binomial model intercept-only
mod.bin <- inla(as.formula(f.bin),
                control.fixed = fe.prec.bin,
                data = mod_data,
                family = "binomial",
                Ntrials = all.cause,
                quantiles = res.quantiles,
                control.predictor=list(compute=TRUE),
                control.compute=list(config = TRUE))

# create long data for poisson trick model
dat.df.agg <- mod_data %>%
    mutate(cause_other = all.cause - get(causename))
dat.df.agg.long <- dat.df.agg %>%
    pivot_longer(cols = c(all_of(causename), 
                          cause_other,
                          all.cause),
                 names_to = "cause",
                 values_to = "deaths") %>%
    filter(cause %in% c(causename, "cause_other")) %>%
    dplyr::select(-contains("cause_"))

# make cause a numeric vector, with cause 1 equal to the other causes 
#   - (so the coefficients from the binomial and Poisson models match)
dat.df.agg.long$cause_factor <- as.factor(ifelse(dat.df.agg.long$cause == causename, 2, 1))

f.ptrick.fe <- "deaths ~ 1 + cause_factor"
f.ptrick.fe <- paste(f.ptrick.fe, 
                     "* year_cent")
f.ptrick.fe <- paste(f.ptrick.fe, 
                     "+ f(year_num, model = 'iid',
                          hyper = list(prec = list(initial = log(0),
                                                   fixed = TRUE)))")

# linear combination without alphas for prediction
mat.pred.nonref <- Diagonal(nrow(dat.df.agg.long), 1)[which(dat.df.agg.long$cause != reference.cause), ]
marg.mat <- -1*model.matrix(~-1 + factor(year_num), data = dat.df.agg.long)[which(dat.df.agg.long$cause != reference.cause), ]
lcs <- inla.make.lincombs(
  "Predictor" = mat.pred.nonref,
  "year_num" = marg.mat
)

mod.ptrick <- inla(as.formula(f.ptrick.fe),
                   control.fixed = fe.prec.pois,
                   family = "poisson",
                   data = dat.df.agg.long,
                   lincomb = lcs,
                   quantiles = res.quantiles,
                   control.predictor=list(compute=TRUE),
                   control.compute=list(config = TRUE))

res.pois <- cbind(dat.df.agg.long[dat.df.agg.long$cause != reference.cause, c("cause", "year", "year_num")],
                  expit(mod.ptrick$summary.lincomb.derived[, paste0(res.quantiles, "quant")]))
res.bin <- cbind(mod_data[, c("year", "year_num")], mod.bin$summary.fitted.values[, paste0(res.quantiles, "quant")])
res.bin <- cbind(cause = rep(causename, nrow(res.bin)), res.bin)
res.pois$model <- "Poisson"
res.bin$model <- "Binomial"
allres <- rbind(res.pois, res.bin)


## ----tablinear-------------------------------------------------------------------------------------------------------------
## compare results
summ_linear <- rbind(mod.bin$summary.fixed,
                     mod.ptrick$summary.fixed)
summ_linear <- cbind(model = c("Multinom", rep("", nrow(mod.bin$summary.fixed) - 1), "Pois", rep("", nrow(mod.ptrick$summary.fixed) - 1)), 
                     coefficient = c("$\\beta_0$","$\\beta_1$","$\\tilde\\beta_0$", "$\\beta_0$", "$\\tilde\\beta_1$", "$\\beta_1$"), 
                     summ_linear)
rownames(summ_linear) <- NULL
kable(summ_linear, 
      caption = "Posterior distributions of shared parameters from Multinomial and Poisson models: fixed intercept, linear slope on time",
      booktabs = TRUE, format = 'markdown', digits = 5) %>% 
    kable_styling(latex_options = "HOLD_position")


## ----get_csmfs_and_plot----------------------------------------------------------------------------------------------------
# ggplot CSMFs for linear trend on time
ggplot(data = allres[allres$cause == causename,], 
       aes(x = year, y = `0.5quant`, color = model)) +
    geom_line(alpha = 0.2) +
    geom_ribbon(aes(x = year, ymin = `0.025quant`, ymax = `0.975quant`, fill = model), alpha = 0.2) +
    labs(x="year", 
         y="Cause fraction",
         subtitle = paste0("Estimated CSMF (",ci.level*100,"% posterior intervals) from yearly national model, linear slope on time"),
         title=unique(allChildCodesBroad$broadCodeName[allChildCodesBroad$cause == causename])) +
    theme_light()


## ----iidtime1, cache = TRUE------------------------------------------------------------------------------------------------
#### IID in time ####
model_name <- "iidyear"

# binomial model formula
f.bin <- paste(causename, "~ 1")
f.bin <- paste(f.bin, 
               "+ f(year, model = 'iid',
                    hyper = list(prec = list(prior = 'pc.prec', 
                                             param = c(0.3, 0.01))))")

# run binomial model
mod.bin <- inla(as.formula(f.bin),
                control.fixed = fe.prec.bin,
                data = mod_data,
                family = "binomial",
                Ntrials = all.cause,
                quantiles = res.quantiles,
                control.predictor=list(compute=TRUE),
                control.compute=list(config = TRUE))

# create long data for poisson trick model
dat.df.agg <- mod_data %>%
    mutate(cause_other = all.cause - get(causename))
dat.df.agg.long <- dat.df.agg %>%
    pivot_longer(cols = c(all_of(causename), 
                          cause_other,
                          all.cause),
                 names_to = "cause",
                 values_to = "deaths") %>%
    filter(cause %in% c(causename, "cause_other")) %>%
    dplyr::select(-contains("cause_"))

# add year variables for each cause (NAs for other cause)
dat.df.agg.long$year_cause <- ifelse(dat.df.agg.long$cause == causename, 
                                     dat.df.agg.long$year, 
                                     NA)
dat.df.agg.long$year_other <- ifelse(dat.df.agg.long$cause != causename, 
                                     dat.df.agg.long$year, 
                                     NA)

# make cause a numeric vector, with cause 1 equal to the other causes 
#   - (so the coefficients from the binomial and Poisson models match)
dat.df.agg.long$cause_factor <- as.factor(ifelse(dat.df.agg.long$cause == causename, 2, 1))

f.ptrick <- "deaths ~ 1 + cause_factor"
f.ptrick <- paste(f.ptrick, 
                  "+ f(year_num, model = 'iid',
                   hyper = list(prec = list(initial = log(0.000001),
                                            fixed = TRUE)))")
f.ptrick <- paste(f.ptrick, 
                  "+ f(year_cause, model = 'iid',
                    hyper = list(prec = list(prior = 'pc.prec', 
                                             param = c(0.3, 0.01))))")

# linear combination without alphas for prediction
mat.pred.nonref <- Diagonal(nrow(dat.df.agg.long), 1)[which(dat.df.agg.long$cause != reference.cause), ]
marg.mat <- -1*model.matrix(~-1 + factor(year_num), data = dat.df.agg.long)[which(dat.df.agg.long$cause != reference.cause), ]
lcs <- inla.make.lincombs(
  "Predictor" = mat.pred.nonref,
  "year_num" = marg.mat
)

mod.ptrick <- inla(as.formula(f.ptrick),
                   control.fixed = fe.prec.pois,
                   family = "poisson",
                   data = dat.df.agg.long,
                   lincomb = lcs,
                   quantiles = res.quantiles,
                   control.predictor=list(compute=TRUE),
                   control.compute=list(config = TRUE))

res.pois <- cbind(dat.df.agg.long[dat.df.agg.long$cause != reference.cause, c("cause", "year", "year_num")],
                  expit(mod.ptrick$summary.lincomb.derived[, paste0(res.quantiles, "quant")]))
res.bin <- cbind(mod_data[, c("year", "year_num")], mod.bin$summary.fitted.values[, paste0(res.quantiles, "quant")])
res.bin <- cbind(cause = rep(causename, nrow(res.bin)), res.bin)
res.pois$model <- "Poisson"
res.bin$model <- "Binomial"
allres <- rbind(res.pois, res.bin)

# make table of FE results
summ_iid_fe <- rbind(mod.bin$summary.fixed,
                     mod.ptrick$summary.fixed)
summ_iid_fe <- cbind(model = c("Multinom", rep("", nrow(mod.bin$summary.fixed) - 1), "Pois", rep("", nrow(mod.ptrick$summary.fixed) - 1)), 
                     coefficient = c("$\\beta_0$","$\\tilde\\beta_0$", "$\\beta_0$"), 
                     summ_iid_fe)
rownames(summ_iid_fe) <- NULL
kable(summ_iid_fe, caption = "Posterior distributions of shared parameters from Multinomial and Poisson models: fixed intercept, IID RE on time",
      booktabs = TRUE, format = 'markdown', digits = 5) %>% kable_styling(latex_options = "HOLD_position")

# make table of RE hyperpar results
binhyper <- mod.bin$summary.hyperpar %>% dplyr::select(-sd, -mode)
poishyper <- mod.ptrick$summary.hyperpar %>% dplyr::select(-sd, -mode)
summ_iid_hyperprec <- rbind(binhyper,
                            poishyper)
summ_iid_hyperprec <- cbind(model = c("Multinom", "Pois"), 
                            hyperpar = c("$\\tau_{\\gamma}$", "$\\tau_{\\gamma}$"),
                            summ_iid_hyperprec)
rownames(summ_iid_hyperprec) <- NULL
kable(summ_iid_hyperprec, caption = "Posterior distributions of IID RE precision from Multinomial and Poisson models: fixed intercept, IID RE on time",
      booktabs = TRUE, format = 'markdown', digits = 5) %>% kable_styling(latex_options = "HOLD_position")

summ_iid_hypersd <- rbind(binhyper^-0.5,
                          poishyper^-0.5)
summ_iid_hypersd <- cbind(model = c("Multinom", "Pois"),  
                          hyperpar = c("$\\sigma_{\\gamma}$", "$\\sigma_{\\gamma}$"),
                          summ_iid_hypersd)
rownames(summ_iid_hypersd) <- NULL
summ_iid_hypersd <- summ_iid_hypersd[, c("model", "hyperpar", "mean", "0.975quant", "0.5quant", "0.025quant")]
names(summ_iid_hypersd) <- c("model", "hyperpar", "mean", "0.025quant", "0.5quant", "0.975quant")
kable(summ_iid_hypersd, caption = "Posterior distributions of IID RE SD from Multinomial and Poisson models: fixed intercept, IID RE on time",
      booktabs = TRUE, format = 'markdown', digits = 5) %>% kable_styling(latex_options = "HOLD_position")

binrand.df <- as.data.frame(mod.bin$summary.random$year)
binrand.df$model <- "Multinomial"
ptrickrand.df <- as.data.frame(mod.ptrick$summary.random$year_cause)
ptrickrand.df$model <- "Poisson"
summ_re <- data.frame(multinomial = binrand.df$`0.5quant`,
                      Poisson = ptrickrand.df$`0.5quant`)
kable(summ_re, caption = "Posterior medians of IID RE from Multinomial and Poisson models: fixed intercept, IID RE on time",
      booktabs = TRUE, format = 'markdown', digits = 5) %>% kable_styling(latex_options = "HOLD_position")
summ_re <- data.frame(multinomial = binrand.df$sd,
                      Poisson = ptrickrand.df$sd)
kable(summ_re, caption = "Posterior SD of IID RE from Multinomial and Poisson models: fixed intercept, IID RE on time",
      booktabs = TRUE, format = 'markdown', digits = 5) %>% kable_styling(latex_options = "HOLD_position")

# ggplot CSMFs
ggplot(data = allres %>% filter(cause == causename), 
       aes(x = year, y = `0.5quant`, color = model)) +
    geom_line(alpha = 0.2) +
    geom_ribbon(aes(x = year, ymin = `0.025quant`, ymax = `0.975quant`, fill = model), alpha = 0.2) +
    labs(x="year", 
         y="Cause fraction",
         subtitle = paste0("Estimated CSMF (",ci.level*100,"% posterior intervals) from yearly national model with IID RE on time"),
         title=unique(allChildCodesBroad$broadCodeName[allChildCodesBroad$cause == causename])) +
    theme_light()



## ----iidtime3, cache = TRUE------------------------------------------------------------------------------------------------
f.ptrick <- "deaths ~ 1 + cause_factor"
f.ptrick <- paste(f.ptrick, 
                  "+ f(year_num, model = 'iid',
                   hyper = list(prec = list(initial = log(0.000001),
                                            fixed = TRUE)))")
f.ptrick <- paste(f.ptrick, 
                  "+ f(year_cause, model = 'iid',
                    hyper = list(prec = list(prior = 'pc.prec', 
                                             param = c(0.3, 0.01))))")
f.ptrick <- paste(f.ptrick, 
                  "+ f(year_other, model = 'iid',
                    hyper = list(prec = list(prior = 'pc.prec', 
                                             param = c(0.3, 0.01))))")
mod.ptrick <- inla(as.formula(f.ptrick),
                   control.fixed = fe.prec.pois,
                   family = "poisson",
                   data = dat.df.agg.long,
                   lincomb = lcs,
                   quantiles = res.quantiles,
                   control.predictor=list(compute=TRUE),
                   control.compute=list(config = TRUE))

res.pois <- cbind(dat.df.agg.long[dat.df.agg.long$cause != reference.cause, c("cause", "year", "year_num")],
                  expit(mod.ptrick$summary.lincomb.derived[, paste0(res.quantiles, "quant")]))
res.pois$model <- "Poisson"
allres <- rbind(res.pois, res.bin)

# make table of FE results
summ_iid_fe <- rbind(mod.bin$summary.fixed,
                     mod.ptrick$summary.fixed)
summ_iid_fe <- cbind(model = c("Multinom", rep("", nrow(mod.bin$summary.fixed) - 1), "Pois", rep("", nrow(mod.ptrick$summary.fixed) - 1)), 
                     coefficient = c("$\\beta_0$","$\\tilde\\beta_0$", "$\\beta_0$"), 
                     summ_iid_fe)
rownames(summ_iid_fe) <- NULL
kable(summ_iid_fe, caption = "Posterior distributions of shared parameters from Multinomial and Poisson models: fixed intercept, IID RE on time",
      booktabs = TRUE, format = 'markdown', digits = 5) %>% kable_styling(latex_options = "HOLD_position")

# make table of RE hyperpar results
binhyper <- mod.bin$summary.hyperpar %>% dplyr::select(-sd, -mode)
poishyper <- mod.ptrick$summary.hyperpar %>% dplyr::select(-sd, -mode)
summ_iid_hyperprec <- rbind(binhyper,
                            poishyper)
summ_iid_hyperprec <- cbind(model = c("Multinom", "Pois", ""), 
                            hyperpar = c("$\\tau_{\\gamma}$", "$\\tau_{\\gamma}$", "$\\tilde\\tau_{\\gamma}$"),
                            summ_iid_hyperprec)
rownames(summ_iid_hyperprec) <- NULL
kable(summ_iid_hyperprec, caption = "Posterior distributions of IID RE precision from Multinomial and Poisson models: fixed intercept, IID RE on time",
      booktabs = TRUE, format = 'markdown', digits = 5) %>% kable_styling(latex_options = "HOLD_position")

summ_iid_hypersd <- rbind(binhyper^-0.5,
                          poishyper^-0.5)
summ_iid_hypersd <- cbind(model = c("Multinom", "Pois", ""),  
                          hyperpar = c("$\\sigma_{\\gamma}$", "$\\sigma_{\\gamma}$", "$\\tilde\\sigma_{\\gamma}$"),
                          summ_iid_hypersd)
summ_iid_hypersd <- summ_iid_hypersd[, c("model", "hyperpar", "mean", "0.975quant", "0.5quant", "0.025quant")]
names(summ_iid_hypersd) <- c("model", "hyperpar", "mean", "0.025quant", "0.5quant", "0.975quant")
rownames(summ_iid_hypersd) <- NULL
kable(summ_iid_hypersd, caption = "Posterior distributions of IID RE SD from Multinomial and Poisson models: fixed intercept, IID RE on time",
      booktabs = TRUE, format = 'markdown', digits = 5) %>% kable_styling(latex_options = "HOLD_position")

binrand.df <- as.data.frame(mod.bin$summary.random$year)
binrand.df$model <- "Multinomial"
ptrickrand.df <- as.data.frame(mod.ptrick$summary.random$year_cause)
ptrickrand.df$model <- "Poisson"
ptrickrandother.df <- as.data.frame(mod.ptrick$summary.random$year_other)
ptrickrandother.df$model <- "Poisson"
summ_re <- data.frame(multinomial = binrand.df$`0.5quant`,
                      Poisson_gamma = ptrickrand.df$`0.5quant`,
                      Poisson_gamma_tilde = ptrickrandother.df$`0.5quant`)
kable(summ_re, caption = "Posterior medians of IID REs from Multinomial and Poisson models: fixed intercept, IID RE on time",
      booktabs = TRUE, format = 'markdown', digits = 5,
      col.names = c("Multinomial $\\gamma_j$", "Poisson $\\gamma_j$", "Poisson $\\tilde\\gamma_j$")) %>% 
  kable_styling(latex_options = "HOLD_position")
summ_re <- data.frame(multinomial = binrand.df$sd,
                      Poisson_gamma = ptrickrand.df$sd,
                      Poisson_gammatilde = ptrickrandother.df$sd)
kable(summ_re, caption = "Posterior SD of IID REs from Multinomial and Poisson models: fixed intercept, IID RE on time",
      booktabs = TRUE, format = 'markdown', digits = 5,
      col.names = c("Multinomial $\\gamma_j$", "Poisson $\\gamma_j$", "Poisson $\\tilde\\gamma_j$")) %>% 
  kable_styling(latex_options = "HOLD_position")

# ggplot CSMFs
ggplot(data = allres %>% filter(cause == causename), 
       aes(x = year, y = `0.5quant`, color = model)) +
    geom_line(alpha = 0.2) +
    geom_ribbon(aes(x = year, ymin = `0.025quant`, ymax = `0.975quant`, fill = model), alpha = 0.2) +
    labs(x="year", 
         y="Cause fraction",
         subtitle = paste0("Estimated CSMF (",ci.level*100,"% posterior intervals) from yearly national model with IID RE on time, one for each cause"),
         title=unique(allChildCodesBroad$broadCodeName[allChildCodesBroad$cause == causename])) +
    theme_light()


## ----rw1, include = FALSE, cache = TRUE------------------------------------------------------------------------------------
#### 1st order RW on time ####
model_name <- "rw1year"

# binomial model formula
if (grepl("rw1year", model_name)) {
  f.bin <- paste(causename, "~ 1")
  f.bin <- paste(f.bin, 
                 "+ f(year_num, model = 'rw1', constr = TRUE,
                        hyper = list(prec = list(prior = 'pc.prec', 
                                                 param = c(0.3, 0.01))))")
  
}

# run binomial model
mod.bin <- inla(as.formula(f.bin),
                control.fixed = fe.prec.bin,
                data = mod_data,
                family = "binomial",
                Ntrials = all.cause,
                quantiles = res.quantiles,
                control.predictor=list(compute=TRUE),
                control.compute=list(config = TRUE))

# create long data for poisson trick model
dat.df.agg <- mod_data %>%
  mutate(cause_other = all.cause - get(causename))
dat.df.agg.long <- dat.df.agg %>%
  pivot_longer(cols = c(all_of(causename), 
                        cause_other,
                        all.cause),
               names_to = "cause",
               values_to = "deaths") %>%
  filter(cause %in% c(causename, "cause_other")) %>%
  dplyr::select(-contains("cause_"))

# add year variables for each cause (NAs for other cause)
dat.df.agg.long$year_cause <- ifelse(dat.df.agg.long$cause == causename, 
                                     dat.df.agg.long$year, 
                                     NA)
dat.df.agg.long$year_other <- ifelse(dat.df.agg.long$cause != causename, 
                                     dat.df.agg.long$year, 
                                     NA)

# make cause a numeric vector, with cause 1 equal to the other causes 
#   - (so the coefficients from the binomial and Poisson models match)
dat.df.agg.long$cause_factor <- as.factor(ifelse(dat.df.agg.long$cause == causename, 2, 1))

f.ptrick <- "deaths ~ 1 + cause_factor"
f.ptrick <- paste(f.ptrick,
                  "+ f(year_num, model = 'iid',
                     hyper = list(prec = list(initial = log(0.000001),
                                              fixed = TRUE)))")
f.ptrick <- paste(f.ptrick, 
                  "+ f(year_cause, model = 'rw1', constr = TRUE,
                      hyper = list(prec = list(prior = 'pc.prec', 
                                               param = c(0.3, 0.01))))")

mod.ptrick <- inla(as.formula(f.ptrick),
                   control.fixed = fe.prec.pois,
                   family = "poisson",
                   data = dat.df.agg.long,
                   lincomb = lcs,
                   quantiles = res.quantiles,
                   control.predictor=list(compute=TRUE),
                   control.compute=list(config = TRUE))

res.pois <- cbind(dat.df.agg.long[dat.df.agg.long$cause != reference.cause, c("cause", "year", "year_num")],
                  expit(mod.ptrick$summary.lincomb.derived[, paste0(res.quantiles, "quant")]))
res.bin <- cbind(mod_data[, c("year", "year_num")], mod.bin$summary.fitted.values[, paste0(res.quantiles, "quant")])
res.bin <- cbind(cause = rep(causename, nrow(res.bin)), res.bin)
res.pois$model <- "Poisson"
res.bin$model <- "Binomial"
allres <- rbind(res.pois, res.bin)


## ----rw1tab1---------------------------------------------------------------------------------------------------------------
# make table of FE results
summ_iid_fe <- rbind(mod.bin$summary.fixed["(Intercept)",],
                     mod.ptrick$summary.fixed["cause_factor2",])
summ_iid_fe <- cbind(coefficient = c("$\\beta_0$",""), 
                     model = c("Multinom", "Pois"), 
                     summ_iid_fe)
rownames(summ_iid_fe) <- NULL
kable(summ_iid_fe, caption = "Posterior distributions of shared parameters from Multinomial and Poisson models: fixed intercept, RW1 on time",
      booktabs = TRUE, format = 'markdown', digits = 5) %>% kable_styling(latex_options = "HOLD_position")

# make table of RE hyperpar results
binhyper <- mod.bin$summary.hyperpar %>% dplyr::select(-sd, -mode)
poishyper <- mod.ptrick$summary.hyperpar %>% dplyr::select(-sd, -mode)
summ_iid_hyperprec <- rbind(binhyper,
                            poishyper)
summ_iid_hyperprec <- cbind(model = c("Multinom", "Pois"), 
                            hyperpar = c("$\\tau_{\\gamma}$", "$\\tau_{\\gamma}$"),
                            summ_iid_hyperprec)
rownames(summ_iid_hyperprec) <- NULL
kable(summ_iid_hyperprec, caption = "Posterior distributions of IID RE precision from Multinomial and Poisson models: fixed intercept, RW1 on time",
      booktabs = TRUE, format = 'markdown', digits = 5) %>% kable_styling(latex_options = "HOLD_position")

summ_iid_hypersd <- rbind(binhyper^-0.5,
                          poishyper^-0.5)
summ_iid_hypersd <- cbind(model = c("Multinom", "Pois"),  
                          hyperpar = c("$\\sigma_{\\gamma}$", "$\\sigma_{\\gamma}$"),
                          summ_iid_hypersd)
summ_iid_hypersd <- summ_iid_hypersd[, c("model", "hyperpar", "mean", "0.975quant", "0.5quant", "0.025quant")]
names(summ_iid_hypersd) <- c("model", "hyperpar", "mean", "0.025quant", "0.5quant", "0.975quant")
rownames(summ_iid_hypersd) <- NULL
kable(summ_iid_hypersd, caption = "Posterior distributions of IID RE SD from Multinomial and Poisson models: fixed intercept, RW1 on time",
      booktabs = TRUE, format = 'markdown', digits = 5) %>% kable_styling(latex_options = "HOLD_position")

binrand.df <- as.data.frame(mod.bin$summary.random$year)
binrand.df$model <- "Multinomial"
ptrickrand.df <- as.data.frame(mod.ptrick$summary.random$year_cause)
ptrickrand.df$model <- "Poisson"
summ_re <- data.frame(multinomial = binrand.df$`0.5quant`,
                      Poisson = ptrickrand.df$`0.5quant`)
kable(summ_re, caption = "Posterior medians of IID RE from Multinomial and Poisson models: fixed intercept, RW1 on time",
      booktabs = TRUE, format = 'markdown', digits = 5) %>% kable_styling(latex_options = "HOLD_position")
summ_re <- data.frame(multinomial = binrand.df$sd,
                      Poisson = ptrickrand.df$sd)
kable(summ_re, caption = "Posterior SD of IID RE from Multinomial and Poisson models: fixed intercept, RW1 on time",
      booktabs = TRUE, format = 'markdown', digits = 5) %>% kable_styling(latex_options = "HOLD_position")

# ggplot CSMFs
ggplot(data = allres %>% filter(cause == causename), 
       aes(x = year, y = `0.5quant`, color = model)) +
  geom_line(alpha = 0.2) +
  geom_ribbon(aes(x = year, ymin = `0.025quant`, ymax = `0.975quant`, fill = model), alpha = 0.2) +
  labs(x="year", 
       y="Cause fraction",
       subtitle = paste0("Estimated CSMF (",ci.level*100,"% posterior intervals) from yearly national model with RW1 on time"),
       title=unique(allChildCodesBroad$broadCodeName[allChildCodesBroad$cause == causename])) +
  theme_light()


## ----rw2slopemodels, include = FALSE, cache = TRUE-------------------------------------------------------------------------
## 2nd order RW
f.bin <- paste(causename, "~ 1 + year_cent")
f.bin <- paste(f.bin, 
               "+ f(year_num, model = 'rw2',
                    constr = TRUE,
                    hyper = list(prec = list(prior = 'pc.prec', 
                                             param = c(0.3, 0.01))))")

# run binomial model
mod.bin <- inla(as.formula(f.bin),
                control.fixed = fe.prec.bin,
                data = mod_data,
                family = "binomial",
                Ntrials = all.cause,
                quantiles = res.quantiles,
                control.predictor=list(compute=TRUE),
                control.compute=list(config = TRUE))

# add year variables for each cause (NAs for other cause)
dat.df.agg.long$year_cause <- ifelse(dat.df.agg.long$cause == causename, 
                                     dat.df.agg.long$year, 
                                     NA)
dat.df.agg.long$year_other <- ifelse(dat.df.agg.long$cause != causename, 
                                     dat.df.agg.long$year, 
                                     NA)
dat.df.agg.long$year2_cause <- ifelse(dat.df.agg.long$cause == causename, 
                                      dat.df.agg.long$year2, 
                                      NA)
dat.df.agg.long$year2_other <- ifelse(dat.df.agg.long$cause != causename, 
                                      dat.df.agg.long$year2, 
                                      NA)
table(dat.df.agg.long$year_cause, dat.df.agg.long$cause)
table(dat.df.agg.long$year_other, dat.df.agg.long$cause)

# Poisson model formula
f.ptrick <- "deaths ~ 1 + cause_factor*year_cent"
f.ptrick <- paste(f.ptrick, 
                  "+ f(year_num, model = 'iid',
                   hyper = list(prec = list(initial = log(0.000001),
                                            fixed = TRUE)))")
f.ptrick <- paste(f.ptrick, 
                  "+ f(year_cause, model = 'rw2',
                    constr = TRUE,
                    hyper = list(prec = list(prior = 'pc.prec', 
                                             param = c(0.3, 0.01))))")

mod.ptrick <- inla(as.formula(f.ptrick),
                   control.fixed = fe.prec.pois,
                   family = "poisson",
                   data = dat.df.agg.long,
                   lincomb = lcs,
                   quantiles = res.quantiles,
                   control.predictor=list(compute=TRUE),
                   control.compute=list(config = TRUE))

res.pois <- cbind(dat.df.agg.long[dat.df.agg.long$cause != reference.cause, c("cause", "year", "year_num")],
                  expit(mod.ptrick$summary.lincomb.derived[, paste0(res.quantiles, "quant")]))
res.bin <- cbind(mod_data[, c("year", "year_num")], mod.bin$summary.fitted.values[, paste0(res.quantiles, "quant")])
res.bin <- cbind(cause = rep(causename, nrow(res.bin)), res.bin)
res.pois$model <- "Poisson"
res.bin$model <- "Binomial"
allres <- rbind(res.pois, res.bin)


## ----rw2slopetabs----------------------------------------------------------------------------------------------------------
# make table of FE results
summ_rw2_fe <- rbind(mod.bin$summary.fixed,
                     mod.ptrick$summary.fixed)
summ_rw2_fe <- cbind(model = c("Multinom", rep("", nrow(mod.bin$summary.fixed) - 1), "Pois", rep("", nrow(mod.ptrick$summary.fixed) - 1)), 
                     coefficient = c("$\\beta_0$","$\\beta_1$","$\\tilde\\beta_0$", "$\\beta_0$", "$\\tilde\\beta_1$", "$\\beta_1$"), 
                     summ_rw2_fe)
rownames(summ_rw2_fe) <- NULL
kable(summ_rw2_fe, caption = "Posterior distributions of shared parameters from Multinomial and Poisson models: fixed intercept, linear slope on time, RW2",
      booktabs = TRUE, format = 'markdown', digits = 5) %>% kable_styling(latex_options = "HOLD_position")

# make table of RE hyperpar results
binhyper <- mod.bin$summary.hyperpar %>% dplyr::select(-sd, -mode)
poishyper <- mod.ptrick$summary.hyperpar %>% dplyr::select(-sd, -mode)
summ_rw2_hyperprec <- rbind(binhyper,
                            poishyper)
summ_rw2_hyperprec <- cbind(model = c("Multinom", "Pois"), 
                            hyperpar = c("$\\tau_{\\gamma}$", "$\\tau_{\\gamma}$"),
                            summ_rw2_hyperprec)
rownames(summ_rw2_hyperprec) <- NULL
kable(summ_rw2_hyperprec, caption = "Posterior distributions of RW2 precision from Multinomial and Poisson models: fixed intercept, linear slope on time, RW2 on time", booktabs = TRUE, format = 'markdown', digits = 5) %>% kable_styling(latex_options = "HOLD_position")

summ_rw2_hypersd <- rbind(binhyper^-0.5,
                          poishyper^-0.5)
summ_rw2_hypersd <- cbind(model = c("Multinom", "Pois"),  
                          hyperpar = c("$\\sigma_{\\gamma}$", "$\\sigma_{\\gamma}$"),
                          summ_rw2_hypersd)
summ_rw2_hypersd <- summ_rw2_hypersd[, c("model", "hyperpar", "mean", "0.975quant", "0.5quant", "0.025quant")]
names(summ_rw2_hypersd) <- c("model", "hyperpar", "mean", "0.025quant", "0.5quant", "0.975quant")
rownames(summ_rw2_hypersd) <- NULL
kable(summ_rw2_hypersd, caption = "Posterior distributions of RW2 SD from Multinomial and Poisson models: fixed intercept, linear slope on time, RW2 on time",
      booktabs = TRUE, format = 'markdown', digits = 5) %>% kable_styling(latex_options = "HOLD_position")


binrand.df <- as.data.frame(mod.bin$summary.random$year)
binrand.df$model <- "Multinomial"
ptrickrand.df <- as.data.frame(mod.ptrick$summary.random$year_cause)
ptrickrand.df$model <- "Poisson"
summ_re <- data.frame(multinomial = binrand.df$`0.5quant`,
                      Poisson = ptrickrand.df$`0.5quant`)
kable(summ_re, caption = "Posterior medians of RW2 RE from Multinomial and Poisson models: fixed intercept, linear slope, RW2 on time",
      booktabs = TRUE, format = 'markdown', digits = 5) %>% kable_styling(latex_options = "HOLD_position")
summ_re <- data.frame(multinomial = binrand.df$sd,
                      Poisson = ptrickrand.df$sd)
kable(summ_re, caption = "Posterior SD of RW2 RE from Multinomial and Poisson models: fixed intercept, linear slope, RW2 on time",
      booktabs = TRUE, format = 'markdown', digits = 5) %>% kable_styling(latex_options = "HOLD_position")



## ----rw2slopecsmfs---------------------------------------------------------------------------------------------------------
# ggplot CSMFs
ggplot(data = allres %>% filter(cause == causename), 
       aes(x = year, y = `0.5quant`, color = model)) +
    geom_line(alpha = 0.2) +
    geom_ribbon(aes(x = year, ymin = `0.025quant`, ymax = `0.975quant`, fill = model), alpha = 0.2) +
    labs(x="year", 
         y="Cause fraction",
         subtitle = paste0("Estimated CSMF (",ci.level*100,"% posterior intervals) from yearly national model with linear slope and RW2 on time"),
         title=unique(allChildCodesBroad$broadCodeName[allChildCodesBroad$cause == causename])) +
    theme_light()



## ----spatial_bym2, echo = FALSE, cache = TRUE------------------------------------------------------------------------------
# set aggregation dimension
aggdim <- "state"

# set data to use
mod_data <- dat.df.agg.list[[aggdim]]

# format data
if (aggdim == "state") {
  mod_data$state_code2 <- mod_data$state_code
}

#### IID in time ####
model_name <- "bym2"

# binomial model formula
if (grepl("bym2", model_name)) {
  f.bin <- paste(causename, "~ 1")
  f.bin <- paste(f.bin, 
                 "+ f(state_code, model = 'bym2',
                                  graph = india.adj, 
                                  scale.model = T, 
                                  constr = T,
                                  hyper=list(phi=list(prior='pc', 
                                                      param=c(0.5, 0.5), 
                                                      initial=1),
                                             prec=list(prior='pc.prec', 
                                                       param=c(0.3,0.01), 
                                                       initial=5)))")
  
}

# run binomial model
mod.bin <- inla(as.formula(f.bin),
                control.fixed = fe.prec.bin,
                data = mod_data,
                family = "binomial",
                Ntrials = all.cause,
                quantiles = res.quantiles,
                control.predictor=list(compute=TRUE),
                control.compute=list(config = TRUE))

# create long data for poisson trick model
dat.df.agg <- mod_data %>%
  mutate(cause_other = all.cause - get(causename))
dat.df.agg.long <- dat.df.agg %>%
  pivot_longer(cols = c(all_of(causename), 
                        cause_other,
                        all.cause),
               names_to = "cause",
               values_to = "deaths") %>%
  filter(cause %in% c(causename, "cause_other")) %>%
  dplyr::select(-contains("cause_"))

# add state_code variables for each cause (NAs for other cause)
dat.df.agg.long$state_code_cause <- ifelse(dat.df.agg.long$cause == causename, 
                                           dat.df.agg.long$state_code, 
                                           NA)
dat.df.agg.long$state_code_other <- ifelse(dat.df.agg.long$cause != causename, 
                                           dat.df.agg.long$state_code, 
                                           NA)
# table(dat.df.agg.long$state_code_cause, dat.df.agg.long$cause)
# table(dat.df.agg.long$state_code_other, dat.df.agg.long$cause)
# make cause a numeric vector, with cause 1 equal to the other causes 
#   - (so the coefficients from the binomial and Poisson models match)
dat.df.agg.long$cause_factor <- as.factor(ifelse(dat.df.agg.long$cause == causename, 2, 1))

if (grepl("bym2", model_name)) {
  f.ptrick <- "deaths ~ 1 + cause_factor"
  f.ptrick <- paste(f.ptrick,
                    "+ f(state_code, model = 'iid',
                       hyper = list(prec = list(initial = log(0.000001),
                                                fixed = TRUE)))")
  f.ptrick <- paste(f.ptrick, 
                    "+ f(state_code_cause, model = 'bym2',
                                     graph=india.adj, 
                                     scale.model=T, 
                                     constr=T,
                                     hyper=list(phi=list(prior='pc', 
                                                         param=c(0.5, 0.5), 
                                                         initial=1),
                                                prec=list(prior='pc.prec', 
                                                          param=c(0.3,0.01), 
                                                          initial=5)))")
}

mat.pred.nonref <- Diagonal(nrow(dat.df.agg.long), 1)[which(dat.df.agg.long$cause != reference.cause), ]
marg.mat <- -1*model.matrix(~-1 + factor(state_code), data = dat.df.agg.long)[which(dat.df.agg.long$cause != reference.cause), ]
lcs <- inla.make.lincombs(
  "Predictor" = mat.pred.nonref,
  "state_code" = marg.mat
)

mod.ptrick <- inla(as.formula(f.ptrick),
                   control.fixed = fe.prec.pois,
                   family = "poisson",
                   data = dat.df.agg.long,
                   lincomb = lcs,
                   quantiles = res.quantiles,
                   control.predictor=list(compute=TRUE),
                   control.compute=list(config = TRUE))

res.pois <- cbind(dat.df.agg.long[dat.df.agg.long$cause != reference.cause, c("cause", "state", "state_code")],
                  expit(mod.ptrick$summary.lincomb.derived[, paste0(res.quantiles, "quant")]))
res.bin <- cbind(mod_data[, c("state" ,"state_code")], mod.bin$summary.fitted.values[, paste0(res.quantiles, "quant")])
res.bin <- cbind(cause = rep(causename, nrow(res.bin)), res.bin)
res.pois$model <- "Poisson"
res.bin$model <- "Binomial"
allres <- rbind(res.pois, res.bin)


## ----spatial_bym2_results, fig.width=8, fig.height=10----------------------------------------------------------------------
# make table of FE results
summ_iid_fe <- rbind(mod.bin$summary.fixed["(Intercept)",],
                     mod.ptrick$summary.fixed["cause_factor2",])
summ_iid_fe <- cbind(coefficient = c("$\\beta_0$","$\\beta_0$"), 
                     model = c("Multinom", "Pois"), 
                     summ_iid_fe)
rownames(summ_iid_fe) <- NULL
kable(summ_iid_fe, caption = "Posterior distributions of shared parameters from Multinomial and Poisson models: fixed intercept, BYM2",
      booktabs = TRUE, format = 'markdown', digits = 5) %>% kable_styling(latex_options = "HOLD_position")

# make table of RE hyperpar results
binhyper <- mod.bin$summary.hyperpar %>% dplyr::select(-sd, -mode)
poishyper <- mod.ptrick$summary.hyperpar %>% dplyr::select(-sd, -mode)
summ_iid_hyperprec <- rbind(binhyper,
                            poishyper)
summ_iid_hyperprec <- cbind(model = c("Multinom", "", "Pois", ""), 
                            hyperpar = c("$\\tau$", "$\\phi$", "$\\tau$", "$\\phi$"),
                            summ_iid_hyperprec)
rownames(summ_iid_hyperprec) <- NULL
kable(summ_iid_hyperprec, caption = "Posterior distributions of BYM2 hyperparameters from Multinomial and Poisson models: fixed intercept, BYM2",
      booktabs = TRUE, format = 'markdown', digits = 5) %>% kable_styling(latex_options = "HOLD_position")

binrand.df <- as.data.frame(mod.bin$summary.random$state_code)
binrand.df$model <- "Multinomial"
ptrickrand.df <- as.data.frame(mod.ptrick$summary.random$state_code_cause)
ptrickrand.df$model <- "Poisson"
summ_re <- data.frame(multinomial = binrand.df$`0.5quant`,
                      Poisson = ptrickrand.df$`0.5quant`)
kable(summ_re, caption = "Posterior medians of BYM2 REs from Multinomial and Poisson models: fixed intercept, BYM2",
      booktabs = TRUE, format = 'markdown', digits = 5) %>% 
  kable_styling(latex_options = "HOLD_position")
summ_re <- data.frame(multinomial = binrand.df$sd,
                      Poisson = ptrickrand.df$sd)
kable(summ_re, caption = "Posterior SD of BYM2 REs from Multinomial and Poisson models: fixed intercept, BYM2",
      booktabs = TRUE, format = 'markdown', digits = 5) %>% 
  kable_styling(latex_options = "HOLD_position")

# ggplot CSMFs
ggplot(data = allres %>% filter(cause == causename), 
       aes(x = model, y = `0.5quant`, color = model)) +
  geom_errorbar(aes(ymin = `0.025quant`, ymax = `0.975quant`), width = .1) +
  geom_point() +
  facet_geo(~ state, grid = india_grid_MDS, scales = "free_y") +
  labs(x="model", 
       y="Cause fraction",
       subtitle = paste0("Estimated CSMF (",ci.level*100,"% posterior intervals) from regional model with BYM2 REs"),
       title=unique(allChildCodesBroad$broadCodeName[allChildCodesBroad$cause == causename])) +
  theme_light()



## ----spatiotemporal--------------------------------------------------------------------------------------------------------
### Spatiotemporal modeling
aggdim_index <- 4

# set aggregation dimension
aggdim <- aggregationDims[aggdim_index]

# set data to use
mod_data <- dat.df.agg.list[[aggdim]]

#### BYM2 in space ####
model_name <- "bym2rw2type1"
# model_name <- "bym2rw2type4"

## Binomial model ##
f.bin <- paste(causename,"~ 1")

# if BYM2 spatial term
if (grepl("bym2", model_name)) {
  f.bin <- paste(f.bin, "+ f(state_code, model = 'bym2',
                                     graph=india.adj, 
                                     scale.model=T, 
                                     constr=T,
                                     hyper=list(phi=list(prior='pc', 
                                                         param=c(0.5, 0.5)),
                                                prec=list(prior='pc.prec', 
                                                          param=c(0.3,0.01))))")
}

# if RW2 temporal trend
if (grepl("rw2", model_name)) {
  f.bin <- paste(f.bin, 
              "+ f(year, model = 'rw2',
                           constr=TRUE,
                           hyper=list(prec=list(prior='pc.prec', 
                                                param=c(0.3,0.01)))) +
                      f(year2, model = 'iid',constr = TRUE,
                        hyper = list(prec = list(prior = 'pc.prec',
                                                 param = c(0.5, 0.01))))")
}

# interaction terms (no term for no interaction)
if (grepl("type1", model_name)  & !grepl("typeIV", model_name)) {
  f.bin <- paste(f.bin,
              "+ f(state.year, model = 'iid',
                           hyper = list(prec = list(prior = 'pc.prec',
                                                    param = c(0.5, 0.01))))")
} else if (grepl("type4", model_name)) {
  f.bin <- paste(f.bin,
              "+ f(state.int, model = 'besag', graph = india.adj, constr = TRUE,
                           group = year.int, control.group = list(model = 'rw2'),
                           hyper = list(prec = list(prior = 'pc.prec',
                                                    param = c(0.5, 0.01))))")
}


# run binomial model
mod.bin <- inla(as.formula(f.bin),
                control.fixed = fe.prec.bin,
                data = mod_data,
                family = "binomial",
                Ntrials = all.cause,
                quantiles = res.quantiles,
                control.predictor=list(compute=TRUE),
                control.compute=list(config = TRUE))

# create long data for poisson trick model
dat.df.agg <- mod_data %>%
  mutate(cause_other = all.cause - get(causename))
dat.df.agg.long <- dat.df.agg %>%
  pivot_longer(cols = c(all_of(causename), 
                        cause_other,
                        all.cause),
               names_to = "cause",
               values_to = "deaths") %>%
  filter(cause %in% c(causename, "cause_other")) %>%
  dplyr::select(-contains("cause_"))

# create an indicator for each cause (with NAs for each other cause) in order to fit interactions for Poisson trick
for (varname in c("year", "year2", "state_code", "state.year", "state.int", "year.int")) {
  vname_cause <- paste(varname, "cause", sep = "_")
  vname_other <- paste(varname, "other", sep = "_")
  dat.df.agg.long[[vname_cause]] <- dat.df.agg.long[[varname]]
  dat.df.agg.long[dat.df.agg.long$cause != causename, vname_cause] <- NA
  dat.df.agg.long[[vname_other]] <- dat.df.agg.long[[varname]]
  dat.df.agg.long[dat.df.agg.long$cause == causename, vname_other] <- NA
}

# make cause a numeric vector, with cause 1 equal to the other causes 
#   - (so the coefficients from the binomial and Poisson models match)
dat.df.agg.long$cause_factor <- as.factor(ifelse(dat.df.agg.long$cause == causename, 2, 1))
dat.df.agg.long$cause_num <- as.numeric(dat.df.agg.long$cause_factor) - 1

# model formula
f.ptrick <- "deaths ~ cause_factor +
          f(state.year, model='iid',
            hyper = list(
              prec = list(
                initial = log(0.000001),
                fixed = TRUE)))"

# if BYM2 spatial term
if (grepl("bym2", model_name)) {
  f.ptrick <- paste(f.ptrick, "+ 
                           f(state_code_cause, model = 'bym2',
                             graph=india.adj, 
                             scale.model=T, 
                             constr=T,
                             hyper=list(phi=list(prior='pc', 
                                                 param=c(0.5, 0.5)),
                                        prec=list(prior='pc.prec', 
                                                  param=c(0.3,0.01))))")
}

# if RW2 temporal trend
if (grepl("rw2", model_name)) {
  f.ptrick <- paste(f.ptrick,"+ 
                          f(year_cause, model = 'rw2',
                            constr = TRUE,
                            hyper = list(prec = list(prior = 'pc.prec', 
                                                     param = c(0.3,0.01)))) +
                          f(year2_cause, model = 'iid',
                            constr = TRUE,
                            hyper = list(prec = list(prior = 'pc.prec',
                                                     param = c(0.5, 0.01))))")
}

# interaction terms (no term for no interaction)
if (grepl("type1", model_name)) {
  f.ptrick <- paste(f.ptrick, "+ 
                           f(state.year_cause, model = 'iid',
                             hyper = list(prec = list(prior = 'pc.prec',
                                                      param = c(0.5, 0.01))))")
} else if (grepl("type4", model_name)) {
  f.ptrick <- paste(f.ptrick, "+ 
                           f(state.int_cause, model = 'besag', graph = india.adj, constr = TRUE,
                             group = year.int_cause, control.group = list(model = 'rw2'),
                             hyper = list(prec = list(prior = 'pc.prec',
                                                      param = c(0.5, 0.01))))")
}

# linear combination without alphas for prediction
mat.pred.nonref <- Diagonal(nrow(dat.df.agg.long), 1)[which(dat.df.agg.long$cause != reference.cause), ]
marg.mat <- -1*model.matrix(~-1 + factor(state.year), data = dat.df.agg.long)[which(dat.df.agg.long$cause != reference.cause), ]
lcs <- inla.make.lincombs(
  "Predictor" = mat.pred.nonref,
  "state.year" = marg.mat
)

mod.ptrick <- inla(as.formula(f.ptrick),
                   control.fixed = fe.prec.pois,
                   family = "poisson",
                   data = dat.df.agg.long,
                   lincomb = lcs,
                   quantiles = res.quantiles,
                   control.predictor=list(compute=TRUE),
                   control.compute=list(config = TRUE))

res.pois <- cbind(dat.df.agg.long[dat.df.agg.long$cause != reference.cause, c("cause", "state", "state_code", "year")],
                  expit(mod.ptrick$summary.lincomb.derived[, paste0(res.quantiles, "quant")]))
res.bin <- cbind(mod_data[, c("state" ,"state_code", "year")], mod.bin$summary.fitted.values[, paste0(res.quantiles, "quant")])
res.bin <- cbind(cause = rep(causename, nrow(res.bin)), res.bin)
res.pois$model <- "Poisson"
res.bin$model <- "Binomial"
allres <- rbind(res.pois, res.bin)


## ----spatiotemporalres, fig.height=10.5, fig.width=8-----------------------------------------------------------------------
# make table of FE results
summ_iid_fe <- rbind(mod.bin$summary.fixed["(Intercept)",],
                     mod.ptrick$summary.fixed["cause_factor2",])
summ_iid_fe <- cbind(coefficient = c("$\\beta_0$","$\\beta_0$"), 
                     model = c("Multinom", "Pois"), 
                     summ_iid_fe)
rownames(summ_iid_fe) <- NULL
kable(summ_iid_fe, caption = "Posterior distributions of shared parameters from Multinomial and Poisson models: Type I",
      booktabs = TRUE, format = 'markdown', digits = 5) %>% kable_styling(latex_options = "HOLD_position")

# make table of RE hyperpar results
binhyper <- mod.bin$summary.hyperpar %>% dplyr::select(-sd, -mode)
poishyper <- mod.ptrick$summary.hyperpar %>% dplyr::select(-sd, -mode)
summ_iid_hyperprec <- rbind(binhyper,
                            poishyper)
summ_iid_hyperprec <- cbind(model = c("Multinom", "", "", "", "", "Pois", "", "", "", ""),
                            parameter = rep(c("$\\tau_{BYM}$", "$\\phi$", "$\\tau_{\\gamma}$", "$\\tau_{\\psi}$", "$\\tau_{\\delta}$"), 2),
                            summ_iid_hyperprec)
rownames(summ_iid_hyperprec) <- NULL
kable(summ_iid_hyperprec, caption = "Posterior distributions of BYM2 hyperparameters from Multinomial and Poisson models: Type I",
      booktabs = TRUE, format = 'markdown', digits = 5) %>% kable_styling(latex_options = "HOLD_position")

# ggplot CSMFs
ggplot(data = allres %>% filter(cause == causename), 
       aes(x = year, y = `0.5quant`, fill = model, color = model, linetype = model)) +
  geom_ribbon(aes(ymin = `0.025quant`, ymax = `0.975quant`), alpha = 0.2) +
  geom_line() +
  facet_geo(~ state, grid = india_grid_MDS, scales = "free_y") +
  labs(x="model", 
       y="Cause fraction",
       subtitle = paste0("Estimated CSMF (",ci.level*100,"% posterior intervals) from state-year spatiotemporal model with type I interaction"),
       title=unique(allChildCodesBroad$broadCodeName[allChildCodesBroad$cause == causename])) +
  theme_light()


## ----spatiotemporal2-------------------------------------------------------------------------------------------------------
### Spatiotemporal modeling
aggdim_index <- 4

# set aggregation dimension
aggdim <- aggregationDims[aggdim_index]

# set data to use
mod_data <- dat.df.agg.list[[aggdim]]

#### BYM2 in space ####
# model_name <- "bym2rw2type1"
model_name <- "bym2rw2type4"

# format data
mod_data$year.int <- as.integer(mod_data$year.int - min(mod_data$year.int) + 1)
mod_data$state.year <- as.integer(as.factor(mod_data$state.year))

## Binomial model ##
f.bin <- paste(causename,"~ 1")

# if BYM2 spatial term
if (grepl("bym2", model_name)) {
  f.bin <- paste(f.bin, "+ f(state_code, model = 'bym2',
                                     graph=india.adj, 
                                     scale.model=T, 
                                     constr=T,
                                     hyper=list(phi=list(prior='pc', 
                                                         param=c(0.5, 0.5)),
                                                prec=list(prior='pc.prec', 
                                                          param=c(0.3,0.01))))")
}

# if RW2 temporal trend
if (grepl("rw2", model_name)) {
  f.bin <- paste(f.bin, 
              "+ f(year, model = 'rw2',
                           constr=TRUE,
                           hyper=list(prec=list(prior='pc.prec', 
                                                param=c(0.3,0.01)))) +
                      f(year2, model = 'iid',constr = TRUE,
                        hyper = list(prec = list(prior = 'pc.prec',
                                                 param = c(0.5, 0.01))))")
}

# interaction terms (no term for no interaction)
if (grepl("type1", model_name)  & !grepl("typeIV", model_name)) {
  f.bin <- paste(f.bin,
              "+ f(state.year, model = 'iid',
                           hyper = list(prec = list(prior = 'pc.prec',
                                                    param = c(0.5, 0.01))))")
} else if (grepl("type4", model_name)) {
  f.bin <- paste(f.bin,
              "+ f(state.int, model = 'besag', graph = india.adj, constr = TRUE,
                           group = year.int, control.group = list(model = 'rw2'),
                           hyper = list(prec = list(prior = 'pc.prec',
                                                    param = c(0.5, 0.01))))")
}


# run binomial model
mod.bin <- inla(as.formula(f.bin),
                control.fixed = fe.prec.bin,
                data = mod_data,
                family = "binomial",
                Ntrials = all.cause,
                quantiles = res.quantiles,
                control.predictor=list(compute=TRUE),
                control.compute=list(config = TRUE))

# create long data for poisson trick model
dat.df.agg <- mod_data %>%
  mutate(cause_other = all.cause - get(causename))
dat.df.agg.long <- dat.df.agg %>%
  pivot_longer(cols = c(all_of(causename), 
                        cause_other,
                        all.cause),
               names_to = "cause",
               values_to = "deaths") %>%
  filter(cause %in% c(causename, "cause_other")) %>%
  dplyr::select(-contains("cause_"))

# create an indicator for each cause (with NAs for each other cause) in order to fit interactions for Poisson trick
for (varname in c("year", "year2", "state_code", "state.year", "state.int", "year.int")) {
  vname_cause <- paste(varname, "cause", sep = "_")
  vname_other <- paste(varname, "other", sep = "_")
  dat.df.agg.long[[vname_cause]] <- dat.df.agg.long[[varname]]
  dat.df.agg.long[dat.df.agg.long$cause != causename, vname_cause] <- NA
  dat.df.agg.long[[vname_other]] <- dat.df.agg.long[[varname]]
  dat.df.agg.long[dat.df.agg.long$cause == causename, vname_other] <- NA
}

# make cause a numeric vector, with cause 1 equal to the other causes 
#   - (so the coefficients from the binomial and Poisson models match)
dat.df.agg.long$cause_factor <- as.factor(ifelse(dat.df.agg.long$cause == causename, 2, 1))
dat.df.agg.long$cause_num <- as.numeric(dat.df.agg.long$cause_factor) - 1

# model formula
f.ptrick <- "deaths ~ cause_factor +
          f(state.year, model='iid',
            hyper = list(
              prec = list(
                initial = log(0.000001),
                fixed = TRUE)))"

# if BYM2 spatial term
if (grepl("bym2", model_name)) {
  f.ptrick <- paste(f.ptrick, "+ 
                           f(state_code_cause, model = 'bym2',
                             graph=india.adj, 
                             scale.model=T, 
                             constr=T,
                             hyper=list(phi=list(prior='pc', 
                                                 param=c(0.5, 0.5)),
                                        prec=list(prior='pc.prec', 
                                                  param=c(0.3,0.01))))")
}

# if RW2 temporal trend
if (grepl("rw2", model_name)) {
  f.ptrick <- paste(f.ptrick,"+ 
                          f(year_cause, model = 'rw2',
                            constr = TRUE,
                            hyper = list(prec = list(prior = 'pc.prec', 
                                                     param = c(0.3,0.01)))) +
                          f(year2_cause, model = 'iid',
                            constr = TRUE,
                            hyper = list(prec = list(prior = 'pc.prec',
                                                     param = c(0.5, 0.01))))")
}

# interaction terms (no term for no interaction)
if (grepl("type1", model_name)) {
  f.ptrick <- paste(f.ptrick, "+ 
                           f(state.year_cause, model = 'iid',
                             hyper = list(prec = list(prior = 'pc.prec',
                                                      param = c(0.5, 0.01))))")
} else if (grepl("type4", model_name)) {
  f.ptrick <- paste(f.ptrick, "+ 
                           f(state.int_cause, model = 'besag', graph = india.adj, constr = TRUE,
                             group = year.int_cause, control.group = list(model = 'rw2'),
                             hyper = list(prec = list(prior = 'pc.prec',
                                                      param = c(0.5, 0.01))))")
}

# linear combination without alphas for prediction
mat.pred.nonref <- Diagonal(nrow(dat.df.agg.long), 1)[which(dat.df.agg.long$cause != reference.cause), ]
marg.mat <- -1*model.matrix(~-1 + factor(state.year), data = dat.df.agg.long)[which(dat.df.agg.long$cause != reference.cause), ]
lcs <- inla.make.lincombs(
  "Predictor" = mat.pred.nonref,
  "state.year" = marg.mat
)

mod.ptrick <- inla(as.formula(f.ptrick),
                   control.fixed = fe.prec.pois,
                   family = "poisson",
                   data = dat.df.agg.long,
                   lincomb = lcs,
                   quantiles = res.quantiles,
                   control.predictor=list(compute=TRUE),
                   control.compute=list(config = TRUE))

res.pois <- cbind(dat.df.agg.long[dat.df.agg.long$cause != reference.cause, c("cause", "state", "state_code", "year")],
                  expit(mod.ptrick$summary.lincomb.derived[, paste0(res.quantiles, "quant")]))
res.bin <- cbind(mod_data[, c("state" ,"state_code", "year")], mod.bin$summary.fitted.values[, paste0(res.quantiles, "quant")])
res.bin <- cbind(cause = rep(causename, nrow(res.bin)), res.bin)
res.pois$model <- "Poisson"
res.bin$model <- "Binomial"
allres <- rbind(res.pois, res.bin)


## ----spatiotemporalres2, fig.height=10.5, fig.width=8----------------------------------------------------------------------
# make table of FE results
summ_iid_fe <- rbind(mod.bin$summary.fixed["(Intercept)",],
                     mod.ptrick$summary.fixed["cause_factor2",])
summ_iid_fe <- cbind(coefficient = c("$\\beta_0$","$\\beta_0$"), 
                     model = c("Multinom", "Pois"), 
                     summ_iid_fe)
rownames(summ_iid_fe) <- NULL
kable(summ_iid_fe, caption = "Posterior distributions of shared parameters from Multinomial and Poisson models: Type I",
      booktabs = TRUE, format = 'markdown', digits = 5) %>% kable_styling(latex_options = "HOLD_position")

# make table of RE hyperpar results
binhyper <- mod.bin$summary.hyperpar %>% dplyr::select(-sd, -mode)
poishyper <- mod.ptrick$summary.hyperpar %>% dplyr::select(-sd, -mode)
summ_iid_hyperprec <- rbind(binhyper,
                            poishyper)
summ_iid_hyperprec <- cbind(model = c("Multinom", "", "", "", "", "Pois", "", "", "", ""),
                            parameter = rep(c("$\\tau_{BYM}$", "$\\phi$", "$\\tau_{\\gamma}$", "$\\tau_{\\psi}$", "$\\tau_{\\delta}$"), 2),
                            summ_iid_hyperprec)
rownames(summ_iid_hyperprec) <- NULL
kable(summ_iid_hyperprec, caption = "Posterior distributions of BYM2 hyperparameters from Multinomial and Poisson models: Type I",
      booktabs = TRUE, format = 'markdown', digits = 5) %>% kable_styling(latex_options = "HOLD_position")

# ggplot CSMFs
ggplot(data = allres %>% filter(cause == causename), 
       aes(x = year, y = `0.5quant`, fill = model, color = model, linetype = model)) +
  geom_ribbon(aes(ymin = `0.025quant`, ymax = `0.975quant`), alpha = 0.2) +
  geom_line() +
  facet_geo(~ state, grid = india_grid_MDS, scales = "free_y") +
  labs(x="model", 
       y="Cause fraction",
       subtitle = paste0("Estimated CSMF (",ci.level*100,"% posterior intervals) from state-year spatiotemporal model with type I interaction"),
       title=unique(allChildCodesBroad$broadCodeName[allChildCodesBroad$cause == causename])) +
  theme_light()

