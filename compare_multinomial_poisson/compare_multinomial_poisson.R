# Austin Schumacher
# 6/12/2020
# Compare binomial and poisson trick models for one cause
#   - building up from intercept to spatial only to temporal only to spatiotemporal

#### clear env and load packages ####
rm(list=ls())
library(SUMMER);
library(milliondeaths); library(tidyverse);
library(INLA);library(mapmisc);library(spdep);
library(ggpubr); library(broom); library(geofacet);
library(knitr); library(kableExtra);

set.seed(530)

#### functions used in this code ####

# function for posterior samples in order to get posterior dist of CSMF from Poisson trick model
getPtrickCSMFposteriorSummary <- function(model, data, causevar, idvar, nsamples, ncause, ci.level) {
    
    # testing
    # model <- mod.ptrick
    # data <- dat.df.agg.long
    # causevar <- "cause"
    # idvar <- "year"
    # nsamples <- 100
    # ncause <- 2
    
    x.samples <- inla.posterior.sample(nsamples, model, use.improved.mean = TRUE, seed = 80085)
    x.contents <- model$misc$configs$contents
    effect <- "Predictor"
    id.effect <- which(x.contents$tag==effect)
    ind.effects <- x.contents$start[id.effect]-1 + (1:x.contents$length[id.effect])
    data <- data %>% dplyr::select(all_of(c(causevar, idvar)))
    results <- matrix(NA, nrow = nrow(data), ncol = nsamples)
    for (i in 1:nsamples) {
        tmp <- data
        tmp$samp <- x.samples[[i]]$latent[ind.effects,]
        csmf.sample <- tmp %>% 
            group_by(get(idvar)) %>%
            mutate(csmf = exp(samp)/sum(exp(samp))) %>%
            ungroup()
        results[, i] <- csmf.sample$csmf
    }
    ci.end <- ci.level
    r1 <- as.data.frame(t(apply(results, 1, quantile, c(0.5, (1-ci.level)/2, 1- (1-ci.level)/2))))
    names(r1) <- c("csmf_med", "csmf_lower", "csmf_upper")
    r1 <- cbind(data, r1)
    return(r1)
}

# function for posterior samples in order to get posterior dist of CSMF from binomial model
getBinCSMFposteriorSummary <- function(model, data, causevar, idvar, nsamples, ci.level) {
    
    # testing
    # model <- mod.bin
    # data <- mod_data
    # causevar <- causename
    # idvar <- "year"
    # nsamples <- 100
    
    x.samples <- inla.posterior.sample(nsamples, model, use.improved.mean = TRUE, seed = 80085)
    x.contents <- model$misc$configs$contents
    effect <- "Predictor"
    id.effect <- which(x.contents$tag==effect)
    ind.effects <- x.contents$start[id.effect]-1 + (1:x.contents$length[id.effect])
    data <- data %>% 
        dplyr::select(c(all_of(idvar)))
    results <- matrix(NA, nrow = nrow(data)*2, ncol = nsamples)
    for (i in 1:nsamples) {
        tmp <- data
        tmp$csmf_cause <- expit(x.samples[[i]]$latent[ind.effects,])
        tmp$csmf_other <- 1 - tmp$csmf_cause
        csmf.sample <- tmp %>% 
            pivot_longer(cols = contains("csmf_"),
                         names_to = "cause",
                         values_to = "csmf",
                         names_prefix = "csmf_")
        results[, i] <- csmf.sample$csmf
    }
    r1 <- as.data.frame(t(apply(results, 1, quantile, c(0.5, (1-ci.level)/2, 1- (1-ci.level)/2))))
    names(r1) <- c("csmf_med", "csmf_lower", "csmf_upper")
    data$csmf_cause <- NA
    data$csmf_other <- NA
    data2 <- data %>% 
        pivot_longer(cols = contains("csmf_"),
                     names_to = "cause",
                     values_to = "csmf",
                     names_prefix = "csmf_")
    data2$cause <- dplyr::recode(data2$cause, cause = causevar, other = "cause_other")
    data2$csmf <- NULL
    r1 <- cbind(data2, r1)
    return(r1)
}

# function to run binomial and Poisson models
# runBinPois <- function(binFormula, poisFormula, model_name, causename, idvar) {
#     
# }

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

# UI level
ci.level <- 0.8

# set fixed effects precisions
fe.prec <- list(prec.intercept = 0,
                prec = 0)


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

### aggdim loop goes here ###

# testing
aggdim_index <- 3

# set aggregation dimension
aggdim <- aggregationDims[aggdim_index]

# set data to use
mod_data <- dat.df.agg.list[[aggdim]]

#######TEMPORARY##########
## add more data in
# mod_data_add <- mod_data
# mod_data_add$year <- mod_data_add$year + 10
# mod_data_add$year2 <- mod_data_add$year2 + 10
# mod_data_add[, 3:10] <- mod_data_add[, 3:10] + rpois(8,90)
# mod_data_add$all.cause <- rowSums(mod_data_add[, 3:10])
# mod_data <- rbind(mod_data, mod_data_add)

# format data
if (aggdim == "year") {
    mod_data$year_num <- mod_data$year - min(mod_data$year) + 1
    mod_data$year_cent <- mod_data$year - mean(mod_data$year)
    mod_data$year2_num <- mod_data$year2 - min(mod_data$year2) + 1
    mod_data$year2_cent <- mod_data$year2 - mean(mod_data$year2)
}

### model loop goes here ###

# testing
mod_dim_index <- 1

# model 
# model_name <- aggDimModels[[aggdim]][mod_dim_index]

### if statement for aggdim goes here ###

### cause loop goes here ###

# testing
causename_index <- 3

# set cause name
causename <- causenames_v[causename_index]

#### linear slope only ####

# set model name
model_name <- "none"

## TRY collapsing data to a single point, fitting models
# mod_data <- mod_data %>% summarize(cause_1C01 = sum(cause_1C01),
#                                    all.cause = sum(all.cause))
# mod_data$id <- 1

# makean id var for binomial model if the agg dim is "none"
if (aggdim == "none") mod_data$id <- 1

# binomial model formula
f.bin <- paste(causename, "~ 1")

if (grepl("linear", model_name)) {
    f.bin <- paste(f.bin, "+ year_cent")
}
    
# run binomial model
mod.bin <- inla(as.formula(f.bin),
                control.fixed = fe.prec,
                data = mod_data,
                family = "binomial",
                Ntrials = all.cause,
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
dat.df.agg.long$id <- 1

# make cause a numeric vector, with cause 1 equal to the other causes 
#   - (so the coefficients from the binomial and Poisson models match)
dat.df.agg.long$cause_factor <- as.factor(ifelse(dat.df.agg.long$cause == causename, 2, 1))
dat.df.agg.long$cause_num <- as.numeric(dat.df.agg.long$cause_factor) - 1

# run poisson trick model intercept-only
# f.ptrick.re.fixedprec <- "deaths ~ -1 + 
#                             f(cause, model='iid',
#                               constr = TRUE,
#                               hyper = list(prec = list(initial = log(0.0001),
#                                                        fixed = TRUE))) +
#                             f(id, model = 'iid', 
#                               hyper = list(prec = list(initial = log(0.000001),
#                                                        fixed = TRUE)))"

f.ptrick.fe <- "deaths ~ 1 + cause_factor"

if (grepl("none", model_name)) {
    f.ptrick.fe <- paste(f.ptrick.fe,
                         "+ f(id, model = 'iid', constr = TRUE,
                              hyper = list(prec = list(initial = log(0.000001),
                                                       fixed = TRUE)))")
}

if (grepl("linear", model_name)) {
    f.ptrick.fe <- paste(f.ptrick.fe, 
                         "* year_cent")
    f.ptrick.fe <- paste(f.ptrick.fe,
                         "+ f(year_num, model = 'iid', constr = TRUE,
                              hyper = list(prec = list(initial = log(0.000001),
                                                       fixed = TRUE)))")
}

mod.ptrick <- inla(as.formula(f.ptrick.fe),
                   control.fixed = fe.prec,
                   family = "poisson",
                   data = dat.df.agg.long,
                   control.predictor=list(compute=TRUE),
                   control.compute=list(config = TRUE))

## compare results
if (grepl(model_name, "none")) {
  summ <- rbind(mod.bin$summary.fixed,
                mod.ptrick$summary.fixed)
  summ <- cbind(model = c("Multinom", rep("", nrow(mod.bin$summary.fixed) - 1), "Pois", rep("", nrow(mod.ptrick$summary.fixed) - 1)), 
                coefficient = c("$\\beta_0$","$\\tilde\\beta_0$", "$\\beta_0$"), 
                summ)
  rownames(summ) <- NULL
  kable(summ, caption = "Posterior distributions of shared parameters from Multinomial and Poisson models: fixed intercept only",
        booktabs = TRUE, digits = 5) %>% kable_styling(latex_options = "HOLD_position")
}
if (grepl(model_name, "linear")) {
  summ_linear <- rbind(mod.bin$summary.fixed,
                       mod.ptrick$summary.fixed)
  summ_linear <- cbind(model = c("Multinom", rep("", nrow(mod.bin$summary.fixed) - 1), "Pois", rep("", nrow(mod.ptrick$summary.fixed) - 1)), 
                       coefficient = c("$\\beta_0$","$\\beta_1$","$\\tilde\\beta_0$", "$\\beta_0$", "$\\tilde\\beta_1$", "$\\beta_1$"), 
                       summ_linear)
  rownames(summ_linear) <- NULL
  kable(summ_linear, caption = "Posterior distributions of shared parameters from Multinomial and Poisson models: fixed intercept, linear slope on time",
        booktabs = TRUE, digits = 5) %>% kable_styling(latex_options = "HOLD_position")
}

# posterior samples of CSMFs
csmf_posterior_bin <- getBinCSMFposteriorSummary(model = mod.bin, 
                                                 data = mod_data, 
                                                 causevar = causename, 
                                                 idvar = "id", 
                                                 nsamples = 1000,
                                                 ci.level = ci.level)
csmf_posterior_ptrick <- getPtrickCSMFposteriorSummary(model = mod.ptrick, 
                                                       data = dat.df.agg.long, 
                                                       causevar = "cause", 
                                                       idvar = "id", 
                                                       nsamples = 1000, 
                                                       ncause = 2,
                                                       ci.level = ci.level)
csmf_posterior_bin$model <- "multinomial"
csmf_posterior_ptrick$model <- "Poisson"
allres <- rbind(csmf_posterior_bin, csmf_posterior_ptrick)

# ggplot CSMFs
ggplot(data = allres %>% filter(cause == causename), 
       aes(x = id, y = csmf_med, color = model)) +
    geom_line(alpha = 0.2) +
    geom_ribbon(aes(x = id, ymin = csmf_lower, ymax = csmf_upper, fill = model), alpha = 0.2) +
    labs(x="year", 
         y="Cause fraction",
         subtitle = paste0("Estimated CSMF (", ci.level*100, "% posterior intervals) from yearly national model with linear slope on time"),
         title=unique(allChildCodesBroad$broadCodeName[allChildCodesBroad$cause == causename])) +
    theme_light()

#### IID in time ####
model_name <- "iidyear"

# binomial model formula
if (grepl("iidyear", model_name)) {
    f.bin <- paste(causename, "~ 1")
    f.bin <- paste(f.bin, 
                   "+ f(year, model = 'iid',
                        hyper = list(prec = list(prior = 'pc.prec', 
                                                 param = c(0.3, 0.01))))")
    
}

# run binomial model
mod.bin <- inla(as.formula(f.bin),
                control.fixed = fe.prec,
                data = mod_data,
                family = "binomial",
                Ntrials = all.cause,
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

if (grepl("iidyear", model_name)) {
    f.ptrick <- "deaths ~ 1 + cause_factor"
    f.ptrick <- paste(f.ptrick,
                      "+ f(year_num, model = 'iid', constr = TRUE,
                       hyper = list(prec = list(initial = log(0.000001),
                                                fixed = TRUE)))")
    f.ptrick <- paste(f.ptrick, 
                      "+ f(year_cause, model = 'iid',
                        hyper = list(prec = list(prior = 'pc.prec', 
                                                 param = c(0.3, 0.01))))")
    
    # f.ptrick <- paste(f.ptrick,
    #                   "+ f(year_other, model = 'iid',
    #                     hyper = list(prec = list(prior = 'pc.prec',
    #                                              param = c(0.3, 0.01))))")
}

mod.ptrick <- inla(as.formula(f.ptrick),
                   control.fixed = fe.prec,
                   family = "poisson",
                   data = dat.df.agg.long,
                   control.predictor=list(compute=TRUE),
                   control.compute=list(config = TRUE))

# posterior samples of CSMFs
csmf_posterior_bin <- getBinCSMFposteriorSummary(model = mod.bin, 
                                                 data = mod_data, 
                                                 causevar = causename, 
                                                 idvar = "year", 
                                                 nsamples = 1000,
                                                 ci.level = ci.level)
csmf_posterior_ptrick <- getPtrickCSMFposteriorSummary(model = mod.ptrick, 
                                                       data = dat.df.agg.long, 
                                                       causevar = "cause", 
                                                       idvar = "year", 
                                                       nsamples = 1000, 
                                                       ncause = 2,
                                                       ci.level = ci.level)
csmf_posterior_bin$model <- "multinomial"
csmf_posterior_ptrick$model <- "Poisson"
allres <- rbind(csmf_posterior_bin, csmf_posterior_ptrick)

# make table of FE results
summ_iid_fe <- rbind(mod.bin$summary.fixed["(Intercept)",],
                     mod.ptrick$summary.fixed["cause_factor2",])
summ_iid_fe <- cbind(coefficient = c("$\\beta_0$",""), 
                     model = c("Multinom", "Pois"), 
                     summ_iid_fe)
rownames(summ_iid_fe) <- NULL
kable(summ_iid_fe, caption = "Posterior distributions of shared parameters from Multinomial and Poisson models: fixed intercept, IID RE on time",
      booktabs = TRUE, digits = 5) %>% kable_styling(latex_options = "HOLD_position")

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
      booktabs = TRUE, digits = 5) %>% kable_styling(latex_options = "HOLD_position")

summ_iid_hypersd <- rbind(binhyper^-0.5,
                          poishyper^-0.5)
summ_iid_hypersd <- cbind(model = c("Multinom", "Pois"),  
                          hyperpar = c("$\\sigma_{\\gamma}$", "$\\sigma_{\\gamma}$"),
                          summ_iid_hypersd)
summ_iid_hypersd <- summ_iid_hypersd[, c("model", "hyperpar", "mean", "0.975quant", "0.5quant", "0.025quant")]
names(summ_iid_hypersd) <- c("model", "hyperpar", "mean", "0.025quant", "0.5quant", "0.975quant")
rownames(summ_iid_hypersd) <- NULL
kable(summ_iid_hypersd, caption = "Posterior distributions of IID RE SD from Multinomial and Poisson models: fixed intercept, IID RE on time",
      booktabs = TRUE, digits = 5) %>% kable_styling(latex_options = "HOLD_position")

binrand.df <- as.data.frame(mod.bin$summary.random$year)
binrand.df$model <- "Multinomial"
ptrickrand.df <- as.data.frame(mod.ptrick$summary.random$year_cause)
ptrickrand.df$model <- "Poisson"
summ_re <- data.frame(multinomial = binrand.df$`0.5quant`,
                      Poisson = ptrickrand.df$`0.5quant`)
kable(summ_re, caption = "Posterior medians of IID RE from Multinomial and Poisson models: fixed intercept, IID RE on time",
      booktabs = TRUE, format = 'markdown', digits = 5) %>% 
  kable_styling(latex_options = "HOLD_position")
summ_re <- data.frame(multinomial = binrand.df$sd,
                      Poisson = ptrickrand.df$sd)
kable(summ_re, caption = "Posterior SD of IID RE from Multinomial and Poisson models: fixed intercept, IID RE on time",
      booktabs = TRUE, format = 'markdown', digits = 5) %>% 
  kable_styling(latex_options = "HOLD_position")

# ggplot CSMFs
ggplot(data = allres %>% filter(cause == causename), 
       aes(x = year, y = csmf_med, color = model)) +
    geom_line(alpha = 0.2) +
    geom_ribbon(aes(x = year, ymin = csmf_lower, ymax = csmf_upper, fill = model), alpha = 0.2) +
    labs(x="year", 
         y="Cause fraction",
         subtitle = paste0("Estimated CSMF (",ci.level*100,"% posterior intervals) from yearly national model with linear slope on time"),
         title=unique(allChildCodesBroad$broadCodeName[allChildCodesBroad$cause == causename])) +
    theme_light()

#### 1st order RW on time ####
model_name <- "rw1year"

# binomial model formula
if (grepl("rw1year", model_name)) {
  f.bin <- paste(causename, "~ 1")
  f.bin <- paste(f.bin, 
                 "+ f(year, model = 'rw1', constr = TRUE,
                        hyper = list(prec = list(prior = 'pc.prec', 
                                                 param = c(0.3, 0.01))))")
  
}

# run binomial model
mod.bin <- inla(as.formula(f.bin),
                control.fixed = fe.prec,
                data = mod_data,
                family = "binomial",
                Ntrials = all.cause,
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

if (grepl("rw1year", model_name)) {
  f.ptrick <- "deaths ~ 1 + cause_factor"
  f.ptrick <- paste(f.ptrick,
                    "+ f(year_num, model = 'iid', constr = TRUE,
                       hyper = list(prec = list(initial = log(0.000001),
                                                fixed = TRUE)))")
  f.ptrick <- paste(f.ptrick, 
                    "+ f(year_cause, model = 'rw1', constr = TRUE,
                        hyper = list(prec = list(prior = 'pc.prec', 
                                                 param = c(0.3, 0.01))))")
  
  # f.ptrick <- paste(f.ptrick,
  #                   "+ f(year_other, model = 'rw1', constr = TRUE,
  #                     hyper = list(prec = list(prior = 'pc.prec',
  #                                              param = c(0.3, 0.01))))")
}

mod.ptrick <- inla(as.formula(f.ptrick),
                   control.fixed = fe.prec,
                   family = "poisson",
                   data = dat.df.agg.long,
                   control.predictor=list(compute=TRUE),
                   control.compute=list(config = TRUE))

# posterior samples of CSMFs
csmf_posterior_bin <- getBinCSMFposteriorSummary(model = mod.bin, 
                                                 data = mod_data, 
                                                 causevar = causename, 
                                                 idvar = "year", 
                                                 nsamples = 1000,
                                                 ci.level = ci.level)
csmf_posterior_ptrick <- getPtrickCSMFposteriorSummary(model = mod.ptrick, 
                                                       data = dat.df.agg.long, 
                                                       causevar = "cause", 
                                                       idvar = "year", 
                                                       nsamples = 1000, 
                                                       ncause = 2,
                                                       ci.level = ci.level)
csmf_posterior_bin$model <- "multinomial"
csmf_posterior_ptrick$model <- "Poisson"
allres <- rbind(csmf_posterior_bin, csmf_posterior_ptrick)

# make table of FE results
summ_iid_fe <- rbind(mod.bin$summary.fixed["(Intercept)",],
                     mod.ptrick$summary.fixed["cause_factor2",])
summ_iid_fe <- cbind(coefficient = c("$\\beta_0$",""), 
                     model = c("Multinom", "Pois"), 
                     summ_iid_fe)
rownames(summ_iid_fe) <- NULL
kable(summ_iid_fe, caption = "Posterior distributions of shared parameters from Multinomial and Poisson models: fixed intercept, RW1 on time",
      booktabs = TRUE, digits = 5) %>% kable_styling(latex_options = "HOLD_position")

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
      booktabs = TRUE, digits = 5) %>% kable_styling(latex_options = "HOLD_position")

summ_iid_hypersd <- rbind(binhyper^-0.5,
                          poishyper^-0.5)
summ_iid_hypersd <- cbind(model = c("Multinom", "Pois"),  
                          hyperpar = c("$\\sigma_{\\gamma}$", "$\\sigma_{\\gamma}$"),
                          summ_iid_hypersd)
summ_iid_hypersd <- summ_iid_hypersd[, c("model", "hyperpar", "mean", "0.975quant", "0.5quant", "0.025quant")]
names(summ_iid_hypersd) <- c("model", "hyperpar", "mean", "0.025quant", "0.5quant", "0.975quant")
rownames(summ_iid_hypersd) <- NULL
kable(summ_iid_hypersd, caption = "Posterior distributions of IID RE SD from Multinomial and Poisson models: fixed intercept, RW1 on time",
      booktabs = TRUE, digits = 5) %>% kable_styling(latex_options = "HOLD_position")

binrand.df <- as.data.frame(mod.bin$summary.random$year)
binrand.df$model <- "Multinomial"
ptrickrand.df <- as.data.frame(mod.ptrick$summary.random$year_cause)
ptrickrand.df$model <- "Poisson"
summ_re <- data.frame(multinomial = binrand.df$`0.5quant`,
                      Poisson = ptrickrand.df$`0.5quant`)
kable(summ_re, caption = "Posterior medians of IID RE from Multinomial and Poisson models: fixed intercept, RW1 on time",
      booktabs = TRUE, digits = 5) %>% kable_styling(latex_options = "HOLD_position")
summ_re <- data.frame(multinomial = binrand.df$sd,
                      Poisson = ptrickrand.df$sd)
kable(summ_re, caption = "Posterior SD of IID RE from Multinomial and Poisson models: fixed intercept, RW1 on time",
      booktabs = TRUE, digits = 5) %>% kable_styling(latex_options = "HOLD_position")

# ggplot CSMFs
ggplot(data = allres %>% filter(cause == causename), 
       aes(x = year, y = csmf_med, color = model)) +
  geom_line(alpha = 0.2) +
  geom_ribbon(aes(x = year, ymin = csmf_lower, ymax = csmf_upper, fill = model), alpha = 0.2) +
  labs(x="year", 
       y="Cause fraction",
       subtitle = paste0("Estimated CSMF (",ci.level*100,"% posterior intervals) from yearly national model with RW1 on time"),
       title=unique(allChildCodesBroad$broadCodeName[allChildCodesBroad$cause == causename])) +
  theme_light()

#### 2nd order RW ####

model_name <- "rw2slope"
# model_name <- "rw2iid"

# binomial model formula
if (grepl("rw2iid", model_name)) {
    f.bin <- paste(causename, "~ 1")
    f.bin <- paste(f.bin, 
                   "+ f(year, model = 'rw2',
                        constr = TRUE,
                        hyper = list(prec = list(prior = 'pc.prec', 
                                                 param = c(0.3, 0.01))))    +
        f(year2, model = 'iid',
          hyper = list(prec = list(prior = 'pc.prec',
                                   param = c(0.5, 0.01))))")

} else if (grepl("rw2slope", model_name)) {
    f.bin <- paste(causename, "~ 1 + year_cent")
    f.bin <- paste(f.bin, 
                   "+ f(year, model = 'rw2',
                        constr = TRUE,
                        hyper = list(prec = list(prior = 'pc.prec', 
                                                 param = c(0.3, 0.01))))")
}

# run binomial model
mod.bin <- inla(as.formula(f.bin),
                control.fixed = fe.prec,
                data = mod_data,
                family = "binomial",
                Ntrials = all.cause,
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
if (grepl("rw2iid", model_name)) {
    f.ptrick <- "deaths ~ 1 + cause_factor"
    f.ptrick <- paste(f.ptrick, 
                      "+ f(year_num, model = 'iid',
                       hyper = list(prec = list(initial = log(0.000001),
                                                fixed = TRUE)))")
    f.ptrick <- paste(f.ptrick, 
                   "+ f(year_cause, model = 'rw2',
                        constr = TRUE,
                        hyper = list(prec = list(prior = 'pc.prec', 
                                                 param = c(0.3, 0.01))))     +
        f(year2_cause, model = 'iid',
          hyper = list(prec = list(prior = 'pc.prec',
                                   param = c(0.5, 0.01))))")

    f.ptrick <- paste(f.ptrick, 
                   "+ f(year_other, model = 'rw2',
                        constr = TRUE,
                        hyper = list(prec = list(prior = 'pc.prec', 
                                                 param = c(0.3, 0.01))))     +
        f(year2_other, model = 'iid',
          hyper = list(prec = list(prior = 'pc.prec',
                                   param = c(0.5, 0.01))))")
} else if (grepl("rw2slope", model_name)) {
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
    # f.ptrick <- paste(f.ptrick, 
    #                   "+ f(year_other, model = 'rw2',
    #                     constr = TRUE,
    #                     hyper = list(prec = list(prior = 'pc.prec', 
    #                                              param = c(0.3, 0.01))))")
}


mod.ptrick <- inla(as.formula(f.ptrick),
                   control.fixed = fe.prec,
                   family = "poisson",
                   data = dat.df.agg.long,
                   control.predictor=list(compute=TRUE),
                   control.compute=list(config = TRUE))

# posterior samples of CSMFs
csmf_posterior_bin <- getBinCSMFposteriorSummary(model = mod.bin, 
                                                 data = mod_data, 
                                                 causevar = causename, 
                                                 idvar = "year", 
                                                 nsamples = 1000,
                                                 ci.level = ci.level)
csmf_posterior_ptrick <- getPtrickCSMFposteriorSummary(model = mod.ptrick, 
                                                       data = dat.df.agg.long, 
                                                       causevar = "cause", 
                                                       idvar = "year", 
                                                       nsamples = 1000, 
                                                       ncause = 2,
                                                       ci.level = ci.level)
csmf_posterior_bin$model <- "multinomial"
csmf_posterior_ptrick$model <- "Poisson"
allres <- rbind(csmf_posterior_bin, csmf_posterior_ptrick)

# make table of FE results
summ_rw2_fe <- rbind(mod.bin$summary.fixed,
                     mod.ptrick$summary.fixed)
summ_rw2_fe <- cbind(model = c("Multinom", rep("", nrow(mod.bin$summary.fixed) - 1), "Pois", rep("", nrow(mod.ptrick$summary.fixed) - 1)), 
                     coefficient = c("$\\beta_0$","$\\beta_1$","$\\tilde\\beta_0$", "$\\beta_0$", "$\\tilde\\beta_1$", "$\\beta_1$"), 
                     summ_rw2_fe)
rownames(summ_rw2_fe) <- NULL
kable(summ_rw2_fe, caption = "Posterior distributions of shared parameters from Multinomial and Poisson models: fixed intercept, linear slope on time, RW2",
      booktabs = TRUE, digits = 5) %>% kable_styling(latex_options = "HOLD_position")

# make table of RE hyperpar results
binhyper <- mod.bin$summary.hyperpar %>% dplyr::select(-sd, -mode)
poishyper <- mod.ptrick$summary.hyperpar %>% dplyr::select(-sd, -mode)
summ_rw2_hyperprec <- rbind(binhyper,
                            poishyper)
summ_rw2_hyperprec <- cbind(model = c("Multinom", "Pois"), 
                            hyperpar = c("$\\tau_1$", "$\\tau_1$"),
                            summ_rw2_hyperprec)
rownames(summ_rw2_hyperprec) <- NULL
kable(summ_rw2_hyperprec, caption = "Posterior distributions of IID RE precision from Multinomial and Poisson models: fixed intercept, linear slope on time, RW2 on time",
      booktabs = TRUE, digits = 5) %>% kable_styling(latex_options = "HOLD_position")

summ_rw2_hypersd <- rbind(binhyper^-0.5,
                          poishyper^-0.5)
summ_rw2_hypersd <- cbind(model = c("Multinom", "Pois"),  
                          hyperpar = c("$\\sigma_1$", "$\\sigma_1$"),
                          summ_rw2_hypersd)
summ_rw2_hypersd <- summ_rw2_hypersd[, c("model", "hyperpar", "mean", "0.975quant", "0.5quant", "0.025quant")]
names(summ_rw2_hypersd) <- c("model", "hyperpar", "mean", "0.025quant", "0.5quant", "0.975quant")
rownames(summ_rw2_hypersd) <- NULL
kable(summ_rw2_hypersd, caption = "Posterior distributions of IID RE SD from Multinomial and Poisson models: fixed intercept, linear slope on time, RW2 on time",
      booktabs = TRUE, digits = 5) %>% kable_styling(latex_options = "HOLD_position")

# ggplot CSMFs
ggplot(data = allres %>% filter(cause == causename), 
       aes(x = year, y = csmf_med, color = model)) +
    geom_line(alpha = 0.2) +
    geom_ribbon(aes(x = year, ymin = csmf_lower, ymax = csmf_upper, fill = model), alpha = 0.2) +
    labs(x="year", 
         y="Cause fraction",
         subtitle = "Estimated CSMF (95% posterior intervals) from yearly national model with linear slope on time",
         title=unique(allChildCodesBroad$broadCodeName[allChildCodesBroad$cause == causename])) +
    theme_light()


##### Spatial modeling #####

aggdim_index <- 2

# set aggregation dimension
aggdim <- aggregationDims[aggdim_index]

# set data to use
mod_data <- dat.df.agg.list[[aggdim]]

#######TEMPORARY##########
## add more data in
# mod_data_add <- mod_data
# mod_data_add$year <- mod_data_add$year + 10
# mod_data_add$year2 <- mod_data_add$year2 + 10
# mod_data_add[, 3:10] <- mod_data_add[, 3:10] + rpois(8,90)
# mod_data_add$all.cause <- rowSums(mod_data_add[, 3:10])
# mod_data <- rbind(mod_data, mod_data_add)

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
                                  graph='india.adj', 
                                  scale.model=T, 
                                  constr=T,
                                  hyper=list(phi=list(prior='pc', 
                                                      param=c(0.5, 0.5), 
                                                      initial=1),
                                             prec=list(prior='pc.prec', 
                                                       param=c(0.3,0.01), 
                                                       initial=5)))")
  
}

# run binomial model
mod.bin <- inla(as.formula(f.bin),
                control.fixed = fe.prec,
                data = mod_data,
                family = "binomial",
                Ntrials = all.cause,
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
table(dat.df.agg.long$state_code_cause, dat.df.agg.long$cause)
table(dat.df.agg.long$state_code_other, dat.df.agg.long$cause)
# make cause a numeric vector, with cause 1 equal to the other causes 
#   - (so the coefficients from the binomial and Poisson models match)
dat.df.agg.long$cause_factor <- as.factor(ifelse(dat.df.agg.long$cause == causename, 2, 1))

if (grepl("bym2", model_name)) {
  f.ptrick <- "deaths ~ 1 + cause_factor"
  f.ptrick <- paste(f.ptrick,
                    "+ f(state_code, model = 'iid', constr = TRUE,
                       hyper = list(prec = list(initial = log(0.000001),
                                                fixed = TRUE)))")
  f.ptrick <- paste(f.ptrick, 
                    "+ f(state_code_cause, model = 'bym2',
                                     graph='india.adj', 
                                     scale.model=T, 
                                     constr=T,
                                     hyper=list(phi=list(prior='pc', 
                                                         param=c(0.5, 0.5), 
                                                         initial=1),
                                                prec=list(prior='pc.prec', 
                                                          param=c(0.3,0.01), 
                                                          initial=5)))")
}

mod.ptrick <- inla(as.formula(f.ptrick),
                   control.fixed = fe.prec,
                   family = "poisson",
                   data = dat.df.agg.long,
                   control.predictor=list(compute=TRUE),
                   control.compute=list(config = TRUE))

# posterior samples of CSMFs
csmf_posterior_bin <- getBinCSMFposteriorSummary(model = mod.bin, 
                                                 data = mod_data, 
                                                 causevar = causename, 
                                                 idvar = "state", 
                                                 nsamples = 1000,
                                                 ci.level = ci.level)
csmf_posterior_ptrick <- getPtrickCSMFposteriorSummary(model = mod.ptrick, 
                                                       data = dat.df.agg.long, 
                                                       causevar = "cause", 
                                                       idvar = "state", 
                                                       nsamples = 1000, 
                                                       ncause = 2,
                                                       ci.level = ci.level)
csmf_posterior_bin$model <- "multinomial"
csmf_posterior_ptrick$model <- "Poisson"
allres <- rbind(csmf_posterior_bin, csmf_posterior_ptrick)

# make table of FE results
summ_iid_fe <- rbind(mod.bin$summary.fixed["(Intercept)",],
                     mod.ptrick$summary.fixed["cause_factor2",])
summ_iid_fe <- cbind(coefficient = c("$\\beta_0$","$\\beta_0$"), 
                     model = c("Multinom", "Pois"), 
                     summ_iid_fe)
rownames(summ_iid_fe) <- NULL
kable(summ_iid_fe, caption = "Posterior distributions of shared parameters from Multinomial and Poisson models: fixed intercept, BYM2",
      booktabs = TRUE, digits = 5) %>% kable_styling(latex_options = "HOLD_position")

# make table of RE hyperpar results
binhyper <- mod.bin$summary.hyperpar %>% dplyr::select(-sd, -mode)
poishyper <- mod.ptrick$summary.hyperpar %>% dplyr::select(-sd, -mode)
summ_iid_hyperprec <- rbind(binhyper,
                            poishyper)
summ_iid_hyperprec <- cbind(model = c("Multinom", "", "Pois", ""), 
                            hyperpar = c("$\\tau_{\\gamma}$", "\\phi", "$\\tau_{\\gamma}$", "\\phi"),
                            summ_iid_hyperprec)
rownames(summ_iid_hyperprec) <- NULL
kable(summ_iid_hyperprec, caption = "Posterior distributions of BYM2 hyperparameters from Multinomial and Poisson models: fixed intercept, BYM2",
      booktabs = TRUE, digits = 5) %>% kable_styling(latex_options = "HOLD_position")

binrand.df <- as.data.frame(mod.bin$summary.random$state_code)
binrand.df$model <- "Multinomial"
ptrickrand.df <- as.data.frame(mod.ptrick$summary.random$state_code_cause)
ptrickrand.df$model <- "Poisson"
summ_re <- data.frame(multinomial = binrand.df$`0.5quant`,
                      Poisson = ptrickrand.df$`0.5quant`)
kable(summ_re, caption = "Posterior medians of BYM2 REs from Multinomial and Poisson models: fixed intercept, BYM2",
      booktabs = TRUE, digits = 5) %>% 
  kable_styling(latex_options = "HOLD_position")
summ_re <- data.frame(multinomial = binrand.df$sd,
                      Poisson = ptrickrand.df$sd)
kable(summ_re, caption = "Posterior SD of BYM2 REs from Multinomial and Poisson models: fixed intercept, BYM2",
      booktabs = TRUE, digits = 5) %>% 
  kable_styling(latex_options = "HOLD_position")

# ggplot CSMFs
ggplot(data = allres %>% filter(cause == causename), 
       aes(x = model, y = csmf_med, color = model)) +
  geom_errorbar(aes(ymin = csmf_lower, ymax = csmf_upper), width = .1) +
  geom_point() +
  facet_geo(~ state, grid = india_grid_MDS, scales = "free_y") +
  labs(x="model", 
       y="Cause fraction",
       subtitle = paste0("Estimated CSMF (",ci.level*100,"% posterior intervals) from regional national model with BYM2 REs"),
       title=unique(allChildCodesBroad$broadCodeName[allChildCodesBroad$cause == causename])) +
  theme_light()



