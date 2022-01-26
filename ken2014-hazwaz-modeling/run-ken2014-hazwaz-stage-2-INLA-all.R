# Austin Schumacher
# 1/21/2022
# Fit all INLA models and save results

# preamble ####
rm(list = ls())

library(tidyverse)
library(spdep)
library(geosphere)
library(haven)
library(knitr)
library(kableExtra)
library(magrittr)
library(rstan)
library(cmdstanr)
library(svyVGAM)
library(mvtnorm)
library(rgdal)
library(bayesplot)
library(INLA)
library(viridis)
library(classInt)
library(gridExtra)
library(ggpubr)
library(loo)


# functions ####
fitINLA <- function(formula, data, lincombs) {
    inla(formula, data = data, 
         family = "gaussian",
         lincomb = lincombs,
         control.family = list(hyper = list(prec = list(initial = 10, fixed=T))), 
         control.predictor = list(compute=T),
         control.compute = list(config=T),
         control.inla = list(lincomb.derived.correlation.matrix=T),
         control.fixed = list(prec = list(default = 0.001), correlation.matrix=T))
}


# read and format data ####
load("/Users/austin/Dropbox/dissertation_2/survey-csmf/data/ken_dhs2014/data/haz-waz-kenDHS2014.rda")
stage_1_list <- read_rds("/Users/austin/Dropbox/dissertation_2/survey-csmf/results/ken2014-hazwaz/ken2014-hazwaz-stage-1.rds")
results <- stage_1_list[["results"]]
V.array <- stage_1_list[["V.array"]]
n_regions <- nrow(poly.adm1)

## reformat data into long form
results.long <- results %>% select(admin1, admin1.name, admin1.char,
                                   meanHAZ.bi, meanWAZ.bi) %>%
    pivot_longer(cols = c(meanHAZ.bi, meanWAZ.bi),
                 names_to = "outcome",
                 names_prefix = "mean",
                 values_to = "value")
results.long$outcome <- ifelse(results.long$outcome == "HAZ.bi", "HAZ", "WAZ")
results.long$admin1.haz <- ifelse(results.long$outcome == "HAZ", results.long$admin1, NA)
results.long$admin1.waz <- ifelse(results.long$outcome == "WAZ", results.long$admin1, NA)
results.long$obs <- 1:nrow(results.long)

# create a list of the data
data <- as.list(results.long)
Nt <- length(unique(results.long$admin1))
N <- nrow(results.long)
data$weights <- rep(1, N)

## Define the index-vectors ii.1 ii.2 etc, which are the
## index's for the iid2d-model at timepoint 1, 2, ...
for(j in 1:Nt) {
    itmp <-  numeric(N)
    itmp[] <-  NA
    itmp[(j-1)*2 + 1] <-  1
    itmp[(j-1)*2 + 2] <-  2
    data <-  c(list(itmp), data)
    names(data)[1] <-  paste("ii.", j, sep="")
}

## add another index for the shared component
data$admin1.haz.2 <- data$admin1.haz

## add indeces for iid components in besag + iid formulation
data$admin1.haz.iid <- data$admin1.haz
data$admin1.waz.iid <- data$admin1.waz

# define model components ####

## priors ####
iid_prior <- list(prec = list(prior = "pc.prec",
                              param = c(1, 0.01)))
bym2_prior <- list(phi=list(prior="logitbeta", param=c(1, 1), initial=0.5), 
                   prec=list(prior="pc.prec", param=c(1, 0.01), initial=5))
lambda_prior <- list(beta = list(prior = 'logtnormal', param = c(0, 1)))
besag_prior <- list(prec=list(prior="pc.prec", param=c(1, 0.01), initial=5))

## Diagonal V formula additions ####

### we now have to add one 'iid2d' model for each observation pair,
### since their cov.matrix is different. we have to do this
### automatically... here I add numbers directly for simplicity

### univariate
add.univariate <- ""
for(j in 1:Nt) {
    init.prec.haz <-  log(1/results$seHAZ.bi[j]^2)
    init.prec.waz <-  log(1/results$seWAZ.bi[j]^2)
    
    add.univariate <-  paste(add.univariate, paste(" + 
                         f(", paste("ii.", j, sep=""), ", model=\"iid2d\", n=2,
                         hyper = list(
                         prec1 = list(
                         initial =", init.prec.haz,", 
                         fixed = TRUE),
                         prec2 = list(
                         initial =", init.prec.waz,", 
                         fixed = TRUE),
                         cor = list(
                         initial = log((1+", 0, ")/(1-", 0, ")), 
                         fixed = TRUE)))"))
    
}

### bivariate
add.bivariate <- ""
for(j in 1:Nt) {
    corr <-  results$corr.bi[j]
    init.prec.haz <-  log(1/results$seHAZ.bi[j]^2)
    init.prec.waz <-  log(1/results$seWAZ.bi[j]^2)
    
    add.bivariate <-  paste(add.bivariate, paste(" + 
                         f(", paste("ii.", j, sep=""), ", model=\"iid2d\", n=2,
                         hyper = list(
                         prec1 = list(
                         initial =", init.prec.haz,", 
                         fixed = TRUE),
                         prec2 = list(
                         initial =", init.prec.waz,", 
                         fixed = TRUE),
                         cor = list(
                         initial = log((1+", corr, ")/(1-", corr, ")), 
                         fixed = TRUE)))"))
    
}

## Model formulas ####
formula.univariate.iid <-  formula(paste("value ~ -1 + outcome + 
                                          f(admin1.haz, model = 'iid', hyper = iid_prior) + 
                                          f(admin1.waz, model = 'iid', hyper = iid_prior)",
                                         paste(add.univariate, collapse = " ")))
formula.univariate.bym <- formula(paste("value ~ -1 + outcome + 
                                         f(admin1.haz, model = 'bym2',
                                           graph = admin1.mat, 
                                           scale.model = T, 
                                           constr = T,
                                           hyper = bym2_prior) +
                                         f(admin1.waz, model = 'bym2',
                                           graph = admin1.mat, 
                                           scale.model = T, 
                                           constr = T,
                                           hyper = bym2_prior)",
                                        paste(add.univariate, collapse = " ")))
formula.bivariate.nonshared.iid <-  formula(paste("value ~ -1 + outcome + 
                                                   f(admin1.haz, model = 'iid', hyper = iid_prior) + 
                                                   f(admin1.waz, model = 'iid', hyper = iid_prior)",
                                                  paste(add.bivariate, collapse = " ")))
formula.bivariate.nonshared.bym <-  formula(paste("value ~ -1 + outcome + 
                                                   f(admin1.haz, model = 'bym2',
                                                     graph = admin1.mat, 
                                                     scale.model = T, 
                                                     constr = T,
                                                     hyper = bym2_prior) +
                                                   f(admin1.waz, model = 'bym2',
                                                     graph = admin1.mat, 
                                                     scale.model = T, 
                                                     constr = T,
                                                     hyper = bym2_prior)",
                                                  paste(add.bivariate, collapse = " ")))
formula.bivariate.shared.iid <-  formula(paste("value ~ -1 + outcome + 
                                                   f(admin1.haz, model = 'iid', hyper = iid_prior) + 
                                                   f(admin1.waz, model = 'iid', hyper = iid_prior) +
                                                   f(admin1.haz.2, copy = \"admin1.waz\", 
                                                     fixed = FALSE, hyper = lambda_prior)",
                                                  paste(add.bivariate, collapse = " ")))
formula.bivariate.shared.bym <-  formula(paste("value ~ -1 + outcome + 
                                                f(admin1.haz, model = 'bym2',
                                                  graph = admin1.mat, 
                                                  scale.model = T, 
                                                  constr = T,
                                                  hyper = bym2_prior) +
                                                f(admin1.waz, model = 'bym2',
                                                  graph = admin1.mat, 
                                                  scale.model = T, 
                                                  constr = T,
                                                  hyper = bym2_prior) +
                                                f(admin1.haz.2, copy = \"admin1.waz\", 
                                                  fixed = FALSE, hyper = lambda_prior)",
                                               paste(add.bivariate, collapse = " ")))
formula.bivariate.shared.bym.alt <-  formula(paste("value ~ -1 + outcome + 
                                                   f(admin1.haz, model = 'besag',
                                                     graph = admin1.mat, 
                                                     scale.model = T, 
                                                     constr = T,
                                                     hyper = besag_prior) +
                                                   f(admin1.haz.iid, model = 'iid', hyper = iid_prior) +
                                                   f(admin1.waz, model = 'besag',
                                                     graph = admin1.mat, 
                                                     scale.model = T, 
                                                     constr = T,
                                                     hyper = besag_prior) +
                                                   f(admin1.waz.iid, model = 'iid', hyper = iid_prior) +
                                                   f(admin1.haz.2, copy = \"admin1.waz\", 
                                                     fixed = FALSE, hyper = lambda_prior)",
                                                  paste(add.bivariate, collapse = " ")))

## linear combinations ####

# linear combination of preds without 2x2 REs
lc.vec.fe <- c(rep(1, n_regions))

diag.na.mat <- matrix(NA, nrow = n_regions, ncol = n_regions)
diag(diag.na.mat) <- 1
lc.vec.admin1.re <- diag.na.mat

## WAZ
lc.all.waz <- inla.make.lincombs(outcomeWAZ = lc.vec.fe, 
                                 admin1.waz = lc.vec.admin1.re)
names(lc.all.waz) <- gsub("lc", "waz.reg.", names(lc.all.waz))

## WAZ if besag + iid
lc.all.waz.alt <- inla.make.lincombs(outcomeWAZ = lc.vec.fe, 
                                     admin1.waz = lc.vec.admin1.re,
                                     admin1.waz.iid = lc.vec.admin1.re)
names(lc.all.waz.alt) <- gsub("lc", "waz.reg.", names(lc.all.waz.alt))

## HAZ nonshared
lc.all.haz.nonshared <- inla.make.lincombs(outcomeHAZ = lc.vec.fe, 
                                           admin1.haz = lc.vec.admin1.re)
names(lc.all.haz.nonshared) <- gsub("lc", "haz.reg.", names(lc.all.haz.nonshared))

## HAZ shared
lc.all.haz.shared <- inla.make.lincombs(outcomeHAZ = lc.vec.fe, 
                                        admin1.haz = lc.vec.admin1.re,
                                        admin1.haz.2 = lc.vec.admin1.re)
names(lc.all.haz.shared) <- gsub("lc", "haz.reg.", names(lc.all.haz.shared))

## HAZ shared if besag + iid
lc.all.haz.shared.alt <- inla.make.lincombs(outcomeHAZ = lc.vec.fe, 
                                            admin1.haz = lc.vec.admin1.re,
                                            admin1.haz.iid = lc.vec.admin1.re,
                                            admin1.haz.2 = lc.vec.admin1.re)
names(lc.all.haz.shared.alt) <- gsub("lc", "haz.reg.", names(lc.all.haz.shared.alt))

# Fit models ####

mod.univariate.iid <- fitINLA(formula = formula.univariate.iid, 
                              data = data, 
                              lincombs = c(lc.all.haz.nonshared, lc.all.waz))
mod.univariate.bym <- fitINLA(formula = formula.univariate.bym, 
                              data = data, 
                              lincombs = c(lc.all.haz.nonshared, lc.all.waz))
mod.bivariate.nonshared.iid <- fitINLA(formula = formula.bivariate.nonshared.iid, 
                                       data = data, 
                                       lincombs = c(lc.all.haz.nonshared, lc.all.waz))
mod.bivariate.nonshared.bym <- fitINLA(formula = formula.bivariate.nonshared.bym, 
                                       data = data, 
                                       lincombs = c(lc.all.haz.nonshared, lc.all.waz))
mod.bivariate.shared.iid <- fitINLA(formula = formula.bivariate.shared.iid, 
                                    data = data, 
                                    lincombs = c(lc.all.haz.shared, lc.all.waz))
mod.bivariate.shared.bym <- fitINLA(formula = formula.bivariate.shared.bym, 
                                    data = data, 
                                    lincombs = c(lc.all.haz.shared, lc.all.waz))
mod.bivariate.shared.bym.alt <- fitINLA(formula = formula.bivariate.shared.bym.alt, 
                                        data = data, 
                                        lincombs = c(lc.all.haz.shared.alt, lc.all.waz.alt))
