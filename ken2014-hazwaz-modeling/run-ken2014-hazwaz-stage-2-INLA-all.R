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

# define model components ####

## priors
iid_prior <- list(prec = list(prior = "pc.prec",
                              param = c(1, 0.01)))
bym2_prior <- list(phi=list(prior="logitbeta", param=c(1, 1), initial=0.5), 
                   prec=list(prior="pc.prec", param=c(1, 0.01), initial=5))
lambda_prior <- list(beta = list(prior = 'logtnormal', param = c(0, 1)))


# "Univariate" diagonal V models

## we now have to add one 'iid2d' model for each observation pair,
## since their cov.matrix is different. we have to do this
## automatically... here I add numbers directly for simplicity
add <- ""
for(j in 1:Nt) {
    corr <-  results$corr.bi[j]
    init.prec.haz <-  log(1/results$seHAZ.bi[j]^2)
    init.prec.waz <-  log(1/results$seWAZ.bi[j]^2)
    
    add <-  paste(add, paste(" + 
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

## Model formula
formula <-  formula(paste("value ~ -1 + outcome + f(admin1.haz, model = 'iid', hyper = iid_prior) + f(admin1.waz, model = 'iid', hyper = iid_prior) + ",
                    paste(add, collapse = " ")))

# linear combination of preds without 2x2 REs
lc.vec.haz.fe <- c(rep(1, n_regions))
lc.vec.waz.fe <- c(rep(1, n_regions))
diag.na.mat <- matrix(NA, nrow = n_regions, ncol = n_regions)
diag(diag.na.mat) <- 1
lc.vec.admin1.haz.re <- diag.na.mat
lc.vec.admin1.waz.re <- diag.na.mat

lc.all.haz <- inla.make.lincombs(outcomeHAZ = lc.vec.haz.fe, 
                                 admin1.haz = lc.vec.admin1.haz.re)
lc.all.waz <- inla.make.lincombs(outcomeWAZ = lc.vec.waz.fe, 
                                 admin1.waz = lc.vec.admin1.waz.re)
names(lc.all.haz) <- gsub("lc", "haz.reg.", names(lc.all.haz))
names(lc.all.waz) <- gsub("lc", "waz.reg.", names(lc.all.waz))

## Fit model
mod.leigh <- inla(formula, data = data, 
                  family = "gaussian",
                  lincomb = c(lc.all.haz, lc.all.waz),
                  control.family = list(hyper = list(prec = list(initial = 10, fixed=T))), 
                  control.predictor = list(compute=T),
                  control.compute = list(config=T),
                  control.inla = list(lincomb.derived.correlation.matrix=T),
                  control.fixed = list(prec = list(default = 0.001), correlation.matrix=T) )

