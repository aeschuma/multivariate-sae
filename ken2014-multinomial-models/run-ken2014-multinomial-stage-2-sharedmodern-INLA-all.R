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
library(svyVGAM)
library(mvtnorm)
library(rgdal)
library(INLA)
library(viridis)
library(classInt)
library(gridExtra)
library(ggpubr)


# functions ####
fitINLA <- function(formula, data, lincombs) {
    inla(formula, data = data, 
         family = "gaussian",
         quantiles = c(0.025, 0.1, 0.5, 0.9, 0.975),
         lincomb = lincombs,
         control.family = list(hyper = list(prec = list(initial = 10, fixed=T))), 
         control.predictor = list(compute=T),
         control.compute = list(config=T, waic = TRUE, dic = TRUE, cpo = TRUE),
         control.inla = list(lincomb.derived.correlation.matrix=T),
         control.fixed = list(prec = list(default = 0.001), correlation.matrix=T))
}

# read and format data ####
load("/Users/austin/Dropbox/dissertation_2/survey-csmf/data/ken_dhs2014/data/haz-waz-kenDHS2014.rda")
stage_1_list <- read_rds("/Users/austin/Dropbox/dissertation_2/survey-csmf/results/ken2014-multinomial/ken2014-multinomial-stage-1.rds")
results <- stage_1_list[["results"]]
V.array <- stage_1_list[["V.array"]]
n_regions <- nrow(poly.adm1)

## reformat data into long form
results.long <- results %>% select(admin1, admin1.name, admin1.char,
                                   mean_beta1, mean_beta2) %>%
    pivot_longer(cols = c(mean_beta1, mean_beta2),
                 names_to = "outcome",
                 names_prefix = "mean_",
                 values_to = "value")
results.long$admin1.beta1 <- ifelse(results.long$outcome == "beta1", results.long$admin1, NA)
results.long$admin1.beta2 <- ifelse(results.long$outcome == "beta2", results.long$admin1, NA)
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
data$admin1.beta2.2 <- data$admin1.beta2

## add indeces for iid components in besag + iid formulation
data$admin1.beta1.iid <- data$admin1.beta1
data$admin1.beta2.iid <- data$admin1.beta2

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
    init.prec.beta1 <-  log(1/results$se_beta1[j]^2)
    init.prec.beta2 <-  log(1/results$se_beta2[j]^2)
    
    add.univariate <-  paste(add.univariate, paste(" + 
                         f(", paste("ii.", j, sep=""), ", model=\"iid2d\", n=2,
                         hyper = list(
                         prec1 = list(
                         initial =", init.prec.beta1,", 
                         fixed = TRUE),
                         prec2 = list(
                         initial =", init.prec.beta2,", 
                         fixed = TRUE),
                         cor = list(
                         initial = log((1+", 0, ")/(1-", 0, ")), 
                         fixed = TRUE)))"))
    
}

### bivariate
add.bivariate <- ""
for(j in 1:Nt) {
    corr <-  results$corr[j]
    init.prec.beta1 <-  log(1/results$se_beta1[j]^2)
    init.prec.beta2 <-  log(1/results$se_beta2[j]^2)
    
    add.bivariate <-  paste(add.bivariate, paste(" + 
                         f(", paste("ii.", j, sep=""), ", model=\"iid2d\", n=2,
                         hyper = list(
                         prec1 = list(
                         initial =", init.prec.beta1,", 
                         fixed = TRUE),
                         prec2 = list(
                         initial =", init.prec.beta2,", 
                         fixed = TRUE),
                         cor = list(
                         initial = log((1+", corr, ")/(1-", corr, ")), 
                         fixed = TRUE)))"))
    
}

## Model formulas ####
formula.univariate.iid <-  formula(paste("value ~ -1 + outcome + 
                                          f(admin1.beta1, model = 'iid', hyper = iid_prior) + 
                                          f(admin1.beta2, model = 'iid', hyper = iid_prior)",
                                         paste(add.univariate, collapse = " ")))
formula.univariate.bym <- formula(paste("value ~ -1 + outcome + 
                                         f(admin1.beta1, model = 'bym2',
                                           graph = admin1.mat, 
                                           scale.model = T, 
                                           constr = T,
                                           hyper = bym2_prior) +
                                         f(admin1.beta2, model = 'bym2',
                                           graph = admin1.mat, 
                                           scale.model = T, 
                                           constr = T,
                                           hyper = bym2_prior)",
                                        paste(add.univariate, collapse = " ")))
formula.bivariate.nonshared.iid <-  formula(paste("value ~ -1 + outcome + 
                                                   f(admin1.beta1, model = 'iid', hyper = iid_prior) + 
                                                   f(admin1.beta2, model = 'iid', hyper = iid_prior)",
                                                  paste(add.bivariate, collapse = " ")))
formula.bivariate.nonshared.bym <-  formula(paste("value ~ -1 + outcome + 
                                                   f(admin1.beta1, model = 'bym2',
                                                     graph = admin1.mat, 
                                                     scale.model = T, 
                                                     constr = T,
                                                     hyper = bym2_prior) +
                                                   f(admin1.beta2, model = 'bym2',
                                                     graph = admin1.mat, 
                                                     scale.model = T, 
                                                     constr = T,
                                                     hyper = bym2_prior)",
                                                  paste(add.bivariate, collapse = " ")))
formula.bivariate.shared.iid <-  formula(paste("value ~ -1 + outcome + 
                                                   f(admin1.beta1, model = 'iid', hyper = iid_prior) + 
                                                   f(admin1.beta2, model = 'iid', hyper = iid_prior) +
                                                   f(admin1.beta2.2, copy = \"admin1.beta1\", 
                                                     fixed = FALSE, hyper = lambda_prior)",
                                                  paste(add.bivariate, collapse = " ")))
formula.bivariate.shared.bym <-  formula(paste("value ~ -1 + outcome + 
                                                f(admin1.beta1, model = 'bym2',
                                                  graph = admin1.mat, 
                                                  scale.model = T, 
                                                  constr = T,
                                                  hyper = bym2_prior) +
                                                f(admin1.beta2, model = 'bym2',
                                                  graph = admin1.mat, 
                                                  scale.model = T, 
                                                  constr = T,
                                                  hyper = bym2_prior) +
                                                f(admin1.beta2.2, copy = \"admin1.beta1\", 
                                                  fixed = FALSE, hyper = lambda_prior)",
                                               paste(add.bivariate, collapse = " ")))
formula.bivariate.shared.bym.alt <-  formula(paste("value ~ -1 + outcome + 
                                                   f(admin1.beta1, model = 'besag',
                                                     graph = admin1.mat, 
                                                     scale.model = T, 
                                                     constr = T,
                                                     hyper = besag_prior) +
                                                   f(admin1.beta1.iid, model = 'iid', hyper = iid_prior) +
                                                   f(admin1.beta2, model = 'besag',
                                                     graph = admin1.mat, 
                                                     scale.model = T, 
                                                     constr = T,
                                                     hyper = besag_prior) +
                                                   f(admin1.beta2.iid, model = 'iid', hyper = iid_prior) +
                                                   f(admin1.beta2.2, copy = \"admin1.beta1\", 
                                                     fixed = FALSE, hyper = lambda_prior)",
                                                  paste(add.bivariate, collapse = " ")))

## linear combinations ####

# linear combination of preds without 2x2 REs
lc.vec.fe <- c(rep(1, n_regions))

diag.na.mat <- matrix(NA, nrow = n_regions, ncol = n_regions)
diag(diag.na.mat) <- 1
lc.vec.admin1.re <- diag.na.mat

## WAZ
lc.all.2 <- inla.make.lincombs(outcomebeta2 = lc.vec.fe, 
                               admin1.beta2 = lc.vec.admin1.re)
names(lc.all.2) <- gsub("lc", "oucome2.reg.", names(lc.all.2))

## WAZ if besag + iid
lc.all.2.alt <- inla.make.lincombs(outcomebeta2 = lc.vec.fe, 
                                   admin1.beta2 = lc.vec.admin1.re,
                                   admin1.beta2.iid = lc.vec.admin1.re)
names(lc.all.2.alt) <- gsub("lc", "outcome2.reg.", names(lc.all.2.alt))

## HAZ nonshared
lc.all.1.nonshared <- inla.make.lincombs(outcomebeta1 = lc.vec.fe, 
                                           admin1.beta1 = lc.vec.admin1.re)
names(lc.all.1.nonshared) <- gsub("lc", "outcome1.reg.", names(lc.all.1.nonshared))

## HAZ shared
lc.all.1.shared <- inla.make.lincombs(outcomebeta1 = lc.vec.fe, 
                                        admin1.beta1 = lc.vec.admin1.re,
                                        admin1.beta2.2 = lc.vec.admin1.re)
names(lc.all.1.shared) <- gsub("lc", "outcome1.reg.", names(lc.all.1.shared))

## HAZ shared if besag + iid
lc.all.1.shared.alt <- inla.make.lincombs(outcomebeta1 = lc.vec.fe, 
                                            admin1.beta1 = lc.vec.admin1.re,
                                            admin1.beta1.iid = lc.vec.admin1.re,
                                            admin1.beta2.2 = lc.vec.admin1.re)
names(lc.all.1.shared.alt) <- gsub("lc", "outcome1.reg.", names(lc.all.1.shared.alt))

# Fit models ####

mod.univariate.iid <- fitINLA(formula = formula.univariate.iid, 
                              data = data, 
                              lincombs = c(lc.all.1.nonshared, lc.all.2))
mod.univariate.bym <- fitINLA(formula = formula.univariate.bym, 
                              data = data, 
                              lincombs = c(lc.all.1.nonshared, lc.all.2))
mod.bivariate.nonshared.iid <- fitINLA(formula = formula.bivariate.nonshared.iid, 
                                       data = data, 
                                       lincombs = c(lc.all.1.nonshared, lc.all.2))
mod.bivariate.nonshared.bym <- fitINLA(formula = formula.bivariate.nonshared.bym, 
                                       data = data, 
                                       lincombs = c(lc.all.1.nonshared, lc.all.2))
mod.bivariate.shared.iid <- fitINLA(formula = formula.bivariate.shared.iid, 
                                    data = data, 
                                    lincombs = c(lc.all.1.shared, lc.all.2))
mod.bivariate.shared.bym <- fitINLA(formula = formula.bivariate.shared.bym, 
                                    data = data, 
                                    lincombs = c(lc.all.1.shared, lc.all.2))
mod.bivariate.shared.bym.alt <- fitINLA(formula = formula.bivariate.shared.bym.alt, 
                                        data = data, 
                                        lincombs = c(lc.all.1.shared.alt, lc.all.2.alt))

# test new lambda prior
# lambda_prior <- list(beta = list(prior = 'normal', param = c(1, 2)))
# mod.bivariate.shared.bym.2 <- fitINLA(formula = formula.bivariate.shared.bym, 
#                                       data = data, 
#                                       lincombs = c(lc.all.haz.shared, lc.all.waz))
# mod.bivariate.shared.bym$summary.hyperpar[5,]
# mod.bivariate.shared.bym.2$summary.hyperpar[5,]

inla.results <- list(mod.univariate.iid,
                     mod.univariate.bym,
                     mod.bivariate.nonshared.iid,
                     mod.bivariate.nonshared.bym,
                     mod.bivariate.shared.iid,
                     mod.bivariate.shared.bym,
                     mod.bivariate.shared.bym.alt)
names(inla.results) <- c("Univariate IID",
                         "Univariate BYM",
                         "Bivariate nonshared IID", 
                         "Bivariate nonshared BYM",
                         "Bivariate shared IID", 
                         "Bivariate shared BYM",
                         "Bivariate shared Besag + IID")

# Save Results ####
write_rds(inla.results, file = "../../../Dropbox/dissertation_2/survey-csmf/results/ken2014-multinomial/ken2014-multinomial-stage-2-sharedmodern-inla-all.rds")
