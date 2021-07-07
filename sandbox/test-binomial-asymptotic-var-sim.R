# Testing asymptotic variance simulation in the binomial case

library(SUMMER)
library(mvtnorm)
library(survey)
library(tidyverse)
library(magrittr)
library(haven)

rm(list = ls())

# parameters
iso <- "bgd"
start_year <- 1990
end_year <- 2020
svy_year <- 2017

# load and format data
load(paste0('../../../Dropbox/dissertation_2/survey-csmf/data/bgd_dhs2017/data/', 'births-file_', iso, '.rda'))

mod.dat.all$years <- as.numeric(as.character(mod.dat.all$time))
dat.years <- sort(unique(mod.dat.all$years))
beg.years <- seq(start_year,end_year,5)
end.years <- beg.years + 4
periods <- paste(beg.years, end.years, sep = "-")
mod.dat.all$period <- as.character(cut(mod.dat.all$years, breaks = c(beg.years, beg.years[length(beg.years)]+5),
                                       include.lowest = T, right = F, labels = periods))
mod.dat.all$v005 <- mod.dat.all$v005/1e6

mod.dat.comp <- mod.dat.all[, c("v001", "v024", "v025", "v005", "age", "strata", "time", 
                                "died")]
mod.dat.comp$total <- 1
formula <- as.formula(".~age + time + strata + v001 + v024 + v025 + v005")
mod.dat.comp <- aggregate(formula, data = mod.dat.comp, FUN = sum, drop = TRUE)

# fit binomial direct 
my.svydesign <- survey::svydesign(ids = ~ v001, 
                                  strata = ~ strata, nest = T, weights = ~v005, 
                                  data = mod.dat.comp)

bin.mod <- svyglm(cbind(died, total) ~ -1 + age,
                  family = quasibinomial, 
                  design = my.svydesign)



betas <- coef(bin.mod)
V <- vcov(bin.mod)

# estimates with delta method var
ns <- c(1, 11, 12, 12, 12, 12)
probs <- expit(betas)

mean.est <- (1 - prod((1 - probs)^ns, na.rm = TRUE))
logit.mean.est <- logitlink((1 - prod((1 - probs)^ns, na.rm = TRUE)))  #*1000  

## partial derivatives ##
gamma <- prod((1 + exp(betas))^ns, na.rm = TRUE)
derivatives <- (gamma)/(gamma - 1) * ns * expit(betas)

## Items to return ##
var.est <- t(derivatives) %*% V %*% derivatives

# simulate asymptotic var
#  simulate betas
expitbetasim <- expit(rmvnorm(100000, mean = betas, sigma = V))
mean.est.sim <- apply(expitbetasim, 1, function(x) {1 - prod((1 - x)^(c(1, 11, 12, 12, 12, 12)), na.rm = TRUE)})
logit.mean.est.sim <- logitlink(mean.est.sim)
var.est.sim <- var(logit.mean.est.sim)

# compare
cat(paste0("delta method var: ", var.est, "\nsimul method var: ", var.est.sim))
