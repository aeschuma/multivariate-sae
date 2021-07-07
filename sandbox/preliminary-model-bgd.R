rm(list = ls())
setwd('~/Desktop/survey-csmf/sandbox')

## Libraries ####
library(SUMMER)
# help(package = "SUMMER", help_type = "html")
library(classInt)
library(RColorBrewer)
library(dplyr)
library(tidyr)
library(rgdal)
library(scales)
library(survey)
library(tidyverse)
library(magrittr)
library(svyVGAM)
library(fields)
library(INLA)
library(mvtnorm)

#### Functions ####

# function for calculating 5q0 from the beta coefficients output from the multnomial regression
get_5q0 <- function(beta, n) {
    ## For testing
    # beta <- coef(mult.mod)
    # n <- c(1, 11, 12, 12, 12, 12)
    
    betas_of_interest <- beta[seq(1,length(beta), by = 2)]
    betas_other <- beta[seq(2,length(beta), by = 2)]
    
    gamma_j <- (1 + exp(betas_of_interest) + exp(betas_other))^-1
    sum_nj_of_gammajs <- rep(0, length(n))
    for (j in 1:length(n)) {
        for (s in 1:n[j]) {
            sum_nj_of_gammajs[j] <- sum_nj_of_gammajs[j] + gamma_j[j]^s
        }
    }
    
    phi_of_interest <- sum(exp(betas_of_interest) * sum_nj_of_gammajs)
    phi_other <- sum(exp(betas_other) * sum_nj_of_gammajs)
    
    return(c(phi_of_interest, phi_other))
}

#### Parameters ####

coi <- "random causes of interest"

iso <- "bgd"
start_year <- 1990
end_year <- 2020
svy_year <- 2017

# xs <- c(0, 1, 12, 24, 36, 48)
# ns <- c(1, 11, 12, 12, 12, 12)
xs <- c(0, 12)
ns <- c(12, 48)
n_age <- length(ns)
n_cause <- 2

## Load data ####

load(paste0('../../../Dropbox/dissertation_2/survey-csmf/data/bgd_dhs2017/data/', 'births-file_', iso, '.rda'))

mod.dat.all$years <- as.numeric(as.character(mod.dat.all$time))
dat.years <- sort(unique(mod.dat.all$years))
beg.years <- seq(start_year,end_year,5)
end.years <- beg.years + 4
periods <- paste(beg.years, end.years, sep = "-")
mod.dat.all$period <- as.character(cut(mod.dat.all$years, breaks = c(beg.years, beg.years[length(beg.years)]+5),
                                   include.lowest = T, right = F, labels = periods))
mod.dat.all$v005 <- mod.dat.all$v005/1e6
mod.dat.all$ones <- 1

# reshape wide by cause of death
mod.dat.wide <- mod.dat.all %>% mutate(alive = ifelse(cause.gp == "alive", 1, 0),
                                       cause_of_interest = ifelse(cause.gp == coi, 1, 0),
                                       cause_other = ifelse(cause.gp == "other", 1, 0)) %>%
    mutate(dead = ifelse(alive == 1, 0, 1))

# compactify
mod.dat.comp <- mod.dat.wide[, c("v001", "v024", "v025", "v005", "age", "strata", "time", 
                                 "dead", "alive", "cause_of_interest", "cause_other")]
mod.dat.comp$total <- 1
formula <- as.formula(".~age + time + strata + v001 + v024 + v025 + v005")
mod.dat.comp <- aggregate(formula, data = mod.dat.comp, FUN = sum, drop = TRUE)
table(mod.dat.comp$age)

# survey design
my.svydesign <- survey::svydesign(ids = ~ v001, 
                                  strata = ~ strata, nest = T, weights = ~v005, 
                                  data = mod.dat.comp)

## national model
mult.mod <- svy_vglm(cbind(cause_of_interest, cause_other, alive) ~ -1 + age, 
                     family = multinomial, 
                     design = my.svydesign)

summary(mult.mod)

# extract betas and sigma-hat to use in calculating 5q0 and its variance
betas <- coef(mult.mod)
V <- stats::vcov(mult.mod)

#  simulate betas
betasim <- rmvnorm(10000, mean = betas, sigma = V)
logit_phisim <- t(apply(betasim, 1, function(x) logitlink(get_5q0(x, ns))))
    
## partial derivatives 
expbetas <- exp(betas)
betamat <- matrix(betas, nrow = n_age, ncol = n_cause, byrow = TRUE)
expbetamat <- exp(betamat)
gammaj <- (1 + apply(expbetamat, 1, sum))^-1
xi <- matrix(0, nrow = n_age, ncol = n_cause)
for (c in 1:n_cause) {
    for (j in 1:n_age) {
        for (s in 1:ns[j]) {
            xi[j, c] <- xi[j, c] + gammaj[j]^s
        }
        xi[j, c] <- xi[j, c] * expbetamat[j, c]
    }
}
phi <- apply(xi, 2, sum)
d_phi_beta_concordant <- matrix(0, nrow = n_age, ncol = n_cause)
for (c in 1:n_cause) {
    for (j in 1:n_age) {
        for (s in 1:ns[j]) {
            d_phi_beta_concordant[j, c] <- d_phi_beta_concordant[j, c] + (gammaj[j]^s * (1 - s * expbetamat[j, c] * gammaj[j]))
        }
        d_phi_beta_concordant[j, c] <- d_phi_beta_concordant[j, c] * expbetamat[j, c]
    }
}
d_phi_beta_concordant_v <- as.vector(d_phi_beta_concordant)

for (c in 1:n_cause) {
    tmpsum <- rep(0, n_age)
    for (j in 1:n_age) {
        for (s in 1:ns[j]) {
            tmpsum[j] <- tmpsum[j] + (-1 * s * gammaj[j]^s)
        }
    }
}
d_phi_beta_discordant <- matrix(rep(apply(expbetamat, 1, prod) * gammaj[j] * tmpsum[j], 2), 
                                nrow = n_age, ncol = n_cause)
d_phi_beta_discordant_v <- as.vector(d_phi_beta_discordant)

d_phi_d_beta_v <- c()
counter <- 0
for (j in 1:n_age) {
    for (c in 1:n_cause) {
        counter <- counter + 1
        d_phi_d_beta_v <- c(d_phi_d_beta_v, d_phi_beta_concordant_v[counter], d_phi_beta_discordant_v[counter])
    }
}
d_phi_d_beta <- matrix(d_phi_d_beta_v, nrow = n_age*n_cause, ncol = n_cause)

d_eta_d_phi <- diag((1/phi) - (1/(1 - phi)), nrow = length(phi), ncol = length(phi))

## Items to return ##
logit.mean.est <- logitlink(get_5q0(betas, ns))
var.est <- var(logit_phisim)
var.est.delta <- t(d_eta_d_phi) %*% t(d_phi_d_beta) %*% V %*% d_phi_d_beta %*% d_eta_d_phi

results_national <- list(logit.mean.est = logit.mean.est,
                         var.est = var.est)

## regional models
regions <- unique(mod.dat.comp$v024)

results <- vector(mode = "list", length = length(regions))

for (r in 1:length(regions)) {
    
    # testing
    # r <- 1
    
    tmp <- subset(my.svydesign, v024 == regions[r])
    
    # stage 1 for each region
    mult.mod <- svy_vglm(cbind(cause_of_interest, cause_other, alive) ~ -1 + age, 
                         family = multinomial, 
                         design = tmp)
    
    summary(mult.mod)
    
    bin.mod <- svyglm(cbind(cause_other, alive+cause_of_interest) ~ -1 + age, 
                      family = binomial, 
                      design = tmp)
    summary(bin.mod)
    betas.bin <- coef(bin.mod)
    V.bin <- vcov(bin.mod)
    
    # extract sigma-hat
    betas <- coef(mult.mod)
    V <- stats::vcov(mult.mod)
    
    #  simulate betas
    betasim <- rmvnorm(10000, mean = betas, sigma = V)
    logit_phisim <- t(apply(betasim, 1, function(x) logitlink(get_5q0(x, ns))))
    
    ## partial derivatives ##TODO
    # var.est <- t(d_eta_d_phi) %*% t(d_phi_d_beta) %*% V %*% d_phi_d_beta %*% d_eta_d_phi 
    
    ## Items to return ##
    logit.mean.est <- logitlink(get_5q0(betas, ns))
    var.est <- var(logit_phisim)
    
    results[[r]] <- list(logit.mean.est = logit.mean.est, 
                         var.est = var.est)
}

# create block diagonal precision matrix and vector of means for use in INLA function to model
# also make a matrix with correlations
diag.list <- vector(mode = "list", length = length(results))
corr.list <- vector(mode = "list", length = length(results))
means <- c()
for (i in 1:length(results))  {
    diag.list[[i]] <- solve(results[[i]]$var.est)
    D <- diag(sqrt(diag(results[[i]]$var.est)))
    DInv <- solve(D)
    corr.list[[i]] <- DInv %*% results[[i]]$var.est %*% DInv
    means <- c(means, results[[i]]$logit.mean.est)
}
Vdes_prec <- bdiag(diag.list)

## second stage model

## load polygon shape files
poly.adm0 <- readOGR(dsn = "/Users/austin/Dropbox/dissertation_2/survey-csmf/data/shapefiles/bgd_adm_bbs_20201113_SHP",
                     layer = "bgd_admbnda_adm0_bbs_20201113")
poly.adm0$region <- "All"
poly.adm1 <- readOGR(dsn = "/Users/austin/Dropbox/dissertation_2/survey-csmf/data/shapefiles/bgd_adm_bbs_20201113_SHP",
                     layer = "bgd_admbnda_adm1_bbs_20201113")
poly.adm1$region <- tolower(poly.adm1$ADM1_EN)

# adjacency matrix
admin1.mat.raw <- poly2nb(SpatialPolygons(poly.adm1@polygons))
admin1.mat <- nb2mat(admin1.mat.raw, style = "B")
nb2INLA("bgd.graph", admin1.mat.raw)

# data
dat <- data.frame(y = means, 
                  reg = rep(as.numeric(regions), each = 2),
                  cause = rep(c(1, 2), length(results)))
dat <- dat[order(dat$reg, dat$cause),]
dat$obs <- 1:nrow(dat)

# priors
fe.prec <- list(prec.intercept = 0,
                prec = 0)
cov_prior <- list(prec = list(prior = "loggamma", param = c(100000, 100000)))
bym2_priors <- list(phi=list(prior="pc", param=c(0.5, 0.5), initial=1), 
                    prec=list(prior="pc.prec", param=c(1, 0.01), initial=5))

# FE on each ccause only model
# model formula: we have a generic RE with a specified C matrix and a tight prior on the precision centered on 1
# m.form.feonly <- y ~ -1 + factor(cause) + f(obs,  model='generic0', Cmatrix = Vdes_prec,
#                                             hyper = cov_prior)
# 
# mod.feonly <- inla(m.form.feonly, 
#                    data = dat,
#                    family = "gaussian",
#                    control.fixed = fe.prec,
#                    control.predictor=list(compute=TRUE),
#                    control.compute=list(config = TRUE),
#                    control.family = list(hyper = list(prec = list(initial = log(1), fixed=TRUE))),
#                    scale = 1000000)
# 
# summary(mod.feonly)

# RE on each cause only model, NO INTERCEPT
# model formula: we have a generic RE with a specified C matrix and a tight prior on the precision centered on 1
# m.form.recause <- y ~ -1 + f(cause, model = "iid") + f(obs,  model='generic0', Cmatrix = Vdes_prec,
#                                             hyper = cov_prior)
# 
# mod.recause <- inla(m.form.recause, 
#                     data = dat,
#                     family = "gaussian",
#                     control.fixed = fe.prec,
#                     control.predictor=list(compute=TRUE),
#                     control.compute=list(config = TRUE),
#                     control.family = list(hyper = list(prec = list(initial = log(1), fixed=TRUE))),
#                     scale = 1000000)
# 
# summary(mod.recause)

# RE on each cause only model, with intercept
# model formula: we have a generic RE with a specified C matrix and a tight prior on the precision centered on 1
# m.form.recause.int <- y ~ 1 + f(cause, model = "iid") + f(obs,  model='generic0', Cmatrix = Vdes_prec,
#                                                        hyper = cov_prior)
# 
# mod.recause.int <- inla(m.form.recause.int, 
#                     data = dat,
#                     family = "gaussian",
#                     control.fixed = fe.prec,
#                     control.predictor=list(compute=TRUE),
#                     control.compute=list(config = TRUE),
#                     control.family = list(hyper = list(prec = list(initial = log(1), fixed=TRUE))),
#                     scale = 1000000)
# 
# summary(mod.recause.int)

# FE on each cause, IID RE on region,  no overall intercept
m.form.iidre <- y ~ -1 + factor(cause) + 
    f(obs,  model='generic0', Cmatrix = Vdes_prec, hyper = cov_prior) +
    f(reg, model = 'iid', constr = TRUE)

mod.iidre <- inla(m.form.iidre, 
                  data = dat,
                  family = "gaussian",
                  control.fixed = fe.prec,
                  control.predictor=list(compute=TRUE),
                  control.compute=list(config = TRUE),
                  control.family = list(hyper = list(prec = list(initial = log(1), fixed=TRUE))),
                  scale = 1000000)

summary(mod.iidre)

# RE on each cause, IID RE on region, overall intercept
# m.form.recause.iidre <- y ~ 1 + f(cause, model = "iid") + f(obs,  model='generic0', Cmatrix = Vdes_prec,
#                                            hyper = cov_prior) + f(reg, model = 'iid')
# 
# mod.recause.iidre <- inla(m.form.recause.iidre, 
#                           data = dat,
#                           family = "gaussian",
#                           control.fixed = fe.prec,
#                           control.predictor=list(compute=TRUE),
#                           control.compute=list(config = TRUE),
#                           control.family = list(hyper = list(prec = list(initial = log(1), fixed=TRUE))),
#                           scale = 1000000)
# 
# summary(mod.recause.iidre)

# FE on each cause, BYM2 on region,  no overall intercept
# m.form.bym2 <- y ~ -1 + factor(cause) + 
#     f(obs,  model='generic0', Cmatrix = Vdes_prec, hyper = cov_prior) + 
#     f(reg, model = 'bym2', 
#       graph = "bgd.graph",
#       scale.model = TRUE,
#       constr = TRUE,
#       hyper = bym2_priors)
# 
# mod.bym2 <- inla(m.form.bym2, 
#                   data = dat,
#                   family = "gaussian",
#                   control.fixed = fe.prec,
#                   control.predictor=list(compute=TRUE),
#                   control.compute=list(config = TRUE),
#                   control.family = list(hyper = list(prec = list(initial = log(1), fixed=TRUE))),
#                   scale = 1000000)
# 
# summary(mod.bym2)

## look at region-specific data amount
# for (i in 1:length(regions)) {
#     cat(paste0("\n REGION: ", regions[i], "\n"))
#     r1 <- mod.dat.comp %>% filter(v024 == regions[i])
#     cat(paste0("cause of interest \n"))
#     print(table(r1$age, r1$cause_of_interest))
#     cat(paste0("\n other cause \n"))
#     print(table(r1$age, r1$cause_other))
#     cat(paste0("\n"))
# }


## testing out direct and smoothed direct

