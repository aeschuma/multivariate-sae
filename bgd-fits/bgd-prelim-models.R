# testing out bgd fits, empirical parameters, and estimated parameters from simpler models

rm(list = ls())

library(rstan); library(cmdstanr); library(tidyverse); library(parallel); library(bayesplot);
library(mvtnorm); library(SUMMER); library(INLA); library(survey); library(svyVGAM); 
library(rgdal); library(spdep);

options(mc.cores = parallel::detectCores())

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

load(paste0('../../data/bgd_dhs2017/data/', 'births-file_', iso, '.rda'))

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

## regional models
regions <- unique(mod.dat.comp$v024)
region_dat <- tibble(reg_num = 1:length(regions),
                     reg = regions)
results <- vector(mode = "list", length = length(regions))
n_regions <- length(regions)

for (r in 1:n_regions) {
    
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

# create array of fixed covs and vector of means
# also make a matrix with correlations
V.array <- array(NA, dim = c(n_regions, n_cause, n_cause))
corr.array <- array(NA, dim = c(n_regions, n_cause, n_cause))
means <- matrix(NA, nrow = n_regions, ncol = n_cause)
for (i in 1:length(results))  {
    Vtmp <- results[[i]]$var.est
    V.array[i, , ] <- Vtmp
    D <- diag(sqrt(diag(Vtmp)))
    DInv <- solve(D)
    corr.array[i, , ] <- DInv %*% Vtmp %*% DInv
    means[i, ] <- results[[i]]$logit.mean.est
}

# get fixed cov values to inform simulations
apply(V.array, c(2, 3), mean)
apply(V.array, c(2, 3), min)
apply(V.array, c(2, 3), max)
apply(corr.array, c(2, 3), mean)
apply(corr.array, c(2, 3), min)
apply(corr.array, c(2, 3), max)
apply(V.array, c(2, 3), function(x) mean(sqrt(x)))
apply(V.array, c(2, 3), function(x) min(sqrt(x)))
apply(V.array, c(2, 3), function(x) max(sqrt(x)))

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

# create node lists
node1 <- c()
node2 <- c()
n_edges <- 0
for (i in 1:n_regions) {
    for (j in 1:i) {
        if(admin1.mat[i,j] == 1) {
            node1 <- c(node1, i)
            node2 <- c(node2, j)
            n_edges <- n_edges + 1
        }
    }
}

# MODEL 1

# hyperparams
tau_gamma_hyper_alpha <- 3.2761
tau_gamma_hyper_beta <- 1.81
tau_phi_hyper_alpha <- 1
tau_phi_hyper_beta <- 1

# create data
datlist <- list(N = nrow(means),
                R = n_regions,
                regions = region_dat$reg_num,
                Sigma = V.array,
                y = means,
                N_edges = n_edges,
                node1 = node1,
                node2 = node2,
                tau_gamma_hyper_alpha = tau_gamma_hyper_alpha,
                tau_gamma_hyper_beta = tau_gamma_hyper_beta,
                tau_phi_hyper_alpha = tau_phi_hyper_alpha,
                tau_phi_hyper_beta = tau_phi_hyper_beta)

# inits_chain1 <- list(beta = c(-3, -3.5))

## second stage models
stan_file <- "../shared-component-models/stan-models/nonshared-icar.stan"

## stan options
niter <- 6000
nchains <- 2
prop_warmup <- 0.5
max_treedepth <- 200
adapt_delta <- 0.75
nthin <- 1

cmd_mod <- cmdstan_model(stan_file = stan_file)
fit <- cmd_mod$sample(data = datlist,
                      iter_warmup = niter*prop_warmup, iter_sampling = niter*(1-prop_warmup),
                      chains = nchains, thin = nthin,
                      adapt_delta = adapt_delta, max_treedepth = max_treedepth
                      # , init = list(inits_chain1)
                      # , refresh = 0
                      )
mod_stan <- rstan::read_stan_csv(fit$output_files())

params_to_extract <- c("beta", "sigma_gamma", "sigma_phi")

# summaries
mod_summary <- summary(mod_stan,
                       pars = params_to_extract,
                       probs = c(0.1, 0.5, 0.9))
mod_summary$summary
stan_diags <- data.frame(pct_divergent = get_num_divergent(mod_stan)/(niter * prop_warmup),
                         pct_max_tree_exceeded = get_num_max_treedepth(mod_stan)/(niter * prop_warmup * nchains),
                         pct_bmfi_low_chains = sum(is.finite(get_low_bfmi_chains(mod_stan)))/nchains)
stan_diags

sample_pars <- c("beta", "tau_gamma", "tau_phi")
stan_trace(mod_stan, pars = sample_pars)
pairs(mod_stan, pars = sample_pars)
mcmc_nuts_energy(nuts_params(mod_stan))

gamma_hat1 <- summary(mod_stan, pars = c("gamma1"), probs = c(0.1, 0.5, 0.9))$summary[, c("10%", "50%" , "90%")]
gamma_hat2 <- summary(mod_stan, pars = c("gamma2"), c(0.1, 0.5, 0.9))$summary[, c("10%", "50%" , "90%")]

phi_hat1 <- summary(mod_stan, pars = c("phi_1"), c(0.1, 0.5, 0.9))$summary[, c("10%", "50%" , "90%")]
phi_hat2 <- summary(mod_stan, pars = c("phi_2"), c(0.1, 0.5, 0.9))$summary[, c("10%", "50%" , "90%")]

summary(gamma_hat1)
summary(gamma_hat2)
summary(phi_hat1)
summary(phi_hat2)


# MODEL 2

# hyperparams
tau_gamma_hyper_alpha <- 20
tau_gamma_hyper_beta <- 1
tau_phi_hyper_alpha <- 0.25
tau_phi_hyper_beta <- 0.25

# create data
datlist <- list(N = nrow(means),
                R = n_regions,
                regions = region_dat$reg_num,
                Sigma = V.array,
                y = means,
                N_edges = n_edges,
                node1 = node1,
                node2 = node2,
                tau_gamma_hyper_alpha = tau_gamma_hyper_alpha,
                tau_gamma_hyper_beta = tau_gamma_hyper_beta,
                tau_phi_hyper_alpha = tau_phi_hyper_alpha,
                tau_phi_hyper_beta = tau_phi_hyper_beta)

# inits_chain1 <- list(beta = c(-3, -3.5))

## second stage models
stan_file <- "../shared-component-models/stan-models/nonshared-icar.stan"

## stan options
niter <- 6000
nchains <- 2
prop_warmup <- 0.5
max_treedepth <- 200
adapt_delta <- 0.75
nthin <- 1

cmd_mod <- cmdstan_model(stan_file = stan_file)
fit <- cmd_mod$sample(data = datlist,
                      iter_warmup = niter*prop_warmup, iter_sampling = niter*(1-prop_warmup),
                      chains = nchains, thin = nthin,
                      adapt_delta = adapt_delta, max_treedepth = max_treedepth
                      # , init = list(inits_chain1)
                      # , refresh = 0
)
mod_stan2 <- rstan::read_stan_csv(fit$output_files())

params_to_extract <- c("beta", "sigma_gamma", "sigma_phi")

# summaries
mod_summary2 <- summary(mod_stan2,
                        pars = params_to_extract,
                        probs = c(0.1, 0.5, 0.9))
mod_summary2$summary
stan_diags2 <- data.frame(pct_divergent = get_num_divergent(mod_stan2)/(niter * prop_warmup),
                          pct_max_tree_exceeded = get_num_max_treedepth(mod_stan2)/(niter * prop_warmup * nchains),
                          pct_bmfi_low_chains = sum(is.finite(get_low_bfmi_chains(mod_stan2)))/nchains)
stan_diags2

sample_pars <- c("beta", "tau_gamma", "tau_phi")
stan_trace(mod_stan2, pars = sample_pars)
pairs(mod_stan2, pars = sample_pars)
mcmc_nuts_energy(nuts_params(mod_stan2))

gamma_hat1a <- summary(mod_stan2, pars = c("gamma1"), probs = c(0.1, 0.5, 0.9))$summary[, c("10%", "50%" , "90%")]
gamma_hat2a <- summary(mod_stan2, pars = c("gamma2"), c(0.1, 0.5, 0.9))$summary[, c("10%", "50%" , "90%")]

phi_hat1a <- summary(mod_stan2, pars = c("phi_1"), c(0.1, 0.5, 0.9))$summary[, c("10%", "50%" , "90%")]
phi_hat2a <- summary(mod_stan2, pars = c("phi_2"), c(0.1, 0.5, 0.9))$summary[, c("10%", "50%" , "90%")]

summary(gamma_hat1a)
summary(gamma_hat2a)
summary(phi_hat1a)
summary(phi_hat2a)

# compare models 1 and 2

mod_summary$summary[, c("10%", "50%", "90%")]
mod_summary2$summary[, c("10%", "50%", "90%")]

par(mfrow = c(2,2))
plot(gamma_hat1[, "50%"] ~ gamma_hat1a[, "50%"], 
     pch = as.character(1:n_regions))
abline(0, 1, col = "darkgreen")
plot(gamma_hat2[, "50%"] ~ gamma_hat2a[, "50%"], 
     pch = as.character(1:n_regions))
abline(0, 1, col = "darkgreen")
plot(phi_hat1[, "50%"] ~ phi_hat1a[, "50%"], 
     pch = as.character(1:n_regions))
abline(0, 1, col = "darkgreen")
plot(phi_hat2[, "50%"] ~ phi_hat2a[, "50%"], 
     pch = as.character(1:n_regions))
abline(0, 1, col = "darkgreen")

plot(gamma_hat1[, "10%"] ~ gamma_hat1a[, "10%"], 
     pch = as.character(1:n_regions))
abline(0, 1, col = "darkgreen")
plot(gamma_hat2[, "10%"] ~ gamma_hat2a[, "10%"], 
     pch = as.character(1:n_regions))
abline(0, 1, col = "darkgreen")
plot(phi_hat1[, "10%"] ~ phi_hat1a[, "10%"], 
     pch = as.character(1:n_regions))
abline(0, 1, col = "darkgreen")
plot(phi_hat2[, "10%"] ~ phi_hat2a[, "10%"], 
     pch = as.character(1:n_regions))
abline(0, 1, col = "darkgreen")

plot(gamma_hat1[, "90%"] ~ gamma_hat1a[, "90%"], 
     pch = as.character(1:n_regions))
abline(0, 1, col = "darkgreen")
plot(gamma_hat2[, "90%"] ~ gamma_hat2a[, "90%"], 
     pch = as.character(1:n_regions))
abline(0, 1, col = "darkgreen")
plot(phi_hat1[, "90%"] ~ phi_hat1a[, "90%"], 
     pch = as.character(1:n_regions))
abline(0, 1, col = "darkgreen")
plot(phi_hat2[, "90%"] ~ phi_hat2a[, "90%"], 
     pch = as.character(1:n_regions))
abline(0, 1, col = "darkgreen")
