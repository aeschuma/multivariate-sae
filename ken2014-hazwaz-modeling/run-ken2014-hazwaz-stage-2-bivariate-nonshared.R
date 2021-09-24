## ----setup, include=FALSE--------------------------------------------------------------------
library(SUMMER)
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

## ----read-data-stage-1-----------------------------------------------------------------------
load("/Users/austin/Dropbox/dissertation_2/survey-csmf/data/ken_dhs2014/data/haz-waz-kenDHS2014.rda")
stage_1_list <- read_rds("/Users/austin/Dropbox/dissertation_2/survey-csmf/results/ken2014-hazwaz/ken2014-hazwaz-stage-1.rds")
results <- stage_1_list[["results"]]
V.array <- stage_1_list[["V.array"]]
n_regions <- nrow(poly.adm1)

## ----bi-nonshared-mod-b11------------------------------------------------------------------------

# hyperpars
rho_beta_1 <- 1
rho_beta_2 <- 1
sigma_normal_sd <- 1

datlist <- list(R = n_regions,
                regions = results$admin1,
                Sigma = V.array,
                y = results[, c("meanHAZ.bi", "meanWAZ.bi")],
                N_edges = node.info$n_edges,
                node1 = node.info$node1,
                node2 = node.info$node2,
                scaling_factor = scaling_factor,
                rho_beta_1 = rho_beta_1,
                rho_beta_2 = rho_beta_2,
                sigma_normal_sd = sigma_normal_sd)


stan_file <- "stan-models/nonshared-bym2.stan"

## stan options
niter <- 10000
nchains <- 3
prop_warmup <- 0.5
max_treedepth <- 250
adapt_delta <- 0.9
nthin <- 1
options(mc.cores = parallel::detectCores())

cmd_mod <- cmdstan_model(stan_file = stan_file)
fit <- cmd_mod$sample(data = datlist,
                      iter_warmup = niter*prop_warmup, iter_sampling = niter*(1-prop_warmup),
                      chains = nchains, thin = nthin,
                      adapt_delta = adapt_delta, max_treedepth = max_treedepth,
                      refresh =  0.1 * niter)
mod.stan.bi.nonshared <- rstan::read_stan_csv(fit$output_files())
params_to_extract <- c("beta", "sigma", "rho")
re_params <- c("v_1", "v_2", "u_1", "u_2")

# summaries
mod.stan.summary.bi.nonshared <- summary(mod.stan.bi.nonshared,
                                         pars = params_to_extract,
                                         probs = c(0.1, 0.5, 0.9))
kable(mod.stan.summary.bi.nonshared$summary, format = "markdown", digits = 3)

stan_diags <- data.frame(pct_divergent = get_num_divergent(mod.stan.bi.nonshared)/(niter * nchains * prop_warmup),
                         pct_max_tree_exceeded = get_num_max_treedepth(mod.stan.bi.nonshared)/(niter * prop_warmup * nchains),
                         pct_bmfi_low_chains = sum(is.finite(get_low_bfmi_chains(mod.stan.bi.nonshared)))/nchains)
print("Stan diagnostics")
stan_diags

# save results
res_list <- list(data = datlist,
                 priors = list(rho_beta_1 = rho_beta_1,
                               rho_beta_2 = rho_beta_2,
                               sigma_normal_sd = sigma_normal_sd),
                 stan_file = stan_file,
                 mod = list(mod.stan.bi.nonshared))

write_rds(res_list, 
          file = "../../../Dropbox/dissertation_2/survey-csmf/results/ken2014-hazwaz/ken2014-hazwaz-stage-2-bivariate-nonshared-b11.rds")

# sample_pars <- c("beta", "sigma", "rho")
# stan_trace(mod.stan.bi.nonshared, pars = sample_pars)
# 
# pairs_pars <- c("beta", "log_sigma", "rho")
# pairs(mod.stan.bi.nonshared, pars = pairs_pars)
# 
# mcmc_nuts_energy(nuts_params(mod.stan.bi.nonshared))

## ----bi-nonshared-mod-b33------------------------------------------------------------------------

# hyperpars
rho_beta_1 <- 3
rho_beta_2 <- 3
sigma_normal_sd <- 1

datlist <- list(R = n_regions,
                regions = results$admin1,
                Sigma = V.array,
                y = results[, c("meanHAZ.bi", "meanWAZ.bi")],
                N_edges = node.info$n_edges,
                node1 = node.info$node1,
                node2 = node.info$node2,
                scaling_factor = scaling_factor,
                rho_beta_1 = rho_beta_1,
                rho_beta_2 = rho_beta_2,
                sigma_normal_sd = sigma_normal_sd)


stan_file <- "stan-models/nonshared-bym2.stan"

## stan options
niter <- 10000
nchains <- 3
prop_warmup <- 0.5
max_treedepth <- 250
adapt_delta <- 0.9
nthin <- 1
options(mc.cores = parallel::detectCores())

cmd_mod <- cmdstan_model(stan_file = stan_file)
fit <- cmd_mod$sample(data = datlist,
                      iter_warmup = niter*prop_warmup, iter_sampling = niter*(1-prop_warmup),
                      chains = nchains, thin = nthin,
                      adapt_delta = adapt_delta, max_treedepth = max_treedepth,
                      refresh =  0.1 * niter)
mod.stan.bi.nonshared <- rstan::read_stan_csv(fit$output_files())
params_to_extract <- c("beta", "sigma", "rho")
re_params <- c("v_1", "v_2", "u_1", "u_2")

# summaries
mod.stan.summary.bi.nonshared <- summary(mod.stan.bi.nonshared,
                                         pars = params_to_extract,
                                         probs = c(0.1, 0.5, 0.9))
kable(mod.stan.summary.bi.nonshared$summary, format = "markdown", digits = 3)

stan_diags <- data.frame(pct_divergent = get_num_divergent(mod.stan.bi.nonshared)/(niter * nchains * prop_warmup),
                         pct_max_tree_exceeded = get_num_max_treedepth(mod.stan.bi.nonshared)/(niter * prop_warmup * nchains),
                         pct_bmfi_low_chains = sum(is.finite(get_low_bfmi_chains(mod.stan.bi.nonshared)))/nchains)
print("Stan diagnostics")
stan_diags

# save results
res_list <- list(data = datlist,
                 priors = list(rho_beta_1 = rho_beta_1,
                               rho_beta_2 = rho_beta_2,
                               sigma_normal_sd = sigma_normal_sd),
                 stan_file = stan_file,
                 mod = list(mod.stan.bi.nonshared))

write_rds(res_list, 
          file = "../../../Dropbox/dissertation_2/survey-csmf/results/ken2014-hazwaz/ken2014-hazwaz-stage-2-bivariate-nonshared-b33.rds")

