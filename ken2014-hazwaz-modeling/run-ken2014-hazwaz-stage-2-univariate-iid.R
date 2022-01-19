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
n_regions <- nrow(poly.adm1)

## ----uni-mod-b11---------------------------------------------------------------------------------
## HAZ

# hyperpars
sigma_normal_sd <- 1

stan_dat_HAZ <- list(R = n_regions,
                     Sigma = results$seHAZ.bi,
                     y = results$meanHAZ.bi,
                     sigma_normal_sd = sigma_normal_sd)


stan_file <- "stan-models/oneD-iid.stan"

## stan options
niter <- 10000
nchains <- 3
prop_warmup <- 0.5
max_treedepth <- 400
adapt_delta <- 0.95
nthin <- 1
options(mc.cores = parallel::detectCores())

cmd_mod <- cmdstan_model(stan_file = stan_file)
fit <- suppressMessages(cmd_mod$sample(data = stan_dat_HAZ,
                                       iter_warmup = niter*prop_warmup, 
                                       iter_sampling = niter*(1-prop_warmup),
                                       chains = nchains, thin = nthin,
                                       adapt_delta = adapt_delta, max_treedepth = max_treedepth,
                                       refresh = 0.1 * niter))
mod.stan.haz.uni <- rstan::read_stan_csv(fit$output_files())
params_to_extract <- c("beta", "sigma", "preds")
re_params <- c("v")

# summaries
mod.stan.summary.haz.uni <- summary(mod.stan.haz.uni,
                                    pars = params_to_extract,
                                    probs = c(0.1, 0.5, 0.9))
kable(mod.stan.summary.haz.uni$summary, format = "markdown", digits = 3)

stan_diags <- data.frame(pct_divergent = get_num_divergent(mod.stan.haz.uni)/(niter * nchains * prop_warmup),
                         pct_max_tree_exceeded = get_num_max_treedepth(mod.stan.haz.uni)/(niter * prop_warmup * nchains),
                         pct_bmfi_low_chains = sum(is.finite(get_low_bfmi_chains(mod.stan.haz.uni)))/nchains)
print("Stan diagnostics")
stan_diags

# sample_pars <- c("beta", "sigma", "rho")
# stan_trace(mod.stan.haz.uni, pars = sample_pars)
# pairs(mod.stan.haz.uni, pars = sample_pars)
# mcmc_nuts_energy(nuts_params(mod.stan.haz.uni))

## WAZ

## hyperpars
sigma_normal_sd <- 1

stan_dat_WAZ <- list(R = n_regions,
                     Sigma = results$seWAZ.bi,
                     y = results$meanWAZ.bi,
                     sigma_normal_sd = sigma_normal_sd)

## stan options
niter <- 10000
nchains <- 3
prop_warmup <- 0.5
max_treedepth <- 400
adapt_delta <- 0.95
nthin <- 1
options(mc.cores = parallel::detectCores())

cmd_mod <- cmdstan_model(stan_file = stan_file)
fit <- suppressMessages(cmd_mod$sample(data = stan_dat_WAZ,
                                       iter_warmup = niter*prop_warmup, 
                                       iter_sampling = niter*(1-prop_warmup),
                                       chains = nchains, thin = nthin,
                                       adapt_delta = adapt_delta, max_treedepth = max_treedepth,
                                       refresh =  0.1 * niter))
mod.stan.waz.uni <- rstan::read_stan_csv(fit$output_files())
params_to_extract <- c("beta", "sigma", "preds")
re_params <- c("v")

# summaries
mod.stan.summary.waz.uni <- summary(mod.stan.waz.uni,
                                    pars = params_to_extract,
                                    probs = c(0.1, 0.5, 0.9))
kable(mod.stan.summary.waz.uni$summary, format = "markdown", digits = 3)

stan_diags <- data.frame(pct_divergent = get_num_divergent(mod.stan.waz.uni)/(niter * nchains * prop_warmup),
                         pct_max_tree_exceeded = get_num_max_treedepth(mod.stan.waz.uni)/(niter * prop_warmup * nchains),
                         pct_bmfi_low_chains = sum(is.finite(get_low_bfmi_chains(mod.stan.waz.uni)))/nchains)
print("Stan diagnostics")
stan_diags

# sample_pars <- c("beta", "sigma", "rho")
# stan_trace(mod.stan.waz.uni, pars = sample_pars)

# pairs_pars <- c("beta", "log_sigma", "rho")
# pairs(mod.stan.waz.uni, pars = pairs_pars)

# mcmc_nuts_energy(nuts_params(mod.stan.waz.uni))

# save results
res_list <- list(data = list(stan_dat_HAZ, 
                             stan_dat_WAZ),
                 priors = list(sigma_normal_sd = sigma_normal_sd),
                 stan_file = stan_file,
                 mod = list(mod.stan.haz.uni,
                            mod.stan.waz.uni))

write_rds(res_list, 
          file = "../../../Dropbox/dissertation_2/survey-csmf/results/ken2014-hazwaz/ken2014-hazwaz-stage-2-univariate-iid.rds")

