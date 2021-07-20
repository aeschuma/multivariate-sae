# testing univariate normal and poisson models with shared IID components

rm(list=ls())

## set the root depending on operating system
root <- ifelse(Sys.info()[1]=="Darwin","~/",
               ifelse(Sys.info()[1]=="Windows","P:/",
                      ifelse(Sys.info()[1]=="Linux","/home/students/aeschuma/",
                             stop("Unknown operating system"))))

## the following code makes rstan work on the Box server and the cluster
if (root == "P:/") {
    Sys.setenv(HOME="C:/Users/aeschuma",
               R_USER="C:/Users/aeschuma",
               R_LIBS_USER="C:/Users/aeschuma/R_libraries")
    .libPaths("C:/Users/aeschuma/R_libraries")
} else if (root == "/home/students/aeschuma/") {
    Sys.setenv(HOME=root,
               R_USER=root,
               R_LIBS_USER=paste0(root,"R/x86_64-pc-linux-gnu-library/3.6"))
    .libPaths(paste0(root,"R/x86_64-pc-linux-gnu-library/3.6"))
}

## load libraries
if (root == "/home/students/aeschuma/") {
    library(Rcpp,lib.loc=paste0(root,"R/x86_64-pc-linux-gnu-library/3.6"));
    library(StanHeaders,lib.loc=paste0(root,"R/x86_64-pc-linux-gnu-library/3.6")); 
    library(BH,lib.loc = paste0(root,"R/x86_64-pc-linux-gnu-library/3.6"))
    library(rstan,lib.loc = paste0(root,"R/x86_64-pc-linux-gnu-library/3.6"));
} else {
    library(Rcpp); library(StanHeaders); library(BH); library(rstan);library(bayesplot);
}
library(mvtnorm); library(MASS);library(gtools); library(parallel);
library(scales); library(RColorBrewer);library(data.table);
library(ggplot2); library(tidyverse); library(magrittr);

if (root == "~/") library(cmdstanr);

# working directory for code
wd <- paste0(root,"Desktop/survey-csmf/stage-2-simulation")

# directory to save results
savedir <- paste0(root,"Dropbox/dissertation_2/survey-csmf/results/stage-2-simulation")

## set directory
setwd(wd)

# set number of cores to the max possible
options(mc.cores = detectCores())

## data generation options
number_of_causes <- 2
number_of_regions <- 20
number_of_replications <- 10

## parameters
beta1 <- 1
beta2 <- 2
beta <- c(beta1, beta2)

sigma_gamma1 <- 1
sigma_gamma2 <- 1

lambda <- 1.5
sigma_delta <- 1

sigmasq_lower <- 0.01
sigmasq_upper <- 0.1

## stan options
niter <- 1000
nchains <- 1
prop_warmup <- 0.5
max_treedepth <- 15
adapt_delta <- 0.8
    
# simulate data
N <- number_of_regions * number_of_replications * number_of_causes
dat_frame <- expand_grid(reg = 1:number_of_regions,
                         rep = 1:number_of_replications,
                         cause = 1:number_of_causes)

sigmasqs <- runif(N, sqrt(sigmasq_lower), sqrt(sigmasq_upper))
epsilons <- rnorm(N, 0, sigmasqs)

dat_frame %<>% mutate(beta = beta[cause],
                      epsilon = epsilons)

gamma1 <- rnorm(number_of_regions, 0, sigma_gamma1)
gamma1_dat <- tibble(reg = 1:number_of_regions,
                     cause = 1,
                     gamma1 = gamma1)
gamma2 <- rnorm(number_of_regions, 0, sigma_gamma2)
gamma2_dat <- tibble(reg = 1:number_of_regions,
                     cause = 2,
                     gamma2 = gamma2)
delta <- rnorm(number_of_regions, 0, sigma_delta)
delta_dat <- tibble(reg = rep(1:number_of_regions, number_of_causes),
                    cause = rep(1:number_of_causes, each = number_of_regions),
                    delta_mod = c(delta * lambda, delta / lambda))

dat_frame %<>% left_join(gamma1_dat) %>% left_join(gamma2_dat) %>% left_join(delta_dat) %>%
    mutate(gamma1 = replace(gamma1, is.na(gamma1), 0),
           gamma2 = replace(gamma2, is.na(gamma2), 0)) %>%
    mutate(eta = beta + gamma1 + gamma2 + delta_mod,
           y_norm = beta + gamma1 + gamma2 + delta_mod + epsilon,
           y_pois = rpois(N, exp(beta + gamma1 + gamma2 + delta_mod)))



# fit STAN model
if (root == "~/") {
    stan_list <- fitSTAN(model_to_run, simulated_data$datlist,
                         niter = niter, nchains = nchains, nthin = 1, prop_warmup = prop_warmup,
                         adapt_delta = adapt_delta, max_treedepth = max_treedepth, 
                         cmd = TRUE)
} else {
    stan_list <- fitSTAN(model_to_run, simulated_data$datlist,
                         niter = niter, nchains = nchains, nthin = 1, prop_warmup = prop_warmup,
                         adapt_delta = adapt_delta, max_treedepth = max_treedepth, 
                         cmd = FALSE)
}

# save parameter names in order to extract and save results
param_names <- names(simulated_data$params)
params_to_extract <- param_names[!(param_names %in% c(grep("gamma_", param_names, value = TRUE),
                                                      grep("delta_", param_names, value = TRUE)))]

# summaries
mod_summary <- summary(stan_list$mod_stan,
                       pars = params_to_extract,
                       probs = c(0.1, 0.5, 0.9))
posterior_qs <- mod_summary$summary[, c("10%", "50%", "90%")]
extracted_params <- rownames(posterior_qs)

# compile results
results <- data.frame(param = rownames(posterior_qs),
                      absolute_bias = NA,
                      relative_bias = NA,
                      coverage = NA,
                      width = NA)
for (i in 1:nrow(results)) {
    tmp_param_name <- rownames(posterior_qs)[i]
    if (tmp_param_name == "beta[1]") {
        tmp_param <- beta1
    } else if (tmp_param_name == "beta[2]") {
        tmp_param <- beta2
    } else if (tmp_param_name == "beta[3]") {
        tmp_param <- beta3
    } else if (tmp_param_name == "sigma_gamma[1]") {
        tmp_param <- sigma_gamma1
    } else if (tmp_param_name == "sigma_gamma[2]") {
        tmp_param <- sigma_gamma2
    } else {
        tmp_param <- get(tmp_param_name)
    }
    results$absolute_bias[i] <- posterior_qs[tmp_param_name, "50%"] - tmp_param
    results$relative_bias[i] <- (posterior_qs[tmp_param_name, "50%"] - tmp_param)/tmp_param
    results$coverage[i] <- (tmp_param > posterior_qs[tmp_param_name, "10%"]) & (tmp_param < posterior_qs[tmp_param_name, "90%"])
    results$width[i] <- posterior_qs[tmp_param_name, "90%"] - posterior_qs[tmp_param_name, "10%"]
}

stan_diags <- data.frame(pct_divergent = get_num_divergent(stan_list$mod_stan)/(niter * prop_warmup),
                         pct_max_tree_exceeded = get_num_max_treedepth(stan_list$mod_stan)/(niter * prop_warmup),
                         pct_bmfi_low_chains = sum(is.finite(get_low_bfmi_chains(stan_list$mod_stan)))/nchains)

if (!testing) {
    # save results
    setwd(savedir)
    saveRDS(results, paste0("tmp/results_run-", run_number, "_sim-", sim, ".rds"))
    saveRDS(stan_diags, paste0("tmp/standiags_run-", run_number, "_sim-", sim, ".rds"))
}

if (testing) {
    
    # make diagnostic plots
    stan_trace(stan_list$mod_stan, pars = c("beta", "sigma_gamma", "lambda", "sigma_delta"))
    pairs(stan_list$mod_stan, pars = c("beta", "sigma_gamma", "lambda", "sigma_delta"))
    # mcmc_areas(as.matrix(stan_list$mod_stan),
    #            pars = c("beta[1]", "beta[2]", 
    #                     "sigma_gamma[1]", "sigma_gamma[2]",
    #                     "lambda", "sigma_delta"),
    #            prob = 0.8)
    mcmc_nuts_energy(nuts_params(stan_list$mod_stan))
    
    # plot gammas vs posterior median gamma-hat
    alpha_hat1 <- summary(stan_list$mod_stan, 
                          pars = c("alpha1"), probs = 0.5)$summary[, "50%"] - 
        summary(stan_list$mod_stan, 
                pars = c("beta[1]"), probs = 0.5)$summary[, "50%"]
    alpha_hat2 <- summary(stan_list$mod_stan, 
                          pars = c("alpha2"), probs = 0.5)$summary[, "50%"] - 
        summary(stan_list$mod_stan, 
                pars = c("beta[2]"), probs = 0.5)$summary[, "50%"]
    delta_hat <- summary(stan_list$mod_stan, 
                         pars = c("delta"), probs = 0.5)$summary[, "50%"]
    
    mean(alpha_hat1)
    mean(alpha_hat2)
    mean(delta_hat)
    par(mfrow = c(1, 3), mar = c(5, 4, 4, 2), oma = c(1, 1, .25, .25))
    # plot((simulated_data$params$gamma_rc[,1] + simulated_data$params$beta[1]) ~ alpha_hat1,
    #      xlab = "Cause 1",
    #      ylab = "True alphas", 
    #      main = "Posterior median of the REs",
    #      pch = 19, col = alpha("dodgerblue", 0.5))
    # abline(0, 1, col = "darkgreen")
    # plot((simulated_data$params$gamma_rc[,2] + simulated_data$params$beta[2]) ~ alpha_hat2,
    #      xlab = "Cause 2",
    #      ylab = "True alphas", 
    #      pch = 19, col = alpha("indianred", 0.5))
    # abline(0, 1, col = "darkgreen")
    # plot(simulated_data$params$delta_rc ~ delta_hat,
    #      xlab = "Shared RE",
    #      ylab = "True deltas", 
    #      pch = 19, col = alpha("orchid3", 0.5))
    # abline(0, 1, col = "darkgreen")
    plot((simulated_data$params$gamma_rc[,1]) ~ alpha_hat1,
         xlab = "Cause 1",
         ylab = "True alphas", 
         main = "Posterior median of the REs",
         pch = 19, col = alpha("dodgerblue", 0.5))
    abline(0, 1, col = "darkgreen")
    plot((simulated_data$params$gamma_rc[,2]) ~ alpha_hat2,
         xlab = "Cause 2",
         ylab = "True alphas", 
         pch = 19, col = alpha("indianred", 0.5))
    abline(0, 1, col = "darkgreen")
    plot(simulated_data$params$delta_rc ~ delta_hat,
         xlab = "Shared RE",
         ylab = "True deltas", 
         pch = 19, col = alpha("orchid3", 0.5))
    abline(0, 1, col = "darkgreen")
}
