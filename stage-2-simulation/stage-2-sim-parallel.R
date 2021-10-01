rm(list=ls())

# Preamble ####

## set the root depending on operating system
root <- ifelse(Sys.info()[1]=="Darwin","~/",
               ifelse(Sys.info()[1]=="Windows","P:/",
                      ifelse(Sys.info()[1]=="Linux","/home/users/aeschuma/",
                             stop("Unknown operating system"))))

## the following code makes rstan work on the Box server and the cluster
if (root == "P:/") {
    Sys.setenv(HOME="C:/Users/aeschuma",
               R_USER="C:/Users/aeschuma",
               R_LIBS_USER="C:/Users/aeschuma/R_libraries")
    .libPaths("C:/Users/aeschuma/R_libraries")
} else if (root == "/home/users/aeschuma/") {
    Sys.setenv(HOME=root,
               R_USER=root,
               R_LIBS_USER=paste0(root,"R/x86_64-pc-linux-gnu-library/3.6"))
    .libPaths(paste0(root,"R/x86_64-pc-linux-gnu-library/3.6"))
}

## load libraries
if (root == "/home/users/aeschuma/") {
    library(Rcpp,lib.loc=paste0(root,"R/x86_64-pc-linux-gnu-library/3.6"));
    library(StanHeaders,lib.loc=paste0(root,"R/x86_64-pc-linux-gnu-library/3.6")); 
    library(BH,lib.loc = paste0(root,"R/x86_64-pc-linux-gnu-library/3.6"))
    library(rstan,lib.loc = paste0(root,"R/x86_64-pc-linux-gnu-library/3.6"));
} else {
    library(Rcpp); library(StanHeaders); library(BH); library(rstan); library(bayesplot); 
}
library(mvtnorm); library(MASS);library(gtools); library(parallel);
library(scales); library(RColorBrewer); library(ggplot2); library(tidyr); library(dplyr);
library(haven); library(knitr); library(INLA); library(readr);

if (root == "~/") library(cmdstanr);

## TESTING THE CODE?
testing <- FALSE

## define directories

# working directory for code
wd <- paste0(root,"Desktop/survey-csmf/stage-2-simulation")

# directory to save results
savedir <- paste0(root,"Dropbox/dissertation_2/survey-csmf/results/stage-2-simulation")

## set directory
setwd(wd)

# set number of cores to the max possible
options(mc.cores = detectCores())

# Load data and functions ####

# load simulation functions:
#   simulateData()
#   fitSTAN()
cat(paste("load sim functions \n"))
source("simulation-functions.R")

# load model csv to load which model we're running
cat(paste("load model info \n"))
models_dat <- read_csv("model-info.csv")
# only if need to rewrite csv to get rid of warning message for incomplete final line
# write.csv(models_dat, file = "model-info.csv", row.names = FALSE)
cat(paste("set parameters from command args \n"))
# Set parameters! ####
if (testing) {
    ## which model to run
    m_number <- 2
    
    ## data generation options
    dgm <- 3
    
    ## stan options
    niter <- 5000
    nchains <- 2
    prop_warmup <- 0.5
    max_treedepth <- 25
    adapt_delta <- 0.8
    
    ## which run
    run_number <- 999
    
    ## which simulation
    sim <- 499
} else {
    ## which model to run
    m_number <- as.numeric(commandArgs(trailingOnly=TRUE)[1])

    ## data generation options
    dgm <- as.numeric(commandArgs(trailingOnly=TRUE)[2])
    
    ## stan options
    niter <- as.numeric(commandArgs(trailingOnly=TRUE)[3])
    nchains <- as.numeric(commandArgs(trailingOnly=TRUE)[4])
    prop_warmup <- as.numeric(commandArgs(trailingOnly=TRUE)[5])
    max_treedepth <- as.numeric(commandArgs(trailingOnly=TRUE)[6])
    adapt_delta <- as.numeric(commandArgs(trailingOnly=TRUE)[7])
    
    ## which run
    run_number <- as.numeric(commandArgs(trailingOnly=TRUE)[8])
    
    ## which simulation
    sim <- as.numeric(commandArgs(trailingOnly=TRUE)[9])
}

cat(paste("Set modeling info \n"))
# set the model to run
my_model <- models_dat %>% dplyr::filter(model_number == m_number) %>% dplyr::select(-1) %>% dplyr::slice(1) %>% unlist()

# load data
dgm_file <- read_csv("dgm-info.csv")
my_dgm <- dgm_file %>% dplyr::filter(dgm_number == dgm) %>% dplyr::select(-dgm_number) %>% dplyr::slice(1)
load(my_dgm %>% pull(geo_data))
random_re <- my_dgm %>% pull(random_re)
my_dgm <- my_dgm %>% dplyr::select(-geo_data, -random_re)

cat(paste("Simulate data \n"))
# Simulate data ####
simulated_data <- simulateData(dgm_specs = my_dgm, 
                               Amat = admin1.mat, 
                               scaling_factor = scaling_factor, 
                               seed_re = ifelse(random_re, sim + 500, 98125), 
                               seed_lik = sim,
                               testing = FALSE)

# update with priors
simulated_data$datlist$rho_beta_1 <- as.numeric(my_model["rho_beta_1"])
simulated_data$datlist$rho_beta_2 <- as.numeric(my_model["rho_beta_2"])
simulated_data$datlist$sigma_normal_sd <- as.numeric(my_model["sigma_normal_sd"])

# Fit STAN model ####
cat(paste("Fit stan model \n"))
if (root == "~/") {
    stan_list <- fitSTAN(stan_file = my_model["stan_model_file"], 
                         data = simulated_data$datlist,
                         niter = niter, nchains = nchains, nthin = 1, prop_warmup = prop_warmup,
                         adapt_delta = adapt_delta, max_treedepth = max_treedepth, 
                         cmd = TRUE)
} else {
    stan_list <- fitSTAN(stan_file = my_model["stan_model_file"], 
                         data = simulated_data$datlist,
                         niter = niter, nchains = nchains, nthin = 1, prop_warmup = prop_warmup,
                         adapt_delta = adapt_delta, max_treedepth = max_treedepth, 
                         cmd = FALSE)
}

# Calculate summary measures ####

# save parameter names in order to extract and save results
params_to_extract <- names(simulated_data$params[["mean_pars"]])
if (!(my_model["stan_model_file"] %in% c("stan-models/shared-bym2-coregionalization.stan"))) {
    params_to_extract <- params_to_extract[which(params_to_extract != "lambda")]
}

cat(paste("Extract summaries for parameters \n"))

# summaries
mod_summary <- summary(stan_list$mod_stan,
                       pars = params_to_extract,
                       probs = c(0.025, 0.1, 0.5, 0.9, 0.975))
posterior_qs <- mod_summary$summary[, c("2.5%", "10%", "50%", "90%", "97.5%")]

cat(paste("Calculate results for parameters \n"))
# compile results
sim_res <- data.frame(param = rownames(posterior_qs),
                      absolute_bias = NA,
                      relative_bias = NA,
                      coverage.80 = NA,
                      width.80 = NA,
                      coverage.95 = NA,
                      width.95 = NA)
for (i in 1:nrow(sim_res)) {
    tmp_param_name <- rownames(posterior_qs)[i]
    tmp_param <- simulated_data$params[["mean_pars"]][tmp_param_name]
    
    sim_res$absolute_bias[i] <- posterior_qs[tmp_param_name, "50%"] - tmp_param
    sim_res$relative_bias[i] <- (posterior_qs[tmp_param_name, "50%"] - tmp_param)/tmp_param
    sim_res$coverage.80[i] <- (tmp_param > posterior_qs[tmp_param_name, "10%"]) & (tmp_param < posterior_qs[tmp_param_name, "90%"])
    sim_res$width.80[i] <- posterior_qs[tmp_param_name, "90%"] - posterior_qs[tmp_param_name, "10%"]
    sim_res$coverage.95[i] <- (tmp_param > posterior_qs[tmp_param_name, "2.5%"]) & (tmp_param < posterior_qs[tmp_param_name, "97.5%"])
    sim_res$width.95[i] <- posterior_qs[tmp_param_name, "97.5%"] - posterior_qs[tmp_param_name, "2.5%"]
}

cat(paste("Extract mean preds \n"))
# posterior observation means
mod_pred_summary <- summary(stan_list$mod_stan,
                            pars = "preds",
                            probs = c(0.025, 0.1, 0.5, 0.9, 0.975))$summary[, c("2.5%", "10%", "50%", "90%", "97.5%")]

cat(paste("Calculate results for mean preds \n"))

mod_pred_50 <- matrix(mod_pred_summary[, "50%"], ncol = 2, byrow = TRUE)
mod_pred_10 <- matrix(mod_pred_summary[, "10%"], ncol = 2, byrow = TRUE)
mod_pred_90 <- matrix(mod_pred_summary[, "90%"], ncol = 2, byrow = TRUE)
mod_pred_025 <- matrix(mod_pred_summary[, "2.5%"], ncol = 2, byrow = TRUE)
mod_pred_975 <- matrix(mod_pred_summary[, "97.5%"], ncol = 2, byrow = TRUE)

real_means <- simulated_data$params$bivariate_means

absolute_bias_pred <- mean(mod_pred_50 - real_means)
rel_bias_pred <- mean((mod_pred_50 - real_means) / real_means)

cov80_pred <- mean( (real_means >= mod_pred_10) & (real_means <= mod_pred_90) )
cov95_pred <- mean( (real_means >= mod_pred_025) & (real_means <= mod_pred_975) )

width80_pred <- mean(mod_pred_90 - mod_pred_10)
width95_pred <- mean(mod_pred_975 - mod_pred_025)

mean_pred_res <- data.frame(param = "mean_preds",
                            absolute_bias = absolute_bias_pred,
                            relative_bias = rel_bias_pred,
                            coverage.80 = cov80_pred,
                            width.80 = width80_pred,
                            coverage.95 = cov95_pred,
                            width.95 = width95_pred)

cat(paste("Create results output \n"))

# add onto results
sim_res <- rbind(sim_res, mean_pred_res)

cat(paste("Calculate stan diagnostics \n"))

# stan diagnostics
stan_diags <- data.frame(pct_divergent = get_num_divergent(stan_list$mod_stan)/(niter * nchains * prop_warmup),
                         pct_max_tree_exceeded = get_num_max_treedepth(stan_list$mod_stan)/(niter * nchains * prop_warmup),
                         pct_bmfi_low_chains = sum(is.finite(get_low_bfmi_chains(stan_list$mod_stan)))/nchains)

# Save results ####
cat(paste("Save results \n"))
if (!testing) {
    # save results
    setwd(savedir)
    saveRDS(sim_res, paste0("tmp/results_run-", run_number, "_sim-", sim, ".rds"))
    saveRDS(stan_diags, paste0("tmp/standiags_run-", run_number, "_sim-", sim, ".rds"))
}

# Make plots if testing ####
if (testing) {
    
    # make diagnostic plots
    # stan_trace(stan_list$mod_stan, pars = params_to_extract)
    # pairs(stan_list$mod_stan, pars = params_to_extract)
    # mcmc_areas(as.matrix(stan_list$mod_stan),
    #            pars = c("beta[1]", "beta[2]", 
    #                     "sigma_gamma[1]", "sigma_gamma[2]",
    #                     "lambda", "sigma_delta"),
    #            prob = 0.8)
    # mcmc_nuts_energy(nuts_params(stan_list$mod_stan))
    
    # plot REs (estimated vs. truth)
    
}
