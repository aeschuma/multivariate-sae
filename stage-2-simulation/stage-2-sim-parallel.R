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
library(ggplot2); 

########
## TESTING THE CODE?
########

testing <- FALSE

## define directories

# working directory for code
wd <- paste0(root,"Desktop/survey-csmf/stage-2-simulation")

# directory to save results
savedir <- paste0(root,"Dropbox/dissertation_2/survey-csmf/results/stage-2-simulation")

## set directory
setwd(wd)

# set number of cores to the max possible
options(mc.cores = parallel::detectCores())

# load simulation functions:
#   simulateData()
#   fitSTAN()
source("simulation-functions.R")

# load model csv to load which model we're running
models_dat <- read.csv("model-info.csv")
# only if need to rewrite csv to get rid of warning message for incomplete final line
# write.csv(models_dat, file = "model-info.csv", row.names = FALSE)

## set parameters!
if (testing) {
    ## which model to run
    model_number <- 4
    model_to_run <- models_dat$model_name[models_dat$model_number == model_number]
    
    ## data generation options
    number_of_causes <- 2
    number_of_regions <- 15
    number_of_replications <- 5
    
    ## parameters
    beta1 <- 1
    beta2 <- 2
    rho_lower <- 0
    rho_upper <- 0.3
    sigmasq_lower <- 0.025
    sigmasq_upper <- 0.5
    sigma_gamma1 <- 1.5
    sigma_gamma2 <- 2.5
    lambda <- 0.5
    sigma_delta <- 1
    rho_gamma <- 0.5
    
    ## stan options
    niter <- 500
    nchains <- 2
    prop_warmup <- 0.5
    max_treedepth <- 25
    adapt_delta <- 0.95
    
    ## which run
    run_number <- 1
    
    ## which simulation
    sim <- 20
} else {
    ## which model to run
    model_number <- as.numeric(commandArgs(trailingOnly=TRUE)[1])
    model_to_run <- models_dat$model_name[models_dat$model_number == model_number]
    
    ## data generation options
    number_of_causes <- as.numeric(commandArgs(trailingOnly=TRUE)[2])
    number_of_regions <- as.numeric(commandArgs(trailingOnly=TRUE)[3])
    number_of_replications <- as.numeric(commandArgs(trailingOnly=TRUE)[4])
    
    ## parameters
    beta1 <- as.numeric(commandArgs(trailingOnly=TRUE)[5])
    beta2 <- as.numeric(commandArgs(trailingOnly=TRUE)[6])
    if (model_to_run == "2 cause FE, shared region IID RE v2") beta3 <- 0
    rho_lower <- as.numeric(commandArgs(trailingOnly=TRUE)[7])
    rho_upper <- as.numeric(commandArgs(trailingOnly=TRUE)[8])
    sigmasq_lower <- as.numeric(commandArgs(trailingOnly=TRUE)[9])
    sigmasq_upper <- as.numeric(commandArgs(trailingOnly=TRUE)[10])
    sigma_gamma1 <- as.numeric(commandArgs(trailingOnly=TRUE)[11])
    sigma_gamma2 <- as.numeric(commandArgs(trailingOnly=TRUE)[12])
    lambda <- as.numeric(commandArgs(trailingOnly=TRUE)[13])
    sigma_delta <- as.numeric(commandArgs(trailingOnly=TRUE)[14])
    rho_gamma <- as.numeric(commandArgs(trailingOnly=TRUE)[15])
    
    ## stan options
    niter <- as.numeric(commandArgs(trailingOnly=TRUE)[16])
    nchains <- as.numeric(commandArgs(trailingOnly=TRUE)[17])
    prop_warmup <- as.numeric(commandArgs(trailingOnly=TRUE)[18])
    max_treedepth <- as.numeric(commandArgs(trailingOnly=TRUE)[19])
    adapt_delta <- as.numeric(commandArgs(trailingOnly=TRUE)[20])
    
    ## which run
    run_number <- as.numeric(commandArgs(trailingOnly=TRUE)[21])
    
    ## which simulation
    sim <- as.numeric(commandArgs(trailingOnly=TRUE)[22])
}

# simulate data
simulated_data <- simulateData(R = number_of_regions, 
                               I = number_of_replications, 
                               C = number_of_causes,
                               beta = c(beta1, beta2), 
                               rho_lower = rho_lower, rho_upper = rho_upper, 
                               sigmasq_lower = sigmasq_lower, sigmasq_upper = sigmasq_upper,
                               sigma_gamma = c(sigma_gamma1, sigma_gamma2),
                               lambda = lambda,
                               sigma_delta = sigma_delta,
                               rho_gamma = rho_gamma,
                               seed = 1 + (sim * 2),
                               dgm = model_to_run)

# fit STAN model
stan_list <- fitSTAN(model_to_run, simulated_data$datlist,
                     niter = niter, nchains = nchains, nthin = 1, prop_warmup = prop_warmup,
                     max_treedepth = max_treedepth, adapt_delta = adapt_delta)

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

stan_diags <- data.frame(n_divergent = get_num_divergent(stan_list$mod_stan),
                         n_max_tree_exceeded = get_num_max_treedepth(stan_list$mod_stan),
                         n_bmfi_low_chains = sum(is.finite(get_low_bfmi_chains(stan_list$mod_stan))))

if (!testing) {
    # save results
    setwd(savedir)
    saveRDS(results, paste0("tmp/results_run-", run_number, "_sim-", sim, ".rds"))
    saveRDS(stan_diags, paste0("tmp/standiags_run-", run_number, "_sim-", sim, ".rds"))
}

if (testing) {
    # posterior estimates
    mod_summary <- summary(stan_list$mod_stan, 
                           pars = c("beta", "sigma_gamma", "lambda", "sigma_delta"),
                           probs = c(0.1, 0.9))
    mod_summary$summary
    
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
                          pars = c("alpha1"), probs = 0.5)$summary[, "50%"]
    alpha_hat2 <- summary(stan_list$mod_stan, 
                          pars = c("alpha2"), probs = 0.5)$summary[, "50%"]
    delta_hat <- summary(stan_list$mod_stan, 
                         pars = c("delta"), probs = 0.5)$summary[, "50%"]

    mean(alpha_hat1)
    mean(alpha_hat2)
    mean(delta_hat)
    par(mfrow = c(1, 3), mar = c(5, 4, 4, 2), oma = c(1, 1, .25, .25))
    plot((simulated_data$params$gamma_rc[,1] + simulated_data$params$beta[1]) ~ alpha_hat1,
         xlab = "Cause 1",
         ylab = "True alphas", 
         main = "Posterior median of the REs",
         pch = 19, col = alpha("dodgerblue", 0.5))
    abline(0, 1, col = "darkgreen")
    plot((simulated_data$params$gamma_rc[,2] + simulated_data$params$beta[2]) ~ alpha_hat2,
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
