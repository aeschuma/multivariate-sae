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
library(ggplot2); library(tidyverse);

########
## TESTING THE CODE?
########

testing <- TRUE

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
models_dat <- read_csv("model-info.csv")

## set parameters!
if (testing) {
    ## which model to run
    model_number <- 1
    model_to_run <- models_dat$model_name[model_number]
    
    ## data generation options
    number_of_causes <- 2
    number_of_regions <- 8
    number_of_replications <- 1
    
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
    rho_gamma <- 0.25
    
    ## stan options
    niter <- 5000
    nchains <- 3
    prop_warmup <- 0.5
    max_treedepth <- 15
    adapt_delta <- 0.8
    
    ## which simulation
    sim <- 3
} else {
    ## which model to run
    model_number <- as.numeric(commandArgs(trailingOnly=TRUE)[1])
    model_to_run <- models_dat$model_name[model_number]
    
    ## data generation options
    number_of_causes <- as.numeric(commandArgs(trailingOnly=TRUE)[2])
    number_of_regions <- as.numeric(commandArgs(trailingOnly=TRUE)[3])
    number_of_replications <- as.numeric(commandArgs(trailingOnly=TRUE)[4])
    
    ## parameters
    beta1 <- as.numeric(commandArgs(trailingOnly=TRUE)[5])
    beta2 <- as.numeric(commandArgs(trailingOnly=TRUE)[6])
    rho_lower <- as.numeric(commandArgs(trailingOnly=TRUE)[7])
    rho_upper <- as.numeric(commandArgs(trailingOnly=TRUE)[8])
    sigmasq_lower <- as.numeric(commandArgs(trailingOnly=TRUE)[9])
    sigmasq_upper <- as.numeric(commandArgs(trailingOnly=TRUE)[10])
    sigma_gamma1 <- as.numeric(commandArgs(trailingOnly=TRUE)[11])
    sigma_gamma2 <- as.numeric(commandArgs(trailingOnly=TRUE)[12])
    lambda <- as.numeric(commandArgs(trailingOnly=TRUE)[13])
    rho_gamma <- as.numeric(commandArgs(trailingOnly=TRUE)[14])
    
    ## stan options
    niter <- as.numeric(commandArgs(trailingOnly=TRUE)[15])
    nchains <- as.numeric(commandArgs(trailingOnly=TRUE)[16])
    prop_warmup <- as.numeric(commandArgs(trailingOnly=TRUE)[17])
    max_treedepth <- as.numeric(commandArgs(trailingOnly=TRUE)[18])
    adapt_delta <- as.numeric(commandArgs(trailingOnly=TRUE)[19])
    
    ## which simulation
    sim <- as.numeric(commandArgs(trailingOnly=TRUE)[20])
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
                               seed = 1 + (sim * 2),
                               dgm = model_to_run)

# fit STAN model
stan_list <- fitSTAN(model_to_run, simulated_data$datlist,
                     niter = niter, nchains = nchains, nthin = 1, prop_warmup = prop_warmup,
                     max_treedepth = max_treedepth, adapt_delta = adapt_delta)

# save parameter names in order to extract and save results
param_names <- names(simulated_data$params)
params_to_extract <- param_names[!(param_names %in% grep("gamma_", param_names, value = TRUE))]

# summaries
mod_summary <- summary(stan_list$mod_stan,
                       pars = params_to_extract,
                       probs = c(0.1, 0.5, 0.9))
posterior_qs <- mod_summary$summary[, c("10%", "50%", "90%")]
extracted_params <- rownames(posterior_qs)

# truth
beta <- c(beta1, beta2)
sigma_gamma <- c(sigma_gamma1, sigma_gamma2)
absolute_bias <- relative_bias <- coverage <- width <- vector(mode = "list", length = length(params_to_extract))
for (i in 1:length(params_to_extract)) {
    tmp_param_name <- params_to_extract[i]
    tmp_param <- get(tmp_param_name)
    if (tmp_param_name == "beta") {
        beta_est <- posterior_qs[c("beta[1]", "beta[2]"), "50%"]
        beta_ci <- posterior_qs[c("beta[1]", "beta[2]"), c("10%", "90%")]
        absolute_bias[[i]] <- beta_est - tmp_param
        relative_bias[[i]] <- (beta_est - tmp_param)/tmp_param
        coverage[[i]]  <- (tmp_param > beta_ci[, "10%"]) & (tmp_param < beta_ci[, "90%"])
        width[[i]] <- beta_ci[, "90%"] - beta_ci[, "10%"]
    } else if (tmp_param_name == "sigma_gamma") {
        sigma_gamma_est <- posterior_qs[c("sigma_gamma[1]", "sigma_gamma[2]"), "50%"]
        sigma_gamma_ci <- posterior_qs[c("sigma_gamma[1]", "sigma_gamma[2]"), c("10%", "90%")]
        absolute_bias[[i]] <- sigma_gamma_est - tmp_param
        relative_bias[[i]] <- (sigma_gamma_est - tmp_param)/tmp_param
        coverage[[i]]  <- (tmp_param > sigma_gamma_ci[, "10%"]) & (tmp_param < sigma_gamma_ci[, "90%"])
        width[[i]] <- sigma_gamma_ci[, "90%"] - sigma_gamma_ci[, "10%"]
    } else {
        absolute_bias[[i]] <- posterior_qs[tmp_param_name, "50%"] - tmp_param
        relative_bias[[i]] <- (posterior_qs[tmp_param_name, "50%"] - tmp_param)/tmp_param
        coverage[[i]]  <- (tmp_param > posterior_qs[tmp_param_name, "10%"]) & (tmp_param < posterior_qs[tmp_param_name, "90%"])
        width[[i]] <- posterior_qs[tmp_param_name, "90%"] - posterior_qs[tmp_param_name, "10%"]
    }
}
results <- list(absolute_bias, relative_bias, coverage, width)
# save results
setwd(savedir)
save_id <- paste(model_number, sep = "-")
saveRDS(results, paste0("tmp/res-", ""))