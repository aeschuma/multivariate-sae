rm(list=ls())

# Preamble ####

## set the root depending on operating system
root <- ifelse(Sys.info()[1]=="Darwin","~/",
               ifelse(Sys.info()[1]=="Windows","P:/",
                      ifelse(Sys.info()[1]=="Linux","/home/users/aeschuma/",
                             stop("Unknown operating system"))))

## load libraries
library(mvtnorm); library(MASS);library(gtools); library(parallel);
library(scales); library(RColorBrewer); library(ggplot2); library(tidyverse);
library(haven); library(knitr); library(INLA); library(readr); library(magrittr);

## TESTING THE CODE?
testing <- FALSE

## define directories

# working directory for code
wd <- paste0(root,"Desktop/survey-csmf/unit-level-simulation")

# directory of the model results to load
resultsdir <- paste0(root,"Dropbox/dissertation_2/survey-csmf/results/ken2014-unit-level/hazwaz")

# directory to save results
savedir <- paste0(root,"Dropbox/dissertation_2/survey-csmf/results/unit-level-simulation")

## set directory
setwd(wd)

# set number of cores to the max possible
options(mc.cores = detectCores())

# start time for run time
start_time <- Sys.time()

# Load data and functions ####

# load simulation functions:
#   simulateData()
#   fitSTAN()
cat(paste("load sim functions \n"))
source("simulation-functions-unitlevel.R")

cat(paste("set parameters from command args \n"))

# Set parameters! ####
if (testing) {
    ## data generation options
    dgm <- 7
    
    ## which run
    run_number <- 15
    
    ## which simulation
    sim <- 2
} else {
    ## data generation options
    dgm <- as.numeric(commandArgs(trailingOnly=TRUE)[1])
    
    ## which run
    run_number <- as.numeric(commandArgs(trailingOnly=TRUE)[2])
    
    ## which simulation
    sim <- as.numeric(commandArgs(trailingOnly=TRUE)[3])
}

cat(paste("Set modeling info \n"))

# load data
dgm_file <- read_csv("dgm-info-unitlevel.csv")
my_dgm <- dgm_file %>% dplyr::filter(dgm_number == dgm) %>% dplyr::select(-dgm_number) %>% dplyr::slice(1)
load(my_dgm %>% pull(geo_data))
if (my_dgm$geo_data != my_dgm$sample_data) load(my_dgm %>% pull(sample_data))
random_re <- my_dgm %>% pull(random_re)
my_dgm <- my_dgm %>% dplyr::select(-geo_data, -sample_data, -random_re)

# get number of individuals to sample
sample_info <- dat %>% select(admin1, rural, cluster) %>% 
    group_by(admin1, rural, cluster) %>%
    summarize(n = n()) %>%
    ungroup() %>%
    arrange(admin1, rural, cluster)

cat(paste("Simulate data \n"))
# Simulate data ####
simulated_data <- simulateData(dgm_specs = my_dgm, 
                               sample_info = sample_info,
                               Amat = admin1.mat, 
                               scaling_factor = scaling_factor, 
                               seed_re = ifelse(random_re, sim + 500*run_number, 98125+run_number), 
                               seed_lik = run_number,
                               testing = FALSE)

# Fit INLA models ####
cat(paste("Fit INLA models \n"))
nsamps <- 100
inla.results <- fitAllINLAmodels(simulated_data = simulated_data, 
                                 nsamps = nsamps,
                                 seed = run_number + 1,
                                 testing = FALSE)

# Calculate summary measures ####
cat(paste("Extract summaries for parameters \n"))

parnames <- c("beta[1]", "beta[2]", 
              "gamma[1]", "gamma[2]", 
              "sigma[1]", "sigma[2]", 
              "rho[1]", "rho[2]", 
              "lambda",
              "sigma_epsilon[1]", "sigma_epsilon[2]",
              "gaussian_sd[1]", "gaussian_sd[2]")

# summaries
mod_summary <- tibble(param = parnames,
                      pct2.5 = NA,
                      pct10 = NA,
                      pct50 = NA,
                      pct90 = NA,
                      pct97.5 = NA)
mod_summary <- rep(list(mod_summary), length(inla.results))
latent_mean_summary <- vector(mode = "list", length = length(inla.results))

for (i in 1:length(inla.results))  {
    # message(paste0("Model ", names(mod_summary)[i]))
    hyperpar_names <- inla.results[[i]]$fit$summary.hyperpar %>% rownames()
    
    # message("FEs")
    mod_summary[[i]][mod_summary[[i]]$param %in% c("beta[1]", "beta[2]", "gamma[1]", "gamma[2]"), c("pct2.5", "pct10", "pct50","pct90","pct97.5")] <- inla.results[[i]]$fit$summary.fixed[, c("0.025quant", "0.1quant", "0.5quant", "0.9quant", "0.975quant")]

    # message("REs")
    for(j in 1:length(hyperpar_names)) {
        tmpname <- hyperpar_names[j]
        if (tmpname == "Precision for admin1.haz") {
            mod_summary[[i]][mod_summary[[i]]$param == "sigma[1]", c("pct2.5", "pct10", "pct50","pct90","pct97.5")] <- rev(inla.results[[i]]$fit$summary.hyperpar[tmpname, c("0.025quant", "0.1quant", "0.5quant", "0.9quant", "0.975quant")]^(-0.5))
        } else if (tmpname == "Precision for admin1.waz") {
            mod_summary[[i]][mod_summary[[i]]$param == "sigma[2]", c("pct2.5", "pct10", "pct50","pct90","pct97.5")] <- rev(inla.results[[i]]$fit$summary.hyperpar[tmpname, c("0.025quant", "0.1quant", "0.5quant", "0.9quant", "0.975quant")]^(-0.5))
        } else if (tmpname == "Phi for admin1.haz") {
            mod_summary[[i]][mod_summary[[i]]$param == "rho[1]", c("pct2.5", "pct10", "pct50","pct90","pct97.5")] <- inla.results[[i]]$fit$summary.hyperpar[tmpname, c("0.025quant", "0.1quant", "0.5quant", "0.9quant", "0.975quant")]
        } else if (tmpname == "Phi for admin1.waz") {
            mod_summary[[i]][mod_summary[[i]]$param == "rho[2]", c("pct2.5", "pct10", "pct50","pct90","pct97.5")] <- inla.results[[i]]$fit$summary.hyperpar[tmpname, c("0.025quant", "0.1quant", "0.5quant", "0.9quant", "0.975quant")]
        } else if (grepl("Beta", tmpname)) {
            mod_summary[[i]][mod_summary[[i]]$param == "lambda", c("pct2.5", "pct10", "pct50","pct90","pct97.5")] <- inla.results[[i]]$fit$summary.hyperpar[tmpname, c("0.025quant", "0.1quant", "0.5quant", "0.9quant", "0.975quant")]
        } else if (tmpname == "Precision for cluster.haz") {
            mod_summary[[i]][mod_summary[[i]]$param == "sigma_epsilon[1]", c("pct2.5", "pct10", "pct50","pct90","pct97.5")] <- rev(inla.results[[i]]$fit$summary.hyperpar[tmpname, c("0.025quant", "0.1quant", "0.5quant", "0.9quant", "0.975quant")]^(-0.5))
        } else if (tmpname == "Precision for cluster.waz") {
            mod_summary[[i]][mod_summary[[i]]$param == "sigma_epsilon[2]", c("pct2.5", "pct10", "pct50","pct90","pct97.5")] <- rev(inla.results[[i]]$fit$summary.hyperpar[tmpname, c("0.025quant", "0.1quant", "0.5quant", "0.9quant", "0.975quant")]^(-0.5))
        } else if (tmpname == "Precision for the Gaussian observations") {
            mod_summary[[i]][mod_summary[[i]]$param == "gaussian_sd[1]", c("pct2.5", "pct10", "pct50","pct90","pct97.5")] <- rev(inla.results[[i]]$fit$summary.hyperpar[tmpname, c("0.025quant", "0.1quant", "0.5quant", "0.9quant", "0.975quant")]^(-0.5))
        } else if (tmpname == "Precision for the Gaussian observations[2]") {
            mod_summary[[i]][mod_summary[[i]]$param == "gaussian_sd[2]", c("pct2.5", "pct10", "pct50","pct90","pct97.5")] <- rev(inla.results[[i]]$fit$summary.hyperpar[tmpname, c("0.025quant", "0.1quant", "0.5quant", "0.9quant", "0.975quant")]^(-0.5))
        }
    }
    latent_mean_summary[[i]] <- inla.results[[i]]$latent_means %>%
        mutate(param = paste0(outcome, "_", admin1, ".",  rural))
}

# true latent means
real_means <- simulated_data$latent_params %>% 
    select(admin1, rural, admin1_rural_HAZ, admin1_rural_WAZ) %>%
    distinct() %>%
    pivot_longer(cols = starts_with("admin1_rural_"),
                 names_prefix = "admin1_rural_",
                 names_to = "outcome",
                 values_to = "real_mean") %>%
    mutate(param = paste0(outcome, "_", admin1, ".",  rural))

cat(paste("Calculate results for parameters \n"))

# compile results
sim_res_latent_means <- tibble(param = real_means$param,
                               truth = real_means$real_mean)
sim_res_latent_means <- rep(list(sim_res_latent_means), length(inla.results))

for (i in 1:length(sim_res_latent_means)) {
    sim_res_latent_means[[i]] %<>% 
        left_join(latent_mean_summary[[i]]) %>%
        mutate(est = `est_50%`,
               coverage.80 = (truth > `est_10%`) & (truth < `est_90%`),
               width.80 = `est_90%` - `est_10%`,
               coverage.95 = (truth > `est_2.5%`) & (truth < `est_97.5%`),
               width.95 = `est_97.5%` - `est_2.5%`) %>%
        select(param, truth, est, 
               coverage.80, width.80, 
               coverage.95, width.95)
}

# model parameters
sim_res_par <- tibble(param = parnames,
                      truth = simulated_data$params,
                      est = NA,
                      coverage.80 = NA,
                      width.80 = NA,
                      coverage.95 = NA,
                      width.95 = NA)
sim_res_par <- rep(list(sim_res_par), length(inla.results))

for (i in 1:length(sim_res_par)) {
    sim_res_par[[i]] %<>% 
        left_join(mod_summary[[i]]) %>%
        mutate(est = `pct50`,
               coverage.80 = (truth > `pct10`) & (truth < `pct90`),
               width.80 = `pct90` - `pct10`,
               coverage.95 = (truth > `pct2.5`) & (truth < `pct97.5`),
               width.95 = `pct97.5` - `pct2.5`) %>%
        select(param, truth, est, 
               coverage.80, width.80, 
               coverage.95, width.95)
}

sim_res <- list()
for (i in 1:length(sim_res_par)) {
    sim_res[[i]] <- rbind(sim_res_par[[i]],
                          sim_res_latent_means[[i]])
}

sim_mod_res$run_time <- Sys.time() - start_time

# Save results ####
cat(paste("Save results \n"))
if (!testing) {
    # save results
    setwd(savedir)
    saveRDS(sim_res, paste0("tmp/results_run-", run_number, "_sim-", sim, ".rds"))
}
