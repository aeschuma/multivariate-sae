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
    dgm <- 1
    
    ## which run
    run_number <- 1
    
    ## which simulation
    sim <- 1
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
                               seed_lik = sim,
                               testing = FALSE)

# Fit INLA models ####
cat(paste("Fit INLA models \n"))
nsamps <- 100
inla.results <- fitAllINLAmodels(simulated_data = simulated_data, 
                                 nsamps = nsamps,
                                 testing = FALSE)

# Calculate summary measures ####
cat(paste("Extract summaries for parameters \n"))

parnames <- c("beta[1]", "beta[2]", 
              "gamma[1]", "gamma[2]", 
              "sigma[1]", "sigma[2]", 
              "rho[1]", "rho[2]", 
              "lambda",
              "sigma_epsilon[1]", "sigma_epsilon[2]")

# summaries
mod_summary <- tibble(params = parnames,
                      pct2.5 = NA,
                      pct10 = NA,
                      pct50 = NA,
                      pct90 = NA,
                      pct97.5 = NA,
                      var = NA)
mod_summary <- rep(list(mod_summary), length(inla.results))
latent_mean_summary <- vector(mode = "list", length = length(inla.results))
for (i in 1:length(inla.results))  {
    mod_summary[[i]][mod_summary[[i]]$params %in% c("beta[1]", "beta[2]", "gamma[1]", "gamma[2]"), c("pct2.5", "pct10", "pct50","pct90","pct97.5")] <- inla.results[[i]]$fit$summary.fixed[, c("0.025quant", "0.1quant", "0.5quant", "0.9quant", "0.975quant")]
    mod_summary[[i]][mod_summary[[i]]$params %in% c("beta[1]", "beta[2]", "gamma[1]", "gamma[2]"), c("var")] <- inla.results[[i]]$fit$summary.fixed[, c("sd")]^2
    hyperpar_names <- rownames(inla.results[[i]]$fit$summary.hyperpar)
    for(j in 1:length(hyperpar_names)) {
        tmpname <- hyperpar_names[j]
        if (str_detect(tmpname, "Precision")) {
            prec_par <- which(names(inla.results[[i]]$fit$internal.marginals.hyperpar) == gsub("Precision", "Log precision", tmpname))
            m.sd <- inla.tmarginal(function(x) sqrt(1/exp(x)),
                                   inla.results[[i]]$fit$internal.marginals.hyperpar[[prec_par]])
            marg <- inla.zmarginal(m.sd, silent = TRUE)
        }
        if (tmpname == "Precision for admin1.haz") {
            mod_summary[[i]][mod_summary[[i]]$params == "sigma[1]", c("pct2.5", "pct10", "pct50","pct90","pct97.5")] <- rev(inla.results[[i]]$fit$summary.hyperpar[tmpname, c("0.025quant", "0.1quant", "0.5quant", "0.9quant", "0.975quant")]^(-0.5))
            mod_summary[[i]][mod_summary[[i]]$params == "sigma[1]", "var"] <- marg$sd^2
        } else if (tmpname == "Precision for admin1.waz") {
            mod_summary[[i]][mod_summary[[i]]$params == "sigma[2]", c("pct2.5", "pct10", "pct50","pct90","pct97.5")] <- rev(inla.results[[i]]$fit$summary.hyperpar[tmpname, c("0.025quant", "0.1quant", "0.5quant", "0.9quant", "0.975quant")]^(-0.5))
            mod_summary[[i]][mod_summary[[i]]$params == "sigma[2]", "var"] <- marg$sd^2
        } else if (tmpname == "Phi for admin1.haz") {
            mod_summary[[i]][mod_summary[[i]]$params == "rho[1]", c("pct2.5", "pct10", "pct50","pct90","pct97.5")] <- inla.results[[i]]$fit$summary.hyperpar[tmpname, c("0.025quant", "0.1quant", "0.5quant", "0.9quant", "0.975quant")]
            mod_summary[[i]][mod_summary[[i]]$params == "rho[1]", "var"] <- inla.results[[i]]$fit$summary.hyperpar[tmpname, "sd"]^2
        } else if (tmpname == "Phi for admin1.waz") {
            mod_summary[[i]][mod_summary[[i]]$params == "rho[2]", c("pct2.5", "pct10", "pct50","pct90","pct97.5")] <- inla.results[[i]]$fit$summary.hyperpar[tmpname, c("0.025quant", "0.1quant", "0.5quant", "0.9quant", "0.975quant")]
            mod_summary[[i]][mod_summary[[i]]$params == "rho[2]", "var"] <- inla.results[[i]]$fit$summary.hyperpar[tmpname, "sd"]^2
        } else if (grepl("Beta", tmpname)) {
            mod_summary[[i]][mod_summary[[i]]$params == "lambda", c("pct2.5", "pct10", "pct50","pct90","pct97.5")] <- inla.results[[i]]$fit$summary.hyperpar[tmpname, c("0.025quant", "0.1quant", "0.5quant", "0.9quant", "0.975quant")]
            mod_summary[[i]][mod_summary[[i]]$params == "lambda", "var"] <- inla.results[[i]]$fit$summary.hyperpar[tmpname, "sd"]^2
        } else if (tmpname == "Precision for cluster.haz") {
            mod_summary[[i]][mod_summary[[i]]$params == "sigma_epsilon[1]", c("pct2.5", "pct10", "pct50","pct90","pct97.5")] <- rev(inla.results[[i]]$fit$summary.hyperpar[tmpname, c("0.025quant", "0.1quant", "0.5quant", "0.9quant", "0.975quant")]^(-0.5))
            mod_summary[[i]][mod_summary[[i]]$params == "sigma_epsilon[1]", "var"] <- marg$sd^2
        } else if (tmpname == "Precision for cluster.waz") {
            mod_summary[[i]][mod_summary[[i]]$params == "sigma_epsilon[2]", c("pct2.5", "pct10", "pct50","pct90","pct97.5")] <- rev(inla.results[[i]]$fit$summary.hyperpar[tmpname, c("0.025quant", "0.1quant", "0.5quant", "0.9quant", "0.975quant")]^(-0.5))
            mod_summary[[i]][mod_summary[[i]]$params == "sigma_epsilon[2]", "var"] <- marg$sd^2
        } else if (tmpname == "Precision for the Gaussian observations") {
            mod_summary[[i]][mod_summary[[i]]$params == "gaussian_sd[1]", c("pct2.5", "pct10", "pct50","pct90","pct97.5")] <- rev(inla.results[[i]]$fit$summary.hyperpar[tmpname, c("0.025quant", "0.1quant", "0.5quant", "0.9quant", "0.975quant")]^(-0.5))
            mod_summary[[i]][mod_summary[[i]]$params == "gaussian_sd[1]", "var"] <- marg$sd^2
        } else if (tmpname == "Precision for the Gaussian observations[2]") {
            mod_summary[[i]][mod_summary[[i]]$params == "gaussian_sd[2]", c("pct2.5", "pct10", "pct50","pct90","pct97.5")] <- rev(inla.results[[i]]$fit$summary.hyperpar[tmpname, c("0.025quant", "0.1quant", "0.5quant", "0.9quant", "0.975quant")]^(-0.5))
            mod_summary[[i]][mod_summary[[i]]$params == "gaussian_sd[2]", "var"] <- marg$sd^2
        } 
    }
    latent_mean_summary[[i]] <- inla.results[[i]]$latent_means
}

cat(paste("Calculate results for parameters \n"))

# compile results
sim_res <- data.frame(param = c(parnames, "HAZ latent means", "WAZ latent means"),
                      bias = NA,
                      absolute_bias = NA,
                      relative_absolute_bias = NA,
                      variance = NA,
                      mse = NA,
                      coverage.80 = NA,
                      width.80 = NA,
                      coverage.95 = NA,
                      width.95 = NA)
sim_res <- rep(list(sim_res), length(inla.results))
names(sim_res) <- names(inla.results)

# true latent means
real_means <- simulated_data$latent_params %>% 
    select(admin1, rural, admin1_rural_HAZ, admin1_rural_WAZ) %>%
    distinct() %>%
    pivot_longer(cols = starts_with("admin1_rural_"),
                 names_prefix = "admin1_rural_",
                 names_to = "outcome",
                 values_to = "real_mean")

for (i in 1:length(sim_res)) {
    for (j in 1:nrow(sim_res[[i]])) {
        tmp_param_name <- sim_res[[i]]$param[j]
        if (grepl("latent means", tmp_param_name)) {
            restmp2 <- real_means %>% 
                left_join(latent_mean_summary[[i]])
            if (grepl("HAZ", tmp_param_name)) {
                restmp <- restmp2 %>% 
                    filter(outcome == "HAZ") 
            } else {
                restmp <- restmp2 %>% 
                    filter(outcome == "WAZ") 
            }    
            sim_res[[i]]$bias[sim_res[[i]]$param == tmp_param_name] <- mean(restmp$`est_50%` - restmp$real_mean)
            sim_res[[i]]$absolute_bias[sim_res[[i]]$param == tmp_param_name] <- mean(abs(restmp$`est_50%` - restmp$real_mean))
            sim_res[[i]]$relative_absolute_bias[sim_res[[i]]$param == tmp_param_name] <- mean(abs(restmp$`est_50%` - restmp$real_mean)/abs(restmp$real_mean))
            sim_res[[i]]$variance[sim_res[[i]]$param == tmp_param_name] <- mean(restmp$var)
            sim_res[[i]]$mse[sim_res[[i]]$param == tmp_param_name] <- (sim_res[[i]]$bias[sim_res[[i]]$param == tmp_param_name])^2 + sim_res[[i]]$variance[sim_res[[i]]$param == tmp_param_name]
            sim_res[[i]]$coverage.80[sim_res[[i]]$param == tmp_param_name] <- mean((restmp$real_mean > restmp$`est_10%`) & (restmp$real_mean < restmp$`est_90%`))
            sim_res[[i]]$width.80[sim_res[[i]]$param == tmp_param_name] <- mean(restmp$`est_90%` - restmp$`est_10%`)
            sim_res[[i]]$coverage.95[sim_res[[i]]$param == tmp_param_name] <- mean((restmp$real_mean > restmp$`est_2.5%`) & (restmp$real_mean < restmp$`est_97.5%`))
            sim_res[[i]]$width.95[sim_res[[i]]$param == tmp_param_name] <- mean(restmp$`est_97.5%` - restmp$`est_2.5%`)
        } else {
            tmp_param <- simulated_data$params[tmp_param_name]
            idx <- mod_summary[[i]]$params == tmp_param_name
            
            sim_res[[i]]$bias[sim_res[[i]]$param == tmp_param_name] <- mod_summary[[i]]$pct50[idx] - tmp_param
            sim_res[[i]]$absolute_bias[sim_res[[i]]$param == tmp_param_name] <- abs(mod_summary[[i]]$pct50[idx] - tmp_param)
            sim_res[[i]]$relative_absolute_bias[sim_res[[i]]$param == tmp_param_name] <- abs(mod_summary[[i]]$pct50[idx] - tmp_param)/abs(tmp_param)
            sim_res[[i]]$variance[sim_res[[i]]$param == tmp_param_name] <- mod_summary[[i]]$var[idx] - tmp_param
            sim_res[[i]]$mse[sim_res[[i]]$param == tmp_param_name] <- (sim_res[[i]]$bias[sim_res[[i]]$param == tmp_param_name])^2 + sim_res[[i]]$var[sim_res[[i]]$param == tmp_param_name]
            sim_res[[i]]$coverage.80[sim_res[[i]]$param == tmp_param_name] <- (tmp_param > mod_summary[[i]]$pct10[idx]) & (tmp_param < mod_summary[[i]]$pct90[idx])
            sim_res[[i]]$width.80[sim_res[[i]]$param == tmp_param_name] <- mod_summary[[i]]$pct90[idx] - mod_summary[[i]]$pct10[idx]
            sim_res[[i]]$coverage.95[sim_res[[i]]$param == tmp_param_name] <- (tmp_param > mod_summary[[i]]$pct2.5[idx]) & (tmp_param < mod_summary[[i]]$pct97.5[idx])
            sim_res[[i]]$width.95[sim_res[[i]]$param == tmp_param_name] <- mod_summary[[i]]$pct97.5[idx] - mod_summary[[i]]$pct2.5[idx]
        }
    }
}

sim_res$run_time <- Sys.time() - start_time

# Save results ####
cat(paste("Save results \n"))
if (!testing) {
    # save results
    setwd(savedir)
    saveRDS(sim_res, paste0("tmp/results_run-", run_number, "_sim-", sim, ".rds"))
}
