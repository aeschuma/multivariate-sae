rm(list=ls())

# Preamble ####

## set the root depending on operating system
root <- ifelse(Sys.info()[1]=="Darwin","~/",
               ifelse(Sys.info()[1]=="Windows","P:/",
                      ifelse(Sys.info()[1]=="Linux","/home/users/aeschuma/",
                             stop("Unknown operating system"))))

## load libraries
library(mvtnorm); library(MASS);library(gtools); library(parallel);
library(scales); library(RColorBrewer); library(ggplot2); library(tidyverse); library(dplyr);
library(haven); library(knitr); library(INLA); library(readr);

## TESTING THE CODE?
testing <- FALSE

## define directories

# working directory for code
wd <- paste0(root,"Desktop/survey-csmf/stage-2-simulation")

# directory of the model results to load
resultsdir <- paste0(root,"Dropbox/dissertation_2/survey-csmf/results/ken2014-hazwaz/")

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

cat(paste("set parameters from command args \n"))

# Set parameters! ####
if (testing) {
    ## data generation options
    dgm <- 10
    
    ## which run
    run_number <- 999
    
    ## which simulation
    sim <- 499
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
dgm_file <- read_csv("dgm-info.csv")
my_dgm <- dgm_file %>% dplyr::filter(dgm_number == dgm) %>% dplyr::select(-dgm_number) %>% dplyr::slice(1)
load(my_dgm %>% pull(geo_data))
random_re <- my_dgm %>% pull(random_re)
my_dgm <- my_dgm %>% dplyr::select(-geo_data, -random_re)

cat(paste("Simulate data \n"))
# Simulate data ####
simulated_data <- simulateData(dgm_specs = my_dgm, 
                               n_r = table(dat$admin1)/(1.5^2),
                               Amat = admin1.mat, 
                               scaling_factor = scaling_factor, 
                               seed_re = ifelse(random_re, sim + 500*run_number, 98125+run_number), 
                               seed_lik = sim,
                               testing = FALSE)

# Fit INLA models ####
cat(paste("Fit INLA models \n"))
inla.results <- fitAllINLAmodels(simulated_data = simulated_data, testing = FALSE)

# Calculate summary measures ####
cat(paste("Extract summaries for parameters \n"))

parnames <- c("beta[1]", "beta[2]", "sigma[1]", "sigma[2]", "rho[1]", "rho[2]", "lambda")

# summaries
mod_summary <- tibble(params = parnames,
                      sd = NA,
                      pct2.5 = NA,
                      pct10 = NA,
                      pct50 = NA,
                      pct90 = NA,
                      pct97.5 = NA)
mod_summary <- rep(list(mod_summary), length(inla.results))
latent_mean_summary <- vector(mode = "list", length = length(inla.results))
for (i in 1:length(inla.results))  {
    mod_summary[[i]][mod_summary[[i]]$params %in% c("beta[1]", "beta[2]"), c("sd", "pct2.5", "pct10", "pct50","pct90","pct97.5")] <- inla.results[[i]]$summary.fixed[, c("sd", "0.025quant", "0.1quant", "0.5quant", "0.9quant", "0.975quant")]
    hyperpar_names <- rownames(inla.results[[i]]$summary.hyperpar)
    for(j in 1:length(hyperpar_names)) {
        tmpname <- hyperpar_names[j]
        if (tmpname == "Precision for admin1.haz") {
            mod_summary[[i]][mod_summary[[i]]$params == "sigma[1]", c("sd", "pct2.5", "pct10", "pct50","pct90","pct97.5")] <- rev(inla.results[[i]]$summary.hyperpar[tmpname, c("sd", "0.025quant", "0.1quant", "0.5quant", "0.9quant", "0.975quant")]^(-0.5))
        } else if (tmpname == "Precision for admin1.waz") {
            mod_summary[[i]][mod_summary[[i]]$params == "sigma[2]", c("sd", "pct2.5", "pct10", "pct50","pct90","pct97.5")] <- rev(inla.results[[i]]$summary.hyperpar[tmpname, c("sd", "0.025quant", "0.1quant", "0.5quant", "0.9quant", "0.975quant")]^(-0.5))
        } else if (tmpname == "Phi for admin1.haz") {
            mod_summary[[i]][mod_summary[[i]]$params == "rho[1]", c("sd", "pct2.5", "pct10", "pct50","pct90","pct97.5")] <- inla.results[[i]]$summary.hyperpar[tmpname, c("sd", "0.025quant", "0.1quant", "0.5quant", "0.9quant", "0.975quant")]
        } else if (tmpname == "Phi for admin1.waz") {
            mod_summary[[i]][mod_summary[[i]]$params == "rho[2]", c("sd", "pct2.5", "pct10", "pct50","pct90","pct97.5")] <- inla.results[[i]]$summary.hyperpar[tmpname, c("sd", "0.025quant", "0.1quant", "0.5quant", "0.9quant", "0.975quant")]
        } else if (grepl("Beta", tmpname)) {
            mod_summary[[i]][mod_summary[[i]]$params == "lambda", c("sd", "pct2.5", "pct10", "pct50","pct90","pct97.5")] <- inla.results[[i]]$summary.hyperpar[tmpname, c("sd", "0.025quant", "0.1quant", "0.5quant", "0.9quant", "0.975quant")]
        } 
    }
    latent_mean_summary[[i]] <- inla.results[[i]]$summary.lincomb.derived[, c("sd", "0.025quant", "0.1quant", "0.5quant", "0.9quant", "0.975quant")]
}

cat(paste("Calculate results for parameters \n"))

true_vals <- c(simulated_data$params$bivariate_means[,1], simulated_data$params$bivariate_means[,2])

# compile results
sim_res_latent_means <- data.frame(param = rownames(latent_mean_summary[[1]])[!grepl("outcome", rownames(latent_mean_summary[[1]]))],
                                   truth = NA,
                                   est = NA,
                                   coverage.80 = NA,
                                   width.80 = NA,
                                   coverage.95 = NA,
                                   width.95 = NA)
sim_res_latent_means <- rep(list(sim_res_latent_means), length(inla.results))

for (i in 1:length(sim_res_latent_means)) {
    lms <- latent_mean_summary[[i]][!grepl("outcome", rownames(latent_mean_summary[[i]])),]
    sim_res_latent_means[[i]]$truth <- true_vals
    sim_res_latent_means[[i]]$est <- lms$`0.5quant`
    sim_res_latent_means[[i]]$coverage.80 <- (true_vals > lms$`0.1quant`) & (true_vals < lms$`0.9quant`)
    sim_res_latent_means[[i]]$width.80 <- lms$`0.9quant` - lms$`0.1quant`
    sim_res_latent_means[[i]]$coverage.95 <- (true_vals > lms$`0.025quant`) & (true_vals < lms$`0.975quant`)
    sim_res_latent_means[[i]]$width.95 <- lms$`0.975quant` - lms$`0.025quant`
}

sim_res_par <- data.frame(param = parnames,
                          truth = NA,
                          est = NA,
                          coverage.80 = NA,
                          width.80 = NA,
                          coverage.95 = NA,
                          width.95 = NA)
sim_res_par <- rep(list(sim_res_par), length(inla.results))

for (i in 1:length(sim_res_par)) {
    for (j in 1:nrow(sim_res_par[[i]])) {
        tmp_param_name <- sim_res_par[[i]]$param[j]
        tmp_param <- simulated_data$params[["mean_pars"]][tmp_param_name]
        idx <- mod_summary[[i]]$params == tmp_param_name
        
        sim_res_par[[i]]$truth[sim_res_par[[i]]$param == tmp_param_name] <- tmp_param
        sim_res_par[[i]]$est[sim_res_par[[i]]$param == tmp_param_name] <- mod_summary[[i]]$pct50[idx]
        sim_res_par[[i]]$coverage.80[sim_res_par[[i]]$param == tmp_param_name] <- (tmp_param > mod_summary[[i]]$pct10[idx]) & (tmp_param < mod_summary[[i]]$pct90[idx])
        sim_res_par[[i]]$width.80[sim_res_par[[i]]$param == tmp_param_name] <- mod_summary[[i]]$pct90[idx] - mod_summary[[i]]$pct10[idx]
        sim_res_par[[i]]$coverage.95[sim_res_par[[i]]$param == tmp_param_name] <- (tmp_param > mod_summary[[i]]$pct2.5[idx]) & (tmp_param < mod_summary[[i]]$pct97.5[idx])
        sim_res_par[[i]]$width.95[sim_res_par[[i]]$param == tmp_param_name] <- mod_summary[[i]]$pct97.5[idx] - mod_summary[[i]]$pct2.5[idx]
    }
}

sim_mod_res <- list()
for (i in 1:length(sim_res_par)) {
    sim_mod_res[[i]] <- rbind(sim_res_par[[i]],
                          sim_res_latent_means[[i]])
}

# Direct estimates ####
direst_values <- simulated_data$datlist$value
std_errs_dir <- c(simulated_data$dattibble$seHAZ.bi, simulated_data$dattibble$seWAZ.bi)
ci_80_direct <- cbind(lower = direst_values - qnorm(0.9)*std_errs_dir, 
                      upper = direst_values + qnorm(0.9)*std_errs_dir) %>% as.data.frame()
ci_95_direct <- cbind(lower = direst_values - qnorm(0.975)*std_errs_dir, 
                      upper = direst_values + qnorm(0.975)*std_errs_dir) %>% as.data.frame()

direst <- data.frame(param = sim_res_latent_means[[1]]$param,
                     truth = true_vals,
                     est = direst_values,
                     coverage.80 = (true_vals > ci_80_direct$lower) & (true_vals < ci_80_direct$upper),
                     width.80 = ci_80_direct$upper - ci_80_direct$lower,
                     coverage.95 = (true_vals > ci_95_direct$lower) & (true_vals < ci_95_direct$upper),
                     width.95 = ci_95_direct$upper - ci_95_direct$lower)

sim_res <- c(list(direst), sim_mod_res)
names(sim_res) <- c("Direct estimates", names(inla.results))

# Save results ####
cat(paste("Save results \n"))
if (!testing) {
    # save results
    setwd(savedir)
    saveRDS(sim_res, paste0("tmp/results_run-", run_number, "_sim-", sim, ".rds"))
}
