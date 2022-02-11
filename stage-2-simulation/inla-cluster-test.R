rm(list=ls())

# Preamble ####

## set the root depending on operating system
root <- ifelse(Sys.info()[1]=="Darwin","~/",
               ifelse(Sys.info()[1]=="Windows","P:/",
                      ifelse(Sys.info()[1]=="Linux","/home/users/aeschuma/",
                             stop("Unknown operating system"))))

library(INLA)
library(tidyverse)

# directory to save results
savedir <- paste0(root,"Dropbox/dissertation_2/survey-csmf/results/stage-2-simulation")

testing <- TRUE 

# Set parameters! ####
if (testing) {
## data generation options
    dgm <- 1
    
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

set.seed(run_number + 10000)

# number of regions
n_regions <- 47

# extract values from data
n_r <- rep(400, n_regions)
d <- 1.5
n_r_star <- n_r/(d^2)
V_r <- rep(0.1, n_regions)

# extract parameters from model for sim
beta_1 <- -1
sigma <- 0.25

# simulate a single set of random effects
set.seed(8008135)
v_r <- rnorm(n_regions, 0, sigma)

# calculate latent means
mu_r <- beta_1 + v_r

# what measures do we want to calculate
measures <- c("mean bias",
              "mean absolute bias",
              "mean relative bias",
              "80% coverage",
              "95% coverage",
              "80% width",
              "95% width")

# start simulations
set.seed(80085)
n_sim <- 2
sim_results <- list(direct_estimates = NA,
                    univariate_iid = NA)
for (i in 1:length(sim_results)) {
    sim_results[[i]] <- array(NA, dim = c(length(measures), n_sim, n_regions))
    # sim_results[[i]] <- tibble(bias = rep(NA, n_sim),
    #                           abs_bias = rep(NA, n_sim),
    #                           rel_abs_bias = rep(NA, n_sim),
    #                           coverage_80 = rep(NA, n_sim),
    #                           width_80 = rep(NA, n_sim),
    #                           coverage_95 = rep(NA, n_sim),
    #                           width_95 = rep(NA, n_sim))
}

for(s in 1:n_sim) {
    # simulate direct estimates
    y_hat_r <- rnorm(n_regions, mu_r, sqrt(V_r))
    
    # simulate V_des
    V_hat_r <- rgamma(n_regions, (n_r_star - 1)/2, rate = (n_r_star - 1)/(2*V_r))
    
    # calculate bias results for direct estimates
    sim_results$direct_estimates[1, s, ] <- y_hat_r - mu_r
    sim_results$direct_estimates[2, s, ] <- abs(y_hat_r - mu_r)
    sim_results$direct_estimates[3, s, ] <- abs(y_hat_r - mu_r)/abs(mu_r)
    
    # calculate CIs
    ci_80_direct <- t(rbind(y_hat_r, y_hat_r) + t(t(c(-1, 1))) %*% t(qt(0.9, n_r_star-1)*sqrt(V_hat_r)))
    ci_95_direct <- t(rbind(y_hat_r, y_hat_r) + t(t(c(-1, 1))) %*% t(qt(0.975, n_r_star-1)*sqrt(V_hat_r)))
    
    # calculate coverage and width results
    sim_results$direct_estimates[4, s, ] <- (mu_r >= ci_80_direct[,1]) & (mu_r <= ci_80_direct[,2])
    sim_results$direct_estimates[5, s, ] <- (mu_r >= ci_95_direct[,1]) & (mu_r <= ci_95_direct[,2])
    sim_results$direct_estimates[6, s, ] <- (ci_80_direct[,2] - ci_80_direct[,1])
    sim_results$direct_estimates[7, s, ] <- (ci_95_direct[,2] - ci_95_direct[,1])
    
    # prep for inla model
    data <- list(y = y_hat_r,
                 admin1 = 1:n_regions)
    formula <- as.formula("y ~ 1 + f(admin1, model = 'iid')")
    scaling <- 1/V_hat_r
    
    # fit smoothing models
    mod1 <- inla(formula,
                 family = "gaussian",
                 data = data,
                 control.family = list(hyper = list(prec = list(initial = log(1), fixed=T))),
                 control.predictor = list(compute=T),
                 control.compute = list(config=T, waic = TRUE, dic = TRUE, cpo = TRUE),
                 control.inla = list(lincomb.derived.correlation.matrix=T),
                 scale = scaling,
                 quantiles = c(0.025, 0.1, 0.5, 0.9, 0.975))
    fitted_res <- mod1$summary.fitted.values
    latent_means <- fitted_res$`0.5quant`
    
    # calculate bias results for direct estimates
    sim_results$univariate_iid[1, s, ] <- latent_means - mu_r
    sim_results$univariate_iid[2, s, ] <- abs(latent_means - mu_r)
    sim_results$univariate_iid[3, s, ] <- abs(latent_means - mu_r)/abs(mu_r)
    
    # calculate coverage and width results
    sim_results$univariate_iid[4, s, ] <- (mu_r >= fitted_res$`0.1quant`) & (mu_r <= fitted_res$`0.9quant`)
    sim_results$univariate_iid[5, s, ] <- (mu_r >= fitted_res$`0.025quant`) & (mu_r <= fitted_res$`0.975quant`)
    sim_results$univariate_iid[6, s, ] <- (fitted_res$`0.9quant` - fitted_res$`0.1quant`)
    sim_results$univariate_iid[7, s, ] <- (fitted_res$`0.975quant` - fitted_res$`0.025quant`)
}

# save results
write_rds(sim_results, paste0(savedir, "/inla-test.rds"))
