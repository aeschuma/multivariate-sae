# Austin Schumacher
# 1/21/2022
# Fit all INLA models and save results

# preamble ####
rm(list = ls())

## set the root depending on operating system
root <- ifelse(Sys.info()[1]=="Darwin","~/",
               ifelse(Sys.info()[1]=="Windows","P:/",
                      ifelse(Sys.info()[1]=="Linux","/home/users/aeschuma/",
                             stop("Unknown operating system"))))

if (root ==  "/home/users/aeschuma/") {
    setwd(paste0(root, "Desktop/survey_csmf/cross-validation"))
}

library(tidyverse)
library(spdep)
library(geosphere)
library(haven)
library(knitr)
library(kableExtra)
library(magrittr)
library(svyVGAM)
library(mvtnorm)
library(rgdal)
library(bayesplot)
library(INLA)
library(viridis)
library(classInt)
library(gridExtra)
library(ggpubr)

# functions ####
fitINLA <- function(formula, data) {
    inla(formula, data = data, 
         family = rep("gaussian", 2),
         control.compute = list(config=T))
}


# read and format data ####
load(paste0(root,"Dropbox/dissertation_2/survey-csmf/data/ken_dhs2014/data/haz-waz-kenDHS2014.rda"))

n_regions <- nrow(poly.adm1)

## reformat data into long form
results.long <- dat %>% select(admin1, admin1.name, admin1.char,
                               cluster, region, strata, rural, weights,
                               HAZ, WAZ) %>%
    pivot_longer(cols = c(HAZ, WAZ),
                 names_to = "outcome",
                 values_to = "value")

results.long$value.haz <- ifelse(results.long$outcome == "HAZ", results.long$value, NA)
results.long$value.waz <- ifelse(results.long$outcome == "WAZ", results.long$value, NA)

results.long$admin1.haz <- ifelse(results.long$outcome == "HAZ", results.long$admin1, NA)
results.long$admin1.waz <- ifelse(results.long$outcome == "WAZ", results.long$admin1, NA)
results.long$obs <- 1:nrow(results.long)

## add outcome specific cluster indicators for nugget terms
results.long$cluster.haz <- ifelse(results.long$outcome == "HAZ", results.long$cluster, NA)
results.long$cluster.waz <- ifelse(results.long$outcome == "WAZ", results.long$cluster, NA)

## add outcome specific urban/rural indicators for strata FEs
results.long$rural.haz <- ifelse(results.long$outcome == "HAZ", results.long$rural, NA)
results.long$rural.waz <- ifelse(results.long$outcome == "WAZ", results.long$rural, NA)

## add another index for the shared component
results.long$admin1.haz.2 <- results.long$admin1.haz

# create a list of the data
data <- as.list(results.long)

# create matrix for outcome in order to have two likelihoods
data$Y <- cbind(data$value.haz, data$value.waz)

# get number of regions
n_regions <- length(unique(data$admin1))

# priors ####
iid_prior <- list(prec = list(prior = "pc.prec",
                              param = c(1, 0.01)))
bym2_prior <- list(phi=list(prior="logitbeta", param=c(1, 1), initial=0.5), 
                   prec=list(prior="pc.prec", param=c(1, 0.01), initial=5))
lambda_prior <- list(beta = list(prior = 'logtnormal', param = c(0, 1)))

# models to run  ####
model_names <- c("IID nonshared", "BYM nonshared",
                 "IID shared", "BYM shared")
n_models <- length(model_names)

# formulas ####
formulas <- vector(mode = "list", length = length(model_names))
names(formulas) <- model_names

formulas[["IID nonshared"]] <- formula("Y ~ -1 + outcome + rural.haz + rural.waz + 
                                        f(admin1.haz, model = 'iid', hyper = iid_prior) +
                                        f(admin1.waz, model = 'iid', hyper = iid_prior) +
                                        f(cluster.haz, model = 'iid', hyper = iid_prior) +
                                        f(cluster.waz, model = 'iid', hyper = iid_prior)")

formulas[["BYM nonshared"]] <- formula("Y ~ -1 + outcome + rural.haz + rural.waz + 
                                        f(admin1.haz, model = 'bym2',
                                             graph = admin1.mat, 
                                             scale.model = T, 
                                             constr = T,
                                             hyper = bym2_prior) +
                                        f(admin1.waz, model = 'bym2',
                                              graph = admin1.mat, 
                                              scale.model = T, 
                                              constr = T,
                                              hyper = bym2_prior) +
                                        f(cluster.haz, model = 'iid', hyper = iid_prior) +
                                        f(cluster.waz, model = 'iid', hyper = iid_prior)")
formulas[["IID shared"]] <- formula("Y ~ -1 + outcome + rural.haz + rural.waz +
                                        f(admin1.haz, model = 'iid', hyper = iid_prior) +
                                        f(admin1.waz, model = 'iid', hyper = iid_prior) +
                                        f(admin1.haz.2, copy = \"admin1.waz\", 
                                          fixed = FALSE, hyper = lambda_prior) +
                                        f(cluster.haz, model = 'iid', hyper = iid_prior) +
                                        f(cluster.waz, model = 'iid', hyper = iid_prior)")
formulas[["BYM shared"]] <- formula("Y ~ -1 + outcome + rural.haz + rural.waz + 
                                        f(admin1.haz, model = 'bym2',
                                             graph = admin1.mat, 
                                             scale.model = T, 
                                             constr = T,
                                             hyper = bym2_prior) +
                                        f(admin1.waz, model = 'bym2',
                                              graph = admin1.mat, 
                                              scale.model = T, 
                                              constr = T,
                                              hyper = bym2_prior) +
                                        f(admin1.haz.2, copy = \"admin1.waz\", 
                                          fixed = FALSE, hyper = lambda_prior) +
                                        f(cluster.haz, model = 'iid', hyper = iid_prior) +
                                        f(cluster.waz, model = 'iid', hyper = iid_prior)")

# Cross validation ####
n_samples <- 250

cv_res_obs_likelihood_areamean_cpo <- tibble(model = rep(model_names, n_regions),
                                            region = rep(1:n_regions, each = n_models),
                                            cpo = NA)
cv_res_directest_bivariate_posterior_areamean_cpo <- tibble(model = rep(model_names, n_regions),
                                               region = rep(1:n_regions, each = n_models),
                                               cpo = NA)
cv_res_obs_areamean_mse_cpo <- tibble(model = rep(model_names, n_regions),
                                      region = rep(1:n_regions, each = n_models),
                                      cpo = NA)
cv_res_directest_areamean_mse_cpo <- tibble(model = rep(model_names, n_regions),
                                            region = rep(1:n_regions, each = n_models),
                                            cpo = NA)
cv_res_directest_areamean_bivariate_posterior <- tibble(model = rep(model_names, n_regions),
                                                        region = rep(1:n_regions, each = n_models),
                                                        cpo = NA)
cv_res_obs_areamean_mse <- tibble(model = rep(model_names, n_regions),
                                        region = rep(1:n_regions, each = n_models),
                                        cpo = NA)
cv_res_directest_areamean_mse <- tibble(model = rep(model_names, n_regions),
                                          region = rep(1:n_regions, each = n_models),
                                          cpo = NA)

# loop through regions to hold out, fit models, and do CV
# for (r in 1:n_regions) {
for (r in c(28, 30)) {
        
    message(paste0("region ", r))
    
    starttime <- Sys.time()
    
    ## Hold out data ####
    tmp <- data
    tmp$value[tmp$admin1 == r] <- NA
    tmp$value.haz <- ifelse(tmp$outcome == "HAZ", tmp$value, NA)
    tmp$value.waz <- ifelse(tmp$outcome == "WAZ", tmp$value, NA)
    
    tmp$Y <- cbind(tmp$value.haz, tmp$value.waz)
    
    haz_obs_urban_true <- data$value.haz[data$admin1 == r & data$rural == 0 & !is.na(data$value.haz)]
    waz_obs_urban_true <- data$value.waz[data$admin1 == r & data$rural == 0 & !is.na(data$value.waz)]
    haz_obs_rural_true <- data$value.haz[data$admin1 == r & data$rural == 1 & !is.na(data$value.haz)] 
    waz_obs_rural_true <- data$value.waz[data$admin1 == r & data$rural == 1 & !is.na(data$value.waz)] 
    
    haz_urban_weights <- data$weights[data$admin1 == r & data$rural == 0 & !is.na(data$value.haz)]
    waz_urban_weights <- data$weights[data$admin1 == r & data$rural == 0 & !is.na(data$value.waz)]
    haz_rural_weights <- data$weights[data$admin1 == r & data$rural == 1 & !is.na(data$value.haz)]
    waz_rural_weights <- data$weights[data$admin1 == r & data$rural == 1 & !is.na(data$value.waz)]
    
    haz_weighted_urban_true <- weighted.mean(haz_obs_urban_true, haz_urban_weights)
    waz_weighted_urban_true <- weighted.mean(waz_obs_urban_true, waz_urban_weights)
    haz_weighted_rural_true <- weighted.mean(haz_obs_rural_true, haz_rural_weights)
    waz_weighted_rural_true <- weighted.mean(waz_obs_rural_true, waz_rural_weights)
    
    ## Fit models ####
    mod_list <- vector(mode = "list", length = n_models)
    names(mod_list) <- model_names
    
    message("Running models...")
    
    mod_list$`IID nonshared` <- fitINLA(formula = formulas$`IID nonshared`, 
                                        data = tmp)
    mod_list$`BYM nonshared` <- fitINLA(formula = formulas$`BYM nonshared`, 
                                        data = tmp)
    mod_list$`IID shared` <- fitINLA(formula = formulas$`IID shared`, 
                                     data = tmp)
    mod_list$`BYM shared` <- fitINLA(formula = formulas$`BYM shared`, 
                                     data = tmp)
    
    message("Posterior sims + calculations...")
    
    # simulating posterior dist and do CV calculations for each model
    for (mm in 1:length(mod_list)) {
        # sample from the posterior
        samp <- inla.posterior.sample(n = n_samples, mod_list[[mm]])
        
        # process fitted values
        haz_fe_idx <- rownames(samp[[1]]$latent) %>% str_detect("outcomeHAZ:1") %>% which()
        waz_fe_idx <- rownames(samp[[1]]$latent) %>% str_detect("outcomeWAZ:1") %>% which()
        haz_fe_mat <- matrix(NA, nrow = length(haz_fe_idx), ncol = n_samples)
        waz_fe_mat <- matrix(NA, nrow = length(waz_fe_idx), ncol = n_samples)
        
        haz_rural_fe_idx <- rownames(samp[[1]]$latent) %>% str_detect("rural\\.haz") %>% which()
        waz_rural_fe_idx <- rownames(samp[[1]]$latent) %>% str_detect("rural\\.waz") %>% which()
        haz_rural_fe_mat <- matrix(NA, nrow = length(haz_rural_fe_idx), ncol = n_samples)
        waz_rural_fe_mat <- matrix(NA, nrow = length(waz_rural_fe_idx), ncol = n_samples)
        
        haz_re_idx <- rownames(samp[[1]]$latent) %>% str_detect("admin1\\.haz:") %>% which()
        waz_re_idx <- rownames(samp[[1]]$latent) %>% str_detect("admin1\\.waz:") %>% which()
        haz_re_mat <- matrix(NA, nrow = length(haz_re_idx), ncol = n_samples)
        waz_re_mat <- matrix(NA, nrow = length(waz_re_idx), ncol = n_samples)
        
        if(!grepl("nonshared", model_names[mm])) {
            shared_re_idx <- rownames(samp[[1]]$latent) %>% str_detect("admin1\\.haz\\.2:") %>% which()
            shared_re_mat <- matrix(NA, nrow = length(shared_re_idx), ncol = n_samples)
        } else {
            shared_re_idx <- integer(0)
            shared_re_mat <- matrix(0, nrow = nrow(haz_re_mat), ncol = n_samples)
        }
        
        # fill in sample matrices
        for (s in 1:n_samples) {
            haz_fe_mat[,s] <- samp[[s]]$latent[haz_fe_idx]
            waz_fe_mat[,s] <- samp[[s]]$latent[waz_fe_idx]
            haz_rural_fe_mat[,s] <- samp[[s]]$latent[haz_rural_fe_idx]
            waz_rural_fe_mat[,s] <- samp[[s]]$latent[waz_rural_fe_idx]
            haz_re_mat[,s] <- samp[[s]]$latent[haz_re_idx]
            waz_re_mat[,s] <- samp[[s]]$latent[waz_re_idx]
            if(!grepl("nonshared", model_names[mm])) shared_re_mat[,s] <- samp[[s]]$latent[shared_re_idx]
        }
        
        # obtain total effect for space - first half of bym2, or all of them for IID
        haz_re_tot_mat <- haz_re_mat[1:n_regions,] + shared_re_mat[1:n_regions,]
        waz_re_tot_mat <- waz_re_mat[1:n_regions,]
        
        # get urban and rural FEs
        haz_urban_fe_tot_mat <- haz_fe_mat
        haz_rural_fe_tot_mat <- haz_fe_mat + haz_rural_fe_mat
        waz_urban_fe_tot_mat <- waz_fe_mat
        waz_rural_fe_tot_mat <- waz_fe_mat + waz_rural_fe_mat
        
        # matrix of fitted estimates
        fitted_haz_urban_mat <- haz_urban_fe_tot_mat[rep(1,n_regions),] + haz_re_tot_mat
        fitted_waz_urban_mat <- waz_urban_fe_tot_mat[rep(1,n_regions),] + waz_re_tot_mat
        fitted_haz_rural_mat <- haz_rural_fe_tot_mat[rep(1,n_regions),] + haz_re_tot_mat
        fitted_waz_rural_mat <- waz_rural_fe_tot_mat[rep(1,n_regions),] + waz_re_tot_mat
        
        # gaussian precisions for observations
        sd_haz <- unlist(lapply(samp,
                                function(s) 1/sqrt(s$hyperpar["Precision for the Gaussian observations"])))
        sd_waz <- unlist(lapply(samp,
                                function(s) 1/sqrt(s$hyperpar["Precision for the Gaussian observations[2]"])))
        
        # bivariate posterior dists of area-strata means
        post_mean_urban <- c(mean(fitted_haz_urban_mat[r,]), mean(fitted_waz_urban_mat[r,]))
        post_sd_urban <- c(sd(fitted_haz_urban_mat[r,]), sd(fitted_waz_urban_mat[r,]))
        post_corr_urban <- cor(fitted_haz_urban_mat[r,], fitted_waz_urban_mat[r,])
        post_sigma_urban <- diag(post_sd_urban) %*% matrix(c(1, post_corr_urban, post_corr_urban, 1), nrow = 2) %*% diag(post_sd_urban) 
        
        post_mean_rural <- c(mean(fitted_haz_rural_mat[r,]), mean(fitted_waz_rural_mat[r,]))
        post_sd_rural <- c(sd(fitted_haz_rural_mat[r,]), sd(fitted_waz_rural_mat[r,]))
        post_corr_rural <- cor(fitted_haz_rural_mat[r,], fitted_waz_rural_mat[r,])
        post_sigma_rural <- diag(post_sd_rural) %*% matrix(c(1, post_corr_rural, post_corr_rural, 1), nrow = 2) %*% diag(post_sd_rural) 
        
        # calculate CV results
        
        # makes sense, but it's not really a posterior predictive measure 
        y_lik_directest_areamean_bivariate_posterior <- c(dmvnorm(c(haz_weighted_urban_true, waz_weighted_urban_true),
                                                                  mean = post_mean_urban, sigma = post_sigma_urban))
        
        # doesn't really make sense
        y_lik_obs_areamean_mse <- c((haz_obs_urban_true - post_mean_urban[1])^2,
                                    (waz_obs_urban_true - post_mean_urban[2])^2)
        
        # makes sense, but no posterior distribution is used
        y_lik_directest_areamean_mse <- c((haz_weighted_urban_true - post_mean_urban[1])^2,
                                          (waz_weighted_urban_true - post_mean_urban[2])^2)
        if (!(r %in% c(28, 30))) {
            y_lik_directest_areamean_bivariate_posterior <- c(y_lik_directest_areamean_bivariate_posterior,
                                                              dmvnorm(c(haz_weighted_rural_true, waz_weighted_rural_true),
                                                                      mean = post_mean_rural, sigma = post_sigma_rural))
            y_lik_obs_areamean_mse <- c(y_lik_obs_areamean_mse,
                                        (haz_obs_rural_true - post_mean_rural[1])^2,
                                        (waz_obs_rural_true - post_mean_rural[2])^2)
            y_lik_directest_areamean_mse <- c(y_lik_directest_areamean_mse,
                                              (haz_weighted_rural_true - post_mean_rural[1])^2,
                                              (waz_weighted_rural_true - post_mean_rural[2])^2)
        }         
        
        # kinda makes sense... if comparing to held out direct estimates is not good
        y_lik_obs_likelihood_areamean_cpo <- c()
        
        # makes sense if direct estimate is ok to use for comparison
        y_lik_directest_bivariate_posterior_areamean_cpo <- c()
        
        # if we want to use posterior distribution along with MSE, then these are good
        y_lik_obs_areamean_mse_cpo <- c()
        y_lik_directest_areamean_mse_cpo <- c()
        
        for (s in 1:n_samples) {
            y_lik_obs_likelihood_areamean_cpo <- c(y_lik_obs_likelihood_areamean_cpo,
                                                   dnorm(haz_obs_urban_true, 
                                                         mean = fitted_haz_urban_mat[r, s],
                                                         sd = sd_haz[s]),
                                                   dnorm(waz_obs_urban_true, 
                                                         mean = fitted_waz_urban_mat[r, s],
                                                         sd = sd_waz[s]))
            y_lik_directest_bivariate_posterior_areamean_cpo <- c(y_lik_directest_bivariate_posterior_areamean_cpo,
                                                                  dmvnorm(c(haz_weighted_urban_true, waz_weighted_urban_true),
                                                                          mean = c(fitted_haz_urban_mat[r, s],
                                                                                   fitted_waz_urban_mat[r, s]), 
                                                                          sigma = post_sigma_urban))
            y_lik_obs_areamean_mse_cpo <- c(y_lik_obs_areamean_mse_cpo,
                                            (haz_obs_urban_true - fitted_haz_urban_mat[r, s])^2,
                                            (waz_obs_urban_true - fitted_waz_urban_mat[r, s])^2)
            y_lik_directest_areamean_mse_cpo <- c(y_lik_directest_areamean_mse_cpo,
                                                  (haz_weighted_urban_true - fitted_haz_urban_mat[r, s])^2,
                                                  (waz_weighted_urban_true - fitted_waz_urban_mat[r, s])^2)
            
            #  completely urban regions can only get urban ests
            if (!(r %in% c(28, 30))) {
                y_lik_obs_likelihood_areamean_cpo <- c(y_lik_obs_likelihood_areamean_cpo,
                                                       dnorm(haz_obs_rural_true, 
                                                             mean = fitted_haz_rural_mat[r, s],
                                                             sd = sd_haz[s]),
                                                       dnorm(waz_obs_rural_true, 
                                                             mean = fitted_waz_rural_mat[r, s],
                                                             sd = sd_waz[s]))
                y_lik_directest_bivariate_posterior_areamean_cpo <- c(y_lik_directest_bivariate_posterior_areamean_cpo,
                                                                      dmvnorm(c(haz_weighted_rural_true, waz_weighted_rural_true),
                                                                              mean = c(fitted_haz_rural_mat[r, s],
                                                                                       fitted_waz_rural_mat[r, s]), 
                                                                              sigma = post_sigma_rural))
                y_lik_obs_areamean_mse_cpo <- c(y_lik_obs_areamean_mse_cpo,
                                                (haz_obs_rural_true - fitted_haz_rural_mat[r, s])^2,
                                                (waz_obs_rural_true - fitted_waz_rural_mat[r, s])^2)
                y_lik_directest_areamean_mse_cpo <- c(y_lik_directest_areamean_mse_cpo,
                                                      (haz_weighted_rural_true - fitted_haz_rural_mat[r, s])^2,
                                                      (waz_weighted_rural_true - fitted_waz_rural_mat[r, s])^2)
            }
        }
        
        cv_res_obs_likelihood_areamean_cpo[cv_res_obs_likelihood_areamean_cpo$model == model_names[mm] & cv_res_obs_likelihood_areamean_cpo$region == r, "cpo"] <- 
            mean(y_lik_obs_likelihood_areamean_cpo)
        cv_res_directest_bivariate_posterior_areamean_cpo[cv_res_directest_bivariate_posterior_areamean_cpo$model == model_names[mm] & cv_res_directest_bivariate_posterior_areamean_cpo$region == r, "cpo"] <- 
            mean(y_lik_directest_bivariate_posterior_areamean_cpo)
        cv_res_obs_areamean_mse_cpo[cv_res_obs_areamean_mse_cpo$model == model_names[mm] & cv_res_obs_areamean_mse_cpo$region == r, "cpo"] <- 
            mean(y_lik_obs_areamean_mse_cpo)
        cv_res_directest_areamean_mse_cpo[cv_res_directest_areamean_mse_cpo$model == model_names[mm] & cv_res_directest_areamean_mse_cpo$region == r, "cpo"] <- 
            mean(y_lik_directest_areamean_mse_cpo)
        cv_res_directest_areamean_bivariate_posterior[cv_res_directest_areamean_bivariate_posterior$model == model_names[mm] & cv_res_directest_areamean_bivariate_posterior$region == r, "cpo"] <- 
            mean(y_lik_directest_areamean_bivariate_posterior)
        cv_res_obs_areamean_mse[cv_res_obs_areamean_mse$model == model_names[mm] & cv_res_obs_areamean_mse$region == r, "cpo"] <- 
            mean(y_lik_obs_areamean_mse)
        cv_res_directest_areamean_mse[cv_res_directest_areamean_mse$model == model_names[mm] & cv_res_directest_areamean_mse$region == r, "cpo"] <- 
            mean(y_lik_directest_areamean_mse)
    }
    
    endtime <- Sys.time()
    totaltime <- endtime - starttime
    message(paste0("DONE! Time: ", totaltime, " min"))
}

# save the ones we want
cv_res <- list(obs_likelihood_areamean_cpo = cv_res_obs_likelihood_areamean_cpo,
               directest_bivariate_posterior_areamean_cpo = cv_res_directest_bivariate_posterior_areamean_cpo, 
               obs_areamean_mse_cpo = cv_res_obs_areamean_mse_cpo,
               directest_areamean_mse_cpo = cv_res_directest_areamean_mse_cpo,
               directest_areamean_bivariate_posterior = cv_res_directest_areamean_bivariate_posterior,
               obs_areamean_mse = cv_res_obs_areamean_mse,
               directest_areamean_mse = cv_res_directest_areamean_mse)

# save results
write_rds(cv_res, file = "../../../Dropbox/dissertation_2/survey-csmf/results/cv/cv_results-hazwaz-unitlevel.rds")

# if on Box, copy to real dropbox
if (root == "P:/") {
    write_rds(cv_res, file = "C:/Users/aeschuma/Dropbox/dissertation_2/survey-csmf/results/cv/cv_results-hazwaz-unitlevel.rds")
}

# format
for (i in 1:length(cv_res)) {
    cv_res[[i]] %<>% mutate(model_factor = factor(model, levels = model_names))
}

logCPOres <- lapply(cv_res, 
                    function(x) {
                        x %<>% mutate(logcpo = log(cpo)) %>% 
                        group_by(model_factor) %>% 
                        summarise(logCPO_num = round(-1 * sum(logcpo), 2)) %>%
                        arrange(desc(logCPO_num))
                        
                        x$logCPO <- ifelse(x$logCPO_num == min(x$logCPO_num), 
                                           paste0("\\textbf{", x$logCPO_num, "}"),
                                           as.character(x$logCPO_num))
                        
                        x
                        })

for (i in 1:length(logCPOres)) {
    logCPOres[[i]] %>% select(model_factor, logCPO) %>%
        kable(format = "markdown", caption = names(logCPOres)[i],
              col.names = c("Model", "$-\\sum\\log(CPO)$")) %>%
        print()
}
