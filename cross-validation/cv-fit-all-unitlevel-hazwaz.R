# Austin Schumacher
# 1/21/2022
# Fit all INLA models and save results

# preamble ####
rm(list = ls())

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
load("/Users/austin/Dropbox/dissertation_2/survey-csmf/data/ken_dhs2014/data/haz-waz-kenDHS2014.rda")

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
n_samples <- 500
cv_res <- tibble(model = rep(model_names, n_regions),
                 region = rep(1:n_regions, each = n_models),
                 cpo = NA)

for (r in 1:n_regions) {
    starttime <- Sys.time()
    
    ## Hold out data ####
    tmp <- data
    tmp$value[tmp$admin1 == r] <- NA
    tmp$value.haz <- ifelse(tmp$outcome == "HAZ", tmp$value, NA)
    tmp$value.waz <- ifelse(tmp$outcome == "WAZ", tmp$value, NA)
    
    tmp$Y <- cbind(tmp$value.haz, tmp$value.waz)
    
    haz_urban_true <- data$value.haz[data$admin1 == r & data$rural == 0 & !is.na(data$value.haz)]
    waz_urban_true <- data$value.waz[data$admin1 == r & data$rural == 0 & !is.na(data$value.waz)]
    haz_rural_true <- data$value.haz[data$admin1 == r & data$rural == 1 & !is.na(data$value.haz)] 
    waz_rural_true <- data$value.waz[data$admin1 == r & data$rural == 1 & !is.na(data$value.waz)] 
    
    haz_urban_weights <- data$weights[data$admin1 == r & data$rural == 0 & !is.na(data$value.haz)]
    waz_urban_weights <- data$weights[data$admin1 == r & data$rural == 0 & !is.na(data$value.waz)]
    haz_rural_weights <- data$weights[data$admin1 == r & data$rural == 1 & !is.na(data$value.haz)]
    waz_rural_weights <- data$weights[data$admin1 == r & data$rural == 1 & !is.na(data$value.waz)]
    
    haz_urban_true <- weighted.mean(haz_urban_true, haz_urban_weights)
    waz_urban_true <- weighted.mean(waz_urban_true, waz_urban_weights)
    haz_rural_true <- weighted.mean(haz_rural_true, haz_rural_weights)
    waz_rural_true <- weighted.mean(waz_rural_true, waz_rural_weights)
    
    ## Fit models ####
    mod_list <- vector(mode = "list", length = n_models)
    names(mod_list) <- model_names
    
    mod_list$`IID nonshared` <- fitINLA(formula = formulas$`IID nonshared`, 
                                        data = tmp)
    mod_list$`BYM nonshared` <- fitINLA(formula = formulas$`BYM nonshared`, 
                                        data = tmp)
    mod_list$`IID shared` <- fitINLA(formula = formulas$`IID shared`, 
                                     data = tmp)
    mod_list$`BYM shared` <- fitINLA(formula = formulas$`BYM shared`, 
                                     data = tmp)
    
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
        
        # densities for heldout region
        gaus_haz_urban <- list(mean = mean(fitted_haz_urban_mat[r,]), 
                            sd = sd(fitted_haz_urban_mat[r,]))
        gaus_waz_urban <- list(mean = mean(fitted_waz_urban_mat[r,]), 
                            sd = sd(fitted_waz_urban_mat[r,]))
        gaus_haz_rural <- list(mean = mean(fitted_haz_rural_mat[r,]), 
                            sd = sd(fitted_haz_rural_mat[r,]))
        gaus_waz_rural <- list(mean = mean(fitted_waz_rural_mat[r,]), 
                            sd = sd(fitted_waz_rural_mat[r,]))
        
        # calculate CV results
        y_lik <- c()
        
        for (i in 1:n_samples) {
            y_lik <- c(y_lik, 
                       dnorm(x = haz_urban_true, 
                             mean = gaus_haz_urban$mean, 
                             sd = gaus_haz_urban$sd),
                       dnorm(x = waz_urban_true, 
                             mean = gaus_waz_urban$mean, 
                             sd = gaus_waz_urban$sd),
                       dnorm(x = haz_rural_true, 
                             mean = gaus_haz_rural$mean, 
                             sd = gaus_haz_rural$sd),
                       dnorm(x = waz_rural_true, 
                             mean = gaus_waz_rural$mean, 
                             sd = gaus_waz_rural$sd))
        }
        cv_res[cv_res$model == model_names[mm] & cv_res$region == r, "cpo"] <- mean(y_lik)
    }
    
    endtime <- Sys.time()
    totaltime <- endtime - starttime
}

# save results
write_rds(cv_res, file = "../../../Dropbox/dissertation_2/survey-csmf/results/cv/cv_results-hazwaz-unitlevel.rds")

# format
cv_res %<>% mutate(model_factor = factor(model, levels = model_names))

logCPOres <- cv_res %>% mutate(logcpo = log(cpo)) %>% 
    group_by(model_factor) %>% 
    summarise(logCPO_num = round(-1 * sum(logcpo), 2))

logCPOres$logCPO <- ifelse(logCPOres$logCPO_num == min(logCPOres$logCPO_num), 
                           paste0("\\textbf{", logCPOres$logCPO_num, "}"),
                           as.character(logCPOres$logCPO_num))

logCPOres %>% select(model_factor, logCPO) %>%
    kable(format = "markdown", caption = "Bivariate $-\\sum\\log(CPO)$ for each model. Bold indicates the best performing model",
          col.names = c("Model", "$-\\sum\\log(CPO)$"))
