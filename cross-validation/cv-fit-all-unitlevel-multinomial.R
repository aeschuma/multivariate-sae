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
         family = "poisson",
         control.compute = list(config=T))
}


# read and format data ####
load("/Users/austin/Dropbox/dissertation_2/survey-csmf/data/ken_dhs2014/data/haz-waz-kenDHS2014.rda")

## reformat data into long form
results.long <- dat_c %>% select(admin1, admin1.name, admin1.char,
                                 cluster, region, strata, rural, weights,
                                 contraceptive) %>%
    group_by(admin1, admin1.name, admin1.char,
             cluster, region, strata, rural, contraceptive) %>%
    summarize(count = n(), weights = unique(weights)) %>%
    ungroup()

results.long$admin1.none<- ifelse(results.long$contraceptive == "none", results.long$admin1, NA)
results.long$admin1.modern <- ifelse(results.long$contraceptive == "modern", results.long$admin1, NA)
results.long$admin1.other <- ifelse(results.long$contraceptive == "other", results.long$admin1, NA)
results.long$obs <- 1:nrow(results.long)

## add outcome specific cluster indicators for nugget terms
results.long$cluster.none<- ifelse(results.long$contraceptive == "none", results.long$cluster, NA)
results.long$cluster.modern <- ifelse(results.long$contraceptive == "modern", results.long$cluster, NA)
results.long$cluster.other <- ifelse(results.long$contraceptive == "other", results.long$cluster, NA)

## add outcome specific urban/rural indicators for strata FEs
results.long$rural.none<- ifelse(results.long$contraceptive == "none", results.long$rural, NA)
results.long$rural.modern <- ifelse(results.long$contraceptive == "modern", results.long$rural, NA)
results.long$rural.other <- ifelse(results.long$contraceptive == "other", results.long$rural, NA)

## add another index for the shared component
results.long$admin1.none.2 <- results.long$admin1.none
results.long$admin1.modern.2 <- results.long$admin1.modern
results.long$admin1.other.2 <- results.long$admin1.other

# create a list of the data
data <- as.list(results.long)
data$contraceptive_factor <- factor(data$contraceptive, levels = c("none", "modern", "other"))

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

formulas[["IID nonshared"]] <- formula("count ~ -1 + contraceptive_factor * rural + 
                                        f(admin1.modern, model = 'iid', hyper = iid_prior) +
                                        f(admin1.none, model = 'iid', hyper = iid_prior) +
                                        f(cluster, model = 'iid', hyper = iid_prior)")

formulas[["BYM nonshared"]] <- formula("count ~ -1 + contraceptive_factor * rural + 
                                        f(admin1.modern, model = 'bym2',
                                             graph = admin1.mat, 
                                             scale.model = T, 
                                             constr = T,
                                             hyper = bym2_prior) +
                                        f(admin1.none, model = 'bym2',
                                              graph = admin1.mat, 
                                              scale.model = T, 
                                              constr = T,
                                              hyper = bym2_prior) +
                                        f(cluster, model = 'iid', hyper = iid_prior)")
formulas[["IID shared"]] <- formula("count ~ -1 + contraceptive_factor * rural +
                                        f(admin1.modern, model = 'iid', hyper = iid_prior) +
                                        f(admin1.none, model = 'iid', hyper = iid_prior) +
                                        f(admin1.none.2, copy = \"admin1.modern\", 
                                          fixed = FALSE, hyper = lambda_prior) +
                                        f(cluster, model = 'iid', hyper = iid_prior)")
formulas[["BYM shared"]] <- formula("count ~ -1 + contraceptive_factor * rural + 
                                        f(admin1.modern, model = 'bym2',
                                             graph = admin1.mat, 
                                             scale.model = T, 
                                             constr = T,
                                             hyper = bym2_prior) +
                                        f(admin1.none, model = 'bym2',
                                              graph = admin1.mat, 
                                              scale.model = T, 
                                              constr = T,
                                              hyper = bym2_prior) +
                                        f(admin1.none.2, copy = \"admin1.modern\", 
                                          fixed = FALSE, hyper = lambda_prior) +
                                        f(cluster, model = 'iid', hyper = iid_prior)")

# Cross validation ####
n_samples <- 500
cv_res <- tibble(model = rep(model_names, n_regions),
                 region = rep(1:n_regions, each = n_models),
                 cpo = NA)

for (r in 1:n_regions) {
    starttime <- Sys.time()
    
    ## Hold out data ####
    tmp <- data
    tmp$count[tmp$admin1 == r] <- NA
    
    # calculate true probabilities
    none_urban_true <- data$count[data$admin1 == r & data$rural == 0 & data$contraceptive == "none"]
    modern_urban_true <- data$count[data$admin1 == r & data$rural == 0 & data$contraceptive == "modern"]
    other_urban_true <- data$count[data$admin1 == r & data$rural == 0 & data$contraceptive == "other"]
    none_rural_true <- data$count[data$admin1 == r & data$rural == 1 & data$contraceptive == "none"]
    modern_rural_true <- data$count[data$admin1 == r & data$rural == 1 & data$contraceptive == "modern"]
    other_rural_true <- data$count[data$admin1 == r & data$rural == 1 & data$contraceptive == "other"]
    
    none_urban_weights_true <- data$weights[data$admin1 == r & data$rural == 0 & data$contraceptive == "none"]
    modern_urban_weights_true <- data$weights[data$admin1 == r & data$rural == 0 & data$contraceptive == "modern"]
    other_urban_weights_true <- data$weights[data$admin1 == r & data$rural == 0 & data$contraceptive == "other"]
    none_rural_weights_true <- data$weights[data$admin1 == r & data$rural == 1 & data$contraceptive == "none"]
    modern_rural_weights_true <- data$weights[data$admin1 == r & data$rural == 1 & data$contraceptive == "modern"]
    other_rural_weights_true <- data$weights[data$admin1 == r & data$rural == 1 & data$contraceptive == "other"]
    
    urban_tot_true <- none_urban_true + modern_urban_true + other_urban_true
    rural_tot_true <- none_rural_true + modern_rural_true + other_rural_true
    
    none_urban_true <- weighted.mean(none_urban_true/urban_tot_true, none_urban_weights_true)
    modern_urban_true <- weighted.mean(modern_urban_true/urban_tot_true, modern_urban_weights_true)
    other_urban_true <- weighted.mean(other_urban_true/urban_tot_true, other_urban_weights_true)
    none_rural_true <- weighted.mean(none_rural_true/rural_tot_true, none_rural_weights_true)
    modern_rural_true <- weighted.mean(modern_rural_true/rural_tot_true, modern_rural_weights_true)
    other_rural_true <- weighted.mean(other_rural_true/rural_tot_true, other_rural_weights_true)
    
    # logit transforms
    logit_m_o_urban_true <- log(modern_urban_true/other_urban_true)
    logit_n_o_urban_true <- log(none_urban_true/other_urban_true)
    logit_m_o_rural_true <- log(modern_rural_true/other_rural_true)
    logit_n_o_rural_true <- log(none_rural_true/other_rural_true)
    
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
        none_fe_idx <- which(rownames(samp[[1]]$latent) %in% c("contraceptive_factornone:1"))
        none_fe_mat <- matrix(NA, nrow = length(none_fe_idx), ncol = nsamps)
        modern_fe_idx <- which(rownames(samp[[1]]$latent) %in% c("contraceptive_factormodern:1"))
        modern_fe_mat <- matrix(NA, nrow = length(modern_fe_idx), ncol = nsamps)
        other_fe_idx <- which(rownames(samp[[1]]$latent) %in% c("contraceptive_factorother:1"))
        other_fe_mat <- matrix(NA, nrow = length(other_fe_idx), ncol = nsamps)
        
        none_rural_fe_idx <- which(rownames(samp[[1]]$latent) %in% c("rural:1"))
        none_rural_fe_mat <- matrix(NA, nrow = length(none_rural_fe_idx), ncol = nsamps)
        modern_rural_fe_idx <- which(rownames(samp[[1]]$latent) %in% c("contraceptive_factormodern:rural:1"))
        modern_rural_fe_mat <- matrix(NA, nrow = length(modern_rural_fe_idx), ncol = nsamps)
        other_rural_fe_idx <- which(rownames(samp[[1]]$latent) %in% c("contraceptive_factorother:rural:1"))
        other_rural_fe_mat <- matrix(NA, nrow = length(other_rural_fe_idx), ncol = nsamps)
        
        modern_re_idx <- rownames(samp[[1]]$latent) %>% str_detect("admin1\\.modern:") %>% which()
        none_re_idx <- rownames(samp[[1]]$latent) %>% str_detect("admin1\\.none:") %>% which()
        modern_re_mat <- matrix(NA, nrow = length(modern_re_idx), ncol = nsamps)
        none_re_mat <- matrix(NA, nrow = length(none_re_idx), ncol = nsamps)
        
        if(!grepl("nonshared", model_names[i])) {
            shared_re_idx <- rownames(samp[[1]]$latent) %>% str_detect("admin1\\.none\\.2:") %>% which()
            shared_re_mat <- matrix(NA, nrow = length(shared_re_idx), ncol = nsamps)
        } else {
            shared_re_idx <- integer(0)
            shared_re_mat <- matrix(0, nrow = length(none_re_idx), ncol = nsamps)
        }
        
        # fill in sample matrices
        for (s in 1:nsamps) {
            none_fe_mat[,s] <- samp[[s]]$latent[none_fe_idx]
            modern_fe_mat[,s] <- samp[[s]]$latent[modern_fe_idx]
            other_fe_mat[,s] <- samp[[s]]$latent[other_fe_idx]
            none_rural_fe_mat[,s] <- samp[[s]]$latent[none_rural_fe_idx]
            modern_rural_fe_mat[,s] <- samp[[s]]$latent[modern_rural_fe_idx]
            other_rural_fe_mat[,s] <- samp[[s]]$latent[other_rural_fe_idx]
            modern_re_mat[,s] <- samp[[s]]$latent[modern_re_idx]
            none_re_mat[,s] <- samp[[s]]$latent[none_re_idx]
            if(!grepl("nonshared", model_names[i])) shared_re_mat[,s] <- samp[[s]]$latent[shared_re_idx]
        }
        
        # # save sample matrices in list for output
        # model_results[[ model_names[i] ]]$samples_all <- list()
        # model_results[[ model_names[i] ]]$samples_all$haz_fe <- haz_fe_mat
        # model_results[[ model_names[i] ]]$samples_all$waz_fe <- waz_fe_mat
        # model_results[[ model_names[i] ]]$samples_all$haz_rural_fe <- haz_rural_fe_mat
        # model_results[[ model_names[i] ]]$samples_all$waz_rural_fe <- waz_rural_fe_mat
        # model_results[[ model_names[i] ]]$samples_all$haz_re <- haz_re_mat
        # model_results[[ model_names[i] ]]$samples_all$waz_re <- waz_re_mat
        # model_results[[ model_names[i] ]]$samples_all$shared_re <- shared_re_mat
        
        # obtain total effect for space - first half of bym2, or all of them for IID
        modern_re_tot_mat <- modern_re_mat[1:n_regions,]
        none_re_tot_mat <- none_re_mat[1:n_regions,] + shared_re_mat[1:n_regions,]
        
        # get urban and rural FEs
        none_urban_fe_tot_mat <- none_fe_mat
        none_rural_fe_tot_mat <- none_fe_mat + none_rural_fe_mat
        modern_urban_fe_tot_mat <- modern_fe_mat
        modern_rural_fe_tot_mat <- modern_fe_mat + modern_rural_fe_mat
        other_urban_fe_tot_mat <- other_fe_mat
        other_rural_fe_tot_mat <- other_fe_mat + other_rural_fe_mat
        
        # matrix of fitted estimates
        fitted_none_urban_mat <- none_urban_fe_tot_mat[rep(1,n_regions),] + none_re_tot_mat
        fitted_none_rural_mat <- none_rural_fe_tot_mat[rep(1,n_regions),] + none_re_tot_mat
        fitted_modern_urban_mat <- modern_urban_fe_tot_mat[rep(1,n_regions),] + modern_re_tot_mat
        fitted_modern_rural_mat <- modern_rural_fe_tot_mat[rep(1,n_regions),] + modern_re_tot_mat
        fitted_other_urban_mat <- other_urban_fe_tot_mat[rep(1,n_regions),]
        fitted_other_rural_mat <- other_rural_fe_tot_mat[rep(1,n_regions),]
        
        # differences in fitted estimates = logits
        fitted_m_o_urban <- fitted_modern_urban_mat - fitted_other_urban_mat
        fitted_n_o_urban <- fitted_none_urban_mat - fitted_other_urban_mat
        fitted_m_o_rural <- fitted_modern_rural_mat - fitted_other_rural_mat
        fitted_n_o_rural <- fitted_none_rural_mat - fitted_other_rural_mat
        
        # densities for heldout region
        gaus_m_o_urban <- list(mean = mean(fitted_m_o_urban[r,]), 
                               sd = sd(fitted_m_o_urban[r,]))
        gaus_n_o_urban <- list(mean = mean(fitted_n_o_urban[r,]), 
                               sd = sd(fitted_n_o_urban[r,]))
        gaus_m_o_rural <- list(mean = mean(fitted_m_o_rural[r,]), 
                               sd = sd(fitted_m_o_rural[r,]))
        gaus_n_o_rural <- list(mean = mean(fitted_n_o_rural[r,]), 
                               sd = sd(fitted_n_o_rural[r,]))
        
        # calculate CV results
        y_lik <- c()
        
        for (i in 1:n_samples) {
            y_lik <- c(y_lik, 
                       dnorm(x = logit_m_o_urban_true, 
                             mean = gaus_m_o_urban$mean, 
                             sd = gaus_m_o_urban$sd),
                       dnorm(x = logit_n_o_urban_true, 
                             mean = gaus_n_o_urban$mean, 
                             sd = gaus_n_o_urban$sd),
                       dnorm(x = logit_m_o_rural_true, 
                             mean = gaus_m_o_rural$mean, 
                             sd = gaus_m_o_rural$sd),
                       dnorm(x = logit_n_o_rural_true, 
                             mean = gaus_n_o_rural$mean, 
                             sd = gaus_n_o_rural$sd))
        }
        cv_res[cv_res$model == model_names[mm] & cv_res$region == r, "cpo"] <- mean(y_lik)
    }
    
    endtime <- Sys.time()
    totaltime <- endtime - starttime
    
    message(paste0("DONE! Time: ", totaltime, " min"))
}

# save results
write_rds(cv_res, file = "../../../Dropbox/dissertation_2/survey-csmf/results/cv/cv_results-multinomial-unitlevel.rds")

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
