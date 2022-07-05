# test unit level models

rm(list = ls())

library(INLA)
library(tidyverse)
library(spdep)
library(geosphere)
library(haven)
library(knitr)
library(kableExtra)
library(magrittr)
library(mvtnorm)
library(rgdal)
library(viridis)
library(classInt)
library(gridExtra)
library(ggpubr)

# load and format data ####
load("/Users/austin/Dropbox/dissertation_2/survey-csmf/data/ken_dhs2014/data/haz-waz-kenDHS2014.rda")
prop_urban <- read_rds("/Users/austin/Dropbox/dissertation_2/survey-csmf/data/ken_dhs2014/data/admin1_2014_urban_frac.rds")

# format prop urban data
admin_df <- dat_c %>% select(admin1.name, admin1) %>% distinct()
prop_urban %<>% left_join(admin_df, 
                          by = c("adm_name" = "admin1.name")) %>%
    select(admin1, urb_frac)

## reformat data into long form
results.long <- dat_c %>% select(admin1, admin1.name, admin1.char,
                                 cluster, region, strata, rural,
                                 contraceptive) %>%
    group_by(admin1, admin1.name, admin1.char,
             cluster, region, strata, rural, contraceptive) %>%
    summarize(count = n()) %>%
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
nugget_prior <- list(prec = list(prior = "gaussian",
                                 param = c(0, 0)))
iid_prior <- list(prec = list(prior = "pc.prec",
                              param = c(1, 0.01)))
bym2_prior <- list(phi=list(prior="logitbeta", param=c(1, 1), initial=0.5), 
                   prec=list(prior="pc.prec", param=c(1, 0.01), initial=5))
lambda_prior <- list(beta = list(prior = 'logtnormal', param = c(0, 1)))

# models to run  ####
model_names <- c("IID nonshared", "BYM nonshared",
                 "IID shared", "BYM shared")

# formulas ####
formulas <- vector(mode = "list", length = length(model_names))
names(formulas) <- model_names

formulas[["IID nonshared"]] <- formula("count ~ -1 + contraceptive_factor + rural.none + rural.modern + rural.other + 
                                        f(admin1.modern, model = 'iid', hyper = iid_prior) +
                                        f(admin1.none, model = 'iid', hyper = iid_prior) +
                                        f(cluster, model = 'iid', hyper = nugget_prior)")

formulas[["BYM nonshared"]] <- formula("count ~ -1 + contraceptive_factor + rural.none + rural.modern + rural.other + 
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
                                        f(cluster, model = 'iid', hyper = nugget_prior)")
formulas[["IID shared"]] <- formula("count ~ -1 + contraceptive_factor + rural.none + rural.modern + rural.other +
                                        f(admin1.modern, model = 'iid', hyper = iid_prior) +
                                        f(admin1.none, model = 'iid', hyper = iid_prior) +
                                        f(admin1.none.2, copy = \"admin1.modern\", 
                                          fixed = FALSE, hyper = lambda_prior) +
                                        f(cluster, model = 'iid', hyper = nugget_prior)")
formulas[["BYM shared"]] <- formula("count ~ -1 + contraceptive_factor + rural.none + rural.modern + rural.other + 
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
                                        f(cluster, model = 'iid', hyper = nugget_prior)")

# modeling ####
model_results <- vector(mode = "list", length = length(model_names))
names(model_results) <- model_names

model_fits <- vector(mode = "list", length = length(model_names))
names(model_fits) <- model_names

nsamps <- 1000

for (i in 1:length(model_names)) {
    message(model_names[i])
    start_time <- Sys.time()
    tmp <- inla(formulas[[ model_names[i] ]], data = data,
                family = "poisson",
                quantiles=c(0.025, 0.1, 0.5, 0.9, 0.975),
                # control.family = list(),
                control.predictor = list(compute=T),
                # control.inla = list(lincomb.derived.correlation.matrix=T),
                # control.fixed = list(prec = list(default = 0.001), correlation.matrix=T),
                control.compute = list(config=T, waic = TRUE, dic = TRUE, cpo = TRUE))
    # print(model_results[[ model_names[i] ]]$summary.fixed)
    # print(model_results[[ model_names[i] ]]$summary.hyperpar)
    
    # extract and format results
    samp <- inla.posterior.sample(n = nsamps, result = tmp)
    
    # save summaries
    model_results[[ model_names[i] ]]$summaries <- list()
    model_results[[ model_names[i] ]]$summaries$summary.fixed <- tmp$summary.fixed
    model_results[[ model_names[i] ]]$summaries$summary.random <- tmp$summary.random
    model_results[[ model_names[i] ]]$summaries$summary.hyperpar <- tmp$summary.hyperpar
    
    # process hyperpars
    hyperpar_names <- names(samp[[1]]$hyperpar)
    hyperpar_mat <- matrix(NA, nrow = length(hyperpar_names), ncol = nsamps)
    for (s in 1:nsamps) {
        hyperpar_mat[,s] <- samp[[s]]$hyperpar
    }
    rownames(hyperpar_mat) <- hyperpar_names
    
    # process fitted values
    none_fe_idx <- which(rownames(samp[[1]]$latent) %in% c("contraceptive_factornone:1"))
    none_fe_mat <- matrix(NA, nrow = length(none_fe_idx), ncol = nsamps)
    modern_fe_idx <- which(rownames(samp[[1]]$latent) %in% c("contraceptive_factormodern:1"))
    modern_fe_mat <- matrix(NA, nrow = length(modern_fe_idx), ncol = nsamps)
    other_fe_idx <- which(rownames(samp[[1]]$latent) %in% c("contraceptive_factorother:1"))
    other_fe_mat <- matrix(NA, nrow = length(other_fe_idx), ncol = nsamps)

    none_rural_fe_idx <- which(rownames(samp[[1]]$latent) %in% c("rural.none:1"))
    none_rural_fe_mat <- matrix(NA, nrow = length(none_rural_fe_idx), ncol = nsamps)
    modern_rural_fe_idx <- which(rownames(samp[[1]]$latent) %in% c("rural.modern:1"))
    modern_rural_fe_mat <- matrix(NA, nrow = length(modern_rural_fe_idx), ncol = nsamps)
    other_rural_fe_idx <- which(rownames(samp[[1]]$latent) %in% c("rural.other:1"))
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
    model_results[[ model_names[i] ]]$samples_all <- list()
    model_results[[ model_names[i] ]]$samples_all$none_fe <- none_fe_mat
    model_results[[ model_names[i] ]]$samples_all$modern_fe <- modern_fe_mat
    model_results[[ model_names[i] ]]$samples_all$other_fe <- other_fe_mat
    model_results[[ model_names[i] ]]$samples_all$none_rural_fe <- none_rural_fe_mat
    model_results[[ model_names[i] ]]$samples_all$modern_rural_fe <- modern_rural_fe_mat
    model_results[[ model_names[i] ]]$samples_all$other_rural_fe <- other_rural_fe_mat
    model_results[[ model_names[i] ]]$samples_all$none_re <- none_re_mat
    model_results[[ model_names[i] ]]$samples_all$modern_re <- modern_re_mat
    model_results[[ model_names[i] ]]$samples_all$shared_re <- shared_re_mat

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
    
    # transform estimates
    fitted_array_urban <- array(c(fitted_none_urban_mat, fitted_modern_urban_mat, fitted_other_urban_mat), 
                                dim = c(n_regions, nsamps, 3))
    fitted_array_rural <- array(c(fitted_none_rural_mat, fitted_modern_rural_mat, fitted_other_rural_mat), 
                                dim = c(n_regions, nsamps, 3))
    
    probs_urban <- array(NA, dim = c(3, n_regions, nsamps))
    probs_rural <- array(NA, dim = c(3, n_regions, nsamps))
    for (s in 1:nsamps) {
        for (r in 1:n_regions) {
            denom_urban <- exp(fitted_none_urban_mat[r, s]) + 
                exp(fitted_modern_urban_mat[r, s]) + 
                exp(fitted_other_urban_mat[r, s])
            probs_urban[1, r, s] <- exp(fitted_none_urban_mat[r, s])/denom_urban
            probs_urban[2, r, s] <- exp(fitted_modern_urban_mat[r, s])/denom_urban
            probs_urban[3, r, s] <- exp(fitted_other_urban_mat[r, s])/denom_urban
            
            denom_rural <- exp(fitted_none_rural_mat[r, s]) + 
                exp(fitted_modern_rural_mat[r, s]) + 
                exp(fitted_other_rural_mat[r, s])
            probs_rural[1, r, s] <- exp(fitted_none_rural_mat[r, s])/denom_rural
            probs_rural[2, r, s] <- exp(fitted_modern_rural_mat[r, s])/denom_rural
            probs_rural[3, r, s] <- exp(fitted_other_rural_mat[r, s])/denom_rural
        }
    }
    
    # output for region-level samples
    model_results[[ model_names[i] ]]$samples_region_ur <- list()
    model_results[[ model_names[i] ]]$samples_region_ur$probs_urban <- probs_urban
    model_results[[ model_names[i] ]]$samples_region_ur$probs_rural <- probs_rural
    
    # aggregate over urban/rural
    fitted_none_mat <- matrix(NA, nrow = n_regions, ncol = nsamps)
    fitted_modern_mat <- matrix(NA, nrow = n_regions, ncol = nsamps)
    fitted_other_mat <- matrix(NA, nrow = n_regions, ncol = nsamps)
    
    for (s in 1:nsamps)  {
        fitted_none_mat[, s] <- (probs_urban[1, , s] * prop_urban$urb_frac) + 
            (probs_rural[1, , s] * (1 - prop_urban$urb_frac))
        fitted_modern_mat[, s] <- (probs_urban[2, , s] * prop_urban$urb_frac) + 
            (probs_rural[2, , s] * (1 - prop_urban$urb_frac))
        fitted_other_mat[, s] <- (probs_urban[3, , s] * prop_urban$urb_frac) + 
            (probs_rural[3, , s] * (1 - prop_urban$urb_frac))
    }
    
    # store results
    model_results[[ model_names[i] ]]$samples_final <- list()
    model_results[[ model_names[i] ]]$samples_final$none <- fitted_none_mat
    model_results[[ model_names[i] ]]$samples_final$modern <- fitted_modern_mat
    model_results[[ model_names[i] ]]$samples_final$other <- fitted_other_mat
    
    # run time
    model_results[[ model_names[i] ]]$run_time <- Sys.time() - start_time
    
    # store full results
    model_fits[[ model_names[i] ]] <- tmp
    rm(tmp)
}

# model results ####
write_rds(model_results, file = "../../../Dropbox/dissertation_2/survey-csmf/results/ken2014-unit-level/multinomial/ken2014-unit-level-multinomial-flatnugget-inla-posteriors.rds")
write_rds(model_fits, file = "../../../Dropbox/dissertation_2/survey-csmf/results/ken2014-unit-level/multinomial/ken2014-unit-level-multinomial-flatnugget-inla-fits.rds")

# model selection
model_selection <- tibble(model = character(),
                          waic = numeric(),
                          dic = numeric(),
                          cpo = numeric())
for(i in 1:length(model_names)) {
    add_row <- tibble(model = model_names[i],
                      waic = model_fits[[ model_names[i] ]]$waic$waic,
                      dic = model_fits[[ model_names[i] ]]$dic$dic,
                      cpo = -sum(log(model_fits[[ model_names[i] ]]$cpo$cpo)))
    model_selection %<>% bind_rows(add_row)
}
model_selection %>% arrange(desc(waic))
