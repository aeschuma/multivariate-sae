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

score_res <- expand_grid(model = model_names, region = 1:n_regions, rural = 0:1) %>% 
    mutate(score = NA)

mse_res <- expand_grid(model = model_names, region = 1:n_regions, rural = 0:1) %>% 
    mutate(mse = NA)

# in order to get direct estimates
my.svydesign <- survey::svydesign(ids = ~ cluster,
                                  strata = ~ strata, nest = T, weights = ~weights,
                                  data = dat)

# loop through regions to hold out, fit models, and do CV
for (r in 1:n_regions) {
# for (r in c(28, 30)) {
        
    message(paste0("region ", r))
    
    starttime <- Sys.time()
    
    ## Hold out data ####
    tmp <- data
    tmp$value[tmp$admin1 == r] <- NA
    tmp$value.haz <- ifelse(tmp$outcome == "HAZ", tmp$value, NA)
    tmp$value.waz <- ifelse(tmp$outcome == "WAZ", tmp$value, NA)
    
    tmp$Y <- cbind(tmp$value.haz, tmp$value.waz)
    
    # direct estimates for held out data
    heldout_urban <- subset(my.svydesign, admin1 == r & rural == 0)
    heldout_rural <- subset(my.svydesign, admin1 == r & rural == 1)
    
    means.svymean_urban <- svymean(~ HAZ + WAZ, heldout_urban)
    means.svymean_rural <- svymean(~ HAZ + WAZ, heldout_rural)
    
    vcov.svymean_urban <- vcov(means.svymean_urban)
    vcov.svymean_rural <- vcov(means.svymean_rural)
    
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
        
        # calculate score and mse results
        y_lik_urban <- c()
        y_lik_rural <- c()
        sq_error_urban <- c()
        sq_error_rural <- c()
        
        for (s in 1:n_samples) {
            
            post_mean_urban <- c(fitted_haz_urban_mat[r, s],
                                 fitted_waz_urban_mat[r, s])
            y_lik_urban <- c(y_lik_urban,
                             dmvnorm(as.numeric(means.svymean_urban),
                                     mean = post_mean_urban,
                                     sigma = vcov.svymean_urban))
            sq_error_urban <- c(sq_error_urban,
                                (means.svymean_urban - post_mean_urban)^2)
            
            #  completely urban regions can only get urban ests
            if (!(r %in% c(28, 30))) {
                post_mean_rural <- c(fitted_haz_rural_mat[r, s],
                                     fitted_waz_rural_mat[r, s])
                y_lik_rural <- c(y_lik_rural,
                                 dmvnorm(as.numeric(means.svymean_rural),
                                         mean = post_mean_rural,
                                         sigma = vcov.svymean_rural))
                sq_error_rural <- c(sq_error_rural,
                                    (means.svymean_rural - post_mean_rural)^2)
            }
        }
        
        score_res[score_res$model == model_names[mm] & score_res$region == r & score_res$rural == 0, "score"] <- -log(mean(y_lik_urban))
        score_res[score_res$model == model_names[mm] & score_res$region == r & score_res$rural == 1, "score"] <- -log(mean(y_lik_rural))
        mse_res[mse_res$model == model_names[mm] & mse_res$region == r & mse_res$rural == 0, "mse"] <- mean(sq_error_urban)
        mse_res[mse_res$model == model_names[mm] & mse_res$region == r & mse_res$rural == 1, "mse"] <- mean(sq_error_rural)
    }
    
    endtime <- Sys.time()
    totaltime <- endtime - starttime
    message(paste0("DONE! Time: ", totaltime, " min"))
}

# save the ones we want
cv_res_list <- list(score = score_res,
                    mse = mse_res)

# save results
write_rds(cv_res_list, file = "../../../Dropbox/dissertation_2/survey-csmf/results/cv/cv_results-hazwaz-unitlevel.rds")

# if on Box, copy to real dropbox
if (root == "P:/") {
    write_rds(cv_res_list, file = "C:/Users/aeschuma/Dropbox/dissertation_2/survey-csmf/results/cv/cv_results-hazwaz-unitlevel.rds")
}

# format
cv_res <- vector(mode = "list", length = length(cv_res_list))
names(cv_res) <- names(cv_res_list)
for (i in 1:length(cv_res)) {
    cv_res[[i]] <- cv_res_list[[i]] %>% mutate(model_factor = factor(model, levels = model_names))
}

score_to_plot <- cv_res$score %>%
    filter(!is.na(score)) %>%
    group_by(model_factor) %>%
    summarise(score = round(mean(score), 3)) %>%
    arrange(desc(score))
score_to_plot$score <- ifelse(score_to_plot$score == min(score_to_plot$score),
                              paste0("\\textbf{", score_to_plot$score, "}"),
                              as.character(score_to_plot$score))

mse_to_plot <- cv_res$mse %>%
    filter(!is.na(mse)) %>%
    group_by(model_factor) %>%
    summarise(mse = round(mean(mse), 3)) %>%
    arrange(desc(mse))
mse_to_plot$mse <- ifelse(mse_to_plot$mse == min(mse_to_plot$mse),
                          paste0("\\textbf{", mse_to_plot$mse, "}"),
                          as.character(mse_to_plot$mse))

mse_to_plot %>%
    select(model_factor, mse) %>%
    kable(format = "markdown", caption = "MSE",
          col.names = c("Model", "MSE"))

score_to_plot %>%
    select(model_factor, score) %>%
    kable(format = "markdown", caption = "LogScore",
          col.names = c("Model", "$LogScore$"))
