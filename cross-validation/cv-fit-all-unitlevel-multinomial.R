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
library(survey)
library(svyVGAM)
library(mvtnorm)
library(rgdal)
library(INLA)
library(viridis)
library(classInt)
library(gridExtra)
library(ggpubr)

# functions ####
fitINLA <- function(formula, data) {
    inla(formula, data = data, 
         family = "poisson",
         control.family = list(link = "log"),
         control.predictor = list(link = 1),
         control.compute = list(config=T))
}


# read and format data ####
load(paste0(root,"Dropbox/dissertation_2/survey-csmf/data/ken_dhs2014/data/haz-waz-kenDHS2014.rda"))

dat_c$cont_factor <- factor(dat_c$contraceptive, 
                            levels = c("none", 
                                       "modern", 
                                       "other"))

# collapse, counting numbers of obs in each cluster
# add in 0 counts for clusters that don't havve obs
results.long <- dat_c %>% select(admin1, admin1.name, admin1.char,
                                 cluster, region, strata, rural, weights,
                                 contraceptive) %>%
    group_by(admin1, admin1.name, admin1.char,
             cluster, region, strata, rural, contraceptive) %>%
    summarize(count = n(), weights = unique(weights)) %>%
    ungroup() %>%
    complete(nesting(strata, cluster), contraceptive,
             fill = list(count = 0)) %>%
    group_by(strata, cluster) %>%
    mutate(admin1 = unique(admin1 %>% na.omit()), 
           admin1.name = unique(admin1.name %>% na.omit()), 
           admin1.char = unique(admin1.char %>% na.omit()), 
           region = unique(region %>% na.omit()), 
           rural = unique(rural %>% na.omit()), 
           weights = unique(weights %>% na.omit())) %>% 
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

# survey object for calculating direct estimates
my.svydesign <- survey::svydesign(ids = ~ cluster,
                                  strata = ~ strata, nest = T, weights = ~weights,
                                  data = dat_c)

score_res <- expand_grid(model = model_names, region = 1:n_regions, rural = 0:1) %>% 
    mutate(score = NA)

mse_res <- expand_grid(model = model_names, region = 1:n_regions, rural = 0:1) %>% 
    mutate(mse = NA)

# Cross validation ####
n_samples <- 250
# for (r in 1:n_regions) {
for (r in 28:n_regions) {

    message(paste0("region ", r))

    starttime <- Sys.time()
    
    ## Hold out data ####
    tmp <- data
    tmp$count[tmp$admin1 == r] <- NA
    
    # direct estimates for held out data
    heldout_urban <- subset(my.svydesign, admin1 == r & rural == 0)
    heldout_rural <- subset(my.svydesign, admin1 == r & rural == 1)
    
    mod_urban <- svy_vglm(cont_factor ~ 1, design = heldout_urban, family = multinomial(refLevel = 1))
    direct_logits_urban <- coef(mod_urban)
    V_urban <- vcov(mod_urban)
    if (!(r %in% c(28, 30))) {
        mod_rural <- svy_vglm(cont_factor ~ 1, design = heldout_rural, family = multinomial(refLevel = 1))
        direct_logits_rural <- coef(mod_rural)
        V_rural <- vcov(mod_rural)
    }

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
        none_fe_idx <- which(rownames(samp[[1]]$latent) %in% c("contraceptive_factornone:1"))
        none_fe_mat <- matrix(NA, nrow = length(none_fe_idx), ncol = n_samples)
        modern_fe_idx <- which(rownames(samp[[1]]$latent) %in% c("contraceptive_factormodern:1"))
        modern_fe_mat <- matrix(NA, nrow = length(modern_fe_idx), ncol = n_samples)
        other_fe_idx <- which(rownames(samp[[1]]$latent) %in% c("contraceptive_factorother:1"))
        other_fe_mat <- matrix(NA, nrow = length(other_fe_idx), ncol = n_samples)
        
        none_rural_fe_idx <- which(rownames(samp[[1]]$latent) %in% c("rural:1"))
        none_rural_fe_mat <- matrix(NA, nrow = length(none_rural_fe_idx), ncol = n_samples)
        modern_rural_fe_idx <- which(rownames(samp[[1]]$latent) %in% c("contraceptive_factormodern:rural:1"))
        modern_rural_fe_mat <- matrix(NA, nrow = length(modern_rural_fe_idx), ncol = n_samples)
        other_rural_fe_idx <- which(rownames(samp[[1]]$latent) %in% c("contraceptive_factorother:rural:1"))
        other_rural_fe_mat <- matrix(NA, nrow = length(other_rural_fe_idx), ncol = n_samples)
        
        modern_re_idx <- rownames(samp[[1]]$latent) %>% str_detect("admin1\\.modern:") %>% which()
        none_re_idx <- rownames(samp[[1]]$latent) %>% str_detect("admin1\\.none:") %>% which()
        modern_re_mat <- matrix(NA, nrow = length(modern_re_idx), ncol = n_samples)
        none_re_mat <- matrix(NA, nrow = length(none_re_idx), ncol = n_samples)
        
        if(!grepl("nonshared", model_names[mm])) {
            shared_re_idx <- rownames(samp[[1]]$latent) %>% str_detect("admin1\\.none\\.2:") %>% which()
            shared_re_mat <- matrix(NA, nrow = length(shared_re_idx), ncol = n_samples)
        } else {
            shared_re_idx <- integer(0)
            shared_re_mat <- matrix(0, nrow = length(none_re_idx), ncol = n_samples)
        }
        
        # fill in sample matrices
        for (s in 1:n_samples) {
            none_fe_mat[,s] <- samp[[s]]$latent[none_fe_idx]
            modern_fe_mat[,s] <- samp[[s]]$latent[modern_fe_idx]
            other_fe_mat[,s] <- samp[[s]]$latent[other_fe_idx]
            none_rural_fe_mat[,s] <- samp[[s]]$latent[none_rural_fe_idx]
            modern_rural_fe_mat[,s] <- samp[[s]]$latent[modern_rural_fe_idx]
            other_rural_fe_mat[,s] <- samp[[s]]$latent[other_rural_fe_idx]
            modern_re_mat[,s] <- samp[[s]]$latent[modern_re_idx]
            none_re_mat[,s] <- samp[[s]]$latent[none_re_idx]
            if(!grepl("nonshared", model_names[mm])) shared_re_mat[,s] <- samp[[s]]$latent[shared_re_idx]
        }
        
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
        fitted_m_n_urban <- fitted_modern_urban_mat - fitted_none_urban_mat
        fitted_o_n_urban <- fitted_other_urban_mat - fitted_none_urban_mat
        fitted_m_n_rural <- fitted_modern_rural_mat - fitted_none_rural_mat
        fitted_o_n_rural <- fitted_other_rural_mat - fitted_none_rural_mat
        
        # calculate score and mse results
        y_lik_urban <- c()
        y_lik_rural <- c()
        sq_error_urban <- c()
        sq_error_rural <- c()
        
        for (s in 1:n_samples) {
            
            if (length(direct_logits_urban) != 1) {
                post_mean_urban <- c(fitted_m_n_urban[r, s],
                                     fitted_o_n_urban[r, s])
                y_lik_urban <- c(y_lik_urban,
                                 dmvnorm(as.numeric(direct_logits_urban),
                                         mean = post_mean_urban,
                                         sigma = V_urban))
                sq_error_urban <- c(sq_error_urban,
                                    (as.numeric(direct_logits_urban) - post_mean_urban)^2)
            }
            #  completely urban regions can only get urban ests
            if (!(r %in% c(28, 30))) {
                if (length(direct_logits_rural) != 1) {
                    post_mean_rural <- c(fitted_m_n_rural[r, s],
                                         fitted_o_n_rural[r, s])
                    y_lik_rural <- c(y_lik_rural,
                                     dmvnorm(as.numeric(direct_logits_rural),
                                             mean = post_mean_rural,
                                             sigma = V_rural))
                    sq_error_rural <- c(sq_error_rural,
                                        (as.numeric(direct_logits_rural) - post_mean_rural)^2)
                }
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
write_rds(cv_res_list, file = "../../../Dropbox/dissertation_2/survey-csmf/results/cv/cv_results-multinomial-unitlevel.rds")

# if on Box, copy to real dropbox
if (root == "P:/") {
    write_rds(cv_res_list, file = "C:/Users/aeschuma/Dropbox/dissertation_2/survey-csmf/results/cv/cv_results-multinomial-unitlevel.rds")
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
