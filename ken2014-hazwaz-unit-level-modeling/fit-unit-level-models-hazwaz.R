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
admin_df <- dat %>% select(admin1.name, admin1) %>% distinct()
prop_urban %<>% left_join(admin_df, 
                          by = c("adm_name" = "admin1.name")) %>%
    select(admin1, urb_frac)

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

# person ID
results.long <- results.long %>% group_by(admin1, cluster) %>% mutate(pid = 1:n()) %>% ungroup()

# create a list of the data
data <- as.list(results.long)
data.haz <- results.long %>% filter(outcome == "HAZ") %>% as.list()
data.waz <- results.long %>% filter(outcome == "WAZ") %>% as.list()

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

# formulas ####
formulas <- list()

# formulas[["IID univariate"]] <- formula("value ~ 1 + rural +
#                                         f(admin1, model = 'iid', hyper = iid_prior) +
#                                         f(cluster, model = 'iid', hyper = iid_prior)")
# 
# formulas[["BYM univariate"]] <- formula("value ~ 1 + rural + 
#                                         f(admin1, model = 'bym2',
#                                              graph = admin1.mat, 
#                                              scale.model = T, 
#                                              constr = T,
#                                              hyper = bym2_prior) +
#                                         f(cluster, model = 'iid', hyper = iid_prior)")
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
model_names <- names(formulas)

# modeling ####
model_results <- vector(mode = "list", length = length(model_names))
names(model_results) <- model_names

model_fits <- vector(mode = "list", length = length(model_names))
names(model_fits) <- model_names

nsamps <- 1000

for (i in 1:length(model_names)) {
    start_time <- Sys.time()
    
    tmp <- inla(formulas[[ model_names[i] ]], data = data,
                family = rep("gaussian",2),
                quantiles=c(0.025, 0.1, 0.5, 0.9, 0.975),
                # control.family = list(),
                control.predictor = list(compute=T),
                # control.inla = list(int.strategy = "grid"),
                # control.fixed = list(prec = list(default = 0.001), correlation.matrix=T),
                control.compute = list(config=T, waic = TRUE, dic = TRUE, cpo = TRUE),
                selection = list(`Precision for the Gaussian observations`, `Precision for cluster.haz`))
    
    # fit and plot univariate models if we want to
    # if (grepl("univariate", model_names[i])) {
    #     tmp.haz <- inla(formulas[[ model_names[i] ]], data = data.haz,
    #                     family = "gaussian",
    #                     quantiles=c(0.025, 0.1, 0.5, 0.9, 0.975),
    #                     # control.family = list(),
    #                     control.predictor = list(compute=T),
    #                     # control.inla = list(lincomb.derived.correlation.matrix=T),
    #                     # control.fixed = list(prec = list(default = 0.001), correlation.matrix=T),
    #                     control.compute = list(config=T, waic = TRUE, dic = TRUE, cpo = TRUE))
    #     tmp.waz <- inla(formulas[[ model_names[i] ]], data = data.waz,
    #                     family = "gaussian",
    #                     quantiles=c(0.025, 0.1, 0.5, 0.9, 0.975),
    #                     # control.family = list(),
    #                     control.predictor = list(compute=T),
    #                     # control.inla = list(lincomb.derived.correlation.matrix=T),
    #                     # control.fixed = list(prec = list(default = 0.001), correlation.matrix=T),
    #                     control.compute = list(config=T, waic = TRUE, dic = TRUE, cpo = TRUE))
    #     tmp.uni.res.haz <- tibble(admin1 = data.haz$admin1,
    #                               cluster = data.haz$cluster,
    #                               pid = data.haz$pid,
    #                               outcome = "HAZ",
    #                               pred.med.uni = tmp.haz$summary.fitted.values$`0.5quant`,
    #                               pred.lwr.uni = tmp.haz$summary.fitted.values$`0.025quant`,
    #                               pred.upr.uni = tmp.haz$summary.fitted.values$`0.975quant`)
    #     tmp.uni.res.waz <- tibble(admin1 = data.waz$admin1,
    #                               cluster = data.waz$cluster,
    #                               pid = data.waz$pid,
    #                               outcome = "WAZ",
    #                               pred.med.uni = tmp.waz$summary.fitted.values$`0.5quant`,
    #                               pred.lwr.uni = tmp.waz$summary.fitted.values$`0.025quant`,
    #                               pred.upr.uni = tmp.waz$summary.fitted.values$`0.975quant`)
    #     tmp.bi.res <- tibble(admin1 = data$admin1,
    #                          cluster = data$cluster,
    #                          pid = data$pid,
    #                          outcome = data$outcome,
    #                          pred.med.bi = tmp$summary.fitted.values$`0.5quant`,
    #                          pred.lwr.bi = tmp$summary.fitted.values$`0.025quant`,
    #                          pred.upr.bi = tmp$summary.fitted.values$`0.975quant`)
    #     
    #     tmp.res <- tmp.uni.res.haz %>% bind_rows(tmp.uni.res.waz) %>% left_join(tmp.bi.res)
    #     
    #     #plot comparison
    #     ggplot(tmp.res, aes(x = pred.med.bi, y = pred.med.uni)) + geom_point(alpha = 0.5) + geom_abline(slope = 1, intercept = 0, col = "darkgreen", alpha = 0.7)
    #     ggplot(tmp.res, aes(x = pred.lwr.bi, y = pred.lwr.uni)) + geom_point(alpha = 0.5) + geom_abline(slope = 1, intercept = 0, col = "darkgreen", alpha = 0.7)
    #     ggplot(tmp.res, aes(x = pred.upr.bi, y = pred.upr.uni)) + geom_point(alpha = 0.5) + geom_abline(slope = 1, intercept = 0, col = "darkgreen", alpha = 0.7)
    #     
    # }
    
    # extract and format results
    samp <- inla.posterior.sample(n = nsamps, result = tmp, use.improved.mean = TRUE)
    samp2 <- inla.hyperpar.sample(n = nsamps, result = tmp)
    
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
    
    # pdf("bivariate-ind-cluster-precision-hyperpar-plot.pdf", width = 8, height = 4)
    # par(mfrow = c(1, 2))
    # plot(samp2[,"Precision for the Gaussian observations"], 
    #      samp2[,"Precision for cluster.haz"],
    #      xlab = "Gaussian precision",
    #      ylab = "Cluster precision", 
    #      main = "HAZ")
    # plot(samp2[,"Precision for the Gaussian observations[2]"], 
    #      samp2[,"Precision for cluster.waz"],
    #      xlab = "Gaussian precision",
    #      ylab = "Cluster precision", 
    #      main = "WAZ")
    # dev.off()
    # 
    # process fitted values
    haz_fe_idx <- rownames(samp[[1]]$latent) %>% str_detect("outcomeHAZ:1") %>% which()
    waz_fe_idx <- rownames(samp[[1]]$latent) %>% str_detect("outcomeWAZ:1") %>% which()
    haz_fe_mat <- matrix(NA, nrow = length(haz_fe_idx), ncol = nsamps)
    waz_fe_mat <- matrix(NA, nrow = length(waz_fe_idx), ncol = nsamps)
    
    haz_rural_fe_idx <- rownames(samp[[1]]$latent) %>% str_detect("rural\\.haz") %>% which()
    waz_rural_fe_idx <- rownames(samp[[1]]$latent) %>% str_detect("rural\\.waz") %>% which()
    haz_rural_fe_mat <- matrix(NA, nrow = length(haz_rural_fe_idx), ncol = nsamps)
    waz_rural_fe_mat <- matrix(NA, nrow = length(waz_rural_fe_idx), ncol = nsamps)
    
    haz_re_idx <- rownames(samp[[1]]$latent) %>% str_detect("admin1\\.haz:") %>% which()
    waz_re_idx <- rownames(samp[[1]]$latent) %>% str_detect("admin1\\.waz:") %>% which()
    haz_re_mat <- matrix(NA, nrow = length(haz_re_idx), ncol = nsamps)
    waz_re_mat <- matrix(NA, nrow = length(waz_re_idx), ncol = nsamps)
    
    if(!grepl("nonshared", model_names[i])) {
        shared_re_idx <- rownames(samp[[1]]$latent) %>% str_detect("admin1\\.haz\\.2:") %>% which()
        shared_re_mat <- matrix(NA, nrow = length(shared_re_idx), ncol = nsamps)
    } else {
        shared_re_idx <- integer(0)
        shared_re_mat <- matrix(0, nrow = nrow(haz_re_mat), ncol = nsamps)
    }

    # fill in sample matrices
    for (s in 1:nsamps) {
        haz_fe_mat[,s] <- samp[[s]]$latent[haz_fe_idx]
        waz_fe_mat[,s] <- samp[[s]]$latent[waz_fe_idx]
        haz_rural_fe_mat[,s] <- samp[[s]]$latent[haz_rural_fe_idx]
        waz_rural_fe_mat[,s] <- samp[[s]]$latent[waz_rural_fe_idx]
        haz_re_mat[,s] <- samp[[s]]$latent[haz_re_idx]
        waz_re_mat[,s] <- samp[[s]]$latent[waz_re_idx]
        if(!grepl("nonshared", model_names[i])) shared_re_mat[,s] <- samp[[s]]$latent[shared_re_idx]
    }
    
    # save sample matrices in list for output
    model_results[[ model_names[i] ]]$samples_all <- list()
    model_results[[ model_names[i] ]]$samples_all$haz_fe <- haz_fe_mat
    model_results[[ model_names[i] ]]$samples_all$waz_fe <- waz_fe_mat
    model_results[[ model_names[i] ]]$samples_all$haz_rural_fe <- haz_rural_fe_mat
    model_results[[ model_names[i] ]]$samples_all$waz_rural_fe <- waz_rural_fe_mat
    model_results[[ model_names[i] ]]$samples_all$haz_re <- haz_re_mat
    model_results[[ model_names[i] ]]$samples_all$waz_re <- waz_re_mat
    model_results[[ model_names[i] ]]$samples_all$shared_re <- shared_re_mat
    
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
    
    # output for region-level samples
    model_results[[ model_names[i] ]]$samples_region_ur <- list()
    model_results[[ model_names[i] ]]$samples_region_ur$haz_urban <- fitted_haz_urban_mat
    model_results[[ model_names[i] ]]$samples_region_ur$waz_urban <- fitted_waz_urban_mat
    model_results[[ model_names[i] ]]$samples_region_ur$haz_rural <- fitted_haz_rural_mat
    model_results[[ model_names[i] ]]$samples_region_ur$waz_rural <- fitted_waz_rural_mat
    
    # aggregate over urban/rural
    fitted_haz_mat <- matrix(NA, nrow = n_regions, ncol = nsamps)
    fitted_waz_mat <- matrix(NA, nrow = n_regions, ncol = nsamps)
    for (s in 1:nsamps)  {
        fitted_haz_mat[, s] <- (fitted_haz_urban_mat[, s] * prop_urban$urb_frac) + 
            (fitted_haz_rural_mat[, s] * (1 - prop_urban$urb_frac))
        fitted_waz_mat[, s] <- (fitted_waz_urban_mat[, s] * prop_urban$urb_frac) + 
            (fitted_waz_rural_mat[, s] * (1 - prop_urban$urb_frac))
    }
    
    # store results
    model_results[[ model_names[i] ]]$samples_final <- list()
    model_results[[ model_names[i] ]]$samples_final$haz <- fitted_haz_mat
    model_results[[ model_names[i] ]]$samples_final$waz <- fitted_waz_mat
    
    # run time
    model_results[[ model_names[i] ]]$run_time <- Sys.time() - start_time
    
    # store full results
    model_fits[[ model_names[i] ]] <- tmp
    rm(tmp)
}

# model results ####
write_rds(model_results, file = "../../../Dropbox/dissertation_2/survey-csmf/results/ken2014-unit-level/hazwaz/ken2014-unit-level-hazwaz-inla-posteriors.rds")
write_rds(model_fits, file = "../../../Dropbox/dissertation_2/survey-csmf/results/ken2014-unit-level/hazwaz/ken2014-unit-level-hazwaz-inla-fits.rds")

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
