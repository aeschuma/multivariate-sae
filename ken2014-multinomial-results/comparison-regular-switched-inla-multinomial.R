# compare regular and switched results

rm(list = ls())

library(SUMMER)
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
library(INLA)
library(viridis)
library(classInt)
library(gridExtra)
library(ggpubr)
library(gridExtra)

savedir <- "~/Dropbox/dissertation_2/survey-csmf/results/proj-2-chapter-results"

# load data ####
regular_results <- read_rds("../../../Dropbox/dissertation_2/survey-csmf/results/ken2014-multinomial/ken2014-multinomial-stage-2-inla-all.rds")
switched_results <- read_rds("../../../Dropbox/dissertation_2/survey-csmf/results/ken2014-multinomial/ken2014-multinomial-stage-2-sharedmodern-inla-all.rds")
n_regions <- 47

# compare latent means ####

## other shared posterior probs ####
fmod <- regular_results$`Bivariate shared BYM`
nsamps <- 1000
samp <- inla.posterior.sample(nsamps, fmod)

# process fitted values
modern_fe_idx <- which(rownames(samp[[1]]$latent) %in% c("outcomebeta1:1"))
modern_fe_mat <- matrix(NA, nrow = length(modern_fe_idx), ncol = nsamps)
other_fe_idx <- which(rownames(samp[[1]]$latent) %in% c("outcomebeta2:1"))
other_fe_mat <- matrix(NA, nrow = length(other_fe_idx), ncol = nsamps)

modern_re_idx <- rownames(samp[[1]]$latent) %>% str_detect("admin1\\.beta1:") %>% which()
other_re_idx <- rownames(samp[[1]]$latent) %>% str_detect("admin1\\.beta2:") %>% which()
modern_re_mat <- matrix(NA, nrow = length(modern_re_idx), ncol = nsamps)
other_re_mat <- matrix(NA, nrow = length(other_re_idx), ncol = nsamps)

shared_re_idx <- rownames(samp[[1]]$latent) %>% str_detect("admin1\\.beta1\\.2:") %>% which()
shared_re_mat <- matrix(NA, nrow = length(shared_re_idx), ncol = nsamps)

# fill in sample matrices
for (s in 1:nsamps) {
    modern_fe_mat[,s] <- samp[[s]]$latent[modern_fe_idx]
    other_fe_mat[,s] <- samp[[s]]$latent[other_fe_idx]
    modern_re_mat[,s] <- samp[[s]]$latent[modern_re_idx]
    other_re_mat[,s] <- samp[[s]]$latent[other_re_idx]
    shared_re_mat[,s] <- samp[[s]]$latent[shared_re_idx]
}

# obtain total effect for space - first half of bym2, or all of them for IID
modern_re_tot_mat <- modern_re_mat[1:n_regions,] + shared_re_mat[1:n_regions,]
other_re_tot_mat <- other_re_mat[1:n_regions,] 

# matrix of fitted estimates
fitted_modern_mat <- modern_fe_mat[rep(1,n_regions),] + modern_re_tot_mat
fitted_other_mat <- other_fe_mat[rep(1,n_regions),] + other_re_tot_mat

# transform to probs
fitted_probs <- array(NA, dim = c(3, n_regions, nsamps))
for (s in 1:nsamps) {
    for (r in 1:n_regions) {
        denom <- 1 + exp(fitted_modern_mat[r, s]) + exp(fitted_other_mat[r, s])
        fitted_probs[1, r, s] <- 1/denom
        fitted_probs[2, r, s] <- exp(fitted_modern_mat[r, s])/denom
        fitted_probs[3, r, s] <- exp(fitted_other_mat[r, s])/denom
    }
}

# summaries of probs
posterior_probs_othershared <- tibble(model = "other.shared",
                                      region = 1:n_regions,
                                      none_med = apply(fitted_probs[1,,], 1, median),
                                      modern_med = apply(fitted_probs[2,,], 1, median),
                                      other_med = apply(fitted_probs[3,,], 1, median),
                                      none_width95 = apply(fitted_probs[1,,], 1, function(x) { sort(x)[975] - sort(x)[25] }),
                                      modern_width95 = apply(fitted_probs[2,,], 1, function(x) { sort(x)[975] - sort(x)[25] }),
                                      other_width95 = apply(fitted_probs[3,,], 1, function(x) { sort(x)[975] - sort(x)[25] }))

## modern shared posterior probs ####
fmod2 <- switched_results$`Bivariate shared BYM`
samp2 <- inla.posterior.sample(nsamps, fmod2)

# process fitted values
modern_fe_idx <- which(rownames(samp2[[1]]$latent) %in% c("outcomebeta1:1"))
modern_fe_mat <- matrix(NA, nrow = length(modern_fe_idx), ncol = nsamps)
other_fe_idx <- which(rownames(samp2[[1]]$latent) %in% c("outcomebeta2:1"))
other_fe_mat <- matrix(NA, nrow = length(other_fe_idx), ncol = nsamps)

modern_re_idx <- rownames(samp2[[1]]$latent) %>% str_detect("admin1\\.beta1:") %>% which()
other_re_idx <- rownames(samp2[[1]]$latent) %>% str_detect("admin1\\.beta2:") %>% which()
modern_re_mat <- matrix(NA, nrow = length(modern_re_idx), ncol = nsamps)
other_re_mat <- matrix(NA, nrow = length(other_re_idx), ncol = nsamps)

shared_re_idx <- rownames(samp2[[1]]$latent) %>% str_detect("admin1\\.beta2\\.2:") %>% which()
shared_re_mat <- matrix(NA, nrow = length(shared_re_idx), ncol = nsamps)

# fill in sample matrices
for (s in 1:nsamps) {
    modern_fe_mat[,s] <- samp2[[s]]$latent[modern_fe_idx]
    other_fe_mat[,s] <- samp2[[s]]$latent[other_fe_idx]
    modern_re_mat[,s] <- samp2[[s]]$latent[modern_re_idx]
    other_re_mat[,s] <- samp2[[s]]$latent[other_re_idx]
    shared_re_mat[,s] <- samp2[[s]]$latent[shared_re_idx]
}

# obtain total effect for space - first half of bym2, or all of them for IID
modern_re_tot_mat <- modern_re_mat[1:n_regions,] 
other_re_tot_mat <- other_re_mat[1:n_regions,] + shared_re_mat[1:n_regions,]

# matrix of fitted estimates
fitted_modern_mat <- modern_fe_mat[rep(1,n_regions),] + modern_re_tot_mat
fitted_other_mat <- other_fe_mat[rep(1,n_regions),] + other_re_tot_mat

# transform to probs
fitted_probs <- array(NA, dim = c(3, n_regions, nsamps))
for (s in 1:nsamps) {
    for (r in 1:n_regions) {
        denom <- 1 + exp(fitted_modern_mat[r, s]) + exp(fitted_other_mat[r, s])
        fitted_probs[1, r, s] <- 1/denom
        fitted_probs[2, r, s] <- exp(fitted_modern_mat[r, s])/denom
        fitted_probs[3, r, s] <- exp(fitted_other_mat[r, s])/denom
    }
}

# summaries of probs
posterior_probs_modernshared <- tibble(model = "modern.shared",
                                       region = 1:n_regions,
                                       none_med = apply(fitted_probs[1,,], 1, median),
                                       modern_med = apply(fitted_probs[2,,], 1, median),
                                       other_med = apply(fitted_probs[3,,], 1, median),
                                       none_width95 = apply(fitted_probs[1,,], 1, function(x) { sort(x)[975] - sort(x)[25] }),
                                       modern_width95 = apply(fitted_probs[2,,], 1, function(x) { sort(x)[975] - sort(x)[25] }),
                                       other_width95 = apply(fitted_probs[3,,], 1, function(x) { sort(x)[975] - sort(x)[25] }))

## compile results ####
posterior_probs <- bind_rows(posterior_probs_othershared,
                             posterior_probs_modernshared) %>%
    pivot_longer(cols = ends_with(c("med", "width95")),
                 names_to = c("category", "measure"),
                 names_sep = "_",
                 values_to = "value") %>%
    pivot_wider(names_from = "model",
                values_from = "value") %>%
    mutate(measure = ifelse(measure == "med", "Median",
                            ifelse(measure == "width95", "95% CI width", NA))) %>%
    mutate(measure = factor(measure, levels = c("Median", "95% CI width")),
           category = factor(category, levels = c("none", "modern", "other"))) 

## plots ####
ggplot(posterior_probs, aes(x = other.shared, y = modern.shared)) +
    geom_abline(slope = 1, intercept = 0, col = "darkgreen", alpha = 0.7) +
    geom_point(alpha = 0.5) +
    facet_wrap(~ measure + category, nrow = 2, ncol = 3, scales = "free") +
    xlab("Other shared (in chapter)") +
    ylab("Modern shared (alternate)") +
    theme_light()
ggsave(paste0(savedir, "/shared_bym_multinomial_different_shared_comparison.pdf"), height = 6, width = 9)
