## ----setup, include=FALSE--------------------------------------------------------------------
library(SUMMER)
library(tidyverse)
library(spdep)
library(geosphere)
library(haven)
library(knitr)
library(kableExtra)
library(magrittr)
library(rstan)
library(cmdstanr)
library(svyVGAM)
library(mvtnorm)
library(rgdal)
library(bayesplot)
library(INLA)
library(viridis)
library(classInt)
library(gridExtra)
library(ggpubr)
library(loo)

## ----read-data-----------------------------------------------------------------------
load("/Users/austin/Dropbox/dissertation_2/survey-csmf/data/ken_dhs2014/data/haz-waz-kenDHS2014.rda")

dat_orig <- dat_c

## format data
dat_c$cont_factor <- factor(dat_c$contraceptive, 
                                 levels = c("none", 
                                            "modern", 
                                            "other"))
outcomes <- unique(dat_c$cont_factor)
noutcomes <- length(outcomes)

## add in a single row with lowest sampling weight in regions with 0 observations
## in order to be able to fit a multinomial model
lowest_weight <- min(dat_c$weights)
admin_outcome_table <- table(dat_c$cont_factor, dat_c$admin1)
admin_outcome_zero_ind <- tibble(admin1 = numeric(),
                                 cont_factor = character())
for (i in 1:n_regions) {
    for (j in 1:noutcomes) {
        if(admin_outcome_table[j, i] == 0) {
            admin_outcome_zero_ind <- bind_rows(admin_outcome_zero_ind, 
                                                tibble(admin1 = i,
                                                       cont_factor = outcomes[j]))
            additional_row <- dat_c %>% filter(admin1 == i) %>% slice(1)
            additional_row$cont_factor[1] <- outcomes[j]
            additional_row$weights[1] <- lowest_weight/2
            dat_c %<>% bind_rows(additional_row)
            }
    }
}
dat_c %<>% select(cluster, region, strata, urban_rural, weights, tmp_admin1,
                  v002, v003, LATNUM, LONGNUM, admin1, admin1.char, admin1.name, cont_factor)

## ----first-stage-----------------------------------------------------------------------------
## Stage 1 direct: univariate and bivariate ####
my.svydesign <- survey::svydesign(ids = ~ cluster,
                                  strata = ~ strata, nest = T, weights = ~weights,
                                  data = dat_c)

n_regions <- length(unique(dat_c$admin1.char))
admin1v <- unique(dat_c$admin1.char)
results <- data.frame(admin1 = unique(dat_c$admin1),
                      admin1.name = unique(dat_c$admin1.name),
                      admin1.char = admin1v,
                      mean_beta1 = rep(NA, n_regions),
                      mean_beta2 = rep(NA, n_regions),
                      mean_beta1_flip = rep(NA, n_regions),
                      mean_beta2_flip = rep(NA, n_regions),
                      se_beta1 = rep(NA, n_regions),
                      se_beta2 = rep(NA, n_regions),
                      corr = rep(NA, n_regions),
                      se_beta1_flip = rep(NA, n_regions),
                      se_beta2_flip = rep(NA, n_regions),
                      corr_flip = rep(NA, n_regions))
V.list <- vector(mode = "list", length = n_regions)
V.array <- array(NA, dim = c(nrow(results), 2, 2))
V.array.flip <- array(NA, dim = c(nrow(results), 2, 2))
names(V.list) <- admin1v
V.list.flip <- V.list

V.array.2 <- array(0, dim = c(nrow(results), 2, 2))

for(i in 1:n_regions) {
    admin1.tmp <- admin1v[i]
    
    tmp <- subset(my.svydesign, admin1.char == admin1.tmp)
    
    # test svy_vglm
    mod <- svy_vglm(cont_factor ~ 1, design = tmp, family = multinomial(refLevel = 1))
    mod.flip <- svy_vglm(cont_factor ~ 1, design = tmp, family = multinomial(refLevel = 2))
    
    index.tmp <- results$admin1.char == admin1.tmp
    
    results$mean_beta1[index.tmp] <- coef(mod)[1]
    results$mean_beta2[index.tmp] <- coef(mod)[2]
    
    results$mean_beta1_flip[index.tmp] <- coef(mod.flip)[1]
    results$mean_beta2_flip[index.tmp] <- coef(mod.flip)[2]
    
    V.tmp <- vcov(mod)
    V.tmp.flip <- vcov(mod.flip)
    V.list[[admin1.tmp]] <- V.tmp
    V.array[i, , ] <- V.tmp
    
    V.array.flip[i, , ] <- V.tmp.flip
    V.list.flip[[admin1.tmp]] <- V.tmp.flip
    
    results$se_beta1[index.tmp] <- V.tmp[1, 1]^0.5
    results$se_beta2[index.tmp] <- V.tmp[2, 2]^0.5
    
    results$se_beta1_flip[index.tmp] <- V.tmp.flip[1, 1]^0.5
    results$se_beta2_flip[index.tmp] <- V.tmp.flip[2, 2]^0.5
    
    D <- diag(sqrt(diag(V.tmp)))
    DInv <- solve(D)
    corr.tmp <- DInv %*% V.tmp %*% DInv
    
    results$corr[index.tmp] <- corr.tmp[1, 2]
    
    D <- diag(sqrt(diag(V.tmp.flip)))
    DInv <- solve(D)
    corr.tmp <- DInv %*% V.tmp.flip %*% DInv
    
    results$corr_flip[index.tmp] <- corr.tmp[1, 2]
}

stage_1_list <- list(results = results,
                     V.array = V.array,
                     V.list = V.list,
                     V.array.flip = V.array.flip,
                     V.list.flip = V.list.flip)

# save stage 1 results
write_rds(stage_1_list, file = "../../../Dropbox/dissertation_2/survey-csmf/results/ken2014-multinomial/ken2014-multinomial-stage-1.rds")
