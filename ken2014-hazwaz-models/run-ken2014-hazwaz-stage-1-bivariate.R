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

## ----first-stage-----------------------------------------------------------------------------
## Stage 1 direct: univariate and bivariate ####
my.svydesign <- survey::svydesign(ids = ~ cluster,
                                  strata = ~ strata, nest = T, weights = ~weights,
                                  data = dat)

n_regions <- length(unique(dat$admin1.char))
admin1v <- unique(dat$admin1.char)
results <- data.frame(admin1 = unique(dat$admin1),
                      admin1.name = unique(dat$admin1.name),
                      admin1.char = admin1v,
                      meanHAZ.bi = rep(NA, n_regions),
                      meanWAZ.bi = rep(NA, n_regions),
                      seHAZ.bi = rep(NA, n_regions),
                      seWAZ.bi = rep(NA, n_regions),
                      corr.bi = rep(NA, n_regions))
V.list <- vector(mode = "list", length = n_regions)
V.array <- array(NA, dim = c(nrow(results), 2, 2))
V.array.flip <- array(NA, dim = c(nrow(results), 2, 2))
names(V.list) <- admin1v
V.list.flip <- V.list


V.array.2 <- array(0, dim = c(nrow(results), 2, 2))

for(i in 1:n_regions) {
    admin1.tmp <- admin1v[i]
    
    tmp <- subset(my.svydesign, admin1.char == admin1.tmp)
    
    means.svymean <- svymean(~HAZ + WAZ, tmp)
    means.svymean.flip <- svymean(~WAZ + HAZ, tmp)
    
    index.tmp <- results$admin1.char == admin1.tmp
    
    results$meanHAZ.bi[index.tmp] <- means.svymean[["HAZ"]]
    results$meanWAZ.bi[index.tmp] <- means.svymean[["WAZ"]]
    
    V.tmp <- vcov(means.svymean)
    V.tmp.flip <- vcov(means.svymean.flip)
    V.list[[admin1.tmp]] <- V.tmp
    V.array[i, , ] <- V.tmp
    
    V.array.flip[i, , ] <- V.tmp.flip
    V.list.flip[[admin1.tmp]] <- V.tmp.flip
    
    results$seHAZ.bi[index.tmp] <- V.tmp[1, 1]^0.5
    results$seWAZ.bi[index.tmp] <- V.tmp[2, 2]^0.5
    
    D <- diag(sqrt(diag(V.tmp)))
    DInv <- solve(D)
    corr.tmp <- DInv %*% V.tmp %*% DInv
    
    results$corr.bi[index.tmp] <- corr.tmp[1, 2]
}

stage_1_list <- list(results = results,
                     V.array = V.array,
                     V.list = V.list,
                     V.array.flip = V.array.flip,
                     V.list.flip = V.list.flip)

# save stage 1 results
write_rds(stage_1_list, file = "../../../Dropbox/dissertation_2/survey-csmf/results/ken2014-hazwaz/ken2014-hazwaz-stage-1.rds")
