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
library(surveillance)

savedir <- "~/Dropbox/dissertation_2/survey-csmf/results/proj-2-chapter-results"

## ----read-data-----------------------------------------------------------------------
load("~/Dropbox/dissertation_2/survey-csmf/data/ken_dhs2014/data/haz-waz-kenDHS2014.rda")

## ----first-stage survey weighted--------------------------------------------------------------------
## Stage 1 direct: univariate and bivariate ####

# unweighted survey design
my.svydesign.un <- survey::svydesign(ids = ~ 1,
                                     data = dat)

# weighted survey design
my.svydesign <- survey::svydesign(ids = ~ cluster,
                                  strata = ~ strata, nest = T, weights = ~weights,
                                  data = dat)

n_regions <- length(unique(dat$admin1.char))
admin1v <- unique(dat$admin1.char)
results <- data.frame(admin1 = unique(dat$admin1),
                      admin1.name = unique(dat$admin1.name),
                      admin1.char = admin1v,
                      meanHAZuni.unweighted = rep(NA, n_regions),
                      meanWAZuni.unweighted = rep(NA, n_regions),
                      meanHAZuni.weighted = rep(NA, n_regions),
                      meanWAZuni.weighted = rep(NA, n_regions),
                      seHAZuni.unweighted = rep(NA, n_regions),
                      seWAZuni.unweighted = rep(NA, n_regions),
                      seHAZuni.weighted = rep(NA, n_regions),
                      seWAZuni.weighted = rep(NA, n_regions),
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

# national
nat_res <- data.frame(meanHAZuni.unweighted = rep(NA, 1),
                      meanWAZuni.unweighted = rep(NA, 1),
                      meanHAZuni.weighted = rep(NA, 1),
                      meanWAZuni.weighted = rep(NA, 1),
                      seHAZuni.unweighted = rep(NA, 1),
                      seWAZuni.unweighted = rep(NA, 1),
                      seHAZuni.weighted = rep(NA, 1),
                      seWAZuni.weighted = rep(NA, 1),
                      meanHAZ.bi = rep(NA, 1),
                      meanWAZ.bi = rep(NA, 1),
                      seHAZ.bi = rep(NA, 1),
                      seWAZ.bi = rep(NA, 1),
                      corr.bi = rep(NA, 1))

# unweighted
means.svymean.un <- svymean(~HAZ + WAZ, my.svydesign.un)

# weighted
means.svymean <- svymean(~HAZ + WAZ, my.svydesign)
means.svymean.flip <- svymean(~WAZ + HAZ, my.svydesign)

nat_res$meanHAZ.bi <- means.svymean[["HAZ"]]
nat_res$meanWAZ.bi <- means.svymean[["WAZ"]]
nat_res$meanHAZuni.weighted <- means.svymean[["HAZ"]]
nat_res$meanWAZuni.weighted <- means.svymean[["WAZ"]]
nat_res$meanHAZuni.unweighted <- means.svymean.un[["HAZ"]]
nat_res$meanWAZuni.unweighted <- means.svymean.un[["WAZ"]]

V.tmp <- vcov(means.svymean)
V.tmp.flip <- vcov(means.svymean.flip)

nat_res$seHAZ.bi <- V.tmp[1, 1]^0.5
nat_res$seWAZ.bi <- V.tmp[2, 2]^0.5
nat_res$seHAZuni.weighted<- sqrt(diag(vcov(means.svymean)))[1]
nat_res$seWAZuni.weighted <- sqrt(diag(vcov(means.svymean)))[2]
nat_res$seHAZuni.unweighted <- sqrt(diag(vcov(means.svymean.un)))[1]
nat_res$seWAZuni.unweighted <- sqrt(diag(vcov(means.svymean.un)))[2]

D <- diag(sqrt(diag(V.tmp)))
DInv <- solve(D)
corr.tmp <- DInv %*% V.tmp %*% DInv

nat_res$corr.bi <- corr.tmp[1, 2]

for(i in 1:n_regions) {
    admin1.tmp <- admin1v[i]
    
    # unweighted
    tmp.un <- subset(my.svydesign.un, admin1.char == admin1.tmp)
    means.svymean.un <- svymean(~HAZ + WAZ, tmp.un)
    
    # weighted
    tmp <- subset(my.svydesign, admin1.char == admin1.tmp)
    means.svymean <- svymean(~HAZ + WAZ, tmp)
    means.svymean.flip <- svymean(~WAZ + HAZ, tmp)
    
    index.tmp <- results$admin1.char == admin1.tmp
    
    results$meanHAZ.bi[index.tmp] <- means.svymean[["HAZ"]]
    results$meanWAZ.bi[index.tmp] <- means.svymean[["WAZ"]]
    results$meanHAZuni.weighted[index.tmp] <- means.svymean[["HAZ"]]
    results$meanWAZuni.weighted[index.tmp] <- means.svymean[["WAZ"]]
    results$meanHAZuni.unweighted[index.tmp] <- means.svymean.un[["HAZ"]]
    results$meanWAZuni.unweighted[index.tmp] <- means.svymean.un[["WAZ"]]
    
    V.tmp <- vcov(means.svymean)
    V.tmp.flip <- vcov(means.svymean.flip)
    V.list[[admin1.tmp]] <- V.tmp
    V.array[i, , ] <- V.tmp
    
    V.array.flip[i, , ] <- V.tmp.flip
    V.list.flip[[admin1.tmp]] <- V.tmp.flip
    
    results$seHAZ.bi[index.tmp] <- V.tmp[1, 1]^0.5
    results$seWAZ.bi[index.tmp] <- V.tmp[2, 2]^0.5
    results$seHAZuni.weighted[index.tmp] <- sqrt(diag(vcov(means.svymean)))[1]
    results$seWAZuni.weighted[index.tmp] <- sqrt(diag(vcov(means.svymean)))[2]
    results$seHAZuni.unweighted[index.tmp] <- sqrt(diag(vcov(means.svymean.un)))[1]
    results$seWAZuni.unweighted[index.tmp] <- sqrt(diag(vcov(means.svymean.un)))[2]
    
    D <- diag(sqrt(diag(V.tmp)))
    DInv <- solve(D)
    corr.tmp <- DInv %*% V.tmp %*% DInv
    
    results$corr.bi[index.tmp] <- corr.tmp[1, 2]
}

# make plots of univariate estimates for data section

# polygon plots
n_cats <- 15
plot.palette <- viridis(n_cats)
outcomes.names <- c("Mean HAZ unweighted", "Mean WAZ unweighted", "Mean HAZ weighted", "Mean WAZ weighted")
outcomes <- c("meanHAZuni.unweighted", "meanWAZuni.unweighted", "meanHAZuni.weighted", "meanWAZuni.weighted")
mean_range <- range(c(results$meanHAZuni.weighted, results$meanHAZuni.unweighted, results$meanWAZuni.weighted, results$meanWAZuni.unweighted))

setwd(savedir)
pdf("weighted-unweighted-uni-mean-compare.pdf", width = 9, height = 9)
tmp <- merge(poly.adm1, results[, c("admin1.name", outcomes)], 
             by.x = "NAME_1", by.y = "admin1.name")
spplot(tmp, outcomes, 
       col.regions = plot.palette, cuts = length(plot.palette) - 1,
       layout = c(2, 2),
       names.attr = outcomes.names)
dev.off()

n_cats <- 15
plot.palette <- viridis(n_cats)
outcomes.names <- c("SE HAZ unweighted", "SE WAZ unweighted", "SE HAZ weighted", "SE WAZ weighted")
outcomes <- c("seHAZuni.unweighted", "seWAZuni.unweighted", "seHAZuni.weighted", "seWAZuni.weighted")
se_range <- range(c(results$seHAZuni.weighted, results$seHAZuni.unweighted, results$senWAZuni.weighted, results$seWAZuni.unweighted))

setwd(savedir)
pdf("weighted-unweighted-uni-se-compare.pdf", width = 9, height = 9)
tmp <- merge(poly.adm1, results[, c("admin1.name", outcomes)], 
             by.x = "NAME_1", by.y = "admin1.name")
spplot(tmp, outcomes, 
       col.regions = plot.palette, cuts = length(plot.palette) - 1,
       layout = c(2, 2),
       names.attr = outcomes.names)
dev.off()

# scatter plot with uncertainty
results$lowerHAZuni.weighted <- results$meanHAZuni.weighted - results$seHAZuni.weighted
results$upperHAZuni.weighted <- results$meanHAZuni.weighted + results$seHAZuni.weighted
results$lowerWAZuni.weighted <- results$meanWAZuni.weighted - results$seWAZuni.weighted
results$upperWAZuni.weighted <- results$meanWAZuni.weighted + results$seWAZuni.weighted

ggplot(results, aes(x = meanWAZuni.weighted, y = meanHAZuni.weighted,
                    ymin = lowerHAZuni.weighted, ymax = upperHAZuni.weighted,
                    xmin = lowerWAZuni.weighted, xmax = upperWAZuni.weighted)) +
    geom_point(cex = 1.5) +
    geom_errorbar() +
    geom_errorbarh() +
    geom_smooth(method = "lm") +
    ylab("HAZ") +
    xlab("WAZ") +
    theme_light()
ggsave("haz-waz-weighted-scatter.pdf", width = 6, height = 6)
