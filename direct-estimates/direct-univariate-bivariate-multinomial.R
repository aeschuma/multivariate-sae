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

stage_1_list <- read_rds("/Users/austin/Dropbox/dissertation_2/survey-csmf/results/ken2014-multinomial/ken2014-multinomial-stage-1.rds")
results <- stage_1_list[["results"]]
V.array <- stage_1_list[["V.array"]]
n_regions <- nrow(poly.adm1)

# make plots of univariate estimates for data section

# polygon plots
n_cats <- 15
plot.palette <- viridis(n_cats)
outcomes.names <- c("Mean beta1", "Mean beta2",
                    "SE beta1", "SE beta2")
outcomes <- c("mean_beta1", "mean_beta2", "se_beta1", "se_beta2")

setwd(savedir)
pdf("direct-mean-se-multinomial.pdf", width = 9, height = 9)
tmp <- merge(poly.adm1, results[, c("admin1.name", outcomes)], 
             by.x = "NAME_1", by.y = "admin1.name")
spplot(tmp, outcomes[grepl("mean", outcomes)], 
       col.regions = plot.palette, cuts = length(plot.palette) - 1,
       layout = c(1, 2),
       names.attr = outcomes.names[grepl("Mean", outcomes.names)])
spplot(tmp, outcomes[grepl("se_", outcomes)], 
       col.regions = plot.palette, cuts = length(plot.palette) - 1,
       layout = c(1, 2),
       names.attr = outcomes.names[grepl("SE", outcomes.names)])
dev.off()

# scatter plot with uncertainty
results$lower_beta1 <- results$mean_beta1 - results$se_beta1
results$upper_beta1 <- results$mean_beta1 + results$se_beta1
results$lower_beta2 <- results$mean_beta2 - results$se_beta2
results$upper_beta2 <- results$mean_beta2 + results$se_beta2

ggplot(results, aes(x = mean_beta2, y = mean_beta1,
                    ymin = lower_beta1, ymax = upper_beta1,
                    xmin = lower_beta2, xmax = upper_beta2)) +
    geom_point(cex = 1.5) +
    geom_errorbar() +
    geom_errorbarh() +
    geom_smooth(method = "lm") +
    ylab("beta1") +
    xlab("beta2") +
    theme_light()
ggsave("multinomial-weighted-scatter.pdf", width = 6, height = 6)
