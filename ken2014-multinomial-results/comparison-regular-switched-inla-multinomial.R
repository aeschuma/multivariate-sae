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
library(gridExtra)


regular_results <- read_rds("../../../Dropbox/dissertation_2/survey-csmf/results/ken2014-hazwaz/ken2014-hazwaz-stage-2-inla-all.rds")
switched_results <- read_rds("../../../Dropbox/dissertation_2/survey-csmf/results/ken2014-hazwaz/ken2014-hazwaz-stage-2-switched-inla-all.rds")

# compare latent means

## haz

reg_res_haz <- regular_results$`Bivariate shared BYM`$summary.lincomb.derived[1:47, "mean"]
switched_res_haz <- switched_results$`Bivariate shared BYM`$summary.lincomb.derived[48:94, "mean"]

plot(reg_res_haz, switched_res_haz)
abline(0, 1, col = "darkgreen")

reg_res_haz95l <- regular_results$`Bivariate shared BYM`$summary.lincomb.derived[1:47, "0.025quant"]
switched_res_haz95l <- switched_results$`Bivariate shared BYM`$summary.lincomb.derived[48:94, "0.025quant"]

plot(reg_res_haz95l, switched_res_haz95l)
abline(0, 1, col = "darkgreen")

reg_res_haz95u <- regular_results$`Bivariate shared BYM`$summary.lincomb.derived[1:47, "0.975quant"]
switched_res_haz95u <- switched_results$`Bivariate shared BYM`$summary.lincomb.derived[48:94, "0.975quant"]

plot(reg_res_haz95u, switched_res_haz95u)
abline(0, 1, col = "darkgreen")

## waz

reg_res_waz <- regular_results$`Bivariate shared BYM`$summary.lincomb.derived[48:94, "mean"]
switched_res_waz <- switched_results$`Bivariate shared BYM`$summary.lincomb.derived[1:47, "mean"]

plot(reg_res_waz, switched_res_waz)
abline(0, 1, col = "darkgreen")

reg_res_waz95l <- regular_results$`Bivariate shared BYM`$summary.lincomb.derived[48:94, "0.025quant"]
switched_res_waz95l <- switched_results$`Bivariate shared BYM`$summary.lincomb.derived[1:47, "0.025quant"]

plot(reg_res_waz95l, switched_res_waz95l)
abline(0, 1, col = "darkgreen")

reg_res_waz95u <- regular_results$`Bivariate shared BYM`$summary.lincomb.derived[48:94, "0.975quant"]
switched_res_waz95u <- switched_results$`Bivariate shared BYM`$summary.lincomb.derived[1:47, "0.975quant"]

plot(reg_res_waz95u, switched_res_waz95u)
abline(0, 1, col = "darkgreen")

# compare REs

reg_convolved_haz <- regular_results$`Bivariate shared BYM`$summary.random$admin1.haz[1:47, "mean"]
reg_convolved_waz <- regular_results$`Bivariate shared BYM`$summary.random$admin1.waz[1:47, "mean"]
reg_hyperpar <- regular_results$`Bivariate shared BYM`$summary.hyperpar
reg_lambda <- reg_hyperpar[grepl("Beta", rownames(reg_hyperpar)), "mean"]

switched_convolved_haz <- switched_results$`Bivariate shared BYM`$summary.random$admin1.haz[1:47, "mean"]
switched_convolved_waz <- switched_results$`Bivariate shared BYM`$summary.random$admin1.waz[1:47, "mean"]
switched_hyperpar <- switched_results$`Bivariate shared BYM`$summary.hyperpar
switched_lambda <- switched_hyperpar[grepl("Beta", rownames(switched_hyperpar)), "mean"]

# compare
comp_waz <- (0.5*reg_convolved_waz) - ((1/(2*reg_lambda)) * reg_convolved_haz)
plot(comp_waz~switched_convolved_waz)
abline(0, 1, col = "darkgreen")

comp_haz <- reg_convolved_haz + (reg_lambda*reg_convolved_waz)
plot(comp_haz~switched_convolved_haz)
abline(0, 1, col = "darkgreen")
