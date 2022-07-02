# Austin Schumacher
# 1/26/2022
# Make comparison plots for INLA models 
#   (some or all to be included in the dissertation chapter)

# preamble ####

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
library(GGally)
library(cowplot)

savedir <- "~/Dropbox/dissertation_2/survey-csmf/results/proj-3-chapter-results/"

# load data and modeling results ####

## raw data 
load("../../../Dropbox/dissertation_2/survey-csmf/data/ken_dhs2014/data/haz-waz-kenDHS2014.rda")

## model results
fits <- read_rds("../../../Dropbox/dissertation_2/survey-csmf/results/ken2014-unit-level/hazwaz/ken2014-unit-level-hazwaz-inla-fits.rds")
posteriors <- read_rds("../../../Dropbox/dissertation_2/survey-csmf/results/ken2014-unit-level/hazwaz/ken2014-unit-level-hazwaz-inla-posteriors.rds")

# model comparisons ####
model_names <- names(posteriors)
n_regions <- nrow(admin1.mat)

# posterior SDs ####
post_sd <- lapply(posteriors,
                  function(x) {
                      lapply(x$samples_final,
                             function(y) {
                                 apply(y, 1, sd)
                             })
                  })

post_sd_res <- tibble(model= rep(model_names, each = n_regions),
                      region = rep(1:n_regions, length(model_names)),
                      haz = NA,
                      waz = NA)
for (i in 1:length(post_sd)) {
    post_sd_res$haz[post_sd_res$model == model_names[i]] <- post_sd[[model_names[i]]]$haz
    post_sd_res$waz[post_sd_res$model == model_names[i]] <- post_sd[[model_names[i]]]$waz
}

post_sd_res %<>% pivot_longer(c("haz", "waz"),
                              names_to = "outcome", 
                              values_to = "sd")

# make posterior sd plots
ggplot(post_sd_res, aes(x = sd, y = reorder(model, sd, na.rm = TRUE), fill = model)) + 
    geom_boxplot(show.legend = FALSE) + 
    scale_fill_viridis(option = "C", discrete = TRUE, alpha = 0.5) +
    xlab("posterior SDs") +
    ylab("Model") +
    theme_light()

ggsave(paste0(savedir, "posterior_sd_compare.pdf"), width = 7, height = 4.5)

# pairwise med comparison ####
post_med <- lapply(posteriors,
                   function(x) {
                       lapply(x$samples_final,
                              function(y) {
                                  apply(y, 1, median)
                              })
                   })

post_med_res <- tibble(model = rep(model_names, each = n_regions),
                       region = rep(1:n_regions, length(model_names)),
                       haz = NA,
                       waz = NA)

for (i in 1:length(post_med)) {
    post_med_res$haz[post_med_res$model == model_names[i]] <- post_med[[model_names[i]]]$haz
    post_med_res$waz[post_med_res$model == model_names[i]] <- post_med[[model_names[i]]]$waz
}

post_med_res %<>% pivot_wider(names_from = model,
                              values_from = c(haz, waz),
                              names_sep = " ")
# names(post_med_res)

# plot
p1 <- ggpairs(post_med_res %>% select(starts_with("haz")))
p2 <- ggpairs(post_med_res %>% select(starts_with("waz")))

plot_grid(
    ggmatrix_gtable(p1),
    ggmatrix_gtable(p2),
    ncol =2
)

ggsave(paste0(savedir, "med_estimate_compare-hazwaz.pdf"), 
       width = 10, height = 5)

# final table of estimates ####
fres <- fits$`BYM shared`

# extract fixed intercepts
fe <- fres$summary.fixed[, c(1, 3:5)]

# extract hyperparameters and lambda
hyper <- fres$summary.hyperpar[, c(1, 3:5)]

# convert precision to st dev
hyper[c(1:3, 5, 7:8),] <- (hyper[c(1:3, 5, 7:8),])^(-0.5)
hyper[c(1:3, 5, 7:8), 2:4] <- hyper[c(1:3, 5, 7:8), 4:2]

# make table
finres <- cbind(parameter =c("$\\beta_1$", "$\\beta_2$",
                             "$\\gamma_1$", "$\\gamma_2$",
                             "$\\omega_1$", "$\\omega_2$",
                             "$\\sigma_1$", "$\\rho_1$",
                             "$\\sigma_2$", "$\\rho_2$",
                             "$\\sigma_{\\epsilon1}$", "$\\sigma_{\\epsilon2}$",
                             "$\\lambda$"),
                rbind(fe, hyper))
rownames(finres) <- NULL
write_rds(finres, paste0(savedir, "final_table.rds"))

# maps

# maps of results ####
n_cats <- 15
plot.palette <- viridis(n_cats)

# plot final estimates ####
names(post_med_res) <- paste0(gsub(" ", "_", names(post_med_res)), "_med")

post_sd_res_2 <- post_sd_res %>% pivot_wider(names_from = c(outcome, model),
                                             values_from = c(sd),
                                             names_sep = "_")
names(post_sd_res_2) <- paste0(gsub(" ", "_", names(post_sd_res_2)), "_sd")

plotres <- left_join(post_med_res, post_sd_res_2, by = c("region_med" = "region_sd"))
plotres$region <- plotres$region_med

poly.adm1$region <- 1:nrow(poly.adm1)
tmp <- merge(poly.adm1, plotres,
             by = "region")
pdf(paste0(savedir, "shared-bym-haz-waz-pred-sd.pdf"), width = 9, height = 9)

grid.arrange(
    sp::spplot(obj = tmp, zcol = c("haz_BYM_shared_med", "waz_BYM_shared_med"),
               col.regions = plot.palette,
               cuts = length(plot.palette) - 1,
               layout = c(1, 2),
               names.attr = c("HAZ post. med", "WAZ post. med")),
    sp::spplot(obj = tmp, zcol = c("haz_BYM_shared_sd", "waz_BYM_shared_sd"),
               col.regions = plot.palette,
               cuts = length(plot.palette) - 1,
               layout = c(1, 2),
               names.attr = c("HAZ post. sd", "WAZ post. sd")),
    ncol = 2
)

dev.off()

# plot iid icar ####
# pairwise med comparison ####
post_meds <- lapply(posteriors,
                   function(x) {
                       lapply(x$samples_all,
                              function(y) {
                                  apply(y, 1, median)
                              })
                   })

post_RE_res <- tibble(region = c(1:n_regions, 1:n_regions),
                      haz_re = post_meds$`BYM shared`$haz_re,
                      waz_re = post_meds$`BYM shared`$waz_re,
                      shared_re = post_meds$`BYM shared`$waz_re)
haz_sigma <- posteriors$`BYM shared`$summaries$summary.hyperpar["Precision for admin1.haz", "0.5quant"]^(-0.5)
haz_phi <- posteriors$`BYM shared`$summaries$summary.hyperpar["Phi for admin1.haz", "0.5quant"]
waz_sigma <- posteriors$`BYM shared`$summaries$summary.hyperpar["Precision for admin1.waz", "0.5quant"]^(-0.5)
waz_phi <- posteriors$`BYM shared`$summaries$summary.hyperpar["Phi for admin1.waz", "0.5quant"]

haz_icar_standardized <- post_RE_res$haz_re[n_regions + (1:n_regions)]
haz_iid_standardized <- ((post_RE_res$haz_re[(1:n_regions)] / haz_sigma) - (sqrt(haz_phi) * haz_icar_standardized))/sqrt(1-haz_phi)
waz_icar_standardized <- post_RE_res$waz_re[n_regions + (1:n_regions)]
waz_iid_standardized <- ((post_RE_res$waz_re[(1:n_regions)] / waz_sigma) - (sqrt(waz_phi) * waz_icar_standardized))/sqrt(1-waz_phi)

RE_dat <- tibble(region = 1:n_regions,
                 haz_icar_standardized = haz_icar_standardized,
                 waz_icar_standardized = waz_icar_standardized,
                 haz_iid_standardized = haz_iid_standardized,
                 waz_iid_standardized = waz_iid_standardized,
                 haz_re_tot = post_RE_res$haz_re[(1:n_regions)],
                 waz_re_tot = post_RE_res$waz_re[(1:n_regions)])

tmp <- merge(poly.adm1, RE_dat,
             by = "region")
pdf(paste0(savedir, "shared-bym-haz-waz-iid-icar.pdf"), width = 12, height = 8)

grid.arrange(
    sp::spplot(obj = tmp, zcol = c("waz_re_tot", "haz_re_tot"),
               col.regions = plot.palette,
               cuts = length(plot.palette) - 1,
               layout = c(1, 2),
               names.attr = c("Total shared (WAZ)) BYM2","Total HAZ BYM2")),
    sp::spplot(obj = tmp, zcol = c("waz_iid_standardized", "waz_icar_standardized",
                                   "haz_iid_standardized", "haz_icar_standardized"),
               col.regions = plot.palette,
               cuts = length(plot.palette) - 1,
               layout = c(2, 2),
               names.attr = c("Shared (WAZ) IID", "Shared (WAZ) ICAR",
                              "HAZ IID", "HAZ ICAR")),
    ncol = 2
)

dev.off()
