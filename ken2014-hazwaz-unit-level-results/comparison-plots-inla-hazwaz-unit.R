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
library(svyVGAM)
library(mvtnorm)
library(rgdal)
library(INLA)
library(viridis)
library(classInt)
library(gridExtra)
library(ggpubr)
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

model_names <- names(posteriors)
n_regions <- nrow(admin1.mat)

# what's our final model?
final_model <- "BYM shared"

# format posterior estimates ####

## SD ####
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
                              values_to = "sd") %>%
    mutate(outcome = ifelse(outcome == "haz", "HAZ", "WAZ"))

## final medians ####
post_med <- lapply(posteriors,
                   function(x) {
                       lapply(x$samples_final,
                              function(y) {
                                  apply(y, 1, median)
                              })
                   })

post_med_res_long <- tibble(model = rep(model_names, each = n_regions),
                            region = rep(1:n_regions, length(model_names)),
                            haz = NA,
                            waz = NA)

for (i in 1:length(post_med)) {
    outs <- names(post_med[[i]])
    for (j in 1:length(outs)) {
        post_med_res_long[post_med_res_long$model == model_names[i], names(post_med[[i]])[j]] <- post_med[[model_names[i]]][[j]]
    }
}

post_med_res <- post_med_res_long %>% 
    pivot_wider(names_from = model,
                values_from = c(haz, waz),
                names_sep = " ")

## CI widths ####
post_95width <- lapply(posteriors,
                  function(x) {
                      lapply(x$samples_final,
                             function(y) {
                                 apply(y, 1, function(z) { z[975] - z[25] })
                             })
                  })
post_95width_res <- tibble(model= rep(model_names, each = n_regions),
                           region = rep(1:n_regions, length(model_names)),
                           haz = NA,
                           waz = NA)
for (i in 1:length(post_95width)) {
    post_95width_res$haz[post_95width_res$model == model_names[i]] <- post_95width[[model_names[i]]]$haz
    post_95width_res$waz[post_95width_res$model == model_names[i]] <- post_95width[[model_names[i]]]$waz
}
post_95width_res %<>% pivot_longer(c("haz", "waz"),
                                   names_to = "outcome", 
                                   values_to = "width.95")

## random effect medians ####
post_RE <- lapply(posteriors,
                    function(x) {
                        lapply(x$samples_all,
                               function(y) {
                                   apply(y, 1, median)
                               })
                    })


post_RE_res <- tibble(region = c(1:n_regions, 1:n_regions),
                      haz_re = post_RE[[final_model]]$haz_re,
                      waz_re = post_RE[[final_model]]$waz_re,
                      shared_re = post_RE[[final_model]]$waz_re)
haz_sigma <- posteriors[[final_model]]$summaries$summary.hyperpar["Precision for admin1.haz", "0.5quant"]^(-0.5)
haz_phi <- posteriors[[final_model]]$summaries$summary.hyperpar["Phi for admin1.haz", "0.5quant"]
waz_sigma <- posteriors[[final_model]]$summaries$summary.hyperpar["Precision for admin1.waz", "0.5quant"]^(-0.5)
waz_phi <- posteriors[[final_model]]$summaries$summary.hyperpar["Phi for admin1.waz", "0.5quant"]

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

# posterior sd plots ####
ggplot(post_sd_res, aes(x = sd, y = reorder(model, sd, na.rm = TRUE), fill = model)) + 
    geom_boxplot(show.legend = FALSE) + 
    facet_wrap(~ outcome, nrow = 1, scales = "free_x") +
    scale_fill_viridis(option = "C", discrete = TRUE, alpha = 0.5) +
    xlab("posterior SDs") +
    ylab("Model") +
    theme_light()

ggsave(paste0(savedir, "posterior_sd_compare-hazwaz.pdf"), width = 7, height = 4)

# pairwise med comparison ####
ggpairs(post_med_res %>% select(starts_with("haz")))

# reshape for comparison plotting
post_med_res_comp <- post_med_res_long %>%
    pivot_longer(cols = c("haz", "waz"),
                 names_to = "outcome",
                 values_to = "value") %>%
    pivot_wider(names_from = "model",
                values_from = "value")

# plot
outcomes <- c("haz", "waz")
for (oo in 1:length(outcomes)) {
    oname <- outcomes[oo]
    tmp <- post_med_res_comp %>% filter(outcome == oname)
    
    mean_abs_diff <- tmp %>% 
        mutate(mad = abs(`IID nonshared` - `BYM nonshared`)) %>%
        pull(mad) %>%
        mean()
    p1 <- ggplot(tmp, aes(x = `IID nonshared`,
                          y = `BYM nonshared`)) +
        geom_abline(slope = 1, intercept = 0, col = "black", alpha = 0.7) +
        geom_point(alpha = 0.5, col = "darkgreen") +
        # xlab(model_names[i]) +
        # ylab(model_names[j]) +
        ggtitle(paste0("Mean absolute diff: ", signif(mean_abs_diff, 2))) +
        theme_light()
    
    mean_abs_diff <- tmp %>% 
        mutate(mad = abs(`IID nonshared` - `IID shared`)) %>%
        pull(mad) %>%
        mean()
    p2 <- ggplot(tmp, aes(x = `IID nonshared`,
                          y = `IID shared`)) +
        geom_abline(slope = 1, intercept = 0, col = "black", alpha = 0.7) +
        geom_point(alpha = 0.5, col = "darkgreen") +
        # xlab(model_names[i]) +
        # ylab(model_names[j]) +
        ggtitle(paste0("Mean absolute diff: ", signif(mean_abs_diff, 2))) +
        theme_light()
    
    mean_abs_diff <- tmp %>% 
        mutate(mad = abs(`IID nonshared` - `BYM shared`)) %>%
        pull(mad) %>%
        mean()
    p3 <- ggplot(tmp, aes(x = `IID nonshared`,
                          y = `BYM shared`)) +
        geom_abline(slope = 1, intercept = 0, col = "black", alpha = 0.7) +
        geom_point(alpha = 0.5, col = "darkgreen") +
        # xlab(model_names[i]) +
        # ylab(model_names[j]) +
        ggtitle(paste0("Mean absolute diff: ", signif(mean_abs_diff, 2))) +
        theme_light()
    
    mean_abs_diff <- tmp %>% 
        mutate(mad = abs(`BYM nonshared` - `IID shared`)) %>%
        pull(mad) %>%
        mean()
    p4 <- ggplot(tmp, aes(x = `BYM nonshared`,
                          y = `IID shared`)) +
        geom_abline(slope = 1, intercept = 0, col = "black", alpha = 0.7) +
        geom_point(alpha = 0.5, col = "darkgreen") +
        # xlab(model_names[i]) +
        # ylab(model_names[j]) +
        ggtitle(paste0("Mean absolute diff: ", signif(mean_abs_diff, 2))) +
        theme_light()
    
    mean_abs_diff <- tmp %>% 
        mutate(mad = abs(`BYM nonshared` - `BYM shared`)) %>%
        pull(mad) %>%
        mean()
    p5 <- ggplot(tmp, aes(x = `BYM nonshared`,
                          y = `BYM shared`)) +
        geom_abline(slope = 1, intercept = 0, col = "black", alpha = 0.7) +
        geom_point(alpha = 0.5, col = "darkgreen") +
        # xlab(model_names[i]) +
        # ylab(model_names[j]) +
        ggtitle(paste0("Mean absolute diff: ", signif(mean_abs_diff, 2))) +
        theme_light()
    
    mean_abs_diff <- tmp %>% 
        mutate(mad = abs(`IID shared` - `BYM shared`)) %>%
        pull(mad) %>%
        mean()
    p6 <- ggplot(tmp, aes(x = `IID shared`,
                          y = `BYM shared`)) +
        geom_abline(slope = 1, intercept = 0, col = "black", alpha = 0.7) +
        geom_point(alpha = 0.5, col = "darkgreen") +
        # xlab(model_names[i]) +
        # ylab(model_names[j]) +
        ggtitle(paste0("Mean absolute diff: ", signif(mean_abs_diff, 2))) +
        theme_light()
    
    ggpubr::ggarrange(p1, p2, p3,  p4, p5, p6, ncol = 3, nrow = 2)
    ggsave(paste0(savedir, paste0("med_estimate_compare-hazwaz-", oname,".pdf")), 
           width = 9, height = 4)
}

# final table of estimates ####
fres <- fits[[final_model]]

# extract fixed intercepts
fe <- fres$summary.fixed[, c(3:7)]

# extract hyperparameters and lambda
hyper <- fres$summary.hyperpar[, c(3:7)]

# convert precision to st dev
hyper[c(1:3, 5, 7:8),] <- (hyper[c(1:3, 5, 7:8),])^(-0.5)
hyper[c(1:3, 5, 7:8), ] <- hyper[c(1:3, 5, 7:8), 5:1]

# make table
finres <- cbind(parameter =c("$\\beta_{\\text{haz}}$", 
                             "$\\beta_{\\text{waz}}$", 
                             "$\\gamma_{\\text{haz}}$", 
                             "$\\gamma_{\\text{waz}}$",
                             "$\\omega_1$", "$\\omega_1$",
                             "$\\sigma_1$", "$\\rho_1$",
                             "$\\sigma_2$", "$\\rho_2$",
                             "$\\sigma_{\\epsilon 1}$",
                             "$\\sigma_{\\epsilon 2}$",
                             "$\\lambda$"),
                rbind(fe, hyper))
rownames(finres) <- NULL
write_rds(finres, paste0(savedir, "final_table_hazwaz.rds"))

# maps

# maps of results ####
n_cats <- 15
plot.palette <- viridis(n_cats)
poly.adm1$region <- 1:nrow(poly.adm1)

## final estimates ####

# format and merge data
names(post_med_res) <- paste0(gsub(" ", "_", names(post_med_res)), "_med")

post_95width_res_2 <- post_95width_res %>% 
    pivot_wider(names_from = c(outcome, model),
                values_from = c(width.95),
                names_sep = "_")
names(post_95width_res_2) <- paste0(gsub(" ", "_", names(post_95width_res_2)), "_width.95")

plotres <- left_join(post_med_res, post_95width_res_2, by = c("region_med" = "region_width.95"))
plotres$region <- plotres$region_med

# plotting
tmp <- merge(poly.adm1, plotres,
             by = "region")
pdf(paste0(savedir, "shared-bym-haz-waz-pred-width95.pdf"), width = 9, height = 9)

grid.arrange(
    sp::spplot(obj = tmp, zcol = c("haz_BYM_shared_med", "waz_BYM_shared_med"),
               col.regions = plot.palette,
               cuts = length(plot.palette) - 1,
               layout = c(2, 1),
               names.attr = c("HAZ", "WAZ"),
               main = "Posterior medians"),
    sp::spplot(obj = tmp, zcol = c("haz_BYM_shared_width.95", "waz_BYM_shared_width.95"),
               col.regions = plot.palette,
               cuts = length(plot.palette) - 1,
               layout = c(2, 1),
               names.attr = c("HAZ", "WAZ"),
               main = "95% CI widths"),
    ncol = 1
)

dev.off()

## random effects ####
tmp <- merge(poly.adm1, RE_dat,
             by = "region")
pdf(paste0(savedir, "shared-bym-haz-waz-iid-icar.pdf"), width = 12, height = 7)

grid.arrange(
    sp::spplot(obj = tmp, zcol = c("waz_re_tot", "haz_re_tot"),
               col.regions = plot.palette,
               cuts = length(plot.palette) - 1,
               layout = c(1, 2),
               names.attr = c("Shared (WAZ)","HAZ"),
               main = "Total BYM2"),
    sp::spplot(obj = tmp, zcol = c("waz_iid_standardized", "waz_icar_standardized",
                                   "haz_iid_standardized", "haz_icar_standardized"),
               col.regions = plot.palette,
               cuts = length(plot.palette) - 1,
               layout = c(2, 2),
               names.attr = c("Shared (WAZ) IID", "Shared (WAZ) ICAR",
                              "HAZ IID", "HAZ ICAR"),
               main = "IID and ICAR separate"),
    ncol = 2
)

dev.off()
