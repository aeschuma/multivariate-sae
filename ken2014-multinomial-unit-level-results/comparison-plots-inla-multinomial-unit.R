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
fits <- read_rds("../../../Dropbox/dissertation_2/survey-csmf/results/ken2014-unit-level/multinomial/ken2014-unit-level-multinomial-inla-fits.rds")
posteriors <- read_rds("../../../Dropbox/dissertation_2/survey-csmf/results/ken2014-unit-level/multinomial/ken2014-unit-level-multinomial-inla-posteriors.rds")

model_names <- names(posteriors)
n_regions <- nrow(admin1.mat)

# what's our final model?
final_model <- "BYM shared"
# format posterior estimates ####

## SDs ####
post_sd <- lapply(posteriors,
                  function(x) {
                      lapply(x$samples_final,
                             function(y) {
                                 apply(y, 1, sd)
                             })
                  })

post_sd_res <- tibble(model= rep(model_names, each = n_regions),
                      region = rep(1:n_regions, length(model_names)),
                      none = NA,
                      modern = NA,
                      other = NA)
for (i in 1:length(post_sd)) {
    post_sd_res$none[post_sd_res$model == model_names[i]] <- post_sd[[model_names[i]]]$none
    post_sd_res$modern[post_sd_res$model == model_names[i]] <- post_sd[[model_names[i]]]$modern
    post_sd_res$other[post_sd_res$model == model_names[i]] <- post_sd[[model_names[i]]]$other
}

post_sd_res %<>% pivot_longer(c("none", "modern", "other"),
                              names_to = "outcome", 
                              values_to = "sd")

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
                            none = NA,
                            modern = NA,
                            other = NA)

for (i in 1:length(post_med)) {
    outs <- names(post_med[[i]])
    for (j in 1:length(outs)) {
        post_med_res_long[post_med_res_long$model == model_names[i], names(post_med[[i]])[j]] <- post_med[[model_names[i]]][[j]]
    }
}

post_med_res <- post_med_res_long %>% 
    pivot_wider(names_from = model,
                values_from = c(none, modern, other),
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
                           none = NA,
                           modern = NA,
                           other = NA)
for (i in 1:length(post_95width)) {
    post_95width_res$none[post_95width_res$model == model_names[i]] <- post_95width[[model_names[i]]]$none
    post_95width_res$modern[post_95width_res$model == model_names[i]] <- post_95width[[model_names[i]]]$modern
    post_95width_res$other[post_95width_res$model == model_names[i]] <- post_95width[[model_names[i]]]$other
}
post_95width_res %<>% pivot_longer(c("none", "modern", "other"),
                                   names_to = "outcome", 
                                   values_to = "width.95")

## RE medians ####
post_RE <- lapply(posteriors,
                  function(x) {
                      lapply(x$samples_all,
                             function(y) {
                                 apply(y, 1, median)
                             })
                  })

post_RE_res <- tibble(region = c(1:n_regions, 1:n_regions),
                      modern_re = post_RE[[final_model]]$modern_re,
                      none_re = post_RE[[final_model]]$none_re)
modern_sigma <- posteriors[[final_model]]$summaries$summary.hyperpar["Precision for admin1.modern", "0.5quant"]^(-0.5)
modern_phi <- posteriors[[final_model]]$summaries$summary.hyperpar["Phi for admin1.modern", "0.5quant"]
none_sigma <- posteriors[[final_model]]$summaries$summary.hyperpar["Precision for admin1.none", "0.5quant"]^(-0.5)
none_phi <- posteriors[[final_model]]$summaries$summary.hyperpar["Phi for admin1.none", "0.5quant"]

modern_icar_standardized <- post_RE_res$modern_re[n_regions + (1:n_regions)]
modern_iid_standardized <- ((post_RE_res$modern_re[(1:n_regions)] / modern_sigma) - (sqrt(modern_phi) * modern_icar_standardized))/sqrt(1-modern_phi)
none_icar_standardized <- post_RE_res$none_re[n_regions + (1:n_regions)]
none_iid_standardized <- ((post_RE_res$none_re[(1:n_regions)] / none_sigma) - (sqrt(none_phi) * none_icar_standardized))/sqrt(1-none_phi)

RE_dat <- tibble(region = 1:n_regions,
                 modern_icar_standardized = modern_icar_standardized,
                 none_icar_standardized = none_icar_standardized,
                 modern_iid_standardized = modern_iid_standardized,
                 none_iid_standardized = none_iid_standardized,
                 modern_re_tot = post_RE_res$modern_re[(1:n_regions)],
                 none_re_tot = post_RE_res$none_re[(1:n_regions)])

# make posterior sd plots ####
ggplot(post_sd_res, aes(x = sd, y = reorder(model, sd, na.rm = TRUE), fill = model)) + 
    geom_boxplot(show.legend = FALSE) + 
    scale_fill_viridis(option = "C", discrete = TRUE, alpha = 0.5) +
    xlab("posterior SDs") +
    ylab("Model") +
    theme_light()

ggsave(paste0(savedir, "posterior_sd_compare_multinomial-overall.pdf"), width = 5, height = 4)

ggplot(post_sd_res, aes(x = sd, y = reorder(model, sd, na.rm = TRUE), fill = model)) + 
    geom_boxplot(show.legend = FALSE) + 
    scale_fill_viridis(option = "C", discrete = TRUE, alpha = 0.5) +
    facet_wrap(~outcome, nrow = 1, scales = "free_x") +
    xlab("posterior SDs") +
    ylab("Model") +
    theme_light()
ggsave(paste0(savedir, "posterior_sd_compare_multinomial-byoutcome.pdf"), width = 8, height = 3)

# pairwise med comparison ####

# reshape for comparison plotting
post_med_res_comp <- post_med_res_long %>%
    pivot_longer(cols = c("none", "modern", "other"),
                 names_to = "outcome",
                 values_to = "value") %>%
    pivot_wider(names_from = "model",
                values_from = "value")

# plot
outcomes <- c("none", "modern", "other")
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
    ggsave(paste0(savedir, paste0("med_estimate_compare-multinomial-", oname,".pdf")), 
           width = 9, height = 4)
}

# final table of estimates ####
fres <- fits[[final_model]]

# extract fixed intercepts
fe <- fres$summary.fixed[1:4, c(3:7)]

# extract hyperparameters and lambda
hyper <- fres$summary.hyperpar[, c(3:7)]

# convert precision to st dev
hyper[c(1, 3),] <- (hyper[c(1, 3),])^(-0.5)
hyper[c(1, 3), ] <- hyper[c(1, 3), 5:1]

# make table
finres <- cbind(parameter =c("$\\beta_{\\text{none}}$", 
                             "$\\beta_{\\text{modern}}$", 
                             "$\\gamma_{\\text{none}}$", 
                             "$\\gamma_{\\text{modern}}$", 
                             "$\\sigma_{\\text{modern}}$", "$\\rho_{\\text{modern}}$",
                             "$\\sigma_{\\text{none}}$", "$\\rho_{\\text{none}}$",
                             "$\\lambda$"),
                rbind(fe, hyper))
rownames(finres) <- NULL
write_rds(finres, paste0(savedir, "final_table_multinomial.rds"))

# maps of results ####

# mapping set up
poly.adm1$region <- 1:nrow(poly.adm1)
n_cats <- 15
plot.palette <- viridis(n_cats)

# plot final estimates ####

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
pdf(paste0(savedir, "shared-bym-multinomial-pred-width95.pdf"), width = 9, height = 6)

grid.arrange(
    sp::spplot(obj = tmp, zcol = c("none_BYM_shared_med"),
               col.regions = plot.palette,
               cuts = length(plot.palette) - 1,
               layout = c(1, 1),
               main = "None post. med"),
    sp::spplot(obj = tmp, zcol = c("modern_BYM_shared_med"),
               col.regions = plot.palette,
               cuts = length(plot.palette) - 1,
               layout = c(1, 1),
               main = "Modern post. med"),
    sp::spplot(obj = tmp, zcol = c("other_BYM_shared_med"),
               col.regions = plot.palette,
               cuts = length(plot.palette) - 1,
               layout = c(1, 1),
               main = "Other post. med"),
    sp::spplot(obj = tmp, zcol = c("none_BYM_shared_width.95"),
               col.regions = plot.palette,
               cuts = length(plot.palette) - 1,
               layout = c(1, 1),
               main = "None 95% width"),
    sp::spplot(obj = tmp, zcol = c("modern_BYM_shared_width.95"),
               col.regions = plot.palette,
               cuts = length(plot.palette) - 1,
               layout = c(1, 1),
               main = "Modern 95% width"),
    sp::spplot(obj = tmp, zcol = c("other_BYM_shared_width.95"),
               col.regions = plot.palette,
               cuts = length(plot.palette) - 1,
               layout = c(1, 1),
               main = "Other 95% width"),
    ncol = 3, 
    nrow = 2
)

dev.off()

# plot iid icar ####

tmp <- merge(poly.adm1, RE_dat,
             by = "region")
pdf(paste0(savedir, "shared-bym-multinomial-iid-icar.pdf"), width = 12, height = 7)

grid.arrange(
    sp::spplot(obj = tmp, zcol = c("none_re_tot", "modern_re_tot"),
               col.regions = plot.palette,
               cuts = length(plot.palette) - 1,
               layout = c(1, 2),
               names.attr = c("None","Modern"),
               main = "Total BYM2"),
    sp::spplot(obj = tmp, zcol = c("none_iid_standardized", "none_icar_standardized",
                                   "modern_iid_standardized", "modern_icar_standardized"),
               col.regions = plot.palette,
               cuts = length(plot.palette) - 1,
               layout = c(2, 2),
               names.attr = c("None IID", "None ICAR",
                              "Modern IID", "Modern ICAR"),
               main = "IID and ICAR separate"),
    ncol = 2
)

dev.off()
