# Austin Schumacher
# 1/26/2022
# Make comparison plots for INLA models 
#   (some or all to be included in the dissertation chapter)

# preamble ####

library(SUMMER);library(tidyverse);library(spdep);library(geosphere);library(haven);
library(knitr);library(kableExtra);library(magrittr);library(svyVGAM);library(mvtnorm);
library(rgdal);library(INLA);library(viridis);library(classInt);library(gridExtra);
library(ggpubr);library(grid);

savedir <- "~/Dropbox/dissertation_2/survey-csmf/results/proj-3-chapter-results/"

# load data and modeling results ####

## raw data 
load("../../../Dropbox/dissertation_2/survey-csmf/data/ken_dhs2014/data/haz-waz-kenDHS2014.rda")
n_regions <- length(unique(dat$admin1))

## area level estimates (from inla)
area_list <- read_rds("../../../Dropbox/dissertation_2/survey-csmf/results/ken2014-hazwaz/ken2014-hazwaz-stage-2-inla-all.rds")
area_list[[7]] <- NULL

## unit level results
unit_fits <- read_rds("../../../Dropbox/dissertation_2/survey-csmf/results/ken2014-unit-level/hazwaz/ken2014-unit-level-hazwaz-inla-fits.rds")
unit_posteriors <- read_rds("../../../Dropbox/dissertation_2/survey-csmf/results/ken2014-unit-level/hazwaz/ken2014-unit-level-hazwaz-inla-posteriors.rds")

## compile
area_mod <- area_list[["Bivariate shared BYM"]]
uni_mod <- unit_fits[["BYM shared"]]
unit_post <- unit_posteriors[["BYM shared"]]
area_list <- NULL
unit_fits <- NULL
unit_posteriors <- NULL

area_post_names <- rownames(area_mod$summary.lincomb.derived)
area_res <- tibble(model = "arealevel",
                   admin1 = rep(1:n_regions, 2),
                   outcome = c(rep("HAZ", n_regions), rep("WAZ", n_regions)),
                   pred = area_mod$summary.lincomb.derived[grepl("reg", area_post_names), "0.5quant"],
                   lower80 = area_mod$summary.lincomb.derived[grepl("reg", area_post_names), "0.1quant"],
                   upper80 = area_mod$summary.lincomb.derived[grepl("reg", area_post_names), "0.9quant"],
                   width80 = area_mod$summary.lincomb.derived[grepl("reg", area_post_names), "0.9quant"] -
                       area_mod$summary.lincomb.derived[grepl("reg", area_post_names), "0.1quant"])
   
unit_res_list <- lapply(unit_post$samples_final,
                       function(y) {
                           tibble(model = "unitlevel",
                                  admin1 = 1:n_regions,
                                  pred = apply(y, 1, median),
                                  lower80 = apply(y, 1, quantile, 0.1),
                                  upper80 = apply(y, 1, quantile, 0.9),
                                  width80 = apply(y, 1, function(x) {quantile(x, 0.9) - quantile(x, 0.1)}))
                       })
unit_res <- bind_rows(unit_res_list$haz %>% mutate(outcome = "HAZ"),
                      unit_res_list$waz %>% mutate(outcome = "WAZ"))

results <- bind_rows(area_res, unit_res) %>% 
    left_join(dat %>% dplyr::select(admin1, admin1.name) %>% distinct())
 
## maps of results ####
n_cats <- 15
plot.palette <- viridis(n_cats)

# format for making maps
all.results.wide <- results %>%
    pivot_wider(names_from = c("model", "outcome"), 
                values_from = c("pred", "lower80", "upper80", "width80"))

tmp <- merge(poly.adm1, all.results.wide,
             by.x = "NAME_1", by.y = "admin1.name")

### comparison of univariate BYM HAZ and WAZ models ####

pdf(paste0(savedir, "area-unit-haz-waz-maps-med-compare.pdf"), width = 9, height = 3)
ggarrange(
    sp::spplot(obj = tmp, zcol = c("pred_arealevel_HAZ", "pred_unitlevel_HAZ"),
               col.regions = plot.palette,
               cuts = length(plot.palette) - 1,
               layout = c(2, 1),
               names.attr = c("HAZ Area-level", "HAZ Unit-level")),
    sp::spplot(obj = tmp, zcol = c("pred_arealevel_WAZ", "pred_unitlevel_WAZ"),
               col.regions = plot.palette,
               cuts = length(plot.palette) - 1,
               layout = c(2, 1),
               names.attr = c("WAZ Area-level", " WAZ Unit-level")),
    ncol = 2, nrow = 1
)
dev.off()

pdf(paste0(savedir, "area-unit-haz-waz-maps-width-compare.pdf"), width = 9, height = 3)

ggarrange(
    sp::spplot(obj = tmp, zcol = c("width80_arealevel_HAZ", "width80_unitlevel_HAZ"),
               col.regions = plot.palette,
               cuts = length(plot.palette) - 1,
               layout = c(2, 1),
               names.attr = c("HAZ Area-level", "HAZ Unit-level")),
    sp::spplot(obj = tmp, zcol = c("width80_arealevel_WAZ", "width80_unitlevel_WAZ"),
               col.regions = plot.palette,
               cuts = length(plot.palette) - 1,
               layout = c(2, 1),
               names.attr = c("WAZ Area-level", "WAZ Unit-level")),
    ncol = 2, nrow = 1
)
dev.off()

## model comparison scatters ####
ggplot(results %>% pivot_wider(names_from = "model",
                               values_from = c("pred", "lower80", "upper80", "width80")),
       aes(x = pred_arealevel, y = pred_unitlevel, 
           xmin = lower80_arealevel, xmax = upper80_arealevel,
           ymin = lower80_unitlevel, ymax = upper80_unitlevel)) +
    geom_abline(slope = 1, intercept = 0, alpha = 0.5, 
                col = "darkgreen", lwd = 1.25) +
    geom_point(alpha = 0.7) +
    geom_errorbar(alpha = 0.7) +
    geom_errorbarh(alpha = 0.7) +
    facet_wrap(~ outcome) +
    xlab("Area-level") +
    ylab("Unit-level") +
    theme_light()
ggsave(paste0(savedir, "area-unit-haz-waz-scatter-compare.pdf"), width = 7, height = 3.5)

# model comparison line plot
results %>% mutate(Model = ifelse(model == "arealevel", "Area level", "Unit level"),
                   admin1.name = fct_reorder(admin1.name, width80, .desc = TRUE)) %>%
    ggplot(aes(x = admin1.name, y = pred, 
               ymin = lower80, ymax = upper80,
               col = Model)) +
    geom_point(aes(pch = Model), position = position_dodge(0.6)) +
    geom_errorbar(aes(lty = Model), position = position_dodge(0.6), width = 0.4) +
    facet_wrap(~ outcome, nrow = 2) +
    scale_color_viridis(option = "D", discrete = TRUE,
                        end = 0.7) +
    ylab("Posterior ests") +
    xlab("County") +
    theme_light() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 6),
          legend.position = "top", 
          legend.margin = margin(t=-0.2, r=0, b=-0.3, l=0, unit="cm"))
ggsave(paste0(savedir, "area-unit-haz-waz-region-lines-compare.pdf"), width = 10, height = 4.5)

