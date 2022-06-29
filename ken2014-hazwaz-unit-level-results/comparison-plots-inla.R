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

savedir <- "~/Dropbox/dissertation_2/survey-csmf/results/proj-2-chapter-results"

# functions ####

# comparison scatterplot function
comparison_scatter <- function(data, measure, x_var, y_var) {
    # testing
    # data = data
    # measure = measures[1]
    # x_var = comp.mod.2[1]
    # y_var = comp.mod.1[1]
    
    data.tmp <- data %>% select("admin1.name", "measure", "model.name", "value")
    
    # object name can't be same name as variable
    my.measure <- measure 
    
    # pivot and subset data for this plot
    tmp <- data.tmp %>% 
        filter(model.name %in% c(x_var,  y_var)) %>%
        pivot_wider(names_from = model.name, values_from = value) %>%
        filter(measure == my.measure)
    
    xyvars <- c(x_var, y_var)
    labs <- c(x_var, y_var)
    
    if (measure == "haz.icar.50") {
        main <- "Median HAZ ICAR RE"
    } else if (measure == "waz.icar.50") {
        main <- "Median WAZ ICAR RE"
    } else if (measure == "haz.iid.50") {
        main <- "Median HAZ IID RE"
    } else if (measure == "waz.iid.50") {
        main <- "Median WAZ IID RE"
    } else if (measure == "haz.pred.50") {
        main <- "Median predicted HAZ"
    } else if (measure == "waz.pred.50") {
        main <- "Median predicted WAZ"
    } else if (measure == "haz.pred.width80") {
        main <- "Width of HAZ posterior 80% interval"
    } else if (measure == "waz.pred.width80") {
        main <- "Width of WAZ posterior 80% interval"
    }
    
    # calculate mean abs diff
    my.x <- tmp %>% pull(x_var)
    my.y <- tmp %>% pull(y_var)
    mean.abs.diff <-  signif(mean(abs(my.x - my.y)), 3)
    
    # make plot
    ggplot(tmp %>% filter(measure == my.measure), aes(x = get(x_var), y = get(y_var))) +
        geom_point(alpha = 0.75, cex = 2) + 
        geom_abline(slope = 1, intercept = 0, col = "darkgreen", size = 1.2) + 
        xlab(labs[1]) +  
        ylab(labs[2]) +
        ggtitle(main,
                subtitle = paste0("Mean absolute diff: ", 
                                  mean.abs.diff)) +
        theme(axis.title = element_text(size=8)) +
        theme_light() 
}

# plot comparisons between 3 models x 2 outcomes (HAZ and WAZ)
twoway_comp_plot <- function(data, measure) {
    
    # testing
    # data <- all.results %>% filter(model.name %in% c("Direct", "Univariate IID", "Univariate BYM"))
    # measure <- "pred"
    
    if (!(measure %in% c("icar", "iid", "pred", "pred.width"))) stop("Measure not supported")
    
    measures <- c(NA, NA)
    if (measure == "icar") {
        measures <- c("haz.icar.50", "waz.icar.50")
    } else if (measure == "iid") {
        measures <- c("haz.iid.50", "waz.iid.50")
    } else if (measure == "pred") {
        measures <- c("haz.pred.50", "waz.pred.50")
    } else if (measure == "pred.width") {
        measures <- c("haz.pred.width95", "waz.pred.width95")
    } 
    
    # storage for list of plots
    models <- unique(data$model.name)
    nmodels <- length(models)
    ncomps <- choose(nmodels,2)
    npages <- ceiling(ncomps/12)
    med.comp.plots <- vector(mode = "list", length = 2*ncomps)
    comp.mod.1 <- c()
    comp.mod.2 <- c()
    for (i in 1:nmodels) {
        for (j in 1:i) {
            if (i == j) next
            comp.mod.1 <- c(comp.mod.1, models[i])
            comp.mod.2 <- c(comp.mod.2, models[j])
        }
    }
    
    ## start plots
    counter <- 0
    for (i in 1:length(measures)) {
      for (j in 1:ncomps) {
        counter <- counter + 1
        if (j == 1) {
          med.comp.plots[[counter]] <- comparison_scatter(data = data, 
                                                          measure = measures[i], 
                                                          x_var = comp.mod.2[j], 
                                                          y_var = comp.mod.1[j])
        } else {
          plt <- comparison_scatter(data = data, 
                                    measure = measures[i], 
                                    x_var = comp.mod.2[j], 
                                    y_var = comp.mod.1[j]) 
          med.comp.plots[[counter]] <- plt + ggtitle("")
        }
      }
    }
    
    # plot them all
    ggarrange(plotlist = med.comp.plots, ncol = ncomps, nrow = length(measures))
}

# load data and modeling results ####

## raw data 
load("../../../Dropbox/dissertation_2/survey-csmf/data/ken_dhs2014/data/haz-waz-kenDHS2014.rda")

## model results
fits <- read_rds("../../../Dropbox/dissertation_2/survey-csmf/results/ken2014-unit-level/hazwaz/ken2014-unit-level-hazwaz-inla-fits.rds")
posteriors <- read_rds("../../../Dropbox/dissertation_2/survey-csmf/results/ken2014-unit-level/hazwaz/ken2014-unit-level-hazwaz-inla-posteriors.rds")

# combine and format data ####

## direct estimates
allests <- results_direct %>% as_tibble() %>%
    mutate(model.name = "Direct",
           haz.pred.width95 = qnorm(0.975) * 2 * seHAZ.bi,
           waz.pred.width95 = qnorm(0.975) * 2 * seWAZ.bi) %>%
    select(model.name, admin1, admin1.name, meanHAZ.bi, meanWAZ.bi, haz.pred.width95, waz.pred.width95) %>%
    rename(haz.pred.50 = meanHAZ.bi,
           waz.pred.50 = meanWAZ.bi) %>%
    mutate(haz.iid.50 = NA,
           waz.iid.50 = NA,
           haz.icar.50 = NA,
           waz.icar.50 = NA)

for (i in 1:length(stage_2_list)) {
    tmpmodname <- names(stage_2_list)[i]
    mod.preds <- stage_2_list[[i]]$summary.lincomb.derived
    haz.pred <- mod.preds[grepl("haz", rownames(mod.preds)),]
    waz.pred <- mod.preds[grepl("waz", rownames(mod.preds)),]
    if (tmpmodname %in% c("Univariate IID", "Bivariate nonshared IID", "Bivariate shared IID")) {
        haz.re.iid <- stage_2_list[[i]]$summary.random$admin1.haz$`0.5quant`
        waz.re.iid <- stage_2_list[[i]]$summary.random$admin1.waz$`0.5quant`
        haz.re.icar <- rep(NA, n_regions)
        waz.re.icar <- rep(NA, n_regions)
    } else if (tmpmodname %in% c("Univariate BYM", "Bivariate nonshared BYM", "Bivariate shared BYM")) {
        u.haz <- stage_2_list[[i]]$summary.random$admin1.haz$`0.5quant`[(n_regions + 1):(n_regions*2)]
        phi.haz <- stage_2_list[[i]]$summary.hyperpar["Phi for admin1.haz", "0.5quant"]
        tau.haz <- stage_2_list[[i]]$summary.hyperpar["Precision for admin1.haz", "0.5quant"]
        x.haz <- stage_2_list[[i]]$summary.random$admin1.haz$`0.5quant`[1:n_regions]
        v.haz <- ((sqrt(tau.haz)*x.haz) - (sqrt(phi.haz)*u.haz))/(sqrt(1-phi.haz))
        
        u.waz <- stage_2_list[[i]]$summary.random$admin1.waz$`0.5quant`[(n_regions + 1):(n_regions*2)]
        phi.waz <- stage_2_list[[i]]$summary.hyperpar["Phi for admin1.waz", "0.5quant"]
        tau.waz <- stage_2_list[[i]]$summary.hyperpar["Precision for admin1.waz", "0.5quant"]
        x.waz <- stage_2_list[[i]]$summary.random$admin1.waz$`0.5quant`[1:n_regions]
        v.waz <- ((sqrt(tau.waz)*x.waz) - (sqrt(phi.waz)*u.waz))/(sqrt(1-phi.waz))
        
        haz.re.iid <- v.haz
        waz.re.iid <- v.waz
        haz.re.icar <- u.haz
        waz.re.icar <- u.waz
    }
    
    
    ## create a data.frame with results
    res.tmp <- tibble(model.name = rep(tmpmodname, n_regions),
                      admin1 = results_direct$admin1,
                      admin1.name = results_direct$admin1.name,
                      haz.iid.50 = haz.re.iid, 
                      waz.iid.50 = waz.re.iid,
                      haz.icar.50 = haz.re.icar, 
                      waz.icar.50 = waz.re.icar,
                      haz.pred.50 = haz.pred$`0.5quant`, 
                      haz.pred.width95 = haz.pred$`0.975quant` - haz.pred$`0.025quant`, 
                      waz.pred.50 = waz.pred$`0.5quant`, 
                      waz.pred.width95 = waz.pred$`0.975quant` - waz.pred$`0.025quant`)
    
    allests <- allests %>% bind_rows(res.tmp)
}

# start making plots ####
setwd(savedir)

# WAIC, DIC, CPO ####

modcomp <- tibble(model = names(stage_2_list),
                  # mlik = rep(NA, n_stage_2_models),
                  waic = rep(NA, n_stage_2_models),
                  dic = rep(NA, n_stage_2_models),
                  cpo = rep(NA, n_stage_2_models))

for (i in 1:n_stage_2_models) {
    tmp <- stage_2_list[[i]]
    # modcomp$mlik[i] <- as.character(round(tmp$mlik[1], digits = 4))
    modcomp$waic[i] <- as.character(round(tmp$waic$waic, digits = 4))
    modcomp$dic[i] <- as.character(round(tmp$dic$dic, digits = 4))
    modcomp$cpo[i] <- as.character(round(-1 * sum(log(tmp$cpo$cpo)), digits = 4))
}

comps <- names(modcomp)[names(modcomp) != "model"]
for (i in 1:length(comps)) {
    if (comps[i] == "mlik") {
        tobold <- as.character(max(modcomp %>% pull(comps[i]) %>% as.numeric()))
    } else {
        tobold <- as.character(min(modcomp %>% pull(comps[i]) %>% as.numeric()))
    }
    modcomp[which(modcomp[, comps[i]] == tobold), comps[i]] <- paste0("\\textbf{", tobold, "}")
}

## save table
write_rds(modcomp, "waic-dic-cpo-table.rds")

kable(modcomp, col.names = c("Model","WAIC", "DIC", "CPO"),
      format = "markdown", escape = FALSE)

# results storage for posterior SDs
model_name <- c("Direct", "Univariate IID", "Univariate BYM",  
           "Bivariate nonshared IID", "Bivariate shared IID",
           "Bivariate nonshared BYM", "Bivariate shared BYM")
outcome <- c("haz", "waz")
reg <- 1:47
posterior_sd_res <- expand_grid(model_name, outcome, reg)
posterior_sd_res$sd <- NA

# store direct estimates
posterior_sd_res$sd[posterior_sd_res$model_name == "Direct" & posterior_sd_res$outcome == "haz"] <- results_direct$seHAZ.bi
posterior_sd_res$sd[posterior_sd_res$model_name == "Direct" & posterior_sd_res$outcome == "waz"] <- results_direct$seWAZ.bi

# save posterior SDs of area-level estimates
for (i in 1:length(stage_2_list)) {
  mod.preds <- stage_2_list[[i]]$summary.lincomb.derived
  haz.pred <- mod.preds[grepl("haz", rownames(mod.preds)),]
  waz.pred <- mod.preds[grepl("waz", rownames(mod.preds)),]
  posterior_sd_res$sd[posterior_sd_res$model_name == names(stage_2_list)[i] & posterior_sd_res$outcome == "haz"] <- haz.pred$sd
  posterior_sd_res$sd[posterior_sd_res$model_name == names(stage_2_list)[i] & posterior_sd_res$outcome == "waz"] <- waz.pred$sd
}

# formatting all results ####
{
# par.names <- c("$\\beta_1$",
#                "$\\beta_2$",
#                "$\\sigma_1$",
#                "$\\sigma_2$",
#                "$\\rho_1$",
#                "$\\rho_2$",
#                "$\\lambda$")
# uni.par.ests.b11 <- as_tibble(rbind(mod.stan.summary.haz.uni.b11$summary["beta",],
#                                     mod.stan.summary.waz.uni.b11$summary["beta",],
#                                     mod.stan.summary.haz.uni.b11$summary["sigma",],
#                                     mod.stan.summary.waz.uni.b11$summary["sigma",],
#                                     mod.stan.summary.haz.uni.b11$summary["rho",],
#                                     mod.stan.summary.waz.uni.b11$summary["rho",],
#                                     rep(NA, ncol(mod.stan.summary.haz.uni.b11$summary))))
# uni.par.ests.b11$parameter <- par.names
# uni.par.ests.b11$model_name <- "univariate bym2"
# uni.par.ests.b11$model <- "univariate.bym2"
# 
# uni.par.ests.iid <- as_tibble(rbind(mod.stan.summary.haz.uni.iid$summary["beta",],
#                                     mod.stan.summary.waz.uni.iid$summary["beta",],
#                                     mod.stan.summary.haz.uni.iid$summary["sigma",],
#                                     mod.stan.summary.waz.uni.iid$summary["sigma",],
#                                     rep(NA, ncol(mod.stan.summary.haz.uni.iid$summary)),
#                                     rep(NA, ncol(mod.stan.summary.haz.uni.iid$summary)),
#                                     rep(NA, ncol(mod.stan.summary.haz.uni.iid$summary))))
# uni.par.ests.iid$parameter <- par.names
# uni.par.ests.iid$model_name <- "univariate iid"
# uni.par.ests.iid$model <- "univariate.iid"
# 
# nonshared.par.ests.b11 <- as_tibble(
#     rbind(mod.stan.summary.bi.nonshared.b11$summary[!grepl("preds", rownames(mod.stan.summary.bi.nonshared.b11$summary)),],
#           rep(NA, ncol(mod.stan.summary.bi.nonshared.b11$summary)))
# )
# nonshared.par.ests.b11$parameter <- par.names
# nonshared.par.ests.b11$model_name <- "bivariate nonshared bym2"
# nonshared.par.ests.b11$model <- "bivariate.nonshared.bym2"
# 
# nonshared.par.ests.iid <- as_tibble(
#     rbind(mod.stan.summary.bi.nonshared.iid$summary[!grepl("preds", rownames(mod.stan.summary.bi.nonshared.iid$summary)),],
#           rep(NA, ncol(mod.stan.summary.bi.nonshared.iid$summary)),
#           rep(NA, ncol(mod.stan.summary.bi.nonshared.iid$summary)),
#           rep(NA, ncol(mod.stan.summary.bi.nonshared.iid$summary)))
# )
# nonshared.par.ests.iid$parameter <- par.names
# nonshared.par.ests.iid$model_name <- "bivariate nonshared iid"
# nonshared.par.ests.iid$model <- "bivariate.nonshared.iid"
# 
# shared.par.ests.b11 <- as_tibble(mod.stan.summary.bi.shared.b11$summary[!grepl("preds", rownames(mod.stan.summary.bi.shared.b11$summary)),])
# shared.par.ests.b11$parameter <- par.names
# shared.par.ests.b11$model_name <- "bivariate shared"
# shared.par.ests.b11$model <- "bivariate.shared"
# 
# shared.par.ests.iid <- as_tibble(rbind(mod.stan.summary.bi.shared.iid$summary[!grepl("preds", rownames(mod.stan.summary.bi.shared.iid$summary)),],
#                                        rep(NA, ncol(mod.stan.summary.bi.shared.iid$summary)),
#                                        rep(NA, ncol(mod.stan.summary.bi.shared.iid$summary))))
# shared.par.ests.iid$parameter <- par.names
# shared.par.ests.iid$model_name <- "bivariate shared iid"
# shared.par.ests.iid$model <- "bivariate.shared.iid"
# 
# all.par.ests <- bind_rows(uni.par.ests.b11,
#                           uni.par.ests.iid,
#                           nonshared.par.ests.b11,
#                           nonshared.par.ests.iid,
#                           shared.par.ests.b11,
#                           shared.par.ests.iid) %>%
#     relocate(c(model, parameter))
}
all.results <- allests %>%
    pivot_longer(col = !c(model.name, admin1, admin1.name),
                 names_to = c("measure"),
                 values_to = "value") %>%
    mutate(model_factor = fct_relevel(as.factor(model.name),
                                      c("Direct", names(stage_2_list))))

posterior_sd_res %<>% mutate(model_factor = fct_relevel(as.factor(model_name),
                                                        c("Direct", names(stage_2_list))))

# make posterior sd plots ####
ggplot(posterior_sd_res %>% filter(outcome == "haz"), aes(x = sd, y = reorder(model_factor, sd, na.rm = TRUE), fill = model_factor)) + 
    geom_boxplot(show.legend = FALSE) + 
    scale_fill_viridis(option = "C", discrete = TRUE, alpha = 0.5) +
    xlab("posterior SDs") +
    ylab("Model") +
    theme_light()

ggsave("posterior_sd_compare.pdf", width = 7, height = 4.5)
    
# maps of results ####
n_cats <- 15
plot.palette <- viridis(n_cats)

# plot final estimates
all.results.wide <- all.results %>%
    select(admin1.name, model.name, measure, value) %>%
    pivot_wider(names_from = c(model.name, measure), values_from = value)
names(all.results.wide) <- gsub(" ", "_", names(all.results.wide))
outcomes <- c("Bivariate_shared_BYM_haz.iid.50", "Bivariate_shared_BYM_haz.icar.50",
              "Bivariate_shared_BYM_waz.iid.50", "Bivariate_shared_BYM_waz.icar.50")
outcomes.names <- c("HAZ IID", "HAZ ICAR",
                    "WAZ IID", "WAZ ICAR")

pdf("shared-bym-haz-waz-iid-icar-compare.pdf", width = 9, height = 9)
tmp <- merge(poly.adm1, all.results.wide[, c("admin1.name", outcomes)],
             by.x = "NAME_1", by.y = "admin1.name")
sp::spplot(obj = tmp, zcol = outcomes,
           col.regions = plot.palette,
           cuts = length(plot.palette) - 1,
           layout = c(2, 2),
           names.attr = outcomes.names)
dev.off()

pdf("shared-bym-haz-waz-pred.pdf", width = 9, height = 9)

outcomes <- c("Bivariate_shared_BYM_haz.pred.50", "Bivariate_shared_BYM_haz.pred.width95",
              "Bivariate_shared_BYM_waz.pred.50", "Bivariate_shared_BYM_waz.pred.width95")
outcomes.names <- c("HAZ posterior median", "HAZ 95% width",
                    "WAZ posterior median", "WAZ 95% width")

tmp <- merge(poly.adm1, all.results.wide[, c("admin1.name", outcomes[c(1, 3)])],
             by.x = "NAME_1", by.y = "admin1.name")
p1 <- sp::spplot(obj = tmp, zcol = outcomes[c(1, 3)],
           col.regions = plot.palette,
           cuts = length(plot.palette) - 1,
           layout = c(2, 1),
           names.attr = outcomes.names[c(1, 3)])

tmp <- merge(poly.adm1, all.results.wide[, c("admin1.name", outcomes[c(2, 4)])],
             by.x = "NAME_1", by.y = "admin1.name")
p2 <- sp::spplot(obj = tmp, zcol = outcomes[c(2, 4)],
           col.regions = plot.palette,
           cuts = length(plot.palette) - 1,
           layout = c(2, 1),
           names.attr = outcomes.names[c(2, 4)])
grid.arrange(p1, p2)
dev.off()

# make model comparison plots ####

# compare median estimate
twoway_comp_plot(all.results %>% filter(model.name %in% c("Direct", "Univariate IID", "Univariate BYM")), "pred")
ggsave("med_estimate_compare_univariate.pdf",
       width = 8.5, height = 3.6)

twoway_comp_plot(all.results %>% filter(model.name %in% c("Direct",  
                                                     "Bivariate nonshared IID", "Bivariate shared IID")), "pred")
ggsave("med_estimate_compare_bivariate_IID.pdf",
       width = 8.5, height = 3.6)

twoway_comp_plot(all.results %>% filter(model.name %in% c("Direct",
                                                          "Bivariate nonshared BYM", "Bivariate shared BYM")), "pred")
ggsave("med_estimate_compare_bivariate_bym2.pdf",
       width = 8.5, height = 3.6)

# IID vs ICAR plots ####
ggplot(all.results %>% 
           filter(measure %in% c("haz.iid.50", "waz.icar.50")) %>%
           filter(model.name %in% c("Bivariate nonshared BYM2", "Bivariate shared BYM2")) %>%
           pivot_wider(names_from = "measure", values_from = "value"),
       aes(x = waz.icar.50, y = haz.iid.50)) +
    geom_point() + 
    geom_smooth(method = "lm", formula = y ~ x) +
    facet_wrap(~ model_name, nrow = 1) +
    ylab("HAZ IID") +
    xlab("WAZ ICAR") +
    theme_light()
ggsave("haz_iid_waz_icar_compare_bivariate_nonshared_vs_shared.pdf")

# compare HAZ and WAZ preds ####
ggplot(all.results %>% 
           filter(measure %in% c("haz.pred.50", "waz.pred.50")) %>%
           pivot_wider(names_from = "measure", values_from = "value"),
       aes(x = waz.pred.50, y = haz.pred.50)) +
    geom_point() + 
    geom_smooth(method = "lm", formula = y ~ x) +
    facet_wrap(~ model, nrow = 1) +
    ylab("HAZ") +
    xlab("WAZ") +
    theme_light()
ggsave("haz_waz_med_estimate_compare.pdf")
