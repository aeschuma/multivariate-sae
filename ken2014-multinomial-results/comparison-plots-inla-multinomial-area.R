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
    
    if (measure == "beta1.icar.50") {
        main <- "Median beta1 ICAR RE"
    } else if (measure == "beta2.icar.50") {
        main <- "Median beta2 ICAR RE"
    } else if (measure == "beta1.iid.50") {
        main <- "Median beta1 IID RE"
    } else if (measure == "beta2.iid.50") {
        main <- "Median beta2 IID RE"
    } else if (measure == "beta1.pred.50") {
        main <- "Med logit(p) modern vs. none"
    } else if (measure == "beta2.pred.50") {
        main <- "Med logit(p) other vs. none"
    } else if (measure == "beta1.pred.width80") {
        main <- "Width of beta1 posterior 80% interval"
    } else if (measure == "beta2.pred.width80") {
        main <- "Width of beta2 posterior 80% interval"
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

# plot comparisons between 3 models x 2 outcomes (beta1 and beta2)
twoway_comp_plot <- function(data, measure) {
    
    # testing
    # data <- all.results %>% filter(model.name %in% c("Direct", "Univariate IID", "Univariate BYM"))
    # measure <- "pred"
    
    if (!(measure %in% c("icar", "iid", "pred", "pred.width"))) stop("Measure not supported")
    
    measures <- c(NA, NA)
    if (measure == "icar") {
        measures <- c("beta1.icar.50", "beta2.icar.50")
    } else if (measure == "iid") {
        measures <- c("beta1.iid.50", "beta2.iid.50")
    } else if (measure == "pred") {
        measures <- c("beta1.pred.50", "beta2.pred.50")
    } else if (measure == "pred.width") {
        measures <- c("beta1.pred.width95", "beta2.pred.width95")
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

## direct estimates (stage 1)
stage_1_list <- read_rds("../../../Dropbox/dissertation_2/survey-csmf/results/ken2014-multinomial/ken2014-multinomial-stage-1.rds")
results_direct <- stage_1_list$results
n_regions <- nrow(results_direct)

## stage 2 estimates (from inla)
stage_2_list <- read_rds("../../../Dropbox/dissertation_2/survey-csmf/results/ken2014-multinomial/ken2014-multinomial-stage-2-inla-all.rds")
stage_2_list[[7]] <- NULL
n_stage_2_models <- length(stage_2_list)

# combine and format data ####

## direct estimates
allests <- results_direct %>% as_tibble() %>%
    mutate(model.name = "Direct",
           beta1.pred.width95 = qnorm(0.975) * 2 * se_beta1,
           beta2.pred.width95 = qnorm(0.975) * 2 * se_beta2) %>%
    select(model.name, admin1, admin1.name, 
           mean_beta1, mean_beta2, 
           beta1.pred.width95, beta2.pred.width95) %>%
    rename(beta1.pred.50 = mean_beta1,
           beta2.pred.50 = mean_beta2) %>%
    mutate(beta1.iid.50 = NA,
           beta2.iid.50 = NA,
           beta1.icar.50 = NA,
           beta2.icar.50 = NA,
           beta1.totre.50 = NA,
           beta2.totre.50 = NA)

for (i in 1:length(stage_2_list)) {
    tmpmodname <- names(stage_2_list)[i]
    mod.preds <- stage_2_list[[i]]$summary.lincomb.derived
    beta1.pred <- mod.preds[grepl("outcome1", rownames(mod.preds)),]
    beta2.pred <- mod.preds[grepl("oucome2", rownames(mod.preds)),]
    if (tmpmodname %in% c("Univariate IID", "Bivariate nonshared IID", "Bivariate shared IID")) {
        beta1.re.iid <- stage_2_list[[i]]$summary.random$admin1.beta1$`0.5quant`
        beta2.re.iid <- stage_2_list[[i]]$summary.random$admin1.beta2$`0.5quant`
        beta1.re.icar <- rep(NA, n_regions)
        beta2.re.icar <- rep(NA, n_regions)
        beta1.re.tot <- rep(NA, n_regions)
        beta2.re.tot <- rep(NA, n_regions)
    } else if (tmpmodname %in% c("Univariate BYM", "Bivariate nonshared BYM", "Bivariate shared BYM")) {
        u.beta1 <- stage_2_list[[i]]$summary.random$admin1.beta1$`0.5quant`[(n_regions + 1):(n_regions*2)]
        phi.beta1 <- stage_2_list[[i]]$summary.hyperpar["Phi for admin1.beta1", "0.5quant"]
        tau.beta1 <- stage_2_list[[i]]$summary.hyperpar["Precision for admin1.beta1", "0.5quant"]
        x.beta1 <- stage_2_list[[i]]$summary.random$admin1.beta1$`0.5quant`[1:n_regions]
        v.beta1 <- ((sqrt(tau.beta1)*x.beta1) - (sqrt(phi.beta1)*u.beta1))/(sqrt(1-phi.beta1))
        
        u.beta2 <- stage_2_list[[i]]$summary.random$admin1.beta2$`0.5quant`[(n_regions + 1):(n_regions*2)]
        phi.beta2 <- stage_2_list[[i]]$summary.hyperpar["Phi for admin1.beta2", "0.5quant"]
        tau.beta2 <- stage_2_list[[i]]$summary.hyperpar["Precision for admin1.beta2", "0.5quant"]
        x.beta2 <- stage_2_list[[i]]$summary.random$admin1.beta2$`0.5quant`[1:n_regions]
        v.beta2 <- ((sqrt(tau.beta2)*x.beta2) - (sqrt(phi.beta2)*u.beta2))/(sqrt(1-phi.beta2))
        
        beta1.re.iid <- v.beta1
        beta2.re.iid <- v.beta2
        beta1.re.icar <- u.beta1
        beta2.re.icar <- u.beta2
        
        beta1.totre <- stage_2_list[[i]]$summary.random$admin1.beta1$`0.5quant`[1:n_regions]
        beta2.totre <- stage_2_list[[i]]$summary.random$admin1.beta2$`0.5quant`[1:n_regions]
        beta1.re.tot <- beta1.totre
        beta2.re.tot <- beta2.totre
    }
    
    
    ## create a data.frame with results
    res.tmp <- tibble(model.name = rep(tmpmodname, n_regions),
                      admin1 = results_direct$admin1,
                      admin1.name = results_direct$admin1.name,
                      beta1.iid.50 = beta1.re.iid, 
                      beta2.iid.50 = beta2.re.iid,
                      beta1.icar.50 = beta1.re.icar, 
                      beta2.icar.50 = beta2.re.icar,
                      beta1.totre.50 = beta1.re.tot, 
                      beta2.totre.50 = beta2.re.tot,
                      beta1.pred.50 = beta1.pred$`0.5quant`, 
                      beta1.pred.width95 = beta1.pred$`0.975quant` - beta1.pred$`0.025quant`, 
                      beta2.pred.50 = beta2.pred$`0.5quant`, 
                      beta2.pred.width95 = beta2.pred$`0.975quant` - beta2.pred$`0.025quant`)
    
    allests <- allests %>% bind_rows(res.tmp)
}

# probability data for only final model
fmod <- stage_2_list$`Bivariate shared BYM`
nsamps <- 1000
samp <- inla.posterior.sample(nsamps, fmod)

# process fitted values
modern_fe_idx <- which(rownames(samp[[1]]$latent) %in% c("outcomebeta1:1"))
modern_fe_mat <- matrix(NA, nrow = length(modern_fe_idx), ncol = nsamps)
other_fe_idx <- which(rownames(samp[[1]]$latent) %in% c("outcomebeta2:1"))
other_fe_mat <- matrix(NA, nrow = length(other_fe_idx), ncol = nsamps)

modern_re_idx <- rownames(samp[[1]]$latent) %>% str_detect("admin1\\.beta1:") %>% which()
other_re_idx <- rownames(samp[[1]]$latent) %>% str_detect("admin1\\.beta2:") %>% which()
modern_re_mat <- matrix(NA, nrow = length(modern_re_idx), ncol = nsamps)
other_re_mat <- matrix(NA, nrow = length(other_re_idx), ncol = nsamps)

shared_re_idx <- rownames(samp[[1]]$latent) %>% str_detect("admin1\\.beta1\\.2:") %>% which()
shared_re_mat <- matrix(NA, nrow = length(shared_re_idx), ncol = nsamps)

# fill in sample matrices
for (s in 1:nsamps) {
    modern_fe_mat[,s] <- samp[[s]]$latent[modern_fe_idx]
    other_fe_mat[,s] <- samp[[s]]$latent[other_fe_idx]
    modern_re_mat[,s] <- samp[[s]]$latent[modern_re_idx]
    other_re_mat[,s] <- samp[[s]]$latent[other_re_idx]
    shared_re_mat[,s] <- samp[[s]]$latent[shared_re_idx]
}

# obtain total effect for space - first half of bym2, or all of them for IID
modern_re_tot_mat <- modern_re_mat[1:n_regions,] + shared_re_mat[1:n_regions,]
other_re_tot_mat <- other_re_mat[1:n_regions,] 

# matrix of fitted estimates
fitted_modern_mat <- modern_fe_mat[rep(1,n_regions),] + modern_re_tot_mat
fitted_other_mat <- other_fe_mat[rep(1,n_regions),] + other_re_tot_mat

# transform to probs
fitted_probs <- array(NA, dim = c(3, n_regions, nsamps))
for (s in 1:nsamps) {
    for (r in 1:n_regions) {
        denom <- 1 + exp(fitted_modern_mat[r, s]) + exp(fitted_other_mat[r, s])
        fitted_probs[1, r, s] <- 1/denom
        fitted_probs[2, r, s] <- exp(fitted_modern_mat[r, s])/denom
        fitted_probs[3, r, s] <- exp(fitted_other_mat[r, s])/denom
    }
}

# summaries of probs
posterior_probs_res <- tibble(region = 1:n_regions,
                              admin1.name = results_direct$admin1.name,
                              prob_none_med = apply(fitted_probs[1,,], 1, median),
                              prob_modern_med = apply(fitted_probs[2,,], 1, median),
                              prob_other_med = apply(fitted_probs[3,,], 1, median),
                              prob_none_width95 = apply(fitted_probs[1,,], 1, function(x) { sort(x)[975] - sort(x)[25] }),
                              prob_modern_width95 = apply(fitted_probs[2,,], 1, function(x) { sort(x)[975] - sort(x)[25] }),
                              prob_other_width95 = apply(fitted_probs[3,,], 1, function(x) { sort(x)[975] - sort(x)[25] }))

# start making plots ####
setwd(savedir)

# WAIC, DIC, CPO ####

# modcomp <- tibble(model = names(stage_2_list),
#                   # mlik = rep(NA, n_stage_2_models),
#                   waic = rep(NA, n_stage_2_models),
#                   dic = rep(NA, n_stage_2_models),
#                   cpo = rep(NA, n_stage_2_models))
# 
# for (i in 1:n_stage_2_models) {
#     tmp <- stage_2_list[[i]]
#     # modcomp$mlik[i] <- as.character(round(tmp$mlik[1], digits = 4))
#     modcomp$waic[i] <- as.character(round(tmp$waic$waic, digits = 4))
#     modcomp$dic[i] <- as.character(round(tmp$dic$dic, digits = 4))
#     modcomp$cpo[i] <- as.character(round(-1 * sum(log(tmp$cpo$cpo)), digits = 4))
# }
# 
# comps <- names(modcomp)[names(modcomp) != "model"]
# for (i in 1:length(comps)) {
#     if (comps[i] == "mlik") {
#         tobold <- as.character(max(modcomp %>% pull(comps[i]) %>% as.numeric()))
#     } else {
#         tobold <- as.character(min(modcomp %>% pull(comps[i]) %>% as.numeric()))
#     }
#     modcomp[which(modcomp[, comps[i]] == tobold), comps[i]] <- paste0("\\textbf{", tobold, "}")
# }
# 
# ## save table
# write_rds(modcomp, "waic-dic-cpo-table.rds")
# 
# kable(modcomp, col.names = c("Model","WAIC", "DIC", "CPO"),
#       format = "markdown", escape = FALSE)

# results storage for posterior SDs
model_name <- c("Direct", "Univariate IID", "Univariate BYM",  
           "Bivariate nonshared IID", "Bivariate shared IID",
           "Bivariate nonshared BYM", "Bivariate shared BYM")
outcome <- c("beta1", "beta2")
reg <- 1:47
posterior_sd_res <- expand_grid(model_name, outcome, reg)
posterior_sd_res$sd <- NA

# store direct estimates
posterior_sd_res$sd[posterior_sd_res$model_name == "Direct" & posterior_sd_res$outcome == "beta1"] <- results_direct$se_beta1
posterior_sd_res$sd[posterior_sd_res$model_name == "Direct" & posterior_sd_res$outcome == "beta2"] <- results_direct$se_beta2

# save posterior SDs of area-level estimates
for (i in 1:length(stage_2_list)) {
  mod.preds <- stage_2_list[[i]]$summary.lincomb.derived
  beta1.pred <- mod.preds[grepl("outcome1", rownames(mod.preds)),]
  beta2.pred <- mod.preds[grepl("oucome2", rownames(mod.preds)),]
  posterior_sd_res$sd[posterior_sd_res$model_name == names(stage_2_list)[i] & posterior_sd_res$outcome == "beta1"] <- beta1.pred$sd
  posterior_sd_res$sd[posterior_sd_res$model_name == names(stage_2_list)[i] & posterior_sd_res$outcome == "beta2"] <- beta2.pred$sd
}

posterior_sd_res %<>% mutate(model_factor = fct_relevel(as.factor(model_name),
                                                        c("Direct", names(stage_2_list)))) %>%
    mutate(outcome = ifelse(outcome == "beta1", "modern vs. none", "other vs. none"))    

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

# make posterior sd plots ####
ggplot(posterior_sd_res, aes(x = sd, y = reorder(model_factor, sd, na.rm = TRUE), fill = model_factor)) + 
    geom_boxplot(show.legend = FALSE) + 
    scale_fill_viridis(option = "C", discrete = TRUE, alpha = 0.5) +
    facet_wrap(~ outcome, nrow = 1, scales = "free_x") +
    xlab("posterior SDs") +
    ylab("Model") +
    theme_light()

ggsave("posterior_sd_compare-multinomial.pdf", width = 7, height = 4)
    
# maps of results ####
n_cats <- 15
plot.palette <- viridis(n_cats)

all.results.wide <- all.results %>%
    select(admin1.name, model.name, measure, value) %>%
    pivot_wider(names_from = c(model.name, measure), values_from = value)
names(all.results.wide) <- gsub(" ", "_", names(all.results.wide))

## RE compare ####
# plot random effects
outcomes <- c("Bivariate_shared_BYM_beta1.iid.50", "Bivariate_shared_BYM_beta1.icar.50", "Bivariate_shared_BYM_beta1.totre.50",
              "Bivariate_shared_BYM_beta2.iid.50", "Bivariate_shared_BYM_beta2.icar.50", "Bivariate_shared_BYM_beta2.totre.50")
outcomes.names <- c("Logit(p) modern vs. none IID", "Logit(p) modern vs. none ICAR", "Logit(p) modern vs. none total BYM",
                    "Shared (WAZ) IID", "Shared (WAZ) ICAR", "Shared (WAZ) total BYM")

tmp <- merge(poly.adm1, all.results.wide[, c("admin1.name", outcomes)],
             by.x = "NAME_1", by.y = "admin1.name")

pdf("shared-bym-multinomial-iid-icar-compare.pdf", width = 12, height = 7)

grid.arrange(
    sp::spplot(obj = tmp, zcol = outcomes[c(6, 3)],
               col.regions = plot.palette,
               cuts = length(plot.palette) - 1,
               layout = c(1, 2),
               names.attr = c("Shared (modern vs. none)","other vs. none"),
               main = "Total BYM2"),
    sp::spplot(obj = tmp, zcol = outcomes[c(4, 5, 1, 2)],
               col.regions = plot.palette,
               cuts = length(plot.palette) - 1,
               layout = c(2, 2),
               names.attr = c("Shared (modern vs. none) IID", "Shared (modern vs. none) ICAR",
                              "other vs. none IID", "other vs. none ICAR"),
               main = "IID and ICAR separate"),
    ncol = 2
)

dev.off()

## post. med and widths ####

pdf("shared-bym-multinomial-pred.pdf", width = 9, height = 9)

outcomes <- c("Bivariate_shared_BYM_beta1.pred.50", "Bivariate_shared_BYM_beta1.pred.width95",
              "Bivariate_shared_BYM_beta2.pred.50", "Bivariate_shared_BYM_beta2.pred.width95")
outcomes.names <- c("Logit(p) of modern vs. none", "Logit(p) of modern vs. none",
                    "Logit(p) of other vs. none", "Logit(p) of other vs. none")

tmp <- merge(poly.adm1, all.results.wide[, c("admin1.name", outcomes[c(1, 3)])],
             by.x = "NAME_1", by.y = "admin1.name")
p1 <- sp::spplot(obj = tmp, zcol = outcomes[c(1, 3)],
           col.regions = plot.palette,
           cuts = length(plot.palette) - 1,
           layout = c(2, 1),
           names.attr = outcomes.names[c(1, 3)],
           main = "Posterior medians")

tmp <- merge(poly.adm1, all.results.wide[, c("admin1.name", outcomes[c(2, 4)])],
             by.x = "NAME_1", by.y = "admin1.name")
p2 <- sp::spplot(obj = tmp, zcol = outcomes[c(2, 4)],
           col.regions = plot.palette,
           cuts = length(plot.palette) - 1,
           layout = c(2, 1),
           names.attr = outcomes.names[c(2, 4)],
           main = "95% CI width")
grid.arrange(p1, p2)
dev.off()

## probabilities meds and widths ####
pdf("shared-bym-multinomial-probs-pred.pdf", width = 9, height = 6)

tmp <- merge(poly.adm1, posterior_probs_res,
             by.x = "NAME_1", by.y = "admin1.name")

ggarrange(
    sp::spplot(obj = tmp, zcol = c("prob_none_med"),
               col.regions = plot.palette,
               cuts = length(plot.palette) - 1,
               main = "None post. med"),
    sp::spplot(obj = tmp, zcol = c("prob_modern_med"),
               col.regions = plot.palette,
               cuts = length(plot.palette) - 1,
               main = "Modern post. med"),
    sp::spplot(obj = tmp, zcol = c("prob_other_med"),
               col.regions = plot.palette,
               cuts = length(plot.palette) - 1,
               main = "Other post. med"),
    sp::spplot(obj = tmp, zcol = c("prob_none_width95"),
               col.regions = plot.palette,
               cuts = length(plot.palette) - 1,
               main = "None 95% width"),
    sp::spplot(obj = tmp, zcol = c("prob_modern_width95"),
               col.regions = plot.palette,
               cuts = length(plot.palette) - 1,
               main = "Modern 95% width"),
    sp::spplot(obj = tmp, zcol = c("prob_other_width95"),
               col.regions = plot.palette,
               cuts = length(plot.palette) - 1,
               main = "Other 95% width"),
    ncol = 3, nrow = 2
)

dev.off()

# model comparison plots ####

## two-way comp ####

# compare median estimate
twoway_comp_plot(all.results %>% filter(model.name %in% c("Direct", "Univariate IID", "Univariate BYM")), "pred")
ggsave("med_estimate_compare_univariate-multinomial.pdf",
       width = 8.5, height = 3.6)

twoway_comp_plot(all.results %>% filter(model.name %in% c("Direct",  
                                                     "Bivariate nonshared IID", "Bivariate shared IID")), "pred")
ggsave("med_estimate_compare_bivariate_IID-multinomial.pdf",
       width = 8.5, height = 3.6)

twoway_comp_plot(all.results %>% filter(model.name %in% c("Direct",
                                                          "Bivariate nonshared BYM", "Bivariate shared BYM")), "pred")
ggsave("med_estimate_compare_bivariate_bym2-multinomial.pdf",
       width = 8.5, height = 3.6)

## IID vs ICAR plots ####
# ggplot(all.results %>% 
#            filter(measure %in% c("beta1.iid.50", "beta2.icar.50")) %>%
#            filter(model_name %in% c("Bivariate nonshared BYM2", "Bivariate shared BYM2")) %>%
#            pivot_wider(names_from = "measure", values_from = "value"),
#        aes(x = beta2.icar.50, y = beta1.iid.50)) +
#     geom_point() + 
#     geom_smooth(method = "lm", formula = y ~ x) +
#     facet_wrap(~ model_name, nrow = 1) +
#     ylab("beta1 IID") +
#     xlab("beta2 ICAR") +
#     theme_light()
# ggsave("beta1_iid_beta2_icar_compare_bivariate_nonshared_vs_shared-multinomial.pdf")

## beta1 vs beta2 preds ####
# ggplot(all.results %>% 
#            filter(measure %in% c("beta1.pred.50", "beta2.pred.50")) %>%
#            pivot_wider(names_from = "measure", values_from = "value"),
#        aes(x = beta2.pred.50, y = beta1.pred.50)) +
#     geom_point() + 
#     geom_smooth(method = "lm", formula = y ~ x) +
#     facet_wrap(~ model_factor, nrow = 1) +
#     ylab("beta1") +
#     xlab("beta2") +
#     theme_light()
# ggsave("multinomial_med_beta1_beta2_compare.pdf", width = 9, height = 5)
