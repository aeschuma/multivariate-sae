# Austin Schumacher
#

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

# functions ####

# plotting maps function
plot_maps <- function(poly, outcome.df, 
                      names.reg.poly, names.reg.outcomes, 
                      outcomes, outcomes.names, 
                      plot.palette) {
    tmp <- merge(poly, outcome.df[, c(names.reg.outcomes, outcomes)], 
                 by.x = names.reg.poly, by.y = names.reg.outcomes)
    plot.grid <- vector(mode = "list", length = length(outcomes))
    for(i in 1:length(outcomes.names)) {
        plot.grid[[i]] <- spplot(tmp, outcomes[i], 
                                 col.regions = plot.palette, cuts = length(plot.palette) - 1,
                                 main = outcomes.names[i])
    }
    do.call(grid.arrange, plot.grid)
}

# comparison scatterplot function
comparison_scatter <- function(data, measure, x_var, y_var) {
    # testing
    # data = data
    # measure = measures[1]
    # x_var = comp.mod.2[1]
    # y_var = comp.mod.1[1]
    
    data <- data %>% select("admin1.name", "measure", "model_name", "value")
    
    # object name can't be same name as variable
    my.measure <- measure 
    
    # pivot and subset data for this plot
    tmp <- data %>% 
        filter(model_name %in% c(x_var,  y_var)) %>%
        pivot_wider(names_from = model_name, values_from = value) %>%
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
        geom_point(alpha = 0.5) + geom_abline(slope = 1, intercept = 0, col = "darkgreen") + 
        xlab(labs[1]) +  ylab(labs[2]) +
        ggtitle(main,
                subtitle = paste0("Mean absolute diff: ", 
                                  mean.abs.diff)) +
        theme(axis.title = element_text(size=8)) +
        theme_light() 
}

# plot comparisons between 3 models x 2 outcomes (HAZ and WAZ)
twoway_comp_plot <- function(data, measure) {
    
    # testing
    # data <- results.b11.b33 %>% filter(model %in% c("univariate_b11", "univariate_b33"))
    # data <- all.results %>% filter(model %in% c("univariate", "bivariate.nonshared"))
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
        measures <- c("haz.pred.width80", "waz.pred.width80")
    } 
    
    # storage for list of plots
    models <- unique(data$model_name)
    nmodels <- length(models)
    ncomps <- choose(nmodels,2)
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
    med.comp.plots[[1]] <- comparison_scatter(data = data, 
                                              measure = measures[1], 
                                              x_var = comp.mod.2[1], 
                                              y_var = comp.mod.1[1])
    med.comp.plots[[2]] <- comparison_scatter(data = data, 
                                              measure = measures[2], 
                                              x_var = comp.mod.2[1], 
                                              y_var = comp.mod.1[1])
    counter <- 2
    if (ncomps >=2) {
        for (i in 2:ncomps) {
            counter <- counter + 1
            med.comp.plots[[counter]] <- comparison_scatter(data = data, 
                                                            measure = measures[1], 
                                                            x_var = comp.mod.2[i], 
                                                            y_var = comp.mod.1[i]) + 
                ggtitle("")
            counter <- counter + 1
            med.comp.plots[[counter]] <- comparison_scatter(data = data, 
                                                            measure = measures[2], 
                                                            x_var = comp.mod.2[i], 
                                                            y_var = comp.mod.1[i]) + 
                ggtitle("")
        }
    }
    
    # plot them all
    ggarrange(plotlist = med.comp.plots, ncol = 2, nrow = ncomps)
}

# load raw data ####
load("../../../../dissertation_2/survey-csmf/data/ken_dhs2014/data/haz-waz-kenDHS2014.rda")

# prelim plots ####
ggplot(dat, aes(x=HAZ)) + 
    geom_histogram() +
    ggtitle("Height for age z-score distribution for BGD DHS 2017") +
    theme_light()

ggplot(dat, aes(x=WAZ)) + 
    geom_histogram() +
    ggtitle("Weight for age z-score distribution for BGD DHS 2017") +
    theme_light()

print("Number of observed children in each Admin 1")
table(dat$admin1.name)

# direct estimates ####
stage_1_list <- read_rds("../../../../dissertation_2/survey-csmf/results/ken2014-hazwaz/ken2014-hazwaz-stage-1.rds")
results_direct <- stage_1_list$results

results <- results_direct
stage_1_list_uni <- read_rds("../../../../dissertation_2/survey-csmf/results/ken2014-hazwaz/ken2014-hazwaz-stage-1-univariate.rds")
results_direct_uni <- stage_1_list_uni$results

# results storage for posterior SDs
model_name <- c("Univariate direct", "Univariate IID", "Univariate BYM2", 
           "Bivariate direct",  
           "Bivariate nonshared IID", "Bivariate shared IID",
           "Bivariate nonshared BYM2", "Bivariate shared BYM2")
outcome <- c("haz", "waz")
reg <- 1:47
posterior_sd_res <- expand_grid(model_name, outcome, reg)
posterior_sd_res$sd <- NA

# store direct estimates
posterior_sd_res$sd[posterior_sd_res$model_name == "Univariate direct" & posterior_sd_res$outcome == "haz"] <- results_direct_uni$seHAZ.uni
posterior_sd_res$sd[posterior_sd_res$model_name == "Univariate direct" & posterior_sd_res$outcome == "waz"] <- results_direct_uni$seWAZ.uni
posterior_sd_res$sd[posterior_sd_res$model_name == "Bivariate direct" & posterior_sd_res$outcome == "haz"] <- results_direct$seHAZ.bi
posterior_sd_res$sd[posterior_sd_res$model_name == "Bivariate direct" & posterior_sd_res$outcome == "waz"] <- results_direct$seWAZ.bi

# univariate models ####

## IID ####

# load results
stage_2_uni_iid_list <- read_rds("../../../../dissertation_2/survey-csmf/results/ken2014-hazwaz/ken2014-hazwaz-stage-2-univariate-iid.rds")
mod.stan.haz.uni.iid <- stage_2_uni_iid_list$mod[[1]]
mod.stan.waz.uni.iid <- stage_2_uni_iid_list$mod[[2]]

params_to_extract <- c("beta", "sigma", "preds")
re_params <- c("v")

# summaries HAZ
mod.stan.summary.haz.uni.iid <- summary(mod.stan.haz.uni.iid,
                                        pars = params_to_extract,
                                        probs = c(0.1, 0.5, 0.9))

# summaries WAZ
mod.stan.summary.waz.uni.iid <- summary(mod.stan.waz.uni.iid,
                                        pars = params_to_extract,
                                        probs = c(0.1, 0.5, 0.9))

# extract REs and posterior widths
haz.re <- summary(mod.stan.haz.uni.iid, pars = re_params, probs = c(0.1, 0.5, 0.9))$summary
waz.re <- summary(mod.stan.waz.uni.iid, pars = re_params, probs = c(0.1, 0.5, 0.9))$summary

haz.re.iid <- haz.re[grep("v", rownames(haz.re)), c("50%")]
waz.re.iid <- waz.re[grep("v", rownames(waz.re)), c("50%")]

haz.re.icar <- rep(NA, length(haz.re.iid))
waz.re.icar <- rep(NA, length(waz.re.iid))

haz.pred <- summary(mod.stan.haz.uni.iid, pars = "preds", probs = c(0.1, 0.5, 0.9))$summary[, c("10%", "50%", "90%")]
waz.pred <- summary(mod.stan.waz.uni.iid, pars = "preds", probs = c(0.1, 0.5, 0.9))$summary[, c("10%", "50%", "90%")]

# create a data.frame with these results
res.df.uni.iid <- as.data.frame(cbind(results$admin1.name,
                                      haz.re.iid, waz.re.iid,
                                      haz.re.icar, waz.re.icar,
                                      haz.pred[, "50%"], haz.pred[, "90%"] - haz.pred[, "10%"], 
                                      waz.pred[, "50%"], waz.pred[, "90%"] - waz.pred[, "10%"]))

names(res.df.uni.iid) <- c("admin1.name", 
                           "haz.iid.50", "waz.iid.50",
                           "haz.icar.50", "waz.icar.50",
                           "haz.pred.50", "haz.pred.width80",
                           "waz.pred.50", "waz.pred.width80")
rownames(res.df.uni.iid) <- NULL

for (vv in 1:ncol(res.df.uni.iid)) {
    if (names(res.df.uni.iid)[vv] == "admin1.name") next
    res.df.uni.iid[, vv] <- as.numeric(res.df.uni.iid[, vv])
}

# save posterior SDs of area-level estimates
posterior_sd_res$sd[posterior_sd_res$model_name == "Univariate IID" & posterior_sd_res$outcome == "haz"] <- mod.stan.summary.haz.uni.iid$summary[grepl("preds", rownames(mod.stan.summary.haz.uni.iid$summary)), "sd"]
posterior_sd_res$sd[posterior_sd_res$model_name == "Univariate IID" & posterior_sd_res$outcome == "waz"] <- mod.stan.summary.waz.uni.iid$summary[grepl("preds", rownames(mod.stan.summary.waz.uni.iid$summary)), "sd"]

## BYM2 ####

# load results
stage_2_uni_b11_list <- read_rds("../../../../dissertation_2/survey-csmf/results/ken2014-hazwaz/ken2014-hazwaz-stage-2-univariate-b11.rds")
mod.stan.haz.uni.b11 <- stage_2_uni_b11_list$mod[[1]]
mod.stan.waz.uni.b11 <- stage_2_uni_b11_list$mod[[2]]

params_to_extract <- c("beta", "sigma", "rho", "preds")
re_params <- c("v", "u")

# summaries HAZ
mod.stan.summary.haz.uni.b11 <- summary(mod.stan.haz.uni.b11,
                                        pars = params_to_extract,
                                        probs = c(0.1, 0.5, 0.9))

# summaries WAZ
mod.stan.summary.waz.uni.b11 <- summary(mod.stan.waz.uni.b11,
                                        pars = params_to_extract,
                                        probs = c(0.1, 0.5, 0.9))

# extract REs and posterior widths
haz.re <- summary(mod.stan.haz.uni.b11, pars = re_params, probs = c(0.1, 0.5, 0.9))$summary
waz.re <- summary(mod.stan.waz.uni.b11, pars = re_params, probs = c(0.1, 0.5, 0.9))$summary

haz.re.iid <- haz.re[grep("v", rownames(haz.re)), c("50%")]
waz.re.iid <- waz.re[grep("v", rownames(waz.re)), c("50%")]

haz.re.icar <- haz.re[grep("u", rownames(haz.re)), c("50%")]
waz.re.icar <- waz.re[grep("u", rownames(waz.re)), c("50%")]

haz.pred <- summary(mod.stan.haz.uni.b11, pars = "preds", probs = c(0.1, 0.5, 0.9))$summary[, c("10%", "50%", "90%")]
waz.pred <- summary(mod.stan.waz.uni.b11, pars = "preds", probs = c(0.1, 0.5, 0.9))$summary[, c("10%", "50%", "90%")]

# create a data.frame with these results
res.df.uni.b11 <- as.data.frame(cbind(results$admin1.name,
                                      haz.re.iid, waz.re.iid,
                                      haz.re.icar, waz.re.icar,
                                      haz.pred[, "50%"], haz.pred[, "90%"] - haz.pred[, "10%"], 
                                      waz.pred[, "50%"], waz.pred[, "90%"] - waz.pred[, "10%"]))

names(res.df.uni.b11) <- c("admin1.name", 
                           "haz.iid.50", "waz.iid.50",
                           "haz.icar.50", "waz.icar.50",
                           "haz.pred.50", "haz.pred.width80",
                           "waz.pred.50", "waz.pred.width80")
rownames(res.df.uni.b11) <- NULL

for (vv in 1:ncol(res.df.uni.b11)) {
    if (names(res.df.uni.b11)[vv] == "admin1.name") next
    res.df.uni.b11[, vv] <- as.numeric(res.df.uni.b11[, vv])
}

# save posterior SDs of area-level estimates
posterior_sd_res$sd[posterior_sd_res$model_name == "Univariate BYM2" & posterior_sd_res$outcome == "haz"] <- mod.stan.summary.haz.uni.b11$summary[grepl("preds", rownames(mod.stan.summary.haz.uni.b11$summary)), "sd"]
posterior_sd_res$sd[posterior_sd_res$model_name == "Univariate BYM2" & posterior_sd_res$outcome == "waz"] <- mod.stan.summary.waz.uni.b11$summary[grepl("preds", rownames(mod.stan.summary.waz.uni.b11$summary)), "sd"]

# bivariate nonshared ####

## IID ####

# load results
stage_2_bi_nonshared_iid_list <- read_rds("../../../../dissertation_2/survey-csmf/results/ken2014-hazwaz/ken2014-hazwaz-stage-2-bivariate-nonshared-iid.rds")
mod.stan.bi.nonshared.iid <- stage_2_bi_nonshared_iid_list$mod[[1]]

params_to_extract <- c("beta", "sigma", "preds")
re_params <- c("v_1", "v_2")

# summaries
mod.stan.summary.bi.nonshared.iid <- summary(mod.stan.bi.nonshared.iid,
                                             pars = params_to_extract,
                                             probs = c(0.1, 0.5, 0.9))

# extract REs and posterior widths
res.re <- summary(mod.stan.bi.nonshared.iid, pars = re_params, probs = c(0.1, 0.5, 0.9))$summary

haz.re.iid <- res.re[grep("v_1", rownames(res.re)), c("50%")]
waz.re.iid <- res.re[grep("v_2", rownames(res.re)), c("50%")]

haz.re.icar <- rep(NA, length(haz.re.iid))
waz.re.icar <- rep(NA, length(waz.re.iid))

res.pred <- summary(mod.stan.bi.nonshared.iid, pars = "preds", probs = c(0.1, 0.5, 0.9))$summary[, c("10%", "50%", "90%")]

haz.pred <- res.pred[grep("preds\\[.*,1", rownames(res.pred)), ]
waz.pred <- res.pred[grep("preds\\[.*,2", rownames(res.pred)), ]

# create a data.frame with these results
res.df.bi.nonshared.iid <- as.data.frame(cbind(results$admin1.name,
                                               haz.re.iid, waz.re.iid,
                                               haz.re.icar, waz.re.icar,
                                               haz.pred[, "50%"], haz.pred[, "90%"] - haz.pred[, "10%"], 
                                               waz.pred[, "50%"], waz.pred[, "90%"] - waz.pred[, "10%"]))

names(res.df.bi.nonshared.iid) <- c("admin1.name", 
                                    "haz.iid.50", "waz.iid.50",
                                    "haz.icar.50", "waz.icar.50",
                                    "haz.pred.50", "haz.pred.width80",
                                    "waz.pred.50", "waz.pred.width80")
rownames(res.df.bi.nonshared.iid) <- NULL

for (vv in 1:ncol(res.df.bi.nonshared.iid)) {
    if (names(res.df.bi.nonshared.iid)[vv] == "admin1.name") next
    res.df.bi.nonshared.iid[, vv] <- as.numeric(res.df.bi.nonshared.iid[, vv])
}

# save posterior SDs of area-level estimates
varnames.tmp <- rownames(mod.stan.summary.bi.nonshared.iid$summary)
tmp.post.sd.haz <- mod.stan.summary.bi.nonshared.iid$summary[grepl("preds", varnames.tmp) & grepl(",1]", varnames.tmp), "sd"]
tmp.post.sd.waz <- mod.stan.summary.bi.nonshared.iid$summary[grepl("preds", varnames.tmp) & grepl(",2]", varnames.tmp), "sd"]
posterior_sd_res$sd[posterior_sd_res$model_name == "Bivariate nonshared IID" & posterior_sd_res$outcome == "haz"] <- tmp.post.sd.haz
posterior_sd_res$sd[posterior_sd_res$model_name == "Bivariate nonshared IID" & posterior_sd_res$outcome == "waz"] <- tmp.post.sd.waz

## BYM2 ####

# load results
stage_2_bi_nonshared_b11_list <- read_rds("../../../../dissertation_2/survey-csmf/results/ken2014-hazwaz/ken2014-hazwaz-stage-2-bivariate-nonshared-b11.rds")
mod.stan.bi.nonshared.b11 <- stage_2_bi_nonshared_b11_list$mod[[1]]

params_to_extract <- c("beta", "sigma", "rho", "preds")
re_params <- c("v_1", "v_2", "u_1", "u_2")

# summaries
mod.stan.summary.bi.nonshared.b11 <- summary(mod.stan.bi.nonshared.b11,
                                             pars = params_to_extract,
                                             probs = c(0.1, 0.5, 0.9))

# extract REs and posterior widths
res.re <- summary(mod.stan.bi.nonshared.b11, pars = re_params, probs = c(0.1, 0.5, 0.9))$summary

haz.re.iid <- res.re[grep("v_1", rownames(res.re)), c("50%")]
waz.re.iid <- res.re[grep("v_2", rownames(res.re)), c("50%")]

haz.re.icar <- res.re[grep("u_1", rownames(res.re)), c("50%")]
waz.re.icar <- res.re[grep("u_2", rownames(res.re)), c("50%")]

res.pred <- summary(mod.stan.bi.nonshared.b11, pars = "preds", probs = c(0.1, 0.5, 0.9))$summary[, c("10%", "50%", "90%")]

haz.pred <- res.pred[grep("preds\\[.*,1", rownames(res.pred)), ]
waz.pred <- res.pred[grep("preds\\[.*,2", rownames(res.pred)), ]

# create a data.frame with these results
res.df.bi.nonshared.b11 <- as.data.frame(cbind(results$admin1.name,
                                               haz.re.iid, waz.re.iid,
                                               haz.re.icar, waz.re.icar,
                                               haz.pred[, "50%"], haz.pred[, "90%"] - haz.pred[, "10%"], 
                                               waz.pred[, "50%"], waz.pred[, "90%"] - waz.pred[, "10%"]))

names(res.df.bi.nonshared.b11) <- c("admin1.name", 
                                    "haz.iid.50", "waz.iid.50",
                                    "haz.icar.50", "waz.icar.50",
                                    "haz.pred.50", "haz.pred.width80",
                                    "waz.pred.50", "waz.pred.width80")
rownames(res.df.bi.nonshared.b11) <- NULL

for (vv in 1:ncol(res.df.bi.nonshared.b11)) {
    if (names(res.df.bi.nonshared.b11)[vv] == "admin1.name") next
    res.df.bi.nonshared.b11[, vv] <- as.numeric(res.df.bi.nonshared.b11[, vv])
}

# save posterior SDs of area-level estimates
varnames.tmp <- rownames(mod.stan.summary.bi.nonshared.b11$summary)
tmp.post.sd.haz <- mod.stan.summary.bi.nonshared.b11$summary[grepl("preds", varnames.tmp) & grepl(",1]", varnames.tmp), "sd"]
tmp.post.sd.waz <- mod.stan.summary.bi.nonshared.b11$summary[grepl("preds", varnames.tmp) & grepl(",2]", varnames.tmp), "sd"]
posterior_sd_res$sd[posterior_sd_res$model_name == "Bivariate nonshared BYM2" & posterior_sd_res$outcome == "haz"] <- tmp.post.sd.haz
posterior_sd_res$sd[posterior_sd_res$model_name == "Bivariate nonshared BYM2" & posterior_sd_res$outcome == "waz"] <- tmp.post.sd.waz

# bivariate shared models ####

## IID ####

# load results
stage_2_bi_shared_iid_list <- read_rds("../../../../dissertation_2/survey-csmf/results/ken2014-hazwaz/ken2014-hazwaz-stage-2-bivariate-shared-iid.rds")
mod.stan.bi.shared.iid <- stage_2_bi_shared_iid_list$mod[[1]]

params_to_extract <- c("beta", "sigma", "lambda", "preds")
re_params <- c("v_1", "v_2")

# summaries
mod.stan.summary.bi.shared.iid <- summary(mod.stan.bi.shared.iid,
                                          pars = params_to_extract,
                                          probs = c(0.1, 0.5, 0.9))

# extract REs and posterior widths
res.re <- summary(mod.stan.bi.shared.iid, pars = re_params, probs = c(0.1, 0.5, 0.9))$summary

haz.re.iid <- res.re[grep("v_1", rownames(res.re)), c("50%")]
waz.re.iid <- res.re[grep("v_2", rownames(res.re)), c("50%")]

haz.re.icar <- rep(NA, length(haz.re.iid))
waz.re.icar <- rep(NA, length(waz.re.iid))

res.pred <- summary(mod.stan.bi.shared.iid, pars = "preds", probs = c(0.1, 0.5, 0.9))$summary[, c("10%", "50%", "90%")]

haz.pred <- res.pred[grep("preds\\[.*,1", rownames(res.pred)), ]
waz.pred <- res.pred[grep("preds\\[.*,2", rownames(res.pred)), ]

# create a data.frame with these results
res.df.bi.shared.iid <- as.data.frame(cbind(results$admin1.name,
                                            haz.re.iid, waz.re.iid,
                                            haz.re.icar, waz.re.icar,
                                            haz.pred[, "50%"], haz.pred[, "90%"] - haz.pred[, "10%"], 
                                            waz.pred[, "50%"], waz.pred[, "90%"] - waz.pred[, "10%"]))

# no longer saving name as "shared.icar.50"
names(res.df.bi.shared.iid) <- c("admin1.name", 
                                 "haz.iid.50", "waz.iid.50",
                                 "haz.icar.50", "waz.icar.50",
                                 "haz.pred.50", "haz.pred.width80",
                                 "waz.pred.50", "waz.pred.width80")
rownames(res.df.bi.shared.iid) <- NULL

for (vv in 1:ncol(res.df.bi.shared.iid)) {
    if (names(res.df.bi.shared.iid)[vv] == "admin1.name") next
    res.df.bi.shared.iid[, vv] <- as.numeric(res.df.bi.shared.iid[, vv])
}

# save posterior SDs of area-level estimates
varnames.tmp <- rownames(mod.stan.summary.bi.shared.iid$summary)
tmp.post.sd.haz <- mod.stan.summary.bi.shared.iid$summary[grepl("preds", varnames.tmp) & grepl(",1]", varnames.tmp), "sd"]
tmp.post.sd.waz <- mod.stan.summary.bi.shared.iid$summary[grepl("preds", varnames.tmp) & grepl(",2]", varnames.tmp), "sd"]
posterior_sd_res$sd[posterior_sd_res$model_name == "Bivariate shared IID" & posterior_sd_res$outcome == "haz"] <- tmp.post.sd.haz
posterior_sd_res$sd[posterior_sd_res$model_name == "Bivariate shared IID" & posterior_sd_res$outcome == "waz"] <- tmp.post.sd.waz

## BYM2 ####

# load results
stage_2_bi_shared_b11_list <- read_rds("../../../../dissertation_2/survey-csmf/results/ken2014-hazwaz/ken2014-hazwaz-stage-2-bivariate-shared-coreg-b11.rds")
mod.stan.bi.shared.b11 <- stage_2_bi_shared_b11_list$mod[[1]]

params_to_extract <- c("beta", "sigma", "rho", "lambda", "preds")
re_params <- c("v_1", "v_2", "u_1", "u_2")

# summaries
mod.stan.summary.bi.shared.b11 <- summary(mod.stan.bi.shared.b11,
                                          pars = params_to_extract,
                                          probs = c(0.1, 0.5, 0.9))

# extract REs and posterior widths
res.re <- summary(mod.stan.bi.shared.b11, pars = re_params, probs = c(0.1, 0.5, 0.9))$summary

haz.re.iid <- res.re[grep("v_1", rownames(res.re)), c("50%")]
waz.re.iid <- res.re[grep("v_2", rownames(res.re)), c("50%")]

haz.re.icar <- res.re[grep("u_1", rownames(res.re)), c("50%")]
waz.re.icar <- res.re[grep("u_2", rownames(res.re)), c("50%")]

res.pred <- summary(mod.stan.bi.shared.b11, pars = "preds", probs = c(0.1, 0.5, 0.9))$summary[, c("10%", "50%", "90%")]

haz.pred <- res.pred[grep("preds\\[.*,1", rownames(res.pred)), ]
waz.pred <- res.pred[grep("preds\\[.*,2", rownames(res.pred)), ]

# create a data.frame with these results
res.df.bi.shared.b11 <- as.data.frame(cbind(results$admin1.name,
                                            haz.re.iid, waz.re.iid,
                                            haz.re.icar, waz.re.icar,
                                            haz.pred[, "50%"], haz.pred[, "90%"] - haz.pred[, "10%"], 
                                            waz.pred[, "50%"], waz.pred[, "90%"] - waz.pred[, "10%"]))

# no longer saving name as "shared.icar.50"
names(res.df.bi.shared.b11) <- c("admin1.name", 
                                 "haz.iid.50", "waz.iid.50",
                                 "haz.icar.50", "waz.icar.50",
                                 "haz.pred.50", "haz.pred.width80",
                                 "waz.pred.50", "waz.pred.width80")
rownames(res.df.bi.shared.b11) <- NULL

for (vv in 1:ncol(res.df.bi.shared.b11)) {
    if (names(res.df.bi.shared.b11)[vv] == "admin1.name") next
    res.df.bi.shared.b11[, vv] <- as.numeric(res.df.bi.shared.b11[, vv])
}

# save posterior SDs of area-level estimates
varnames.tmp <- rownames(mod.stan.summary.bi.shared.b11$summary)
tmp.post.sd.haz <- mod.stan.summary.bi.shared.b11$summary[grepl("preds", varnames.tmp) & grepl(",1]", varnames.tmp), "sd"]
tmp.post.sd.waz <- mod.stan.summary.bi.shared.b11$summary[grepl("preds", varnames.tmp) & grepl(",2]", varnames.tmp), "sd"]
posterior_sd_res$sd[posterior_sd_res$model_name == "Bivariate shared BYM2" & posterior_sd_res$outcome == "haz"] <- tmp.post.sd.haz
posterior_sd_res$sd[posterior_sd_res$model_name == "Bivariate shared BYM2" & posterior_sd_res$outcome == "waz"] <- tmp.post.sd.waz

# formatting all results ####

par.names <- c("$\\alpha_1$",
               "$\\alpha_2$",
               "$\\sigma_1$",
               "$\\sigma_2$",
               "$\\rho_1$",
               "$\\rho_2$",
               "$\\lambda$")
uni.par.ests.b11 <- as_tibble(rbind(mod.stan.summary.haz.uni.b11$summary["beta",],
                                    mod.stan.summary.waz.uni.b11$summary["beta",],
                                    mod.stan.summary.haz.uni.b11$summary["sigma",],
                                    mod.stan.summary.waz.uni.b11$summary["sigma",],
                                    mod.stan.summary.haz.uni.b11$summary["rho",],
                                    mod.stan.summary.waz.uni.b11$summary["rho",],
                                    rep(NA, ncol(mod.stan.summary.haz.uni.b11$summary))))
uni.par.ests.b11$parameter <- par.names
uni.par.ests.b11$model_name <- "univariate bym2"
uni.par.ests.b11$model <- "univariate.bym2"

uni.par.ests.iid <- as_tibble(rbind(mod.stan.summary.haz.uni.iid$summary["beta",],
                                    mod.stan.summary.waz.uni.iid$summary["beta",],
                                    mod.stan.summary.haz.uni.iid$summary["sigma",],
                                    mod.stan.summary.waz.uni.iid$summary["sigma",],
                                    rep(NA, ncol(mod.stan.summary.haz.uni.iid$summary)),
                                    rep(NA, ncol(mod.stan.summary.haz.uni.iid$summary)),
                                    rep(NA, ncol(mod.stan.summary.haz.uni.iid$summary))))
uni.par.ests.iid$parameter <- par.names
uni.par.ests.iid$model_name <- "univariate iid"
uni.par.ests.iid$model <- "univariate.iid"

nonshared.par.ests.b11 <- as_tibble(
    rbind(mod.stan.summary.bi.nonshared.b11$summary[!grepl("preds", rownames(mod.stan.summary.bi.nonshared.b11$summary)),],
          rep(NA, ncol(mod.stan.summary.bi.nonshared.b11$summary)))
)
nonshared.par.ests.b11$parameter <- par.names
nonshared.par.ests.b11$model_name <- "bivariate nonshared bym2"
nonshared.par.ests.b11$model <- "bivariate.nonshared.bym2"

nonshared.par.ests.iid <- as_tibble(
    rbind(mod.stan.summary.bi.nonshared.iid$summary[!grepl("preds", rownames(mod.stan.summary.bi.nonshared.iid$summary)),],
          rep(NA, ncol(mod.stan.summary.bi.nonshared.iid$summary)),
          rep(NA, ncol(mod.stan.summary.bi.nonshared.iid$summary)),
          rep(NA, ncol(mod.stan.summary.bi.nonshared.iid$summary)))
)
nonshared.par.ests.iid$parameter <- par.names
nonshared.par.ests.iid$model_name <- "bivariate nonshared iid"
nonshared.par.ests.iid$model <- "bivariate.nonshared.iid"

shared.par.ests.b11 <- as_tibble(mod.stan.summary.bi.shared.b11$summary[!grepl("preds", rownames(mod.stan.summary.bi.shared.b11$summary)),])
shared.par.ests.b11$parameter <- par.names
shared.par.ests.b11$model_name <- "bivariate shared"
shared.par.ests.b11$model <- "bivariate.shared"

shared.par.ests.iid <- as_tibble(rbind(mod.stan.summary.bi.shared.iid$summary[!grepl("preds", rownames(mod.stan.summary.bi.shared.iid$summary)),],
                                       rep(NA, ncol(mod.stan.summary.bi.shared.iid$summary)),
                                       rep(NA, ncol(mod.stan.summary.bi.shared.iid$summary))))
shared.par.ests.iid$parameter <- par.names
shared.par.ests.iid$model_name <- "bivariate shared iid"
shared.par.ests.iid$model <- "bivariate.shared.iid"

all.par.ests <- bind_rows(uni.par.ests.b11,
                          uni.par.ests.iid,
                          nonshared.par.ests.b11,
                          nonshared.par.ests.iid,
                          shared.par.ests.b11,
                          shared.par.ests.iid) %>%
    relocate(c(model, parameter))

model_name_df <- all.par.ests %>% select(model, model_name) %>% unique()

# format direct estimates (first stage)
res.direct <- results_direct %>% mutate(haz.pred.width80_direct = qnorm(0.9) * 2 * seHAZ.bi,
                                 waz.pred.width80_direct = qnorm(0.9) * 2 * seWAZ.bi) %>%
    select(admin1.name, meanHAZ.bi, meanWAZ.bi, haz.pred.width80_direct, waz.pred.width80_direct) %>%
    rename(haz.pred.50_direct = meanHAZ.bi,
           waz.pred.50_direct = meanWAZ.bi)

res.direct.uni <- results_direct_uni %>% mutate(haz.pred.width80_direct.uni = qnorm(0.9) * 2 * seHAZ.uni,
                                 waz.pred.width80_direct.uni = qnorm(0.9) * 2 * seWAZ.uni) %>%
    select(admin1.name, meanHAZ.uni, meanWAZ.uni, haz.pred.width80_direct.uni, waz.pred.width80_direct.uni) %>%
    rename(haz.pred.50_direct.uni = meanHAZ.uni,
           waz.pred.50_direct.uni = meanWAZ.uni)

# combine results into one dataset
names(res.df.uni.iid) <- c("admin1.name", 
                           paste0(names(res.df.uni.iid)[-1], 
                                  "_univariate.iid"))
names(res.df.uni.b11) <- c("admin1.name", 
                           paste0(names(res.df.uni.b11)[-1], 
                                  "_univariate.bym2"))
names(res.df.bi.nonshared.iid) <- c("admin1.name", 
                                    paste0(names(res.df.bi.nonshared.iid)[-1], 
                                           "_bivariate.nonshared.iid"))
names(res.df.bi.nonshared.b11) <- c("admin1.name", 
                                    paste0(names(res.df.bi.nonshared.b11)[-1], 
                                           "_bivariate.nonshared.bym2"))
names(res.df.bi.shared.iid) <- c("admin1.name", 
                                 paste0(names(res.df.bi.shared.iid)[-1], 
                                        "_bivariate.shared.iid"))
names(res.df.bi.shared.b11) <- c("admin1.name", 
                                 paste0(names(res.df.bi.shared.b11)[-1], 
                                        "_bivariate.shared.bym2"))

# reference for model names
model_name_dat <- tibble(model = c("direct.uni", "univariate.iid", "univariate.bym2", 
                                   "direct",  
                                   "bivariate.nonshared.iid", "bivariate.shared.iid",
                                   "bivariate.nonshared.bym2", "bivariate.shared.bym2"),
                         model_name = c("Univariate direct", "Univariate IID", "Univariate BYM2", 
                                        "Bivariate direct",  
                                        "Bivariate nonshared IID", "Bivariate shared IID",
                                        "Bivariate nonshared BYM2", "Bivariate shared BYM2"))

all.results <- res.direct %>% 
    left_join(res.direct.uni) %>%
    left_join(res.df.uni.iid) %>% 
    left_join(res.df.uni.b11) %>% 
    left_join(res.df.bi.nonshared.iid) %>%
    left_join(res.df.bi.nonshared.b11) %>%
    left_join(res.df.bi.shared.iid) %>%
    left_join(res.df.bi.shared.b11) %>%
    pivot_longer(col = !admin1.name,
                 names_to = c("measure", "model"),
                 names_sep = "_",
                 values_to = "value") %>%
    mutate(model_factor = fct_relevel(as.factor(model),
                                      c("direct.uni", "univariate.iid", "univariate.bym2", 
                                        "direct",  
                                        "bivariate.nonshared.iid", "bivariate.shared.iid",
                                        "bivariate.nonshared.bym2", "bivariate.shared.bym2"))) %>%
    left_join(model_name_dat) %>%
    mutate(model_name_factor = fct_relevel(as.factor(model_name),
                                           c("Univariate direct", "Univariate IID", "Univariate BYM2", 
                                             "Bivariate direct",  
                                             "Bivariate nonshared IID", "Bivariate shared IID",
                                             "Bivariate nonshared BYM2", "Bivariate shared BYM2")))

posterior_sd_res %<>% mutate(model_name_factor = fct_relevel(as.factor(model_name),
                                                             c("Univariate direct", "Univariate IID", "Univariate BYM2", 
                                                               "Bivariate direct",  
                                                               "Bivariate nonshared IID", "Bivariate shared IID",
                                                               "Bivariate nonshared BYM2", "Bivariate shared BYM2")))

# make posterior sd plots ####
ggplot(posterior_sd_res %>% filter(outcome == "haz"), aes(x = sd, y = reorder(model_name_factor, sd, na.rm = TRUE), fill = model_name_factor)) + 
    geom_boxplot(show.legend = FALSE) + 
    scale_fill_viridis(option = "C", discrete = TRUE, alpha = 0.5) +
    xlab("posterior SDs") +
    ylab("Model") +
    theme_light()

ggsave("posterior_sd_compare.pdf")
    
# maps of results ####
n_cats <- 15
plot.palette <- viridis(n_cats)

# plot smoothed direct IID and BYM2
tmpdat <- all.results %>% filter(model %in% c("univariate.iid", "univariate.bym2")) %>% 
    select(admin1.name, model, measure, value) %>%
    pivot_wider(names_from = c(model, measure), values_from = value)
outcomes <- c("univariate.iid_haz.pred.50", "univariate.bym2_haz.pred.50",
              "univariate.iid_waz.pred.50", "univariate.bym2_waz.pred.50")
outcomes.names <- c("Univariate IID: HAZ", "Univariate BYM2: HAZ",
                    "Univariate IID: WAZ", "Univariate BYM2: WAZ")

pdf("uni-haz-waz-iid-bym2-med-compare.pdf", width = 9, height = 9)
plot_maps(poly = poly.adm1, outcome.df = tmpdat, 
          names.reg.poly = "NAME_1", names.reg.outcomes = "admin1.name", 
          outcomes = outcomes, outcomes.names = outcomes.names, 
          plot.palette = plot.palette)
dev.off()

# make model comparison plots ####

# compare median estimate
twoway_comp_plot(all.results %>% filter(model %in% c("direct.uni", "univariate.iid", "univariate.bym2")), "pred")
ggsave("med_estimate_compare_univariate.pdf",
       width = 8, height = 10)

twoway_comp_plot(all.results %>% filter(model %in% c("direct.uni", "univariate.iid", "univariate.bym2")), "pred")
ggsave("med_estimate_compare_univariate.pdf",
       width = 8, height = 10)

twoway_comp_plot(all.results %>% filter(model_name %in% c("Bivariate direct",  
                                                     "Bivariate nonshared IID", "Bivariate shared IID")), "pred")
ggsave("med_estimate_compare_bivariate_IID.pdf",
       width = 8, height = 10)

twoway_comp_plot(all.results %>% filter(model_name %in% c("Bivariate direct",
                                                          "Bivariate nonshared BYM2", "Bivariate shared BYM2")), "pred")
ggsave("med_estimate_compare_bivariate_bym2.pdf",
       width = 8, height = 10)

# IID vs ICAR plots ####
ggplot(all.results %>% 
           filter(measure %in% c("haz.iid.50", "waz.icar.50")) %>%
           filter(model_name %in% c("Bivariate nonshared BYM2", "Bivariate shared BYM2")) %>%
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
