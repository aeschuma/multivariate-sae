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


## ----functions, cache = FALSE----------------------------------------------------------------
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
    # object name can't be same name as variable
    my.measure <- measure 
    
    # pivot and subset data for this plot
    tmp <- data %>% 
        filter(model %in% c(x_var,  y_var)) %>%
        pivot_wider(names_from = model, values_from = value) %>%
        filter(measure == my.measure)
    
    xyvars <- c(x_var, y_var)
    labs <- c(NA, NA)
    for (i in 1:2) {
        if (xyvars[i] == "univariate") {
            labs[i] <- "Univariate"
        } else if (xyvars[i] == "bivariate.nonshared") {
            labs[i] <- "Bivariate Nonshared"
        } else if (xyvars[i] == "bivariate.shared") {
            labs[i] <- "Bivariate Shared (WAZ)"
        } else if (xyvars[i] == "bivariate.shared.2") {
            labs[i] <- "Bivariate Shared (HAZ)"
        } else if (xyvars[i] == "direct") {
            labs[i] <- "Direct"
        }
    } 
    
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
    models <- unique(data$model)
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


## ----read-data-stage-1-----------------------------------------------------------------------
load("/Users/austin/Dropbox/dissertation_2/survey-csmf/data/ken_dhs2014/data/haz-waz-kenDHS2014.rda")
stage_1_list <- read_rds("/Users/austin/Dropbox/dissertation_2/survey-csmf/results/ken2014-hazwaz/ken2014-hazwaz-stage-1.rds")
results <- stage_1_list[["results"]]
V.array <- stage_1_list[["V.array"]]
n_regions <- nrow(poly.adm1)

## ----bi-nonshared-mod------------------------------------------------------------------------

# hyperpars
rho_beta_1 <- 1
rho_beta_2 <- 1
sigma_normal_sd <- 1

datlist <- list(R = n_regions,
                regions = results$admin1,
                Sigma = V.array,
                y = results[, c("meanHAZ.bi", "meanWAZ.bi")],
                N_edges = node.info$n_edges,
                node1 = node.info$node1,
                node2 = node.info$node2,
                scaling_factor = scaling_factor,
                rho_beta_1 = rho_beta_1,
                rho_beta_2 = rho_beta_2,
                sigma_normal_sd = sigma_normal_sd)


stan_file <- "stan-models/nonshared-bym2.stan"

## stan options
niter <- 10000
nchains <- 3
prop_warmup <- 0.5
max_treedepth <- 250
adapt_delta <- 0.9
nthin <- 1
options(mc.cores = parallel::detectCores())

cmd_mod <- cmdstan_model(stan_file = stan_file)
fit <- cmd_mod$sample(data = datlist,
                      iter_warmup = niter*prop_warmup, iter_sampling = niter*(1-prop_warmup),
                      chains = nchains, thin = nthin,
                      adapt_delta = adapt_delta, max_treedepth = max_treedepth,
                      refresh =  0.1 * niter)
mod.stan.bi.nonshared <- rstan::read_stan_csv(fit$output_files())
params_to_extract <- c("beta", "sigma", "rho")
re_params <- c("v_1", "v_2", "u_1", "u_2")

# summaries
mod.stan.summary.bi.nonshared <- summary(mod.stan.bi.nonshared,
                                         pars = params_to_extract,
                                         probs = c(0.1, 0.5, 0.9))
kable(mod.stan.summary.bi.nonshared$summary, format = "markdown", digits = 3)

stan_diags <- data.frame(pct_divergent = get_num_divergent(mod.stan.bi.nonshared)/(niter * nchains * prop_warmup),
                         pct_max_tree_exceeded = get_num_max_treedepth(mod.stan.bi.nonshared)/(niter * prop_warmup * nchains),
                         pct_bmfi_low_chains = sum(is.finite(get_low_bfmi_chains(mod.stan.bi.nonshared)))/nchains)
print("Stan diagnostics")
stan_diags

# save results
res_list <- list(data = datlist,
                 priors = list(rho_beta_1 = rho_beta_1,
                               rho_beta_2 = rho_beta_2,
                               sigma_normal_sd = sigma_normal_sd),
                 stan_file = stan_file,
                 mod = list(mod.stan.bi.nonshared))

write_rds(res_list, 
          file = "../../../Dropbox/dissertation_2/survey-csmf/results/ken2014-hazwaz/ken2014-hazwaz-stage-2-bivariate-nonshared.rds")

# sample_pars <- c("beta", "sigma", "rho")
# stan_trace(mod.stan.bi.nonshared, pars = sample_pars)
# 
# pairs_pars <- c("beta", "log_sigma", "rho")
# pairs(mod.stan.bi.nonshared, pars = pairs_pars)
# 
# mcmc_nuts_energy(nuts_params(mod.stan.bi.nonshared))

# extract REs and posterior widths
res.re <- summary(mod.stan.bi.nonshared, pars = re_params, probs = c(0.1, 0.5, 0.9))$summary

haz.re.iid <- res.re[grep("v_1", rownames(res.re)), c("50%")]
waz.re.iid <- res.re[grep("v_2", rownames(res.re)), c("50%")]

haz.re.icar <- res.re[grep("u_1", rownames(res.re)), c("50%")]
waz.re.icar <- res.re[grep("u_2", rownames(res.re)), c("50%")]

res.pred <- summary(mod.stan.bi.nonshared, pars = "preds", probs = c(0.1, 0.5, 0.9))$summary[, c("10%", "50%", "90%")]

haz.pred <- res.pred[grep("preds\\[.*,1", rownames(res.pred)), ]
waz.pred <- res.pred[grep("preds\\[.*,2", rownames(res.pred)), ]

# create a data.frame with these results
res.df.bi.nonshared <- as.data.frame(cbind(results$admin1.name,
                                           haz.re.iid, waz.re.iid,
                                           haz.re.icar, waz.re.icar,
                                           haz.pred[, "50%"], haz.pred[, "90%"] - haz.pred[, "10%"], 
                                           waz.pred[, "50%"], waz.pred[, "90%"] - waz.pred[, "10%"]))

names(res.df.bi.nonshared) <- c("admin1.name", 
                                "haz.iid.50", "waz.iid.50",
                                "haz.icar.50", "waz.icar.50",
                                "haz.pred.50", "haz.pred.width80",
                                "waz.pred.50", "waz.pred.width80")
rownames(res.df.bi.nonshared) <- NULL

for (vv in 1:ncol(res.df.bi.nonshared)) {
    if (names(res.df.bi.nonshared)[vv] == "admin1.name") next
    res.df.bi.nonshared[, vv] <- as.numeric(res.df.bi.nonshared[, vv])
}


## ----bi-nonshared-plots, fig.height = 6, fig.width = 6, cache = FALSE------------------------
n_cats <- 15
plot.palette <- viridis(n_cats)

# plot REs (spatial and nonspatial)
outcomes <- c("haz.iid.50", "waz.iid.50",
              "haz.icar.50", "waz.icar.50")
outcomes.names <- c("HAZ IID RE median", "WAZ IID RE median",
                    "HAZ ICAR RE median", "WAZ ICAR RE median")

plot_maps(poly = poly.adm1, outcome.df = res.df.bi.nonshared, 
          names.reg.poly = "NAME_1", names.reg.outcomes = "admin1.name",
          outcomes = outcomes, outcomes.names = outcomes.names, 
          plot.palette = plot.palette)


## ----bi-nonshared-plots-2, fig.height = 6, fig.width = 6, cache = FALSE----------------------
# plot predictions
outcomes <- c("haz.pred.50", "haz.pred.width80",
              "waz.pred.50", "waz.pred.width80")
outcomes.names <- c("HAZ pred median", "HAZ pred 80% UI width",
                    "WAZ pred median", "WAZ pred 80% UI width")

plot_maps(poly = poly.adm1, outcome.df = res.df.bi.nonshared, 
          names.reg.poly = "NAME_1", names.reg.outcomes = "admin1.name",
          outcomes = outcomes, outcomes.names = outcomes.names, 
          plot.palette = plot.palette)
