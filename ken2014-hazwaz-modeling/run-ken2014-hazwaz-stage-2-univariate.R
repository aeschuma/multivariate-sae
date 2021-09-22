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
n_regions <- nrow(poly.adm1)

## ----uni-mod---------------------------------------------------------------------------------
## HAZ

# hyperpars
rho_beta_1 <- 1
rho_beta_2 <- 1
sigma_normal_sd <- 1

stan_dat_HAZ <- list(R = n_regions,
                     Sigma = results$seHAZ.bi,
                     y = results$meanHAZ.bi,
                     N_edges = node.info$n_edges,
                     node1 = node.info$node1,
                     node2 = node.info$node2,
                     scaling_factor = scaling_factor,
                     rho_beta_1 = rho_beta_1,
                     rho_beta_2 = rho_beta_2,
                     sigma_normal_sd = sigma_normal_sd)


stan_file <- "stan-models/oneD-bym2.stan"

## stan options
niter <- 10000
nchains <- 3
prop_warmup <- 0.5
max_treedepth <- 400
adapt_delta <- 0.95
nthin <- 1
options(mc.cores = parallel::detectCores())

cmd_mod <- cmdstan_model(stan_file = stan_file)
fit <- suppressMessages(cmd_mod$sample(data = stan_dat_HAZ,
                                       iter_warmup = niter*prop_warmup, 
                                       iter_sampling = niter*(1-prop_warmup),
                                       chains = nchains, thin = nthin,
                                       adapt_delta = adapt_delta, max_treedepth = max_treedepth,
                                       refresh = 0.1 * niter))
mod.stan.haz.uni <- rstan::read_stan_csv(fit$output_files())
params_to_extract <- c("beta", "sigma", "rho")
re_params <- c("v", "u")

# summaries
mod.stan.summary.haz.uni <- summary(mod.stan.haz.uni,
                                    pars = params_to_extract,
                                    probs = c(0.1, 0.5, 0.9))
kable(mod.stan.summary.haz.uni$summary, format = "markdown", digits = 3)

stan_diags <- data.frame(pct_divergent = get_num_divergent(mod.stan.haz.uni)/(niter * nchains * prop_warmup),
                         pct_max_tree_exceeded = get_num_max_treedepth(mod.stan.haz.uni)/(niter * prop_warmup * nchains),
                         pct_bmfi_low_chains = sum(is.finite(get_low_bfmi_chains(mod.stan.haz.uni)))/nchains)
print("Stan diagnostics")
stan_diags

# sample_pars <- c("beta", "sigma", "rho")
# stan_trace(mod.stan.haz.uni, pars = sample_pars)
# pairs(mod.stan.haz.uni, pars = sample_pars)
# mcmc_nuts_energy(nuts_params(mod.stan.haz.uni))

## WAZ

## hyperpars
rho_beta_1 <- 1
rho_beta_2 <- 1
sigma_normal_sd <- 1

stan_dat_WAZ <- list(R = n_regions,
                     Sigma = results$seWAZ.bi,
                     y = results$meanWAZ.bi,
                     N_edges = node.info$n_edges,
                     node1 = node.info$node1,
                     node2 = node.info$node2,
                     scaling_factor = scaling_factor,
                     rho_beta_1 = rho_beta_1,
                     rho_beta_2 = rho_beta_2,
                     sigma_normal_sd = sigma_normal_sd)

## stan options
niter <- 10000
nchains <- 3
prop_warmup <- 0.5
max_treedepth <- 400
adapt_delta <- 0.95
nthin <- 1
options(mc.cores = parallel::detectCores())

cmd_mod <- cmdstan_model(stan_file = stan_file)
fit <- suppressMessages(cmd_mod$sample(data = stan_dat_WAZ,
                                       iter_warmup = niter*prop_warmup, 
                                       iter_sampling = niter*(1-prop_warmup),
                                       chains = nchains, thin = nthin,
                                       adapt_delta = adapt_delta, max_treedepth = max_treedepth,
                                       refresh =  0.1 * niter))
mod.stan.waz.uni <- rstan::read_stan_csv(fit$output_files())
params_to_extract <- c("beta", "sigma", "rho")
re_params <- c("v", "u")

# summaries
mod.stan.summary.waz.uni <- summary(mod.stan.waz.uni,
                                    pars = params_to_extract,
                                    probs = c(0.1, 0.5, 0.9))
kable(mod.stan.summary.waz.uni$summary, format = "markdown", digits = 3)

stan_diags <- data.frame(pct_divergent = get_num_divergent(mod.stan.waz.uni)/(niter * nchains * prop_warmup),
                         pct_max_tree_exceeded = get_num_max_treedepth(mod.stan.waz.uni)/(niter * prop_warmup * nchains),
                         pct_bmfi_low_chains = sum(is.finite(get_low_bfmi_chains(mod.stan.waz.uni)))/nchains)
print("Stan diagnostics")
stan_diags

# sample_pars <- c("beta", "sigma", "rho")
# stan_trace(mod.stan.waz.uni, pars = sample_pars)

# pairs_pars <- c("beta", "log_sigma", "rho")
# pairs(mod.stan.waz.uni, pars = pairs_pars)

# mcmc_nuts_energy(nuts_params(mod.stan.waz.uni))

# save results
res_list <- list(data = list(stan_dat_HAZ, 
                             stan_dat_WAZ),
                 priors = list(rho_beta_1 = rho_beta_1,
                               rho_beta_2 = rho_beta_2,
                               sigma_normal_sd = sigma_normal_sd),
                 stan_file = stan_file,
                 mod = list(mod.stan.haz.uni,
                            mod.stan.waz.uni))

write_rds(res_list, 
          file = "../../../Dropbox/dissertation_2/survey-csmf/results/ken2014-hazwaz/ken2014-hazwaz-stage-2-univariate.rds")

# plotting!

# extract REs and posterior widths
haz.re <- summary(mod.stan.haz.uni, pars = re_params, probs = c(0.1, 0.5, 0.9))$summary
waz.re <- summary(mod.stan.waz.uni, pars = re_params, probs = c(0.1, 0.5, 0.9))$summary

haz.re.iid <- haz.re[grep("v", rownames(haz.re)), c("50%")]
waz.re.iid <- waz.re[grep("v", rownames(waz.re)), c("50%")]

haz.re.icar <- haz.re[grep("u", rownames(haz.re)), c("50%")]
waz.re.icar <- waz.re[grep("u", rownames(waz.re)), c("50%")]

haz.pred <- summary(mod.stan.haz.uni, pars = "preds", probs = c(0.1, 0.5, 0.9))$summary[, c("10%", "50%", "90%")]
waz.pred <- summary(mod.stan.waz.uni, pars = "preds", probs = c(0.1, 0.5, 0.9))$summary[, c("10%", "50%", "90%")]

# create a data.frame with these results
res.df.uni <- as.data.frame(cbind(results$admin1.name,
                                  haz.re.iid, waz.re.iid,
                                  haz.re.icar, waz.re.icar,
                                  haz.pred[, "50%"], haz.pred[, "90%"] - haz.pred[, "10%"], 
                                  waz.pred[, "50%"], waz.pred[, "90%"] - waz.pred[, "10%"]))

names(res.df.uni) <- c("admin1.name", 
                       "haz.iid.50", "waz.iid.50",
                       "haz.icar.50", "waz.icar.50",
                       "haz.pred.50", "haz.pred.width80",
                       "waz.pred.50", "waz.pred.width80")
rownames(res.df.uni) <- NULL

for (vv in 1:ncol(res.df.uni)) {
    if (names(res.df.uni)[vv] == "admin1.name") next
    res.df.uni[, vv] <- as.numeric(res.df.uni[, vv])
}


## ----uni-plots-1, fig.height = 6, fig.width = 6, cache = FALSE-------------------------------
n_cats <- 15
plot.palette <- viridis(n_cats)

# plot REs (spatial and nonspatial)
outcomes <- c("haz.iid.50", "waz.iid.50",
              "haz.icar.50", "waz.icar.50")
outcomes.names <- c("HAZ IID RE median", "WAZ IID RE median",
                    "HAZ ICAR RE median", "WAZ ICAR RE median")

plot_maps(poly = poly.adm1, outcome.df = res.df.uni, 
          names.reg.poly = "NAME_1", names.reg.outcomes = "admin1.name",
          outcomes = outcomes, outcomes.names = outcomes.names, 
          plot.palette = plot.palette)


## ----uni-plots-2, fig.height = 6, fig.width = 6, cache = FALSE-------------------------------
# plot predictions
outcomes <- c("haz.pred.50", "haz.pred.width80",
              "waz.pred.50", "waz.pred.width80")
outcomes.names <- c("HAZ pred median", "HAZ pred 80% UI width",
                    "WAZ pred median", "WAZ pred 80% UI width")

plot_maps(poly = poly.adm1, outcome.df = res.df.uni, 
          names.reg.poly = "NAME_1", names.reg.outcomes = "admin1.name",
          outcomes = outcomes, outcomes.names = outcomes.names, 
          plot.palette = plot.palette)

# extract convolved REs
haz.conv <- summary(mod.stan.haz.uni, pars = "convolved_re", probs = c(0.5))$summary[, "50%"]
waz.conv <- summary(mod.stan.waz.uni, pars = "convolved_re", probs = c(0.5))$summary[, "50%"]

# plot convolved REs, IID, and ICARs against each other to see if there is indeed a shared component
plot(haz.re.iid ~ waz.re.iid, main = "IID REs", xlab = "WAZ", ylab = "HAZ")
abline(lm(haz.re.iid ~ waz.re.iid), col = "red")

plot(haz.re.icar ~ waz.re.icar, main = "ICAR REs", xlab = "WAZ", ylab = "HAZ")
abline(lm(haz.re.icar ~ waz.re.icar), col = "red")

plot(haz.conv ~ waz.conv, main = "Convolved BYM2 REs", xlab = "WAZ", ylab = "HAZ")
abline(lm(haz.conv ~ waz.conv), col = "red")

