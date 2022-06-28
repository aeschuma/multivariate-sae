# Austin Schumacher
# 1/21/2022
# Fit all INLA models and save results

# preamble ####
rm(list = ls())

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
fitINLA <- function(formula, data, lincombs) {
    inla(formula, data = data, 
         family = "gaussian",
         lincomb = lincombs,
         control.family = list(hyper = list(prec = list(initial = 10, fixed=T))), 
         control.compute = list(config=T),
         control.fixed = list(prec = list(default = 0.001)))
}


# read and format data ####
load("/Users/austin/Dropbox/dissertation_2/survey-csmf/data/ken_dhs2014/data/haz-waz-kenDHS2014.rda")
subsample <- FALSE
if (subsample) {
    samplesize <- 0.05
    stage_1_list <- read_rds(paste0("../../../Dropbox/dissertation_2/survey-csmf/results/ken2014-multinomial/ken2014-subsample-",gsub("\\.","_",samplesize),"-multinomial-stage-1.rds"))
} else {
    stage_1_list <- read_rds(paste0("../../../Dropbox/dissertation_2/survey-csmf/results/ken2014-multinomial/ken2014-multinomial-stage-1.rds"))
}
results <- stage_1_list[["results"]]
V.array <- stage_1_list[["V.array"]]
n_regions <- nrow(poly.adm1)

## reformat data into long form
results.long <- results %>% select(admin1, admin1.name, admin1.char,
                                   mean_beta1, mean_beta2) %>%
    pivot_longer(cols = c(mean_beta1, mean_beta2),
                 names_to = "outcome",
                 names_prefix = "mean_",
                 values_to = "value")
results.long$admin1.beta1 <- ifelse(results.long$outcome == "beta1", results.long$admin1, NA)
results.long$admin1.beta2 <- ifelse(results.long$outcome == "beta2", results.long$admin1, NA)
results.long$obs <- 1:nrow(results.long)

# create a list of the data
data <- as.list(results.long)
Nt <- length(unique(results.long$admin1))
N <- nrow(results.long)
data$weights <- rep(1, N)

## Define the index-vectors ii.1 ii.2 etc, which are the
## index's for the iid2d-model at timepoint 1, 2, ...
for(j in 1:Nt) {
    itmp <-  numeric(N)
    itmp[] <-  NA
    itmp[(j-1)*2 + 1] <-  1
    itmp[(j-1)*2 + 2] <-  2
    data <-  c(list(itmp), data)
    names(data)[1] <-  paste("ii.", j, sep="")
}

## add another index for the shared component
data$admin1.beta1.2 <- data$admin1.beta1

## add indeces for iid components in besag + iid formulation
data$admin1.beta1.iid <- data$admin1.beta1
data$admin1.beta2.iid <- data$admin1.beta2

# define model components ####

## priors ####
iid_prior <- list(prec = list(prior = "pc.prec",
                              param = c(1, 0.01)))
bym2_prior <- list(phi=list(prior="logitbeta", param=c(1, 1), initial=0.5), 
                   prec=list(prior="pc.prec", param=c(1, 0.01), initial=5))
lambda_prior <- list(beta = list(prior = 'logtnormal', param = c(0, 1)))
besag_prior <- list(prec=list(prior="pc.prec", param=c(1, 0.01), initial=5))

## Diagonal V formula additions ####

### we now have to add one 'iid2d' model for each observation pair,
### since their cov.matrix is different. we have to do this
### automatically... here I add numbers directly for simplicity

### univariate
add.univariate <- ""
for(j in 1:Nt) {
    init.prec.beta1 <-  log(1/results$se_beta1[j]^2)
    init.prec.beta2 <-  log(1/results$se_beta2[j]^2)
    
    add.univariate <-  paste(add.univariate, paste(" + 
                         f(", paste("ii.", j, sep=""), ", model=\"iid2d\", n=2,
                         hyper = list(
                         prec1 = list(
                         initial =", init.prec.beta1,", 
                         fixed = TRUE),
                         prec2 = list(
                         initial =", init.prec.beta2,", 
                         fixed = TRUE),
                         cor = list(
                         initial = log((1+", 0, ")/(1-", 0, ")), 
                         fixed = TRUE)))"))
    
}

### bivariate
add.bivariate <- ""
for(j in 1:Nt) {
    corr <-  results$corr[j]
    init.prec.beta1 <-  log(1/results$se_beta1[j]^2)
    init.prec.beta2 <-  log(1/results$se_beta2[j]^2)
    
    add.bivariate <-  paste(add.bivariate, paste(" + 
                         f(", paste("ii.", j, sep=""), ", model=\"iid2d\", n=2,
                         hyper = list(
                         prec1 = list(
                         initial =", init.prec.beta1,", 
                         fixed = TRUE),
                         prec2 = list(
                         initial =", init.prec.beta2,", 
                         fixed = TRUE),
                         cor = list(
                         initial = log((1+", corr, ")/(1-", corr, ")), 
                         fixed = TRUE)))"))
    
}

## Model formulas ####
formula.univariate.iid <-  formula(paste("value ~ -1 + outcome + 
                                          f(admin1.beta1, model = 'iid', hyper = iid_prior) + 
                                          f(admin1.beta2, model = 'iid', hyper = iid_prior)",
                                         paste(add.univariate, collapse = " ")))
formula.univariate.bym <- formula(paste("value ~ -1 + outcome + 
                                         f(admin1.beta1, model = 'bym2',
                                           graph = admin1.mat, 
                                           scale.model = T, 
                                           constr = T,
                                           hyper = bym2_prior) +
                                         f(admin1.beta2, model = 'bym2',
                                           graph = admin1.mat, 
                                           scale.model = T, 
                                           constr = T,
                                           hyper = bym2_prior)",
                                        paste(add.univariate, collapse = " ")))
formula.bivariate.nonshared.iid <-  formula(paste("value ~ -1 + outcome + 
                                                   f(admin1.beta1, model = 'iid', hyper = iid_prior) + 
                                                   f(admin1.beta2, model = 'iid', hyper = iid_prior)",
                                                  paste(add.bivariate, collapse = " ")))
formula.bivariate.nonshared.bym <-  formula(paste("value ~ -1 + outcome + 
                                                   f(admin1.beta1, model = 'bym2',
                                                     graph = admin1.mat, 
                                                     scale.model = T, 
                                                     constr = T,
                                                     hyper = bym2_prior) +
                                                   f(admin1.beta2, model = 'bym2',
                                                     graph = admin1.mat, 
                                                     scale.model = T, 
                                                     constr = T,
                                                     hyper = bym2_prior)",
                                                  paste(add.bivariate, collapse = " ")))
formula.bivariate.shared.iid <-  formula(paste("value ~ -1 + outcome + 
                                                   f(admin1.beta1, model = 'iid', hyper = iid_prior) + 
                                                   f(admin1.beta2, model = 'iid', hyper = iid_prior) +
                                                   f(admin1.beta1.2, copy = \"admin1.beta2\", 
                                                     fixed = FALSE, hyper = lambda_prior)",
                                               paste(add.bivariate, collapse = " ")))
formula.bivariate.shared.bym <-  formula(paste("value ~ -1 + outcome + 
                                                f(admin1.beta1, model = 'bym2',
                                                  graph = admin1.mat, 
                                                  scale.model = T, 
                                                  constr = T,
                                                  hyper = bym2_prior) +
                                                f(admin1.beta2, model = 'bym2',
                                                  graph = admin1.mat, 
                                                  scale.model = T, 
                                                  constr = T,
                                                  hyper = bym2_prior) +
                                                f(admin1.beta1.2, copy = \"admin1.beta2\", 
                                                  fixed = FALSE, hyper = lambda_prior)",
                                               paste(add.bivariate, collapse = " ")))
formula.bivariate.shared.bym.alt <-  formula(paste("value ~ -1 + outcome + 
                                                   f(admin1.beta1, model = 'besag',
                                                     graph = admin1.mat, 
                                                     scale.model = T, 
                                                     constr = T,
                                                     hyper = besag_prior) +
                                                   f(admin1.beta1.iid, model = 'iid', hyper = iid_prior) +
                                                   f(admin1.beta2, model = 'besag',
                                                     graph = admin1.mat, 
                                                     scale.model = T, 
                                                     constr = T,
                                                     hyper = besag_prior) +
                                                   f(admin1.beta2.iid, model = 'iid', hyper = iid_prior) +
                                                   f(admin1.beta1.2, copy = \"admin1.beta2\", 
                                                     fixed = FALSE, hyper = lambda_prior)",
                                                   paste(add.bivariate, collapse = " ")))

## linear combinations ####

# linear combination of preds without 2x2 REs
lc.vec.fe <- c(rep(1, n_regions))

diag.na.mat <- matrix(NA, nrow = n_regions, ncol = n_regions)
diag(diag.na.mat) <- 1
lc.vec.admin1.re <- diag.na.mat

## beta2
lc.all.2 <- inla.make.lincombs(outcomebeta2 = lc.vec.fe, 
                               admin1.beta2 = lc.vec.admin1.re)
names(lc.all.2) <- gsub("lc", "oucome2.reg.", names(lc.all.2))

## beta2 if besag + iid
lc.all.2.alt <- inla.make.lincombs(outcomebeta2 = lc.vec.fe, 
                                   admin1.beta2 = lc.vec.admin1.re,
                                   admin1.beta2.iid = lc.vec.admin1.re)
names(lc.all.2.alt) <- gsub("lc", "outcome2.reg.", names(lc.all.2.alt))

## beta1 nonshared
lc.all.1.nonshared <- inla.make.lincombs(outcomebeta1 = lc.vec.fe, 
                                         admin1.beta1 = lc.vec.admin1.re)
names(lc.all.1.nonshared) <- gsub("lc", "outcome1.reg.", names(lc.all.1.nonshared))

## beta1 shared
lc.all.1.shared <- inla.make.lincombs(outcomebeta1 = lc.vec.fe, 
                                      admin1.beta1 = lc.vec.admin1.re,
                                      admin1.beta1.2 = lc.vec.admin1.re)
names(lc.all.1.shared) <- gsub("lc", "outcome1.reg.", names(lc.all.1.shared))

## beta1 shared if besag + iid
lc.all.1.shared.alt <- inla.make.lincombs(outcomebeta1 = lc.vec.fe, 
                                          admin1.beta1 = lc.vec.admin1.re,
                                          admin1.beta1.iid = lc.vec.admin1.re,
                                          admin1.beta1.2 = lc.vec.admin1.re)
names(lc.all.1.shared.alt) <- gsub("lc", "outcome1.reg.", names(lc.all.1.shared.alt))

# model info
model_names <- c("Univariate IID",
                 "Univariate BYM",
                 "Bivariate nonshared IID", 
                 "Bivariate nonshared BYM",
                 "Bivariate shared IID", 
                 "Bivariate shared BYM")
n_models <- length(model_names)

# Cross validation ####
n_samples <- 1000
cv_res <- tibble(model = rep(model_names, n_regions),
                 region = rep(1:n_regions, each = n_models),
                 cpo = NA,
                 pit_beta1 = NA,
                 pit_beta2 = NA)

for (r in 1:n_regions) {
    
    # region index
    idx <- ((r)*2-1):(r*2)
    
    ## Hold out data ####
    tmp <- data
    tmp$value[idx] <- NA
    y_true <- data$value[idx]
    
    ## Fit models ####
    mod_list <- vector(mode = "list", length = n_models)
    names(mod_list) <- model_names
    
    mod_list$`Univariate IID` <- fitINLA(formula = formula.univariate.iid, 
                                         data = tmp, 
                                         lincombs = c(lc.all.1.nonshared, lc.all.2))
    mod_list$`Univariate BYM` <- fitINLA(formula = formula.univariate.bym, 
                                         data = data, 
                                         lincombs = c(lc.all.1.nonshared, lc.all.2))
    mod_list$`Bivariate nonshared IID` <- fitINLA(formula = formula.bivariate.nonshared.iid, 
                                                  data = data, 
                                                  lincombs = c(lc.all.1.nonshared, lc.all.2))
    mod_list$`Bivariate nonshared BYM` <- fitINLA(formula = formula.bivariate.nonshared.bym, 
                                                  data = data, 
                                                  lincombs = c(lc.all.1.nonshared, lc.all.2))
    mod_list$`Bivariate shared IID` <- fitINLA(formula = formula.bivariate.shared.iid, 
                                               data = data, 
                                               lincombs = c(lc.all.1.shared, lc.all.2))
    mod_list$`Bivariate shared BYM` <- fitINLA(formula = formula.bivariate.shared.bym, 
                                               data = data, 
                                               lincombs = c(lc.all.1.shared, lc.all.2))
    starttime <- Sys.time()
    # simulating posterior dist and do CV calculations for each model
    for (mm in 1:length(mod_list)) {
        # sample from the posterior
        samp <- inla.posterior.sample(n = n_samples, mod_list[[mm]])
        
        if (names(mod_list)[mm] %in% c("Univariate IID", "Univariate BYM", 
                                       "Bivariate nonshared IID", "Bivariate nonshared BYM")) {
            # fixed effects
            beta1_int_idx <- rownames(samp[[1]]$latent) %>% str_detect("outcomebeta1:1") %>% which()
            beta2_int_idx <- rownames(samp[[1]]$latent) %>% str_detect("outcomebeta2:1") %>% which()
            beta1_int_mat <- matrix(0, nrow = length(beta1_int_idx), ncol = n_samples)
            beta2_int_mat <- matrix(0, nrow = length(beta2_int_idx), ncol = n_samples)
            
            # random effects
            beta1_re_idx <- which(rownames(samp[[1]]$latent) == paste0("admin1.beta1:",r))
            beta2_re_idx <- which(rownames(samp[[1]]$latent) == paste0("admin1.beta2:",r))
            beta1_re_mat <- matrix(0, nrow = length(beta1_re_idx), ncol = n_samples)
            beta2_re_mat <- matrix(0, nrow = length(beta2_re_idx), ncol = n_samples)
            
            # fill matrices for each sample
            for (i in 1:n_samples) {
                beta1_int_mat[,i] <- samp[[i]]$latent[beta1_int_idx]
                beta2_int_mat[,i] <- samp[[i]]$latent[beta2_int_idx]
                beta1_re_mat[,i] <- samp[[i]]$latent[beta1_re_idx]
                beta2_re_mat[,i] <- samp[[i]]$latent[beta2_re_idx]
            }
            
            # set up fitted mat
            beta1_pred_mat <- beta1_int_mat + beta1_re_mat
            beta2_pred_mat <- beta2_int_mat + beta2_re_mat
            posterior_pred <- rbind(beta1_pred_mat, beta2_pred_mat)
            
        } else if (names(mod_list)[mm] %in% c("Bivariate shared IID", "Bivariate shared BYM")) {
            # lambda (it's a hyperparameter in our INLA function)
            hyperpar_lambda_idx <- names(samp[[1]]$hyperpar) %>% str_detect("Beta") %>% which()
            lambda_mat <- matrix(0, nrow = length(hyperpar_lambda_idx), ncol = n_samples)
            
            # fixed effects
            beta1_int_idx <- rownames(samp[[1]]$latent) %>% str_detect("outcomebeta1:1") %>% which()
            beta2_int_idx <- rownames(samp[[1]]$latent) %>% str_detect("outcomebeta2:1") %>% which()
            beta1_int_mat <- matrix(0, nrow = length(beta1_int_idx), ncol = n_samples)
            beta2_int_mat <- matrix(0, nrow = length(beta2_int_idx), ncol = n_samples)
            
            # random effects
            beta1_re_idx <- which(rownames(samp[[1]]$latent) == paste0("admin1.beta1:",r))
            beta2_re_idx <- which(rownames(samp[[1]]$latent) == paste0("admin1.beta2:",r))
            beta1_re_mat <- matrix(0, nrow = length(beta1_re_idx), ncol = n_samples)
            beta2_re_mat <- matrix(0, nrow = length(beta2_re_idx), ncol = n_samples)
            
            # fill in matrices for each sample
            for (i in 1:n_samples) {
                lambda_mat[,i] <- samp[[i]]$hyperpar[hyperpar_lambda_idx]
                beta1_int_mat[,i] <- samp[[i]]$latent[beta1_int_idx]
                beta2_int_mat[,i] <- samp[[i]]$latent[beta2_int_idx]
                beta1_re_mat[,i] <- samp[[i]]$latent[beta1_re_idx]
                beta2_re_mat[,i] <- samp[[i]]$latent[beta2_re_idx]
            }
            
            # set up fitted mat
            beta1_pred_mat <- beta1_int_mat + beta1_re_mat + (lambda_mat * beta2_re_mat)
            beta2_pred_mat <- beta2_int_mat + beta2_re_mat
            posterior_pred <- rbind(beta1_pred_mat, beta2_pred_mat)
            
        } 
        
        # calculate CV results
        y_lik <- rep(NA, n_samples)
        beta1_prob <- rep(NA, n_samples)
        beta2_prob <- rep(NA, n_samples)
        
        for (i in 1:n_samples) {
            y_lik[i] <- dmvnorm(x = y_true, mean = posterior_pred[,i], sigma = V.array[r,,])
            beta1_prob[i] <- pnorm(y_true[1], mean = posterior_pred[1,i], sd = sqrt(V.array[r,1,1]))
            beta2_prob[i] <- pnorm(y_true[2], mean = posterior_pred[2,i], sd = sqrt(V.array[r,2,2]))
        }
        cv_res[cv_res$model == model_names[mm] & cv_res$region == r, "cpo"] <- mean(y_lik)
        cv_res[cv_res$model == model_names[mm] & cv_res$region == r, "pit_beta1"] <- mean(beta1_prob)
        cv_res[cv_res$model == model_names[mm] & cv_res$region == r, "pit_beta2"] <- mean(beta2_prob)
    }
    
    endtime <- Sys.time()
    totaltime <- endtime - starttime
}

# save results
write_rds(cv_res, file = "../../../Dropbox/dissertation_2/survey-csmf/results/cv/cv_results-arealevel-multinomial.rds")

# format
cv_res %<>% mutate(model_factor = factor(model, levels = model_names))

# make summary tables and plots
pit_beta1_plot <- ggplot(cv_res, aes(x = pit_beta1)) +
    geom_histogram() + 
    facet_wrap(~ model_factor, ncol = 2) +
    theme_light()

pit_beta2_plot <- ggplot(cv_res, aes(x = pit_beta2)) +
    geom_histogram() + 
    facet_wrap(~ model_factor, ncol = 2) +
    theme_light()

ggarrange(plotlist = list(pit_beta1_plot, pit_beta2_plot), nrow = 1, labels = c("beta1", "beta2"))

logCPOres <- cv_res %>% mutate(logcpo = log(cpo)) %>% 
    group_by(model_factor) %>% 
    summarise(logCPO_num = round(-1 * sum(logcpo), 2))

logCPOres$logCPO <- ifelse(logCPOres$logCPO_num == min(logCPOres$logCPO_num), 
                           paste0("\\textbf{", logCPOres$logCPO_num, "}"),
                           as.character(logCPOres$logCPO_num))

logCPOres %>% select(model_factor, logCPO) %>%
    kable(format = "markdown", caption = "Bivariate $-\\sum\\log(CPO)$ for each model. Bold indicates the best performing model",
          col.names = c("Model", "$-\\sum\\log(CPO)$"))
