# Austin Schumacher
# 1/21/2022
# Fit all INLA models and save results

# preamble ####
rm(list = ls())

## set the root depending on operating system
root <- ifelse(Sys.info()[1]=="Darwin","~/",
               ifelse(Sys.info()[1]=="Windows","P:/",
                      ifelse(Sys.info()[1]=="Linux","/home/users/aeschuma/",
                             stop("Unknown operating system"))))

if (root ==  "/home/users/aeschuma/") {
    setwd(paste0(root, "Desktop/survey_csmf/cross-validation"))
}

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
load(paste0(root,"Dropbox/dissertation_2/survey-csmf/data/ken_dhs2014/data/haz-waz-kenDHS2014.rda"))
subsample <- FALSE
if (subsample) {
    samplesize <- 0.05
    stage_1_list <- read_rds(paste0("../../../Dropbox/dissertation_2/survey-csmf/results/ken2014-hazwaz/ken2014-subsample-",gsub("\\.","_",samplesize),"-hazwaz-stage-1.rds"))
} else {
    stage_1_list <- read_rds(paste0("../../../Dropbox/dissertation_2/survey-csmf/results/ken2014-hazwaz/ken2014-hazwaz-stage-1.rds"))
}
results <- stage_1_list[["results"]]
V.array <- stage_1_list[["V.array"]]
n_regions <- nrow(poly.adm1)

## reformat data into long form
results.long <- results %>% select(admin1, admin1.name, admin1.char,
                                   meanHAZ.bi, meanWAZ.bi) %>%
    pivot_longer(cols = c(meanHAZ.bi, meanWAZ.bi),
                 names_to = "outcome",
                 names_prefix = "mean",
                 values_to = "value")
results.long$outcome <- ifelse(results.long$outcome == "HAZ.bi", "HAZ", "WAZ")
results.long$admin1.haz <- ifelse(results.long$outcome == "HAZ", results.long$admin1, NA)
results.long$admin1.waz <- ifelse(results.long$outcome == "WAZ", results.long$admin1, NA)
results.long$obs <- 1:nrow(results.long)

# create a list of the data
data <- as.list(results.long)
Nt <- length(unique(results.long$admin1))
N <- nrow(results.long)
data$weights <- rep(1, N)

## add another index for the shared component
data$admin1.haz.2 <- data$admin1.haz

## add indeces for iid components in besag + iid formulation
data$admin1.haz.iid <- data$admin1.haz
data$admin1.waz.iid <- data$admin1.waz

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

# define model components ####

## priors ####
iid_prior <- list(prec = list(prior = "pc.prec",
                              param = c(1, 0.01)))
bym2_prior <- list(phi=list(prior="logitbeta", param=c(1, 1), initial=0.5), 
                   prec=list(prior="pc.prec", param=c(1, 0.01), initial=5))
lambda_prior <- list(beta = list(prior = 'logtnormal', param = c(0, 1)))

## Diagonal V formula additions ####

### we now have to add one 'iid2d' model for each observation pair,
### since their cov.matrix is different. we have to do this
### automatically... here I add numbers directly for simplicity

### univariate
add.univariate <- ""
for(j in 1:Nt) {
    init.prec.haz <-  log(1/results$seHAZ.bi[j]^2)
    init.prec.waz <-  log(1/results$seWAZ.bi[j]^2)
    
    add.univariate <-  paste(add.univariate, paste(" + 
                         f(", paste("ii.", j, sep=""), ", model=\"iid2d\", n=2,
                         hyper = list(
                         prec1 = list(
                         initial =", init.prec.haz,", 
                         fixed = TRUE),
                         prec2 = list(
                         initial =", init.prec.waz,", 
                         fixed = TRUE),
                         cor = list(
                         initial = log((1+", 0, ")/(1-", 0, ")), 
                         fixed = TRUE)))"))
    
}

### bivariate
add.bivariate <- ""
for(j in 1:Nt) {
    corr <-  results$corr.bi[j]
    init.prec.haz <-  log(1/results$seHAZ.bi[j]^2)
    init.prec.waz <-  log(1/results$seWAZ.bi[j]^2)
    
    add.bivariate <-  paste(add.bivariate, paste(" + 
                         f(", paste("ii.", j, sep=""), ", model=\"iid2d\", n=2,
                         hyper = list(
                         prec1 = list(
                         initial =", init.prec.haz,", 
                         fixed = TRUE),
                         prec2 = list(
                         initial =", init.prec.waz,", 
                         fixed = TRUE),
                         cor = list(
                         initial = log((1+", corr, ")/(1-", corr, ")), 
                         fixed = TRUE)))"))
    
}

## Model formulas ####
formula.univariate.iid <-  formula(paste("value ~ -1 + outcome + 
                                          f(admin1.haz, model = 'iid', hyper = iid_prior) + 
                                          f(admin1.waz, model = 'iid', hyper = iid_prior)",
                                         paste(add.univariate, collapse = " ")))
formula.univariate.bym <- formula(paste("value ~ -1 + outcome + 
                                         f(admin1.haz, model = 'bym2',
                                           graph = admin1.mat, 
                                           scale.model = T, 
                                           constr = T,
                                           hyper = bym2_prior) +
                                         f(admin1.waz, model = 'bym2',
                                           graph = admin1.mat, 
                                           scale.model = T, 
                                           constr = T,
                                           hyper = bym2_prior)",
                                        paste(add.univariate, collapse = " ")))
formula.bivariate.nonshared.iid <-  formula(paste("value ~ -1 + outcome + 
                                                   f(admin1.haz, model = 'iid', hyper = iid_prior) + 
                                                   f(admin1.waz, model = 'iid', hyper = iid_prior)",
                                                  paste(add.bivariate, collapse = " ")))
formula.bivariate.nonshared.bym <-  formula(paste("value ~ -1 + outcome + 
                                                   f(admin1.haz, model = 'bym2',
                                                     graph = admin1.mat, 
                                                     scale.model = T, 
                                                     constr = T,
                                                     hyper = bym2_prior) +
                                                   f(admin1.waz, model = 'bym2',
                                                     graph = admin1.mat, 
                                                     scale.model = T, 
                                                     constr = T,
                                                     hyper = bym2_prior)",
                                                  paste(add.bivariate, collapse = " ")))
formula.bivariate.shared.iid <-  formula(paste("value ~ -1 + outcome + 
                                                   f(admin1.haz, model = 'iid', hyper = iid_prior) + 
                                                   f(admin1.waz, model = 'iid', hyper = iid_prior) +
                                                   f(admin1.haz.2, copy = \"admin1.waz\", 
                                                     fixed = FALSE, hyper = lambda_prior)",
                                               paste(add.bivariate, collapse = " ")))
formula.bivariate.shared.bym <-  formula(paste("value ~ -1 + outcome + 
                                                f(admin1.haz, model = 'bym2',
                                                  graph = admin1.mat, 
                                                  scale.model = T, 
                                                  constr = T,
                                                  hyper = bym2_prior) +
                                                f(admin1.waz, model = 'bym2',
                                                  graph = admin1.mat, 
                                                  scale.model = T, 
                                                  constr = T,
                                                  hyper = bym2_prior) +
                                                f(admin1.haz.2, copy = \"admin1.waz\", 
                                                  fixed = FALSE, hyper = lambda_prior)",
                                               paste(add.bivariate, collapse = " ")))

## linear combinations ####

# linear combination of preds without 2x2 REs
lc.vec.fe <- c(rep(1, n_regions))

diag.na.mat <- matrix(NA, nrow = n_regions, ncol = n_regions)
diag(diag.na.mat) <- 1
lc.vec.admin1.re <- diag.na.mat

## WAZ
lc.all.waz <- inla.make.lincombs(outcomeWAZ = lc.vec.fe, 
                                 admin1.waz = lc.vec.admin1.re)
names(lc.all.waz) <- gsub("lc", "waz.reg.", names(lc.all.waz))

## WAZ if besag + iid
lc.all.waz.alt <- inla.make.lincombs(outcomeWAZ = lc.vec.fe, 
                                     admin1.waz = lc.vec.admin1.re,
                                     admin1.waz.iid = lc.vec.admin1.re)
names(lc.all.waz.alt) <- gsub("lc", "waz.reg.", names(lc.all.waz.alt))

## HAZ nonshared
lc.all.haz.nonshared <- inla.make.lincombs(outcomeHAZ = lc.vec.fe, 
                                           admin1.haz = lc.vec.admin1.re)
names(lc.all.haz.nonshared) <- gsub("lc", "haz.reg.", names(lc.all.haz.nonshared))

## HAZ shared
lc.all.haz.shared <- inla.make.lincombs(outcomeHAZ = lc.vec.fe, 
                                        admin1.haz = lc.vec.admin1.re,
                                        admin1.haz.2 = lc.vec.admin1.re)
names(lc.all.haz.shared) <- gsub("lc", "haz.reg.", names(lc.all.haz.shared))

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
                 pit_haz = NA,
                 pit_waz = NA)

for (r in 1:n_regions) {
    starttime <- Sys.time()
    
    message(paste0("region ", r))
    
    # region index
    idx <- data$admin1 == r
    
    ## Hold out data ####
    tmp <- data
    tmp$value[idx] <- NA
    y_true <- data$value[idx]
    
    message("Running models...")
    
    ## Fit models ####
    mod_list <- vector(mode = "list", length = n_models)
    names(mod_list) <- model_names
    
    mod_list$`Univariate IID` <- fitINLA(formula = formula.univariate.iid, 
                                         data = tmp, 
                                         lincombs = c(lc.all.haz.nonshared, lc.all.waz))
    mod_list$`Univariate BYM` <- fitINLA(formula = formula.univariate.bym, 
                                         data = tmp, 
                                         lincombs = c(lc.all.haz.nonshared, lc.all.waz))
    mod_list$`Bivariate nonshared IID` <- fitINLA(formula = formula.bivariate.nonshared.iid, 
                                                  data = tmp, 
                                                  lincombs = c(lc.all.haz.nonshared, lc.all.waz))
    mod_list$`Bivariate nonshared BYM` <- fitINLA(formula = formula.bivariate.nonshared.bym, 
                                                  data = tmp, 
                                                  lincombs = c(lc.all.haz.nonshared, lc.all.waz))
    mod_list$`Bivariate shared IID` <- fitINLA(formula = formula.bivariate.shared.iid, 
                                               data = tmp, 
                                               lincombs = c(lc.all.haz.shared, lc.all.waz))
    mod_list$`Bivariate shared BYM` <- fitINLA(formula = formula.bivariate.shared.bym, 
                                               data = tmp, 
                                               lincombs = c(lc.all.haz.shared, lc.all.waz))
    
    message("Posterior sims + calculations...")
    
    # simulating posterior dist and do CV calculations for each model
    for (mm in 1:length(mod_list)) {
        # sample from the posterior
        samp <- inla.posterior.sample(n = n_samples, mod_list[[mm]])
        
        if (names(mod_list)[mm] %in% c("Univariate IID", "Univariate BYM", 
                                       "Bivariate nonshared IID", "Bivariate nonshared BYM")) {
            # fixed effects
            haz_int_idx <- rownames(samp[[1]]$latent) %>% str_detect("outcomeHAZ:1") %>% which()
            waz_int_idx <- rownames(samp[[1]]$latent) %>% str_detect("outcomeWAZ:1") %>% which()
            haz_int_mat <- matrix(0, nrow = length(haz_int_idx), ncol = n_samples)
            waz_int_mat <- matrix(0, nrow = length(waz_int_idx), ncol = n_samples)
            
            # random effects
            haz_re_idx <- which(rownames(samp[[1]]$latent) == paste0("admin1.haz:",r))
            waz_re_idx <- which(rownames(samp[[1]]$latent) == paste0("admin1.waz:",r))
            haz_re_mat <- matrix(0, nrow = length(haz_re_idx), ncol = n_samples)
            waz_re_mat <- matrix(0, nrow = length(waz_re_idx), ncol = n_samples)
            
            # fill matrices for each sample
            for (i in 1:n_samples) {
                haz_int_mat[,i] <- samp[[i]]$latent[haz_int_idx]
                waz_int_mat[,i] <- samp[[i]]$latent[waz_int_idx]
                haz_re_mat[,i] <- samp[[i]]$latent[haz_re_idx]
                waz_re_mat[,i] <- samp[[i]]$latent[waz_re_idx]
            }
            
            # set up fitted mat
            haz_pred_mat <- haz_int_mat + haz_re_mat
            waz_pred_mat <- waz_int_mat + waz_re_mat
            posterior_pred <- rbind(haz_pred_mat, waz_pred_mat)
            
        } else if (names(mod_list)[mm] %in% c("Bivariate shared IID", "Bivariate shared BYM")) {
            # lambda (it's a hyperparameter in our INLA function)
            hyperpar_lambda_idx <- names(samp[[1]]$hyperpar) %>% str_detect("Beta") %>% which()
            lambda_mat <- matrix(0, nrow = length(hyperpar_lambda_idx), ncol = n_samples)
            
            # fixed effects
            haz_int_idx <- rownames(samp[[1]]$latent) %>% str_detect("outcomeHAZ:1") %>% which()
            waz_int_idx <- rownames(samp[[1]]$latent) %>% str_detect("outcomeWAZ:1") %>% which()
            haz_int_mat <- matrix(0, nrow = length(haz_int_idx), ncol = n_samples)
            waz_int_mat <- matrix(0, nrow = length(waz_int_idx), ncol = n_samples)
            
            # random effects
            haz_re_idx <- which(rownames(samp[[1]]$latent) == paste0("admin1.haz:",r))
            waz_re_idx <- which(rownames(samp[[1]]$latent) == paste0("admin1.waz:",r))
            haz_re_mat <- matrix(0, nrow = length(haz_re_idx), ncol = n_samples)
            waz_re_mat <- matrix(0, nrow = length(waz_re_idx), ncol = n_samples)
            
            # fill in matrices for each sample
            for (i in 1:n_samples) {
                lambda_mat[,i] <- samp[[i]]$hyperpar[hyperpar_lambda_idx]
                haz_int_mat[,i] <- samp[[i]]$latent[haz_int_idx]
                waz_int_mat[,i] <- samp[[i]]$latent[waz_int_idx]
                haz_re_mat[,i] <- samp[[i]]$latent[haz_re_idx]
                waz_re_mat[,i] <- samp[[i]]$latent[waz_re_idx]
            }
            
            # set up fitted mat
            haz_pred_mat <- haz_int_mat + haz_re_mat + (lambda_mat * waz_re_mat)
            waz_pred_mat <- waz_int_mat + waz_re_mat
            posterior_pred <- rbind(haz_pred_mat, waz_pred_mat)
            
        } 
        
        # calculate CV results
        y_lik <- rep(NA, n_samples)
        haz_prob <- rep(NA, n_samples)
        waz_prob <- rep(NA, n_samples)
        
        for (i in 1:n_samples) {
            y_lik[i] <- dmvnorm(x = y_true, mean = posterior_pred[,i], sigma = V.array[r,,])
            haz_prob[i] <- pnorm(y_true[1], mean = posterior_pred[1,i], sd = sqrt(V.array[r,1,1]))
            waz_prob[i] <- pnorm(y_true[2], mean = posterior_pred[2,i], sd = sqrt(V.array[r,2,2]))
        }
        cv_res[cv_res$model == model_names[mm] & cv_res$region == r, "cpo"] <- mean(y_lik)
        cv_res[cv_res$model == model_names[mm] & cv_res$region == r, "pit_haz"] <- mean(haz_prob)
        cv_res[cv_res$model == model_names[mm] & cv_res$region == r, "pit_waz"] <- mean(waz_prob)
    }
    
    endtime <- Sys.time()
    totaltime <- endtime - starttime
    
    message(paste0("DONE! Time: ", totaltime, " min"))
}

# save results
write_rds(cv_res, file = "../../../Dropbox/dissertation_2/survey-csmf/results/cv/cv_results-hazwaz-arealevel.rds")

# if on Box, copy to real dropbox
if (root == "P:/") {
    write_rds(cv_res, file = "C:/Users/aeschuma/Dropbox/dissertation_2/survey-csmf/results/cv/cv_results-hazwaz-arealevel.rds")
}

# format
cv_res %<>% mutate(model_factor = factor(model, levels = model_names))

# make summary tables and plots
pit_haz_plot <- ggplot(cv_res, aes(x = pit_haz)) +
    geom_histogram(bins = 7) + 
    facet_wrap(~ model_factor, ncol = 2) +
    theme_light()

pit_waz_plot <- ggplot(cv_res, aes(x = pit_waz)) +
    geom_histogram(bins = 7) + 
    facet_wrap(~ model_factor, ncol = 2) +
    theme_light()

ggarrange(plotlist = list(pit_haz_plot, pit_waz_plot), nrow = 1, labels = c("HAZ", "WAZ"))

logCPOres <- cv_res %>% mutate(logcpo = log(cpo)) %>% 
    group_by(model_factor) %>% 
    summarise(logCPO_num = round(-1 * sum(logcpo), 2)) %>%
    arrange(desc(logCPO_num))

logCPOres$logCPO <- ifelse(logCPOres$logCPO_num == min(logCPOres$logCPO_num), 
                           paste0("\\textbf{", logCPOres$logCPO_num, "}"),
                           as.character(logCPOres$logCPO_num))

logCPOres %>% select(model_factor, logCPO) %>%
    kable(format = "markdown", caption = "Bivariate $-\\sum\\log(CPO)$ for each model. Bold indicates the best performing model",
          col.names = c("Model", "$-\\sum\\log(CPO)$"))
