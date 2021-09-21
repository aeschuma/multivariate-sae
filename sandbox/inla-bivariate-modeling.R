# Austin Schumacher
# 9/21/2021
# Trying INLA modeling for shared component and nonshared bivariate normal model

# Preamble ####
rm(list = ls())

library(SUMMER); library(tidyverse); library(spdep); library(geosphere); library(haven);
library(knitr); library(kableExtra); library(magrittr); library(rstan); library(cmdstanr);
library(svyVGAM); library(mvtnorm); library(rgdal); library(bayesplot); library(INLA); 
library(viridis); library(classInt);

# Functions ####

# first stage bivariate model
# In: data: data frame with the following variables
#           - "cluster"     "region"      "strata"      "weights"
#           - "admin1"      "admin1.char"       "admin1.name"
#           - the 2 outcomes
#     outcomes: character vector of length 2 with names of variables in `data` with the 2 outcomes
# Out: a list with 5 elements
#       - results: data frame with 
            # admin1
            # admin1.name
            # admin1.char
            # `outcomes[1]`.mean
            # `outcomes[2]`.mean
            # `outcomes[1]`.se
            # `outcomes[2]`.se
            # corr
#       - results.long: long version (by outcome) of the results df for use with INLA
#       - V.array: 2x2 fixed cov matrices as an array
#       - V.list: 2x2 fixed cov matrices as a list
#       - Vdes_prec: block diagonal precision matrix with the 2x2 fixed cov matrices on the diagonal

bivariate_first_stage_model <- function(data, outcomes) {
    # testing
    # data <- dat
    # outcomes <- c("HAZ", "WAZ")
    
    if (length(outcomes) != 2) stop("Need exactly 2 outcomes")
    
    my.svydesign <- survey::svydesign(ids = ~ cluster,
                                      strata = ~ strata, nest = T, weights = ~weights,
                                      data = data)
    
    n_regions <- length(unique(data$admin1.char))
    admin1v <- unique(data$admin1.char)
    results <- data.frame(admin1 = unique(data$admin1),
                          admin1.name = unique(data$admin1.name),
                          admin1.char = admin1v,
                          mean1 = rep(NA, n_regions),
                          mean2 = rep(NA, n_regions),
                          se1 = rep(NA, n_regions),
                          se2 = rep(NA, n_regions),
                          corr = rep(NA, n_regions))
    V.list <- vector(mode = "list", length = n_regions)
    V.array <- array(NA, dim = c(nrow(results), 2, 2))
    names(V.list) <- admin1v
    diag.list <- vector(mode = "list", length = length(V.list))
    
    for(i in 1:n_regions) {
        admin1.tmp <- admin1v[i]
        
        tmp <- subset(my.svydesign, admin1.char == admin1.tmp)
        fmla <- as.formula(paste0("~ ", outcomes[1], " + ", outcomes[2]))
        means.svymean <- svymean(fmla, tmp)
        
        index.tmp <- results$admin1.char == admin1.tmp
        
        results$mean1[index.tmp] <- means.svymean[[outcomes[1]]]
        results$mean2[index.tmp] <- means.svymean[[outcomes[2]]]
        
        V.tmp <- vcov(means.svymean)
        V.list[[admin1.tmp]] <- V.tmp
        V.array[i, , ] <- V.tmp
        
        results$se1[index.tmp] <- V.tmp[1, 1]^0.5
        results$se2[index.tmp] <- V.tmp[2, 2]^0.5
        
        D <- diag(sqrt(diag(V.tmp)))
        DInv <- solve(D)
        corr.tmp <- DInv %*% V.tmp %*% DInv
        
        results$corr[index.tmp] <- corr.tmp[1, 2]
        
        # create block diagonal matrix from list of fixed covariances
        diag.list[[i]] <- solve(V.tmp)
    }
    
    # reformat data
    results.long <- results %>% select(admin1, admin1.name, admin1.char,
                                       mean1, mean2) %>%
        pivot_longer(cols = c(mean1, mean2),
                     names_to = "outcome",
                     names_prefix = "mean",
                     values_to = "value")
    results.long$admin1.1 <- ifelse(results.long$outcome == 1, results.long$admin1, NA)
    results.long$admin1.2 <- ifelse(results.long$outcome == 2, results.long$admin1, NA)
    results.long$outcome <- ifelse(results.long$outcome == 1, outcomes[1], outcomes[2])
    results.long$obs <- 1:nrow(results.long)
    
    # update naems for results
    names(results) <- c("admin1", "admin1.name", "admin1.char", 
                        paste0(outcomes[1], ".mean"), paste0(outcomes[2], ".mean"), 
                        paste0(outcomes[1], ".se"), paste0(outcomes[2], ".se"),      
                        "corr")
    
    Vdes_prec <- bdiag(diag.list)
    
    return(list(results = results,
                results.long = results.long,
                V.array = V.array,
                V.list = V.list,
                Vdes_prec = Vdes_prec))
}

bivariate_second_stage_model_inla <- function(data, formula, 
                                              fe.prior, bym2.prior, cov.prior, 
                                              lambda.prior = NULL,
                                              scale) {

    
    # fit model
    smooth.direct.bi <- inla(formula,
                             data = res.list$results.long,
                             family = "gaussian",
                             control.fixed = fe.prior,
                             control.predictor=list(compute=TRUE),
                             control.compute=list(config = TRUE, dic = TRUE, cpo = TRUE, waic = TRUE),
                             control.family = list(hyper = list(prec = list(initial = log(1), fixed=TRUE))),
                             scale = scale)
    
    return(smooth.direct.bi)
}
# Load data ####
# loads: "admin1.mat", "bivariate_first_stage_model", "dat", "node.info", "poly.adm1"  
load("/Users/austin/Dropbox/dissertation_2/survey-csmf/data/ken_dhs2014/data/haz-waz-kenDHS2014.rda")

# Stage 1 direct: univariate and bivariate ####
res.list <- bivariate_first_stage_model(data = dat, outcomes = c("HAZ", "WAZ"))

# Stage 2 smoothing with INLA: nonshared ####

## model formula
m.form <- value ~ -1 + outcome + 
    f(admin1.1, model = 'bym2',
      graph = admin1.mat, 
      scale.model = T, 
      constr = T,
      hyper = bym2.prior) +
    f(admin1.2, model = 'bym2',
      graph = admin1.mat, 
      scale.model = T, 
      constr = T,
      hyper = bym2.prior) +
    f(obs,  model='generic0', Cmatrix = res.list$Vdes_prec,
      hyper = cov.prior)

# priors
fe_prec_prior <- list(prec.intercept = 0,
                      prec = 0)
bym2_prior <- list(phi=list(prior="pc", param=c(0.5, 0.25)), 
                   prec=list(prior="pc.prec", param=c(0.5/0.31, 0.01)))
cov_prec_prior <- list(prec = list(prior = "loggamma", param = c(1000000, 100000)))

# scale for fixing covariance
scale <- 10000000

# fit model
mod.inla.nonshared <- bivariate_second_stage_model_inla(data = res.list$results.long, 
                                                        formula = m.form, 
                                                        fe.prior = fe_prec_prior, 
                                                        bym2.prior = bym2_prior, 
                                                        cov.prior = cov_prec_prior,
                                                        scale = scale)
mod.inla.nonshared$summary.fixed
mod.inla.nonshared$summary.hyperpar
