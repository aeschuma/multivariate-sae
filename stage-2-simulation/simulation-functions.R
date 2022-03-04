# function:
#   sim.Q()
# purpose:
#   simulate from the precision matrix of the GMRF
# Gratuitously stolen from the rst() function in the SUMMER package
sim.Q <- function(Q){
    eigenQ <- eigen(Q)
    rankQ <- qr(Q)$rank
    sim <- as.vector(eigenQ$vectors[,1:rankQ] %*% 
                         matrix(
                             stats::rnorm(rep(1, rankQ), rep(0, rankQ), 1/sqrt(eigenQ$values[1:rankQ])),
                             ncol = 1))
    sim
}

# function: 
#   simulateData()
# purpose: 
#   simulate data to fit a STAN model on 
#   for a given set of parameters
# input:
#   node_list: or poly_file?
#   dgm: data generating mechanism, parameters in file "dgm-info.csv"
#   testing: logical, TRUE if you are just testing the function for troubleshooting
# output:
#   datlist: list of data to be used as input to various STAN models
#   params: values of parameters used in the simulation
simulateData <- function(dgm_specs, n_r, Amat, scaling_factor, seed_re, seed_lik, testing = FALSE) {
    
    # testing
    if (testing == TRUE) {
        dgm_specs <- my_dgm
        n_r <- table(dat$admin1)
        Amat <- admin1.mat
        scaling_factor <- scaling_factor
        seed_re <- 80085
        seed_lik <- 8008135
    }
    
    n_regions <- nrow(Amat)
    n_cause <- 2
    n_obs <- n_regions * n_cause
    
    # load dgm data
    my_dgm <- dgm_specs %>% unlist()
    one_precise_one_not <- my_dgm["one_precise_one_not"]
    my_dgm <- my_dgm[!(names(my_dgm) == "one_precise_one_not")]
    
    # load first stage model results
    mod_lists <- list(readRDS("~/Dropbox/dissertation_2/survey-csmf/results/ken2014-hazwaz/ken2014-hazwaz-stage-1.rds"))
    names(mod_lists) <- "First stage"
    
    # load all 2nd stage model results
    mod_lists_2 <- readRDS("~/Dropbox/dissertation_2/survey-csmf/results/ken2014-hazwaz/ken2014-hazwaz-stage-2-inla-all.rds")
    
    # combine lists
    mod_lists <- c(mod_lists, mod_lists_2)
    
    # set parameters
    all_par_names <- names(my_dgm)
    
    V_pars <- my_dgm[grep("V_", all_par_names, value = TRUE)]
    V_pars_numeric <- suppressWarnings(as.numeric(V_pars))
    if (length(unique(V_pars)) == 1 & sum(is.na(V_pars_numeric)) == length(V_pars_numeric)) { # if all parameters are from the same model
        this.V.mod <- unique(V_pars)
        V.array <- mod_lists[[this.V.mod]]$V.array
        num_corrs <- dim(V.array)[1]
        my.corrs <- rep(NA, num_corrs)
        for (j in 1:num_corrs) {
            Vtmp <- V.array[j,,]
            D <- diag(sqrt(diag(Vtmp)))
            DInv <- solve(D)
            corr.tmp <- DInv %*% Vtmp %*% DInv
            my.corrs[j] <- corr.tmp[1,2]
        }
    } else {
        # V correlation
        corr_pars <- V_pars[grepl("V_corr", names(V_pars))]
        corr_pars_num <- suppressWarnings(as.numeric(corr_pars))
        names(corr_pars_num) <- names(corr_pars)

        if (sum(is.na(corr_pars_num)) == 0) { # if correlation parameters are numeric
            
            if (corr_pars_num["V_corr_lower"] > corr_pars_num["V_corr_upper"]) {
                stop("Correlation lower bound is larger than upper bound!")
            }
            
            # correlation parameters are random uniform between bounds
            num_corrs <- n_regions
            
            # set seed for these
            set.seed(1)
            my.corrs <- runif(num_corrs, corr_pars_num["V_corr_lower"], corr_pars_num["V_corr_upper"])
            
        } else if (sum(is.na(corr_pars_num)) == length(corr_pars_num)) { # if correlation parameters are from a model
            
            if (length(unique(corr_pars)) != 1) {
                stop("If based on a model, correlation parameters must be from the SAME model")
            }
            corr_par <- unique(corr_pars)
            V.array.tmp <- mod_lists[[corr_par]]$V.array
            num_corrs <- dim(V.array.tmp)[1]
            my.corrs <- rep(NA, num_corrs)
            for (j in 1:num_corrs) {
                Vtmp <- V.array.tmp[j,,]
                D <- diag(sqrt(diag(Vtmp)))
                DInv <- solve(D)
                corr.tmp <- DInv %*% Vtmp %*% DInv
                my.corrs[j] <- corr.tmp[1,2]
            }
        } else {
            stop("Need to either have correlation lower and upper bounds both from a model, or both numeric")
        }
        
        # V diags
        vdiag_pars <- V_pars[grepl("V_diag", names(V_pars))]
        vdiag_pars_num <- suppressWarnings(as.numeric(vdiag_pars))
        names(vdiag_pars_num) <- names(vdiag_pars)
        
        if (one_precise_one_not) {
            if (vdiag_pars_num["V_diag_lower"] > vdiag_pars_num["V_diag_upper"]) {
                stop("V diagonal lower bound (precise) is larger than upper bound (imprecise)!")
            }
            
            my.diag.mat <- array(NA, dim = c(n_regions, n_cause, n_cause))
            for (jj in 1:dim(my.diag.mat)[1]) {
                my.diag.mat[jj,,] <- diag(vdiag_pars_num, 
                                          nrow = 2, ncol = 2)
            }
        } else if (sum(is.na(vdiag_pars_num)) == 0) { # if V diagonal parameters are numeric
            
            if (vdiag_pars_num["V_diag_lower"] > vdiag_pars_num["V_diag_upper"]) {
                stop("V diagonal lower bound is larger than upper bound!")
            }

            # set seed for these
            set.seed(2)
            my.diag.mat <- array(NA, dim = c(n_regions, n_cause, n_cause))
            for (jj in 1:dim(my.diag.mat)[1]) {
                my.diag.mat[jj,,] <- diag(runif(2, vdiag_pars_num["V_diag_lower"], vdiag_pars_num["V_diag_upper"]), 
                                          nrow = 2, ncol = 2)
            }
        } else if (sum(is.na(vdiag_pars_num)) == length(vdiag_pars_num)) { # if V diagonal parameters are from a model
            
            if (length(unique(vdiag_pars)) != 1) {
                stop("If based on a model, V diagonal parameters must be from the SAME model")
            }
            vdiag_par <- unique(vdiag_pars)
            V.array.tmp <- mod_lists[[vdiag_par]]$V.array
            num_vdiags <- dim(V.array.tmp)[1]
            my.diag.mat <- V.array.tmp
            for (jj in 1:num_vdiags) {
                my.diag.mat[jj,1,2] <- 0
                my.diag.mat[jj,2,1] <- 0
            }
        } else {
            stop("Need to either have V diagonal lower and upper bounds both from a model, or both numeric")
        }
        
        if (length(my.corrs) != dim(my.diag.mat)[1]) {
            stop("Number of correlations is not equal to the number of V diagonals")
        }
        
        # combine corr and V diagonal
        V.array <- array(NA, dim = dim(my.diag.mat))
        for (jj in 1:dim(my.diag.mat)[1]) {
            Omega <- matrix(c(1, my.corrs[jj], my.corrs[jj], 1), ncol = 2, nrow = 2)
            tmp.mat <- sqrt(my.diag.mat[jj,,]) %*% Omega %*% sqrt(my.diag.mat[jj,,])
            V.array[jj,,] <- tmp.mat
        }
    }
    
    mean_pars_names <- all_par_names[!(all_par_names %in% names(V_pars))]
    
    mean_pars.tmp <- my_dgm[mean_pars_names]
    mean_pars_num.tmp <- suppressWarnings(as.numeric(mean_pars.tmp))
    
    mean_pars <- rep(NA, length(mean_pars.tmp))
    
    # mean_pars_from_model <- mean_pars[which(is.na(mean_pars_num))]
    
    for (i in 1:length(mean_pars)) {
        if (!is.na(mean_pars_num.tmp[i])) {
            mean_pars[i] <- as.numeric(mean_pars.tmp[i])
        } else {
            par_mod.tmp <- mean_pars.tmp[i]
            if (mean_pars_names[i] == "lambda" & grepl("nonshared", par_mod.tmp)) {
                mean_pars[i] <- 0
            } else {
                if (names(par_mod.tmp) == "beta[1]") {
                    mean_pars[i] <- mod_lists[[par_mod.tmp]]$summary.fixed["outcomeHAZ", "0.5quant"]
                } else if (names(par_mod.tmp) == "beta[2]") {
                    mean_pars[i] <- mod_lists[[par_mod.tmp]]$summary.fixed["outcomeWAZ", "0.5quant"]
                } else if (names(par_mod.tmp) == "sigma[1]") {
                    mean_pars[i] <- mod_lists[[par_mod.tmp]]$summary.hyperpar["Precision for admin1.haz", "0.5quant"]^(-0.5)
                } else if (names(par_mod.tmp) == "sigma[2]") {
                    mean_pars[i] <- mod_lists[[par_mod.tmp]]$summary.hyperpar["Precision for admin1.waz", "0.5quant"]^(-0.5)
                } else if (names(par_mod.tmp) == "rho[1]") {
                    mean_pars[i] <- mod_lists[[par_mod.tmp]]$summary.hyperpar["Phi for admin1.haz", "0.5quant"]
                } else if (names(par_mod.tmp) == "rho[2]") {
                    mean_pars[i] <- mod_lists[[par_mod.tmp]]$summary.hyperpar["Phi for admin1.waz", "0.5quant"]
                } else if (names(par_mod.tmp) == "lambda") {
                    mean_pars[i] <- mod_lists[[par_mod.tmp]]$summary.hyperpar["Beta for admin1.haz.2", "0.5quant"]
                } 
            }
        }
    }
    names(mean_pars) <- mean_pars_names
    
    # simulate data!
    
    # set seed for random effects
    set.seed(seed_re)
    
    ## IID normal(0, 1)
    v_1 <- rnorm(n_regions, 0, 1)
    v_2 <- rnorm(n_regions, 0, 1)
    
    ## ICAR 
    if(ncol(Amat) != nrow(Amat)){
        stop("Amat does not have the same number of rows and columns.")
    }
    if(sum(Amat %in% c(0, 1)) < length(Amat)){
        stop("Amat contains values other than 0 and 1.")
    }
    Q <- Amat * -1
    diag(Q) <- 0
    diag(Q) <- -1 * apply(Q, 2, sum)

    u_1 <- sim.Q(Q)
    u_2 <- sim.Q(Q) 
    
    ## Convolved REs
    convolved_re_1 <- (sqrt(1 - mean_pars["rho[1]"]) * v_1) + (sqrt(mean_pars["rho[1]"] / scaling_factor) * u_1)
    convolved_re_2 <- (sqrt(1 - mean_pars["rho[2]"]) * v_2) + (sqrt(mean_pars["rho[2]"] / scaling_factor) * u_2)
    
    # set seed for the likelihood
    set.seed(seed_lik)
    
    ## simulate from likelihood with fixed covariance
    y <- matrix(NA, nrow = n_regions, ncol = 2)
    mu <- matrix(NA, nrow = n_regions, ncol = 2)
    for (i in 1:n_regions) {
        mu[i, ] <- c(mean_pars["beta[1]"] + (convolved_re_1[i] * mean_pars["sigma[1]"]) + (mean_pars["lambda"] * convolved_re_2[i]), 
                     mean_pars["beta[2]"] + (convolved_re_2[i] * mean_pars["sigma[2]"]))
        y[i, ] <- rmvnorm(1, mu[i, ], V.array[i,,])
    }
    
    # simulate V_des
    V_hat_array <- V.array
    for (i in 1:n_regions) {
        V_hat_array[i,,] <- rWishart(1, n_r[i]-1, V.array[i,,])/(n_r[i]-1)
    }
    num_corrs <- dim(V_hat_array)[1]
    my.corrs_hat <- rep(NA, num_corrs)
    for (j in 1:num_corrs) {
        Vtmp <- V_hat_array[j,,]
        D <- diag(sqrt(diag(Vtmp)))
        DInv <- solve(D)
        corr.tmp <- DInv %*% Vtmp %*% DInv
        my.corrs_hat[j] <- corr.tmp[1,2]
    }

    ## reformat data into long form
    results <- tibble(admin1 = 1:n_regions,
                      meanHAZ.bi = y[, 1],
                      meanWAZ.bi = y[, 2],
                      seHAZ.bi = V_hat_array[,1,1],
                      seWAZ.bi = V_hat_array[,2,2],
                      corr.bi = my.corrs_hat)
    results.long <- results %>% select(admin1,
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
    data$weights <- rep(1, n_regions)
    
    params <- list(V_pars = V_pars,
                   V.array = V.array,
                   mean_pars = mean_pars,
                   REs = list(v_1, v_2, u_1, u_2, convolved_re_1, convolved_re_2),
                   bivariate_means = mu)
    
    # output
    output <- list(datlist = data,
                   params = params,
                   dattibble = results)
    
    # return the data
    return(output)
}

# function: 
#   fitSTAN()
# purpose: 
#   fit a specified model in STAN with certain STAN fitting parameters and data
# input:
#   stan_file: model to be fit in STAN, one of:
#        c("nonshared-bym2.stan", 
#          "shared-bym2-coregionalization.stan")
#   data: list of input data for the specific STAN model
#   niter: number of TOTAL iterations for HMC
#   nchains: number of chains to run
#   nthin: number of iterations to thin
#   prop_warmup: percentage of the total iterations to use for warming up the HMC
#   max_treedepth: max treedepth for the HMC
#   inits: function defining initial values in STAN, if desired
# output:
#   list with the following components:
#       stan_file: the file for which STAN model was run
#       mod_stan: output from STAN model
#       elapsed_time: total time it took to run the STAN model
fitSTAN <- function(stan_file, data,
                    niter, nchains, nthin, prop_warmup,
                    max_treedepth = NULL, adapt_delta = NULL,
                    inits = NULL,
                    cmd = TRUE) {
    # checks
    if (!(stan_file %in% c("stan-models/nonshared-bym2.stan", 
                           "stan-models/shared-bym2-coregionalization.stan"))) {
        stop("this STAN model is not supported")
    }
    
    start.time <- proc.time() 
    cat(paste("Starting stan modeling! \n"))
    if (cmd == FALSE) {
        if (is.null(inits)) {
            mod_stan <- stan(file = stan_file,
                             data = data,
                             iter = niter, chains = nchains, thin = nthin, 
                             warmup = niter*prop_warmup,
                             control = list(max_treedepth = max_treedepth,
                                            adapt_delta = adapt_delta),
                             seed = 530)
        } else {
            mod_stan <- stan(file = stan_file,
                             data = data,
                             iter = niter, chains = nchains, thin = nthin, 
                             warmup = niter*prop_warmup,
                             init = inits,
                             control = list(max_treedepth = max_treedepth,
                                            adapt_delta = adapt_delta),
                             seed = 530)
        }
    } else {
        if (is.null(inits)) {
            cmd_mod <- cmdstan_model(stan_file = stan_file)
            fit <- cmd_mod$sample(data = data,
                                  iter_warmup = niter*prop_warmup, iter_sampling = niter*(1-prop_warmup),
                                  chains = nchains, thin = nthin,
                                  adapt_delta = adapt_delta, max_treedepth = max_treedepth,
                                  refresh = 0,
                                  seed = 530)
            mod_stan <- rstan::read_stan_csv(fit$output_files())
        } else {
            cmd_mod <- cmdstan_model(stan_file = stan_file)
            fit <- cmd_mod$sample(data = data,
                                  iter_warmup = niter*prop_warmup, iter_sampling = niter*(1-prop_warmup),
                                  chains = nchains, thin = nthin,
                                  adapt_delta = adapt_delta, max_treedepth = max_treedepth,
                                  init = inits,
                                  refresh = 0,
                                  seed = 530)
            mod_stan <- rstan::read_stan_csv(fit$output_files())
        }
    }
    stop.time <- proc.time()
    
    cat(paste("Creating output! \n"))
    output <- list(stan_file = stan_file,
                   mod_stan = mod_stan,
                   elapsed_time = stop.time[3] - start.time[3])
    cat(paste("Done! \n"))
    return(output)
}

fitAllINLAmodels <- function(simulated_data, testing = FALSE) {
    
    if (testing) {
        simulated_data <- simulated_data
    }
    
    data <- simulated_data$datlist
    Nt <- length(unique(data$admin1))
    N <- length(data$admin1)
    n_regions <- Nt
    
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
    data$admin1.haz.2 <- data$admin1.haz
    
    ## add indeces for iid components in besag + iid formulation
    data$admin1.haz.iid <- data$admin1.haz
    data$admin1.waz.iid <- data$admin1.waz
    
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
        init.prec.haz <-  log(simulated_data$dattibble$seHAZ.bi[j]^(-2))
        init.prec.waz <-  log(simulated_data$dattibble$seWAZ.bi[j]^(-2))
        
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
        corr <-  simulated_data$dattibble$corr.bi[j]
        init.prec.haz <-  log(simulated_data$dattibble$seHAZ.bi[j]^(-2))
        init.prec.waz <-  log(simulated_data$dattibble$seWAZ.bi[j]^(-2))
        
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
    formula.bivariate.shared.bym.alt <-  formula(paste("value ~ -1 + outcome + 
                                                   f(admin1.haz, model = 'besag',
                                                     graph = admin1.mat, 
                                                     scale.model = T, 
                                                     constr = T,
                                                     hyper = besag_prior) +
                                                   f(admin1.haz.iid, model = 'iid', hyper = iid_prior) +
                                                   f(admin1.waz, model = 'besag',
                                                     graph = admin1.mat, 
                                                     scale.model = T, 
                                                     constr = T,
                                                     hyper = besag_prior) +
                                                   f(admin1.waz.iid, model = 'iid', hyper = iid_prior) +
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
    
    ## HAZ shared if besag + iid
    lc.all.haz.shared.alt <- inla.make.lincombs(outcomeHAZ = lc.vec.fe, 
                                                admin1.haz = lc.vec.admin1.re,
                                                admin1.haz.iid = lc.vec.admin1.re,
                                                admin1.haz.2 = lc.vec.admin1.re)
    names(lc.all.haz.shared.alt) <- gsub("lc", "haz.reg.", names(lc.all.haz.shared.alt))
    
    # Fit models ####
    fitINLA <- function(formula, data, lincombs) {
        inla(formula, data = data, 
             family = "gaussian",
             lincomb = lincombs,
             control.family = list(hyper = list(prec = list(initial = 10, fixed=T))), 
             control.predictor = list(compute=T),
             control.compute = list(config=T, waic = TRUE, dic = TRUE, cpo = TRUE),
             control.inla = list(lincomb.derived.correlation.matrix=T),
             control.fixed = list(prec = list(default = 0.001), correlation.matrix=T),
             quantiles=c(0.025, 0.1, 0.5, 0.9, 0.975))
    }
    
    cat(paste("Univariate IID \n"))
    mod.univariate.iid <- fitINLA(formula = formula.univariate.iid, 
                                  data = data, 
                                  lincombs = c(lc.all.haz.nonshared, lc.all.waz) )
    
    cat(paste("Univariate BYM \n"))
    mod.univariate.bym <- fitINLA(formula = formula.univariate.bym, 
                                  data = data, 
                                  lincombs = c(lc.all.haz.nonshared, lc.all.waz))
    
    cat(paste("Bivariate nonshared IID \n"))
    mod.bivariate.nonshared.iid <- fitINLA(formula = formula.bivariate.nonshared.iid, 
                                           data = data, 
                                           lincombs = c(lc.all.haz.nonshared, lc.all.waz))
    
    cat(paste("Bivariate nonshared BYM \n"))
    mod.bivariate.nonshared.bym <- fitINLA(formula = formula.bivariate.nonshared.bym, 
                                           data = data, 
                                           lincombs = c(lc.all.haz.nonshared, lc.all.waz))
    
    cat(paste("Bivariate shared IID \n"))
    mod.bivariate.shared.iid <- fitINLA(formula = formula.bivariate.shared.iid, 
                                        data = data, 
                                        lincombs = c(lc.all.haz.shared, lc.all.waz))
    
    cat(paste("Bivariate shared BYM \n"))
    mod.bivariate.shared.bym <- fitINLA(formula = formula.bivariate.shared.bym, 
                                        data = data, 
                                        lincombs = c(lc.all.haz.shared, lc.all.waz))

    inla.results <- list(mod.univariate.iid,
                         mod.univariate.bym,
                         mod.bivariate.nonshared.iid,
                         mod.bivariate.nonshared.bym,
                         mod.bivariate.shared.iid,
                         mod.bivariate.shared.bym)
    names(inla.results) <- c("Univariate IID",
                             "Univariate BYM",
                             "Bivariate nonshared IID", 
                             "Bivariate nonshared BYM",
                             "Bivariate shared IID", 
                             "Bivariate shared BYM")
    
    return(inla.results)
}
