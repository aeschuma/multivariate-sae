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
simulateData <- function(dgm_specs, Amat, scaling_factor, seed_re, seed_lik, testing = FALSE) {
    
    # testing
    if (testing == TRUE) {
        dgm_specs <- my_dgm
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
    
    # load models to pull from 
    numeric_pars <- suppressWarnings(as.numeric(my_dgm))
    if (any(is.na(numeric_pars))) {
        models_to_load <- unique(my_dgm[is.na(numeric_pars)])
        mod_lists <- vector(mode = "list", length = length(models_to_load))
        names(mod_lists) <- models_to_load
        for (i in 1:length(mod_lists)) {
            # check par model is legit
            if (!(models_to_load[i] %in% c("ken2014-hazwaz-stage-2-bivariate-shared-coreg-b11",
                                           "ken2014-hazwaz-stage-2-bivariate-nonshared-b11",
                                           "ken2014-hazwaz-stage-2-bivariate-shared-coreg-b33",
                                           "ken2014-hazwaz-stage-2-bivariate-nonshared-b33"))) {
                stop(paste0("Model ",models_to_load[i] ," not yet supported"))
            }
            
            # load model
            mod_lists[[i]] <- read_rds(paste0("~/Dropbox/dissertation_2/survey-csmf/results/ken2014-hazwaz/", models_to_load[i], ".rds"))
        }
    }
    
    # set parameters
    all_par_names <- names(my_dgm)
    
    V_pars <- my_dgm[grep("V_", all_par_names, value = TRUE)]
    V_pars_numeric <- suppressWarnings(as.numeric(V_pars))
    if (length(unique(V_pars)) == 1 & sum(is.na(V_pars_numeric)) == length(V_pars_numeric)) { # if all parameters are from the same model
        this.V.mod <- unique(V_pars)
        V.array <- mod_lists[[this.V.mod]]$data$Sigma
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
            V.array.tmp <- mod_lists[[corr_par]]$data$Sigma
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
            V.array.tmp <- mod_lists[[vdiag_par]]$data$Sigma
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
                mean_pars[i] <- summary(mod_lists[[par_mod.tmp]][["mod"]][[1]], par = mean_pars_names[i])$summary[,"50%"]
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
        mu[i, ] <- c(mean_pars["beta[1]"] + (convolved_re_1[i] * mean_pars["sigma[1]"]) + (mean_pars["lambda"] * u_2[i]), 
                     mean_pars["beta[2]"] + (convolved_re_2[i] * mean_pars["sigma[2]"]))
        y[i, ] <- rmvnorm(1, mu[i, ], V.array[i,,])
    }
        
    datlist <- list(R = n_regions,
                    regions = 1:n_regions,
                    Sigma = V.array,
                    y = y,
                    N_edges = node.info$n_edges,
                    node1 = node.info$node1,
                    node2 = node.info$node2,
                    scaling_factor = scaling_factor,
                    rho_beta_1 = NA,
                    rho_beta_2 = NA,
                    sigma_normal_sd = NA)
    
    RE_list <- list(u_1 = u_1,
                    u_2 = u_2,
                    v_1 = v_1,
                    v_2 = v_2,
                    convolved_re_1 = convolved_re_1,
                    convolved_re_2 = convolved_re_2)
    
    params <- list(V_pars = V_pars,
                   mean_pars = mean_pars,
                   REs = RE_list,
                   bivariate_means = mu)
    # output
    output <- list(datlist = datlist,
                   params = params)
    
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