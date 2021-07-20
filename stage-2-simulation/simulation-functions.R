# function: 
#   simulateData()
# purpose: 
#   simulate data to fit a STAN model on 
#   for a given set of dimensions, parameters, 
#   and data generating mechanism
# input:
#   R: number of regions
#   I: number of replicationsk
#   C: number of causes
#   beta: fixed intercepts on cause
#   rho_lower, rho_upper: lower and upper bounds for correlation
#                         of fixed covariance (generated form runif())
#   sigmasq_lower, sigmasq_upper: lower and upper bounds for variance components
#                                 of fixed covariance (generated form runif())
#   sigma_gamma: SD of IID random effects (on region or on region x cause)
#   lambda: scaling parameter for shared component specification
#   sigma_delta: SD of shared IID RE (on region) for the shared RE specification
#   rho_gamma: correlation for bivariate RE specification
#   dgm: data generating mechanism, one of:
#        c("cause FE only",
#          "cause FE, region IID RE sum to zero",
#          "cause FE, region IID RE",
#          "cause FE, region x cause IID RE sum to zero",
#          "cause FE, region x cause IID RE",
#          "one cause, region IID RE sum to zero",
#          "one cause, region IID RE",
#          "cause FE, separate region IID RE by cause",
#          "cause FE, separate region IID RE by cause sum to zero",
#          "2 cause FE, shared region IID RE vX", # different versions, v1, v2, etc i.e. X = 1, 2, ...
#          "2 cause FE, bivariate region IID RE",
#          "2 cause FE, bivariate region IID RE noncentered")
#   testing: logical, TRUE if you are just testing the function for troubleshooting
# output:
#   datlist: list of data to be used as input to various STAN models;
#            includes some or all of the following:
#       y: multivariate outcome as a matrix
#       N: total number of bivariate observations
#       C: number of causes
#       R: number of regions
#       I: number of replications
#       Sigma: N x C x C array with all fixed covariance matrices
#       sigma: length N vector with fixed variances in the case with only 1 cause
#       regions: vector of length N specifying which region 
#               each multivariate observation belongs to
#   params: values of parameters used in the simulation
simulateData <- function(R, I, C, 
                         beta = 1:C, 
                         rho_lower = -0.2, rho_upper = 0.2, 
                         sigmasq_lower = 0.05, sigmasq_upper = 0.5,
                         sigma_gamma = 2,
                         lambda = NULL,
                         sigma_delta = NULL,
                         rho_gamma = NULL,
                         dgm,
                         seed = 8008135,
                         testing = FALSE) {
    # testing
    if (testing == TRUE) {
        R = 10
        I = 5
        C = 2
        beta = 1:C 
        rho_lower = -0.2
        rho_upper = 0.2
        sigmasq_lower = 0.05
        sigmasq_upper = 0.5
        sigma_gamma = c(1.5, 2)
        lambda = 0.5
        sigma_delta = 1
        rho_gamma = 0.5
        seed = 8008135
        dgm = "2 cause FE, shared region IID RE"
    }
    
    # set seed
    set.seed(seed)
    
    # checks
    if (length(beta) != C & !(dgm %in% c("2 cause FE, shared region IID RE v3"))) {
        stop("number of betas (cause FEs) must be equal to number of causes for this dgm")
    }
    if (length(beta) != 1 & (dgm %in% c("2 cause FE, shared region IID RE v3"))) {
        stop("only one beta must be specified for this dgm")
    }
    if ((!C %in% c(1, 2))) stop("Only 1 or 2 causes currently supported")
    if (!(dgm %in% c("cause FE only", 
                     "cause FE, region IID RE sum to zero",
                     "cause FE, region IID RE",
                     "cause FE, region x cause IID RE sum to zero",
                     "cause FE, region x cause IID RE",
                     "one cause, region IID RE sum to zero",
                     "one cause, region IID RE",
                     "cause FE, separate region IID RE by cause",
                     "cause FE, separate region IID RE by cause sum to zero",
                     "2 cause FE, shared region IID RE v1",
                     "2 cause FE, shared region IID RE v2",
                     "2 cause FE, shared region IID RE v3",
                     "2 cause FE, shared region IID RE v4",
                     "2 cause FE, shared region IID RE v5",
                     "2 cause FE, bivariate region IID RE",
                     "2 cause FE, bivariate region IID RE noncentered"))) {
        stop("this data generating mechanism not supported")
    }
    if (length(sigma_gamma) > 1 & !(dgm %in% c("cause FE, separate region IID RE by cause",
                                               "cause FE, separate region IID RE by cause sum to zero",
                                               "2 cause FE, shared region IID RE v1",
                                               "2 cause FE, shared region IID RE v2",
                                               "2 cause FE, shared region IID RE v3",
                                               "2 cause FE, shared region IID RE v4",
                                               "2 cause FE, shared region IID RE v5",
                                               "2 cause FE, bivariate region IID RE",
                                               "2 cause FE, bivariate region IID RE noncentered"))) {
        stop(paste0("multiple sigma_gamma parameters specified which is incompatible with this data generating mechanism: ", dgm))
    }
    if (C != 1 & dgm %in% c("one cause, region IID RE sum to zero",
                            "one cause, region IID RE")) {
        stop("This dgm can only support one cause")
    }
    if (C < 2 & dgm %in% c("cause FE, separate region IID RE by cause",
                           "cause FE, separate region IID RE by cause sum to zero")) {
        stop("This dgm can only support 2 or more causes")
    } 
    if (C != 2 & dgm %in% c("2 cause FE, shared region IID RE v1",
                            "2 cause FE, shared region IID RE v2",
                            "2 cause FE, shared region IID RE v3",
                            "2 cause FE, shared region IID RE v4",
                            "2 cause FE, shared region IID RE v5",
                            "2 cause FE, bivariate region IID RE",
                            "2 cause FE, bivariate region IID RE noncentered")) {
        stop(paste0("Only 2 causes currently supported for dgm: ", dgm))
    }
    if (length(sigma_gamma) != C & dgm %in% c("cause FE, separate region IID RE by cause",
                                              "cause FE, separate region IID RE by cause sum to zero",
                                              "2 cause FE, shared region IID RE v1",
                                              "2 cause FE, shared region IID RE v2",
                                              "2 cause FE, shared region IID RE v3",
                                              "2 cause FE, shared region IID RE v4",
                                              "2 cause FE, shared region IID RE v5",
                                              "2 cause FE, bivariate region IID RE",
                                              "2 cause FE, bivariate region IID RE noncentered")) {
        stop("This dgm needs separate sigma_gamma parameters for each cause")
    }
    if (is.null(lambda) & dgm %in% c("2 cause FE, shared region IID RE v1",
                                     "2 cause FE, shared region IID RE v2",
                                     "2 cause FE, shared region IID RE v3",
                                     "2 cause FE, shared region IID RE v4",
                                     "2 cause FE, shared region IID RE v5")) {
        stop(paste0("need to specify lambda for dgm: ", dgm))
    }
    if (is.null(sigma_delta) & dgm %in% c("2 cause FE, shared region IID RE v1",
                                          "2 cause FE, shared region IID RE v2",
                                          "2 cause FE, shared region IID RE v3",
                                          "2 cause FE, shared region IID RE v4",
                                          "2 cause FE, shared region IID RE v5")) {
        stop(paste0("need to specify sigma_delta for dgm: ", dgm))
    }
    if (is.null(rho_gamma) & dgm %in% c("2 cause FE, bivariate region IID RE",
                                        "2 cause FE, bivariate region IID RE noncentered")) {
        stop(paste0("need to specify rho_gamma for dgm: ", dgm))
    }
    
    # parameters
    N <- R * I
    rhos <- runif(N, rho_lower, rho_upper)
    sigmas <- matrix(sqrt(runif(N * C, sigmasq_lower, sigmasq_upper)), 
                     nrow = N, ncol = C)
    
    # simulate data
    gamma_r <- rnorm(R, 0, sigma_gamma)
    if (dgm %in% c("cause FE, region IID RE sum to zero",
                   "one cause, region IID RE sum to zero")) {
        gamma_r <- gamma_r - mean(gamma_r)
    }
    gamma_rc <- rnorm(R * C, 0, sigma_gamma)
    if (dgm == "cause FE, region x cause IID RE sum to zero") {
        gamma_rc <- gamma_rc - mean(gamma_r)
    }
    regions <- rep(1:R, I)
    gamma_r_mat <- matrix(rep(rep(gamma_r, each = C), I), nrow = N, ncol = C, byrow = TRUE)
    gamma_rc_mat <- matrix(rep(gamma_rc, I), nrow = N, ncol = C, byrow = TRUE)
    if (dgm %in% c("cause FE, separate region IID RE by cause",
                   "cause FE, separate region IID RE by cause sum to zero")) {
        gamma_rc <- matrix(NA, nrow = R, ncol = C)
        for (c in 1:C) {
            gamma_rc[, c] <- rnorm(R, 0, sigma_gamma[c])
            if (dgm == "cause FE, separate region IID RE by cause sum to zero") {
                gamma_rc[, c] <- gamma_rc[, c] - mean(gamma_rc[, c])
            }
        }
        gamma_rc_mat <- gamma_rc[rep(1:R, I),]
    }
    if (dgm %in% c("2 cause FE, shared region IID RE v1",
                   "2 cause FE, shared region IID RE v2")) {
        gamma_rc <- matrix(NA, nrow = R, ncol = C)
        for (c in 1:C) {
            gamma_rc[, c] <- rnorm(R, 0, sigma_gamma[c])
            gamma_rc[, c] <- gamma_rc[, c] - mean(gamma_rc[, c])
        }
        gamma_rc_mat <- gamma_rc[rep(1:R, I),]
        delta <- rnorm(R, 0, sigma_delta)
        delta <- delta - mean(delta)
        delta_mat <- delta[rep(1:R, I)]
    }
    if (dgm %in% c("2 cause FE, shared region IID RE v3")) {
        gamma_rc <- matrix(NA, nrow = R, ncol = C)
        for (c in 1:C) {
            gamma_rc[, c] <- rnorm(R, 0, sigma_gamma[c])
            # gamma_rc[, c] <- gamma_rc[, c] - mean(gamma_rc[, c])
        }
        gamma_rc_mat <- gamma_rc[rep(1:R, I),]
        delta <- rnorm(R, beta, sigma_delta)
        # delta <- delta - mean(delta)
        delta_mat <- delta[rep(1:R, I)]
    }
    if (dgm %in% c("2 cause FE, shared region IID RE v4",
                   "2 cause FE, shared region IID RE v5")) {
        gamma_rc <- matrix(NA, nrow = R, ncol = C)
        for (c in 1:C) {
            gamma_rc[, c] <- rnorm(R, 0, sigma_gamma[c])
        }
        gamma_rc_mat <- gamma_rc[rep(1:R, I),]
        delta <- rnorm(R, 0, sigma_delta)
        delta_mat <- delta[rep(1:R, I)]
    }
    if (dgm %in% c("2 cause FE, bivariate region IID RE",
                   "2 cause FE, bivariate region IID RE noncentered")) {
        Omega <- matrix(c(1, rho_gamma, rho_gamma, 1), ncol = 2, nrow = 2)
        V <- diag(sigma_gamma, nrow = 2, ncol = 2) %*% Omega %*% diag(sigma_gamma, nrow = 2, ncol = 2)
        gamma_rc <- rmvnorm(R, rep(0, 2), V)
        gamma_rc_mat <- gamma_rc[rep(1:R, I),]
    }
    if (C == 1) {
        sigmas <- as.vector(sigmas)
        gamma_r_vec <- as.vector(gamma_r_mat)
        y <- rnorm(N, beta + gamma_r_vec, sigmas)
        
        # create data
        if (dgm %in% c("one cause, region IID RE",
                       "one cause, region IID RE sum to zero")) {
            datlist <- list(N = N,
                            R = R,
                            regions = regions,
                            sigma = sigmas,
                            y = y)
            
            # parameters to output
            params <- list(beta = beta,
                           sigma_gamma = sigma_gamma,
                           gamma_r = gamma_r)
        } else {
            stop("this data generating mechanism not yet supported for 1 cause")
        }
        
        # output
        output <- list(datlist = datlist,
                       params = params)
        
        # return the data
        return(output)
    } else {
        Sigma.array <- array(NA, dim = c(N, C, C))
        y <- matrix(NA, nrow = N, ncol = C)
        for (i in 1:N) {
            Rho <- matrix(c(1, rhos[i], rhos[i], 1), ncol = C)
            Sigma.array[i,,] <- diag(sigmas[i, ], nrow = C, ncol = C) %*% 
                Rho %*% 
                diag(sigmas[i, ], nrow = C, ncol = C)
            if (dgm == "cause FE only") {
                y[i, ] <- rmvnorm(1, beta, Sigma.array[i,,])
            } else if (dgm %in% c("cause FE, region IID RE",
                                  "cause FE, region IID RE sum to zero")) {
                y[i, ] <- rmvnorm(1, beta + gamma_r_mat[i, ], Sigma.array[i,,])
            } else if (dgm %in% c("cause FE, region x cause IID RE",
                                  "cause FE, region x cause IID RE sum to zero",
                                  "cause FE, separate region IID RE by cause",
                                  "cause FE, separate region IID RE by cause sum to zero")) {
                y[i, ] <- rmvnorm(1, beta + gamma_rc_mat[i, ], Sigma.array[i,,])
            } else if (dgm %in% c("2 cause FE, shared region IID RE v1",
                                  "2 cause FE, shared region IID RE v2")) {
                mu <- c(beta[1] + gamma_rc_mat[i, 1] + (lambda * delta_mat[i]),
                        beta[2] + gamma_rc_mat[i, 2] + (1/lambda * delta_mat[i]))
                y[i, ] <- rmvnorm(1, mu, Sigma.array[i,,])
            } else if (dgm %in% c("2 cause FE, shared region IID RE v3")) {
                mu <- c(gamma_rc_mat[i, 1] + (delta_mat[i] * lambda),
                        gamma_rc_mat[i, 2] + (delta_mat[i]) / lambda)
                y[i, ] <- rmvnorm(1, mu, Sigma.array[i,,])
            } else if (dgm %in% c("2 cause FE, shared region IID RE v4",
                                  "2 cause FE, shared region IID RE v5")) {
                mu <- c(beta[1] + gamma_rc_mat[i, 1] + (lambda * delta_mat[i]),
                        beta[2] + gamma_rc_mat[i, 2] + (1/lambda * delta_mat[i]))
                y[i, ] <- rmvnorm(1, mu, Sigma.array[i,,])
            } else if (dgm %in% c("2 cause FE, bivariate region IID RE",
                                  "2 cause FE, bivariate region IID RE noncentered")) {
                mu <- c(beta[1] + gamma_rc_mat[i, 1],
                        beta[2] + gamma_rc_mat[i, 2])
                y[i, ] <- rmvnorm(1, mu, Sigma.array[i,,])
            } else {
                stop("this data generating mechanism not yet supported for multiple causes")
            }
        }
        
        # create data
        if (dgm == "cause FE only") {
            # data to be used in STAN
            datlist <- list(N = N,
                            C = C,
                            Sigma = Sigma.array,
                            y = y)
            
            # parameters to output
            params <- list(beta = beta)
            
        } else if (dgm %in% c("cause FE, region IID RE",
                              "cause FE, region IID RE sum to zero")) {
            datlist <- list(N = N,
                            R = R,
                            C = C,
                            regions = regions,
                            Sigma = Sigma.array,
                            y = y)
            
            # parameters to output
            params <- list(beta = beta,
                           sigma_gamma = sigma_gamma,
                           gamma_r = gamma_r)
        } else if (dgm %in% c("cause FE, region x cause IID RE",
                              "cause FE, region x cause IID RE sum to zero",
                              "cause FE, separate region IID RE by cause",
                              "cause FE, separate region IID RE by cause sum to zero")) {
            datlist <- list(N = N,
                            R = R,
                            C = C,
                            regions = regions,
                            Sigma = Sigma.array,
                            y = y)
            
            # parameters to output
            params <- list(beta = beta,
                           sigma_gamma = sigma_gamma,
                           gamma_rc = gamma_rc,
                           gamma_rc_mat = gamma_rc_mat)
        } else if (dgm %in% c("2 cause FE, shared region IID RE v1",
                              "2 cause FE, shared region IID RE v2",
                              "2 cause FE, shared region IID RE v3",
                              "2 cause FE, shared region IID RE v4",
                              "2 cause FE, shared region IID RE v5")) {
            datlist <- list(N = N,
                            R = R,
                            regions = regions,
                            Sigma = Sigma.array,
                            y = y)
            
            # parameters to output
            params <- list(beta = beta,
                           sigma_gamma = sigma_gamma,
                           gamma_rc = gamma_rc,
                           gamma_rc_mat = gamma_rc_mat,
                           lambda = lambda,
                           sigma_delta = sigma_delta,
                           delta_rc = delta,
                           delta_mat = delta_mat)
        } else if (dgm %in% c("2 cause FE, bivariate region IID RE",
                              "2 cause FE, bivariate region IID RE noncentered")) {
            datlist <- list(N = N,
                            R = R,
                            regions = regions,
                            Sigma = Sigma.array,
                            y = y)
            
            # parameters to output
            params <- list(beta = beta,
                           sigma_gamma = sigma_gamma,
                           gamma_rc = gamma_rc,
                           gamma_rc_mat = gamma_rc_mat,
                           rho_gamma = rho_gamma)
        } else {
            stop("this data generating mechanism not yet supported for multiple causes")
        }
        
        # output
        output <- list(datlist = datlist,
                       params = params)
        
        # return the data
        return(output)
    }
}

# function: 
#   fitSTAN()
# purpose: 
#   fit a specified model in STAN with certain STAN fitting parameters and data
# input:
#   stan_model: model to be fit in STAN (matches data generating mechanism), one of:
#        c("cause FE only",
#          "cause FE, region IID RE sum to zero",
#          "cause FE, region IID RE",
#          "cause FE, region x cause IID RE sum to zero",
#          "cause FE, region x cause IID RE",
#          "one cause, region IID RE sum to zero",
#          "one cause, region IID RE",
#          "cause FE, separate region IID RE by cause",
#          "cause FE, separate region IID RE by cause sum to zero",
#          "2 cause FE, shared region IID RE vX", # different versions, v1, v2, etc i.e. X = 1, 2, ...
#          "2 cause FE, bivariate region IID RE",
#          "2 cause FE, bivariate region IID RE noncentered")
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
fitSTAN <- function(stan_model, data,
                    niter, nchains, nthin, prop_warmup,
                    max_treedepth = NULL, adapt_delta = NULL,
                    inits = NULL,
                    cmd = TRUE) {
    # checks
    if (!(stan_model %in% c("cause FE only", 
                            "cause FE, region IID RE sum to zero",
                            "cause FE, region IID RE",
                            "cause FE, region x cause IID RE sum to zero",
                            "cause FE, region x cause IID RE",
                            "one cause, region IID RE sum to zero",
                            "one cause, region IID RE",
                            "cause FE, separate region IID RE by cause",
                            "cause FE, separate region IID RE by cause sum to zero",
                            "2 cause FE, shared region IID RE v1",
                            "2 cause FE, shared region IID RE v2",
                            "2 cause FE, shared region IID RE v3",
                            "2 cause FE, shared region IID RE v4",
                            "2 cause FE, shared region IID RE v5",
                            "2 cause FE, bivariate region IID RE",
                            "2 cause FE, bivariate region IID RE noncentered"))) {
        stop("this STAN model is not supported")
    }
    
    start.time <- proc.time() 
    if (stan_model == "cause FE only") {
        stan_file <- "stan-models/fixed-cov-fe-cause.stan"
    } else if(stan_model == "cause FE, region IID RE") {
        stan_file <- "stan-models/fixed-cov-fe-cause-re-region.stan"
    } else if(stan_model == "cause FE, region IID RE sum to zero") {
        stan_file <- "stan-models/fixed-cov-fe-cause-re-regionS2Z.stan"
    } else if(stan_model == "cause FE, region x cause IID RE") {
        stan_file <- "stan-models/fixed-cov-fe-cause-re-regionXcause.stan"
    } else if(stan_model == "cause FE, region x cause IID RE sum to zero") {
        stan_file <- "stan-models/fixed-cov-fe-cause-re-regionXcauseS2Z.stan"
    } else if(stan_model == "one cause, region IID RE sum to zero") {
        stan_file <- "stan-models/fixed-var1d-fixedintercept-re-regionS2Z.stan"
    } else if(stan_model == "one cause, region IID RE") {
        stan_file <- "stan-models/fixed-var1d-fixedintercept-re-region.stan"
    } else if(stan_model == "cause FE, separate region IID RE by cause") {
        stan_file <- "stan-models/fixed-cov-fe-cause-separate-re-region.stan"
    } else if(stan_model == "cause FE, separate region IID RE by cause sum to zero") {
        stan_file <- "stan-models/fixed-cov-fe-cause-separate-re-regionS2Z.stan"
    } else if(stan_model == "2 cause FE, shared region IID RE v1") {
        stan_file <- "stan-models/fixed-cov-2cause-fe-shared-iid-re-region-v1.stan"
    } else if(stan_model == "2 cause FE, shared region IID RE v2") {
        stan_file <- "stan-models/fixed-cov-2cause-fe-shared-iid-re-region-v2.stan"
    } else if(stan_model == "2 cause FE, shared region IID RE v3") {
        stan_file <- "stan-models/fixed-cov-2cause-fe-shared-iid-re-region-v3.stan"
    } else if(stan_model == "2 cause FE, shared region IID RE v4") {
        stan_file <- "stan-models/fixed-cov-2cause-fe-shared-iid-re-region-v4.stan"
    } else if(stan_model == "2 cause FE, shared region IID RE v5") {
        stan_file <- "stan-models/fixed-cov-2cause-fe-shared-iid-re-region-v5.stan"
    } else if(stan_model == "2 cause FE, bivariate region IID RE") {
        stan_file <- "stan-models/fixed-cov-2cause-fe-bivariate-iid-re-region.stan"
    } else if(stan_model == "2 cause FE, bivariate region IID RE noncentered") {
        stan_file <- "stan-models/fixed-cov-2cause-fe-bivariate-iid-re-region-noncentered.stan"
    } else {
        stop("this STAN model is not yet supported")
    }
    
    if (cmd == FALSE) {
        if (is.null(inits)) {
            mod_stan <- stan(file = stan_file,
                             data = data,
                             iter = niter, chains = nchains, thin = nthin, 
                             warmup = niter*prop_warmup,
                             control = list(max_treedepth = max_treedepth,
                                            adapt_delta = adapt_delta))
        } else {
            mod_stan <- stan(file = stan_file,
                             data = data,
                             iter = niter, chains = nchains, thin = nthin, 
                             warmup = niter*prop_warmup,
                             init = inits,
                             control = list(max_treedepth = max_treedepth,
                                            adapt_delta = adapt_delta))
        }
    } else {
        if (is.null(inits)) {
            cmd_mod <- cmdstan_model(stan_file = stan_file)
            fit <- cmd_mod$sample(data = data,
                                  iter_warmup = niter*prop_warmup, iter_sampling = niter*(1-prop_warmup),
                                  chains = nchains, thin = nthin,
                                  adapt_delta = adapt_delta, max_treedepth = max_treedepth)
            mod_stan <- rstan::read_stan_csv(fit$output_files())
        } else {
            cmd_mod <- cmdstan_model(stan_file = stan_file)
            fit <- cmd_mod$sample(data = data,
                                  iter_warmup = niter*prop_warmup, iter_sampling = niter*(1-prop_warmup),
                                  chains = nchains, thin = nthin,
                                  adapt_delta = adapt_delta, max_treedepth = max_treedepth,
                                  init = inits)
            mod_stan <- rstan::read_stan_csv(fit$output_files())
        }
    }
    stop.time <- proc.time()
    output <- list(stan_file = stan_file,
                   mod_stan = mod_stan,
                   elapsed_time = stop.time[3] - start.time[3])
    return(output)
}