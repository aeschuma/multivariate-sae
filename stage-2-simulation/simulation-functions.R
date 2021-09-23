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
simulateData <- function(dgm, Amat, scaling_factor, seed, testing = FALSE) {
    
    # testing
    if (testing == TRUE) {
        Amat <- admin1.mat
        scaling_factor <- scaling_factor
        dgm <- 1
        seed <- 80085
    }
    
    n_regions <- nrow(Amat)
    
    # load dgm data
    dgm_file <- read_csv("dgm-info.csv")
    my_dgm <- dgm_file %>% filter(dgm_number == dgm) %>% select(-1) %>% slice(1) %>% unlist()
    data_based <- as.numeric(my_dgm["data_based"])
    my_dgm <- my_dgm[-1:-2]
    
    # check parameter model and set
    if (data_based == 1 & length(unique(my_dgm)) != 1) stop("Not all parameters are from the same data-based model")
    par_mod <- NA
    if (data_based == 1) {
        par_mod <- unique(my_dgm)
        mod_list <- read_rds(paste0("~/Dropbox/dissertation_2/survey-csmf/results/ken2014-hazwaz/", par_mod, ".rds"))
        
    } else {
        mod_list <- NA
    }
    
    # check par model is legit
    if (data_based == 1 & !(par_mod %in% c("ken2014-hazwaz-stage-2-bivariate-shared-coreg",
                                           "ken2014-hazwaz-stage-2-bivariate-nonshared"))) {
        stop("Model for parameter values not yet supported")
    }
    
    # set seed
    set.seed(seed)
    
    # set parameters
    all_par_names <- names(my_dgm)
    
    V_pars <- my_dgm[grep("V_", all_par_names, value = TRUE)]
    if (data_based == 1) {
        V.array <- mod_list$data$Sigma
    } else (
        stop("Creating fixed covariance matrix for non-data-based parameters not yet implemented")
    )
    
    mean_pars_names <- all_par_names[!(all_par_names %in% names(V_pars))]
    mean_pars <- suppressWarnings(as.numeric(my_dgm[mean_pars_names]))
    if (data_based == 1) {
        for (i in 1:length(mean_pars)) {
            if (mean_pars_names[i] == "lambda" & par_mod == "ken2014-hazwaz-stage-2-bivariate-nonshared") {
                mean_pars[i] <- 0
            } else {
                mean_pars[i] <- summary(mod_list[["mod"]][[1]], par = mean_pars_names[i])$summary[,"50%"]
            }
        }
    }
    names(mean_pars) <- mean_pars_names
    
    # simulate data!
    
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
    
    ## simulate from likelihood with fixed covariance
    y <- matrix(NA, nrow = n_regions, ncol = 2)
    for (i in 1:n_regions) {
        mu <- c(mean_pars["beta[1]"] + (convolved_re_1[i] * mean_pars["sigma[1]"]) + (mean_pars["lambda"] * u_2[i]), 
                mean_pars["beta[2]"] + (convolved_re_2[i] * mean_pars["sigma[2]"]))
        y[i, ] <- rmvnorm(1, mu, V.array[i,,])
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
                   REs = RE_list)
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
                                  adapt_delta = adapt_delta, max_treedepth = max_treedepth,
                                  refresh = 0)
            mod_stan <- rstan::read_stan_csv(fit$output_files())
        } else {
            cmd_mod <- cmdstan_model(stan_file = stan_file)
            fit <- cmd_mod$sample(data = data,
                                  iter_warmup = niter*prop_warmup, iter_sampling = niter*(1-prop_warmup),
                                  chains = nchains, thin = nthin,
                                  adapt_delta = adapt_delta, max_treedepth = max_treedepth,
                                  init = inits,
                                  refresh = 0)
            mod_stan <- rstan::read_stan_csv(fit$output_files())
        }
    }
    stop.time <- proc.time()
    output <- list(stan_file = stan_file,
                   mod_stan = mod_stan,
                   elapsed_time = stop.time[3] - start.time[3])
    return(output)
}