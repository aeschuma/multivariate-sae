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
simulateData <- function(dgm_specs, sample_info, Amat, scaling_factor, seed_re, seed_lik, testing = FALSE) {
    
    # testing
    if (testing == TRUE) {
        dgm_specs <- my_dgm
        sample_info <- sample_info
        Amat <- admin1.mat
        scaling_factor <- scaling_factor
        seed_re <- 80085
        seed_lik <- 8008135
    }
    
    # sample information
    n_regions <- nrow(Amat)
    n_cause <- 2
    n_outcomes <- n_regions * n_cause
    sample_data <- sample_info %>% uncount(n)
    
    # load dgm data
    my_dgm <- dgm_specs %>% unlist()
    dgm_model <- as.character(my_dgm["mech"])
    my_dgm <- my_dgm[which(names(my_dgm) != "mech")]
    
    # load model results
    mod_lists <- readRDS("~/Dropbox/dissertation_2/survey-csmf/results/ken2014-unit-level/hazwaz/ken2014-unit-level-hazwaz-inla-posteriors.rds")

    # set parameters
    all_par_names <- names(my_dgm)

    all_pars.tmp <- my_dgm[all_par_names]
    all_pars_num.tmp <- suppressWarnings(as.numeric(all_pars.tmp))
    
    all_pars <- rep(NA, length(all_pars.tmp))
    
    for (i in 1:length(all_pars)) {
        if (!is.na(all_pars_num.tmp[i])) {
            all_pars[i] <- all_pars_num.tmp[i]
        } else {
            par_mod.tmp <- all_pars.tmp[i]
            if (all_par_names[i] == "lambda" & grepl("nonshared", par_mod.tmp)) {
                all_pars[i] <- 0
            } else {
                if (names(par_mod.tmp) == "beta[1]") {
                    all_pars[i] <- mod_lists[[par_mod.tmp]]$summaries$summary.fixed["outcomeHAZ", "0.5quant"]
                } else if (names(par_mod.tmp) == "beta[2]") {
                    all_pars[i] <- mod_lists[[par_mod.tmp]]$summaries$summary.fixed["outcomeWAZ", "0.5quant"]
                } else if (names(par_mod.tmp) == "sigma[1]") {
                    all_pars[i] <- mod_lists[[par_mod.tmp]]$summaries$summary.hyperpar["Precision for admin1.haz", "0.5quant"]^(-0.5)
                } else if (names(par_mod.tmp) == "sigma[2]") {
                    all_pars[i] <- mod_lists[[par_mod.tmp]]$summaries$summary.hyperpar["Precision for admin1.waz", "0.5quant"]^(-0.5)
                } else if (names(par_mod.tmp) == "rho[1]") {
                    all_pars[i] <- mod_lists[[par_mod.tmp]]$summaries$summary.hyperpar["Phi for admin1.haz", "0.5quant"]
                } else if (names(par_mod.tmp) == "rho[2]") {
                    all_pars[i] <- mod_lists[[par_mod.tmp]]$summaries$summary.hyperpar["Phi for admin1.waz", "0.5quant"]
                } else if (names(par_mod.tmp) == "lambda") {
                    all_pars[i] <- mod_lists[[par_mod.tmp]]$summaries$summary.hyperpar["Beta for admin1.haz.2", "0.5quant"]
                } else if (names(par_mod.tmp) == "gamma[1]") {
                    all_pars[i] <- mod_lists[[par_mod.tmp]]$summaries$summary.fixed["rural.haz", "0.5quant"]
                } else if (names(par_mod.tmp) == "gamma[2]") {
                    all_pars[i] <- mod_lists[[par_mod.tmp]]$summaries$summary.fixed["rural.waz", "0.5quant"]
                } else if (names(par_mod.tmp) == "sigma_epsilon[1]") {
                    all_pars[i] <- mod_lists[[par_mod.tmp]]$summaries$summary.hyperpar["Precision for cluster.haz", "0.5quant"]^(-0.5)
                } else if (names(par_mod.tmp) == "sigma_epsilon[2]") {
                    all_pars[i] <- mod_lists[[par_mod.tmp]]$summaries$summary.hyperpar["Precision for cluster.waz", "0.5quant"]^(-0.5)
                } else if (names(par_mod.tmp) == "gaussian_sd[1]") {
                    all_pars[i] <- mod_lists[[par_mod.tmp]]$summaries$summary.hyperpar["Precision for the Gaussian observations", "0.5quant"]^(-0.5)
                } else if (names(par_mod.tmp) == "gaussian_sd[2]") {
                    all_pars[i] <- mod_lists[[par_mod.tmp]]$summaries$summary.hyperpar["Precision for the Gaussian observations[2]", "0.5quant"]^(-0.5)
                } 
            }
        }
    }
    names(all_pars) <- all_par_names
    
    # set seed for random effects
    set.seed(seed_re)
    
    ## IID normal(0, 1)
    v_haz <- rnorm(n_regions, 0, 1)
    v_waz <- rnorm(n_regions, 0, 1)
    
    if (grepl("BYM", dgm_model)) {
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
        
        u_haz <- sim.Q(Q)
        u_waz <- sim.Q(Q) 
        
        ## Convolved REs
        convolved_re_haz <- all_pars["sigma[1]"]  * ((sqrt(1 - all_pars["rho[1]"]) * v_haz) + (sqrt(all_pars["rho[1]"] / scaling_factor) * u_haz))
        convolved_re_waz <- all_pars["sigma[2]"]  * ((sqrt(1 - all_pars["rho[2]"]) * v_waz) + (sqrt(all_pars["rho[2]"] / scaling_factor) * u_waz))
    } else {
        u_haz <- u_waz <- convolved_re_haz <- convolved_re_waz <- rep(0, n_regions)
    }
    
    # merge on REs
    re_data <- tibble(admin1 = 1:n_regions,
                      v_haz = v_haz * all_pars["sigma[1]"],
                      v_waz = v_waz * all_pars["sigma[2]"],
                      u_haz = u_haz  / sqrt(scaling_factor),
                      u_waz = u_waz / sqrt(scaling_factor),
                      convolved_re_haz = convolved_re_haz,
                      convolved_re_waz = convolved_re_waz)
    if (grepl("BYM", dgm_model)) {
        re_data %<>% mutate(spatial_re_haz = convolved_re_haz,
                            spatial_re_waz = convolved_re_waz)
    } else if (grepl("IID", dgm_model)) {
        re_data %<>% mutate(spatial_re_haz = v_haz,
                            spatial_re_waz = v_waz)
    } else {
        stop("Need either a BYM or IID model!")
    }
    sample_data %<>% left_join(re_data)
    
    # merge on FEs
    fe_data <- tibble(rural = c(0,1),
                      beta_haz = rep(all_pars["beta[1]"], 2),
                      beta_waz = rep(all_pars["beta[2]"], 2),
                      gamma_haz = c(0, all_pars["gamma[1]"]),
                      gamma_waz = c(0, all_pars["gamma[2]"])) %>%
        mutate(final_fe_haz = beta_haz + gamma_haz,
               final_fe_waz = beta_waz + gamma_waz)
    sample_data %<>% left_join(fe_data)
    
    # merge on cluster effects
    sample_data %<>% group_by(cluster) %>%
        mutate(epsilon_haz = rnorm(1, 0, all_pars["sigma_epsilon[1]"]),
               epsilon_waz = rnorm(1, 0, all_pars["sigma_epsilon[2]"])) %>%
        ungroup()
    
    # calculate latent means
    sample_data %<>% mutate(mu_haz = final_fe_haz + spatial_re_haz + (all_pars["lambda"] * spatial_re_waz) + epsilon_haz,
                            mu_waz = final_fe_waz + spatial_re_waz + epsilon_waz,
                            admin1_rural_HAZ = final_fe_haz + spatial_re_haz + (all_pars["lambda"] * spatial_re_waz),
                            admin1_rural_WAZ = final_fe_waz + spatial_re_waz)
    
    # parameter data to save
    latent_params <- sample_data
    
    # set seed for the likelihood
    set.seed(seed_lik)
    
    ## simulate individual values from likelihood 
    sample_data <- latent_params %>% mutate(y_haz = rnorm(n(), mu_haz, all_pars["gaussian_sd[1]"]),
                                            y_waz = rnorm(n(), mu_waz, all_pars["gaussian_sd[2]"])) %>% 
        select(admin1, rural, cluster, y_haz, y_waz) %>%
        pivot_longer(cols = c("y_haz", "y_waz"),
                     names_prefix = "y_",
                     names_to = "outcome",
                     values_to = "value") %>%
        mutate(outcome = toupper(outcome))
    
    sample_data$value.haz <- ifelse(sample_data$outcome == "HAZ", sample_data$value, NA)
    sample_data$value.waz <- ifelse(sample_data$outcome == "WAZ", sample_data$value, NA)
    
    sample_data$admin1.haz <- ifelse(sample_data$outcome == "HAZ", sample_data$admin1, NA)
    sample_data$admin1.waz <- ifelse(sample_data$outcome == "WAZ", sample_data$admin1, NA)
    sample_data$obs <- 1:nrow(sample_data)
    
    ## add outcome specific cluster indicators for nugget terms
    sample_data$cluster.haz <- ifelse(sample_data$outcome == "HAZ", sample_data$cluster, NA)
    sample_data$cluster.waz <- ifelse(sample_data$outcome == "WAZ", sample_data$cluster, NA)
    
    ## add outcome specific urban/rural indicators for strata FEs
    sample_data$rural.haz <- ifelse(sample_data$outcome == "HAZ", sample_data$rural, NA)
    sample_data$rural.waz <- ifelse(sample_data$outcome == "WAZ", sample_data$rural, NA)
    
    ## add another index for the shared component
    sample_data$admin1.haz.2 <- sample_data$admin1.haz
    
    # create a list of the data
    data <- as.list(sample_data)
    
    # create matrix for outcome in order to have two likelihoods
    data$Y <- cbind(data$value.haz, data$value.waz)
    
    # output
    output <- list(data = data,
                   params = all_pars,
                   latent_params = latent_params,
                   Amat = Amat)
    
    # return the data
    return(output)
}

fitAllINLAmodels <- function(simulated_data, nsamps, seed = 1, testing = FALSE) {
    
    if (testing) {
        simulated_data <- simulated_data
        seed <- run_number + 1
        nsamps <- 100
    }
    
    message("format data and set up model")
    
    data <- simulated_data$data
    admin1.mat <- simulated_data$Amat
    
    # define model components ####
    
    # get number of regions
    n_regions <- length(unique(data$admin1))
    
    # priors ####
    iid_prior <- list(prec = list(prior = "pc.prec",
                                  param = c(1, 0.01)))
    bym2_prior <- list(phi=list(prior="logitbeta", param=c(1, 1), initial=0.5), 
                       prec=list(prior="pc.prec", param=c(1, 0.01), initial=5))
    lambda_prior <- list(beta = list(prior = 'logtnormal', param = c(0, 1)))
    
    # models to run  ####
    model_names <- c("IID nonshared", "BYM nonshared",
                     "IID shared", "BYM shared")
    
    # formulas ####
    formulas <- vector(mode = "list", length = length(model_names))
    names(formulas) <- model_names
    
    formulas[["IID nonshared"]] <- formula("Y ~ -1 + outcome + rural.haz + rural.waz + 
                                        f(admin1.haz, model = 'iid', hyper = iid_prior) +
                                        f(admin1.waz, model = 'iid', hyper = iid_prior) +
                                        f(cluster.haz, model = 'iid', hyper = iid_prior) +
                                        f(cluster.waz, model = 'iid', hyper = iid_prior)")
    
    formulas[["BYM nonshared"]] <- formula("Y ~ -1 + outcome + rural.haz + rural.waz + 
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
                                        f(cluster.haz, model = 'iid', hyper = iid_prior) +
                                        f(cluster.waz, model = 'iid', hyper = iid_prior)")
    formulas[["IID shared"]] <- formula("Y ~ -1 + outcome + rural.haz + rural.waz +
                                        f(admin1.haz, model = 'iid', hyper = iid_prior) +
                                        f(admin1.waz, model = 'iid', hyper = iid_prior) +
                                        f(admin1.haz.2, copy = \"admin1.waz\", 
                                          fixed = FALSE, hyper = lambda_prior) +
                                        f(cluster.haz, model = 'iid', hyper = iid_prior) +
                                        f(cluster.waz, model = 'iid', hyper = iid_prior)")
    formulas[["BYM shared"]] <- formula("Y ~ -1 + outcome + rural.haz + rural.waz + 
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
                                          fixed = FALSE, hyper = lambda_prior) +
                                        f(cluster.haz, model = 'iid', hyper = iid_prior) +
                                        f(cluster.waz, model = 'iid', hyper = iid_prior)")
    
    # results storage
    results <- vector(mode = "list", length = length(model_names))
    names(results) <- model_names
    
    message("fitting models")
    
    for (i in 1:length(model_names)) {
        
        message(paste0("Fitting ", model_names[i]))
        tmp <- inla(formulas[[ model_names[i] ]], data = data,
                    family = rep("gaussian", 2),
                    quantiles=c(0.025, 0.1, 0.5, 0.9, 0.975),
                    # control.family = list(),
                    control.predictor = list(compute=T),
                    # control.inla = list(lincomb.derived.correlation.matrix=T),
                    # control.fixed = list(prec = list(default = 0.001), correlation.matrix=T),
                    control.compute = list(config=T, waic = F, dic = F, cpo = F))
        
        message(paste0("Posterior sampling ", model_names[i]))
        
        samp <- inla.posterior.sample(n = nsamps, result = tmp, seed = seed)
        
        message(paste0("Posterior summaries ", model_names[i]))
        
        # process hyperpars
        hyperpar_names <- names(samp[[1]]$hyperpar)
        hyperpar_mat <- matrix(NA, nrow = length(hyperpar_names), ncol = nsamps)
        for (s in 1:nsamps) {
            hyperpar_mat[,s] <- samp[[s]]$hyperpar
        }
        rownames(hyperpar_mat) <- hyperpar_names
        
        # process fitted values
        haz_fe_idx <- rownames(samp[[1]]$latent) %>% str_detect("outcomeHAZ:1") %>% which()
        waz_fe_idx <- rownames(samp[[1]]$latent) %>% str_detect("outcomeWAZ:1") %>% which()
        haz_fe_mat <- matrix(NA, nrow = length(haz_fe_idx), ncol = nsamps)
        waz_fe_mat <- matrix(NA, nrow = length(waz_fe_idx), ncol = nsamps)
        
        haz_rural_fe_idx <- rownames(samp[[1]]$latent) %>% str_detect("rural\\.haz") %>% which()
        waz_rural_fe_idx <- rownames(samp[[1]]$latent) %>% str_detect("rural\\.waz") %>% which()
        haz_rural_fe_mat <- matrix(NA, nrow = length(haz_rural_fe_idx), ncol = nsamps)
        waz_rural_fe_mat <- matrix(NA, nrow = length(waz_rural_fe_idx), ncol = nsamps)
        
        haz_re_idx <- rownames(samp[[1]]$latent) %>% str_detect("admin1\\.haz:") %>% which()
        waz_re_idx <- rownames(samp[[1]]$latent) %>% str_detect("admin1\\.waz:") %>% which()
        haz_re_mat <- matrix(NA, nrow = length(haz_re_idx), ncol = nsamps)
        waz_re_mat <- matrix(NA, nrow = length(waz_re_idx), ncol = nsamps)
        
        if(!grepl("nonshared", model_names[i])) {
            shared_re_idx <- rownames(samp[[1]]$latent) %>% str_detect("admin1\\.haz\\.2:") %>% which()
            shared_re_mat <- matrix(NA, nrow = length(shared_re_idx), ncol = nsamps)
        } else {
            shared_re_idx <- integer(0)
            shared_re_mat <- matrix(0, nrow = nrow(haz_re_mat), ncol = nsamps)
        }
        
        # fill in sample matrices
        for (s in 1:nsamps) {
            haz_fe_mat[,s] <- samp[[s]]$latent[haz_fe_idx]
            waz_fe_mat[,s] <- samp[[s]]$latent[waz_fe_idx]
            haz_rural_fe_mat[,s] <- samp[[s]]$latent[haz_rural_fe_idx]
            waz_rural_fe_mat[,s] <- samp[[s]]$latent[waz_rural_fe_idx]
            haz_re_mat[,s] <- samp[[s]]$latent[haz_re_idx]
            waz_re_mat[,s] <- samp[[s]]$latent[waz_re_idx]
            if(!grepl("nonshared", model_names[i])) shared_re_mat[,s] <- samp[[s]]$latent[shared_re_idx]
        }
        
        # obtain total effect for space - first half of bym2, or all of them for IID
        haz_re_tot_mat <- haz_re_mat[1:n_regions,] + shared_re_mat[1:n_regions,]
        waz_re_tot_mat <- waz_re_mat[1:n_regions,]
        
        # get urban and rural FEs
        haz_urban_fe_tot_mat <- haz_fe_mat
        haz_rural_fe_tot_mat <- haz_fe_mat + haz_rural_fe_mat
        waz_urban_fe_tot_mat <- waz_fe_mat
        waz_rural_fe_tot_mat <- waz_fe_mat + waz_rural_fe_mat
        
        # matrix of fitted estimates
        fitted_haz_urban_mat <- haz_urban_fe_tot_mat[rep(1,n_regions),] + haz_re_tot_mat
        fitted_waz_urban_mat <- waz_urban_fe_tot_mat[rep(1,n_regions),] + waz_re_tot_mat
        fitted_haz_rural_mat <- haz_rural_fe_tot_mat[rep(1,n_regions),] + haz_re_tot_mat
        fitted_waz_rural_mat <- waz_rural_fe_tot_mat[rep(1,n_regions),] + waz_re_tot_mat
        
        # summaries
        latent_means <- do.call(rbind, 
                                list(haz_urban = t(apply(fitted_haz_urban_mat, 1, quantile, c(0.025, 0.1, 0.5, 0.9, 0.975))),
                                     haz_rural = t(apply(fitted_haz_rural_mat, 1, quantile, c(0.025, 0.1, 0.5, 0.9, 0.975))),
                                     waz_urban = t(apply(fitted_waz_urban_mat, 1, quantile, c(0.025, 0.1, 0.5, 0.9, 0.975))),
                                     waz_rural = t(apply(fitted_waz_rural_mat, 1, quantile, c(0.025, 0.1, 0.5, 0.9, 0.975)))))
        latent_means <- as_tibble(latent_means) %>% setNames(paste0("est_", names(.)))

        latent_mean_vars <- c(apply(fitted_haz_urban_mat, 1, var),
                              apply(fitted_haz_rural_mat, 1, var),
                              apply(fitted_waz_urban_mat, 1, var),
                              apply(fitted_waz_rural_mat, 1, var))

        latent_means <- bind_cols(admin1 = rep(1:n_regions, 4),
                                  rural = rep(rep(0:1, each = n_regions), 2),
                                  outcome = c(rep("HAZ", n_regions * 2), rep("WAZ", n_regions * 2)),
                                  latent_means, 
                                  var = latent_mean_vars)
        
        # save results
        results[[ model_names[i] ]]$fit <- tmp
        results[[ model_names[i] ]]$latent_means <- latent_means
        
        message("Done!")
        
    }
    
    return(results)
}
