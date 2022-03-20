#' Function to fit smoothCluster() (as in SUMMER) for a binomial family
#' without the cluster-level random effect or lono correction.
#' Used to compare results to tmb code, and also for 
#' checking to make sure I know what the model is actually doing.
#' 
smooth_cluster <- function(dat, 
                           Amat,
                           graph_path,
                           family = "binomial",
                           cluster = "cluster",
                           y = "died",
                           Ntrials = "total",
                           region = "admin1") {
  
  
  N <- length(unique(binom_df[,region]))
  
  # re-assign column names
  dat$y <- dat[,y]
  dat$cluster <- dat[,cluster]
  dat$Ntrials <- dat[,Ntrials]
  dat$region <- dat[,region]
  
  # make dataset for inla call
  inla_dat <- list(space_struct = as.numeric(factor(dat$region)),
                   Ntrials = dat$Ntrials,
                   y = dat$y)
  
  # make formula
  formula <- y ~ 1 +
  f(space_struct, model = "bym2", graph = graph_path,
    hyper = list(prec = list(prior = "pc.prec", param = c(1, 0.01)), 
                 phi = list(prior = "pc", param = c(0.5, 2/3))),
    constr = TRUE, rankdef = 1, scale.model = TRUE) 
  
  message("Fitting INLA model...")
  
  begin_time <- Sys.time()
  
  # fit INLA model
  result <- inla(formula,
                 data = inla_dat,
                 family = "binomial",
                 Ntrials = Ntrials,
                 control.predictor = list(compute = TRUE),
                 control.compute = list(config = TRUE))
  
  message("Getting posterior samples...")
  
  # get posterior samples
  mod_list <- vector(mode = "list", length = 6)
  mod_list[[6]] <- result
  samp <- inla.posterior.sample(n = 1000, result = mod_list[[6]])
  
  end_time <- Sys.time()
  
  # process hyperpars
  hyperpar_names <- names(samp[[1]]$hyperpar)
  hyperpar_mat <- matrix(0, nrow = length(hyperpar_names), ncol = 1000)
  for (i in 1:1000) {
    hyperpar_mat[,i] <- samp[[i]]$hyperpar
  }
  rownames(hyperpar_mat) <- hyperpar_names
  
  # process fitted values
  int_idx <- rownames(samp[[1]]$latent) %>% str_detect("outcomeHAZ:1") %>% which()
  space_struct_idx <- rownames(samp[[1]]$latent) %>% str_detect("space_struct") %>% which()
  
  int_mat <- matrix(0, nrow = length(age_idx), ncol = 1000)
  space_struct_mat <- matrix(0, nrow = length(space_struct_idx), ncol = 1000)
  
  # fill in sample matrices
  for (i in 1:1000) {
    int_mat[,i] <- samp[[i]]$latent[int_idx]
    space_struct_mat[,i] <- samp[[i]]$latent[space_struct_idx]
  }
  
  # save sample matrices in list for output
  re_list <- list()
  re_list$intercept <- int_mat
  re_list$space_struct <- space_struct_mat
  
  # obtain total effect for space - first half of bym2
  space_tot_mat <- space_struct_mat[1:length(unique(inla_dat$space_struct)),]
  
  # set up fitted mat
  # add in time_unstruct, space_tot, space_time, age
  s_mat <- int_mat[rep(1,length(unique(inla_dat$space_struct))),] +
    space_tot_mat
  
  # return objects
  return(list(fitted_mat = expit(s_mat),
              re_list = re_list,
              param_mat = hyperpar_mat,
              fit = result,
              runtime = end_time - begin_time))
  
}

