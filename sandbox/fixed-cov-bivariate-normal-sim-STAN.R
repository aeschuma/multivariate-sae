# STAN modeling bivariate normal
rm(list = ls())

# libraries
library(rstan)
library(tidyverse)
library(haven)
library(magrittr)
library(mvtnorm)
library(Matrix)
library(bayesplot)

# directories
savedir <- "../../../Dropbox/dissertation_2/survey-csmf/sandbox"

# stan options
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

# set seed
set.seed(8008135)

# parameters
nreg <- 8
nt <- 5
ncause <- 2
N <- nreg * nt * ncause
alpha <- 1:ncause
rhos <- runif(nreg*nt, -0.2, 0.2)
Sigma <- vector(mode = "list", length = nreg*nt)
for (i in 1:(nreg*nt)) {
    R <- matrix(c(1, rhos[i], rhos[i], 1), ncol = ncause)
    sigs <- sqrt(runif(ncause, 0.05, 0.1))
    Sigma[[i]] <- diag(sigs) %*% R %*% diag(sigs)
}
gamma_var <- 1

# simulation settings and results storage
nsim <- 1
alpha.res <- matrix(NA, nrow = nsim, ncol = ncause)
coverage.res <- matrix(NA, nrow = nsim, ncol = ncause)
covprec.res <- rep(NA, nsim)
reprec.res <- rep(NA, nsim)
coverage.re.sd.res <- rep(NA, nsim)

# make block diagnonal covariance matrix
bdiag_matlist_prec <- list()
bdiag_matlist <- list()
Sigma.array <- array(NA, dim = c(nt*nreg, ncause, ncause))
for(i in 1:(nreg*nt)) { 
    bdiag_matlist_prec[[i]] <- solve(Sigma[[i]]) 
    bdiag_matlist[[i]] <- Sigma[[i]]
    Sigma.array[i,,] <- Sigma[[i]]
}
block_prec <- bdiag(bdiag_matlist_prec)
block_var <- bdiag(bdiag_matlist)

# simulate epsilon
epsilon <- rmvnorm(nsim, mean = rep(0, ncause * nreg*nt), sigma = as.matrix(block_var))

# simulate gammas
gamma <- rmvnorm(nsim, mean = rep(0, nreg), sigma = diag(rep(gamma_var, nreg)))
# for (s in 1:nsim) {
#     gamma[s, ] <- gamma[s, ] - mean(gamma[s, ]) # sum to zero
# }

# for (s in 1:nsim) {
    # testing 
    s <- 1
    
    cat(paste0("Starting sim ", s, "....\n"))
    # simulate outcome
    all_alpha <- rep(alpha, nreg*nt)
    all_gamma <- rep(rep(gamma[s, ], each = ncause), nt)
    all_epsilon <- epsilon[s, ]
    y <- all_alpha + all_gamma + all_epsilon
    
    dat <- data.frame(reg = rep(rep(1:nreg, each = ncause), nt),
                      cause = rep(1:ncause, nreg*nt),
                      time = rep(1:nt, each = ncause*nreg),
                      y = y)
    
    # create dataframe
    datlist <- list(nreg = nreg,
                    ncause = ncause,
                    N = nreg * nt,
                    region = dat$reg[seq(1, length(y), by = ncause)],
                    Sigma = Sigma.array,
                    y = cbind(dat$y[which(dat$cause == 1)],
                              dat$y[which(dat$cause == 2)]))
    niter <- 1000
    nchains <- 2
    nthin <- 1
    prop_warmup <- 0.5
    inits <- NULL
    start.time <- proc.time() 
    mod_stan <- stan(file = "fixed-cov-bivariate-normal.stan",
                     data = datlist,
                     iter = niter, chains = nchains, thin=nthin,
                     warmup = niter*prop_warmup,
                     # init = inits,
                     control = list(max_treedepth = 15))
    stop.time <- proc.time()
    elapsed.time <- stop.time[3] - start.time[3]
    
    # save model
    write_rds(mod_stan, paste0(savedir, "/fixed_cov_bivariate_normal_stan.rds"))
    
    # diagnostic plots
    stan_trace(mod_stan, pars = c("sigma", "gamma", "beta"))
    mcmc_scatter(
        as.matrix(mod_stan),
        pars = c("beta[1]", "beta[2]"),
        np = nuts_params(mod_stan),
        np_style = scatter_style_np(div_color = "green", div_alpha = 0.8)
    )
    pairs(mod_stan, pars = c("sigma", "beta"))
    mcmc_areas(as.matrix(mod_stan),
               pars = c("beta[1]", "beta[2]", "sigma"),
               prob = 0.8)
    mcmc_nuts_energy(nuts_params(mod_stan))
# }
    
    
    
    
    
# testing normal (not bivariate)
    # sigma_ep <- exp(rnorm(N, log(0.1), 0.5))
    # ep <- rnorm(N, rep(0, N), sigma_ep)
    # dat2 <- data.frame(reg = rep(rep(1:nreg, each = ncause), nt),
    #                    cause = rep(1:ncause, nreg*nt),
    #                    time = rep(1:nt, each = ncause*nreg),
    #                    y = all_alpha + all_gamma + ep)
    # X_beta <- model.matrix(~ -1 + factor(cause), data = dat2)
    # Z_gamma <- model.matrix(~ -1 + factor(reg), data = dat2)
    # 
    # # create data frame
    # datlist2 <- list(N = N,
    #                  nreg = nreg,
    #                  ncause = ncause,
    #                  Sigma = sigma_ep,
    #                  y = dat2$y,
    #                  X = X_beta,
    #                  Z = Z_gamma)
    # niter <- 2000
    # nchains <- 2
    # nthin <- 1
    # prop_warmup <- 0.5
    # inits <- NULL
    # start.time <- proc.time() 
    # mod_stan <- stan(file = "fixed-var-normal.stan",
    #                  data = datlist2,
    #                  iter = niter, chains = nchains, thin=nthin,
    #                  warmup = niter*prop_warmup,
    #                  # init = inits,
    #                  control = list(max_treedepth = 15))
    # stop.time <- proc.time()
    # elapsed.time <- stop.time[3] - start.time[3]
    # 
    # # diagnostic plots
    # stan_trace(mod_stan, pars = c("sigma", "gamma", "beta"))
    # mcmc_scatter(
    #     as.matrix(mod_stan),
    #     pars = c("beta[1]", "beta[2]"),
    #     np = nuts_params(mod_stan),
    #     np_style = scatter_style_np(div_color = "green", div_alpha = 0.8)
    # )
    # pairs(mod_stan, pars = c("sigma", "beta"))
    # mcmc_areas(as.matrix(mod_stan),
    #            pars = c("beta[1]", "beta[2]", "sigma"),
    #            prob = 0.8)
    # mcmc_nuts_energy(nuts_params(mod_stan))
    # 