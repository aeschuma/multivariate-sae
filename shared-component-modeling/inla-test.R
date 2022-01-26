# Austin Schumacher
# 1/1/2022
# Testing INLA code for shared component modeling from Leigh's paper with Jon

rm(list = ls())

## Libraries ####
library(SUMMER)
# help(package = "SUMMER", help_type = "html")
library(classInt)
library(RColorBrewer)
library(dplyr)
library(tidyr)
library(rgdal)
library(scales)
library(survey)
library(tidyverse)
library(spdep)
library(INLA)
library(knitr)
library(rstan)
library(cmdstanr)
library(GGally)

set.seed(858)

# Load data ####
load("/Users/austin/Dropbox/dissertation_2/survey-csmf/data/ken_dhs2014/data/haz-waz-kenDHS2014.rda")

## Stage 1 direct: univariate and bivariate ####

my.svydesign <- survey::svydesign(ids = ~ cluster,
                                  strata = ~ strata, nest = T, weights = ~weights,
                                  data = dat)

n_regions <- length(unique(dat$admin1.char))
admin1v <- unique(dat$admin1.char)
results <- data.frame(admin1 = unique(dat$admin1),
                      admin1.name = unique(dat$admin1.name),
                      admin1.char = admin1v,
                      meanHAZ.bi = rep(NA, n_regions),
                      meanWAZ.bi = rep(NA, n_regions),
                      seHAZ.bi = rep(NA, n_regions),
                      seWAZ.bi = rep(NA, n_regions),
                      corr.bi = rep(NA, n_regions))
V.list <- vector(mode = "list", length = n_regions)
V.array <- array(NA, dim = c(nrow(results), 2, 2))
names(V.list) <- admin1v

V.array.2 <- array(0, dim = c(nrow(results), 2, 2))

for(i in 1:n_regions) {
    admin1.tmp <- admin1v[i]
    
    tmp <- subset(my.svydesign, admin1.char == admin1.tmp)
    means.svymean <- svymean(~HAZ + WAZ, tmp)
    
    index.tmp <- results$admin1.char == admin1.tmp
    
    results$meanHAZ.bi[index.tmp] <- means.svymean[["HAZ"]]
    results$meanWAZ.bi[index.tmp] <- means.svymean[["WAZ"]]
    
    V.tmp <- vcov(means.svymean)
    V.list[[admin1.tmp]] <- V.tmp
    V.array[i, , ] <- V.tmp
    tmat <- matrix(0, nrow = 2, ncol = 2)
    diag(tmat) <- diag(V.tmp)
    V.array.2[i, , ] <- tmat
    results$seHAZ.bi[index.tmp] <- V.tmp[1, 1]^0.5
    results$seWAZ.bi[index.tmp] <- V.tmp[2, 2]^0.5
    
    D <- diag(sqrt(diag(V.tmp)))
    DInv <- solve(D)
    corr.tmp <- DInv %*% V.tmp %*% DInv
    
    results$corr.bi[index.tmp] <- corr.tmp[1, 2]
    
    # test recovery of vcov
    # omega <- matrix(c(1, corr.tmp[1, 2], corr.tmp[1, 2], 1), nrow = 2, ncol = 2)
    # D <- diag(c(V.tmp[1, 1]^0.5, V.tmp[2, 2]^0.5))
    # V.test <- D %*% omega %*% D
}

# reformat data
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

# create block diagonal matrix from list of fixed covariances
diag.list <- vector(mode = "list", length = length(V.list))
for (i in 1:length(V.list))  {
    diag.list[[i]] <- solve(V.list[[i]])
    D <- diag(sqrt(diag(V.list[[i]])))
    DInv <- solve(D)
}
Vdes_prec <- bdiag(diag.list)

# "Univariate" models ####
# diagonal V matrix!

# Old way; nonshared iid ####

# priors
fe.prec <- list(prec.intercept = 0,
                prec = 5)
iid.prior <- list(theta = list(prior = "loggamma",
                               param = c(1, 5e-05)))
cov_prior <- list(prec = list(prior = "loggamma", param = c(10000, 10000)))

# model formula
m.form <- value ~ -1 + outcome + 
    f(admin1.haz, model = 'iid',
      hyper = iid.prior) +
    f(admin1.waz, model = 'iid',
      hyper = iid.prior) +
    f(obs,  model='generic0', Cmatrix = Vdes_prec,
      hyper = cov_prior)

# fit model
smooth.direct.bi <- inla(m.form,
                         data = results.long,
                         family = "gaussian",
                         control.fixed = fe.prec,
                         control.predictor=list(compute=TRUE),
                         control.compute=list(config = TRUE),
                         control.family = list(hyper = list(prec = list(initial = log(1), fixed=TRUE))),
                         scale = 1000000)
smooth.direct.bi$summary.fixed
smooth.direct.bi$summary.hyperpar

# leigh's way: nonshared iid ####
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

## we now have to add one 'iid2d' model for each observation pair,
## since their cov.matrix is different. we have to do this
## automatically... here I add numbers directly for simplicity
add <- ""
for(j in 1:Nt) {
    corr <-  results$corr.bi[j]
    init.prec.haz <-  log(1/results$seHAZ.bi[j]^2)
    init.prec.waz <-  log(1/results$seWAZ.bi[j]^2)
    
    add <-  paste(add, paste(" + 
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

## Model formula
formula.iid2d <-  formula(paste("value ~ -1 + outcome + f(admin1.haz, model = 'iid', hyper = iid.prior) + f(admin1.waz, model = 'iid', hyper = iid.prior) + ",
                        paste(add, collapse = " ")))

# linear combinatin of preds without 2x2 REs
lc.vec.haz.fe <- c(rep(1, n_regions))
lc.vec.waz.fe <- c(rep(1, n_regions))
diag.na.mat <- matrix(NA, nrow = n_regions, ncol = n_regions)
diag(diag.na.mat) <- 1
lc.vec.admin1.haz.re <- diag.na.mat
lc.vec.admin1.waz.re <- diag.na.mat

lc.all.haz <- inla.make.lincombs(outcomeHAZ = lc.vec.haz.fe, 
                                 admin1.haz = lc.vec.admin1.haz.re)
lc.all.waz <- inla.make.lincombs(outcomeWAZ = lc.vec.waz.fe, 
                                 admin1.waz = lc.vec.admin1.waz.re)
names(lc.all.haz) <- gsub("lc", "haz.reg.", names(lc.all.haz))
names(lc.all.waz) <- gsub("lc", "waz.reg.", names(lc.all.waz))

## Fit model
mod.leigh <- inla(formula.iid2d, data = data, 
          family = "gaussian",
          lincomb = c(lc.all.haz, lc.all.waz),
          control.family = list(hyper = list(prec = list(initial = 10, fixed=T))), 
          control.predictor = list(compute=T),
          control.compute = list(config=T),
          control.inla = list(lincomb.derived.correlation.matrix=T),
          control.fixed = list(prec = list(default = 0.001), correlation.matrix=T) )

mod.leigh.2 <- inla(formula.iid2d, data = data, 
                    family = "gaussian",
                    lincomb = c(lc.all.haz, lc.all.waz),
                    control.family = list(hyper = list(prec = list(initial = log(1), fixed=T))), 
                    control.predictor = list(compute=TRUE),
                    control.compute = list(config=TRUE),
                    control.inla = list(lincomb.derived.correlation.matrix=T),
                    control.fixed = fe.prec,
                    scale = 1000000)

# stan nonshared iid ####
sigma_normal_sd <- 1

datlist <- list(R = n_regions,
                regions = results$admin1,
                Sigma = V.array,
                y = results[, c("meanHAZ.bi", "meanWAZ.bi")],
                sigma_normal_sd = sigma_normal_sd)

stan_file <- "../ken2014-hazwaz-modeling/stan-models/nonshared-iid.stan"

## stan options
niter <- 10000
nchains <- 3
prop_warmup <- 0.5
max_treedepth <- 100
adapt_delta <- 0.8
nthin <- 1
options(mc.cores = parallel::detectCores())

cmd_mod <- cmdstan_model(stan_file = stan_file)
fit <- cmd_mod$sample(data = datlist,
                      iter_warmup = niter*prop_warmup, iter_sampling = niter*(1-prop_warmup),
                      chains = nchains, thin = nthin,
                      adapt_delta = adapt_delta, max_treedepth = max_treedepth,
                      refresh =  0.1 * niter)
mod.stan.bi.nonshared.iid <- rstan::read_stan_csv(fit$output_files())

stan_diag(mod.stan.bi.nonshared.iid)
stan_diags <- data.frame(pct_divergent = get_num_divergent(mod.stan.bi.nonshared.iid)/(niter * nchains * prop_warmup),
                         pct_max_tree_exceeded = get_num_max_treedepth(mod.stan.bi.nonshared.iid)/(niter * prop_warmup * nchains),
                         pct_bmfi_low_chains = sum(is.finite(get_low_bfmi_chains(mod.stan.bi.nonshared.iid)))/nchains)
stan_diags
traceplot(mod.stan.bi.nonshared.iid, pars = c("beta", "sigma"))

params_to_extract <- c("beta", "sigma", "preds")
re_params <- c("v_1", "v_2", "u_1", "u_2")

# compare results
mod.stan.summary.bi.nonshared.iid <- summary(mod.stan.bi.nonshared.iid,
                                         pars = params_to_extract[params_to_extract != "preds"],
                                         probs = c(0.025, 0.5, 0.975))

kable(mod.stan.summary.bi.nonshared.iid$summary, format = "markdown", digits = 3)

# compare
mod.leigh$summary.fixed
1/sqrt(mod.leigh$summary.hyperpar[, c("0.025quant", "0.5quant", "0.975quant")])

mod.leigh.2$summary.fixed
mod.leigh.2$summary.hyperpar

smooth.direct.bi$summary.fixed
smooth.direct.bi$summary.hyperpar

# compare random effects
stan.re.summary <- summary(mod.stan.bi.nonshared.iid, pars = c("gamma_1", "gamma_2"), probs = c(0.025, 0.5, 0.975))$summary[,c("2.5%", "50%", "97.5%")]
inla.re.summary <- rbind(mod.leigh$summary.random[c("admin1.haz")][[1]],
                         mod.leigh$summary.random[c("admin1.waz")][[1]])
inla2.re.summary <- rbind(mod.leigh.2$summary.random[c("admin1.haz")][[1]],
                          mod.leigh.2$summary.random[c("admin1.waz")][[1]])

res.median.re <- tibble(inla = inla.re.summary$`0.5quant`,
                        inla.2 = inla2.re.summary$`0.5quant`,
                        stan = stan.re.summary[,"50%"])

pairs.med <- ggpairs(res.median.re) +
  theme_light()
for (i in 2:pairs.med$nrow) {
  for (j in 1:(i-1)) {
    pairs.med[i,j] <- pairs.med[i,j] + geom_abline(intercept = 0,slope = 1, col = "darkgreen")
  }
}
pairs.med

res.width <- tibble(inla = inla.re.summary$`0.975quant` - inla.re.summary$`0.025quant`,
                    inla.2 = inla2.re.summary$`0.975quant` - inla2.re.summary$`0.025quant`,
                    stan = stan.re.summary[,"97.5%"] - stan.re.summary[,"2.5%"])

pairs.width <- ggpairs(res.width) +
  theme_light()
for (i in 2:pairs.width$nrow) {
  for (j in 1:(i-1)) {
    pairs.width[i,j] <- pairs.width[i,j] + geom_abline(intercept = 0,slope = 1, col = "darkgreen")
  }
}
pairs.width

# compare predictions
stan.preds.summary <- summary(mod.stan.bi.nonshared.iid, pars = "preds", probs = c(0.025, 0.5, 0.975))$summary[,c("2.5%", "50%", "97.5%")]
stan.preds.summary <- stan.preds.summary[c(seq(1, (n_regions*2), by = 2), seq(2, (n_regions*2), by = 2)),]
mod.leigh.preds <- mod.leigh$summary.lincomb.derived[c(names(lc.all.haz), names(lc.all.waz)),]
mod.leigh.2.preds <- mod.leigh.2$summary.lincomb.derived[c(names(lc.all.haz), names(lc.all.waz)),]

res.median <- tibble(inla = mod.leigh.preds$`0.5quant`,
                     inla.2 = mod.leigh.2.preds$`0.5quant`,
                     stan = stan.preds.summary[, "50%"])

pairs.med <- ggpairs(res.median) +
  theme_light()
for (i in 2:pairs.med$nrow) {
  for (j in 1:(i-1)) {
    pairs.med[i,j] <- pairs.med[i,j] + geom_abline(intercept = 0,slope = 1, col = "darkgreen")
  }
}
pairs.med

res.width <- tibble(inla = mod.leigh.preds$`0.975quant` - mod.leigh.preds$`0.025quant`,
                    inla.2 = mod.leigh.2.preds$`0.975quant` - mod.leigh.2.preds$`0.025quant`,
                    stan = stan.preds.summary[, "97.5%"] - stan.preds.summary[, "2.5%"])

pairs.width <- ggpairs(res.width) +
  theme_light()
for (i in 2:pairs.width$nrow) {
  for (j in 1:(i-1)) {
    pairs.width[i,j] <- pairs.width[i,j] + geom_abline(intercept = 0,slope = 1, col = "darkgreen")
  }
}
pairs.width

# nonshared BYM2 old way ####

# priors
bym2_prior <- list(phi=list(prior="logitbeta", param=c(1, 1), initial=0.5), 
                   prec=list(prior="pc.prec", param=c(1, 0.01), initial=5))

# model formula
m.form <- value ~ -1 + outcome + 
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
  f(obs,  model='generic0', Cmatrix = Vdes_prec,
    hyper = cov_prior)

# fit model
smooth.direct.bi.nonshared <- inla(m.form,
                                data = results.long,
                                family = "gaussian",
                                control.fixed = fe.prec,
                                control.predictor=list(compute=TRUE),
                                control.compute=list(config = TRUE),
                                control.family = list(hyper = list(prec = list(initial = log(1), fixed=TRUE))),
                                scale = 1000000)
smooth.direct.bi.nonshared$summary.fixed
smooth.direct.bi.nonshared$summary.hyperpar

# nonshared bym2 leigh's way ####
data$admin1.haz.iid <- data$admin1.haz
data$admin1.waz.iid <- data$admin1.waz

## Model formula
formula.iid2d <-  formula(paste("value ~ -1 + outcome + f(admin1.haz, model = 'bym2',
    graph = admin1.mat, 
    scale.model = T, 
    constr = T,
    hyper = bym2_prior) +
  f(admin1.waz, model = 'bym2',
    graph = admin1.mat, 
    scale.model = T, 
    constr = T,
    hyper = bym2_prior)",
                                paste(add, collapse = " ")))

# linear combinatin of preds without 2x2 REs
lc.vec.haz.fe <- c(rep(1, n_regions))
lc.vec.waz.fe <- c(rep(1, n_regions))
diag.na.mat <- matrix(NA, nrow = n_regions, ncol = n_regions)
diag(diag.na.mat) <- 1
lc.vec.admin1.haz.re <- diag.na.mat
lc.vec.admin1.waz.re <- diag.na.mat

lc.all.haz <- inla.make.lincombs(outcomeHAZ = lc.vec.haz.fe, 
                                 admin1.haz = lc.vec.admin1.haz.re)
lc.all.waz <- inla.make.lincombs(outcomeWAZ = lc.vec.waz.fe, 
                                 admin1.waz = lc.vec.admin1.waz.re)
names(lc.all.haz) <- gsub("lc", "haz.reg.", names(lc.all.haz))
names(lc.all.waz) <- gsub("lc", "waz.reg.", names(lc.all.waz))


## Fit model
mod.leigh.nonshared <- inla(formula.iid2d, data = data, 
                  family = "gaussian",
                  lincomb = c(lc.all.haz, lc.all.waz),
                  control.inla = list(lincomb.derived.correlation.matrix=T),
                  control.family = list(hyper = list(prec = list(initial = 10, fixed=T))), 
                  control.predictor = list(compute=T),
                  control.compute = list(config=T),
                  control.fixed = list(prec = list(default = 0.001), correlation.matrix=T) )

mod.leigh.2.nonshared <- inla(formula.iid2d, data = data, 
                    family = "gaussian",
                    lincomb = c(lc.all.haz, lc.all.waz),
                    control.inla = list(lincomb.derived.correlation.matrix=T),
                    control.family = list(hyper = list(prec = list(initial = log(1), fixed=T))), 
                    control.predictor = list(compute=TRUE),
                    control.compute = list(config=TRUE),
                    control.fixed = fe.prec,
                    scale = 1000000)

besag_prior <- list(prec=list(prior="pc.prec", param=c(1, 0.01), initial=5))

formula.iid2d.alt <- formula(paste("value ~ -1 + outcome + 
  f(admin1.haz, model = 'besag',
    graph = admin1.mat, 
    scale.model = T, 
    constr = T,
    hyper = besag_prior) +
  f(admin1.haz.iid, model = 'iid') +
  f(admin1.waz, model = 'besag',
    graph = admin1.mat, 
    scale.model = T, 
    constr = T,
    hyper = besag_prior) +
  f(admin1.waz.iid, model = 'iid') ",
                                   paste(add, collapse = " ")))

# adapted lincombs
lc.all.haz.2 <- inla.make.lincombs(outcomeHAZ = lc.vec.haz.fe, 
                                   admin1.haz = lc.vec.admin1.haz.re,
                                   admin1.haz.iid = lc.vec.admin1.haz.re)
lc.all.waz.2 <- inla.make.lincombs(outcomeWAZ = lc.vec.waz.fe, 
                                   admin1.waz = lc.vec.admin1.waz.re,
                                   admin1.waz.iid = lc.vec.admin1.waz.re)
names(lc.all.haz.2) <- gsub("lc", "haz.reg.", names(lc.all.haz))
names(lc.all.waz.2) <- gsub("lc", "waz.reg.", names(lc.all.waz))

mod.leigh.3.nonshared <- inla(formula.iid2d.alt, data = data, 
                           family = "gaussian",
                           lincomb = c(lc.all.haz, lc.all.waz),
                           control.family = list(hyper = list(prec = list(initial = 10, fixed=T))), 
                           control.predictor = list(compute=T),
                           control.compute = list(config=T),
                           control.inla = list(lincomb.derived.correlation.matrix=T),
                           control.fixed = list(prec = list(default = 0.001), correlation.matrix=T) )

# STAN nonshared bym2 to compare ####
rho_beta_1 <- 1
rho_beta_2 <- 1
sigma_normal_sd <- 1

datlist <- list(R = n_regions,
                regions = results$admin1,
                Sigma = V.array,
                y = results[, c("meanHAZ.bi", "meanWAZ.bi")],
                N_edges = node.info$n_edges,
                node1 = node.info$node1,
                node2 = node.info$node2,
                scaling_factor = scaling_factor,
                rho_beta_1 = rho_beta_1,
                rho_beta_2 = rho_beta_2,
                sigma_normal_sd = sigma_normal_sd)

stan_file <- "../ken2014-hazwaz-modeling/stan-models/nonshared-bym2.stan"

## stan options
niter <- 10000
nchains <- 3
prop_warmup <- 0.5
max_treedepth <- 100
adapt_delta <- 0.8
nthin <- 1
options(mc.cores = parallel::detectCores())

cmd_mod <- cmdstan_model(stan_file = stan_file)
fit <- cmd_mod$sample(data = datlist,
                      iter_warmup = niter*prop_warmup, iter_sampling = niter*(1-prop_warmup),
                      chains = nchains, thin = nthin,
                      adapt_delta = adapt_delta, max_treedepth = max_treedepth,
                      refresh =  0.1 * niter)
mod.stan.bi.nonshared <- rstan::read_stan_csv(fit$output_files())
params_to_extract <- c("beta", "sigma", "rho", "preds")
re_params <- c("v_1", "v_2", "u_1", "u_2")

# compare results
mod.stan.summary.bi.nonshared <- summary(mod.stan.bi.nonshared,
                                      pars = params_to_extract[params_to_extract != "preds"],
                                      probs = c(0.025, 0.5, 0.975))
kable(mod.stan.summary.bi.nonshared$summary, format = "markdown", digits = 3)

mod.leigh.nonshared$summary.fixed
mod.leigh.nonshared$summary.hyperpar

mod.leigh.2.nonshared$summary.fixed
mod.leigh.2.nonshared$summary.hyperpar

mod.leigh.3.nonshared$summary.fixed
mod.leigh.3.nonshared$summary.hyperpar

# compare random effects
# stan.re.summary <- summary(mod.stan.bi.nonshared, pars = c("convolved_re_1", "convolved_re_2"), probs = c(0.025, 0.5, 0.975))$summary[,c("2.5%", "50%", "97.5%")]
# inla.re.summary <- rbind(mod.leigh.nonshared$summary.random[c("admin1.haz")][[1]],
#                          mod.leigh.nonshared$summary.random[c("admin1.waz")][[1]])
# inla2.re.summary <- rbind(mod.leigh.2.nonshared$summary.random[c("admin1.haz")][[1]],
#                           mod.leigh.2.nonshared$summary.random[c("admin1.waz")][[1]])
# 
# res.median.re <- tibble(inla = inla.re.summary$`0.5quant`,
#                         inla.2 = inla2.re.summary$`0.5quant`,
#                         stan = stan.re.summary[,"50%"])
# 
# pairs.med <- ggpairs(res.median.re) +
#   theme_light()
# for (i in 2:pairs.med$nrow) {
#   for (j in 1:(i-1)) {
#     pairs.med[i,j] <- pairs.med[i,j] + geom_abline(intercept = 0,slope = 1, col = "darkgreen")
#   }
# }
# pairs.med
# 
# res.width <- tibble(inla = inla.re.summary$`0.975quant` - inla.re.summary$`0.025quant`,
#                     inla.2 = inla2.re.summary$`0.975quant` - inla2.re.summary$`0.025quant`,
#                     stan = stan.re.summary[,"97.5%"] - stan.re.summary[,"2.5%"])
# 
# pairs.width <- ggpairs(res.width) +
#   theme_light()
# for (i in 2:pairs.width$nrow) {
#   for (j in 1:(i-1)) {
#     pairs.width[i,j] <- pairs.width[i,j] + geom_abline(intercept = 0,slope = 1, col = "darkgreen")
#   }
# }
# pairs.width

# compare predictions
stan.preds.summary <- summary(mod.stan.bi.nonshared, pars = "preds", probs = c(0.025, 0.5, 0.975))$summary[,c("2.5%", "50%", "97.5%")]
stan.preds.summary <- stan.preds.summary[c(seq(1, (n_regions*2), by = 2), seq(2, (n_regions*2), by = 2)),]
mod.leigh.preds <- mod.leigh.nonshared$summary.lincomb.derived[c(names(lc.all.haz), names(lc.all.waz)),]
mod.leigh.2.preds <- mod.leigh.2.nonshared$summary.lincomb.derived[c(names(lc.all.haz), names(lc.all.waz)),]

res.median <- tibble(inla = mod.leigh.preds$`0.5quant`,
                     inla.2 = mod.leigh.2.preds$`0.5quant`,
                     stan = stan.preds.summary[, "50%"])

pairs.med <- ggpairs(res.median) +
  theme_light()
for (i in 2:pairs.med$nrow) {
  for (j in 1:(i-1)) {
    pairs.med[i,j] <- pairs.med[i,j] + geom_abline(intercept = 0,slope = 1, col = "darkgreen")
  }
}
pairs.med

res.width <- tibble(inla = mod.leigh.preds$`0.975quant` - mod.leigh.preds$`0.025quant`,
                    inla.2 = mod.leigh.2.preds$`0.975quant` - mod.leigh.2.preds$`0.025quant`,
                    stan = stan.preds.summary[, "97.5%"] - stan.preds.summary[, "2.5%"])

pairs.width <- ggpairs(res.width) +
  theme_light()
for (i in 2:pairs.width$nrow) {
  for (j in 1:(i-1)) {
    pairs.width[i,j] <- pairs.width[i,j] + geom_abline(intercept = 0,slope = 1, col = "darkgreen")
  }
}
pairs.width

# shared bym2 leigh ####
data$admin1.haz.2 <- data$admin1.haz

lambda_prior <- list(beta = list(prior = 'logtnormal', param = c(0, 1)))
iid_prior <- list(prec=list(prior="pc.prec", param=c(0.1, 0.1), initial=5))

## Model formula
formula.iid2d <-  formula(paste("value ~ -1 + outcome + f(admin1.haz, model = 'bym2',
    graph = admin1.mat, 
    scale.model = T, 
    constr = T,
    hyper = bym2_prior) +
  f(admin1.waz, model = 'bym2',
    graph = admin1.mat, 
    scale.model = T, 
    constr = T,
    hyper = bym2_prior) +
  f(admin1.haz.2, copy = \"admin1.waz\", fixed = FALSE, hyper = lambda_prior)",
                                paste(add, collapse = " ")))

# linear combinatin of preds without 2x2 REs
lc.vec.haz.fe <- c(rep(1, n_regions))
lc.vec.waz.fe <- c(rep(1, n_regions))
diag.na.mat <- matrix(NA, nrow = n_regions, ncol = n_regions)
diag(diag.na.mat) <- 1
lc.vec.admin1.haz.re <- diag.na.mat
lc.vec.admin1.waz.re <- diag.na.mat

lc.all.haz <- inla.make.lincombs(outcomeHAZ = lc.vec.haz.fe, 
                                 admin1.haz = lc.vec.admin1.haz.re,
                                 admin1.haz.2 = lc.vec.admin1.haz.re)
lc.all.waz <- inla.make.lincombs(outcomeWAZ = lc.vec.waz.fe, 
                                 admin1.waz = lc.vec.admin1.waz.re)

names(lc.all.haz) <- gsub("lc", "haz.reg.", names(lc.all.haz))
names(lc.all.waz) <- gsub("lc", "waz.reg.", names(lc.all.waz))

## Fit model
mod.leigh.shared <- inla(formula.iid2d, data = data, 
                         family = "gaussian",
                         lincomb = c(lc.all.haz, lc.all.waz),
                         control.inla = list(lincomb.derived.correlation.matrix=T,
                                             strategy = "laplace", npoints = 21),
                         control.family = list(hyper = list(prec = list(initial = 10, fixed=T))), 
                         control.predictor = list(compute=T),
                         control.compute = list(config=T, 
                                                waic = TRUE, 
                                                cpo = TRUE, 
                                                dic = TRUE),
                         control.fixed = list(prec = list(default = 0.001), correlation.matrix=T) )

mod.leigh.2.shared <- inla(formula.iid2d, data = data, 
                           family = "gaussian",
                           lincomb = c(lc.all.haz, lc.all.waz),
                           control.inla = list(lincomb.derived.correlation.matrix=T),
                           control.family = list(hyper = list(prec = list(initial = log(1), fixed=T))), 
                           control.predictor = list(compute=TRUE),
                           control.compute = list(config=TRUE, waic = TRUE, cpo = TRUE, dic = TRUE),
                           control.fixed = fe.prec,
                           scale = 1000000)

formula.iid2d.alt <- formula(paste("value ~ -1 + outcome + 
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
  f(admin1.haz.2, copy = \"admin1.waz\", fixed = FALSE, hyper = lambda_prior)" ,
                                   paste(add, collapse = " ")))

# adapted lincombs
lc.all.haz.2 <- inla.make.lincombs(outcomeHAZ = lc.vec.haz.fe, 
                                   admin1.haz = lc.vec.admin1.haz.re,
                                   admin1.haz.iid = lc.vec.admin1.haz.re,
                                   admin1.haz.2 = lc.vec.admin1.haz.re)
lc.all.waz.2 <- inla.make.lincombs(outcomeWAZ = lc.vec.waz.fe, 
                                   admin1.waz = lc.vec.admin1.waz.re,
                                   admin1.waz.iid = lc.vec.admin1.waz.re)
names(lc.all.haz.2) <- gsub("lc", "haz.reg.", names(lc.all.haz))
names(lc.all.waz.2) <- gsub("lc", "waz.reg.", names(lc.all.waz))

mod.leigh.3.shared <- inla(formula.iid2d.alt, data = data, 
                           family = "gaussian",
                           lincomb = c(lc.all.haz.2, lc.all.waz.2),
                           control.inla = list(lincomb.derived.correlation.matrix=T),
                           control.family = list(hyper = list(prec = list(initial = 10, fixed=T))), 
                           control.predictor = list(compute=T),
                           control.compute = list(config=T),
                           control.fixed = list(prec = list(default = 0.001), correlation.matrix=T) )

# shared bym2 STAN JUST SHARED SPATIAL####
datlist <- list(R = n_regions,
                regions = results$admin1,
                Sigma = V.array,
                y = results[, c("meanHAZ.bi", "meanWAZ.bi")],
                N_edges = node.info$n_edges,
                node1 = node.info$node1,
                node2 = node.info$node2,
                scaling_factor = scaling_factor,
                rho_beta_1 = rho_beta_1,
                rho_beta_2 = rho_beta_2,
                sigma_normal_sd = sigma_normal_sd)

stan_file <- "../ken2014-hazwaz-modeling/stan-models/shared-bym2-coregionalization.stan"

## stan options
niter <- 10000
nchains <- 3
prop_warmup <- 0.5
max_treedepth <- 100
adapt_delta <- 0.8
nthin <- 1
options(mc.cores = parallel::detectCores())

cmd_mod <- cmdstan_model(stan_file = stan_file)
fit <- cmd_mod$sample(data = datlist,
                      iter_warmup = niter*prop_warmup, iter_sampling = niter*(1-prop_warmup),
                      chains = nchains, thin = nthin,
                      adapt_delta = adapt_delta, max_treedepth = max_treedepth,
                      refresh =  0.1 * niter)
mod.stan.bi.shared <- rstan::read_stan_csv(fit$output_files())
params_to_extract <- c("beta", "sigma", "rho", "lambda")
re_params <- c("v_1", "v_2", "u_1", "u_2")
mod.stan.summary.bi.shared <- summary(mod.stan.bi.shared,
                                      pars = params_to_extract,
                                      probs = c(0.025, 0.5, 0.975))

# shared bym2 STAN SHARED SPATIAL + IID####
stan_file <- "../ken2014-hazwaz-modeling/stan-models/shared-bym2-coregionalization-iidplusbesagshared.stan"

cmd_mod <- cmdstan_model(stan_file = stan_file)
fit <- cmd_mod$sample(data = datlist,
                      iter_warmup = niter*prop_warmup, iter_sampling = niter*(1-prop_warmup),
                      chains = nchains, thin = nthin,
                      adapt_delta = adapt_delta, max_treedepth = max_treedepth,
                      refresh =  0.1 * niter)
mod.stan.bi.shared.iisplusbesag <- rstan::read_stan_csv(fit$output_files())
mod.stan.summary.bi.shared.iisplusbesag <- summary(mod.stan.bi.shared.iisplusbesag,
                                      pars = params_to_extract,
                                      probs = c(0.025, 0.5, 0.975))
# compare shared results ####
kable(mod.stan.summary.bi.shared$summary, format = "markdown", digits = 3)
kable(mod.stan.summary.bi.shared.iisplusbesag$summary, format = "markdown", digits = 3)

mod.leigh.shared$summary.fixed
ml1.hyper <- mod.leigh.shared$summary.hyperpar[, c("0.025quant", "0.5quant", "0.975quant")]
ml1.hyper[grepl("Precision", rownames(ml1.hyper)),] <- 1/sqrt(ml1.hyper[grepl("Precision", rownames(ml1.hyper)),])[3:1]
rownames(ml1.hyper) <- gsub("Precision", "SD", rownames(ml1.hyper))
ml1.hyper

mod.leigh.2.shared$summary.fixed
ml2.hyper <- mod.leigh.2.shared$summary.hyperpar[, c("0.025quant", "0.5quant", "0.975quant")]
ml2.hyper[grepl("Precision", rownames(ml2.hyper)),] <- 1/sqrt(ml2.hyper[grepl("Precision", rownames(ml2.hyper)),])[3:1]
rownames(ml2.hyper) <- gsub("Precision", "SD", rownames(ml2.hyper))
ml2.hyper

mod.leigh.3.shared$summary.fixed
ml3.hyper <- mod.leigh.3.shared$summary.hyperpar[, c("0.025quant", "0.5quant", "0.975quant")]
ml3.hyper[grepl("Precision", rownames(ml3.hyper)),] <- 1/sqrt(ml3.hyper[grepl("Precision", rownames(ml3.hyper)),])[3:1]
rownames(ml3.hyper) <- gsub("Precision", "SD", rownames(ml3.hyper))
ml3.hyper

# compare predictions
stan.preds.summary <- summary(mod.stan.bi.shared, pars = "preds", probs = c(0.025, 0.5, 0.975))$summary[,c("2.5%", "50%", "97.5%")]
stan.preds.summary <- stan.preds.summary[c(seq(1, (n_regions*2), by = 2), seq(2, (n_regions*2), by = 2)),]

stan.iidplusbesag.preds.summary <- summary(mod.stan.bi.shared.iisplusbesag, pars = "preds", probs = c(0.025, 0.5, 0.975))$summary[,c("2.5%", "50%", "97.5%")]
stan.iidplusbesag.preds.summary <- stan.iidplusbesag.preds.summary[c(seq(1, (n_regions*2), by = 2), seq(2, (n_regions*2), by = 2)),]

mod.leigh.preds <- mod.leigh.shared$summary.lincomb.derived[c(names(lc.all.haz), names(lc.all.waz)),]
mod.leigh.2.preds <- mod.leigh.2.shared$summary.lincomb.derived[c(names(lc.all.haz), names(lc.all.waz)),]
mod.leigh.3.preds <- mod.leigh.3.shared$summary.lincomb.derived[c(names(lc.all.haz), names(lc.all.waz)),]

res.median <- tibble(inla.bym2.shared.iid.spatial = mod.leigh.preds$`0.5quant`,
                     inla.iid.besag.shared.spatial = mod.leigh.3.preds$`0.5quant`,
                     stan.bym2.shared.spatial = stan.preds.summary[, "50%"],
                     stan.bym2.shared.iid.spatial = stan.iidplusbesag.preds.summary[, "50%"])

pairs.med <- ggpairs(res.median) +
  theme_light()
for (i in 2:pairs.med$nrow) {
  for (j in 1:(i-1)) {
    pairs.med[i,j] <- pairs.med[i,j] + geom_abline(intercept = 0,slope = 1, col = "darkgreen")
  }
}
pairs.med

res.width <- tibble(inla.bym2.shared.iid.spatial = mod.leigh.preds$`0.975quant` - mod.leigh.preds$`0.025quant`,
                    inla.iid.besag.shared.spatial = mod.leigh.3.preds$`0.975quant` - mod.leigh.3.preds$`0.025quant`,
                    stan.bym2.shared.spatial = stan.preds.summary[, "97.5%"] - stan.preds.summary[, "2.5%"],
                    stan.bym2.shared.iid.spatial = stan.iidplusbesag.preds.summary[, "97.5%"] - stan.iidplusbesag.preds.summary[, "2.5%"])

pairs.width <- ggpairs(res.width) +
  theme_light()
for (i in 2:pairs.width$nrow) {
  for (j in 1:(i-1)) {
    pairs.width[i,j] <- pairs.width[i,j] + geom_abline(intercept = 0,slope = 1, col = "darkgreen")
  }
}
pairs.width

mean(abs(res.width$inla.bym2.shared.iid.spatial - res.width$inla.iid.besag.shared.spatial))
sum(res.width$inla.bym2.shared.iid.spatial - res.width$inla.iid.besag.shared.spatial)

# cross validation measures
mod.leigh.shared$cpo$pit
