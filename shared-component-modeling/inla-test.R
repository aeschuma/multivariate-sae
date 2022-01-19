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

set.seed(858)

# Load data ####
load("/Users/austin/Dropbox/dissertation_2/survey-csmf/data/ken_dhs2014/data/haz-waz-kenDHS2014.rda")

# Univariate models ####

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
    
    # FOR DOING THIS BY HAND:
    {
        ## STEP 0: subtract the overall (over the entire region) survey-weighted mean HAZ and WAZ scores
        x.all <- tmp$variables[, c("HAZ", "WAZ")]
        pweights <- tmp$variables$weights
        psum <- sum(pweights)
        average<-colSums(x.all*pweights/psum)
        centered.x.all <-sweep(x.all, 2, average)
        
        ## STEP 1: multiply centered HAZ and WAZ scores by survey weight, divide by sum of all weights in the region (for all observations)
        centered.weighted.x.all <- centered.x.all * pweights/psum
        
        # Loop through strata (it's a tapply that's summed at the end in the original code)
        strata <- unique(tmp$strata$strata)
        output <- vector(mode = "list", length = length(strata))
        for (ss in 1:length(strata)) {
            index <- which(tmp$strata$strata == strata[ss])
            x <- centered.weighted.x.all[index,]
            cluster <- tmp$cluster$cluster[index]
            nPSU <- as.vector(tmp$fpc$sampsize)[index]
            fpc <- NULL
            f<-rep(1,NROW(x))
            scale <- ifelse(nPSU>1, f*nPSU/(nPSU-1), f)
            
            scale<-scale[!duplicated(cluster)]
            
            ## STEP 2a: sum all of the x's (which are HAZ and WAZ multiplied by survey weights and divided by the sum of all) 
            ##              within each cluster cluster.
            ##              This gives one measure of summed HAZ and WAZ for each cluster.
            ## STEP 2b: Then, subtract the stratum-specific mean HAZ and mean WAZ
            
            # 2a
            x<-rowsum(x,cluster)
            
            # 2b
            x <- sweep(x, 2, colMeans(x), "-")
            
            ## STEP 3: takes the crossproduct, which multiplies (off diagonal) or squares (diagonal) the summed and centered 
            ##            weight adjusted HAZ and WAZ score contributions
            ##            Then multiplies each by a scale factor
            ##         The scale factor is: with n = number of sampled clusters in this strata, n/(n-1) 
            ##         -- given no finite population correction (fpc). If there's an fpc, also multiply by (popsize-nPSU)/popsize
            output[[ss]] <- crossprod(as.matrix(x)*sqrt(scale))
            
        }
        nstrat<-length(unique(strata))
        nokstrat<-sum(sapply(output,function(m) !any(is.na(m))))
        V.tmp.2 <- apply(array(unlist(output),c(2,2,length(output))),1:2,sum,na.rm=TRUE)*nstrat/nokstrat
        
        V.array.2[i, , ] <- V.tmp.2
    }
    
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

# Old way; nonshared

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

# leigh's way: nonshared
data <- as.list(results.long)
Nt <- length(unique(results.long$admin1))
N <- nrow(results.long)
## Define the index-vectors ii.1 ii.2 etc, which are the
## index's for the iid2d-model at timepoint 1, 2, ...
for(j in 1:Nt) {
    itmp = numeric(N)
    itmp[] = NA
    itmp[j] = 1
    itmp[j+Nt] = 2
    data = c(list(itmp), data)
    names(data)[1] = paste("ii.", j, sep="")
}

## we now have to add one 'iid2d' model for each observation pair,
## since their cov.matrix is different. we have to do this
## automatically... here I add numbers directly for simplicity
add=""
for(j in 1:Nt) {
    corr = yGt[j, "covlECt.hat"]/sqrt(prod(yGt[j, c("var.lEt.hat", "var.lCt.hat")]))
    init.precE = log(1/yGt[j, "var.lEt.hat"])
    init.precC = log(1/yGt[j, "var.lCt.hat"])
    
    add = paste(add, paste(" + 
                         f(", paste("ii.", j, sep=""), ", weights, model=\"iid2d\", n=2,
                         hyper = list(
                         prec1 = list(
                         initial =", init.precE,", 
                         fixed = TRUE),
                         prec2 = list(
                         initial =", init.precC,", 
                         fixed = TRUE),
                         cor = list(
                         initial = log((1+", corr, ")/(1-", corr, ")), 
                         fixed = TRUE)))"))
    
}



## define hyper parameters
hypers=c(list(prec = list(param = c(10, 1e-2) )), 
         list(prec = list(param = c(10, 1e-2) ))  )

## Model
S1=matrix(1:Nt, ncol=Nt, nrow=1)
formula.lGt = lGt ~ -1 + interE + interC + f(wkE, model="rw2", hyper=hypers[1], scale.model=F, constr=T, extraconstr=list(A=S1, e=0)) + f(wkC, model="rw2", hyper=hypers[2], scale.model=F, constr=T, extraconstr=list(A=S1, e=0)) + wkE2 + wkC2 + tempE + tempC

## Add iid2d pair to formula
formula.lGt=update(formula.lGt, as.formula(paste(". ~ . ", add)))


## Fit model
mod= inla(formula.lGt, data=data, 
          family = "gaussian",
          control.family = list(hyper = list(prec = list(initial = 10, fixed=T))), 
          control.predictor=list(compute=T),
          control.compute=list(config=T),
          control.inla=list(lincomb.derived.correlation.matrix=T),
          control.fixed=list(prec=list(default=0.001), correlation.matrix=T) )

