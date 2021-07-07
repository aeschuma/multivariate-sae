rm(list=ls())

## set the root depending on operating system
root <- ifelse(Sys.info()[1]=="Darwin","~/",
               ifelse(Sys.info()[1]=="Windows","P:/",
                      ifelse(Sys.info()[1]=="Linux","/home/students/aeschuma/",
                             stop("Unknown operating system"))))

## the following code makes rstan work on the Box server and the cluster
if (root == "P:/") {
    Sys.setenv(HOME="C:/Users/aeschuma",
               R_USER="C:/Users/aeschuma",
               R_LIBS_USER="C:/Users/aeschuma/R_libraries")
    .libPaths("C:/Users/aeschuma/R_libraries")
} else if (root == "/home/students/aeschuma/") {
    Sys.setenv(HOME=root,
               R_USER=root,
               R_LIBS_USER=paste0(root,"R/x86_64-pc-linux-gnu-library/3.6"))
    .libPaths(paste0(root,"R/x86_64-pc-linux-gnu-library/3.6"))
}

## load libraries
if (root == "/home/students/aeschuma/") {
    library(Rcpp,lib.loc=paste0(root,"R/x86_64-pc-linux-gnu-library/3.6"));
    library(StanHeaders,lib.loc=paste0(root,"R/x86_64-pc-linux-gnu-library/3.6")); 
    library(BH,lib.loc = paste0(root,"R/x86_64-pc-linux-gnu-library/3.6"))
    library(rstan,lib.loc = paste0(root,"R/x86_64-pc-linux-gnu-library/3.6"));
} else {
    library(Rcpp); library(StanHeaders); library(BH); library(rstan);library(bayesplot);
}
library(mvtnorm); library(MASS);library(gtools); library(parallel);
library(scales); library(RColorBrewer);library(data.table);
library(ggplot2);

## define directories

# working directory for code
wd <- paste0(root,"Desktop/survey-csmf/stage-2-simulation")

# directory to save results
savedir <- paste0(root,"Dropbox/dissertation_2/survey-csmf/stage-2-simulation")

## set directory
setwd(wd)
setwd(savedir)

# set number of cores to the max possible
options(mc.cores = parallel::detectCores())

number_of_causes <- 2
number_of_regions <- 8
number_of_replications <- 1

## parameters
beta1 <- 1
beta2 <- 2
rho_lower <- -0.2
rho_upper <- 0.2
sigmasq_lower <- 0.05
sigmasq_upper <- 0.5
sigma_gamma1 <- 1.5
sigma_gamma2 <- 2.5
lambda <- 0.5

## STAN options
niter <- 10000
nchains <- 3
prop_warmup <- 0.5
max_treedepth <- 15
adapt_delta <- 0.8

## Simulation number
number_of_sims <- 5

## which model to run
print("models_to_run")
print(as.logical(commandArgs(trailingOnly=TRUE)[1]))

## data generation options
print("number_of_causes")
print(as.logical(commandArgs(trailingOnly=TRUE)[2]))
print("number_of_regions")
print(as.logical(commandArgs(trailingOnly=TRUE)[3]))
print("number_of_replications")
print(as.logical(commandArgs(trailingOnly=TRUE)[4]))

## parameters
print("beta1")
print(as.logical(commandArgs(trailingOnly=TRUE)[5]))
print("beta2")
print(as.logical(commandArgs(trailingOnly=TRUE)[6]))
print("rho_lower")
print(as.logical(commandArgs(trailingOnly=TRUE)[7]))
print("rho_upper")
print(as.logical(commandArgs(trailingOnly=TRUE)[8]))
print("sigmasq_lower")
print(as.logical(commandArgs(trailingOnly=TRUE)[9]))
print("sigmasq_upper")
print(as.numeric(commandArgs(trailingOnly=TRUE)[10]))
print("sigma_gamma1")
print(as.numeric(commandArgs(trailingOnly=TRUE)[11]))
print("sigma_gamma2")
print(as.numeric(commandArgs(trailingOnly=TRUE)[12]))
print("lambda")
print(as.numeric(commandArgs(trailingOnly=TRUE)[13]))

## stan options
print("niter")
print(strsplit(as.character(commandArgs(trailingOnly=TRUE)[14]),"_")[[1]])
print("nchains")
print(strsplit(as.character(commandArgs(trailingOnly=TRUE)[15]),"_")[[1]])
print("prop_warmup")
print(as.numeric(commandArgs(trailingOnly=TRUE)[16]))
print("max_treedepth")
print(as.numeric(commandArgs(trailingOnly=TRUE)[17]))
print("adapt_delta")
print(as.numeric(commandArgs(trailingOnly=TRUE)[18]))

## which simulation
print("simulation number")
print(as.numeric(commandArgs(trailingOnly=TRUE)[19]))
