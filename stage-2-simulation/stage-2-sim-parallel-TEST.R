rm(list=ls())

# Preamble ####

## set the root depending on operating system
root <- ifelse(Sys.info()[1]=="Darwin","~/",
               ifelse(Sys.info()[1]=="Windows","P:/",
                      ifelse(Sys.info()[1]=="Linux","/home/users/aeschuma/",
                             stop("Unknown operating system"))))

## the following code makes rstan work on the Box server and the cluster
if (root == "P:/") {
    Sys.setenv(HOME="C:/Users/aeschuma",
               R_USER="C:/Users/aeschuma",
               R_LIBS_USER="C:/Users/aeschuma/R_libraries")
    .libPaths("C:/Users/aeschuma/R_libraries")
} else if (root == "/home/users/aeschuma/") {
    Sys.setenv(HOME=root,
               R_USER=root,
               R_LIBS_USER=paste0(root,"R/x86_64-pc-linux-gnu-library/3.6"))
    .libPaths(paste0(root,"R/x86_64-pc-linux-gnu-library/3.6"))
}

## load libraries
if (root == "/home/users/aeschuma/") {
    library(Rcpp,lib.loc=paste0(root,"R/x86_64-pc-linux-gnu-library/3.6"));
    library(StanHeaders,lib.loc=paste0(root,"R/x86_64-pc-linux-gnu-library/3.6")); 
    library(BH,lib.loc = paste0(root,"R/x86_64-pc-linux-gnu-library/3.6"))
    library(rstan,lib.loc = paste0(root,"R/x86_64-pc-linux-gnu-library/3.6"));
} else {
    library(Rcpp); library(StanHeaders); library(BH); library(rstan); library(bayesplot); 
}
library(mvtnorm); library(MASS);library(gtools); library(parallel);
library(scales); library(RColorBrewer); library(ggplot2); library(tidyr);
library(haven); library(knitr); library(rgdal); library(INLA);

if (root == "~/") library(cmdstanr);

## TESTING THE CODE?
testing <- TRUE

## define directories

# working directory for code
wd <- paste0(root,"Desktop/survey-csmf/stage-2-simulation")

# directory to save results
savedir <- paste0(root,"Dropbox/dissertation_2/survey-csmf/results/stage-2-simulation")

## set directory
setwd(wd)

# set number of cores to the max possible
options(mc.cores = detectCores())

# Load data and functions ####

# load simulation functions:
#   simulateData()
#   fitSTAN()
source("simulation-functions.R")

# load model csv to load which model we're running
models_dat <- read_csv("model-info.csv")
# only if need to rewrite csv to get rid of warning message for incomplete final line
# write.csv(models_dat, file = "model-info.csv", row.names = FALSE)

## set directory
setwd(wd)
setwd(savedir)

# set number of cores to the max possible
options(mc.cores = parallel::detectCores())

## which model to run
print("model_to_run")
print(as.numeric(commandArgs(trailingOnly=TRUE)[1]))

## data generation options
print("dgm_to_run")
print(as.numeric(commandArgs(trailingOnly=TRUE)[2]))

## stan options
print("niter")
print(as.numeric(commandArgs(trailingOnly=TRUE)[3]))
print("nchains")
print(as.numeric(commandArgs(trailingOnly=TRUE)[4]))
print("prop_warmup")
print(as.numeric(commandArgs(trailingOnly=TRUE)[5]))
print("max_treedepth")
print(as.numeric(commandArgs(trailingOnly=TRUE)[6]))
print("adapt_delta")
print(as.numeric(commandArgs(trailingOnly=TRUE)[7]))

## run number
print("run number")
print(as.numeric(commandArgs(trailingOnly=TRUE)[8]))

## which simulation
print("simulation number")
print(as.numeric(commandArgs(trailingOnly=TRUE)[9]))
