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
wd <- paste(root,"Desktop/dissertation/model_development/recover_parameters",sep="")

# directory to save results
savedir <- paste(root,"Dropbox/dissertation_2/cause_specific_child_mort/model_development/recover_parameters",sep="")

## set directory
setwd(wd)
setwd(savedir)

# set number of cores to the max possible
options(mc.cores = parallel::detectCores())

## data generation options
print("cause_re")
print(as.logical(commandArgs(trailingOnly=TRUE)[1]))
print("few_corrs")
print(as.logical(commandArgs(trailingOnly=TRUE)[2]))
print("time_rw")
print(as.logical(commandArgs(trailingOnly=TRUE)[3]))
print("no_corr")
print(as.logical(commandArgs(trailingOnly=TRUE)[4]))
print("FE_interactions")
print(as.logical(commandArgs(trailingOnly=TRUE)[5]))

## modeling options
print("cause_re_model")
print(as.logical(commandArgs(trailingOnly=TRUE)[6]))
print("time_rw_model")
print(as.logical(commandArgs(trailingOnly=TRUE)[7]))
print("no_corr_model")
print(as.logical(commandArgs(trailingOnly=TRUE)[8])) # IMPLIES CAUSE RE AND TIME RW BUT NO CORRELATIONS
print("FE_interactions_model")
print(as.logical(commandArgs(trailingOnly=TRUE)[9]))

## dimensions
print("nyear")
print(as.numeric(commandArgs(trailingOnly=TRUE)[10]))
print("nreg")
print(as.numeric(commandArgs(trailingOnly=TRUE)[11]))
print("nageg")
print(as.numeric(commandArgs(trailingOnly=TRUE)[12]))
print("ncause")
print(as.numeric(commandArgs(trailingOnly=TRUE)[13]))
print("rw_vars")
print(strsplit(as.character(commandArgs(trailingOnly=TRUE)[14]),"_")[[1]])
print("re_vars")
print(strsplit(as.character(commandArgs(trailingOnly=TRUE)[15]),"_")[[1]])

## parameters/hyperpars
print("sigma_rw")
print(as.numeric(commandArgs(trailingOnly=TRUE)[16]))
print("sigma_re")
print(as.numeric(commandArgs(trailingOnly=TRUE)[17]))
print("rho")
print(as.numeric(commandArgs(trailingOnly=TRUE)[18]))
print("exposure")
print(as.numeric(commandArgs(trailingOnly=TRUE)[19]))
print("lkj_hyperpar")
print(as.numeric(commandArgs(trailingOnly=TRUE)[20]))

## STAN options
print("niter")
print(as.numeric(commandArgs(trailingOnly=TRUE)[21]))
print("nthin")
print(as.numeric(commandArgs(trailingOnly=TRUE)[22]))
print("prop_warmup")
print(as.numeric(commandArgs(trailingOnly=TRUE)[23]))

# what chain is this?
print("chain_number")
print(as.numeric(commandArgs(trailingOnly=TRUE)[24]))
