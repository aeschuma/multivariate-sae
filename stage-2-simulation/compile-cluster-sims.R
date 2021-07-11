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

########
## TESTING THE CODE?
########

testing <- FALSE
