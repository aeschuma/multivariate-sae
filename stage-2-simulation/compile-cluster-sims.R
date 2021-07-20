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
library(parallel);library(scales); library(RColorBrewer);library(data.table);library(ggplot2); 

########
## TESTING THE CODE?
########

testing <- FALSE

## define directories

# working directory for code
wd <- paste0(root,"Desktop/survey-csmf/stage-2-simulation")

# directory to save results
savedir <- paste0(root,"Dropbox/dissertation_2/survey-csmf/results/stage-2-simulation")

## set directory
setwd(paste0(savedir,"/tmp"))

## files
files <- list.files()
resfiles <- grep("results_", files, value = TRUE)
diagfiles <- grep("standiags_", files, value = TRUE)

## extract run number to name saved results
run_number <- as.numeric(regmatches(resfiles[1], regexec("results_run-\\s*(.*?)\\s*_", resfiles[1]))[[1]][2])

# compile results
results <- readRDS(resfiles[1])
paramnames <- as.character(results$param)
measurenames <- names(results)[which(names(results) != "param")]
results <- array(NA, c(nrow(results), ncol(results) - 1, length(resfiles)))
for (i in 1:length(resfiles)) {
    results[,,i] <- as.matrix(readRDS(resfiles[i])[, measurenames])
}
results_comp <- cbind(paramnames, as.data.frame(apply(results, c(1, 2), mean)))
names(results_comp) <- c("param", paste0("mean_",measurenames))

# compile diagnosics
diags <- readRDS(diagfiles[1])
stan_diag_names <- names(diags)
diags <- matrix(NA, nrow = length(diagfiles), ncol = ncol(diags))
for (i in 1:length(diagfiles)) {
    diags[i,] <- as.matrix(readRDS(diagfiles[i]))
}
diags_comp <- as.data.frame(t(apply(diags, 2, mean)))
names(diags_comp) <- paste0("mean_", stan_diag_names)

if (!testing) {
    # save results
    setwd(savedir)
    save(results_comp, diags_comp, file = paste0("results-diags_run-", run_number,".Rdata"))
    
    # delete tmp files
    setwd(paste0(savedir,"/tmp"))
    sapply(c(resfiles, diagfiles), unlink)
    setwd(paste0(savedir,"/../../out/stage-2-simulation"))
    sapply(list.files(), unlink)
    
    setwd(root)
    system("rm -f /home/users/aeschuma/qsub-compile-sims*")
}
