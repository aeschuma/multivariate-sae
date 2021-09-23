rm(list=ls())

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
    library(Rcpp); library(StanHeaders); library(BH); library(rstan);library(bayesplot);
}
library(scales); library(RColorBrewer); library(ggplot2); 

########
## TESTING THE CODE?
########

testing <- TRUE

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
run_numbers <- c()
for (i in 1:length(resfiles)) {
    run_numbers <- c(run_numbers,
                     as.numeric(regmatches(resfiles[i], regexec("results_run-\\s*(.*?)\\s*_", resfiles[i]))[[1]][2]))
}

run_numbers <- unique(run_numbers)
results_comp <- diags_comp <- vector(mode = "list", length = length(run_numbers))
for (rn in 1:length(run_numbers)) {
    run_number <- run_numbers[rn]
    
    # compile results
    results <- readRDS(resfiles[1])
    paramnames <- as.character(results$param)
    measurenames <- names(results)[which(names(results) != "param")]
    results <- array(NA, c(nrow(results), ncol(results) - 1, length(resfiles)))
    for (i in 1:length(resfiles)) {
        results[,,i] <- as.matrix(readRDS(resfiles[i])[, measurenames])
    }
    results_comp[[rn]] <- cbind(paramnames, as.data.frame(apply(results, c(1, 2), mean)))
    names(results_comp[[rn]]) <- c("param", paste0("mean_",measurenames))
    
    # compile diagnosics
    diags <- readRDS(diagfiles[1])
    stan_diag_names <- names(diags)
    diags <- matrix(NA, nrow = length(diagfiles), ncol = ncol(diags))
    for (i in 1:length(diagfiles)) {
        diags[i,] <- as.matrix(readRDS(diagfiles[i]))
    }
    diags_comp[[rn]] <- as.data.frame(t(apply(diags, 2, mean)))
    names(diags_comp[[rn]]) <- paste0("mean_", stan_diag_names)
    
    if (!testing) {
        # save results
        setwd(savedir)
        save(results_comp[[rn]], diags_comp[[rn]], file = paste0("results-diags_run-", run_number,".Rdata"))
        
        # move results to the real dropbox if on Box
        if (root == "P:/") {
            setwd("C:/Users/aeschuma/Dropbox/dissertation_2/survey-csmf/results/stage-2-simulation")
            save(results_comp[[rn]], diags_comp[[rn]], file = paste0("results-diags_run-", run_number,".Rdata"))
        } 
    }
}

if (!testing) {
    # delete tmp files
    setwd(paste0(savedir,"/tmp"))
    sapply(c(resfiles, diagfiles), unlink)
    setwd(paste0(savedir,"/../../out/stage-2-simulation"))
    sapply(list.files(), unlink)
    
    # delete .pe and .po files
    if (root == "P:/") {
        setwd(root)
        pe_po_files <- grep("s2sim_", list.files(), value = TRUE)
        sapply(pe_po_files, unlink)
    } else if (root == "/home/users/aeschuma/") {
        setwd(root)
        system("rm -f /home/users/aeschuma/s2sim_*")
    }
}

