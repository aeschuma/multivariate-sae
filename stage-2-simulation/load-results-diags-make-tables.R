rm(list = ls())

library(tidyverse)

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
library(knitr)

code_dir <- paste0(root,"Desktop/survey-csmf/stage-2-simulation")
dropbox_dir <- paste0(root,"Dropbox/dissertation_2/survey-csmf/results/stage-2-simulation")

# load results run info
setwd(code_dir)
res_run <- read_csv("results-run-info.csv")
model_info <- read_csv("model-info.csv")
dgm_info <- read_csv("dgm-info.csv")

# set dir to files dir
setwd(dropbox_dir)

# file list
fs <- grep("results", list.files(), value = TRUE)

## extract run number to name saved results
run_numbers <- c()
for (i in 1:length(fs)) {
    run_numbers <- c(run_numbers,
                     as.numeric(regmatches(fs[i], regexec("results-diags_run-\\s*(.*?)\\s*\\.Rdata", fs[i]))[[1]][2]))
}

run_numbers <- sort(run_numbers)

# read in results and display in tables
tmps <- vector(mode = "list", length = length(run_numbers))
for (i in 1:length(run_numbers)) {
    run_number <- run_numbers[i]
    myfile <- paste0("results-diags_run-",run_number,".Rdata")
    tmps[[i]] <- vector(mode = "list", length = 2)
    load(myfile)
    tmps[[i]][[1]] <- results_comp
    tmps[[i]][[2]] <- diags_comp
    rm(results_comp)
    rm(diags_comp)
    cat(paste0("model: ", res_run$model_number[res_run$run_number == run_number], 
               "; dgm: ", res_run$dgm_number[res_run$run_number == run_number], "\n"))
    print(kable(tmps[[i]][[1]], format = "markdown", caption = "Simulation results", digits = 3))
    print(kable(tmps[[i]][[2]], format = "markdown", caption = "Stan diagnostics"))
}
