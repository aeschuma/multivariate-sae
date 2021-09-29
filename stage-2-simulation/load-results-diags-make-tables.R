rm(list = ls())

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

setwd(paste0(root,"Dropbox/dissertation_2/survey-csmf/results/stage-2-simulation"))

# read in results, display in tables
fs <- grep("results", list.files(), value = TRUE)

tmps <- vector(mode = "list", length = length(fs))

for (i in 1:length(tmps)) {
    tmps[[i]] <- vector(mode = "list", length = 2)
    load(fs[i])
    tmps[[1]] <- results_comp
    tmps[[2]] <- diags_comp
    rm(results_comp)
    rm(diags_comp)
    print(paste0("results run ", i))
    print(kable(tmps[[1]], format = "markdown"))
    print(kable(tmps[[2]], format = "markdown"))
}
