rm(list=ls())

## set the root depending on operating system
root <- ifelse(Sys.info()[1]=="Darwin","~/",
               ifelse(Sys.info()[1]=="Windows","P:/",
                      ifelse(Sys.info()[1]=="Linux","/home/users/aeschuma/",
                             stop("Unknown operating system"))))

library(tidyverse); library(magrittr);

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

## extract run number to name saved results
run_numbers <- c()
for (i in 1:length(resfiles)) {
    run_numbers <- c(run_numbers,
                     as.numeric(regmatches(resfiles[i], regexec("results_run-\\s*(.*?)\\s*_", resfiles[i]))[[1]][2]))
}

run_numbers <- sort(unique(run_numbers))

for (rn in 1:length(run_numbers)) {
    setwd(paste0(savedir,"/tmp"))
    
    cat(paste0("Compiling run: ", run_numbers[rn], "\n")); flush.console()
    
    run_number <- run_numbers[rn]
    
    tmp.resfiles <- grep(paste0("results_run-", run_number, "_sim-"), resfiles, value = TRUE)

    # storage for results
    one_tmp <- readRDS(tmp.resfiles[1])
    results <- vector(mode = "list", length = length(one_tmp))
    for (i in 1:length(results)) {
        results[[i]] <- vector(mode = "list", length = 8)
        names(results[[i]]) <- c("bias", "relative bias", "variance", "mse", 
                                 "coverage80", "width80", "coverage95", "width95")
        for (j in 1:length(results[[i]])) {
            results[[i]][[j]] <- matrix(NA, nrow = nrow(one_tmp[[1]]), ncol = length(tmp.resfiles))
        }
    }
    
    # loop through files and store results
    for (j in 1:length(tmp.resfiles)) {
        tmpfile <- readRDS(tmp.resfiles[j])
        for (i in 1:length(results)) {
            results[[i]]$`bias`[, j] <- tmpfile[[i]]$est - tmpfile[[i]]$truth
            results[[i]]$`relative bias`[, j] <- (tmpfile[[i]]$est - tmpfile[[i]]$truth)/tmpfile[[i]]$truth
            results[[i]]$`variance`[, j] <- var(tmpfile[[i]]$est)
            results[[i]]$`mse`[, j] <- mean((tmpfile[[i]]$est - tmpfile[[i]]$truth)^2)
            results[[i]]$`coverage80`[, j] <- tmpfile[[i]]$coverage.80
            results[[i]]$`coverage95`[, j] <- tmpfile[[i]]$coverage.95
            
        }
    }
   
   
    if (!testing) {
        
        # save results
        setwd(savedir)
        save(results_comp, file = paste0("results_run-", run_number,".Rdata"))
        
        # move results to the real dropbox if on Box
        if (root == "P:/") {
            setwd("C:/Users/aeschuma/Dropbox/dissertation_2/survey-csmf/results/stage-2-simulation")
            save(results_comp, file = paste0("results_run-", run_number,".Rdata"))
        } 
    }
}

if (!testing) {

    if (root == "P:/") {
        # delete tmp files
        setwd(paste0(savedir,"/tmp"))
        sapply(c(resfiles), unlink)
        setwd(paste0(savedir,"/../../out/stage-2-simulation"))
        sapply(list.files(), unlink)
        
        # delete .pe and .po files
        setwd(root)
        pe_po_files <- grep("s2sim_", list.files(), value = TRUE)
        sapply(pe_po_files, unlink)
    } else if (root == "/home/users/aeschuma/") {
        
        # delete tmp files
        tmp_del_res <- paste0(savedir,"/tmp/results*")
        system(paste0("rm -f ", tmp_del_res))

        # delete .pe and .po files
        system("rm -f /home/users/aeschuma/s2sim_*")
    }
}

