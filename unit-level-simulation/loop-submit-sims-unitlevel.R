## Austin Schumacher
## 2021-07-05
## Loop and submit simulations to the cluster that will fit stage 2 models
## NOTE:  -l h="b34|b35|b36|b37" qsub option for just the newer nodes so STAN will run

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

## testing?
testing_loop <- FALSE
testing_script <- FALSE

## define directories

# working directory for code
wd <- paste0(root,"Desktop/survey-csmf/unit-level-simulation")

# directory to save results
savedir <- paste0(root,"Dropbox/dissertation_2/survey-csmf")

setwd(wd)

#######################
## SIMULATION SETTINGS
#######################

## which dgm to run?
# dgm_to_run <- 1:7
dgm_to_run <- c(7)

## Simulation number
number_of_sims <- 2

## loop and submit jobs
running_total <- 0
for (a in dgm_to_run) {
                        
    ## set run number from run info CSV
    run_info <- suppressWarnings(read.csv("results-run-info-unitlevel.csv"))
    run_number <- max(run_info$run_number) + 1
    
    # save run info
    if (!testing_loop & !testing_script) {
        new_info <- c(run_number, a, number_of_sims)
        run_info <- rbind(run_info, new_info)
        write.csv(run_info, "results-run-info-unitlevel.csv", row.names = FALSE)
    }
    
    for (s in  1:number_of_sims) {
        running_total <- running_total + 1
        if (testing_script) {
            script <- "qsub-sims-unitlevel-TEST.sh"
        } else {
            script <- "qsub-sims-unitlevel.sh" 
        }
        sub <- paste0("qsub -l h=\"biostat-b34|biostat-b35|biostat-b36|biostat-b37\" -pe local 2 -v ",
                      "a=",a,
                      ",rr=",run_number,
                      ",s=",s,
                      ",t=",running_total,
                      " -N s2sim_",running_total,
                      " ",
                      script)
        
        if (testing_loop) {
            print(sub)
        } else {
            system(sub)
        }
    }
}
