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
wd <- paste0(root,"Desktop/survey-csmf/stage-2-simulation")

# directory to save results
savedir <- paste0(root,"Dropbox/dissertation_2/survey-csmf")

setwd(wd)

#######################
## SIMULATION SETTINGS
#######################

## which model and dgm to run?
# see model-info.csv
model_to_run <- c(1, 2)
# model_to_run <- 3:4
# dgm_to_run <- c(1, 2)
dgm_to_run <- c(3, 8, 15, 16)

## STAN options
niter <- 10000
nchains <- 2
prop_warmup <- 0.5
max_treedepth <- 25
adapt_delta <- 0.8

## Simulation number
number_of_sims <- 500

## loop and submit jobs
running_total <- 0
for (a in model_to_run) {
    for (b in dgm_to_run) {
        for (c in  niter) {
            for (d in  nchains) {
                for (e in  prop_warmup) {
                    for (f in  max_treedepth) {
                        for (g in  adapt_delta) {
                            
                            ## set run number from run info CSV
                            run_info <- suppressWarnings(read.csv("results-run-info.csv"))
                            run_number <- max(run_info$run_number) + 1
                            
                            # save run info
                            if (!testing_loop & !testing_script) {
                                new_info <- c(run_number, a,
                                              b, c, d,
                                              e, f, g,
                                              number_of_sims)
                                run_info <- rbind(run_info, new_info)
                                write.csv(run_info, "results-run-info.csv", row.names = FALSE)
                            }
                            
                            for (s in  1:number_of_sims) {
                                running_total <- running_total + 1
                                if (testing_script) {
                                    script <- "qsub-stage-2-sims-TEST.sh"
                                } else {
                                    script <- "qsub-stage-2-sims.sh" 
                                }
                                sub <- paste0("qsub -l h=\"biostat-b34|biostat-b35|biostat-b36|biostat-b37\" -pe local ", d, " -v ",
                                              "a=",a,
                                              ",b=",b,
                                              ",c=",c,
                                              ",d=",d,
                                              ",e=",e,
                                              ",f=",f,
                                              ",g=",g,
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
                    }
                }
            }
        }
    }
}

