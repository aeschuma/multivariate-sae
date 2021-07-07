## Austin Schumacher
## 2021-07-05
## Loop and submit simulations to the cluster that will fit stage 2 models
## NOTE:  -l h="b34|b35|b36|b37" qsub option for just the newer nodes so STAN will run

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

## testing?
testing_loop <- FALSE
testing_script <- TRUE

## define directories

# working directory for code
wd <- paste0(root,"Desktop/survey-csmf/stage-2-simulation")

# directory to save results
savedir <- paste0(root,"Dropbox/dissertation_2/survey-csmf")

#######################
## SIMULATION SETTINGS
#######################

## which models to run?
# see model-info.csv
models_to_run <- c(1)

## data generating options
number_of_causes <- 2
number_of_regions <- 8
number_of_replications <- 1

## parameters
# beta1 <- 1
# beta2 <- 2
# rho_lower <- -0.2
# rho_upper <- 0.2
# sigmasq_lower <- 0.05
# sigmasq_upper <- 0.5
# sigma_gamma1 <- 1.5
# sigma_gamma2 <- 2.5
# lambda <- 0.5

## STAN options
# niter <- 10000
# nchains <- 3
# prop_warmup <- 0.5
# max_treedepth <- 15
# adapt_delta <- 0.8

## Simulation number
number_of_sims <- 5

## loop and submit jobs
setwd(wd)
for (a in models_to_run) {
    for (b in number_of_causes) {
        for (c in number_of_regions) {
            for (d in number_of_replications) {
                for (e in  1:number_of_sims) {
                    if (testing_script) {
                        script <- "qsub-stage-2-sims-TEST.sh"
                    } else {
                        script <- "qsub-stage-2-sims.sh" 
                    }
                    sub <- paste0("qsub -l h=\"b34|b35|b36|b37\" -pe local 4 -v ",
                                  "a=",a,
                                  ",b=",b,
                                  ",c=",c,
                                  ",d=",d,
                                  ",e=",e,
                                  " -N s2sim_",
                                  a,"_",
                                  b,"_",
                                  c,"_",
                                  d,"_",
                                  e,
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
