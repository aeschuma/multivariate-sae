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

## which models to run?
# see model-info.csv
model_to_run <- 6

## data generating options
number_of_causes <- 2
number_of_regions <- 50
number_of_replications <- 50

## parameters (TODO: currently must specify all of them, but not all models run will use all parameters)
## perhaps want to manually set the parameters in the run info CSV ahead of time and then just pull them from there?

beta1 <- 1
beta2 <- 2
rho_lower <- -0.2
rho_upper <- 0.2
sigmasq_lower <- 0.05
sigmasq_upper <- 0.5
sigma_gamma1 <- 1.5
sigma_gamma2 <- 2.5
lambda <- 0.5
sigma_delta <- 1
rho_gamma <- 0.5

## STAN options
niter <- 10000
nchains <- 2
prop_warmup <- 0.5
max_treedepth <- 25
adapt_delta <- 0.95

## Simulation number
number_of_sims <- 1000

## loop and submit jobs
for (a in model_to_run) {
    for (b in number_of_causes) {
        for (c in number_of_regions) {
            for (d in number_of_replications) {
                for (e in  beta1) {
                    for (f in  beta2) {
                        for (g in rho_lower) {
                            for (h in rho_upper) {
                                for (i in  sigmasq_lower) {
                                    for (j in  sigmasq_upper) {
                                        for (k in  sigma_gamma1) {
                                            for (l in  sigma_gamma2) {
                                                for (m in  lambda) {
                                                    for (mmm in sigma_delta) {
                                                        for (mm in  rho_gamma) {
                                                            for (n in  niter) {
                                                                for (o in  nchains) {
                                                                    for (p in  prop_warmup) {
                                                                        for (q in  max_treedepth) {
                                                                            for (r in  adapt_delta) {
                                                                                ## set run number from run info CSV
                                                                                run_info <- read.csv("results-run-info.csv")
                                                                                run_number <- max(run_info$run_number) + 1
                                                                                
                                                                                # save run info
                                                                                if (!testing_loop & !testing_script) {
                                                                                    new_info <- c(run_number, a,
                                                                                                  b, c, d,
                                                                                                  e, f, g, h, i, j,
                                                                                                  k, l, m, mmm, mm,
                                                                                                  n, o, p, q, r,
                                                                                                  number_of_sims)
                                                                                    run_info <- rbind(run_info, new_info)
                                                                                    write.csv(run_info, "results-run-info.csv", row.names = FALSE)
                                                                                }
                                                                                
                                                                                for (s in  1:number_of_sims) {
                                                                                    
                                                                                    if (testing_script) {
                                                                                        script <- "qsub-stage-2-sims-TEST.sh"
                                                                                    } else {
                                                                                        script <- "qsub-stage-2-sims.sh" 
                                                                                    }
                                                                                    sub <- paste0("qsub -l h=\"b34|b35|b36|b37\" -pe local ", o, " -v ",
                                                                                                  "a=",a,
                                                                                                  ",b=",b,
                                                                                                  ",c=",c,
                                                                                                  ",d=",d,
                                                                                                  ",e=",e,
                                                                                                  ",f=",f,
                                                                                                  ",g=",g,
                                                                                                  ",h=",h,
                                                                                                  ",ii=",i,
                                                                                                  ",j=",j,
                                                                                                  ",k=",k,
                                                                                                  ",l=",l,
                                                                                                  ",m=",m,
                                                                                                  ",mmm=",mmm,
                                                                                                  ",mm=",mm,
                                                                                                  ",n=",n,
                                                                                                  ",o=",o,
                                                                                                  ",p=",p,
                                                                                                  ",q=",q,
                                                                                                  ",r=",r,
                                                                                                  ",rr=",run_number,
                                                                                                  ",s=",s,
                                                                                                  " -N s2sim_",
                                                                                                  a,"_",
                                                                                                  b,"_",
                                                                                                  c,"_",
                                                                                                  d,"_",
                                                                                                  e, "_",
                                                                                                  f, "_",
                                                                                                  g,"_",
                                                                                                  h,"_",
                                                                                                  i,"_",
                                                                                                  j,"_",
                                                                                                  k, "_",
                                                                                                  l, "_",
                                                                                                  m,"_",
                                                                                                  mmm,"_",
                                                                                                  mm,"_",
                                                                                                  n,"_",
                                                                                                  o,"_",
                                                                                                  p,"_",
                                                                                                  q, "_",
                                                                                                  r, "_",
                                                                                                  run_number, "_",
                                                                                                  s,
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
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}
