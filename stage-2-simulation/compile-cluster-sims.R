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
    params <- one_tmp[[2]]$param
    results <- vector(mode = "list", length = length(one_tmp))
    for (i in 1:length(results)) {
        results[[i]] <- vector(mode = "list", length = 8)
        names(results[[i]]) <- c("bias", "relative bias", "est", "squared error", 
                                 "coverage80", "width80", "coverage95", "width95")
        for (j in 1:length(results[[i]])) {
            results[[i]][[j]] <- matrix(NA, nrow = nrow(one_tmp[[2]]), ncol = length(tmp.resfiles))
        }
    }
    
    # loop through files and store results
    for (j in 1:length(tmp.resfiles)) {
        tmpfile <- readRDS(tmp.resfiles[j])
        for (i in 1:length(results)) {
            if (i == 1) {
                nadat <- data.frame(param = c("beta[1]", "beta[2]", 
                                              "sigma[1]", "sigma[2]", 
                                              "rho[1]", "rho[2]",
                                              "lambda"),
                                    truth = NA, est = NA, 
                                    coverage.80 = NA, width.80  = NA, 
                                    coverage.95 = NA, width.95 = NA)
                tmpfile[[i]] <- rbind(nadat,
                                      tmpfile[[i]])
            }
            results[[i]]$`bias`[, j] <- tmpfile[[i]]$est - tmpfile[[i]]$truth
            results[[i]]$`relative bias`[, j] <- (tmpfile[[i]]$est - tmpfile[[i]]$truth)/tmpfile[[i]]$truth
            results[[i]]$`est`[, j] <- tmpfile[[i]]$est
            results[[i]]$`squared error`[, j] <- (tmpfile[[i]]$est - tmpfile[[i]]$truth)^2
            results[[i]]$`coverage80`[, j] <- tmpfile[[i]]$coverage.80
            results[[i]]$`coverage95`[, j] <- tmpfile[[i]]$coverage.95
            results[[i]]$`width80`[, j] <- tmpfile[[i]]$width.80
            results[[i]]$`width95`[, j] <- tmpfile[[i]]$width.95
        }
    }
    
    results_comp <- list()
    length(results_comp) <- length(results)
    
    # compile measures across sims
    for (i in 1:length(results)) {
        results_comp[[i]] <- data.frame(param = params,
                                        bias = apply(results[[i]]$bias, 1, mean),
                                        relbias = apply(results[[i]]$`relative bias`, 1, mean),
                                        var = apply(results[[i]]$est, 1, var),
                                        mse = apply(results[[i]]$`squared error`, 1, mean),
                                        cov80 = apply(results[[i]]$coverage80, 1, mean),
                                        width80 = apply(results[[i]]$width80, 1, mean),
                                        cov95 = apply(results[[i]]$coverage95, 1, mean),
                                        width95 = apply(results[[i]]$width95, 1, mean))
        # collapse latent means
        lm_haz_params <- grep("haz.reg", params, value = TRUE)
        lm_waz_params <- grep("waz.reg", params, value = TRUE)
        lm_haz <- apply(results_comp[[i]] %>% 
                            filter(param %in% lm_haz_params) %>% 
                            select(-param),
                        2, mean)
        lm_waz <- apply(results_comp[[i]] %>% 
                            filter(param %in% lm_waz_params) %>% 
                            select(-param),
                        2, mean)
        
        lm_res <- as.data.frame(rbind(lm_haz, lm_waz))
        lm_res$param <- rownames(lm_res)
        rownames(lm_res) <- NULL
        results_comp[[i]] <- results_comp[[i]] %>% 
            filter(param %in% c("beta[1]", "beta[2]", 
                                "sigma[1]", "sigma[2]", 
                                "rho[1]", "rho[2]",
                                "lambda")) %>%
            bind_rows(lm_res)
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

