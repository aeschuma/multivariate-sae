rm(list=ls())

## set the root depending on operating system
root <- ifelse(Sys.info()[1]=="Darwin","~/",
               ifelse(Sys.info()[1]=="Windows","P:/",
                      ifelse(Sys.info()[1]=="Linux","/home/users/aeschuma/",
                             stop("Unknown operating system"))))

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

    # compile results
    results <- readRDS(tmp.resfiles[1])
    paramnames <- vector(mode = "list", length = length(results))
    measurenames <- vector(mode = "list", length = length(results))
    results_comp <- vector(mode = "list", length = length(results))
    for (j in 1:length(results)) {
        paramnames[[j]] <- as.character(results[[j]]$param)
        measurenames[[j]] <- names(results[[j]])[which(names(results[[j]]) != "param")]
        results[[j]] <- array(NA, c(nrow(results[[j]]), ncol(results[[j]]) - 1, length(tmp.resfiles)))
    }    
    for (i in 1:length(tmp.resfiles)) {
        tmpres <- readRDS(tmp.resfiles[i])
        for (j in 1:length(tmpres)) {
            results[[j]][,,i] <- as.matrix(tmpres[[j]][, measurenames[[j]]])
        }
    }
    for (j in 1:length(results)) {
        results_comp[[j]] <- cbind(paramnames[[j]], as.data.frame(apply(results[[j]], c(1, 2), mean)))
        names(results_comp[[j]]) <- c("param", paste0("mean_",measurenames[[j]]))
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

