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

# which runs do we want for the paper?
paper_dgms <- c(2, 4, 15)
paper_models <- c(1,2)
which_runs <- res_run %>% filter(dgm_number %in% paper_dgms & model_number %in% paper_models) %>% pull(run_number)
shared_mods <- res_run %>% filter(dgm_number %in% paper_dgms & model_number == 2) %>% pull(run_number)
nonshared_mods <- res_run %>% filter(dgm_number %in% paper_dgms & model_number == 1) %>% pull(run_number)

model_results <- vector(mode = "list", length = length(which_runs))

# read in results and display in tables
tmps <- vector(mode = "list", length = length(run_numbers))
for (i in 1:length(which_runs)) {
    run_number <- which_runs[i]
    myfile <- paste0("results-diags_run-",run_number,".Rdata")
    load(myfile)
    model_results[[i]] <- results_comp
    
    cat(paste0("model: ", res_run$model_number[res_run$run_number == run_number], 
               "; dgm: ", res_run$dgm_number[res_run$run_number == run_number], "\n"))
    print(kable(diags_comp, format = "markdown", caption = "Stan diagnostics"))
    
    rm(results_comp)
    rm(diags_comp)
}

model_results

# parameter name crosswalk
parnames <- tibble(param = model_results[which_runs %in% shared_mods][[1]][, c("param")],
                   parname = c("$\\beta_1$", "$\\beta_2$",
                               "$\\sigma_1$", "$\\sigma_2$",
                               "$\\rho_1$", "$\\rho_2$", 
                               "$\\lambda$", "latent means"))

# create big table
tables.list <- vector(mode = "list", length = length(paper_dgms))

for (j in 1:length(tables.list)) {
    s1runs <- res_run %>% filter(dgm_number == paper_dgms[j] & model_number %in% paper_models) %>% pull(run_number)
    tmp1 <- model_results[which_runs %in% s1runs & which_runs %in% nonshared_mods][[1]][, c("param",
                                                                                            "mean_absolute_bias",
                                                                                            "mean_coverage.95",
                                                                                            "mean_width.95")]
    tmp2 <- model_results[which_runs %in% s1runs & which_runs %in% shared_mods][[1]][, c("param",
                                                                                         "mean_absolute_bias",
                                                                                         "mean_coverage.95",
                                                                                         "mean_width.95")]
    tmp2$order <- 1:nrow(tmp2)
    tmp2 <- merge(tmp2, parnames)
    tmp2 <- tmp2[order(tmp2$order),]
    
    tmp1$order <- 1:nrow(tmp1)
    tmp1 <- merge(tmp1, parnames)
    tmp1 <- tmp1[order(tmp1$order),]
    
    tmptabmods <- rep(c("shared", "nonshared"), nrow(tmp2))
    sctab <- c()
    for (i in 1:nrow(tmp2)) {
        parname <- tmp2$param[i]
        
        if (parname == "lambda") {
            tt <- rbind(tmp2[tmp2$param == parname,c("parname", 
                                                     "mean_absolute_bias",
                                                     "mean_coverage.95",
                                                     "mean_width.95")], 
                        c("", rep("---", ncol(tmp1) - 1)))
        } else {
            tmp1$parname[tmp1$param == parname] <- ""
            tt <- rbind(tmp2[tmp2$param == parname,c("parname", 
                                                     "mean_absolute_bias",
                                                     "mean_coverage.95",
                                                     "mean_width.95")], 
                        tmp1[tmp1$param == parname,c("parname", 
                                                     "mean_absolute_bias",
                                                     "mean_coverage.95",
                                                     "mean_width.95")])
        }
        tmptab <- cbind(model = tmptabmods[((i-1)*2 + 1):((i-1)*2 + 2)],
                        tt)
        sctab <- rbind(sctab, tmptab)
        sctab <- sctab[,c("parname", 
                          "model", 
                          "mean_absolute_bias",
                          "mean_coverage.95",
                          "mean_width.95")]
        for (vv in c("mean_absolute_bias", "mean_coverage.95", "mean_width.95")) {
            sctab[,vv] <- as.numeric(sctab[, vv])
        }
    }
    rownames(sctab) <- NULL
    tables.list[[j]] <- sctab
}

kable(tables.list, booktabs = TRUE, format = "markdown", caption = "Simulation results", valign = 't',
      col.names = c("Model", "Parameter", "Mean bias", "Mean 95\\% coverage", "Mean 95\\% width"),
      digits = 4)

saveRDS(tables.list, file = "../proj-2-chapter-results/simulation-results-table.rds")
