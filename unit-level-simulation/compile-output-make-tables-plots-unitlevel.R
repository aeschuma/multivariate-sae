# Austin
# create results from all simulations (stage 2) to make a table in the chapter
## 3/20/2022

rm(list = ls())

library(tidyverse)
library(spdep)
library(geosphere)
library(haven)
library(knitr)
library(kableExtra)
library(magrittr)
library(rgdal)
library(viridis)
library(classInt)
library(gridExtra)
library(ggpubr)

# load results run info
res_run <- read_csv("results-run-info-unitlevel.csv")
dgm_info <- read_csv("dgm-info-unitlevel.csv")
model_names <- c("I",
                 "II",
                 "III", 
                 "IV")

dropbox_dir <- "~/Dropbox/dissertation_2/survey-csmf/results/unit-level-simulation"

# file list
fs <- grep("results", list.files(dropbox_dir), value = TRUE)

## extract run number to name saved results
run_numbers <- c()
for (i in 1:length(fs)) {
    tmprn <- regmatches(fs[i], regexec("results_run-\\s*(.*?)\\s*\\.Rdata", fs[i]))[[1]][2]
    if (is.na(as.numeric(tmprn))) tmprn <- str_split(tmprn, "_")[[1]][1]
    run_numbers <- c(run_numbers,
                     as.numeric(tmprn))
}

run_numbers <- sort(run_numbers)

# read in results and display in tables
tmps <- vector(mode = "list", length = length(run_numbers))
which_runs <- c(1:7) # remember to put in the order of the scenarios in the sims in the chapter

run_results <- vector(mode = "list", length = length(which_runs))
for (i in 1:length(which_runs)) {
    # cat(paste(i, "\n"))
    run_number <- run_numbers[which_runs[i]]
    myfile <- grep(paste0("results_run-",run_number,"\\s*(.*?)\\s*"), fs, value = TRUE)
    load(paste0(dropbox_dir, "/", myfile))
    tmps[[i]] <- results_comp
    rm(results_comp)
    
    dgmnum <- res_run$dgm_number[res_run$run_number == run_number]
    
    run_results[[i]] <- cbind(model = c(model_names[1], ""), 
                              tmps[[i]][[1]][tmps[[i]][[1]]$param %in% c("HAZ latent means", "WAZ latent means"),])
    if (length(tmps[[i]]) > 1) {
        for (j in 2:length(tmps[[i]])) {
            # extract latent mean results
            tmp_res <- cbind(model = c(model_names[j], ""), 
                             tmps[[i]][[j]][tmps[[i]][[j]]$param %in% c("HAZ latent means", "WAZ latent means"),])
            run_results[[i]] <- rbind(run_results[[i]], tmp_res)
        }   
    }
    run_results[[i]]$param <- ifelse(run_results[[i]]$param == "HAZ latent means", "outcome 1", "outcome 2")
    run_results[[i]]$dgm <- dgmnum
    rownames(run_results[[i]]) <- NULL
}

write_rds(run_results, paste0(dropbox_dir, "/", "simulation-results-tables.rds"))
write_rds(run_results, paste0(dropbox_dir, "/../proj-3-chapter-results/simulation-results-tables.rds"))

# load tables
for (i in 1:length(run_results)) {
    run_results[[i]]$scenario <- run_results[[i]]$dgm
    run_results[[i]]$Model <- rep(c("I", "II", "III", "IV"), each = 2)
    run_results[[i]]$Model <- factor(run_results[[i]]$Model, levels = c("I", "II", "III", "IV"))
}
all_sim_tables <- do.call(rbind, run_results)

plotResults <- function(data, scenarios, measure, measure_name) {
    scenarios_xaxis <- 1:length(scenarios)
    scenario_xwalk <- tibble(scenario = scenarios, 
                             scenario_xaxis = scenarios_xaxis)
    mydata <- data %>% filter(scenario %in% scenarios)
    mydata <- mydata %>% left_join(scenario_xwalk)
    key_pch <- as.character(unique(mydata$Model))
    mydata$scenario_xaxis <- mydata$scenario_xaxis + (as.numeric(mydata$Model) + 0.4)/(max(as.numeric(mydata$Model)) + 2)
    ggplot(mydata, 
           aes(x = scenario_xaxis, y = get(measure), label = Model, color = Model)) +
        geom_text(size = 3.5, fontface = 'bold') +
        geom_vline(xintercept = 1:(length(scenarios) + 1), alpha = 0.65) +
        scale_fill_viridis(discrete = TRUE, option = "C", alpha = 0.75) +
        scale_x_continuous(breaks = (1:length(scenarios)) + 0.5, labels = scenarios) +
        facet_wrap(~ param, scales = "fixed") +
        xlab("Scenario") +
        ylab(measure_name) +
        ggtitle(measure_name) +
        theme_light() +
        theme(panel.grid.major.x = element_blank() ,
            panel.grid.major.y = element_line(size=.1, color="gray") 
        )
}

# graph scenarios 1-4
p1 <- plotResults(data = all_sim_tables, scenarios = 1:4, measure = "mean_bias", measure_name = "Bias")
p2 <- plotResults(data = all_sim_tables, scenarios = 1:4, measure = "mean_absolute_bias", measure_name = "Absolute bias")
p3 <- plotResults(data = all_sim_tables, scenarios = 1:4, measure = "mean_variance", measure_name = "Variance")
p4 <- plotResults(data = all_sim_tables, scenarios = 1:4, measure = "mean_mse", measure_name = "MSE")
p5 <- plotResults(data = all_sim_tables, scenarios = 1:4, measure = "mean_coverage.95", measure_name = "95% coverage")
p6 <- plotResults(data = all_sim_tables, scenarios = 1:4, measure = "mean_width.95", measure_name = "95% width")

ggarrange(p1, p2, p3, p4, p5, p6, 
          ncol = 2, nrow = 3, legend = "none")
ggsave(paste0(dropbox_dir, "/../proj-3-chapter-results/simulation-results-graph-scenarios1thru4.pdf"), 
       width = 8, height = 10)

# graph scenarios 5-7
p1b <- plotResults(data = all_sim_tables, scenarios = 5:7, measure = "mean_bias", measure_name = "Bias")
p2b <- plotResults(data = all_sim_tables, scenarios = 5:7, measure = "mean_absolute_bias", measure_name = "Absolute bias")
p3b <- plotResults(data = all_sim_tables, scenarios = 5:7, measure = "mean_variance", measure_name = "Variance")
p4b <- plotResults(data = all_sim_tables, scenarios = 5:7, measure = "mean_mse", measure_name = "MSE")
p5b <- plotResults(data = all_sim_tables, scenarios = 5:7, measure = "mean_coverage.95", measure_name = "95% coverage")
p6b <- plotResults(data = all_sim_tables, scenarios = 5:7, measure = "mean_width.95", measure_name = "95% width")

ggarrange(p1b, p2b, p3b, p4b, p5b, p6b, 
          ncol = 2, nrow = 3, legend = "none")
ggsave(paste0(dropbox_dir, "/../proj-3-chapter-results/simulation-results-graph-scenarios5thru7.pdf"), 
       width = 8, height = 10)
