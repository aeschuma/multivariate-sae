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
res_run <- read_csv("results-run-info.csv")
dgm_info <- read_csv("dgm-info.csv")
model_names <- c("O", 
                 "I",
                 "II",
                 "III", 
                 "IV",
                 "V", 
                 "VI")

dropbox_dir <- "~/Dropbox/dissertation_2/survey-csmf/results/stage-2-simulation"

# file list
fs <- grep("results", list.files(dropbox_dir), value = TRUE)

## extract run number to name saved results
run_numbers <- c()
for (i in 1:length(fs)) {
    run_numbers <- c(run_numbers,
                     as.numeric(regmatches(fs[i], regexec("results_run-\\s*(.*?)\\s*\\.Rdata", fs[i]))[[1]][2]))
}

run_numbers <- sort(run_numbers)

# read in results and display in tables
tmps <- vector(mode = "list", length = length(run_numbers))
which_runs <- c(1, 6:10, 2:3, 5) # remember to put in the order of the scenarios in the sims in the chapter

run_results <- vector(mode = "list", length = length(which_runs))
for (i in 1:length(which_runs)) {
    # cat(paste(i, "\n"))
    run_number <- run_numbers[which_runs[i]]
    myfile <- paste0("results_run-",run_number,".Rdata")
    load(paste0(dropbox_dir, "/", myfile))
    tmps[[i]] <- results_comp
    rm(results_comp)
    
    dgmnum <- res_run$dgm_number[res_run$run_number == run_number]
    
    run_results[[i]] <- cbind(model = c(model_names[1], ""), tmps[[i]][[1]])
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
write_rds(run_results, paste0(dropbox_dir, "/../proj-2-chapter-results/simulation-results-tables.rds"))

# load tables
for (i in 1:length(run_results)) {
    run_results[[i]]$scenario <- i
    run_results[[i]]$Model <- rep(c("O", "I", "II", "III", "IV", "V", "VI"), each = 2)
    run_results[[i]]$Model <- factor(run_results[[i]]$Model, levels = c("O", "I", "II", "III", "IV", "V", "VI"))
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

# graph scenarios 1-6
p1 <- plotResults(data = all_sim_tables, scenarios = 1:6, measure = "mean_bias", measure_name = "Mean bias")
p2 <- plotResults(data = all_sim_tables, scenarios = 1:6, measure = "mean_absolute_bias", measure_name = "Mean absolute bias")
p3 <- plotResults(data = all_sim_tables, scenarios = 1:6, measure = "mean_coverage.95", measure_name = "95% coverage")
p4 <- plotResults(data = all_sim_tables, scenarios = 1:6, measure = "mean_width.95", measure_name = "95% width")

ggarrange(p1, p2, p3, p4, ncol = 1, nrow = 4, legend = "none")
ggsave(paste0(dropbox_dir, "/../proj-2-chapter-results/simulation-results-graph-scenarios1thru6.pdf"), 
       width = 8.5, height = 11)

# graph scenarios 7-9
p5 <- plotResults(data = all_sim_tables, scenarios = 7:9, measure = "mean_bias", measure_name = "Mean bias")
p6 <- plotResults(data = all_sim_tables, scenarios = 7:9, measure = "mean_absolute_bias", measure_name = "Mean absolute bias")
p7 <- plotResults(data = all_sim_tables, scenarios = 7:9, measure = "mean_coverage.95", measure_name = "95% coverage")
p8 <- plotResults(data = all_sim_tables, scenarios = 7:9, measure = "mean_width.95", measure_name = "95% width")

ggarrange(p5, p6, p7, p8, ncol = 1, nrow = 4, legend = "none")
ggsave(paste0(dropbox_dir, "/../proj-2-chapter-results/simulation-results-graph-scenarios7thru9.pdf"), 
       width = 8.5, height = 11)
