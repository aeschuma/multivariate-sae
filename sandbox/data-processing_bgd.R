#  to do


rm(list = ls())
setwd('~/Desktop/survey-csmf/sandbox')

#### Libraries ####
library(SUMMER)
library(dplyr)
library(tidyr)
library(spdep)
library(geosphere)
library(readstata13)
library(knitr)
library(kableExtra)

help(package = "SUMMER", help_type = "html")
# utils::browseVignettes(package = "SUMMER")

#### Parameters ####

iso <- "bgd"
start_year <- 1990
end_year <- 2020
svy_year <- 2017

#### Load Data and format births####

mod.dat <- getBirths(filepath = "/Users/austin/Dropbox/dissertation_2/survey-csmf/data/bgd_dhs2017/data/BDBR7RDT/BDBR7RFL.DTA",
                     surveyyear = svy_year,
                     year.cut = seq(start_year, end_year, 1),
                     strata = c("v024", "v025"), compact = T)

#### Save frame ####

save(mod.dat, file = paste0('/Users/austin/Dropbox/dissertation_2/survey-csmf/data/bgd_dhs2017/data/', 'births-file_', iso, '.rda'))



## explore data

mod.dat %>% group_by(strata) %>%
  summarise(deaths = sum(died),
            person_months = sum(total)) %>%
  kable(format = 'latex')

mod.dat %>% group_by(time) %>%
  summarise(deaths = sum(died),
            person_months = sum(total)) %>%
  kable(format = 'latex', booktabs = TRUE, linesep = "") 

plotdat <- mod.dat %>% group_by(strata, time) %>%
  summarise(deaths = sum(died),
            person_months = sum(total)) %>% ungroup() 

ggplot(plotdat, aes(x = time, y = deaths, group = strata)) + geom_line(aes(color = strata))
