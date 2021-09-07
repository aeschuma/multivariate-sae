#  to do

rm(list = ls())

#### Libraries ####
library(SUMMER)
library(tidyverse)
library(tidyr)
library(spdep)
library(geosphere)
library(haven)
library(knitr)
library(kableExtra)
library(magrittr)

# help(package = "SUMMER", help_type = "html")
# utils::browseVignettes(package = "SUMMER")

#### Parameters ####

iso <- "bgd"
start_year <- 1990
end_year <- 2020
svy_year <- 2017

#### Load Data and format births ####



# rawbr <- read_stata("/Users/austin/Dropbox/dissertation_2/survey-csmf/data/bgd_dhs2017/data/BDBR7RDT/BDBR7RFL.DTA")
# rawbr %<>% mutate(mom.id = paste(v001, v002, v003, sep = "_"))

mod.dat <- getBirths(filepath = "/Users/austin/Dropbox/dissertation_2/survey-csmf/data/bgd_dhs2017/data/BDBR7RDT/BDBR7RFL.DTA",
                     variables = c("caseid", "v001", "v002",  "v003", "v004", "v005", "v008", "v021", "v022", "v023",
                                   "v024", "v025", "v139", "bidx", "b1", "b2"),
                     surveyyear = svy_year,
                     year.cut = seq(start_year, end_year, 1),
                     strata = c("v024", "v025"), compact = FALSE)

mod.dat %<>% mutate(mom.id = paste(v001, v002, v003, sep = "_"),
                    b1 = as.numeric(as.character(b1)),
                    b2 = as.numeric(as.character(b2))) %>% 
  as_tibble()

#### Merge on cause-identifiers ####

# load VA data
vadat <- read_stata("/Users/austin/Dropbox/dissertation_2/survey-csmf/data/bgd_dhs2017/data/BDVA7RDT/BDVA7RFL.DTA")

# formatting
vadat %<>% mutate(mom.id = paste(qncluster, qnhnumber, qnmother, sep = "_"),
                  birth.year = ifelse(is.na(qn303y), qc303y, qn303y),
                  birth.month = ifelse(is.na(qn303m), qc303m, qn303m)) %>%
  filter(!is.na(birth.year)) %>%
  select(mom.id, birth.year, birth.month, qfinicd)

## NOTE: need to get a better ID to merge on (currently dropping duplicates)

# merge
mod.dat.all <- mod.dat %>% full_join(vadat, by = c("mom.id" = "mom.id", 
                                                   "b1" = "birth.month", 
                                                   "b2" = "birth.year"))

summary(mod.dat.all$time[!is.na(mod.dat.all$qfinicd)])

# View(mod.dat.all %>% filter(is.na(v001)))

## 16 VA deaths are not in the births file:
mod.dat.all %>% filter(is.na(v001)) %>% pull(mom.id)

# only want those in 60 months prior to survey
mod.dat.all %<>% filter(obsmonth >= v008 - 60)

## Collapse age groups for multinomial model (so there are no zero death age-cause combos)
mod.dat.all$age_old <- mod.dat.all$age
mod.dat.all$age <- ifelse(mod.dat.all$age_old %in% c("12-23", "24-35", "36-47", "48-59"), "12-59", "0-11")
table(mod.dat.all$age)   

# cause groups
table(vadat$qfinicd)
attr(vadat$qfinicd, "labels")
table(mod.dat.all$age, mod.dat.all$qfinicd)

mod.dat.all$cause <- mod.dat.all$qfinicd
mod.dat.all$cause <- ifelse(mod.dat.all$died == 0,  "alive", mod.dat.all$cause)
table(mod.dat.all$cause, mod.dat.all$age)
table(mod.dat.all$cause, mod.dat.all$v024)
table(mod.dat.all$cause, mod.dat.all$age, mod.dat.all$v024)

mod.dat.all$cause.gp <- ifelse(mod.dat.all$qfinicd %in% c(6,7,13,3,14), "random causes of interest",
                               ifelse(!is.na(mod.dat.all$qfinicd), "other", "alive"))
mod.dat.all$cause.gp <- ifelse(mod.dat.all$died == 0,  "alive", mod.dat.all$cause.gp)

table(mod.dat.all$cause.gp, mod.dat.all$age, mod.dat.all$v024)


#### Save frame ####

save(mod.dat.all, file = paste0('/Users/austin/Dropbox/dissertation_2/survey-csmf/data/bgd_dhs2017/data/', 'births-file_', iso, '.rda'))

## explore data

# mod.dat %>% group_by(strata) %>%
#   summarise(deaths = sum(died),
#             person_months = sum(total)) %>%
#   kable(format = 'latex')
# 
# mod.dat %>% group_by(time) %>%
#   summarise(deaths = sum(died),
#             person_months = sum(total)) %>%
#   kable(format = 'latex', booktabs = TRUE, linesep = "") 
# 
# plotdat <- mod.dat %>% group_by(strata, time) %>%
#   summarise(deaths = sum(died),
#             person_months = sum(total)) %>% ungroup() 
# 
# ggplot(plotdat, aes(x = time, y = deaths, group = strata)) + geom_line(aes(color = strata))
