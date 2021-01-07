# Austin Schumacher
# 2/28/2020
# load sample USA data

# libraries
library(tidyverse); library(haven);

# directories
wonder_dir <- "../../../Dropbox/Age specific cod_shared/US-WONDER"

# set directory of Wonder data
setwd(wonder_dir)

# load one year of data
wonder <- list()
years <- c(2015)
i <- 1
year <- years[i]
wonder[[i]] <- read_stata(paste0("mort",year,".dta"))
