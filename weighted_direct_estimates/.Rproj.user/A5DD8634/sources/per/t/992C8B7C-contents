# Austin
# 12/17/2019
# Calculate weighted direct estimates of child mortality from MDS data
# - all cause U5MR and separately NMR
#   - national yearly
#   - each state yearly
# - some common causes yearly
#   - birth trauma NMR 
#     - national over time
#     - each state over time

# load libraries
library(SUMMER); library(INLA); library(survey);
library(tidyverse); 
library(sp); library(maptools); library(rgdal);
library(raster); library(milliondeaths);  library(RPostgreSQL);


# which data are we using?
dat <- allcause_dat
# dat <- allData

# look at the data
dim(dat$data)
head(dat$data)
summary(as.data.frame(dat$data))
summary(as.data.frame(dat$data[is.na(dat$data$pop),]))

# collapse data into U5, both sexes (for now, we're only doing one cause)


# IHME data
ihme <- read.csv("/Users/austin/Dropbox/csmr_india/data/ihme/gbd2017_india_national_cause/IHME-GBD_2017_DATA-cd1311cc-1.csv",
                 header = TRUE, stringsAsFactors = FALSE)

#------------------------
# birth trauma infant national over time
#------------------------

# my.svydesign <- survey::svydesign(ids = ~srs_id, 
#                                   strata = ~strata, nest = T, weights = ~weight_alan, data = tmp)

# glm.ob <- survey::svyglm(Birth.asphyxia.birth.trauma ~ 1, 
#                          design = my.svydesign, family = quasipoisson, 
#                          maxit = 50, offset = log(my.svydesign$variables$pop))
# 


