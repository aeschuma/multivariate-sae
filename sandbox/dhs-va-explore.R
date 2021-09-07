# Austin Schumacher
# 9/7/2020
# test bangladesh dhs va data

rm(list = ls())
library(tidyverse); library(haven); library(INLA); library(sjmisc); library(svyVGAM); library(SUMMER);

setwd("/Users/austin/Desktop/india_mds_csmr/afg")
datadir <- "../../../Dropbox/csmr_india/afghanistan/data"

# read in data
datVA1 <- read_dta(paste0(datadir, "/AFVA66DT/VAFORM1.DTA"))
datVA2 <- read_dta(paste0(datadir, "/AFVA66DT/VAFORM2.DTA"))
datCHu5 <- read_dta(paste0(datadir, "/AFCH66DT/UNDER5.DTA"))
# datCHchildren <- read_dta(paste0(datadir, "/AFCH66DT/CHILDREN.DTA"))

## format u5 alive children data
datCHu5$nmonths <- datCHu5$qintc - datCHu5$q317c
datCHtmp <- datCHu5[datCHu5$nmonths < 60 & !is.na(datCHu5$nmonths),]
datCHtmp$agecat <- ifelse(datCHtmp$nmonths < 1, "0-28d", "1m-5y")
datCHtmp$icd10_condensed <- "alive"
datCHtmp$qhprov <- datCHtmp$qprov
datCHtmp$yod <- 1389
datCHtmp$qhweight <- datCHtmp$qweight
datCHtmp$qhtype <- datCHtmp$qtype
datCHtmp <- datCHtmp[, c("nmonths", "agecat", "yod", "icd10_condensed", "qhprov", "qhtype", "qhweight", "qhclust", "qhnumber")]

# formatting deaths data
datVA1$agedays <- datVA1$qv304o
datVA1$agedays[datVA1$qv304o > 28] <- NA
datVA1$agecat <- "0-28d"

datVA1$nmonths <- 1

# THIS IS TEMPORARY --- NEED TO DO THIS MORE ACCURATELY USING DAY OF SURVEY AND DAY OF DEATH
datVA2$agemonthsappx <- datVA2$qv304o
datVA2$agemonthsappx[datVA2$qv304v == 1] <- datVA2$agemonthsappx[datVA2$qv304v == 1] * 12
datVA2$agemonthsappx[datVA2$qv304v == 9] <- NA
datVA2$agemonthsappx[datVA2$qv304o %in% c(29, 40, 60, 98, 99)] <- NA
datVA2$agecat <- ifelse(datVA2$agemonthsappx < 60, "1m-5y", ">=5y")

datVA2$nmonths <- datVA2$agemonthsappx

# 1105 neonatal deaths (0-28 days)
# 831 1 - 60m deaths
# 157 5 - 12y deaths

# format VA death data
datVA1$yod <- datVA1$qv305z
datVA2$yod <- datVA2$qv305z

varnames <- c("nmonths", "agecat", "yod", "icd10_condensed", "qhprov", "qhtype", "qhweight", "qhclust", "qhnumber")
tmpdat <- rbind(datVA1 %>% dplyr::select(all_of(varnames)),
                datVA2 %>% dplyr::select(all_of(varnames)))

table(tmpdat$agecat)
table(tmpdat$yod)

# 4 years of data

table(tmpdat$icd10_condensed)

# label define ICD10_CONDENSED
# 1 "Certain infectious and parasitic diseases"
# 2 "Neoplasms"
# 3 "Diseases of the blood and blood-forming organs and certain disorders involving the immune mechanism"
# 4 "Endocrine, nutritional and metabolic diseases"
# 5 "Mental and behavioural disorders"
# 6 "Diseases of the nervous system"
# 7 "Diseases of the eye and adnexa"
# 8 "Diseases of the ear and mastoid process"
# 9 "Diseases of the circulatory system"
# 10 "Diseases of the respiratory system"
# 11 "Diseases of the digestive system"
# 12 "Diseases of the skin and subcutaneous tissue"
# 13 "Diseases of the musculoskeletal system and connective tissue"
# 14 "Diseases of the genitourinary system"
# 15 "Pregnancy, childbirth and the puerperium"
# 16 "Certain conditions originating in the perinatal period"
# 17 "Congenital malformations, deformations and chromosomal abnormalities"
# 18 "Symptoms, signs and abnormal clinical and laboratory findings, not elsewhere classified"
# 19 "External causes of morbidity and mortality"
# 20 "Other"
# 99 "Missing"

table(tmpdat$qhprov)

# crosstabs
table(tmpdat$agecat, tmpdat$yod)
table(tmpdat$agecat, tmpdat$icd10_condensed)
table(tmpdat$agecat, tmpdat$qhprov)
table(tmpdat$yod, tmpdat$icd10_condensed)
table(tmpdat$yod, tmpdat$qhprov)
table(tmpdat$icd10_condensed, tmpdat$qhprov)

# combine deaths and alives
dat <- rbind(datCHtmp, tmpdat)
dat$pweight <- 1/dat$qhweight
dat <- dat[!is.na(dat$nmonths), ]

#--------------
# individual level data
#--------------

# make a long dataset with an observation for each person-month
dat$pid <- 1:nrow(dat)
dat_long_list <- list()
length(dat_long_list) <- nrow(dat)
for (i in 1:nrow(dat)) {
    dati <- dat[dat$pid == i,]
    tmp <- dati[rep(1,dati$nmonths),]
    tmp$outcome <- "alive"
    tmp$outcome[length(tmp$outcome)] <- dati$icd10_condensed
    dat_long_list[[i]] <- tmp
}
dat_ind <- do.call(rbind, dat_long_list)

# TEMP: keep deaths all cause
# dat_ind$outcome_short <- ifelse(dat_ind$outcome == "16", "16",
#                                 ifelse(dat_ind$outcome == "alive", "alive", "other"))
dat_ind$died <- ifelse(dat_ind$outcome == "alive", 0, 1)

# reshape wide by cause
dat_ind_wide <- dat_ind %>% 
    to_dummy(died, suffix = "numeric") %>% 
    bind_cols(dat_ind)
dat_ind_wide$alive <- dat_ind_wide$`...1`
dat_ind_wide$dead <- dat_ind_wide$`...2`
dat_ind_wide$`...1` <- NULL
dat_ind_wide$`...2` <- NULL
# age groups


# fit preliminary model
my.svydesign <- survey::svydesign(ids = ~ qhclust,
                                  strata = ~ qhtype,
                                  nest = T, 
                                  weights = ~ qhweight, 
                                  data = dat_ind_wide)

mult.prem_noall0 <- svy_vglm(cbind(alive, dead) ~ 1, 
                             family = multinomial, 
                             design = my.svydesign)


# summer model
dat_ind_wide$year <- as.character(dat_ind_wide$yod)
dat_ind_wide$strata <- as.character(dat_ind_wide$qhtype)
# summer.mod <- getDirect(births = dat_ind_wide,
#                         years = c("1386","1387","1388","1389"),
#                         regionVar = NULL,
#                         weightsVar = "qhweight", 
#                         clusterVar = "~qhclust", 
#                         ageVar = NULL,
#                         timeVar = NULL)

#-------------
# Aggregated data
#-------------

# naive aggregation
dat_agg <- dat %>% group_by(qhprov, icd10_condensed) %>% summarise(death = n())

# variables
dat_agg$intercept1 <- ifelse(dat_agg$icd10_condensed == "alive", 1, 0)
dat_agg$intercept2 <- ifelse(dat_agg$icd10_condensed == "16", 1, 0)
dat_agg$intercept3 <- ifelse(dat_agg$icd10_condensed == "10", 1, 0)

dat_agg$prov_num <- as.numeric(as.factor(dat_agg$qhprov))
dat_agg$s1 <- dat_agg$s2 <- dat_agg$s3 <- dat_agg$s12 <- dat_agg$s13 <- dat_agg$s23 <- dat_agg$prov_num
tmp <- dat_agg %>% filter(icd10_condensed %in% c("alive", "16", "10"))

# testing coregionalization model
f1 <- death ~ icd10_condensed + 
    f(s1, model = "iid") + f(s2, model = "iid") + f(s3, model = "iid") + 
    f(s12, copy = "s1", fixed = FALSE) + 
    f(s13, copy = "s1", fixed = FALSE) + 
    f(s23, copy = "s2", fixed = FALSE) 

# model
result <- inla(f1, "poisson", 
               data = tmp)
