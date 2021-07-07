#  to do

rm(list = ls())
setwd('~/Desktop/survey-csmf/sandbox')

## Libraries ####
library(SUMMER)
# help(package = "SUMMER", help_type = "html")
library(classInt)
library(RColorBrewer)
library(dplyr)
library(tidyr)
library(rgdal)
library(scales)
library(survey)
library(tidyverse)

#### Parameters ####

iso <- "bgd"
start_year <- 1990
end_year <- 2020
svy_year <- 2017

## Load data ####

load(paste0('../../../Dropbox/dissertation_2/survey-csmf/data/bgd_dhs2017/data/', 'births-file_', iso, '.rda'))

mod.dat$years <- as.numeric(as.character(mod.dat$time))
dat.years <- sort(unique(mod.dat$years))
beg.years <- seq(start_year,end_year,5)
end.years <- beg.years + 4
periods <- paste(beg.years, end.years, sep = "-")
mod.dat$period <- as.character(cut(mod.dat$years, breaks = c(beg.years, beg.years[length(beg.years)]+5),
                                   include.lowest = T, right = F, labels = periods))

# IGME

igme <- read_csv("../../../Dropbox/dissertation_2/survey-csmf/data/un-igme/UN IGME Child Mortality and Stillbirth Estimates.csv") %>%
    filter(`Geographic area` == "Bangladesh" & Indicator == "Under-five mortality rate" & 
               `Series Name` %in% c("Demographic and Health Survey 2017-2018 (Direct)",
                                    "UN IGME estimate") &
               Sex == "Total") %>%
    select(`Series Name`, `Series Year`, `TIME_PERIOD`, `OBS_VALUE`) %>%
    mutate(year = as.numeric(substr(TIME_PERIOD, 1, 4)))

#### National Direct ####

mod.dat$v005 <- mod.dat$v005/1e6


direct.natl <-  getDirect(mod.dat, periods,
                          regionVar = "v024",
                          timeVar = "period", 
                          clusterVar =  "~v001",
                          ageVar = "age", Ntrials = "total",
                          weightsVar = "v005",national.only = T)
direct.natl$survey <- 1
direct.natl$surveyYears <- svy_year

direct.natl.yearly <-  getDirect(mod.dat, start_year:end_year,
                                 regionVar = "v024",
                                 timeVar = "years", 
                                 clusterVar =  "~v001",
                                 ageVar = "age", Ntrials = "total",
                                 weightsVar = "v005",national.only = T)
direct.natl.yearly$survey <- 1
direct.natl.yearly$surveyYears <- svy_year


save(direct.natl, file = paste0('/Users/austin/Dropbox/dissertation_2/survey-csmf/results/bgd/', 'natl-direct_', iso, '.rda'))
save(direct.natl.yearly, file = paste0('/Users/austin/Dropbox/dissertation_2/survey-csmf/results/bgd/', 'natl-direct-yearly_', iso, '.rda'))

## Admin1 Direct ####

direct.admin1 <-  getDirect(mod.dat, periods,
                            regionVar = "v024",
                            timeVar = "period", 
                            clusterVar =  "~v001",
                            ageVar = "age", Ntrials = "total",
                            weightsVar = "v005",national.only = F)
direct.admin1$survey <- 1
direct.admin1$surveyYears <- svy_year

save(direct.admin1, file = paste0('/Users/austin/Dropbox/dissertation_2/survey-csmf/results/bgd/', 'admin1-direct_', iso, '.rda'))

#### Polygon Plots ####

## load polygon shape files
poly.adm0 <- readOGR(dsn = "/Users/austin/Dropbox/dissertation_2/survey-csmf/data/shapefiles/bgd_adm_bbs_20201113_SHP",
                     layer = "bgd_admbnda_adm0_bbs_20201113")
poly.adm0$region <- "All"
poly.adm1 <- readOGR(dsn = "/Users/austin/Dropbox/dissertation_2/survey-csmf/data/shapefiles/bgd_adm_bbs_20201113_SHP",
                     layer = "bgd_admbnda_adm1_bbs_20201113")
poly.adm1$region <- tolower(poly.adm1$ADM1_EN)

## National ##
med.palette <- brewer.pal(5, name = "Purples")
med.int <- classIntervals(round(direct.natl$logit.est, 2),
                          n = 5, style = 'jenks')
med.col <- findColours(med.int, med.palette)

pdf(paste0('/Users/austin/Dropbox/dissertation_2/survey-csmf/plots/bgd/','direct-natl-poly_', iso, '.pdf'), height = 9, width = 16)

par(mfrow = c(2, 4))
for(year in periods){
  idx <- which(direct.natl$years == year)
  plot(poly.adm0, border = F, col = med.col[idx],
       axes = F, main = year)
}
plot(NA, xlim = c(0,1), ylim = c(0,1), axes = F, xlab = "", ylab = "")
legend(x = "center",inset = 0,
       legend = names(attr(med.col, 'table')),
       fill = med.palette, cex= 1, horiz = FALSE, bty = 'n')

dev.off()

## Admin1 ##

med.palette <- brewer.pal(n = 5, name = "Purples")
med.int <- classIntervals(round(direct.admin1$logit.est, 2),
                          n = 5, style = 'jenks')
med.col <- findColours(med.int, med.palette)

pdf(paste0('/Users/austin/Dropbox/dissertation_2/survey-csmf/plots/bgd/','direct-admin1-poly_', iso, '.pdf'), height = 9, width = 16)

par(mfrow = c(2, 4))
for(year in periods){
  tmp <- poly.adm1
  for (reg in tmp$region) {
      tmp$col[tmp$region == reg] <- med.col[direct.admin1$years == year & direct.admin1$region == reg]
  }
  plot(tmp, border = F, col = tmp$col,
       axes = F, main = year)
}
plot(NA, xlim = c(0,1), ylim = c(0,1), axes = F, xlab = "", ylab = "")
legend(x = "center",inset = 0,
       legend = names(attr(med.col, 'table')),
       fill = med.palette, cex= 1, horiz = FALSE, bty = 'n')

dev.off()

## Spaghetti Plots ####
plot.years <- seq(start_year, end_year, 5)

## National ##
pdf(paste0('/Users/austin/Dropbox/dissertation_2/survey-csmf/plots/bgd/','direct-natl-spaghetti_', iso, '.pdf'), width = 12, height = 7)
plot.max <- max(direct.natl$upper+.025, na.rm = T)
plot.min <- min(direct.natl$lower - 0.025, na.rm = T)

pane.years <- plot.years
tmp <- direct.natl
tmp$width <- tmp$upper - tmp$lower
tmp$cex2 <- median(tmp$width, na.rm = T)/tmp$width
tmp$cex2[tmp$cex2 > 6] <- 6

par(mfrow = c(1,1))
plot(NA,
     xlab = "Year", ylab = "U5MR",
     ylim = c(plot.min, plot.max),
     xlim = c(start_year, end_year),
     type = 'n', col = 'indianred', lwd = 2,
     main = paste0(iso, " national direct"))

lines(pane.years, tmp$mean, col = 'indianred',
      lwd = 2)
points(pane.years, tmp$mean,
       col = alpha('indianred', 0.35),
       cex = tmp$cex2, pch = 19)
lines(igme$year[igme$`Series Name` == "UN IGME estimate" & igme$year <= 2020 & igme$year >= 1990], 
      igme$OBS_VALUE[igme$`Series Name` == "UN IGME estimate" & igme$year <= 2020 & igme$year >= 1990]/1000,
      col = "orchid4", lty = 2, lwd = 2)
points(igme$year[igme$`Series Name` == "UN IGME estimate" & igme$year <= 2020 & igme$year >= 1990], 
      igme$OBS_VALUE[igme$`Series Name` == "UN IGME estimate" & igme$year <= 2020 & igme$year >= 1990]/1000,
      col = alpha("orchid4", 0.5), pch = 17)
lines(igme$year[igme$`Series Name` == "Demographic and Health Survey 2017-2018 (Direct)"], 
      igme$OBS_VALUE[igme$`Series Name` == "Demographic and Health Survey 2017-2018 (Direct)"]/1000,
      col = "dodgerblue", lty = 4, lwd = 2)
points(igme$year[igme$`Series Name` == "Demographic and Health Survey 2017-2018 (Direct)"], 
      igme$OBS_VALUE[igme$`Series Name` == "Demographic and Health Survey 2017-2018 (Direct)"]/1000,
      col = alpha("dodgerblue", 0.5), pch = 15)
legend("topright", c("Direct", "UN-IGME DHS 2017", "UN-IGME final est"),
       lty = c(1, 4, 2), lwd = 2, col = c("indianred", "dodgerblue", "orchid4"))

dev.off()

pdf(paste0('/Users/austin/Dropbox/dissertation_2/survey-csmf/plots/bgd/','direct-natl-yearly-spaghetti_', iso, '.pdf'), width = 12, height = 7)
plot.max <- max(direct.natl.yearly$upper+.025, na.rm = T)
plot.min <- min(direct.natl.yearly$lower - 0.025, na.rm = T)

pane.years <- jitter(start_year:end_year)
tmp <- direct.natl.yearly
tmp$width <- tmp$upper - tmp$lower
tmp$cex2 <- median(tmp$width, na.rm = T)/tmp$width
tmp$cex2[tmp$cex2 > 6] <- 6

par(mfrow = c(1,1))
plot(NA,
     xlab = "Year", ylab = "U5MR",
     ylim = c(plot.min, plot.max),
     xlim = c(start_year, end_year),
     type = 'n', col = 'orchid3', lwd = 2,
     main = paste0(iso, " national direct, yearly"))

lines(pane.years, tmp$mean, col = 'orchid3',
      lwd = 2)

points(pane.years, tmp$mean,
       col = alpha('orchid3', 0.35),
       cex = tmp$cex2, pch = 19)

dev.off()


## Admin1 ##

pdf(paste0('/Users/austin/Dropbox/dissertation_2/survey-csmf/plots/bgd/','direct-admin1-spaghetti_', iso, '.pdf'), width = 12, height = 7)
par(mfrow = c(2,4))
cols <- brewer.pal(8, 'Set3')
areas <- unique(direct.admin1$region[direct.admin1$region != "All"])
for(area in 1:length(areas)){
  tmp.area <- direct.admin1[direct.admin1$region == areas[area],]
  tmp.area$width <- tmp.area$upper - tmp.area$lower
  tmp.area$cex2 <- median(tmp.area$width, na.rm = T)/tmp.area$width
  tmp.area$cex2[tmp.area$cex2 > 6] <- 6
  
  plot.max <- max(tmp.area$upper+.025, na.rm = T)
  
  pane.years <- jitter(plot.years)
  
  plot(NA,
       xlab = "Year", ylab = "U5MR",
       ylim = c(0, plot.max),
       xlim = c(start_year, end_year),
       type = 'l', col = cols[area], lwd = 2,
       main = paste0(areas[area], ": direct admin1"))
      
  lines(pane.years, tmp.area$mean, cex = tmp.area$cex2,
        type = 'l', col = cols[area], lwd = 2)
  
  points(pane.years, tmp.area$mean, pch = 19,
         col = alpha(cols[area], 0.35),
         cex = tmp$cex2)
}
dev.off()

## gg spaghetti
direct.admin1$plot.year <- as.numeric(substr(direct.admin1$years, 1, 4)) + 2.5
ggplot(direct.admin1 %>% filter(region != "All"), aes(x = plot.year, y = logit.est, col = region)) + geom_line()
ggsave(paste0('/Users/austin/Dropbox/dissertation_2/survey-csmf/plots/bgd/','direct-admin1-spaghetti-all_', iso, '.pdf'),
       width = 10, height = 6)
