#  to do

rm(list = ls())
setwd('~/Desktop/survey-csmf/sandbox')

## Libraries ####
library(SUMMER)
help(package = "SUMMER", help_type = "html")
library(classInt)
library(RColorBrewer)
library(dplyr)
library(tidyr)
library(rgdal)
library(scales)
library(survey)
library(tidyverse)
library(spdep)

#### Parameters ####

iso <- "bgd"
start_year <- 1990
end_year <- 2020
svy_year <- 2017

## Load data ####

load(paste0('/Users/austin/Dropbox/dissertation_2/survey-csmf/data/bgd_dhs2017/data/', 'births-file_', iso, '.rda'))

mod.dat$years <- as.numeric(as.character(mod.dat$time))
dat.years <- sort(unique(mod.dat$years))
beg.years <- seq(start_year,end_year,5)
end.years <- beg.years + 4
periods <- paste(beg.years, end.years, sep = "-")
mod.dat$period <- as.character(cut(mod.dat$years, breaks = c(beg.years, beg.years[length(beg.years)]+5),
                                   include.lowest = T, right = F, labels = periods))

## direct.natl
load(paste0('/Users/austin/Dropbox/dissertation_2/survey-csmf/results/bgd/', 'natl-direct_', iso, '.rda'))

# direct.natl.yearly
load(paste0('/Users/austin/Dropbox/dissertation_2/survey-csmf/results/bgd/', 'natl-direct-yearly_', iso, '.rda'))

# direct.admin1
load(paste0('/Users/austin/Dropbox/dissertation_2/survey-csmf/results/bgd/', 'admin1-direct_', iso, '.rda'))

## load polygon shape files
poly.adm0 <- readOGR(dsn = "/Users/austin/Dropbox/dissertation_2/survey-csmf/data/shapefiles/bgd_adm_bbs_20201113_SHP",
                     layer = "bgd_admbnda_adm0_bbs_20201113")
poly.adm0$region <- "All"
poly.adm1 <- readOGR(dsn = "/Users/austin/Dropbox/dissertation_2/survey-csmf/data/shapefiles/bgd_adm_bbs_20201113_SHP",
                     layer = "bgd_admbnda_adm1_bbs_20201113")
poly.adm1$region <- tolower(poly.adm1$ADM1_EN)

# adjacency matrix
admin1.mat <- poly2nb(SpatialPolygons(poly.adm1@polygons))
admin1.mat <- nb2mat(admin1.mat, zero.policy = TRUE)
colnames(admin1.mat) <- rownames(admin1.mat) <- poly.adm1$region

#### National Model ####
fit.natl <- fitINLA(direct.natl, geo = NULL, Amat = NULL,
                    year_label = c(periods),
                    year_range = c(start_year, end_year + 5), is.yearly = F)
res.natl <- getSmoothed(fit.natl, year_range = c(beg.year, end.year+5),
                        year_label = c(periods, proj.per))
res.natl$years.num <- seq(start_year+2, end_year+5, 5)
save(res.natl, file = paste0('/Users/austin/Dropbox/dissertation_2/survey-csmf/results/bgd/', 'natl-smoothed-direct_', iso, '.rda'))

fit.natl.yearly <- fitINLA(direct.natl.yearly, geo = NULL, Amat = NULL,
                           year_label = as.character(start_year:(end_year + 5)),
                           year_range = c(start_year, end_year + 5), is.yearly = F)
res.natl.yearly <- getSmoothed(fit.natl.yearly, year_range = c(start_year, end_year + 5),
                               year_label = as.character(start_year:(end_year + 5)))
res.natl.yearly$years.num <- start_year:(end_year + 5)
save(res.natl.yearly, file = paste0('/Users/austin/Dropbox/dissertation_2/survey-csmf/results/bgd/', 
                                    'natl-smoothed-direct-yearly_', iso, '.rda'))

#### Admin1 Model ####
fit.admin1 <- fitINLA(direct.admin1, geo = poly.adm1, Amat = admin1.mat,
                      year_label = periods,
                      year_range = c(1990, 2024), is.yearly = F)
res.admin1 <- getSmoothed(fit.admin1, Amat = admin1.mat,
                          year_range = c(1990, 2024),
                          year_label = periods)
res.admin1$years.num <- seq(start_year+2, end_year+5, 5)[match(res.admin1$years, periods)]
save(res.admin1, file = paste0('/Users/austin/Dropbox/dissertation_2/survey-csmf/results/bgd/', 
                               'admin1-smoothed-direct_', iso, '.rda'))

#### Polygon Plots ####

## National ##
med.palette <- brewer.pal(5, name = "Purples")
med.int <- classIntervals(round(res.natl$logit.median, 2),
                          n = 5, style = 'jenks')
med.col <- findColours(med.int, med.palette)

pdf(paste0('/Users/austin/Dropbox/dissertation_2/survey-csmf/plots/bgd/','smoothed-direct-natl-poly_', iso, '.pdf'), height = 9, width = 16)

par(mfrow = c(2, 4))
for(year in periods){
  idx <- which(res.natl$years == year)
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
med.int <- classIntervals(round(res.admin1$logit.median, 2),
                          n = 5, style = 'jenks')
med.col <- findColours(med.int, med.palette)

pdf(paste0('/Users/austin/Dropbox/dissertation_2/survey-csmf/plots/bgd/','smoothed-direct-admin1-poly_', iso, '.pdf'), height = 9, width = 16)

par(mfrow = c(2, 4))
for(year in periods){
  tmp <- poly.adm1
  for (reg in tmp$region) {
      tmp$col[tmp$region == reg] <- med.col[res.admin1$years == year & res.admin1$region == reg]
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
plot.years <- seq(start_year + 2.5, end_year + 2.5, 5)

## National ##
pdf(paste0('/Users/austin/Dropbox/dissertation_2/survey-csmf/plots/bgd/','smoothed-direct-natl-spaghetti_', iso, '.pdf'), width = 12, height = 7)
plot.max <- max(res.natl$logit.upper+.025, na.rm = T)
plot.min <- min(res.natl$logit.lower - 0.025, na.rm = T)

pane.years <- jitter(plot.years)
tmp <- res.natl
tmp$width <- tmp$logit.upper - tmp$logit.lower
tmp$cex2 <- median(tmp$width, na.rm = T)/tmp$width
tmp$cex2[tmp$cex2 > 6] <- 6

par(mfrow = c(1,1))
plot(NA,
     xlab = "Year", ylab = "logit U5MR",
     ylim = c(plot.min, plot.max),
     xlim = range(plot.years),
     type = 'n', col = 'indianred', lwd = 2,
     main = paste0(iso, " national smoothed direct"))

lines(pane.years, tmp$logit.median, col = 'indianred',
      lwd = 2)
points(pane.years, tmp$logit.median,
       col = alpha('indianred', 0.35),
       cex = tmp$cex2, pch = 19)

dev.off()

pdf(paste0('/Users/austin/Dropbox/dissertation_2/survey-csmf/plots/bgd/','smoothed-direct-natl-yearly-spaghetti_', iso, '.pdf'), width = 12, height = 7)
plot.max <- max(res.natl.yearly$logit.upper+.025, na.rm = T)
plot.min <- min(res.natl.yearly$logit.lower - 0.025, na.rm = T)

pane.years <- jitter(res.natl.yearly$years.num)
tmp <- res.natl.yearly
tmp$width <- tmp$logit.upper - tmp$logit.lower
tmp$cex2 <- median(tmp$width, na.rm = T)/tmp$width
tmp$cex2[tmp$cex2 > 6] <- 6

par(mfrow = c(1,1))
plot(NA,
     xlab = "Year", ylab = "logit U5MR",
     ylim = c(plot.min, plot.max),
     xlim = range(pane.years),
     type = 'n', col = 'orchid3', lwd = 2,
     main = paste0(iso, " national smoothed direct, yearly"))

lines(pane.years, tmp$logit.median, col = 'orchid3',
      lwd = 2)

points(pane.years, tmp$logit.median,
       col = alpha('orchid3', 0.35),
       cex = tmp$cex2, pch = 19)

dev.off()


## Admin1 ##

pdf(paste0('/Users/austin/Dropbox/dissertation_2/survey-csmf/plots/bgd/','smoothed-direct-admin1-spaghetti_', iso, '.pdf'), width = 12, height = 7)
par(mfrow = c(2,4))
cols <- brewer.pal(8, 'Set3')
areas <- unique(res.admin1$region[res.admin1$region != "All"])
for(area in 1:length(areas)){
  tmp.area <- res.admin1[res.admin1$region == areas[area],]
  tmp.area$width <- tmp.area$logit.upper - tmp.area$logit.lower
  tmp.area$cex2 <- median(tmp.area$width, na.rm = T)/tmp.area$width
  tmp.area$cex2[tmp.area$cex2 > 6] <- 6
  
  plot.max <- max(tmp.area$logit.upper+.025, na.rm = T)
  plot.min <- min(tmp.area$logit.lower-.025, na.rm = T)
  
  pane.years <- jitter(plot.years)
  
  plot(NA,
       xlab = "Year", ylab = "logit U5MR",
       ylim = c(plot.min, plot.max),
       xlim = range(plot.years),
       type = 'l', col = cols[area], lwd = 2,
       main = paste0(areas[area], ": direct admin1"))
      
  lines(pane.years, tmp.area$logit.median, cex = tmp.area$cex2,
        type = 'l', col = cols[area], lwd = 2)
  
  points(pane.years, tmp.area$logit.median, pch = 19,
         col = alpha(cols[area], 0.35),
         cex = tmp$cex2)
}
dev.off()

## gg spaghetti
ggplot(res.admin1 %>% filter(region != "All"), aes(x = years.num, y = logit.median, col = region)) + geom_line()
ggsave(paste0('/Users/austin/Dropbox/dissertation_2/survey-csmf/plots/bgd/','smoothed-direct-admin1-spaghetti-all_', iso, '.pdf'),
       width = 10, height = 6)

