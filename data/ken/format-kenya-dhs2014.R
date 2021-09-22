# Format Kenya DHS 2014 data for modeling HAZ and WAZ

rm(list = ls())

library(SUMMER)
library(tidyverse)
library(spdep)
library(geosphere)
library(haven)
library(knitr)
library(kableExtra)
library(magrittr)
library(rgdal)

# Read in data

svyfile <- "/Users/austin/Dropbox/dissertation_2/survey-csmf/data/ken_dhs2014/data/KE_2014_DHS_08252021_2356_143252/KEBR72DT/KEBR72FL.DTA"

dat <- read_stata(svyfile)

dat %<>% filter(!(hw5 %in% c(9998, 9999)) & !(hw8 %in% c(9998, 9999))) %>%
    mutate(HAZ = hw5/100, WAZ = hw8/100, 
           region = v024, urban_rural = v025, strata = v023,
           weights = v005,  cluster = v001, tmp_admin1 = scounty) %>% 
    select(HAZ, WAZ,
           cluster, region, strata, urban_rural, weights, tmp_admin1) %>%
    filter(!is.na(HAZ) & !is.na(WAZ))

# data plotting
ggplot(dat, aes(x=HAZ)) + 
    geom_histogram() +
    ggtitle("Height for age z-score distribution for KEN DHS 2014") +
    theme_light()

ggplot(dat, aes(x=WAZ)) + 
    geom_histogram() +
    ggtitle("Weight for age z-score distribution for KEN DHS 2014") +
    theme_light()

table(dat$tmp_admin1)

# geo data
poly.layer.adm0 <- paste('gadm36', 'KEN',
                         '0', sep = "_")
poly.layer.adm1 <- paste('gadm36', 'KEN',
                         '1', sep = "_")
poly.layer.adm2 <- paste('gadm36', 'KEN',
                         '2', sep = "_")

poly.path <- "/Users/austin/Dropbox/dissertation_2/survey-csmf/data/shapefiles/shapeFiles_gadm/ken"
poly.adm1 <- readOGR(dsn = poly.path,
                     layer = as.character(poly.layer.adm1))
proj4string(poly.adm1)

# adjacency matrix
admin1.mat <- poly2nb(SpatialPolygons(poly.adm1@polygons))
admin1.mat <- nb2mat(admin1.mat, zero.policy = TRUE, style = "B")
colnames(admin1.mat) <- 
    rownames(admin1.mat) <- paste0("admin1_",
                                   1:dim(admin1.mat)[1])
admin1.names <- data.frame(GADM = poly.adm1@data$NAME_1,
                           Internal = rownames(admin1.mat))

# assign lat and long
points.path <- "/Users/austin/Dropbox/dissertation_2/survey-csmf/data/ken_dhs2014/data/KE_2014_DHS_08252021_2356_143252/KEGE71FL"
points <- readOGR(dsn = points.path,
                  layer = "KEGE71FL")

wrong.points <- which(points@data$LATNUM == 0.0 &
                          points@data$LONGNUM == 0.0)
plot(points[-wrong.points,])

dat.tmp <- dat
dat.tmp <- dat.tmp[!(dat.tmp$cluster %in% points@data$DHSCLUST[wrong.points]),]
if (sum(points@data$DHSCLUST[wrong.points] %in% unique(dat.tmp$cluster) != 0)) stop("deletion of wrong points unsuccessful")

dat.tmp$LONGNUM <- dat.tmp$LATNUM <- NA
for(i in 1:dim(points)[1]){
    dat.tmp$LATNUM[dat.tmp$cluster == points@data$DHSCLUST[i]] <- points@data$LATNUM[i]
    dat.tmp$LONGNUM[dat.tmp$cluster == points@data$DHSCLUST[i]] <- points@data$LONGNUM[i]
}

miss <- which(dat.tmp$LATNUM == 0 & dat.tmp$LONGNUM == 0)
if(length(miss != 0)){
    dat.tmp <- dat.tmp[-miss,]
}

points.frame <- as.data.frame(dat.tmp[,c("LONGNUM", "LATNUM")])
points.frame <- SpatialPoints(points.frame)
poly.over.adm1 <- SpatialPolygons(poly.adm1@polygons)
proj4string(points.frame) <- proj4string(poly.over.adm1) <- 
    proj4string(poly.adm1) 
admin1.key <- over(points.frame, poly.over.adm1)
miss.frame.adm1 <- unique(points.frame@coords[which(is.na(admin1.key)),])

if(dim(miss.frame.adm1)[1] != 0){
    miss.poly.adm1 <- dist2Line( miss.frame.adm1, poly.over.adm1)
    
    for(i in 1:dim(miss.poly.adm1)[1]){
        long.ids <- which(points.frame@coords[,c("LONGNUM")] %in% miss.frame.adm1[i,1])
        lat.ids <- which(points.frame@coords[,c("LATNUM")] %in% miss.frame.adm1[i,2])
        ids <- intersect(long.ids, lat.ids)
        
        # ids[(length(ids)/2 + 1):length(ids)] <- ids[(length(ids)/2 + 1):length(ids)] - dim(points.frame@coords)[1]
        # ids <- unique(ids)
        admin1.key[ids] <- rep(miss.poly.adm1[i, 'ID'], length(ids))
    }
}

dat.tmp$admin1 <- admin1.key
dat.tmp$admin1.char <- paste0("admin1_", admin1.key)
dat.tmp$admin1.name <- as.character(poly.adm1@data$NAME_1)[admin1.key]

# create node lists
n_regions <- length(unique(dat.tmp$admin1))
node1 <- c()
node2 <- c()
n_edges <- 0
for (i in 1:n_regions) {
    for (j in 1:i) {
        if(admin1.mat[i,j] != 0) {
            node1 <- c(node1, i)
            node2 <- c(node2, j)
            n_edges <- n_edges + 1
        }
    }
}
node.info <- list(node1 = node1, node2 = node2, n_edges = n_edges)

# scaling factor for ICAR models

## make INLA adjacency matrix
adj.matrix <- sparseMatrix(i=node.info$node1,j=node.info$node2,x=1,symmetric=TRUE)

## The ICAR precision matrix (note! This is singular)
Q <- Diagonal(n_regions, rowSums(adj.matrix)) - adj.matrix

## Add a small jitter to the diagonal for numerical stability (optional but recommended)
Q_pert <- Q + Diagonal(n_regions) * max(diag(Q)) * sqrt(.Machine$double.eps)

## Compute the diagonal elements of the covariance matrix subject to the 
## constraint that the entries of the ICAR sum to zero.
## See the inla.qinv function help for further details.
Q_inv <- inla.qinv(Q_pert, constr=list(A = matrix(1,1,n_regions),e=0))

## Compute the geometric mean of the variances, which are on the diagonal of Q_inv
scaling_factor <- exp(mean(log(diag(Q_inv))))

# format data
dat <- dat.tmp
dat <- dat[order(dat$admin1),]
colnames(admin1.mat) <- 1:ncol(admin1.mat)
rownames(admin1.mat) <- 1:nrow(admin1.mat)

# save data
save(dat, poly.adm1, admin1.mat, node.info, scaling_factor,
     file = "/Users/austin/Dropbox/dissertation_2/survey-csmf/data/ken_dhs2014/data/haz-waz-kenDHS2014.rda")
