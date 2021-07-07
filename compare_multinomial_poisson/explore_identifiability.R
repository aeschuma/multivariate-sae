# Austin Schumacher
# 7/6/2020
# Check identifiability of Poisson trick formulations

rm(list = ls())

# if we only have two causes, two time points, 
# and we fit a model with just an intercept, 
# then the Poisson trick model has an RE for each time point, and a FE for each cause

X <- as.matrix(rbind(c(1,0,1,0),
                     c(1,0,0,1),
                     c(0,1,1,0),
                     c(0,1,0,1)))
solve(t(X)%*%X)

# if we have 3 time points, then
X <- as.matrix(rbind(c(1,0,0,1,0),
                     c(1,0,0,0,1),
                     c(0,1,0,1,0),
                     c(0,1,0,0,1),
                     c(0,0,1,1,0),
                     c(0,0,1,0,1)))
solve(t(X)%*%X)

# if we have 4 time points, then
X <- as.matrix(rbind(c(1,0,0,0,1,0),
                     c(1,0,0,0,0,1),
                     c(0,1,0,0,1,0),
                     c(0,1,0,0,0,1),
                     c(0,0,1,0,1,0),
                     c(0,0,1,0,0,1),
                     c(0,0,0,1,1,0),
                     c(0,0,0,1,0,1)))
solve(t(X)%*%X)