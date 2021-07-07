# Austin
# 10/16/2020
# Simulation example for coregionalization model

rm(list=ls())

library(tidyverse); library(mvtnorm); library(INLA);

# USE ICAR SPATIAL FOR THE Z1 and Z2

# simulate data
lambda <- 0.25
alpha <- c(1, 2)
n <- 10000

z1 <- rnorm(n , 0, 1)
z2 <- rnorm(n , 0, 1)

e1 <- rnorm(n , 0, 1)
e2 <- rnorm(n , 0, 1)

y1 <- alpha[1] + z1 
y2 <- alpha[2] + lambda*z1 + z2 

dat <- data.frame(y = c(y1, y2),
                  int1 = c(rep(1, n), rep(0, n)),
                  int2 = c(rep(0, n), rep(1, n)),
                  s1 = c(1:n, rep(NA, n)),
                  s2 = c(rep(NA, n), 1:n),
                  s12 = c(rep(NA, n), 1:n))

# fit model

# testing coregionalization model
f1 <- y ~ 0 + int1 + int2 +
  f(s1, model = "iid") + f(s2, model = "iid") + 
  f(s12, copy = "s1", fixed = FALSE) 

# model
result <- inla(f1, family = "gaussian", 
               data = dat)
summary(result)
result$summary.hyperpar
