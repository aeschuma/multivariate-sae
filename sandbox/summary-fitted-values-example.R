rm(list = ls())
library(SUMMER); library(INLA);

## test #1: linear trend data, model fit with intercept only
n <- 25
b <- 1
x <- (1:n)/n
N <- 100
y <- rbinom(n, N, expit(x*b))

mod.bin.test <- inla(y ~ 1,
                     control.fixed = list(prec.intercept = 0),
                     family = "binomial",
                     Ntrials = N,
                     data = list(y = y,
                                 N = N),
                     control.predictor = list(compute = TRUE),
                     control.compute = list(config = TRUE))

# the following two plots show that the median summary.fitted.values have exactly the same relative distribution as 
# the empirical data points
#   i.e. the plots look exactly the same except for the scale of the y-axis
plot(mod.bin.test$summary.fitted.values$`0.5quant`)
plot(y/N)

# but with data and predictions plotted on the same graph, we get this:
plot(y/N)
points(mod.bin.test$summary.fitted.values$`0.5quant`, col = "red", pch = 19)


## test #2: completely random data, model fit with intercept only
test_n <- 1000
test_pi <- 0.1
test_N <- 500
test_data2 <- data.frame(ydat = rbinom(test_n, test_N, test_pi),
                         Ndat = rep(test_N, test_n))
mod.bin.test2 <- inla(ydat ~ 1,
                      Ntrials = Ndat,
                      control.fixed = list(prec.intercept = 0),
                      data = test_data2,
                      family = "binomial",
                      control.predictor = list(compute = TRUE),
                      control.compute = list(config = TRUE))

# the same issue as above is here too
plot(mod.bin.test2$summary.fitted.values$`0.5quant`)
plot(test_data2$ydat/test_data2$Ndat)

# data and predictions plotted on the same graph
plot(test_data2$ydat/test_data2$Ndat)
points(mod.bin.test2$summary.fitted.values$`0.5quant`, col = "red", pch = 19)
