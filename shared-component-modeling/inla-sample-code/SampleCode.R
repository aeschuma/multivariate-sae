setwd("./inla-sample-code")

#### Estimate unobserved disease counts ####
dat=read.table("SampleData.csv", header=T, sep=",")

# Obtain Empirical Bayes prior estimates 
alphaMat=do.call(rbind, by(dat, dat$strat, function(subdat){
  tmp=data.frame(strat=subdat$strat[1])
  tmp[, c("aES", "aCS", "aOS")]=colSums(subdat[, c("zES", "zCS", "zOS")])/colSums(dat[, c("zES", "zCS", "zOS")])
  tmp[, c("aEM", "aCM", "aOM")]=colSums(subdat[, c("zEM", "zCM", "zOM")])/colSums(dat[, c("zEM", "zCM", "zOM")])
  return(tmp)
}))
alphaMat$a.S=rowSums(alphaMat[, c("aES", "aCS", "aOS")])
alphaMat$a.M=rowSums(alphaMat[, c("aEM", "aCM", "aOM")])

# Estimate disease counts
dat$yES.hat=with(dat, zES+(yS-kS)*(alphaMat[dat$strat, "aES"]+zES)/(alphaMat[dat$strat, "a.S"]+kS) )
dat$yEM.hat=with(dat, zEM+(yM-kM)*(alphaMat[dat$strat, "aEM"]+zEM)/(alphaMat[dat$strat, "a.M"]+kM) )
dat$yCS.hat=with(dat, zCS+(yS-kS)*(alphaMat[dat$strat, "aCS"]+zCS)/(alphaMat[dat$strat, "a.S"]+kS) )
dat$yCM.hat=with(dat, zCM+(yM-kM)*(alphaMat[dat$strat, "aCM"]+zCM)/(alphaMat[dat$strat, "a.M"]+kM) )
dat$yOS.hat=with(dat, zOS+(yS-kS)*(alphaMat[dat$strat, "aOS"]+zOS)/(alphaMat[dat$strat, "a.S"]+kS) )
dat$yOM.hat=with(dat, zOM+(yM-kM)*(alphaMat[dat$strat, "aOM"]+zOM)/(alphaMat[dat$strat, "a.M"]+kM) )

dat$yE.hat=dat$yES.hat+dat$yEM.hat
dat$yC.hat=dat$yCS.hat+dat$yCM.hat
dat$yO.hat=dat$yOS.hat+dat$yOM.hat

# Estimate variances
dat$var.yES=with(dat, (yS-kS)*(alphaMat[dat$strat, "aES"]+zES)*(yS+alphaMat[dat$strat, "a.S"])*(alphaMat[dat$strat, "a.S"]+kS-alphaMat[dat$strat, "aES"]-zES)/((kS+alphaMat[dat$strat, "a.S"]+1)*(kS+alphaMat[dat$strat, "a.S"])^2) )
dat$var.yEM=with(dat, (yM-kM)*(alphaMat[dat$strat, "aEM"]+zEM)*(yM+alphaMat[dat$strat, "a.M"])*(alphaMat[dat$strat, "a.M"]+kM-alphaMat[dat$strat, "aEM"]-zEM)/((kM+alphaMat[dat$strat, "a.M"]+1)*(kM+alphaMat[dat$strat, "a.M"])^2) )
dat$var.yCS=with(dat, (yS-kS)*(alphaMat[dat$strat, "aCS"]+zCS)*(yS+alphaMat[dat$strat, "a.S"])*(alphaMat[dat$strat, "a.S"]+kS-alphaMat[dat$strat, "aCS"]-zCS)/((kS+alphaMat[dat$strat, "a.S"]+1)*(kS+alphaMat[dat$strat, "a.S"])^2) )
dat$var.yCM=with(dat, (yM-kM)*(alphaMat[dat$strat, "aCM"]+zCM)*(yM+alphaMat[dat$strat, "a.M"])*(alphaMat[dat$strat, "a.M"]+kM-alphaMat[dat$strat, "aCM"]-zCM)/((kM+alphaMat[dat$strat, "a.M"]+1)*(kM+alphaMat[dat$strat, "a.M"])^2) )
dat$var.yOS=with(dat, (yS-kS)*(alphaMat[dat$strat, "aOS"]+zOS)*(yS+alphaMat[dat$strat, "a.S"])*(alphaMat[dat$strat, "a.S"]+kS-alphaMat[dat$strat, "aOS"]-zOS)/((kS+alphaMat[dat$strat, "a.S"]+1)*(kS+alphaMat[dat$strat, "a.S"])^2) )
dat$var.yOM=with(dat, (yM-kM)*(alphaMat[dat$strat, "aOM"]+zOM)*(yM+alphaMat[dat$strat, "a.M"])*(alphaMat[dat$strat, "a.M"]+kM-alphaMat[dat$strat, "aOM"]-zOM)/((kM+alphaMat[dat$strat, "a.M"]+1)*(kM+alphaMat[dat$strat, "a.M"])^2) )

dat$var.yE=dat$var.yES+dat$var.yEM
dat$var.yC=dat$var.yCS+dat$var.yCM
dat$var.yO=dat$var.yOS+dat$var.yOM

## Compute covariances ##
dat$covECS.hat=with(dat, -(yS-kS)*(alphaMat[dat$strat, "aES"]+zES)*(alphaMat[dat$strat, "aCS"]+zCS)*(yS+alphaMat[dat$strat, "a.S"])/((kS+alphaMat[dat$strat, "a.S"]+1)*(kS+alphaMat[dat$strat, "a.S"])^2) )
dat$covECM.hat=with(dat, -(yM-kM)*(alphaMat[dat$strat, "aEM"]+zEM)*(alphaMat[dat$strat, "aCM"]+zCM)*(yM+alphaMat[dat$strat, "a.M"])/((kM+alphaMat[dat$strat, "a.M"]+1)*(kM+alphaMat[dat$strat, "a.M"])^2) )

## Estimate Reference Probabilities ##
pEj.hat=coef(glm(yE.hat~offset(log(Nj))-1+as.factor(strat), family=poisson, data=dat))
pCj.hat=coef(glm(yC.hat~offset(log(Nj))-1+as.factor(strat), family=poisson, data=dat))
pOj.hat=coef(glm(yO.hat~offset(log(Nj))-1+as.factor(strat), family=poisson, data=dat))


## Collapse for weekly data ##
yGt=0
yGt=data.frame(do.call(rbind, by(dat, dat$wk, function(wkdat){
  c(wk=wkdat$wk[1], colSums(wkdat[, -c(1:2)]))
})))

## Estimate Expected Numbers ##
Nj=dat[dat$wk==1, "Nj"]
yGt$EE.hat=rep(crossprod(exp(pEj.hat), Nj), dim(yGt)[1])
yGt$EC.hat=rep(crossprod(exp(pCj.hat), Nj), dim(yGt)[1])
yGt$EO.hat=rep(crossprod(exp(pOj.hat), Nj), dim(yGt)[1])

## Estimate thetaGts ##
yGt$thetaEt.hat=with(yGt, yE.hat/EE.hat)
yGt$thetaCt.hat=with(yGt, yC.hat/EC.hat)
yGt$thetaOt.hat=with(yGt, yO.hat/EO.hat)

## Estimate variances of thetaGts ## 
yGt$varEt=with(yGt, var.yE/(EE.hat^2))
yGt$varCt=with(yGt, var.yC/(EC.hat^2))
yGt$varOt=with(yGt, var.yO/(EO.hat^2))

## Estimate covariance of thetaEt and thetaCt ##
yGt$covECt.hat=with(yGt, (covECS.hat+covECM.hat)/(EE.hat*EC.hat) )

# Get log thetas and variances #
yGt$lEt.hat=with(yGt, log(yE.hat/EE.hat))
yGt$lCt.hat=with(yGt, log(yC.hat/EC.hat))
yGt$lOt.hat=with(yGt, log(yO.hat/EO.hat))

yGt$var.lEt.hat=with(yGt, varEt/(thetaEt.hat^2))
yGt$var.lCt.hat=with(yGt, varCt/(thetaCt.hat^2))
yGt$var.lOt.hat=with(yGt, varOt/(thetaOt.hat^2))

## Estimate covariance of log thetaEt and log thetaCt ##
yGt$covlECt.hat=with(yGt, covECt.hat/(thetaEt.hat*thetaCt.hat) )


#### Fit joint model using INLA ####
library(INLA)
Nt=156 ## Number of times
temp=dat$temp[1:Nt] ## get temperature covariate

# set up 2D model for INLA
N <- 2*Nt
weights = rep(1, N)

# define data in long form
data=list(lGt=c(yGt$lEt.hat, yGt$lCt.hat), 
          wkE=c(1:Nt, rep(NA, Nt)),
          wkE2=c(1:Nt, rep(NA, Nt)),
          wkC=c(rep(NA, Nt), 1:Nt),
          wkC2=c(rep(NA, Nt), 1:Nt),
          interE=c(rep(1,Nt), rep(NA,Nt)), 
          interC=c(rep(NA,Nt),rep(1,Nt)), 
          tempE=c(temp, rep(NA, Nt)),
          tempC=c(rep(NA, Nt), temp)
)

## Define the index-vectors ii.1 ii.2 etc, which are the
## index's for the iid2d-model at timepoint 1, 2, ...
for(j in 1:Nt) {
  itmp = numeric(N)
  itmp[] = NA
  itmp[j] = 1
  itmp[j+Nt] = 2
  data = c(list(itmp), data)
  names(data)[1] = paste("ii.", j, sep="")
}

## we now have to add one 'iid2d' model for each observation pair,
## since their cov.matrix is different. we have to do this
## automatically... here I add numbers directly for simplicity
add=""
for(j in 1:Nt) {
  corr = yGt[j, "covlECt.hat"]/sqrt(prod(yGt[j, c("var.lEt.hat", "var.lCt.hat")]))
  init.precE = log(1/yGt[j, "var.lEt.hat"])
  init.precC = log(1/yGt[j, "var.lCt.hat"])
  
  add = paste(add, paste(" + 
                         f(", paste("ii.", j, sep=""), ", weights, model=\"iid2d\", n=2,
                         hyper = list(
                         prec1 = list(
                         initial =", init.precE,", 
                         fixed = TRUE),
                         prec2 = list(
                         initial =", init.precC,", 
                         fixed = TRUE),
                         cor = list(
                         initial = log((1+", corr, ")/(1-", corr, ")), 
                         fixed = TRUE)))"))
  
}



## define hyper parameters
hypers=c(list(prec = list(param = c(10, 1e-2) )), 
         list(prec = list(param = c(10, 1e-2) ))  )

## Model
S1=matrix(1:Nt, ncol=Nt, nrow=1)
formula.lGt = lGt ~ -1 + interE + interC + f(wkE, model="rw2", hyper=hypers[1], scale.model=F, constr=T, extraconstr=list(A=S1, e=0)) + f(wkC, model="rw2", hyper=hypers[2], scale.model=F, constr=T, extraconstr=list(A=S1, e=0)) + wkE2 + wkC2 + tempE + tempC

## Add iid2d pair to formula
formula.lGt=update(formula.lGt, as.formula(paste(". ~ . ", add)))


## Fit model
mod= inla(formula.lGt, data=data, 
          family = "gaussian",
          control.family = list(hyper = list(prec = list(initial = 10, fixed=T))), 
          control.predictor=list(compute=T),
          control.compute=list(config=T),
          control.inla=list(lincomb.derived.correlation.matrix=T),
          control.fixed=list(prec=list(default=0.001), correlation.matrix=T) )

