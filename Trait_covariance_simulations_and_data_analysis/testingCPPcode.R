library(Rcpp)
library(RcppEigen)

bivariateNormalR <-function(traitmeans, traitsds, corr, randnorm) {
  p <- length(traitmeans);
  mu <- matrix(0,p,1);
  sigma <- matrix(0,p,p);
  ## Generate Lognormal random variables
  ## To produce a distribution with mean equal to traitmean
  ## and variance equal to traitsd^2, the parameters of
  ## the corresponding lognormal distribution are:
  for (i in 1:p) {
    mu[i] = log(traitmeans[i]^2 / sqrt(traitsds[i]^2+traitmeans[i]^2))
    sigma[i,i] = log(1 + traitsds[i]^2/traitmeans[i]^2);
  }
  ## generate the correlation matrix
  corr <- matrix(c(1,corr,corr,1),nrow=2,byrow=TRUE)
  ## generate the covariance matrix
  Cov <- sqrt(sigma)%*%corr%*%sqrt(sigma);
  Cov.Eigen = eigen(Cov);
  Cov.EigenVal <- Cov.Eigen$values*(abs(Cov.Eigen$values) >= 1e-12);
  RootCov = (Cov.Eigen$vectors)%*%diag(sqrt(Cov.EigenVal))%*%t(Cov.Eigen$vectors);
  X <- matrix(randnorm, nrow = N, ncol = p);
  for (i in 1:N){
    X[i,] <- exp(mu + RootCov%*%X[i,])
  }
  colnames(X) <- names(traitmeans) ## give the columns names
  return(X);
}

## RcppEigen function
sourceCpp("multivariateLognormal.cpp")
## Rcpp function
sourceCpp("SIRcov.cpp")

## Confirm the R and RcppEigen code agree
traitmeans = c(1,1)
traitsds = c(1, 1)
randnorm = rnorm(2)
bivariateNormalR(traitmeans, traitsds, 0.5, randnorm)
bivariateNormalEigen(traitmeans, traitsds, 0.5, randnorm)
bivariateNormal(traitmeans, traitsds, 0.5)

## Test for speed: RcppEigen code is ~35x faster, Cpp code is ~40x faster!
library(microbenchmark)
microbenchmark(R = bivariateNormalR(traitmeans, traitsds, 0.5, rnorm(2)), 
               Eigen = bivariateNormalEigen(traitmeans, traitsds, 0.5, rnorm(2)),
               Cpp = bivariateNormal(traitmeans, traitsds, 0.5),
               times=10000)

sirGEM(traitmeans, traitsds, 0.5)


SIRdet <- function(t, y, params) {
  S = y[1];
  I = y[2];
  R = y[3];
  b = params["b"];
  bs = params["bs"];
  d = params["d"]
  c = params["c"]
  s = params["s"]
  a = params["a"]
  g = params["g"]
  
  dSdt = (b-bs*(S+I+R))*(S+I+R) - c*s/(1+s)*S*I - d*S
  dIdt = c*s/(1+s)*S*I - (d+a+g)*I
  dRdt = g*I - d*R
  
  return(list(c(dSdt,dIdt,dRdt)))
}
library(deSolve)  

## No variation simulation
params = c(b=2.5, bs=0.01, d=0.1, 
           c=0.1, s=1/9, a=0.1, g=0.1,
           cSD=0.0001, sSD=0.0001,
           corr=0.5)
x = c(S=230, I=10, R=0)

## Deterministic run
det = ode(x, seq(0,100,0.1), SIRdet, params)
  
## Stochastic run in Rcpp
sourceCpp("SIRcov.cpp")
out = SIRcovCS(params, x, 0.1)

## Stochastic run in R
source("GEM_SIR_cov_storage.R")
c(  c=0.1,    shed=1/9,      alpha=0.1,    gamma=0.1, 
    sd_c=0.0001, sd_shed=0.0001, sd_alpha=0.0001, sd_gamma=0.0001, 
    b=2.5, d=.1, bs=.01) -> paramsR
gillespie.SIR.cov_storage(tmax=100, 
                          params=paramsR, 
                          corr=0.5, 
                          x=x, 
                          covParams=c("c","shed"),
                          seed=101347) -> outR

## Compare dynamics to confirm that the stochastic simulation works!
par(mfrow=c(1,3), oma=rep(0.5,4), mar=c(3, 3, 0.5, 0.5))
plot(out[[1]], out[[2]][,1], type='l')
lines(det[,c(1,2)], col=2, lwd=2)
lines(outR[[1]][,c(1,2)], col=4, lwd=2, lty=2)

plot(out[[1]], out[[2]][,2], type='l')
lines(det[,c(1,3)], col=2, lwd=2)
lines(outR[[1]][,c(1,3)], col=4, lwd=2, lty=2)

plot(out[[1]], out[[2]][,3], type='l')
lines(det[,c(1,4)], col=2, lwd=2)
lines(outR[[1]][,c(1,4)], col=4, lwd=2, lty=2)

## Compare the speed of the Rcpp code compared to the raw R code
tic <- Sys.time()
for (i in 1:10) 
  gillespie.SIR.cov_storage(tmax=100, 
                            params=paramsR, 
                            corr=0.5, 
                            x=x, 
                            covParams=c("c","shed"),
                            seed=101347) -> outR
toc <- Sys.time()
toc - tic

paramsCS = c(b=2.5, bs=0.01, d=0.1, 
             c=0.1, s=1/9, a=0.1, g=0.1,
             cSD=0.0001, sSD=0.0001,
             corr=0.5)
tic <- Sys.time()
for (i in 1:10) 
  out = SIRcovCS(params, c(S=230, I=10, R=0), 50)
toc <- Sys.time()
toc - tic

## Stochastic run in Rcpp
sourceCpp("SIRcov.cpp")
x = c(S=230, I=10, R=0)

## No variation simulations (should all be very similar)
paramsCS = c(b=2.5, bs=0.01, d=0.1, 
             c=0.1, s=1/9, a=0.1, g=0.1,
             cSD=0.0001, sSD=0.0001,
             corr=0.5)
tic <- Sys.time()
for (i in 1:10) 
  outCS = SIRcovCS(paramsCS, x, 100)
toc <- Sys.time()
toc - tic

tic <- Sys.time()
for (i in 1:10) 
  outCS = SIRcovCS2(paramsCS, x, 100)
toc <- Sys.time()
toc - tic


paramsCA = c(b=2.5, bs=0.01, d=0.1, 
             c=0.1, s=1/9, a=0.1, g=0.1,
             cSD=0.0001, aSD=0.0001,
             corr=0.5)
outCA = SIRcovCA(paramsCA, x, 50)

paramsCG = c(b=2.5, bs=0.01, d=0.1, 
             c=0.1, s=1/9, a=0.1, g=0.1,
             cSD=0.0001, gSD=0.0001,
             corr=0.5)
outCG = SIRcovCG(paramsCG, x, 50)

paramsSA = c(b=2.5, bs=0.01, d=0.1, 
             c=0.1, s=1/9, a=0.1, g=0.1,
             sSD=0.0001, aSD=0.0001,
             corr=0.5)
outSA = SIRcovSA(paramsSA, x, 50)

paramsSG = c(b=2.5, bs=0.01, d=0.1, 
             c=0.1, s=1/9, a=0.1, g=0.1,
             sSD=0.0001, gSD=0.0001,
             corr=0.5)
outSG = SIRcovSG(paramsSG, x, 50)

paramsAG = c(b=2.5, bs=0.01, d=0.1, 
             c=0.1, s=1/9, a=0.1, g=0.1,
             aSD=0.0001, gSD=0.0001,
             corr=0.5)
outAG = SIRcovAG(paramsAG, x, 50)

paramsC = c(b=2.5, bs=0.01, d=0.1, 
             c=0.1, s=1/9, a=0.1, g=0.1,
             cSD=0.0001)
outC = SIRvarC(paramsC, x, 50)

paramsS = c(b=2.5, bs=0.01, d=0.1, 
            c=0.1, s=1/9, a=0.1, g=0.1,
            sSD=0.0001)
outS = SIRvarS(paramsS, x, 50)

paramsA = c(b=2.5, bs=0.01, d=0.1, 
            c=0.1, s=1/9, a=0.1, g=0.1,
            aSD=0.0001)
outA = SIRvarA(paramsA, x, 50)

paramsG = c(b=2.5, bs=0.01, d=0.1, 
            c=0.1, s=1/9, a=0.1, g=0.1,
            gSD=0.0001)
outG = SIRvarG(paramsG, x, 50)

params0 = c(b=2.5, bs=0.01, d=0.1, 
            c=0.1, s=1/9, a=0.1, g=0.1)
out0 = SIRnovar(params0, x, 50)


par(mfrow=c(1,1), mar=c(4,4,0.5,0.5), oma=rep(0.5,4))
plot(outCS[[1]], outCS[[2]][,2], type='l', xlab="Time", ylab="No. infected")
lines(outCA[[1]], outCA[[2]][,2], col=2)
lines(outCG[[1]], outCG[[2]][,2], col=3)
lines(outSA[[1]], outSA[[2]][,2], col=4)
lines(outSG[[1]], outSG[[2]][,2], col=5)
lines(outAG[[1]], outAG[[2]][,2], col=6)
lines(outC[[1]], outC[[2]][,2], col=7)
lines(outS[[1]], outS[[2]][,2], col=8)
lines(outA[[1]], outA[[2]][,2], col=9)
lines(outG[[1]], outG[[2]][,2], col=10)
lines(out0[[1]], out0[[2]][,2], col=11)


