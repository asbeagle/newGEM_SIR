#### SHEDDING AND GAMMA COVARIATION ANALYSIS

library(tidyverse)
library(parallel)

## set up
seeds <- floor(runif(20,1,1e5)) # set seeds
tmax <- 150

nocorr <- matrix(c(1,0,0,1), nrow=2, byrow=T)
negcorr <- matrix(c(1,-.5,-.5,1), nrow=2, byrow=T)
poscorr <- matrix(c(1,.5,.5,1), nrow=2, byrow=T)

gamma_shed_pars = c(c=.1, shed=.05, alpha=.1, gamma=.1, beta=.25, d=.1, 
                    b=2.5, bs=.01, varS=0.0099, sd_s=.1, sd_g=0.1) # R0=3.8


initial_state <- floor(c(S =unname(((contact_shed_pars["b"]-contact_shed_pars["d"])/contact_shed_pars["bs"]))-5, I=5, R=0))

###### RUN MULTIPLE SIMULATIONS ######
## no variation
source("GEM_SIR_noVar.R")
mclapply(seeds,
         function(s) gillespie.SIR.noVar(tmax, analytical_parms_shed, initial_state),
         mc.cores=4) -> out_no_var

## no corr
source("GEM_SIR_tau&gamma.cov.R")
mclapply(seeds,
         function(s) gillespie.SIR.cov_shedgamma(tmax, gamma_shed_pars, nocorr, initial_state),
         mc.cores=4) -> out_nocov_shed_gamma

## positive corr
source("GEM_SIR_tau&gamma.cov.R")
mclapply(seeds,
         function(s) gillespie.SIR.cov_shedgamma(tmax, gamma_shed_pars, poscorr, initial_state),
         mc.cores=4) -> out_poscov_shed_gamma

## negative corr
source("GEM_SIR_tau&gamma.cov.R")
mclapply(seeds,
         function(s) gillespie.SIR.cov_shedalpha(tmax, gamma_shed_pars, negcorr, initial_state),
         mc.cores=4) -> out_negcov_shed_gamma

###### FORMAT OUTPUT ######
## no cov
## infected
timeSeq <- 0:150
storeMatrix.cov.SG.I <- array(NA, dim=c(length(timeSeq),length(out_nocov_shed_gamma)))
for (j in 1:length(out_nocov_shed_gamma)) {
  o <- out_nocov_shed_gamma[[j]]
  storeMatrix.cov.SG.I[,j] <- o[,3] 
}

## susceptible
storeMatrix.cov.SG.S <- array(NA, dim=c(length(timeSeq),length(out_nocov_shed_gamma)))
for (j in 1:length(out_nocov_shed_gamma)) {
  o <- out_nocov_shed_gamma[[j]]
  storeMatrix.cov.SG.S[,j] <- o[,2] 
}

## recovered
storeMatrix.cov.SG.R <- array(NA, dim=c(length(timeSeq),length(out_nocov_shed_gamma)))
for (j in 1:length(out_nocov_shed_gamma)) {
  o <- out_nocov_shed_gamma[[j]]
  storeMatrix.cov.SG.R[,j] <- o[,4] 
}

## positive cov
## infected
timeSeq <- 0:150
storeMatrix.poscov.SG.I <- array(NA, dim=c(length(timeSeq),length(out_poscov_shed_gamma)))
for (j in 1:length(out_poscov_shed_gamma)) {
  o <- out_poscov_shed_gamma[[j]]
  storeMatrix.poscov.SG.I[,j] <- o[,3] 
}

## susceptible
storeMatrix.poscov.SG.S <- array(NA, dim=c(length(timeSeq),length(out_poscov_shed_gamma)))
for (j in 1:length(out_poscov_shed_gamma)) {
  o <- out_poscov_shed_gamma[[j]]
  storeMatrix.poscov.SG.S[,j] <- o[,2] 
}

## recovered
storeMatrix.poscov.SG.R <- array(NA, dim=c(length(timeSeq),length(out_poscov_shed_gamma)))
for (j in 1:length(out_poscov_shed_gamma)) {
  o <- out_poscov_shed_gamma[[j]]
  storeMatrix.poscov.SG.R[,j] <- o[,4]
}

## neg cov
## infected
timeSeq <- 0:150
storeMatrix.negcov.SG.I <- array(NA, dim=c(length(timeSeq),length(out_negcov_shed_gamma)))
for (j in 1:length(out_negcov_shed_gamma)) {
  o <- out_negcov_shed_gamma[[j]]
  storeMatrix.negcov.SG.I[,j] <- o[,3] 
}

## susceptible
storeMatrix.negcov.SG.S <- array(NA, dim=c(length(timeSeq),length(out_negcov_shed_gamma)))
for (j in 1:length(out_negcov_shed_gamma)) {
  o <- out_negcov_shed_gamma[[j]]
  storeMatrix.negcov.SG.S[,j] <- o[,2] 
}

## recovered
storeMatrix.negcov.SG.R <- array(NA, dim=c(length(timeSeq),length(out_negcov_shed_gamma)))
for (j in 1:length(out_negcov_shed_gamma)) {
  o <- out_negcov_shed_gamma[[j]]
  storeMatrix.negcov.SG.R[,j] <- o[,4] 
}

###### PLOT ######
par(mfrow=c(1,3))

## no cov
plot(0:150, apply(storeMatrix.cov.SG.I, 1, mean), col="red", lwd=1.75, type="l", ylim=c(0,250), ylab="N", xlab="Time", main="No Cov")
lines(0:150, apply(storeMatrix.cov.SG.S, 1, mean), col="blue", lwd=1.75, type="l", ylim=c(0,150), ylab="N", xlab="Time")
lines(0:150, apply(storeMatrix.cov.SG.R, 1, mean), col="green", lwd=1.75, type="l", ylim=c(0,150), ylab="N", xlab="Time")
#abline(h=63, lty=2)
legend("topright",legend=c("S","I","R"),fill=c("blue","red","green"), cex=0.25)

## pos cov
plot(0:150, apply(storeMatrix.poscov.SG.I, 1, mean), col="red", lwd=1.75, type="l", ylim=c(0,250), ylab="N", xlab="Time", main="Pos Cov")
lines(0:150, apply(storeMatrix.poscov.SG.S, 1, mean), col="blue", lwd=1.75, type="l", ylim=c(0,150), ylab="N", xlab="Time")
lines(0:150, apply(storeMatrix.poscov.SG.R, 1, mean), col="green", lwd=1.75, type="l", ylim=c(0,150), ylab="N", xlab="Time")
#abline(h=63, lty=2)
legend("topright",legend=c("S","I","R"),fill=c("blue","red","green"), cex=0.25)

## neg cov
plot(0:150, apply(storeMatrix.negcov.SG.I, 1, mean), col="red", lwd=1.75, type="l", ylim=c(0,250), ylab="N", xlab="Time", main="Neg Cov")
lines(0:150, apply(storeMatrix.negcov.SG.S, 1, mean), col="blue", lwd=1.75, type="l", ylim=c(0,150), ylab="N", xlab="Time")
lines(0:150, apply(storeMatrix.negcov.SG.R, 1, mean), col="green", lwd=1.75, type="l", ylim=c(0,150), ylab="N", xlab="Time")
#abline(h=63, lty=2)
legend("topright",legend=c("S","I","R"),fill=c("blue","red","green"), cex=0.25)