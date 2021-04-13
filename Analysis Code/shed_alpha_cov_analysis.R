#### SHEDDING AND ALPHA COVARIATION ANALYSIS

library(tidyverse)
library(parallel)

## set up
seeds <- floor(runif(20,1,1e5)) # set seeds
tmax <- 150

nocorr <- matrix(c(1,0,0,1), nrow=2, byrow=T)
negcorr <- matrix(c(1,-.5,-.5,1), nrow=2, byrow=T)
poscorr <- matrix(c(1,.5,.5,1), nrow=2, byrow=T)

alpha_shed_pars = c(c=.1, shed=.05, alpha=.1, gamma=.1, beta=.25, d=.1, 
                      b=2.5, bs=.01, varS=0.0099, sd_s=.1, sd_a=0.1) # R0=3.8


initial_state <- floor(c(S =unname(((contact_shed_pars["b"]-contact_shed_pars["d"])/contact_shed_pars["bs"]))-5, I=5, R=0))

###### RUN MULTIPLE SIMULATIONS ######
## no variation
source("GEM_SIR_noVar.R")
mclapply(seeds,
         function(s) gillespie.SIR.noVar(tmax, analytical_parms_shed, initial_state),
         mc.cores=4) -> out_no_var

## no corr
source("GEM_SIR_tau&alpha.cov.R")
mclapply(seeds,
         function(s) gillespie.SIR.cov_shedalpha(tmax, alpha_shed_pars, nocorr, initial_state),
         mc.cores=4) -> out_nocov_shed_alpha

## positive corr
source("GEM_SIR_tau&alpha.cov.R")
mclapply(seeds,
         function(s) gillespie.SIR.cov_shedalpha(tmax, alpha_shed_pars, poscorr, initial_state),
         mc.cores=4) -> out_poscov_shed_alpha

## negative corr
source("GEM_SIR_tau&alpha.cov.R")
mclapply(seeds,
         function(s) gillespie.SIR.cov_shedalpha(tmax, alpha_shed_pars, negcorr, initial_state),
         mc.cores=4) -> out_negcov_shed_alpha

###### FORMAT OUTPUT ######
## no cov
## infected
timeSeq <- 0:150
storeMatrix.cov.SA.I <- array(NA, dim=c(length(timeSeq),length(out_nocov_shed_alpha)))
for (j in 1:length(out_nocov_shed_alpha)) {
  o <- out_nocov_shed_alpha[[j]]
  storeMatrix.cov.SA.I[,j] <- o[,3] 
}

## susceptible
storeMatrix.cov.SA.S <- array(NA, dim=c(length(timeSeq),length(out_nocov_shed_alpha)))
for (j in 1:length(out_nocov_shed_alpha)) {
  o <- out_nocov_shed_alpha[[j]]
  storeMatrix.cov.SA.S[,j] <- o[,2] 
}

## recovered
storeMatrix.cov.SA.R <- array(NA, dim=c(length(timeSeq),length(out_nocov_shed_alpha)))
for (j in 1:length(out_nocov_shed_alpha)) {
  o <- out_nocov_shed_alpha[[j]]
  storeMatrix.cov.SA.R[,j] <- o[,4] 
}

## positive cov
## infected
timeSeq <- 0:150
storeMatrix.poscov.SA.I <- array(NA, dim=c(length(timeSeq),length(out_poscov_shed_alpha)))
for (j in 1:length(out_poscov_shed_alpha)) {
  o <- out_poscov_shed_alpha[[j]]
  storeMatrix.poscov.SA.I[,j] <- o[,3] 
}

## susceptible
storeMatrix.poscov.SA.S <- array(NA, dim=c(length(timeSeq),length(out_poscov_shed_alpha)))
for (j in 1:length(out_poscov_shed_alpha)) {
  o <- out_poscov_shed_alpha[[j]]
  storeMatrix.poscov.SA.S[,j] <- o[,2] 
}

## recovered
storeMatrix.poscov.SA.R <- array(NA, dim=c(length(timeSeq),length(out_poscov_shed_alpha)))
for (j in 1:length(out_poscov_shed_alpha)) {
  o <- out_poscov_shed_alpha[[j]]
  storeMatrix.poscov.SA.R[,j] <- o[,4]
}

## neg cov
## infected
timeSeq <- 0:150
storeMatrix.negcov.SA.I <- array(NA, dim=c(length(timeSeq),length(out_negcov_shed_alpha)))
for (j in 1:length(out_negcov_shed_alpha)) {
  o <- out_negcov_shed_alpha[[j]]
  storeMatrix.negcov.SA.I[,j] <- o[,3] 
}

## susceptible
storeMatrix.negcov.SA.S <- array(NA, dim=c(length(timeSeq),length(out_negcov_shed_alpha)))
for (j in 1:length(out_negcov_shed_alpha)) {
  o <- out_negcov_shed_alpha[[j]]
  storeMatrix.negcov.SA.S[,j] <- o[,2] 
}

## recovered
storeMatrix.negcov.SA.R <- array(NA, dim=c(length(timeSeq),length(out_negcov_shed_alpha)))
for (j in 1:length(out_negcov_shed_alpha)) {
  o <- out_negcov_shed_alpha[[j]]
  storeMatrix.negcov.SA.R[,j] <- o[,4] 
}

###### PLOT ######
par(mfrow=c(1,3))

## no cov
plot(0:150, apply(storeMatrix.cov.SA.I, 1, mean), col="red", lwd=1.75, type="l", ylim=c(0,250), ylab="N", xlab="Time", main="No Cov")
lines(0:150, apply(storeMatrix.cov.SA.S, 1, mean), col="blue", lwd=1.75, type="l", ylim=c(0,150), ylab="N", xlab="Time")
lines(0:150, apply(storeMatrix.cov.SA.R, 1, mean), col="green", lwd=1.75, type="l", ylim=c(0,150), ylab="N", xlab="Time")
#abline(h=63, lty=2)
legend("topright",legend=c("S","I","R"),fill=c("blue","red","green"), cex=0.25)

## pos cov
plot(0:150, apply(storeMatrix.poscov.SA.I, 1, mean), col="red", lwd=1.75, type="l", ylim=c(0,250), ylab="N", xlab="Time", main="Pos Cov")
lines(0:150, apply(storeMatrix.poscov.SA.S, 1, mean), col="blue", lwd=1.75, type="l", ylim=c(0,150), ylab="N", xlab="Time")
lines(0:150, apply(storeMatrix.poscov.SA.R, 1, mean), col="green", lwd=1.75, type="l", ylim=c(0,150), ylab="N", xlab="Time")
#abline(h=63, lty=2)
legend("topright",legend=c("S","I","R"),fill=c("blue","red","green"), cex=0.25)

## neg cov
plot(0:150, apply(storeMatrix.negcov.SA.I, 1, mean), col="red", lwd=1.75, type="l", ylim=c(0,250), ylab="N", xlab="Time", main="Neg Cov")
lines(0:150, apply(storeMatrix.negcov.SA.S, 1, mean), col="blue", lwd=1.75, type="l", ylim=c(0,150), ylab="N", xlab="Time")
lines(0:150, apply(storeMatrix.negcov.SA.R, 1, mean), col="green", lwd=1.75, type="l", ylim=c(0,150), ylab="N", xlab="Time")
#abline(h=63, lty=2)
legend("topright",legend=c("S","I","R"),fill=c("blue","red","green"), cex=0.25)

