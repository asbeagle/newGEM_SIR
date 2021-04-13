#### SHEDDING AND CONTACT COVARIATION ANALYSIS CODE

library(tidyverse)
library(parallel)

## set up
seeds <- floor(runif(20,1,1e5)) # set seeds
tmax <- 150

nocorr <- matrix(c(1,0,0,1), nrow=2, byrow=T)
negcorr <- matrix(c(1,-.5,-.5,1), nrow=2, byrow=T)
poscorr <- matrix(c(1,.5,.5,1), nrow=2, byrow=T)


contact_shed_pars = c(c=.1, shed=.05, alpha=.1, gamma=.1, beta=.25, d=.1, 
                      b=2.5, bs=.01, varS=0.0099, epsilon=0.01, sd_s=.1, sd_c=0.01) # R0=3.8


initial_state <- floor(c(S =unname(((contact_shed_pars["b"]-contact_shed_pars["d"])/contact_shed_pars["bs"]))-5, I=5, R=0))

###### RUN MULTIPLE SIMULATIONS ######
## no variation
source("GEM_SIR_noVar.R")
mclapply(seeds,
         function(s) gillespie.SIR.noVar(tmax, analytical_parms_shed, initial_state),
         mc.cores=4) -> out_no_var

## no corr
source("GEM_SIR_contact&tau.cov.R")
mclapply(seeds,
         function(s) gillespie.SIR.cov_cshed(tmax, contact_shed_pars, nocorr, initial_state),
         mc.cores=4) -> out_nocov_contact_shed

## positive corr
source("GEM_SIR_contact&tau.cov.R")
mclapply(seeds,
         function(s) gillespie.SIR.cov_cshed(tmax, contact_shed_pars, poscorr, initial_state),
         mc.cores=4) -> out_poscov_contact_shed

## negative corr
source("GEM_SIR_contact&tau.cov.R")
mclapply(seeds,
         function(s) gillespie.SIR.cov_cshed(tmax, contact_shed_pars, negcorr, initial_state),
         mc.cores=4) -> out_negcov_contact_shed

###### FORMAT OUTPUT ######
## no correlation
## infected
timeSeq <- 0:150
storeMatrix.cov.CS.I <- array(NA, dim=c(length(timeSeq),length(out_nocov_contact_shed)))
for (j in 1:length(out_nocov_contact_shed)) {
  o <- out_nocov_contact_shed[[j]]
  storeMatrix.cov.CS.I[,j] <- o[,3] 
}

## susceptible
storeMatrix.cov.CS.S <- array(NA, dim=c(length(timeSeq),length(out_nocov_contact_shed)))
for (j in 1:length(out_nocov_contact_shed)) {
  o <- out_nocov_contact_shed[[j]]
  storeMatrix.cov.CS.S[,j] <- o[,2] 
}

## recovered
storeMatrix.cov.CS.R <- array(NA, dim=c(length(timeSeq),length(out_nocov_contact_shed)))
for (j in 1:length(out_nocov_contact_shed)) {
  o <- out_nocov_contact_shed[[j]]
  storeMatrix.cov.CS.R[,j] <- o[,4] 
}

## positive cov
## infected
timeSeq <- 0:150
storeMatrix.poscov.CS.I <- array(NA, dim=c(length(timeSeq),length(out_poscov_contact_shed)))
for (j in 1:length(out_poscov_contact_shed)) {
  o <- out_poscov_contact_shed[[j]]
  storeMatrix.poscov.CS.I[,j] <- o[,3] 
}

## susceptible
storeMatrix.poscov.CS.S <- array(NA, dim=c(length(timeSeq),length(out_poscov_contact_shed)))
for (j in 1:length(out_poscov_contact_shed)) {
  o <- out_poscov_contact_shed[[j]]
  storeMatrix.poscov.CS.S[,j] <- o[,2] 
}

## recovered
storeMatrix.poscov.CS.R <- array(NA, dim=c(length(timeSeq),length(out_poscov_contact_shed)))
for (j in 1:length(out_poscov_contact_shed)) {
  o <- out_poscov_contact_shed[[j]]
  storeMatrix.poscov.CS.R[,j] <- o[,4]
}

## neg cov
## infected
timeSeq <- 0:150
storeMatrix.negcov.CS.I <- array(NA, dim=c(length(timeSeq),length(out_negcov_contact_shed)))
for (j in 1:length(out_negcov_contact_shed)) {
  o <- out_negcov_contact_shed[[j]]
  storeMatrix.negcov.CS.I[,j] <- o[,3] 
}

## susceptible
storeMatrix.negcov.CS.S <- array(NA, dim=c(length(timeSeq),length(out_negcov_contact_shed)))
for (j in 1:length(out_negcov_contact_shed)) {
  o <- out_negcov_contact_shed[[j]]
  storeMatrix.negcov.CS.S[,j] <- o[,2] 
}

## recovered
storeMatrix.negcov.CS.R <- array(NA, dim=c(length(timeSeq),length(out_negcov_contact_shed)))
for (j in 1:length(out_negcov_contact_shed)) {
  o <- out_negcov_contact_shed[[j]]
  storeMatrix.negcov.CS.R[,j] <- o[,4] 
}

###### PLOT ######
par(mfrow=c(1,3))

## no cov
plot(0:150, apply(storeMatrix.cov.CS.I, 1, mean), col="red", lwd=1.75, type="l", ylim=c(0,250), ylab="N", xlab="Time", main="No Cov")
lines(0:150, apply(storeMatrix.cov.CS.S, 1, mean), col="blue", lwd=1.75, type="l", ylim=c(0,150), ylab="N", xlab="Time")
lines(0:150, apply(storeMatrix.cov.CS.R, 1, mean), col="green", lwd=1.75, type="l", ylim=c(0,150), ylab="N", xlab="Time")
#abline(h=63, lty=2)
legend("topright",legend=c("S","I","R"),fill=c("blue","red","green"), cex=0.25)

## pos cov
plot(0:150, apply(storeMatrix.poscov.CS.I, 1, mean), col="red", lwd=1.75, type="l", ylim=c(0,250), ylab="N", xlab="Time", main="Pos Cov")
lines(0:150, apply(storeMatrix.poscov.CS.S, 1, mean), col="blue", lwd=1.75, type="l", ylim=c(0,150), ylab="N", xlab="Time")
lines(0:150, apply(storeMatrix.poscov.CS.R, 1, mean), col="green", lwd=1.75, type="l", ylim=c(0,150), ylab="N", xlab="Time")
#abline(h=63, lty=2)
legend("topright",legend=c("S","I","R"),fill=c("blue","red","green"), cex=0.25)

## neg cov
plot(0:150, apply(storeMatrix.negcov.CS.I, 1, mean), col="red", lwd=1.75, type="l", ylim=c(0,250), ylab="N", xlab="Time", main="Neg Cov")
lines(0:150, apply(storeMatrix.negcov.CS.S, 1, mean), col="blue", lwd=1.75, type="l", ylim=c(0,150), ylab="N", xlab="Time")
lines(0:150, apply(storeMatrix.negcov.CS.R, 1, mean), col="green", lwd=1.75, type="l", ylim=c(0,150), ylab="N", xlab="Time")
#abline(h=63, lty=2)
legend("topright",legend=c("S","I","R"),fill=c("blue","red","green"), cex=0.25)

