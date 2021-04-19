#### CONTACT AND GAMMA COVARIATION ANALYSIS CODE
library(tidyverse)
library(parallel)

## set up
seeds <- floor(runif(20,1,1e5)) # set seeds
tmax <- 150

nocorr <- matrix(c(1,0,0,1), nrow=2, byrow=T)
negcorr <- matrix(c(1,-.5,-.5,1), nrow=2, byrow=T)
poscorr <- matrix(c(1,.5,.5,1), nrow=2, byrow=T)


contact_gamma_pars = c(c=.1, shed=.05, alpha=.1, gamma=.1, beta=.25, d=.1, 
                      b=2.5, bs=.01, sd_c=1e-6, sd_g=1e-6) # R0=3.8


initial_state <- floor(c(S =unname(((contact_gamma_pars["b"]-contact_gamma_pars["d"])/contact_gamma_pars["bs"]))-5, I=5, R=0))

###### RUN MULTIPLE SIMULATIONS ######
## no variation
source("GEM_SIR_noVar.R")
mclapply(seeds,
         function(s) gillespie.SIR.noVar(tmax, contact_gamma_pars, initial_state),
         mc.cores=4) -> out_no_var

## no corr
source("GEM_SIR_contact&gamma.cov.R")
mclapply(seeds,
         function(s) gillespie.SIR.cov_cgamma(tmax, contact_gamma_pars, nocorr, initial_state),
         mc.cores=4) -> out_nocov_contact_gamma

## positive corr
source("GEM_SIR_contact&gamma.cov.R")
mclapply(seeds,
         function(s) gillespie.SIR.cov_cgamma(tmax, contact_gamma_pars, poscorr, initial_state),
         mc.cores=4) -> out_poscov_contact_gamma

## negative corr
source("GEM_SIR_contact&gamma.cov.R")
mclapply(seeds,
         function(s) gillespie.SIR.cov_cgamma(tmax, contact_gamma_pars, negcorr, initial_state),
         mc.cores=4) -> out_negcov_contact_gamma


###### FORMAT OUTPUT ######
## no correlation
## infected
timeSeq <- 0:150
storeMatrix.cov.CG.I <- array(NA, dim=c(length(timeSeq),length(out_nocov_contact_gamma)))
for (j in 1:length(out_nocov_contact_gamma)) {
  o <- out_nocov_contact_gamma[[j]]
  storeMatrix.cov.CG.I[,j] <- o[,3] 
}

## susceptible
storeMatrix.cov.CG.S <- array(NA, dim=c(length(timeSeq),length(out_nocov_contact_gamma)))
for (j in 1:length(out_nocov_contact_gamma)) {
  o <- out_nocov_contact_gamma[[j]]
  storeMatrix.cov.CG.S[,j] <- o[,2] 
}

## recovered
storeMatrix.cov.CG.R <- array(NA, dim=c(length(timeSeq),length(out_nocov_contact_gamma)))
for (j in 1:length(out_nocov_contact_gamma)) {
  o <- out_nocov_contact_gamma[[j]]
  storeMatrix.cov.CG.R[,j] <- o[,4] 
}

## positive cov
## infected
timeSeq <- 0:150
storeMatrix.poscov.CG.I <- array(NA, dim=c(length(timeSeq),length(out_poscov_contact_gamma)))
for (j in 1:length(out_poscov_contact_gamma)) {
  o <- out_poscov_contact_gamma[[j]]
  storeMatrix.poscov.CG.I[,j] <- o[,3] 
}

## susceptible
storeMatrix.poscov.CG.S <- array(NA, dim=c(length(timeSeq),length(out_poscov_contact_gamma)))
for (j in 1:length(out_poscov_contact_gamma)) {
  o <- out_poscov_contact_gamma[[j]]
  storeMatrix.poscov.CG.S[,j] <- o[,2] 
}

## recovered
storeMatrix.poscov.CG.R <- array(NA, dim=c(length(timeSeq),length(out_poscov_contact_gamma)))
for (j in 1:length(out_poscov_contact_gamma)) {
  o <- out_poscov_contact_gamma[[j]]
  storeMatrix.poscov.CG.R[,j] <- o[,4]
}

## neg cov
## infected
timeSeq <- 0:150
storeMatrix.negcov.CG.I <- array(NA, dim=c(length(timeSeq),length(out_negcov_contact_gamma)))
for (j in 1:length(out_negcov_contact_gamma)) {
  o <- out_negcov_contact_gamma[[j]]
  storeMatrix.negcov.CG.I[,j] <- o[,3] 
}

## susceptible
storeMatrix.negcov.CG.S <- array(NA, dim=c(length(timeSeq),length(out_negcov_contact_gamma)))
for (j in 1:length(out_negcov_contact_gamma)) {
  o <- out_negcov_contact_gamma[[j]]
  storeMatrix.negcov.CG.S[,j] <- o[,2] 
}

## recovered
storeMatrix.negcov.CG.R <- array(NA, dim=c(length(timeSeq),length(out_negcov_contact_gamma)))
for (j in 1:length(out_negcov_contact_gamma)) {
  o <- out_negcov_contact_gamma[[j]]
  storeMatrix.negcov.CG.R[,j] <- o[,4] 
}

###### PLOT ######
par(mfrow=c(1,3))

## no cov
plot(0:150, apply(storeMatrix.cov.CG.I, 1, mean), col="red", lwd=1.75, type="l", ylim=c(0,250), ylab="N", xlab="Time", main="No Cov")
lines(0:150, apply(storeMatrix.cov.CG.S, 1, mean), col="blue", lwd=1.75, type="l", ylim=c(0,150), ylab="N", xlab="Time")
lines(0:150, apply(storeMatrix.cov.CG.R, 1, mean), col="green", lwd=1.75, type="l", ylim=c(0,150), ylab="N", xlab="Time")
#abline(h=63, lty=2)
legend("topright",legend=c("S","I","R"),fill=c("blue","red","green"), cex=0.25)

## pos cov
plot(0:150, apply(storeMatrix.poscov.CG.I, 1, mean), col="red", lwd=1.75, type="l", ylim=c(0,250), ylab="N", xlab="Time", main="Pos Cov")
lines(0:150, apply(storeMatrix.poscov.CG.S, 1, mean), col="blue", lwd=1.75, type="l", ylim=c(0,150), ylab="N", xlab="Time")
lines(0:150, apply(storeMatrix.poscov.CG.R, 1, mean), col="green", lwd=1.75, type="l", ylim=c(0,150), ylab="N", xlab="Time")
#abline(h=63, lty=2)
legend("topright",legend=c("S","I","R"),fill=c("blue","red","green"), cex=0.25)

## neg cov
plot(0:150, apply(storeMatrix.negcov.CG.I, 1, mean), col="red", lwd=1.75, type="l", ylim=c(0,250), ylab="N", xlab="Time", main="Neg Cov")
lines(0:150, apply(storeMatrix.negcov.CG.S, 1, mean), col="blue", lwd=1.75, type="l", ylim=c(0,150), ylab="N", xlab="Time")
lines(0:150, apply(storeMatrix.negcov.CG.R, 1, mean), col="green", lwd=1.75, type="l", ylim=c(0,150), ylab="N", xlab="Time")
#abline(h=63, lty=2)
legend("topright",legend=c("S","I","R"),fill=c("blue","red","green"), cex=0.25)

