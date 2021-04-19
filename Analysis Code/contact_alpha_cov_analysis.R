#### CONTACT AND ALPHA COVARIATION ANALYSIS CODE

library(tidyverse)
library(parallel)

## set up
seeds <- floor(runif(20,1,1e5)) # set seeds
tmax <- 150

nocorr <- matrix(c(1,0,0,1), nrow=2, byrow=T)
negcorr <- matrix(c(1,-.5,-.5,1), nrow=2, byrow=T)
poscorr <- matrix(c(1,.5,.5,1), nrow=2, byrow=T)


contact_alpha_pars = c(c=.1, shed=.05, alpha=.1, gamma=.1, beta=.25, d=.1, 
                       b=2.5, bs=.01, sd_c=1e-6, sd_a=1e-6) # R0=3.8


initial_state <- floor(c(S =unname(((contact_alpha_pars["b"]-contact_alpha_pars["d"])/contact_alpha_pars["bs"]))-5, I=5, R=0))

###### RUN MULTIPLE SIMULATIONS ######
## no variation
source("GEM_SIR_noVar.R")
mclapply(seeds,
         function(s) gillespie.SIR.noVar(tmax, contact_alpha_pars, initial_state),
         mc.cores=4) -> out_no_var

## no corr
source("GEM_SIR_contact&alpha.cov.R")
mclapply(seeds,
         function(s) gillespie.SIR.cov_calpha(tmax, contact_alpha_pars, nocorr, initial_state),
         mc.cores=4) -> out_nocov_contact_alpha

## positive corr
source("GEM_SIR_contact&alpha.cov.R")
mclapply(seeds,
         function(s) gillespie.SIR.cov_calpha(tmax, contact_alpha_pars, poscorr, initial_state),
         mc.cores=4) -> out_poscov_contact_alpha

## negative corr
source("GEM_SIR_contact&alpha.cov.R")
mclapply(seeds,
         function(s) gillespie.SIR.cov_calpha(tmax, contact_alpha_pars, negcorr, initial_state),
         mc.cores=4) -> out_negcov_contact_alpha

###### FORMAT OUTPUT ######
## no correlation
## infected
timeSeq <- 0:150
storeMatrix.cov.CA.I <- array(NA, dim=c(length(timeSeq),length(out_nocov_contact_alpha)))
for (j in 1:length(out_nocov_contact_alpha)) {
  o <- out_nocov_contact_alpha[[j]]
  storeMatrix.cov.CA.I[,j] <- o[,3] 
}

## susceptible
storeMatrix.cov.CA.S <- array(NA, dim=c(length(timeSeq),length(out_nocov_contact_alpha)))
for (j in 1:length(out_nocov_contact_alpha)) {
  o <- out_nocov_contact_alpha[[j]]
  storeMatrix.cov.CA.S[,j] <- o[,2] 
}

## recovered
storeMatrix.cov.CA.R <- array(NA, dim=c(length(timeSeq),length(out_nocov_contact_alpha)))
for (j in 1:length(out_nocov_contact_alpha)) {
  o <- out_nocov_contact_alpha[[j]]
  storeMatrix.cov.CA.R[,j] <- o[,4] 
}

## positive cov
## infected
timeSeq <- 0:150
storeMatrix.poscov.CA.I <- array(NA, dim=c(length(timeSeq),length(out_poscov_contact_alpha)))
for (j in 1:length(out_poscov_contact_alpha)) {
  o <- out_poscov_contact_alpha[[j]]
  storeMatrix.poscov.CA.I[,j] <- o[,3] 
}

## susceptible
storeMatrix.poscov.CA.S <- array(NA, dim=c(length(timeSeq),length(out_poscov_contact_alpha)))
for (j in 1:length(out_poscov_contact_alpha)) {
  o <- out_poscov_contact_alpha[[j]]
  storeMatrix.poscov.CA.S[,j] <- o[,2] 
}

## recovered
storeMatrix.poscov.CA.R <- array(NA, dim=c(length(timeSeq),length(out_poscov_contact_alpha)))
for (j in 1:length(out_poscov_contact_alpha)) {
  o <- out_poscov_contact_alpha[[j]]
  storeMatrix.poscov.CA.R[,j] <- o[,4]
}

## neg cov
## infected
timeSeq <- 0:150
storeMatrix.negcov.CA.I <- array(NA, dim=c(length(timeSeq),length(out_negcov_contact_alpha)))
for (j in 1:length(out_negcov_contact_alpha)) {
  o <- out_negcov_contact_alpha[[j]]
  storeMatrix.negcov.CA.I[,j] <- o[,3] 
}

## susceptible
storeMatrix.negcov.CA.S <- array(NA, dim=c(length(timeSeq),length(out_negcov_contact_alpha)))
for (j in 1:length(out_negcov_contact_alpha)) {
  o <- out_negcov_contact_alpha[[j]]
  storeMatrix.negcov.CA.S[,j] <- o[,2] 
}

## recovered
storeMatrix.negcov.CA.R <- array(NA, dim=c(length(timeSeq),length(out_negcov_contact_alpha)))
for (j in 1:length(out_negcov_contact_alpha)) {
  o <- out_negcov_contact_alpha[[j]]
  storeMatrix.negcov.CA.R[,j] <- o[,4] 
}

###### PLOT ######
par(mfrow=c(1,3))

## no cov
plot(0:150, apply(storeMatrix.cov.CA.I, 1, mean), col="red", lwd=1.75, type="l", ylim=c(0,250), ylab="N", xlab="Time", main="No Cov")
lines(0:150, apply(storeMatrix.cov.CA.S, 1, mean), col="blue", lwd=1.75, type="l", ylim=c(0,150), ylab="N", xlab="Time")
lines(0:150, apply(storeMatrix.cov.CA.R, 1, mean), col="green", lwd=1.75, type="l", ylim=c(0,150), ylab="N", xlab="Time")
#abline(h=63, lty=2)
legend("topright",legend=c("S","I","R"),fill=c("blue","red","green"), cex=0.25)

## pos cov
plot(0:150, apply(storeMatrix.poscov.CA.I, 1, mean), col="red", lwd=1.75, type="l", ylim=c(0,250), ylab="N", xlab="Time", main="Pos Cov")
lines(0:150, apply(storeMatrix.poscov.CA.S, 1, mean), col="blue", lwd=1.75, type="l", ylim=c(0,150), ylab="N", xlab="Time")
lines(0:150, apply(storeMatrix.poscov.CA.R, 1, mean), col="green", lwd=1.75, type="l", ylim=c(0,150), ylab="N", xlab="Time")
#abline(h=63, lty=2)
legend("topright",legend=c("S","I","R"),fill=c("blue","red","green"), cex=0.25)

## neg cov
plot(0:150, apply(storeMatrix.negcov.CA.I, 1, mean), col="red", lwd=1.75, type="l", ylim=c(0,250), ylab="N", xlab="Time", main="Neg Cov")
lines(0:150, apply(storeMatrix.negcov.CA.S, 1, mean), col="blue", lwd=1.75, type="l", ylim=c(0,150), ylab="N", xlab="Time")
lines(0:150, apply(storeMatrix.negcov.CA.R, 1, mean), col="green", lwd=1.75, type="l", ylim=c(0,150), ylab="N", xlab="Time")
#abline(h=63, lty=2)
legend("topright",legend=c("S","I","R"),fill=c("blue","red","green"), cex=0.25)



