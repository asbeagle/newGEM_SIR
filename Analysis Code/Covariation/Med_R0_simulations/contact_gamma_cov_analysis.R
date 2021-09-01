#### CONTACT AND GAMMA COVARIATION ANALYSIS CODE
library(tidyverse)
library(parallel)

## set up
seeds <- floor(runif(20,1,1e5)) # set seeds
tmax <- 100

nocorr <- matrix(c(1,0,0,1), nrow=2, byrow=T)
negcorr <- matrix(c(1,-.5,-.5,1), nrow=2, byrow=T)
poscorr <- matrix(c(1,.5,.5,1), nrow=2, byrow=T)


contact_gamma_pars_og = c(c=.1, shed=.05, alpha=.1, gamma=.1, beta=.25, d=.1, 
                      b=2.5, bs=.01, sd_c=.5, sd_g=.5) # R0=3.8

contact_gamma_pars = c(c=.1, shed=.05, alpha=.1, gamma=.1, sd_c=0.05, sd_shed=0.025, 
                       sd_alpha=0.5, sd_gamma=0.5, b=2.5, d=.1, bs=.01)


initial_state <- floor(c(S =unname(((contact_gamma_pars["b"]-contact_gamma_pars["d"])/contact_gamma_pars["bs"]))-5, I=5, R=0))

###### RUN MULTIPLE SIMULATIONS ######
## no variation
source("GEM_SIR_noVar.R")
mclapply(seeds,
         function(s) gillespie.SIR.noVar(tmax, contact_gamma_pars, initial_state),
         mc.cores=4) -> out_no_var

## no covariation
mclapply(seeds, 
         function(i) gillespie.SIR.cov_storage(tmax=100, 
                                               params=contact_gamma_pars, 
                                               corr=nocorr, 
                                               initial_state, 
                                               covParams=c('c','gamma')),
         mc.cores=4) -> out_nocov_gamma_contact

## negative covariation
mclapply(seeds, 
         function(i) gillespie.SIR.cov_storage(tmax=100, 
                                               params=contact_gamma_pars, 
                                               corr=negcorr, 
                                               x=c(S=235,I=5,R=0), 
                                               covParams=c('c','gamma')),
         mc.cores=4) -> out_negcov_gamma_contact

## positive covariation
mclapply(seeds, 
         function(i) gillespie.SIR.cov_storage(tmax=100, 
                                               params=contact_gamma_pars, 
                                               corr=poscorr, 
                                               initial_state, 
                                               covParams=c('c','gamma')),
         mc.cores=4) -> out_poscov_gamma_contact

### Format Output ###

### no covariation
# susceptible
timeSeq <- 0:100
contact_gamma_S_nocov <- array(NA, dim=c(length(timeSeq), length(out_nocov_contact_gamma)+1))
contact_gamma_S_nocov[,1] <- timeSeq
for (i in 1:length(contact_gamma_S_nocov)) contact_gamma_S_nocov[,i+1] <- out_nocov_contact_gamma[[i]][[1]]$S

# infected
contact_gamma_I_nocov <- array(NA, dim=c(length(timeSeq), length(out_nocov_contact_gamma)+1))
contact_gamma_I_nocov[,1] <- timeSeq
for (i in 1:length(contact_gamma_I_nocov)) contact_gamma_I_nocov[,i+1] <- out_nocov_contaact_gamma[[i]][[1]]$I

# recovered
contact_gamma_R_nocov <- array(NA, dim=c(length(timeSeq), length(out_nocov_contact_gamma)+1))
contact_gamma_R_nocov[,1] <- timeSeq
for (i in 1:length(contact_gamma_R_nocov)) contact_gamma_R_nocov[,i+1] <- out_nocov_contact_gamma[[i]][[1]]$R

### negative covariation
# susceptible
contact_gamma_S_negcov <- array(NA, dim=c(length(timeSeq), length(out_negcov_contact_gamma)+1))
contact_gamma_S_negcov[,1] <- timeSeq
for (i in 1:length(contact_gamma_S_negcov)) contact_gamma_S_negcov[,i+1] <- out_negcov_contact_gamma[[i]][[1]]$S

# infected
contact_gamma_I_negcov <- array(NA, dim=c(length(timeSeq), length(out_negcov_contact_gamma)+1))
contact_gamma_I_negcov[,1] <- timeSeq
for (i in 1:length(contact_gamma_I_negcov)) contact_gamma_I_negcov[,i+1] <- out_negcov_contact_gamma[[i]][[1]]$I

# recovered
contact_gamma_R_negcov <- array(NA, dim=c(length(timeSeq), length(out_negcov_contact_gamma)+1))
contact_gamma_R_negcov[,1] <- timeSeq
for (i in 1:length(contact_gamma_R_negcov)) contact_gamma_R_negcov[,i+1] <- out_negcov_contact_gamma[[i]][[1]]$R

### positive covariation
# susceptible
contact_gamma_S_poscov <- array(NA, dim=c(length(timeSeq), length(out_poscov_contact_gamma)+1))
contact_gamma_S_poscov[,1] <- timeSeq
for (i in 1:length(contact_gamma_S_poscov)) contact_gamma_S_poscov[,i+1] <- out_poscov_contact_gamma[[i]][[1]]$S

# infected
contact_gamma_I_poscov <- array(NA, dim=c(length(timeSeq), length(out_poscov_contact_gamma)+1))
contact_gamma_I_poscov[,1] <- timeSeq
for (i in 1:length(contact_gamma_I_poscov)) contact_gamma_I_poscov[,i+1] <- out_poscov_contact_gamma[[i]][[1]]$I

# recovered
contact_gamma_R_poscov <- array(NA, dim=c(length(timeSeq), length(out_poscov_contact_gamma)+1))
contact_gamma_R_poscov[,1] <- timeSeq
for (i in 1:length(contact_gamma_R_poscov)) contact_gamma_R_poscov[,i+1] <- out_poscov_contact_gamma[[i]][[1]]$R

### Plot 
par(mfrow=c(1,3))
## No Cov
plot(0:100, apply(contact_gamma_S, 1, mean), col="blue", lwd=1.75, type="l", ylim=c(0,250), ylab="S", xlab="Time", main="No Cov")
lines(0:100, apply(contact_gamma_I, 1, mean), col="red", lwd=1.75, type="l", ylim=c(0,250), ylab="S", xlab="Time")
lines(0:100, apply(contact_gamma_R, 1, mean), col="green", lwd=1.75, type="l", ylim=c(0,250), ylab="S", xlab="Time")

## Neg Cov
plot(0:100, apply(contact_gamma_S_negcov, 1, mean), col="blue", lwd=1.75, type="l", ylim=c(0,250), ylab="S", xlab="Time", main="Neg Cov")
lines(0:100, apply(contact_gamma_I_negcov, 1, mean), col="red", lwd=1.75, type="l", ylim=c(0,250), ylab="S", xlab="Time")
lines(0:100, apply(contact_gamma_R_negcov, 1, mean), col="green", lwd=1.75, type="l", ylim=c(0,250), ylab="S", xlab="Time")

## Pos Cov
plot(0:100, apply(contact_gamma_S_poscov, 1, mean), col="blue", lwd=1.75, type="l", ylim=c(0,250), ylab="S", xlab="Time", main="Pos Cov")
lines(0:100, apply(contact_gamma_I_poscov, 1, mean), col="red", lwd=1.75, type="l", ylim=c(0,250), ylab="S", xlab="Time")
lines(0:100, apply(contact_gamma_R_poscov, 1, mean), col="green", lwd=1.75, type="l", ylim=c(0,250), ylab="S", xlab="Time")




######################################### ORIGINAL GEM CODE #############################################

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
par(mfrow=c(1,3), mar=c(4,4,1,1), oma=c(1,1,2,1))

## no cov
plot(0:150, apply(storeMatrix.cov.CG.I, 1, mean), col="red", lwd=1.75, type="l", ylim=c(0,250), ylab="N", xlab="Time", main="No Cov")
lines(0:150, apply(storeMatrix.cov.CG.S, 1, mean), col="blue", lwd=1.75, type="l", ylim=c(0,150), ylab="N", xlab="Time")
lines(0:150, apply(storeMatrix.cov.CG.R, 1, mean), col="green", lwd=1.75, type="l", ylim=c(0,150), ylab="N", xlab="Time")
abline(h=63, lty=2)
legend("topright",legend=c("S","I","R"),fill=c("blue","red","green"), cex=0.25)

## pos cov - if you get NA
plot(0:150, apply(storeMatrix.poscov.CG.I[,-which(apply(storeMatrix.poscov.CG.I, 2, function(col) any(is.na(col))))], 1, mean), 
     col="red", lwd=1.75, type="l", ylim=c(0,250), ylab="N", xlab="Time", main="Pos Cov")
lines(0:150, apply(storeMatrix.poscov.CG.S[,-which(apply(storeMatrix.poscov.CG.S, 2, function(col) any(is.na(col))))], 1, mean), 
      col="blue", lwd=1.75, type="l", ylim=c(0,250), ylab="N", xlab="Time")
lines(0:150, apply(storeMatrix.poscov.CG.R[,-which(apply(storeMatrix.poscov.CG.R, 2, function(col) any(is.na(col))))], 1, mean), 
      col="green", lwd=1.75, type="l", ylim=c(0,250), ylab="N", xlab="Time")
legend("topright",legend=c("S","I","R"),fill=c("blue","red","green"), cex=0.25)
abline(h=63, lty=2)

## neg cov - if you get NA
plot(0:150, apply(storeMatrix.negcov.CG.I[,-which(apply(storeMatrix.negcov.CG.I, 2, function(col) any(is.na(col))))], 1, mean), 
     col="red", lwd=1.75, type="l", ylim=c(0,250), ylab="N", xlab="Time", main="Neg Cov")
lines(0:150, apply(storeMatrix.negcov.CG.S[,-which(apply(storeMatrix.negcov.CG.S, 2, function(col) any(is.na(col))))], 1, mean), 
      col="blue", lwd=1.75, type="l", ylim=c(0,250), ylab="N", xlab="Time")
lines(0:150, apply(storeMatrix.negcov.CG.R[,-which(apply(storeMatrix.negcov.CG.R, 2, function(col) any(is.na(col))))], 1, mean), 
      col="green", lwd=1.75, type="l", ylim=c(0,250), ylab="N", xlab="Time")
legend("topright",legend=c("S","I","R"),fill=c("blue","red","green"), cex=0.25)
abline(h=63, lty=2)

mtext(side=3, outer=TRUE, "Covariance between contact rates and recovery")


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

### all in one graphs
par(mfrow=c(1,3))
## no cov
plot(0:150, apply(storeMatrix.cov.CG.S, 1, mean), col="blue", lwd=1.75, type="l", ylim=c(0,200), ylab="S", xlab="Time", main="No Cov")
lines(0:150, apply(storeMatrix.cont.varG.S, 1, mean), col=alpha("black",.4), lwd=1.75, type="l", ylim=c(0,150), ylab="N", xlab="Time")
lines(0:150, apply(storeMatrix.cont.varB.S, 1, mean), col=alpha("black",.4), lwd=1.75, type="l", ylim=c(0,150), ylab="N", xlab="Time")
lines(0:150, apply(storeMatrix.no.var.S, 1, mean), col="black", lwd=1.75, type="l", ylim=c(0,150), ylab="N", xlab="Time")
abline(h=63, lty=2)
#legend("topright",legend=c("S","I","R"),fill=c("blue","red","green"), cex=0.25)

## pos cov
plot(0:150, apply(storeMatrix.poscov.CG.S[,-which(apply(storeMatrix.poscov.CG.S, 2, function(col) any(is.na(col))))], 1, mean), 
     col="blue", lwd=1.75, type="l", ylim=c(0,200), ylab="S", xlab="Time", main="Pos Cov")
lines(0:150, apply(storeMatrix.cont.varG.S, 1, mean), col=alpha("black",.4), lwd=1.75, type="l", ylim=c(0,150), ylab="N", xlab="Time")
lines(0:150, apply(storeMatrix.cont.varB.S, 1, mean), col=alpha("black",.4), lwd=1.75, type="l", ylim=c(0,150), ylab="N", xlab="Time")
lines(0:150, apply(storeMatrix.no.var.S, 1, mean), col="black", lwd=1.75, type="l", ylim=c(0,150), ylab="N", xlab="Time")
abline(h=63, lty=2)
#legend("topright",legend=c("S","I","R"),fill=c("blue","red","green"), cex=0.25)

## neg cov
plot(0:150, apply(storeMatrix.negcov.CG.S[,-which(apply(storeMatrix.negcov.CG.S, 2, function(col) any(is.na(col))))], 1, mean), 
     col="blue", lwd=1.75, type="l", ylim=c(0,200), ylab="S", xlab="Time", main="Neg Cov")
lines(0:150, apply(storeMatrix.cont.varG.S, 1, mean), col=alpha("black",.4), lwd=1.75, type="l", ylim=c(0,150), ylab="N", xlab="Time")
lines(0:150, apply(storeMatrix.cont.varB.S, 1, mean), col=alpha("black",.4), lwd=1.75, type="l", ylim=c(0,150), ylab="N", xlab="Time")
lines(0:150, apply(storeMatrix.no.var.S, 1, mean), col="black", lwd=1.75, type="l", ylim=c(0,150), ylab="N", xlab="Time")
abline(h=63, lty=2)


