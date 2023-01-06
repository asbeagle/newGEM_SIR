#### SHEDDING AND CONTACT COVARIATION ANALYSIS CODE

library(tidyverse)
library(parallel)

## set up
seeds <- floor(runif(20,1,1e5)) # set seeds
tmax <- 150

nocorr <- matrix(c(1,0,0,1), nrow=2, byrow=T)
negcorr <- matrix(c(1,-.5,-.5,1), nrow=2, byrow=T)
poscorr <- matrix(c(1,.5,.5,1), nrow=2, byrow=T)


contact_shed_pars_og = c(c=.1, shed=.05, alpha=.1, gamma=.1, beta=.25, d=.1, 
                      b=2.5, bs=.01, sd_s=.1, sd_c=.1) # R0=3.8

contact_shed_pars = c(c=.1, shed=.05, alpha=.1, gamma=.1, sd_c=0.1, sd_shed=0.1, 
                       sd_alpha=0.5, sd_gamma=0.5, b=2.5, d=.1, bs=.01)


initial_state <- floor(c(S =unname(((contact_shed_pars["b"]-contact_shed_pars["d"])/contact_shed_pars["bs"]))-5, I=5, R=0))

###### RUN MULTIPLE SIMULATIONS ######
## no variation
source("GEM_SIR_noVar.R")
mclapply(seeds,
         function(s) gillespie.SIR.noVar(tmax, analytical_parms_shed, initial_state),
         mc.cores=4) -> out_no_var

## no covariation
mclapply(seeds, 
         function(i) gillespie.SIR.cov_storage(tmax=100, 
                                               params=contact_shed_pars, 
                                               corr=nocorr, 
                                               initial_state, 
                                               covParams=c('shed','c')),
         mc.cores=4) -> out_nocov_shed_contact

## negative covariation
mclapply(seeds, 
         function(i) gillespie.SIR.cov_storage(tmax=100, 
                                               params=contact_shed_pars, 
                                               corr=negcorr, 
                                               x=c(S=235,I=5,R=0), 
                                               covParams=c('shed','c')),
         mc.cores=4) -> out_negcov_shed_contact

## positive covariation
mclapply(seeds, 
         function(i) gillespie.SIR.cov_storage(tmax=100, 
                                               params=contact_shed_pars, 
                                               corr=poscorr, 
                                               initial_state, 
                                               covParams=c('shed','c')),
         mc.cores=4) -> out_poscov_shed_contact

### Format Output ###

### no covariation
# susceptible
timeSeq <- 0:100
SC_S_nocov <- array(NA, dim=c(length(timeSeq), length(out_nocov_shed_contact)+1))
SC_S_nocov[,1] <- timeSeq
for (i in 1:length(SC_S_nocov)) SC_S_nocov[,i+1] <- out_nocov_shed_contact[[i]][[1]]$S

# infected
SC_I_nocov <- array(NA, dim=c(length(timeSeq), length(out_nocov_shed_contact)+1))
SC_I_nocov[,1] <- timeSeq
for (i in 1:length(SC_I_nocov)) SC_I_nocov[,i+1] <- out_nocov_shed_contact[[i]][[1]]$I

# recovered
SC_R_nocov <- array(NA, dim=c(length(timeSeq), length(out_nocov_shed_contact)+1))
SC_R_nocov[,1] <- timeSeq
for (i in 1:length(SC_R_nocov)) SC_R_nocov[,i+1] <- out_nocov_shed_contact[[i]][[1]]$R

### negative covariation
# susceptible
timeSeq <- 0:100
SC_S_nocov <- array(NA, dim=c(length(timeSeq), length(out_negcov_shed_contact)+1))
SC_S_nocov[,1] <- timeSeq
for (i in 1:length(SC_S_nocov)) SC_S_nocov[,i+1] <- out_negcov_shed_contact[[i]][[1]]$S

# infected
SC_I_negcov <- array(NA, dim=c(length(timeSeq), length(out_negcov_shed_contact)+1))
SC_I_negcov[,1] <- timeSeq
for (i in 1:length(SC_I_negcov)) SC_I_negcov[,i+1] <- out_negcov_shed_contact[[i]][[1]]$I

# recovered
SC_R_negcov <- array(NA, dim=c(length(timeSeq), length(out_negcov_shed_contact)+1))
SC_R_negcov[,1] <- timeSeq
for (i in 1:length(SC_R_negcov)) SC_R_negcov[,i+1] <- out_negcov_shed_contact[[i]][[1]]$R

### positive covariation
# susceptible
timeSeq <- 0:100
SC_S_poscov <- array(NA, dim=c(length(timeSeq), length(out_poscov_shed_contact)+1))
SC_S_poscov[,1] <- timeSeq
for (i in 1:length(SC_S_poscov)) SC_S_poscov[,i+1] <- out_poscov_shed_contact[[i]][[1]]$S

# infected
SC_I_poscov <- array(NA, dim=c(length(timeSeq), length(out_poscov_shed_contact)+1))
SC_I_poscov[,1] <- timeSeq
for (i in 1:length(SC_I_poscov)) SC_I_poscov[,i+1] <- out_poscov_shed_contact[[i]][[1]]$I

# recovered
SC_R_poscov <- array(NA, dim=c(length(timeSeq), length(out_poscov_shed_contact)+1))
SC_R_poscov[,1] <- timeSeq
for (i in 1:length(SC_R_poscov)) SC_R_poscov[,i+1] <- out_poscov_shed_contact[[i]][[1]]$R

### Plot 
par(mfrow=c(1,3))
## No Cov
plot(0:100, apply(SC_S_nocov, 1, mean), col="blue", lwd=1.75, type="l", ylim=c(0,250), ylab="S", xlab="Time", main="No Cov")
lines(0:100, apply(SC_I_nocov, 1, mean), col="red", lwd=1.75, type="l", ylim=c(0,250), ylab="S", xlab="Time")
lines(0:100, apply(SC_R_nocov, 1, mean), col="green", lwd=1.75, type="l", ylim=c(0,250), ylab="S", xlab="Time")

## Neg Cov
plot(0:100, apply(SC_S_negcov, 1, mean), col="blue", lwd=1.75, type="l", ylim=c(0,250), ylab="S", xlab="Time", main="Neg Cov")
lines(0:100, apply(SC_I_negcov, 1, mean), col="red", lwd=1.75, type="l", ylim=c(0,250), ylab="S", xlab="Time")
lines(0:100, apply(SC_R_negcov, 1, mean), col="green", lwd=1.75, type="l", ylim=c(0,250), ylab="S", xlab="Time")

## Pos Cov
plot(0:100, apply(SC_S_poscov, 1, mean), col="blue", lwd=1.75, type="l", ylim=c(0,250), ylab="S", xlab="Time", main="Pos Cov")
lines(0:100, apply(SC_I_poscov, 1, mean), col="red", lwd=1.75, type="l", ylim=c(0,250), ylab="S", xlab="Time")
lines(0:100, apply(SC_R_poscov, 1, mean), col="green", lwd=1.75, type="l", ylim=c(0,250), ylab="S", xlab="Time")






###################################### ORIGINAL GEM CODE ################################################
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
par(mfrow=c(1,3), mar=c(4,4,1,1), oma=c(1,1,5,1))

## no cov
plot(0:150, apply(storeMatrix.cov.CS.I, 1, mean), col="red", lwd=1.75, type="l", ylim=c(0,250), ylab="N", xlab="Time", main="No Cov")
lines(0:150, apply(storeMatrix.cov.CS.S, 1, mean), col="blue", lwd=1.75, type="l", ylim=c(0,150), ylab="N", xlab="Time")
lines(0:150, apply(storeMatrix.cov.CS.R, 1, mean), col="green", lwd=1.75, type="l", ylim=c(0,150), ylab="N", xlab="Time")
abline(h=63, lty=2)
legend("topright",legend=c("S","I","R"),fill=c("blue","red","green"), cex=0.25)

## pos cov
plot(0:150, apply(storeMatrix.poscov.CS.I, 1, mean), col="red", lwd=1.75, type="l", ylim=c(0,250), ylab="N", xlab="Time", main="Pos Cov")
lines(0:150, apply(storeMatrix.poscov.CS.S, 1, mean), col="blue", lwd=1.75, type="l", ylim=c(0,150), ylab="N", xlab="Time")
lines(0:150, apply(storeMatrix.poscov.CS.R, 1, mean), col="green", lwd=1.75, type="l", ylim=c(0,150), ylab="N", xlab="Time")
abline(h=63, lty=2)
legend("topright",legend=c("S","I","R"),fill=c("blue","red","green"), cex=0.25)



## neg cov
plot(0:150, apply(storeMatrix.negcov.CS.I, 1, mean), col="red", lwd=1.75, type="l", ylim=c(0,250), ylab="N", xlab="Time", main="Neg Cov")
lines(0:150, apply(storeMatrix.negcov.CS.S, 1, mean), col="blue", lwd=1.75, type="l", ylim=c(0,150), ylab="N", xlab="Time")
lines(0:150, apply(storeMatrix.negcov.CS.R, 1, mean), col="green", lwd=1.75, type="l", ylim=c(0,150), ylab="N", xlab="Time")
#abline(h=63, lty=2)
legend("topright",legend=c("S","I","R"),fill=c("blue","red","green"), cex=0.25)


mtext(side=3, outer=TRUE, "Covariance between shedding and contact")
par(mfrow=c(1,3), mar=c(4,4,1,1), oma=c(1,1,2,1))

## no cov - if you get NA
plot(0:150, apply(storeMatrix.cov.CS.I[,-which(apply(storeMatrix.cov.CS.I, 2, function(col) any(is.na(col))))], 1, mean), 
     col="red", lwd=1.75, type="l", ylim=c(0,250), ylab="N", xlab="Time", main="No Cov")
lines(0:150, apply(storeMatrix.cov.CS.S[,-which(apply(storeMatrix.cov.CS.S, 2, function(col) any(is.na(col))))], 1, mean), 
      col="blue", lwd=1.75, type="l", ylim=c(0,250), ylab="N", xlab="Time")
lines(0:150, apply(storeMatrix.cov.CS.R[,-which(apply(storeMatrix.cov.CS.R, 2, function(col) any(is.na(col))))], 1, mean), 
      col="green", lwd=1.75, type="l", ylim=c(0,250), ylab="N", xlab="Time")
legend("topright",legend=c("S","I","R"),fill=c("blue","red","green"), cex=0.25)
abline(h=63, lty=2)

## pos cov - if you get NA
plot(0:150, apply(storeMatrix.poscov.CS.I[,-which(apply(storeMatrix.poscov.CS.I, 2, function(col) any(is.na(col))))], 1, mean), 
     col="red", lwd=1.75, type="l", ylim=c(0,250), ylab="N", xlab="Time", main="Pos Cov")
lines(0:150, apply(storeMatrix.poscov.CS.S[,-which(apply(storeMatrix.poscov.CS.S, 2, function(col) any(is.na(col))))], 1, mean), 
      col="blue", lwd=1.75, type="l", ylim=c(0,250), ylab="N", xlab="Time")
lines(0:150, apply(storeMatrix.poscov.CS.R[,-which(apply(storeMatrix.poscov.CS.R, 2, function(col) any(is.na(col))))], 1, mean), 
      col="green", lwd=1.75, type="l", ylim=c(0,250), ylab="N", xlab="Time")
legend("topright",legend=c("S","I","R"),fill=c("blue","red","green"), cex=0.25)
abline(h=63, lty=2)

## neg cov - if you get NA
plot(0:150, apply(storeMatrix.negcov.CS.I[,-which(apply(storeMatrix.negcov.CS.I, 2, function(col) any(is.na(col))))], 1, mean), 
     col="red", lwd=1.75, type="l", ylim=c(0,250), ylab="N", xlab="Time", main="Neg Cov")
lines(0:150, apply(storeMatrix.negcov.CS.S[,-which(apply(storeMatrix.negcov.CS.S, 2, function(col) any(is.na(col))))], 1, mean), 
      col="blue", lwd=1.75, type="l", ylim=c(0,250), ylab="N", xlab="Time")
lines(0:150, apply(storeMatrix.negcov.CS.R[,-which(apply(storeMatrix.negcov.CS.R, 2, function(col) any(is.na(col))))], 1, mean), 
      col="green", lwd=1.75, type="l", ylim=c(0,250), ylab="N", xlab="Time")
legend("topright",legend=c("S","I","R"),fill=c("blue","red","green"), cex=0.25)
abline(h=63, lty=2)

mtext(side=3, outer=TRUE, "Covariance between shedding and contact")

### all in one graphs

## no cov
plot(0:150, apply(storeMatrix.cov.CS.S[,-which(apply(storeMatrix.cov.CS.S, 2, function(col) any(is.na(col))))], 1, mean), 
     col="blue", lwd=1.75, type="l", ylim=c(0,200), ylab="S", xlab="Time", main="No Cov")
lines(0:150, apply(storeMatrix.cont.varB.S, 1, mean), col=alpha("black",.4), lwd=1.75, type="l", ylim=c(0,150), ylab="N", xlab="Time")
lines(0:150, apply(storeMatrix.cont.varS.S[,-which(apply(storeMatrix.cont.varS.S, 2, function(col) any(is.na(col))))], 1, mean), 
      col=alpha("black",.4), lwd=1.75, type="l", ylim=c(0,250), ylab="N", xlab="Time")
lines(0:150, apply(storeMatrix.no.var.S, 1, mean), col="black", lwd=1.75, type="l", ylim=c(0,150), ylab="N", xlab="Time")
abline(h=63, lty=2)
#legend("topright",legend=c("S","I","R"),fill=c("blue","red","green"), cex=0.25)

## pos cov
plot(0:150, apply(storeMatrix.poscov.CS.S[,-which(apply(storeMatrix.poscov.CS.S, 2, function(col) any(is.na(col))))], 1, mean), 
     col="blue", lwd=1.75, type="l", ylim=c(0,200), ylab="S", xlab="Time", main="Pos Cov")
lines(0:150, apply(storeMatrix.cont.varS.S[,-which(apply(storeMatrix.cont.varS.S, 2, function(col) any(is.na(col))))], 1, mean), 
      col=alpha("black",.4), lwd=1.75, type="l", ylim=c(0,250), ylab="N", xlab="Time")
lines(0:150, apply(storeMatrix.cont.varB.S, 1, mean), col=alpha("black",.4), lwd=1.75, type="l", ylim=c(0,150), ylab="N", xlab="Time")
lines(0:150, apply(storeMatrix.no.var.S, 1, mean), col="black", lwd=1.75, type="l", ylim=c(0,150), ylab="N", xlab="Time")
abline(h=63, lty=2)
#legend("topright",legend=c("S","I","R"),fill=c("blue","red","green"), cex=0.25)

## neg cov
plot(0:150, apply(storeMatrix.negcov.CS.S[,-which(apply(storeMatrix.negcov.CS.S, 2, function(col) any(is.na(col))))], 1, mean), 
     col="blue", lwd=1.75, type="l", ylim=c(0,200), ylab="S", xlab="Time", main="Neg Cov")
lines(0:150, apply(storeMatrix.cont.varB.S, 1, mean), col=alpha("black",.4), lwd=1.75, type="l", ylim=c(0,150), ylab="N", xlab="Time")
lines(0:150, apply(storeMatrix.cont.varS.S[,-which(apply(storeMatrix.cont.varS.S, 2, function(col) any(is.na(col))))], 1, mean), 
      col=alpha("black",.4), lwd=1.75, type="l", ylim=c(0,250), ylab="N", xlab="Time")
lines(0:150, apply(storeMatrix.no.var.S, 1, mean), col="black", lwd=1.75, type="l", ylim=c(0,150), ylab="N", xlab="Time")
abline(h=63, lty=2)
