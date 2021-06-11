#### CONTACT AND ALPHA COVARIATION ANALYSIS CODE

library(tidyverse)
library(parallel)

## set up
seeds <- floor(runif(20,1,1e5)) # set seeds
tmax <- 150

nocorr <- matrix(c(1,0,0,1), nrow=2, byrow=T)
negcorr <- matrix(c(1,-.5,-.5,1), nrow=2, byrow=T)
poscorr <- matrix(c(1,.5,.5,1), nrow=2, byrow=T)


contact_alpha_pars_og = c(c=.1, shed=.05, alpha=.1, gamma=.1, beta=.25, d=.1, 
                       b=2.5, bs=.01, sd_c=.5, sd_a=.5) # R0=3.8

contact_alpha_pars = c(c=.1, shed=.05, alpha=.1, gamma=.1, sd_c=0.5, sd_shed=0.025, 
                     sd_alpha=0.5, sd_gamma=0.5, b=2.5, d=.1, bs=.01)



initial_state <- floor(c(S =unname(((contact_alpha_pars["b"]-contact_alpha_pars["d"])/contact_alpha_pars["bs"]))-5, I=5, R=0))

###### RUN MULTIPLE SIMULATIONS ######
## no variation
source("GEM_SIR_noVar.R")
mclapply(seeds,
         function(s) gillespie.SIR.noVar(tmax, contact_alpha_pars, initial_state),
         mc.cores=4) -> out_no_var

## no covariation
source("GEM_SIR_cov_storage.R")
mclapply(seeds, 
         function(i) gillespie.SIR.cov_storage(tmax=20, 
                                               params=contact_alpha_pars, 
                                               corr=nocorr, 
                                               initial_state, 
                                               covParams=c('alpha','c'),
                                               seed=i),
         mc.cores=4) -> example

## negative covariation
mclapply(seeds, 
         function(i) gillespie.SIR.cov_storage(tmax=100, 
                                               params=contact_alpha_pars, 
                                               corr=negcorr, 
                                               x=c(S=235,I=5,R=0), 
                                               covParams=c('alpha','c')),
         mc.cores=4) -> out_negcov_alpha_contact

## positive covariation
mclapply(seeds, 
         function(i) gillespie.SIR.cov_storage(tmax=100, 
                                               params=contact_alpha_pars, 
                                               corr=poscorr, 
                                               initial_state, 
                                               covParams=c('alpha','c')),
         mc.cores=4) -> out_poscov_alpha_contact

### Format Output ###

### no covariation
# susceptible
timeSeq <- 0:100
alpha_contact_S_nocov <- array(NA, dim=c(length(timeSeq), length(out_nocov_alpha_contact)+1))
alpha_contact_S_nocov[,1] <- timeSeq
for (i in 1:length(alpha_contact_S_nocov)) alpha_contact_S_nocov[,i+1] <- out_nocov_alpha_contact[[i]][[1]]$S

# infected
alpha_contact_I_nocov <- array(NA, dim=c(length(timeSeq), length(out_nocov_alpha_contact)+1))
alpha_contact_I_nocov[,1] <- timeSeq
for (i in 1:length(alpha_contact_I_nocov)) alpha_contact_I_nocov[,i+1] <- out_nocov_alpha_contact[[i]][[1]]$I

# recovered
alpha_contact_R_nocov <- array(NA, dim=c(length(timeSeq), length(out_nocov_alpha_contact)+1))
alpha_contact_R_nocov[,1] <- timeSeq
for (i in 1:length(alpha_contact_R_nocov)) alpha_contact_R_nocov[,i+1] <- out_nocov_alpha_contact[[i]][[1]]$R

### negative covariation
# susceptible
timeSeq <- 0:100
alpha_contact_S_negcov <- array(NA, dim=c(length(timeSeq), length(out_negcov_alpha_contact)+1))
alpha_contact_S_negcov[,1] <- timeSeq
for (i in 1:length(alpha_contact_S_negcov)) alpha_contact_S_negcov[,i+1] <- out_negcov_alpha_contact[[i]][[1]]$S

# infected
alpha_contact_I_negcov <- array(NA, dim=c(length(timeSeq), length(out_negcov_alpha_contact)+1))
alpha_contact_I_negcov[,1] <- timeSeq
for (i in 1:length(alpha_contact_I_negcov)) alpha_contact_I_negcov[,i+1] <- out_negcov_alpha_contact[[i]][[1]]$I

# recovered
alpha_contact_R_negcov <- array(NA, dim=c(length(timeSeq), length(out_negcov_alpha_contact)+1))
alpha_contact_R_negcov[,1] <- timeSeq
for (i in 1:length(alpha_contact_R_negcov)) alpha_contact_R_negcov[,i+1] <- out_negcov_alpha_contact[[i]][[1]]$R

### positive covariation
# susceptible
timeSeq <- 0:100
alpha_contact_S_poscov <- array(NA, dim=c(length(timeSeq), length(out_poscov_alpha_contact)+1))
alpha_contact_S_poscov[,1] <- timeSeq
for (i in 1:length(alpha_contact_S_poscov)) alpha_contact_S_poscov[,i+1] <- out_poscov_alpha_contact[[i]][[1]]$S

# infected
alpha_contact_I_poscov <- array(NA, dim=c(length(timeSeq), length(out_poscov_alpha_contact)+1))
alpha_contact_I_poscov[,1] <- timeSeq
for (i in 1:length(alpha_contact_I_poscov)) alpha_contact_I_poscov[,i+1] <- out_poscov_alpha_contact[[i]][[1]]$I

# recovered
alpha_contact_R_poscov <- array(NA, dim=c(length(timeSeq), length(out_poscov_alpha_contact)+1))
alpha_contact_R_poscov[,1] <- timeSeq
for (i in 1:length(alpha_contact_R_poscov)) alpha_contact_R_poscov[,i+1] <- out_poscov_alpha_contact[[i]][[1]]$R

### Plot 
par(mfrow=c(1,3))
## No Cov
plot(0:100, apply(alpha_contact_S_nocov, 1, mean), col="blue", lwd=1.75, type="l", ylim=c(0,250), ylab="S", xlab="Time", main="No Cov")
lines(0:100, apply(alpha_contact_I_nocov, 1, mean), col="red", lwd=1.75, type="l", ylim=c(0,250), ylab="S", xlab="Time")
lines(0:100, apply(alpha_contact_R_nocov, 1, mean), col="green", lwd=1.75, type="l", ylim=c(0,250), ylab="S", xlab="Time")

## Neg Cov
plot(0:100, apply(alpha_contact_S_negcov, 1, mean), col="blue", lwd=1.75, type="l", ylim=c(0,250), ylab="S", xlab="Time", main="Neg Cov")
lines(0:100, apply(alpha_contact_I_negcov, 1, mean), col="red", lwd=1.75, type="l", ylim=c(0,250), ylab="S", xlab="Time")
lines(0:100, apply(alpha_contact_R_negcov, 1, mean), col="green", lwd=1.75, type="l", ylim=c(0,250), ylab="S", xlab="Time")

## Pos Cov
plot(0:100, apply(alpha_contact_S_poscov, 1, mean), col="blue", lwd=1.75, type="l", ylim=c(0,250), ylab="S", xlab="Time", main="Pos Cov")
lines(0:100, apply(alpha_contact_I_poscov, 1, mean), col="red", lwd=1.75, type="l", ylim=c(0,250), ylab="S", xlab="Time")
lines(0:100, apply(alpha_contact_R_poscov, 1, mean), col="green", lwd=1.75, type="l", ylim=c(0,250), ylab="S", xlab="Time")





################################### ORIGINAL GEM CODE ##################################################
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
par(mfrow=c(1,3), mar=c(4,4,1,1), oma=c(1,1,2,1))

## no cov
plot(0:150, apply(storeMatrix.cov.CA.I, 1, mean), col="red", lwd=1.75, type="l", ylim=c(0,250), ylab="N", xlab="Time", main="No Cov")
lines(0:150, apply(storeMatrix.cov.CA.S, 1, mean), col="blue", lwd=1.75, type="l", ylim=c(0,150), ylab="N", xlab="Time")
lines(0:150, apply(storeMatrix.cov.CA.R, 1, mean), col="green", lwd=1.75, type="l", ylim=c(0,150), ylab="N", xlab="Time")
#abline(h=80, lty=2)
legend("topright",legend=c("S","I","R"),fill=c("blue","red","green"), cex=0.25)

## no cov - if you get NA
plot(0:150, apply(storeMatrix.cov.CA.I[,-which(apply(storeMatrix.cov.CA.I, 2, function(col) any(is.na(col))))], 1, mean), 
     col="red", lwd=1.75, type="l", ylim=c(0,250), ylab="N", xlab="Time", main="No Cov")
lines(0:150, apply(storeMatrix.cov.CA.S[,-which(apply(storeMatrix.cov.CA.S, 2, function(col) any(is.na(col))))], 1, mean), 
      col="blue", lwd=1.75, type="l", ylim=c(0,250), ylab="N", xlab="Time")
lines(0:150, apply(storeMatrix.cov.CA.R[,-which(apply(storeMatrix.cov.CA.R, 2, function(col) any(is.na(col))))], 1, mean), 
      col="green", lwd=1.75, type="l", ylim=c(0,250), ylab="N", xlab="Time")
legend("topright",legend=c("S","I","R"),fill=c("blue","red","green"), cex=0.25)
abline(h=63, lty=2)

## pos cov - if you get NA
plot(0:150, apply(storeMatrix.poscov.CA.I[,-which(apply(storeMatrix.poscov.CA.I, 2, function(col) any(is.na(col))))], 1, mean), 
     col="red", lwd=1.75, type="l", ylim=c(0,250), ylab="N", xlab="Time", main="Pos Cov")
lines(0:150, apply(storeMatrix.poscov.CA.S[,-which(apply(storeMatrix.poscov.CA.S, 2, function(col) any(is.na(col))))], 1, mean), 
      col="blue", lwd=1.75, type="l", ylim=c(0,250), ylab="N", xlab="Time")
lines(0:150, apply(storeMatrix.poscov.CA.R[,-which(apply(storeMatrix.poscov.CA.R, 2, function(col) any(is.na(col))))], 1, mean), 
      col="green", lwd=1.75, type="l", ylim=c(0,250), ylab="N", xlab="Time")
legend("topright",legend=c("S","I","R"),fill=c("blue","red","green"), cex=0.25)
abline(h=63, lty=2)

## neg cov - if you get NA
plot(0:150, apply(storeMatrix.negcov.CA.I[,-which(apply(storeMatrix.negcov.CA.I, 2, function(col) any(is.na(col))))], 1, mean), 
     col="red", lwd=1.75, type="l", ylim=c(0,250), ylab="N", xlab="Time", main="Neg Cov")
lines(0:150, apply(storeMatrix.negcov.CA.S[,-which(apply(storeMatrix.negcov.CA.S, 2, function(col) any(is.na(col))))], 1, mean), 
      col="blue", lwd=1.75, type="l", ylim=c(0,250), ylab="N", xlab="Time")
lines(0:150, apply(storeMatrix.negcov.CA.R[,-which(apply(storeMatrix.negcov.CA.R, 2, function(col) any(is.na(col))))], 1, mean), 
      col="green", lwd=1.75, type="l", ylim=c(0,250), ylab="N", xlab="Time")
legend("topright",legend=c("S","I","R"),fill=c("blue","red","green"), cex=0.25)
abline(h=63, lty=2)

mtext(side=3, outer=TRUE, "Covariance between contact rates and virulence")

## pos cov
plot(0:150, apply(storeMatrix.poscov.CA.I, 1, mean), col="red", lwd=1.75, type="l", ylim=c(0,250), ylab="N", xlab="Time", main="Pos Cov")
lines(0:150, apply(storeMatrix.poscov.CA.S, 1, mean), col="blue", lwd=1.75, type="l", ylim=c(0,150), ylab="N", xlab="Time")
lines(0:150, apply(storeMatrix.poscov.CA.R, 1, mean), col="green", lwd=1.75, type="l", ylim=c(0,150), ylab="N", xlab="Time")
abline(h=80, lty=2)
legend("topright",legend=c("S","I","R"),fill=c("blue","red","green"), cex=0.25)

## neg cov
plot(0:150, apply(storeMatrix.negcov.CA.I, 1, mean), col="red", lwd=1.75, type="l", ylim=c(0,250), ylab="N", xlab="Time", main="Neg Cov")
lines(0:150, apply(storeMatrix.negcov.CA.S, 1, mean), col="blue", lwd=1.75, type="l", ylim=c(0,150), ylab="N", xlab="Time")
lines(0:150, apply(storeMatrix.negcov.CA.R, 1, mean), col="green", lwd=1.75, type="l", ylim=c(0,150), ylab="N", xlab="Time")
abline(h=80, lty=2)
legend("topright",legend=c("S","I","R"),fill=c("blue","red","green"), cex=0.25)

mtext(side=3, outer=TRUE, "Covariance between recovery and virulence")


### all in one graphs
par(mfrow=c(1,3))
## no cov
plot(0:150, apply(storeMatrix.cov.CA.S[,-which(apply(storeMatrix.cov.CA.S, 2, function(col) any(is.na(col))))], 1, mean), 
     col="blue", lwd=1.75, type="l", ylim=c(0,200), ylab="S", xlab="Time", main="No Cov")
lines(0:150, apply(storeMatrix.cont.varA.S, 1, mean), col=alpha("black",.4), lwd=1.75, type="l", ylim=c(0,150), ylab="N", xlab="Time")
lines(0:150, apply(storeMatrix.cont.varB.S, 1, mean), col=alpha("black",.4), lwd=1.75, type="l", ylim=c(0,150), ylab="N", xlab="Time")
lines(0:150, apply(storeMatrix.no.var.S, 1, mean), col="black", lwd=1.75, type="l", ylim=c(0,150), ylab="N", xlab="Time")
abline(h=63, lty=2)
#legend("topright",legend=c("S","I","R"),fill=c("blue","red","green"), cex=0.25)

## pos cov
plot(0:150, apply(storeMatrix.poscov.CA.S[,-which(apply(storeMatrix.poscov.CA.S, 2, function(col) any(is.na(col))))], 1, mean), 
     col="blue", lwd=1.75, type="l", ylim=c(0,200), ylab="S", xlab="Time", main="Pos Cov")
lines(0:150, apply(storeMatrix.cont.varA.S, 1, mean), col=alpha("black",.4), lwd=1.75, type="l", ylim=c(0,150), ylab="N", xlab="Time")
lines(0:150, apply(storeMatrix.cont.varB.S, 1, mean), col=alpha("black",.4), lwd=1.75, type="l", ylim=c(0,150), ylab="N", xlab="Time")
lines(0:150, apply(storeMatrix.no.var.S, 1, mean), col="black", lwd=1.75, type="l", ylim=c(0,150), ylab="N", xlab="Time")
abline(h=63, lty=2)
#legend("topright",legend=c("S","I","R"),fill=c("blue","red","green"), cex=0.25)

## neg cov
plot(0:150, apply(storeMatrix.negcov.CA.S[,-which(apply(storeMatrix.negcov.CA.S, 2, function(col) any(is.na(col))))], 1, mean), 
     col="blue", lwd=1.75, type="l", ylim=c(0,200), ylab="S", xlab="Time", main="Neg Cov")
lines(0:150, apply(storeMatrix.cont.varA.S, 1, mean), col=alpha("black",.4), lwd=1.75, type="l", ylim=c(0,150), ylab="N", xlab="Time")
lines(0:150, apply(storeMatrix.cont.varB.S, 1, mean), col=alpha("black",.4), lwd=1.75, type="l", ylim=c(0,150), ylab="N", xlab="Time")
lines(0:150, apply(storeMatrix.no.var.S, 1, mean), col="black", lwd=1.75, type="l", ylim=c(0,150), ylab="N", xlab="Time")
abline(h=63, lty=2)

