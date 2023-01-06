#### SHEDDING AND ALPHA COVARIATION ANALYSIS

library(tidyverse)
library(parallel)

## set up
seeds <- floor(runif(20,1,1e5)) # set seeds
tmax <- 150

nocorr <- matrix(c(1,0,0,1), nrow=2, byrow=T)
negcorr <- matrix(c(1,-.5,-.5,1), nrow=2, byrow=T)
poscorr <- matrix(c(1,.5,.5,1), nrow=2, byrow=T)

alpha_shed_pars_og = c(c=.1, shed=.05, alpha=.1, gamma=.1, beta=.25, d=.1, 
                      b=2.5, bs=.01, varS=0.0099, sd_s=.1, sd_a=0.5) # R0=3.8

alpha_shed_pars = c(c=.1, shed=.05, alpha=.1, gamma=.1, sd_c=0.1, sd_shed=0.05, 
                      sd_alpha=0.5, sd_gamma=0.5, b=2.5, d=.1, bs=.01)


initial_state <- floor(c(S =unname(((alpha_shed_pars["b"]-alpha_shed_pars["d"])/alpha_shed_pars["bs"]))-5, I=5, R=0))

###### RUN MULTIPLE SIMULATIONS ######
## no variation
source("GEM_SIR_noVar.R")
mclapply(seeds,
         function(s) gillespie.SIR.noVar(tmax, analytical_parms_shed, initial_state),
         mc.cores=4) -> out_no_var

## no covariation
mclapply(seeds, 
         function(i) gillespie.SIR.cov_storage(tmax=100, 
                                               params=alpha_shed_pars, 
                                               corr=nocorr, 
                                               initial_state, 
                                               covParams=c('shed','alpha')),
         mc.cores=4) -> out_nocov_shed_alpha

## negative covariation
mclapply(seeds, 
         function(i) gillespie.SIR.cov_storage(tmax=100, 
                                               params=alpha_shed_pars, 
                                               corr=matrix(c(1,-.5,-.5,1), nrow=2, byrow=T), 
                                               x=c(S=235,I=5,R=0), 
                                               covParams=c('shed','alpha')),
         mc.cores=4) -> out_negcov_shed_alpha

## positive covariation
mclapply(seeds, 
         function(i) gillespie.SIR.cov_storage(tmax=100, 
                                               params=alpha_shed_pars, 
                                               corr=poscorr, 
                                               initial_state, 
                                               covParams=c('shed','alpha')),
         mc.cores=4) -> out_poscov_shed_alpha

### Format Output ###

### no covariation
# susceptible
timeSeq <- 0:100
SA_S_nocov <- array(NA, dim=c(length(timeSeq), length(out_nocov_shed_alpha)+1))
SA_S_nocov[,1] <- timeSeq
for (i in 1:length(SA_S_nocov)) SA_S_nocov[,i+1] <- out_nocov_shed_alpha[[i]][[1]]$S

# infected
SA_I_nocov <- array(NA, dim=c(length(timeSeq), length(out_nocov_shed_alpha)+1))
SA_I_nocov[,1] <- timeSeq
for (i in 1:length(SA_I_nocov)) SA_I_nocov[,i+1] <- out_nocov_shed_alpha[[i]][[1]]$I

# recovered
SA_R_nocov <- array(NA, dim=c(length(timeSeq), length(out_nocov_shed_alpha)+1))
SA_R_nocov[,1] <- timeSeq
for (i in 1:length(SA_R_nocov)) SA_R_nocov[,i+1] <- out_nocov_shed_alpha[[i]][[1]]$R

### negative covariation
# susceptible
timeSeq <- 0:100
SA_S_negcov <- array(NA, dim=c(length(timeSeq), length(out_negcov_shed_alpha)+1))
SA_S_negcov[,1] <- timeSeq
for (i in 1:length(SA_S_nocov)) SA_S_negcov[,i+1] <- out_negcov_shed_alpha[[i]][[1]]$S

# infected
SA_I_negcov <- array(NA, dim=c(length(timeSeq), length(out_negcov_shed_alpha)+1))
SA_I_negcov[,1] <- timeSeq
for (i in 1:length(SA_I_negcov)) SA_I_negcov[,i+1] <- out_negcov_shed_alpha[[i]][[1]]$I

# recovered
SA_R_negcov <- array(NA, dim=c(length(timeSeq), length(out_negcov_shed_alpha)+1))
SA_R_negcov[,1] <- timeSeq
for (i in 1:length(SA_R_negcov)) SA_R_negcov[,i+1] <- out_negcov_shed_alpha[[i]][[1]]$R

### positive covariation
# susceptible
timeSeq <- 0:100
SA_S_poscov <- array(NA, dim=c(length(timeSeq), length(out_poscov_shed_alpha)+1))
SA_S_poscov[,1] <- timeSeq
for (i in 1:length(SA_S_poscov)) SA_S_poscov[,i+1] <- out_poscov_shed_alpha[[i]][[1]]$S

# infected
SA_I_poscov <- array(NA, dim=c(length(timeSeq), length(out_poscov_shed_alpha)+1))
SA_I_poscov[,1] <- timeSeq
for (i in 1:length(SA_I_poscov)) SA_I_poscov[,i+1] <- out_poscov_shed_alpha[[i]][[1]]$I

# recovered
SA_R_poscov <- array(NA, dim=c(length(timeSeq), length(out_poscov_shed_alpha)+1))
SA_R_poscov[,1] <- timeSeq
for (i in 1:length(SA_R_poscov)) SA_R_poscov[,i+1] <- out_poscov_shed_alpha[[i]][[1]]$R

### Plot 
par(mfrow=c(1,3))
## No Cov
plot(0:100, apply(SA_S_nocov, 1, mean), col="blue", lwd=1.75, type="l", ylim=c(0,250), ylab="S", xlab="Time", main="No Cov")
lines(0:100, apply(SA_I_nocov, 1, mean), col="red", lwd=1.75, type="l", ylim=c(0,250), ylab="S", xlab="Time")
lines(0:100, apply(SA_R_nocov, 1, mean), col="green", lwd=1.75, type="l", ylim=c(0,250), ylab="S", xlab="Time")
lines(0:150, apply(storeMatrix.novar.S, 1, mean), col="black", lwd=1.75, type="l", ylim=c(0,250), ylab="S", xlab="Time")
abline(h=63, lty=2)

## Pos Cov
plot(0:100, apply(SA_S_poscov, 1, mean), col="blue", lwd=1.75, type="l", ylim=c(0,250), ylab="S", xlab="Time", main="Pos Cov")
lines(0:100, apply(SA_I_poscov, 1, mean), col="red", lwd=1.75, type="l", ylim=c(0,250), ylab="S", xlab="Time")
lines(0:100, apply(SA_R_poscov, 1, mean), col="green", lwd=1.75, type="l", ylim=c(0,250), ylab="S", xlab="Time")
lines(0:150, apply(storeMatrix.novar.S, 1, mean), col="black", lwd=1.75, type="l", ylim=c(0,250), ylab="S", xlab="Time")
abline(h=63, lty=2)

## Neg Cov
plot(0:100, apply(SA_S_negcov, 1, mean), col="blue", lwd=1.75, type="l", ylim=c(0,250), ylab="S", xlab="Time", main="Neg Cov")
lines(0:100, apply(SA_I_negcov, 1, mean), col="red", lwd=1.75, type="l", ylim=c(0,250), ylab="S", xlab="Time")
lines(0:100, apply(SA_R_negcov, 1, mean), col="green", lwd=1.75, type="l", ylim=c(0,250), ylab="S", xlab="Time")
lines(0:150, apply(storeMatrix.novar.S, 1, mean), col="black", lwd=1.75, type="l", ylim=c(0,250), ylab="S", xlab="Time")
abline(h=63, lty=2)





############################################## ORIGINAL GEM CODE ########################################

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
par(mfrow=c(1,3), mar=c(4,4,1,1), oma=c(1,1,2,1))

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
abline(h=63, lty=2)
legend("topright",legend=c("S","I","R"),fill=c("blue","red","green"), cex=0.25)

## neg cov
plot(0:150, apply(storeMatrix.negcov.SA.I, 1, mean), col="red", lwd=1.75, type="l", ylim=c(0,250), ylab="N", xlab="Time", main="Neg Cov")
lines(0:150, apply(storeMatrix.negcov.SA.S, 1, mean), col="blue", lwd=1.75, type="l", ylim=c(0,150), ylab="N", xlab="Time")
lines(0:150, apply(storeMatrix.negcov.SA.R, 1, mean), col="green", lwd=1.75, type="l", ylim=c(0,150), ylab="N", xlab="Time")
abline(h=63, lty=2)
legend("topright",legend=c("S","I","R"),fill=c("blue","red","green"), cex=0.25)

## no cov - if you get NA
plot(0:150, apply(storeMatrix.cov.SA.I[,-which(apply(storeMatrix.cov.SA.I, 2, function(col) any(is.na(col))))], 1, mean), 
     col="red", lwd=1.75, type="l", ylim=c(0,250), ylab="N", xlab="Time", main="No Cov")
lines(0:150, apply(storeMatrix.cov.SA.S[,-which(apply(storeMatrix.cov.SA.S, 2, function(col) any(is.na(col))))], 1, mean), 
      col="blue", lwd=1.75, type="l", ylim=c(0,250), ylab="N", xlab="Time")
lines(0:150, apply(storeMatrix.cov.SA.R[,-which(apply(storeMatrix.cov.SA.R, 2, function(col) any(is.na(col))))], 1, mean), 
      col="green", lwd=1.75, type="l", ylim=c(0,250), ylab="N", xlab="Time")
legend("topright",legend=c("S","I","R"),fill=c("blue","red","green"), cex=0.25)
abline(h=63, lty=2)

## pos cov - if you get NA
plot(0:150, apply(storeMatrix.poscov.SA.I[,-which(apply(storeMatrix.poscov.SA.I, 2, function(col) any(is.na(col))))], 1, mean), 
     col="red", lwd=1.75, type="l", ylim=c(0,250), ylab="N", xlab="Time", main="Pos Cov")
lines(0:150, apply(storeMatrix.poscov.SA.S[,-which(apply(storeMatrix.poscov.SA.S, 2, function(col) any(is.na(col))))], 1, mean), 
      col="blue", lwd=1.75, type="l", ylim=c(0,250), ylab="N", xlab="Time")
lines(0:150, apply(storeMatrix.poscov.SA.R[,-which(apply(storeMatrix.poscov.SA.R, 2, function(col) any(is.na(col))))], 1, mean), 
      col="green", lwd=1.75, type="l", ylim=c(0,250), ylab="N", xlab="Time")
legend("topright",legend=c("S","I","R"),fill=c("blue","red","green"), cex=0.25)
abline(h=63, lty=2)

# plot code if you get NA
## neg cov - if you get NA
plot(0:150, apply(storeMatrix.negcov.SA.I[,-which(apply(storeMatrix.negcov.SA.I, 2, function(col) any(is.na(col))))], 1, mean), 
     col="red", lwd=1.75, type="l", ylim=c(0,250), ylab="N", xlab="Time", main="Neg Cov")
lines(0:150, apply(storeMatrix.negcov.SA.S[,-which(apply(storeMatrix.negcov.SA.S, 2, function(col) any(is.na(col))))], 1, mean), 
      col="blue", lwd=1.75, type="l", ylim=c(0,250), ylab="N", xlab="Time")
lines(0:150, apply(storeMatrix.negcov.SA.R[,-which(apply(storeMatrix.negcov.SA.R, 2, function(col) any(is.na(col))))], 1, mean), 
      col="green", lwd=1.75, type="l", ylim=c(0,250), ylab="N", xlab="Time")
legend("topright",legend=c("S","I","R"),fill=c("blue","red","green"), cex=0.25)
abline(h=63, lty=2)

mtext(side=3, outer=TRUE, "Covariance between shedding and virulence")

### all in one graphs
par(mfrow=c(1,3))
## no cov
plot(0:150, apply(storeMatrix.cov.SA.S[,-which(apply(storeMatrix.cov.SA.S, 2, function(col) any(is.na(col))))], 1, mean), 
     col="blue", lwd=1.75, type="l", ylim=c(0,200), ylab="S", xlab="Time", main="No Cov")
lines(0:150, apply(storeMatrix.cont.varA.S, 1, mean), col=alpha("black",.4), lwd=1.75, type="l", ylim=c(0,150), ylab="N", xlab="Time")
lines(0:150, apply(storeMatrix.cont.varS.S[,-which(apply(storeMatrix.cont.varS.S, 2, function(col) any(is.na(col))))], 1, mean), 
      col=alpha("black",.4), lwd=1.75, type="l", ylim=c(0,250), ylab="N", xlab="Time")
lines(0:150, apply(storeMatrix.no.var.S, 1, mean), col="black", lwd=1.75, type="l", ylim=c(0,150), ylab="N", xlab="Time")
abline(h=63, lty=2)
#legend("topright",legend=c("S","I","R"),fill=c("blue","red","green"), cex=0.25)

## pos cov
plot(0:150, apply(storeMatrix.poscov.SA.S[,-which(apply(storeMatrix.poscov.SA.S, 2, function(col) any(is.na(col))))], 1, mean), 
     col="blue", lwd=1.75, type="l", ylim=c(0,200), ylab="S", xlab="Time", main="Pos Cov")
lines(0:150, apply(storeMatrix.cont.varA.S, 1, mean), col=alpha("black",.4), lwd=1.75, type="l", ylim=c(0,150), ylab="N", xlab="Time")
lines(0:150, apply(storeMatrix.cont.varS.S[,-which(apply(storeMatrix.cont.varS.S, 2, function(col) any(is.na(col))))], 1, mean), 
      col=alpha("black",.4), lwd=1.75, type="l", ylim=c(0,250), ylab="N", xlab="Time")
lines(0:150, apply(storeMatrix.no.var.S, 1, mean), col="black", lwd=1.75, type="l", ylim=c(0,150), ylab="N", xlab="Time")
abline(h=63, lty=2)
#legend("topright",legend=c("S","I","R"),fill=c("blue","red","green"), cex=0.25)

## neg cov
plot(0:150, apply(storeMatrix.negcov.SA.S, 1, mean), col="blue", lwd=1.75, type="l", ylim=c(0,200), ylab="S", xlab="Time", main="Neg Cov")
lines(0:150, apply(storeMatrix.cont.varA.S, 1, mean), col=alpha("black",.4), lwd=1.75, type="l", ylim=c(0,200), ylab="N", xlab="Time")
lines(0:150, apply(storeMatrix.cont.varS.S[,-which(apply(storeMatrix.cont.varS.S, 2, function(col) any(is.na(col))))], 1, mean), 
      col=alpha("black",.4), lwd=1.75, type="l", ylim=c(0,250), ylab="N", xlab="Time")
lines(0:150, apply(storeMatrix.no.var.S, 1, mean), col="black", lwd=1.75, type="l", ylim=c(0,150), ylab="N", xlab="Time")
abline(h=63, lty=2)

