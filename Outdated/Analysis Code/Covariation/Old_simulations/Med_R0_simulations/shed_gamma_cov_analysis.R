#### SHEDDING AND GAMMA COVARIATION ANALYSIS

library(tidyverse)
library(parallel)

## set up
seeds <- floor(runif(20,1,1e5)) # set seeds
tmax <- 150

nocorr <- matrix(c(1,0,0,1), nrow=2, byrow=T)
negcorr <- matrix(c(1,-.5,-.5,1), nrow=2, byrow=T)
poscorr <- matrix(c(1,.5,.5,1), nrow=2, byrow=T)

gamma_shed_pars_og = c(c=.1, shed=.05, alpha=.1, gamma=.1, beta=.25, d=.1, 
                    b=2.5, bs=.01, varS=0.0099, sd_s=.5, sd_g=.5) # R0=3.8

gamma_shed_pars = c(c=.1, shed=.05, alpha=.1, gamma=.1, sd_c=0.1, sd_shed=0.05, 
                    sd_alpha=0.5, sd_gamma=0.5, b=2.5, d=.1, bs=.01)


initial_state <- floor(c(S =unname(((gamma_shed_pars["b"]-gamma_shed_pars["d"])/gamma_shed_pars["bs"]))-5, I=5, R=0))

###### RUN MULTIPLE SIMULATIONS ######
## no variation
source("GEM_SIR_noVar.R")
mclapply(seeds,
         function(s) gillespie.SIR.noVar(tmax, analytical_parms_shed, initial_state),
         mc.cores=4) -> out_no_var

## no covariation
mclapply(seeds, 
         function(i) gillespie.SIR.cov_storage(tmax=100, 
                                               params=gamma_shed_pars, 
                                               corr=nocorr, 
                                               initial_state, 
                                               covParams=c('shed','gamma')),
         mc.cores=4) -> out_nocov_shed_gamma

## negative covariation
mclapply(seeds, 
         function(i) gillespie.SIR.cov_storage(tmax=100, 
                                               params=gamma_shed_pars, 
                                               corr=negcorr, 
                                               initial_state, 
                                               covParams=c('shed','gamma')),
         mc.cores=4) -> out_negcov_shed_gamma

## positive covariation
mclapply(seeds, 
         function(i) gillespie.SIR.cov_storage(tmax=100, 
                                               params=gamma_shed_pars, 
                                               corr=poscorr, 
                                               initial_state, 
                                               covParams=c('shed','gamma')),
         mc.cores=4) -> out_poscov_shed_gamma

### Format Output ###

### no covariation
# susceptible
timeSeq <- 0:100
SG_S_nocov <- array(NA, dim=c(length(timeSeq), length(out_nocov_shed_gamma)+1))
SG_S_nocov[,1] <- timeSeq
for (i in 1:length(SG_S_nocov)) SG_S_nocov[,i+1] <- out_nocov_shed_gamma[[i]][[1]]$S

# infected
SG_I_nocov <- array(NA, dim=c(length(timeSeq), length(out_nocov_shed_gamma)+1))
SG_I_nocov[,1] <- timeSeq
for (i in 1:length(SG_I_nocov)) SG_I_nocov[,i+1] <- out_nocov_shed_gamma[[i]][[1]]$I

# recovered
SG_R_nocov <- array(NA, dim=c(length(timeSeq), length(out_nocov_shed_gamma)+1))
SG_R_nocov[,1] <- timeSeq
for (i in 1:length(SG_R_nocov)) SG_R_nocov[,i+1] <- out_nocov_shed_gamma[[i]][[1]]$R

### negative covariation
# susceptible
timeSeq <- 0:100
SG_S_negcov <- array(NA, dim=c(length(timeSeq), length(out_negcov_shed_gamma)+1))
SG_S_negcov[,1] <- timeSeq
for (i in 1:length(SG_S_negcov)) SG_S_negcov[,i+1] <- out_negcov_shed_gamma[[i]][[1]]$S

# infected
SG_I_negcov <- array(NA, dim=c(length(timeSeq), length(out_negcov_shed_gamma)+1))
SG_I_negcov[,1] <- timeSeq
for (i in 1:length(SG_I_negcov)) SG_I_negcov[,i+1] <- out_negcov_shed_gamma[[i]][[1]]$I

# recovered
SG_R_negcov <- array(NA, dim=c(length(timeSeq), length(out_negcov_shed_gamma)+1))
SG_R_negcov[,1] <- timeSeq
for (i in 1:length(SG_R_negcov)) SG_R_negcov[,i+1] <- out_negcov_shed_gamma[[i]][[1]]$R

### positive covariation
# susceptible
timeSeq <- 0:100
SG_S_poscov <- array(NA, dim=c(length(timeSeq), length(out_poscov_shed_gamma)+1))
SG_S_poscov[,1] <- timeSeq
for (i in 1:length(SG_S_poscov)) SG_S_poscov[,i+1] <- out_poscov_shed_gamma[[i]][[1]]$S

# infected
SG_I_poscov <- array(NA, dim=c(length(timeSeq), length(out_poscov_shed_gamma)+1))
SG_I_poscov[,1] <- timeSeq
for (i in 1:length(SG_I_poscov)) SG_I_poscov[,i+1] <- out_poscov_shed_gamma[[i]][[1]]$I

# recovered
SG_R_poscov <- array(NA, dim=c(length(timeSeq), length(out_poscov_shed_gamma)+1))
SG_R_poscov[,1] <- timeSeq
for (i in 1:length(SG_R_poscov)) SG_R_poscov[,i+1] <- out_poscov_shed_gamma[[i]][[1]]$R

### Plot 
par(mfrow=c(1,3))
## No Cov
plot(0:100, apply(SG_S_nocov, 1, mean), col="blue", lwd=1.75, type="l", ylim=c(0,250), ylab="S", xlab="Time", main="No Cov")
lines(0:100, apply(SG_I_nocov, 1, mean), col="red", lwd=1.75, type="l", ylim=c(0,250), ylab="S", xlab="Time")
lines(0:100, apply(SG_R_nocov, 1, mean), col="green", lwd=1.75, type="l", ylim=c(0,250), ylab="S", xlab="Time")
lines(0:150, apply(storeMatrix.novar.S, 1, mean), col="black", lwd=1.75, type="l", ylim=c(0,250), ylab="S", xlab="Time")
abline(h=63, lty=2)

## Neg Cov
plot(0:100, apply(SG_S_negcov, 1, mean), col="blue", lwd=1.75, type="l", ylim=c(0,250), ylab="S", xlab="Time", main="Neg Cov")
lines(0:100, apply(SG_I_negcov, 1, mean), col="red", lwd=1.75, type="l", ylim=c(0,250), ylab="S", xlab="Time")
lines(0:100, apply(SG_R_negcov, 1, mean), col="green", lwd=1.75, type="l", ylim=c(0,250), ylab="S", xlab="Time")
lines(0:150, apply(storeMatrix.novar.S, 1, mean), col="black", lwd=1.75, type="l", ylim=c(0,250), ylab="S", xlab="Time")
abline(h=63, lty=2)

## Pos Cov
plot(0:100, apply(SG_S_poscov, 1, mean), col="blue", lwd=1.75, type="l", ylim=c(0,250), ylab="S", xlab="Time", main="Pos Cov")
lines(0:100, apply(SG_I_poscov, 1, mean), col="red", lwd=1.75, type="l", ylim=c(0,250), ylab="S", xlab="Time")
lines(0:100, apply(SG_R_poscov, 1, mean), col="green", lwd=1.75, type="l", ylim=c(0,250), ylab="S", xlab="Time")
lines(0:150, apply(storeMatrix.novar.S, 1, mean), col="black", lwd=1.75, type="l", ylim=c(0,250), ylab="S", xlab="Time")
abline(h=63, lty=2)


###################################### ORIGINAL GEM CODE ################################################

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
         function(s) gillespie.SIR.cov_shedgamma(tmax, gamma_shed_pars, negcorr, initial_state),
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
par(mfrow=c(1,3), mar=c(4,4,1,1), oma=c(1,1,2,1))

## no cov - if you get NA
plot(0:150, apply(storeMatrix.cov.SG.I[,-which(apply(storeMatrix.cov.SG.I, 2, function(col) any(is.na(col))))], 1, mean), 
     col="red", lwd=1.75, type="l", ylim=c(0,250), ylab="N", xlab="Time", main="No Cov")
lines(0:150, apply(storeMatrix.cov.SG.S[,-which(apply(storeMatrix.cov.SG.S, 2, function(col) any(is.na(col))))], 1, mean), 
      col="blue", lwd=1.75, type="l", ylim=c(0,250), ylab="N", xlab="Time")
lines(0:150, apply(storeMatrix.cov.SG.R[,-which(apply(storeMatrix.cov.SG.R, 2, function(col) any(is.na(col))))], 1, mean), 
      col="green", lwd=1.75, type="l", ylim=c(0,250), ylab="N", xlab="Time")
legend("topright",legend=c("S","I","R"),fill=c("blue","red","green"), cex=0.25)
abline(h=63, lty=2)

## pos cov - if you get NA
plot(0:150, apply(storeMatrix.poscov.SG.I[,-which(apply(storeMatrix.poscov.SG.I, 2, function(col) any(is.na(col))))], 1, mean), 
     col="red", lwd=1.75, type="l", ylim=c(0,250), ylab="N", xlab="Time", main="Pos Cov")
lines(0:150, apply(storeMatrix.poscov.SG.S[,-which(apply(storeMatrix.poscov.SG.S, 2, function(col) any(is.na(col))))], 1, mean), 
      col="blue", lwd=1.75, type="l", ylim=c(0,250), ylab="N", xlab="Time")
lines(0:150, apply(storeMatrix.poscov.SG.R[,-which(apply(storeMatrix.poscov.SG.R, 2, function(col) any(is.na(col))))], 1, mean), 
      col="green", lwd=1.75, type="l", ylim=c(0,250), ylab="N", xlab="Time")
legend("topright",legend=c("S","I","R"),fill=c("blue","red","green"), cex=0.25)
abline(h=63, lty=2)

# plot code if you get NA
## neg cov - if you get NA
plot(0:150, apply(storeMatrix.negcov.SG.I[,-which(apply(storeMatrix.negcov.SG.I, 2, function(col) any(is.na(col))))], 1, mean), 
     col="red", lwd=1.75, type="l", ylim=c(0,250), ylab="N", xlab="Time", main="Neg Cov")
lines(0:150, apply(storeMatrix.negcov.SG.S[,-which(apply(storeMatrix.negcov.SG.S, 2, function(col) any(is.na(col))))], 1, mean), 
      col="blue", lwd=1.75, type="l", ylim=c(0,250), ylab="N", xlab="Time")
lines(0:150, apply(storeMatrix.negcov.SG.R[,-which(apply(storeMatrix.negcov.SG.R, 2, function(col) any(is.na(col))))], 1, mean), 
      col="green", lwd=1.75, type="l", ylim=c(0,250), ylab="N", xlab="Time")
legend("topright",legend=c("S","I","R"),fill=c("blue","red","green"), cex=0.25)
abline(h=63, lty=2)

mtext(side=3, outer=TRUE, "Covariance between shedding and recovery")


## no cov
plot(0:150, apply(storeMatrix.cov.SG.I[,-which(apply(storeMatrix.cov.SG.I, 2, function(col) any(is.na(col))))], 1, mean), 
     col="blue", lwd=1.75, type="l", ylim=c(0,250), ylab="N", xlab="Time", main="No Cov")
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



par(mfrow=c(1,3), mar=c(4,4,1,1), oma=c(1,1,2,1))

## no cov
plot(0:150, apply(storeMatrix.cov.SG.S[,-which(apply(storeMatrix.cov.SG.S, 2, function(col) any(is.na(col))))], 1, mean), 
     col="blue", lwd=1.75, type="l", ylim=c(0,200), ylab="S", xlab="Time", main="No Cov")
lines(0:150, apply(storeMatrix.cont.varG.S, 1, mean), col=alpha("black",.3), lwd=1.75, type="l", ylim=c(0,150), ylab="N", xlab="Time")
lines(0:150, apply(storeMatrix.cont.varS.S[,-which(apply(storeMatrix.cont.varS.S, 2, function(col) any(is.na(col))))], 1, mean), 
      col=alpha("black",.4), lwd=1.75, type="l", ylim=c(0,250), ylab="N", xlab="Time")
lines(0:150, apply(storeMatrix.no.var.S, 1, mean), col="black", lwd=1.75, type="l", ylim=c(0,150), ylab="N", xlab="Time")
abline(h=63, lty=2)
#legend("topright",legend=c("S","I","R"),fill=c("blue","red","green"), cex=0.25)

## pos cov
plot(0:150, apply(storeMatrix.poscov.SG.S[,-which(apply(storeMatrix.poscov.SG.S, 2, function(col) any(is.na(col))))], 1, mean), 
     col="blue", lwd=1.75, type="l", ylim=c(0,200), ylab="S", xlab="Time", main="Pos Cov")
lines(0:150, apply(storeMatrix.cont.varS.S[,-which(apply(storeMatrix.cont.varS.S, 2, function(col) any(is.na(col))))], 1, mean), 
      col=alpha("black",.4), lwd=1.75, type="l", ylim=c(0,250), ylab="N", xlab="Time")
lines(0:150, apply(storeMatrix.cont.varG.S, 1, mean), col=alpha("black",.3), lwd=1.75, type="l", ylim=c(0,150), ylab="N", xlab="Time")
lines(0:150, apply(storeMatrix.no.var.S, 1, mean), col="black", lwd=1.75, type="l", ylim=c(0,150), ylab="N", xlab="Time")
abline(h=63, lty=2)
#legend("topright",legend=c("S","I","R"),fill=c("blue","red","green"), cex=0.25)

## neg cov
plot(0:150, apply(storeMatrix.negcov.SG.S[,-which(apply(storeMatrix.negcov.SG.S, 2, function(col) any(is.na(col))))], 1, mean), 
     col="blue", lwd=1.75, type="l", ylim=c(0,200), ylab="S", xlab="Time", main="Neg Cov")
lines(0:150, apply(storeMatrix.cont.varG.S, 1, mean), col=alpha("black",.3), lwd=1.75, type="l", ylim=c(0,150), ylab="N", xlab="Time")
lines(0:150, apply(storeMatrix.cont.varS.S[,-which(apply(storeMatrix.cont.varS.S, 2, function(col) any(is.na(col))))], 1, mean), 
      col=alpha("black",.4), lwd=1.75, type="l", ylim=c(0,250), ylab="N", xlab="Time")
lines(0:150, apply(storeMatrix.no.var.S, 1, mean), col="black", lwd=1.75, type="l", ylim=c(0,150), ylab="N", xlab="Time")
abline(h=63, lty=2)