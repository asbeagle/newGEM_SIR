#### ALPHA AND GAMMA COVARIATION ANALYSIS

library(tidyverse)
library(parallel)

## set up
seeds <- floor(runif(20,1,1e5)) # set seeds
tmax <- 150

nocorr <- matrix(c(1,0,0,1), nrow=2, byrow=T)
negcorr <- matrix(c(1,-.5,-.5,1), nrow=2, byrow=T)
poscorr <- matrix(c(1,.5,.5,1), nrow=2, byrow=T)

alpha_gamma_pars = c(c=.1, shed=.05, alpha=.1, gamma=.1, beta=.25, d=.1, 
                    b=2.5, bs=.01, sd_g=.5, sd_a=.5) # R0=3.8


initial_state <- floor(c(S =unname(((alpha_gamma_pars["b"]-alpha_gamma_pars["d"])/alpha_gamma_pars["bs"]))-5, I=5, R=0))

###### RUN MULTIPLE SIMULATIONS ######
## no variation
source("GEM_SIR_noVar.R")
mclapply(seeds,
         function(s) gillespie.SIR.noVar(tmax, analytical_parms_shed, initial_state),
         mc.cores=4) -> out_no_var

## no corr
source("GEM_SIR_alpha&gamma.cov.R")
mclapply(seeds,
         function(s) gillespie.SIR.cov_alphagamma(tmax, alpha_gamma_pars, nocorr, initial_state),
         mc.cores=4) -> out_nocov_alpha_gamma

## positive corr
source("GEM_SIR_alpha&gamma.cov.R")
mclapply(seeds,
         function(s) gillespie.SIR.cov_alphagamma(tmax, alpha_gamma_pars, poscorr, initial_state),
         mc.cores=4) -> out_poscov_alpha_gamma

## negative corr
source("GEM_SIR_alpha&gamma.cov.R")
mclapply(seeds,
         function(s) gillespie.SIR.cov_alphagamma(tmax, alpha_gamma_pars, negcorr, initial_state),
         mc.cores=4) -> out_negcov_alpha_gamma

### running new model clay made
source("GEM_SIR_cov_storage.R")
mclapply(seeds,
         function(s) gillespie.SIR.cov_storage(tmax=100, 
                          params=params, 
                          corr=matrix(c(1,0,0,1), nrow=2, byrow=T), 
                          x=c(S=140,I=10,R=0), 
                          covParams=c('alpha','gamma'))) -> alpha_gamma_cov_out




###### FORMAT OUTPUT ######
## no correlation
## infected
timeSeq <- 0:10094
storeMatrix.cov.AG.I <- array(NA, dim=c(length(timeSeq),length(alpha_gamma_cov_out)))
for (j in 1:length(alpha_gamma_cov_out)) {
  for (i in 1:length(alpha_gamma_cov_out)){
  o <- alpha_gamma_cov_out[[j]] # pull out one seed run
  p <- o[[i]][3] # pull out I from each individual seed
  storeMatrix.cov.AG.I[,j] <- p[,1] 
  }
}



## susceptible
storeMatrix.cov.AG.S <- array(NA, dim=c(length(timeSeq),length(out_nocov_alpha_gamma)))
for (j in 1:length(out_nocov_alpha_gamma)) {
  o <- out_nocov_alpha_gamma[[j]]
  storeMatrix.cov.AG.S[,j] <- o[,2]
}

## recovered
storeMatrix.cov.AG.R <- array(NA, dim=c(length(timeSeq),length(out_nocov_alpha_gamma)))
for (j in 1:length(out_nocov_alpha_gamma)) {
  o <- out_nocov_alpha_gamma[[j]]
  storeMatrix.cov.AG.R[,j] <- o[,4] 
}

## positive cov
## infected
timeSeq <- 0:150
storeMatrix.poscov.AG.I <- array(NA, dim=c(length(timeSeq),length(out_poscov_alpha_gamma)))
for (j in 1:length(out_poscov_alpha_gamma)) {
  o <- out_poscov_alpha_gamma[[j]]
  storeMatrix.poscov.AG.I[,j] <- o[,3] 
}

## susceptible
storeMatrix.poscov.AG.S <- array(NA, dim=c(length(timeSeq),length(out_poscov_alpha_gamma)))
for (j in 1:length(out_poscov_alpha_gamma)) {
  o <- out_poscov_alpha_gamma[[j]]
  storeMatrix.poscov.AG.S[,j] <- o[,2] 
}

## recovered
storeMatrix.poscov.AG.R <- array(NA, dim=c(length(timeSeq),length(out_poscov_alpha_gamma)))
for (j in 1:length(out_poscov_alpha_gamma)) {
  o <- out_poscov_alpha_gamma[[j]]
  storeMatrix.poscov.AG.R[,j] <- o[,4]
}

## neg cov
## infected
timeSeq <- 0:150
storeMatrix.negcov.AG.I <- array(NA, dim=c(length(timeSeq),length(out_negcov_alpha_gamma)))
for (j in 1:length(out_negcov_alpha_gamma)) {
  o <- out_negcov_alpha_gamma[[j]]
  storeMatrix.negcov.AG.I[,j] <- o[,3] 
}

## susceptible
storeMatrix.negcov.AG.S <- array(NA, dim=c(length(timeSeq),length(out_negcov_alpha_gamma)))
for (j in 1:length(out_negcov_alpha_gamma)) {
  o <- out_negcov_alpha_gamma[[j]]
  storeMatrix.negcov.AG.S[,j] <- o[,2] 
}

## recovered
storeMatrix.negcov.AG.R <- array(NA, dim=c(length(timeSeq),length(out_negcov_alpha_gamma)))
for (j in 1:length(out_negcov_alpha_gamma)) {
  o <- out_negcov_alpha_gamma[[j]]
  storeMatrix.negcov.AG.R[,j] <- o[,4] 
}

###### PLOT ######
par(mfrow=c(1,3), mar=c(4,4,1,1), oma=c(1,1,2,1))

## no cov
plot(0:150, apply(storeMatrix.cov.AG.S, 1, mean), col="blue", lwd=1.75, type="l", ylim=c(0,200), ylab="S", xlab="Time", main="No Cov")
#lines(0:150, apply(storeMatrix.cov.AG.S, 1, mean), col="blue", lwd=1.75, type="l", ylim=c(0,150), ylab="N", xlab="Time")
#lines(0:150, apply(storeMatrix.cov.AG.R, 1, mean), col="green", lwd=1.75, type="l", ylim=c(0,150), ylab="N", xlab="Time")
lines(0:150, apply(storeMatrix.cont.varG.S, 1, mean), col=alpha("black",.4), lwd=1.75, type="l", ylim=c(0,150), ylab="N", xlab="Time")
lines(0:150, apply(storeMatrix.cont.varA.S, 1, mean), col=alpha("black",.4), lwd=1.75, type="l", ylim=c(0,150), ylab="N", xlab="Time")
lines(0:150, apply(storeMatrix.no.var.S, 1, mean), col="black", lwd=1.75, type="l", ylim=c(0,150), ylab="N", xlab="Time")
abline(h=63, lty=2)
#legend("topright",legend=c("S","I","R"),fill=c("blue","red","green"), cex=0.25)

## pos cov
plot(0:150, apply(storeMatrix.poscov.AG.S, 1, mean), col="blue", lwd=1.75, type="l", ylim=c(0,200), ylab="S", xlab="Time", main="Pos Cov")
#lines(0:150, apply(storeMatrix.poscov.AG.S, 1, mean), col="blue", lwd=1.75, type="l", ylim=c(0,150), ylab="N", xlab="Time")
#lines(0:150, apply(storeMatrix.poscov.AG.R, 1, mean), col="green", lwd=1.75, type="l", ylim=c(0,150), ylab="N", xlab="Time")
lines(0:150, apply(storeMatrix.cont.varG.S, 1, mean), col=alpha("black",.4), lwd=1.75, type="l", ylim=c(0,150), ylab="N", xlab="Time")
lines(0:150, apply(storeMatrix.cont.varA.S, 1, mean), col=alpha("black",.4), lwd=1.75, type="l", ylim=c(0,150), ylab="N", xlab="Time")
lines(0:150, apply(storeMatrix.no.var.S, 1, mean), col="black", lwd=1.75, type="l", ylim=c(0,150), ylab="N", xlab="Time")
abline(h=63, lty=2)
#legend("topright",legend=c("S","I","R"),fill=c("blue","red","green"), cex=0.25)

## neg cov
plot(0:150, apply(storeMatrix.negcov.AG.S, 1, mean), col="blue", lwd=1.75, type="l", ylim=c(0,200), ylab="S", xlab="Time", main="Neg Cov")
#lines(0:150, apply(storeMatrix.negcov.AG.S, 1, mean), col="blue", lwd=1.75, type="l",lty=3, ylim=c(0,150), ylab="N", xlab="Time")
#lines(0:150, apply(storeMatrix.negcov.AG.R, 1, mean), col="green", lwd=1.75, type="l", lty=4,ylim=c(0,150), ylab="N", xlab="Time")
lines(0:150, apply(storeMatrix.cont.varG.S, 1, mean), col=alpha("black",.4), lwd=1.75, type="l", ylim=c(0,150), ylab="N", xlab="Time")
lines(0:150, apply(storeMatrix.cont.varA.S, 1, mean), col=alpha("black",.4), lwd=1.75, type="l", ylim=c(0,150), ylab="N", xlab="Time")
lines(0:150, apply(storeMatrix.no.var.S, 1, mean), col="black", lwd=1.75, type="l", ylim=c(0,150), ylab="N", xlab="Time")
abline(h=63, lty=2)
legend("topright",legend=c("No Var","Var","Cov"),fill=c("black","grey","blue"), cex=.95)

mtext(side=3, outer=TRUE, "Covariance between recovery and virulence")


### all in one graphs
plot(0:150, apply(storeMatrix.cov.AG.S, 1, mean), col="green", lwd=1.75, type="l", ylim=c(0,250), ylab="N", xlab="Time")
lines(0:150, apply(storeMatrix.poscov.AG.S, 1, mean), col="blue", lwd=1.75, type="l", ylim=c(0,150), ylab="N", xlab="Time")
lines(0:150, apply(storeMatrix.negcov.AG.S, 1, mean), col="lightblue", lwd=1.75, type="l", ylim=c(0,150), ylab="N", xlab="Time")
lines(0:150, apply(storeMatrix.no.var.S, 1, mean), col="black", lwd=1.75, type="l", ylim=c(0,250), ylab="S", xlab="Time")
lines(0:150, apply(storeMatrix.cont.varG.S, 1, mean), col="red", lwd=1.75, type="l", ylim=c(0,150), ylab="N", xlab="Time")
abline(h=63, lty=2)

