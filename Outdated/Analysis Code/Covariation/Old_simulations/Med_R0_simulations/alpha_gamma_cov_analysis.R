#### ALPHA AND GAMMA COVARIATION ANALYSIS

library(tidyverse)
library(parallel)

## set up
seeds <- floor(runif(20,1,1e5)) # set seeds
tmax <- 150

nocorr <- matrix(c(1,0,0,1), nrow=2, byrow=T)
negcorr <- matrix(c(1,-.5,-.5,1), nrow=2, byrow=T)
poscorr <- matrix(c(1,.5,.5,1), nrow=2, byrow=T)

alpha_gamma_pars_og = c(c=.1, shed=.05, alpha=.1, gamma=.1, beta=.25, d=.1, 
                    b=2.5, bs=.01, sd_g=.5, sd_a=.5) # R0=3.8

alpha_gamma_pars = c(c=.1, shed=.05, alpha=.1, gamma=.1, sd_c=0.05, sd_shed=0.025, 
                    sd_alpha=0.5, sd_gamma=0.5, b=2.5, d=.1, bs=.01)



initial_state <- floor(c(S =unname(((alpha_gamma_pars["b"]-alpha_gamma_pars["d"])/alpha_gamma_pars["bs"]))-5, I=5, R=0))

### Run Multiple Simulations
## no variation
source("GEM_SIR_noVar.R")
mclapply(seeds,
         function(s) gillespie.SIR.noVar(tmax, alpha_gamma_pars, initial_state),
         mc.cores=4) -> out_no_var

## no covariation
mclapply(seeds, 
         function(i) gillespie.SIR.cov_storage(tmax=100, 
                                               params=alpha_gamma_pars, 
                                               corr=nocorr, 
                                               initial_state, 
                                               covParams=c('alpha','gamma')),
         mc.cores=4) -> out_nocov_alpha_gamma

## negative covariation
mclapply(seeds, 
         function(i) gillespie.SIR.cov_storage(tmax=100, 
                                               params=alpha_gamma_pars, 
                                               corr=matrix(c(1,-.5,-.5,1), nrow=2, byrow=T), 
                                               x=c(S=235,I=5,R=0), 
                                               covParams=c('alpha','gamma')),
         mc.cores=4) -> out_negcov_alpha_gamma

## positive covariation
mclapply(seeds, 
         function(i) gillespie.SIR.cov_storage(tmax=100, 
                                               params=alpha_gamma_pars, 
                                               corr=poscorr, 
                                               initial_state, 
                                               covParams=c('alpha','gamma')),
         mc.cores=4) -> out_poscov_alpha_gamma

### Format Output ###

### no covariation
# susceptible
timeSeq <- 0:100
alpha_gamma_S <- array(NA, dim=c(length(timeSeq), length(out_nocov_alpha_gamma)+1))
alpha_gamma_S[,1] <- timeSeq
for (i in 1:length(alpha_gamma_S)) alpha_gamma_S[,i+1] <- out_nocov_alpha_gamma[[i]][[1]]$S

# infected
alpha_gamma_I <- array(NA, dim=c(length(timeSeq), length(out_nocov_alpha_gamma)+1))
alpha_gamma_I[,1] <- timeSeq
for (i in 1:length(alpha_gamma_I)) alpha_gamma_I[,i+1] <- out_nocov_alpha_gamma[[i]][[1]]$I

# recovered
alpha_gamma_R <- array(NA, dim=c(length(timeSeq), length(out_nocov_alpha_gamma)+1))
alpha_gamma_R[,1] <- timeSeq
for (i in 1:length(alpha_gamma_R)) alpha_gamma_R[,i+1] <- out_nocov_alpha_gamma[[i]][[1]]$R

### negative covariation
# susceptible
alpha_gamma_S_negcov <- array(NA, dim=c(length(timeSeq), length(out_negcov_alpha_gamma)+1))
alpha_gamma_S_negcov[,1] <- timeSeq
for (i in 1:length(alpha_gamma_S_negcov)) alpha_gamma_S_negcov[,i+1] <- out_negcov_alpha_gamma[[i]][[1]]$S

# infected
alpha_gamma_I_negcov <- array(NA, dim=c(length(timeSeq), length(out_negcov_alpha_gamma)+1))
alpha_gamma_I_negcov[,1] <- timeSeq
for (i in 1:length(alpha_gamma_I_negcov)) alpha_gamma_I_negcov[,i+1] <- out_negcov_alpha_gamma[[i]][[1]]$I

# recovered
alpha_gamma_R_negcov <- array(NA, dim=c(length(timeSeq), length(out_negcov_alpha_gamma)+1))
alpha_gamma_R_negcov[,1] <- timeSeq
for (i in 1:length(alpha_gamma_R_negcov)) alpha_gamma_R_negcov[,i+1] <- out_negcov_alpha_gamma[[i]][[1]]$R

### positive covariation
# susceptible
alpha_gamma_S_poscov <- array(NA, dim=c(length(timeSeq), length(out_poscov_alpha_gamma)+1))
alpha_gamma_S_poscov[,1] <- timeSeq
for (i in 1:length(alpha_gamma_S_poscov)) alpha_gamma_S_poscov[,i+1] <- out_poscov_alpha_gamma[[i]][[1]]$S

# infected
alpha_gamma_I_poscov <- array(NA, dim=c(length(timeSeq), length(out_poscov_alpha_gamma)+1))
alpha_gamma_I_poscov[,1] <- timeSeq
for (i in 1:length(alpha_gamma_I_poscov)) alpha_gamma_I_poscov[,i+1] <- out_poscov_alpha_gamma[[i]][[1]]$I

# recovered
alpha_gamma_R_poscov <- array(NA, dim=c(length(timeSeq), length(out_poscov_alpha_gamma)+1))
alpha_gamma_R_poscov[,1] <- timeSeq
for (i in 1:length(alpha_gamma_R_poscov)) alpha_gamma_R_poscov[,i+1] <- out_poscov_alpha_gamma[[i]][[1]]$R

### Plot 
par(mfrow=c(1,3))
## No Cov
plot(0:100, apply(alpha_gamma_S, 1, mean), col="blue", lwd=1.75, type="l", ylim=c(0,250), ylab="S", xlab="Time", main="No Cov")
lines(0:100, apply(alpha_gamma_I, 1, mean), col="red", lwd=1.75, type="l", ylim=c(0,250), ylab="S", xlab="Time")
lines(0:100, apply(alpha_gamma_R, 1, mean), col="green", lwd=1.75, type="l", ylim=c(0,250), ylab="S", xlab="Time")
lines(0:150, apply(storeMatrix.novar.S, 1, mean), col="black", lwd=1.75, type="l", ylim=c(0,250), ylab="S", xlab="Time")
abline(h=63, lty=2)

## Neg Cov
plot(0:100, apply(alpha_gamma_S_negcov, 1, mean), col="blue", lwd=1.75, type="l", ylim=c(0,250), ylab="S", xlab="Time", main="Neg Cov")
lines(0:100, apply(alpha_gamma_I_negcov, 1, mean), col="red", lwd=1.75, type="l", ylim=c(0,250), ylab="S", xlab="Time")
lines(0:100, apply(alpha_gamma_R_negcov, 1, mean), col="green", lwd=1.75, type="l", ylim=c(0,250), ylab="S", xlab="Time")
lines(0:150, apply(storeMatrix.novar.S, 1, mean), col="black", lwd=1.75, type="l", ylim=c(0,250), ylab="S", xlab="Time")
abline(h=63, lty=2)

## Pos Cov
plot(0:100, apply(alpha_gamma_S_poscov, 1, mean), col="blue", lwd=1.75, type="l", ylim=c(0,250), ylab="S", xlab="Time", main="Pos Cov")
lines(0:100, apply(alpha_gamma_I_poscov, 1, mean), col="red", lwd=1.75, type="l", ylim=c(0,250), ylab="S", xlab="Time")
lines(0:100, apply(alpha_gamma_R_poscov, 1, mean), col="green", lwd=1.75, type="l", ylim=c(0,250), ylab="S", xlab="Time")
lines(0:150, apply(storeMatrix.novar.S, 1, mean), col="black", lwd=1.75, type="l", ylim=c(0,250), ylab="S", xlab="Time")
abline(h=63, lty=2)


## No Variation Format Output
## no correlation
## infected
timeSeq1<-1:151
storeMatrix.novar.I <- array(NA, dim=c(length(timeSeq1),length(out_no_var)))
for (j in 1:length(out_no_var)) {
  o <- out_no_var[[j]]
  storeMatrix.novar.I[,j] <- o[,2]
}

## susceptible
storeMatrix.novar.S <- array(NA, dim=c(length(timeSeq1),length(out_no_var)))
for (j in 1:length(out_no_var)) {
  o <- out_no_var[[j]]
  storeMatrix.novar.S[,j] <- o[,2]
}

## recovered
storeMatrix.novar.R <- array(NA, dim=c(length(timeSeq1),length(out_no_var)))
for (j in 1:length(out_no_var)) {
  o <- out_no_var[[j]]
  storeMatrix.novar.R[,j] <- o[,4] 
}



########################################### ORIGINAL GEM CODE ###########################################
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



###### FORMAT OUTPUT ######
## no correlation
## infected
storeMatrix.cov.AG.I <- array(NA, dim=c(length(timeSeq),length(out_nocov_alpha_gamma)))
for (j in 1:length(out_nocov_alpha_gamma)) {
  o <- out_nocov_alpha_gamma[[j]]
  storeMatrix.cov.AG.I[,j] <- o[,2]
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
storeMatrix.novar.I <- array(NA, dim=c(length(timeSeq),length(out_no_var)))
for (j in 1:length(out_no_var)) {
  o <- out_no_var[[j]]
  storeMatrix.novar.I[,j] <- o[,3] 
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

