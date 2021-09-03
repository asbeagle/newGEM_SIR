#### BETA VARIATION COMPARISON CODE

library(tidyverse)
library(parallel)

## set up
seeds <- floor(runif(20,1,1e5)) # set seeds
tmax <- 150

betaparams =  c(alpha = .01, gamma = .01, beta = .1, d = 0.001, b = 2.5, bs = .01, varB = .01, 
                 epsilon_b = 0.01)

betaparams2 =  c(alpha = .03, gamma = .03, beta = .2, d = 0.001, b = 2.5, bs = .01, varB = .01, 
                epsilon_b = 0.01)

analytical_parms_beta = c(c=.1, shed=.05, alpha=.1, gamma=.1, beta=.25, d=.1, 
                           b=2.5, bs=.01, varC=0.01, epsilon_c=0.01) # R0=3.8


initial_state <- floor(c(S =unname(((analytical_parms_beta["b"]-analytical_parms_beta["d"])/analytical_parms_beta["bs"]))-5, I=5, R=0))

###### RUN MULTIPLE SIMULATIONS ######
## no variation
source("GEM_SIR_noVar.R")
mclapply(seeds,
         function(s) gillespie.SIR.noVar(tmax, analytical_parms_beta, initial_state, seed=s),
         mc.cores=4) -> out_no_var

## stratified variation
source("GEM_SIR_beta.variation.R")
mclapply(seeds,
         function(s) gillespie.SIR.strat.varB(tmax, analytical_parms_beta, initial_state, seed=s),
         mc.cores=4) -> out_strat_var_beta

## continuous variation
source("GEM_SIR_beta.variation.R")
mclapply(seeds,
         function(s) gillespie.SIR.varB(tmax, analytical_parms_beta, initial_state,seed=s),
         mc.cores=4) -> out_cont_var_beta

###### FORMAT OUTPUT ######
## continuous variation
## infected
timeSeq <- 0:150
storeMatrix.cont.varB.I <- array(NA, dim=c(length(timeSeq),length(out_cont_var_beta)))
for (j in 1:length(out_cont_var_beta)) {
  o <- out_cont_var_beta[[j]]
  storeMatrix.cont.varB.I[,j] <- o[,3] 
}

## susceptible
storeMatrix.cont.varB.S <- array(NA, dim=c(length(timeSeq),length(out_cont_var_beta)))
for (j in 1:length(out_cont_var_beta)) {
  o <- out_cont_var_beta[[j]]
  storeMatrix.cont.varB.S[,j] <- o[,2] 
}

## recovered
storeMatrix.cont.varB.R <- array(NA, dim=c(length(timeSeq),length(out_cont_var_beta)))
for (j in 1:length(out_cont_var_beta)) {
  o <- out_cont_var_beta[[j]]
  storeMatrix.cont.varB.R[,j] <- o[,4] 
}

## stratified variation
## infected
timeSeq <- 0:150
storeMatrix.beta.I <- array(NA, dim=c(length(timeSeq),length(out_strat_var_beta)))
for (j in 1:length(out_strat_var_beta)) {
  o <- out_strat_var_beta[[j]]
  storeMatrix.beta.I[,j] <- o[,3] 
}

## susceptible
storeMatrix.beta.S <- array(NA, dim=c(length(timeSeq),length(out_strat_var_beta)))
for (j in 1:length(out_strat_var_beta)) {
  o <- out_strat_var_beta[[j]]
  storeMatrix.beta.S[,j] <- o[,2]
}

## recovered
storeMatrix.beta.R <- array(NA, dim=c(length(timeSeq),length(out_strat_var_beta)))
for (j in 1:length(out_strat_var_beta)) {
  o <- out_strat_var_beta[[j]]
  storeMatrix.beta.R[,j] <- o[,4] 
}

## no variation
## infected
timeSeq <- 0:150
storeMatrix.no.var.I <- array(NA, dim=c(length(timeSeq),length(out_no_var)))
for (j in 1:length(out_no_var)) {
  o <- out_no_var[[j]]
  storeMatrix.no.var.I[,j] <- o[,3] 
}

## susceptible
storeMatrix.no.var.S <- array(NA, dim=c(length(timeSeq),length(out_no_var)))
for (j in 1:length(out_no_var)) {
  o <- out_no_var[[j]]
  storeMatrix.no.var.S[,j] <- o[,2] 
}

## recovered
storeMatrix.no.var.R <- array(NA, dim=c(length(timeSeq),length(out_no_var)))
for (j in 1:length(out_no_var)) {
  o <- out_no_var[[j]]
  storeMatrix.no.var.R[,j] <- o[,4] 
}

###### PLOT ######
par(mfrow=c(1,1))

par(mfrow=c(1,1), mar=c(4,4,2,2), oma=c(1,1,2,1))

## no variation
plot(0:150, apply(storeMatrix.no.var.I, 1, mean), col="red", lwd=1.75, type="l", ylim=c(0,250), ylab="N", xlab="Time", main="No Var")
lines(0:150, apply(storeMatrix.no.var.S, 1, mean), col="blue", lwd=1.75, type="l", ylim=c(0,150), ylab="N", xlab="Time")
lines(0:150, apply(storeMatrix.no.var.R, 1, mean), col="green", lwd=1.75, type="l", ylim=c(0,150), ylab="N", xlab="Time")
abline(h=63, lty=2)
#legend("topright",legend=c("S","I","R"),fill=c("blue","red","green"), cex=0.5)

plot(0:150, apply(storeMatrix.no.var.S, 1, mean), col="black", lwd=1.75, type="l", ylim=c(0,200), ylab="S", xlab="Time")
lines(0:150, apply(storeMatrix.cont.varB.S, 1, mean), col="blue", lwd=1.75, type="l", ylim=c(0,150), ylab="N", xlab="Time")
abline(h=63, lty=2)
legend("topright",legend=c("No Var","Var"),fill=c("blue","green"), cex=0.95)

## stratified variation
plot(0:150, apply(storeMatrix.beta.I, 1, mean), col="red", lwd=1.75, type="l", ylim=c(0,300), ylab="Susceptible Individuals", xlab="Time", main="Strat Var")
lines(0:150, apply(storeMatrix.beta.S, 1, mean), col="blue", lwd=1.75, type="l", ylim=c(0,150), ylab="N", xlab="Time")
lines(0:150, apply(storeMatrix.beta.R, 1, mean), col="green", lwd=1.75, type="l", ylim=c(0,150), ylab="N", xlab="Time")
legend("topright",legend=c("S","I","R"),fill=c("blue","red","green"), cex=0.25)

## continuous variation
plot(0:150, apply(storeMatrix.cont.varB.I, 1, mean), col="red", lwd=1.75, type="l", ylim=c(0,250), ylab="N", xlab="Time", main="Cont Var")
lines(0:150, apply(storeMatrix.cont.varB.S, 1, mean), col="blue", lwd=1.75, type="l", ylim=c(0,150), ylab="N", xlab="Time")
lines(0:150, apply(storeMatrix.cont.varB.R, 1, mean), col="green", lwd=1.75, type="l", ylim=c(0,150), ylab="N", xlab="Time")
abline(h=63, lty=2)
#legend("topright",legend=c("S","I","R"),fill=c("blue","red","green"), cex=0.5)

## continuous variation
plot(0:150, apply(storeMatrix.cont.varB.I[,-which(apply(storeMatrix.cont.varB.I, 2, function(col) any(is.na(col))))], 1, mean), 
     col="red", lwd=1.75, type="l", ylim=c(0,250), ylab="N", xlab="Time", main="Cont. Var")
lines(0:150, apply(storeMatrix.cont.varB.S[,-which(apply(storeMatrix.cont.varB.S, 2, function(col) any(is.na(col))))], 1, mean), 
      col="blue", lwd=1.75, type="l", ylim=c(0,250), ylab="N", xlab="Time")
lines(0:150, apply(storeMatrix.cont.varB.R[,-which(apply(storeMatrix.cont.varB.R, 2, function(col) any(is.na(col))))], 1, mean), 
      col="green", lwd=1.75, type="l", ylim=c(0,150), ylab="N", xlab="Time")
legend("topright",legend=c("S","I","R"),fill=c("blue","red","green"), cex=0.25)



