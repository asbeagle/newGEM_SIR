#### ALPHA VARIATION COMPARISON CODE

library(tidyverse)
library(parallel)

## set up
seeds <- floor(runif(20,1,1e5)) # set seeds
tmax <- 150

analytical_parms_alpha = c(c=.1, shed=.05, alpha=.1, gamma=.1, beta=.25, d=.1, 
                           b=2.5, bs=.01, varA=0.1, epsilon=0.01) # R0=3.8


initial_state <- floor(c(S =unname(((analytical_parms_alpha["b"]-analytical_parms_alpha["d"])/analytical_parms_alpha["bs"]))-5, I=5, R=0))

###### RUN MULTIPLE SIMULATIONS ######
## no variation
source("GEM_SIR_noVar.R")
mclapply(seeds,
         function(s) gillespie.SIR.noVar(tmax, analytical_parms_alpha, initial_state, seed=s),
         mc.cores=4) -> out_no_var

## stratified variation
source("GEM_SIR_alpha.variation.R")
mclapply(seeds,
         function(s) gillespie.SIR.strat.varA(tmax, analytical_parms_alpha, initial_state,seed=s),
         mc.cores=4) -> out_strat_var_alpha

## continuous variation
source("GEM_SIR_alpha.variation.R")
mclapply(seeds,
         function(s) gillespie.SIR.varA(tmax, analytical_parms_alpha, initial_state,seed=s),
         mc.cores=4) -> out_cont_var_alpha

###### FORMAT OUTPUT ######
## continuous variation
## infected
timeSeq <- 0:150
storeMatrix.cont.varA.I <- array(NA, dim=c(length(timeSeq),length(out_cont_var_alpha)))
for (j in 1:length(out_cont_var_alpha)) {
  o <- out_cont_var_alpha[[j]]
  storeMatrix.cont.varA.I[,j] <- o[,3] 
}

## susceptible
storeMatrix.cont.varA.S <- array(NA, dim=c(length(timeSeq),length(out_cont_var_alpha)))
for (j in 1:length(out_cont_var_alpha)) {
  o <- out_cont_var_alpha[[j]]
  storeMatrix.cont.varA.S[,j] <- o[,2] 
}

## recovered
storeMatrix.cont.varA.R <- array(NA, dim=c(length(timeSeq),length(out_cont_var_alpha)))
for (j in 1:length(out_cont_var_alpha)) {
  o <- out_cont_var_alpha[[j]]
  storeMatrix.cont.varA.R[,j] <- o[,4] 
}

## stratified variation
## infected
timeSeq <- 0:150
storeMatrix.alpha.I <- array(NA, dim=c(length(timeSeq),length(out_strat_var_alpha)))
for (j in 1:length(out_strat_var_alpha)) {
  o <- out_strat_var_alpha[[j]]
  storeMatrix.alpha.I[,j] <- o[,3] 
}

## susceptible
storeMatrix.alpha.S <- array(NA, dim=c(length(timeSeq),length(out_strat_var_alpha)))
for (j in 1:length(out_strat_var_alpha)) {
  o <- out_strat_var_alpha[[j]]
  storeMatrix.alpha.S[,j] <- o[,2] 
}

## recovered
storeMatrix.alpha.R <- array(NA, dim=c(length(timeSeq),length(out_strat_var_alpha)))
for (j in 1:length(out_strat_var_alpha)) {
  o <- out_strat_var_alpha[[j]]
  storeMatrix.alpha.R[,j] <- o[,4] 
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
par(mfrow=c(1,3))

## no variation
plot(0:150, apply(storeMatrix.no.var.I, 1, mean), col="red", lwd=1.75, type="l", ylim=c(0,300), ylab="N", xlab="Time", main="No Var")
lines(0:150, apply(storeMatrix.no.var.S, 1, mean), col="blue", lwd=1.75, type="l", ylim=c(0,150), ylab="N", xlab="Time")
lines(0:150, apply(storeMatrix.no.var.R, 1, mean), col="green", lwd=1.75, type="l", ylim=c(0,150), ylab="N", xlab="Time")
abline(h=56, lty=2)
legend("topright",legend=c("S","I","R"),fill=c("blue","red","green"), cex=0.25)

## stratified variation
plot(0:150, apply(storeMatrix.alpha.I, 1, mean), col="red", lwd=1.75, type="l", ylim=c(0,300), ylab="N", xlab="Time", main="Strat Var")
lines(0:150, apply(storeMatrix.alpha.S, 1, mean), col="blue", lwd=1.75, type="l", ylim=c(0,150), ylab="N", xlab="Time")
lines(0:150, apply(storeMatrix.alpha.R, 1, mean), col="green", lwd=1.75, type="l", ylim=c(0,150), ylab="N", xlab="Time")
abline(h=56, lty=2)
legend("topright",legend=c("S","I","R"),fill=c("blue","red","green"), cex=0.25)

## continous variation
plot(0:150, apply(storeMatrix.cont.varA.I, 1, mean), col="red", lwd=1.75, type="l", ylim=c(0,300), ylab="N", xlab="Time", main="Cont Var")
lines(0:150, apply(storeMatrix.cont.varA.S, 1, mean), col="blue", lwd=1.75, type="l", ylim=c(0,150), ylab="N", xlab="Time")
lines(0:150, apply(storeMatrix.cont.varA.R, 1, mean), col="green", lwd=1.75, type="l", ylim=c(0,150), ylab="N", xlab="Time")
abline(h=56, lty=2)
legend("topright",legend=c("S","I","R"),fill=c("blue","red","green"), cex=0.25)



