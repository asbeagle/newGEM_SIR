#### BETA VARIATION COMPARISON CODE

library(tidyverse)
library(parallel)

## set up
seeds <- floor(runif(20,1,1e5)) # set seeds
tmax <- 150

baselineparams_high = c(c=.5, shed=3.5, alpha=.15, gamma=.15, beta=.01, d=.001, 
                   b=2.5, bs=.01, varA=.15,varB=.001, varG=.15, varS=.15,
                   epsilon=0.1, epsilon_b=.001)

baselineparams_low = c(c=.5, shed=.5, alpha=.01, gamma=.01, beta=.0001, d=.001, 
                       b=2.5, bs=.01, varA=.01,varB=.0001, varG=.01, varS=.01,
                       epsilon=0.05, epsilon_b=.0001)

baselineparams_int = c(c=.5, shed=1, alpha=.1, gamma=.1, beta=.005, d=.001, 
                       b=2.5, bs=.01, varA=.1,varB=.0001, varG=.1, varS=.1,
                       epsilon=0.01, epsilon_b=.0001)



betaparams =  c(alpha = .001, gamma = .001, beta = .01, d = 0.001, b = 2.5, bs = .01, varB = .001, 
                 epsilon_b = 0.001)

initial_state <- floor(c(S =unname(((baselineparams["b"]-baselineparams["d"])/baselineparams["bs"]))-5, I=5, R=0))

###### RUN MULTIPLE SIMULATIONS ######
## no variation
source("GEM_SIR_noVar.R")
mclapply(seeds,
         function(s) gillespie.SIR.noVar(tmax, baselineparams_int, initial_state, seed=s),
         mc.cores=4) -> out_no_var

## stratified variation
source("GEM_SIR_beta.variation.R")
mclapply(seeds,
         function(s) gillespie.SIR.strat.varB(tmax, baselineparams_int, initial_state, seed=s),
         mc.cores=4) -> out_strat_var_beta

## continuous variation
source("GEM_SIR_beta.variation.R")
mclapply(seeds,
         function(s) gillespie.SIR.varB(tmax, baselineparams_int, initial_state,seed=s),
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
par(mfrow=c(1,3))

## no variation
plot(0:150, apply(storeMatrix.no.var.I, 1, mean), col="red", lwd=1.75, type="l", ylim=c(0,300), ylab="N", xlab="Time", main="No Var")
lines(0:150, apply(storeMatrix.no.var.S, 1, mean), col="blue", lwd=1.75, type="l", ylim=c(0,150), ylab="N", xlab="Time")
lines(0:150, apply(storeMatrix.no.var.R, 1, mean), col="green", lwd=1.75, type="l", ylim=c(0,150), ylab="N", xlab="Time")
legend("topright",legend=c("S","I","R"),fill=c("blue","red","green"), cex=0.25)

## stratified variation
plot(0:150, apply(storeMatrix.beta.I, 1, mean), col="red", lwd=1.75, type="l", ylim=c(0,300), ylab="N", xlab="Time", main="Strat Var")
lines(0:150, apply(storeMatrix.beta.S, 1, mean), col="blue", lwd=1.75, type="l", ylim=c(0,150), ylab="N", xlab="Time")
lines(0:150, apply(storeMatrix.beta.R, 1, mean), col="green", lwd=1.75, type="l", ylim=c(0,150), ylab="N", xlab="Time")
legend("topright",legend=c("S","I","R"),fill=c("blue","red","green"), cex=0.25)

## continuous variation
plot(0:150, apply(storeMatrix.cont.varB.I, 1, mean), col="red", lwd=1.75, type="l", ylim=c(0,300), ylab="N", xlab="Time", main="Cont Var")
lines(0:150, apply(storeMatrix.cont.varB.S, 1, mean), col="blue", lwd=1.75, type="l", ylim=c(0,150), ylab="N", xlab="Time")
lines(0:150, apply(storeMatrix.cont.varB.R, 1, mean), col="green", lwd=1.75, type="l", ylim=c(0,150), ylab="N", xlab="Time")
legend("topright",legend=c("S","I","R"),fill=c("blue","red","green"), cex=0.25)

## continuous variation
plot(0:150, apply(storeMatrix.cont.varB.I[,-which(apply(storeMatrix.cont.varB.I, 2, function(col) any(is.na(col))))], 1, mean), 
     col="red", lwd=1.75, type="l", ylim=c(0,300), ylab="N", xlab="Time", main="Cont. Var")
lines(0:150, apply(storeMatrix.cont.varB.S[,-which(apply(storeMatrix.cont.varB.S, 2, function(col) any(is.na(col))))], 1, mean), 
      col="blue", lwd=1.75, type="l", ylim=c(0,150), ylab="N", xlab="Time")
lines(0:150, apply(storeMatrix.cont.varB.R[,-which(apply(storeMatrix.cont.varB.R, 2, function(col) any(is.na(col))))], 1, mean), 
      col="green", lwd=1.75, type="l", ylim=c(0,150), ylab="N", xlab="Time")
legend("topright",legend=c("S","I","R"),fill=c("blue","red","green"), cex=0.25)



