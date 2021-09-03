#### SHEDDING VARIATION COMPARISON CODE

library(tidyverse)
library(parallel)

## set up
seeds <- floor(runif(20,1,1e5)) # set seeds
tmax <- 150

analytical_parms_shed = c(c=.1, shed=.05, alpha=.1, gamma=.1, beta=.25, d=.1, 
                          b=2.5, bs=.01, varS=0.001, epsilon=0.01) # R0=3.8

analytical_parms_shed = c(c=.1, shed=.05, alpha=.1, gamma=.1, beta=.25, d=.1, 
                          b=2.5, bs=.01, varS=0.0099, epsilon=0.01) # R0=3.8



initial_state <- floor(c(S =unname(((analytical_parms_shed["b"]-analytical_parms_shed["d"])/analytical_parms_shed["bs"]))-5, I=5, R=0))

###### RUN MULTIPLE SIMULATIONS ######
## no variation
source("GEM_SIR_noVar.R")
mclapply(seeds,
         function(s) gillespie.SIR.noVar(tmax, analytical_parms_shed, initial_state),
         mc.cores=4) -> out_no_var

## stratified variation
source("GEM_SIR_shedding.variation.R")
mclapply(seeds,
         function(s) gillespie.SIR.strat.varS(tmax, analytical_parms_shed, initial_state),
         mc.cores=4) -> out_strat_var_shed

## continuous variation
source("GEM_SIR_shedding.variation.R")
mclapply(seeds,
         function(s) gillespie.SIR.varS(tmax, analytical_parms_shed, initial_state),
         mc.cores=4) -> out_cont_var_shed

###### FORMAT OUTPUT ######
## continuous variation
## infected
timeSeq <- 0:150
storeMatrix.cont.varS.I <- array(NA, dim=c(length(timeSeq),length(out_cont_var_shed)))
for (j in 1:length(out_cont_var_shed)) {
  o <- out_cont_var_shed[[j]]
  storeMatrix.cont.varS.I[,j] <- o[,3] 
}

## susceptible
storeMatrix.cont.varS.S <- array(NA, dim=c(length(timeSeq),length(out_cont_var_shed)))
for (j in 1:length(out_cont_var_shed)) {
  o <- out_cont_var_shed[[j]]
  storeMatrix.cont.varS.S[,j] <- o[,2] 
}

## recovered
storeMatrix.cont.varS.R <- array(NA, dim=c(length(timeSeq),length(out_cont_var_shed)))
for (j in 1:length(out_cont_var_shed)) {
  o <- out_cont_var_shed[[j]]
  storeMatrix.cont.varS.R[,j] <- o[,4] 
}

## stratified variation
## infected
timeSeq <- 0:150
storeMatrix.strat.varS.I <- array(NA, dim=c(length(timeSeq),length(out_strat_var_shed)))
for (j in 1:length(out_strat_var_shed)) {
  o <- out_strat_var_shed[[j]]
  storeMatrix.strat.varS.I[,j] <- o[,3] 
}

## susceptible
storeMatrix.strat.varS.S <- array(NA, dim=c(length(timeSeq),length(out_strat_var_shed)))
for (j in 1:length(out_strat_var_shed)) {
  o <- out_strat_var_shed[[j]]
  storeMatrix.strat.varS.S[,j] <- o[,2] 
}

## recovered
storeMatrix.strat.varS.R <- array(NA, dim=c(length(timeSeq),length(out_strat_var_shed)))
for (j in 1:length(out_strat_var_shed)) {
  o <- out_strat_var_shed[[j]]
  storeMatrix.strat.varS.R[,j] <- o[,4]
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
par(mfrow=c(1,2))

## no variation
plot(0:150, apply(storeMatrix.no.var.I, 1, mean), col="red", lwd=1.75, type="l", ylim=c(0,250), ylab="N", xlab="Time", main="No Var")
lines(0:150, apply(storeMatrix.no.var.S, 1, mean), col="blue", lwd=1.75, type="l", ylim=c(0,150), ylab="N", xlab="Time")
lines(0:150, apply(storeMatrix.no.var.R, 1, mean), col="green", lwd=1.75, type="l", ylim=c(0,150), ylab="N", xlab="Time")
abline(h=63, lty=2)
legend("topright",legend=c("S","I","R"),fill=c("blue","red","green"), cex=0.25)

## stratified variation
plot(0:150, apply(storeMatrix.strat.varS.I, 1, mean), col="red", lwd=1.75, type="l", ylim=c(0,250), ylab="N", xlab="Time", main="Strat Var")
lines(0:150, apply(storeMatrix.strat.varS.S, 1, mean), col="blue", lwd=1.75, type="l", ylim=c(0,150), ylab="N", xlab="Time")
lines(0:150, apply(storeMatrix.strat.varS.R, 1, mean), col="green", lwd=1.75, type="l", ylim=c(0,150), ylab="N", xlab="Time")
abline(h=63, lty=2)
legend("topright",legend=c("S","I","R"),fill=c("blue","red","green"), cex=0.25)

plot(0:150, apply(storeMatrix.no.var.S, 1, mean), col="black", lwd=1.75, type="l", ylim=c(0,200), ylab="S", xlab="Time")
lines(0:150, apply(storeMatrix.cont.varS.S[,-which(apply(storeMatrix.cont.varS.S, 2, function(col) any(is.na(col))))], 1, mean), 
      col="blue", lwd=1.75, type="l", ylim=c(0,250), ylab="N", xlab="Time")
abline(h=63, lty=2)

## continuous variation
plot(0:150, apply(storeMatrix.cont.varS.I, 1, mean), col="red", lwd=1.75, type="l", ylim=c(0,250), ylab="N", xlab="Time", main="Cont Var")
lines(0:150, apply(storeMatrix.cont.varS.S, 1, mean), col="blue", lwd=1.75, type="l", ylim=c(0,150), ylab="N", xlab="Time")
lines(0:150, apply(storeMatrix.cont.varS.R, 1, mean), col="green", lwd=1.75, type="l", ylim=c(0,150), ylab="N", xlab="Time")
abline(h=63, lty=2)
legend("topright",legend=c("S","I","R"),fill=c("blue","red","green"), cex=0.25)

## continuous variation - if you get NA
plot(0:150, apply(storeMatrix.cont.varS.I[,-which(apply(storeMatrix.cont.varS.I, 2, function(col) any(is.na(col))))], 1, mean), 
     col="red", lwd=1.75, type="l", ylim=c(0,250), ylab="N", xlab="Time", main="Cont. Var")
lines(0:150, apply(storeMatrix.cont.varS.S[,-which(apply(storeMatrix.cont.varS.S, 2, function(col) any(is.na(col))))], 1, mean), 
      col="blue", lwd=1.75, type="l", ylim=c(0,250), ylab="N", xlab="Time")
lines(0:150, apply(storeMatrix.cont.varS.R[,-which(apply(storeMatrix.cont.varS.R, 2, function(col) any(is.na(col))))], 1, mean), 
      col="green", lwd=1.75, type="l", ylim=c(0,250), ylab="N", xlab="Time")
legend("topright",legend=c("S","I","R"),fill=c("blue","red","green"), cex=0.25)
abline(h=63, lty=2)



