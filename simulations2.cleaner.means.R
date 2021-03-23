###### PREVIOUS PARAMS ######
#params = c(c=.2, shed=.2, sd_s=.1, sd_a=.1, sd_c=.01, sd_g=.01, 
           #h=.1, alpha=.01, gamma=.3, b=2.5, d=.4, bs=.01) # R0 = 13
#params2 = c(c=.2, shed=.2, sd_s=.25, sd_a=.25, sd_c=.25, sd_g=.25, 
            #h=.1, alpha=.01, gamma=.3, b=2.5, d=.4, bs=.01) # R0 = 13
#params3 = c(c=.4, shed=.4, sd_s=.25, sd_a=.25, sd_c=.25, sd_g=.25, 
            #h=.2, alpha=.1, gamma=.15, b=2.5, d=.4, bs=.01) # R0 = 28

#params6 = c(c=.05, shed=.05, sd_s=.05, sd_a=.1, sd_c=.05, sd_g=.15, 
            #h=.15, alpha=.1, gamma=.15, b=2.5, d=.4, bs=.01) # R0 = 1.326 
#params7 = c(c=.07, shed=.07, sd_s=.07, sd_a=.1, sd_c=.07, sd_g=.15, 
            #h=.13, alpha=.1, gamma=.15, b=2.5, d=.4, bs=.01) # R0 = 2.522
#params8 = c(c=.085, shed=.089, sd_s=.089, sd_a=.13, sd_c=.085, sd_g=.15, 
            #h=.13, alpha=.13, gamma=.15, b=2.5, d=.4, bs=.01) # R0 = 3.55

### covariation baseline params
#baselineparams = c(c=.035, shed=.05, h=.15, alpha=.15, gamma=.15, beta=.25, d=.2, 
#                   sd_c=.035, sd_s=.05, sd_a=.15, sd_g=.15, b=2.5, bs=.01, varA=1e-3, 
#                   varB=1e-3, varG=1e-3) # R0 = 1.225

### old baseline, changed 3/23/21
### baselineparams.old = c(c=.5, shed=.5, h=.15, alpha=.15, gamma=.15, beta=.25, d=.2, 
#                       sd_c=.035, sd_s=.05, sd_a=.15, sd_g=.15, b=2.5, bs=.01, varA=1e-3, 
#                       varB=1e-3, varG=1e-3, varS=1e-3)
#x = c(S=70, I=10, R=0)

## R0 CALCULATION
c=.1
shed=.05
h=.15
alpha=.15
gamma=.15
d=.2
R0<- (((c*(shed/(shed+h)))*70)/(alpha+gamma+d))

# packages
library(tidyverse)
library(parallel)

# seeds and current parameter values
seeds <- floor(runif(20,1,1e5)) # set seeds
tmax <- 150

## params_novar <- c(beta=0.0025, alpha=0.15, gamma=0.001, varG=1e-3, b=2.5, d=0.001, bs=0.01, ds=0.01, varA=0)
baselineparams = c(c=.5, shed=1, alpha=.15, gamma=.001, beta=.0025, d=.001, 
                   sd_c=.035, sd_s=.05, sd_a=.15, sd_g=.15, b=2.5, bs=.01, varA=.15, 
                   varB=.001, varG=.15, varS=.15,epsilon=0.1, epsilon_b=.001)

initial_state <- floor(c(S =unname(((baselineparams["b"]-baselineparams["d"])/baselineparams["bs"]))-5, I=5, R=0))


contact.tau.params1 = c(c=.035, shed=.05, h=.15, alpha=.12, gamma=.05, d=.07,
                        sd_c=.035, sd_s=.05, sd_a=.12, sd_g=.05,b=2.5, bs=.01) # R0 = 2.55
contact.tau.params2 = c(c=.035, shed=.05, h=.15, alpha=.1, gamma=.005, d=.07,
                        sd_c=.035, sd_s=.05, sd_a=.1, sd_g=.005,b=2.5, bs=.01) # R0 = 3.5

## Note: you reduced alpha from the baseline here, and these two parameter sets are identical
contact.alpha.params1 = c(c=.035, shed=.22, h=.2, alpha=.15, gamma=.15, d=.2,
                          sd_c=.035, sd_s=.22, sd_a=.15, sd_g=.15, b=2.5, bs=.01) # R0 = 2.566
contact.alpha.params2 = c(c=.035, shed=.3, h=.12, alpha=.15, gamma=.15, d=.2,
                          sd_c=.035, sd_s=.3, sd_a=.15, sd_g=.15, b=2.5, bs=.01) # R0 = 3.5

tau.alpha.params1 = c(c=.074, shed=.05, h=.15, alpha=.15, gamma=.15, d=.2,
                      sd_c=.074, sd_s=.05, sd_a=.15, sd_g=.15, b=2.5, bs=.01) # R0 = 2.59
tau.alpha.params2 = c(c=.1, shed=.05, h=.15, alpha=.15, gamma=.15, d=.2,
                      sd_c=.074, sd_s=.05, sd_a=.15, sd_g=.15, b=2.5, bs=.01) # R0 = 3.5


nocorr <- matrix(c(1,0,0,1), nrow=2, byrow=T)
negcorr <- matrix(c(1,-.5,-.5,1), nrow=2, byrow=T)
poscorr <- matrix(c(1,.5,.5,1), nrow=2, byrow=T)

#  model simulations
## no covariance
mclapply(seeds,
         function(s) gillespie.SIR.cov_taualpha(tmax, tau.alpha.params2, nocorr, x),
         mc.cores=4) -> out_nocov_taualpha

#mclapply(seeds,
         #function(s) gillespie.SIR.cov_taugamma(tmax, baselineparams, nocorr, x),
         #mc.cores=4) -> out_nocov_taugamma

#mclapply(seeds,
         #function(s) gillespie.SIR.cov_cgamma(tmax, baselineparams, nocorr, x),
         #mc.cores=4) -> out_nocov_cgamma

mclapply(seeds,
         function(s) gillespie.SIR.cov_calpha(tmax, contact.alpha.params2, nocorr, x),
         mc.cores=4) -> out_nocov_calpha

mclapply(seeds,
         function(s) gillespie.SIR.cov_ctau(tmax, contact.tau.params2, nocorr, x),
         mc.cores=4) -> out_nocov_ctau

## negative covariance
mclapply(seeds,
         function(s) gillespie.SIR.cov_taualpha(tmax, tau.alpha.params2, negcorr, x),
         mc.cores=4) -> out_negcov_taualpha

#mclapply(seeds,
         #function(s) gillespie.SIR.cov_taugamma(tmax, params8, negcorr, x),
         #mc.cores=4) -> out_negcov_taugamma

#mclapply(seeds,
         #function(s) gillespie.SIR.cov_cgamma(tmax, params8, negcorr, x),
         #mc.cores=4) -> out_negcov_cgamma

mclapply(seeds,
         function(s) gillespie.SIR.cov_calpha(tmax, contact.alpha.params2, negcorr, x),
         mc.cores=4) -> out_negcov_calpha

mclapply(seeds,
         function(s) gillespie.SIR.cov_ctau(tmax, contact.tau.params2, negcorr, x),
         mc.cores=4) -> out_negcov_ctau

## positive cov
mclapply(seeds,
         function(s) gillespie.SIR.cov_taualpha(tmax, tau.alpha.params2, poscorr, x),
         mc.cores=4) -> out_poscov_taualpha

#mclapply(seeds,
         #function(s) gillespie.SIR.cov_taugamma(tmax, params8, poscorr, x),
         #mc.cores=4) -> out_poscov_taugamma

#mclapply(seeds,
         #function(s) gillespie.SIR.cov_cgamma(tmax, params8, poscorr, x),
         #mc.cores=4) -> out_poscov_cgamma

mclapply(seeds,
         function(s) gillespie.SIR.cov_calpha(tmax, contact.alpha.params, poscorr, x),
         mc.cores=4) -> out_poscov_calpha

mclapply(seeds,
         function(s) gillespie.SIR.cov_ctau(tmax, contact.tau.params, poscorr, x),
         mc.cores=4) -> out_poscov_ctau

### SIMS WITH STRATIFIED VARIATION
source("GEM_SIR_alpha.variation.R")
mclapply(seeds,
         function(s) gillespie.SIR.strat.varA(tmax, baselineparams, initial_state),
         mc.cores=4) -> out_strat_var_alpha

mclapply(seeds,
         function(s) gillespie.SIR.strat.varB(tmax, baselineparams, initial_state),
         mc.cores=4) -> out_strat_var_beta

mclapply(seeds,
         function(s) gillespie.SIR.strat.varG(tmax, baselineparams, initial_state),
         mc.cores=4) -> out_strat_var_gamma

mclapply(seeds,
         function(s) gillespie.SIR.strat.varS(tmax, baselineparams, initial_state),
         mc.cores=4) -> out_strat_var_shed

### SIM WITH NO VARIATION
mclapply(seeds,
         function(s) gillespie.SIR.noVar(tmax, baselineparams, initial_state),
         mc.cores=4) -> out_no_var

### SIMs WITH CONTINUOUS VARIATION
mclapply(seeds,
         function(s) gillespie.SIR.varB(tmax, baselineparams, initial_state),
         mc.cores=4) -> out_cont_var_beta

mclapply(seeds,
         function(s) gillespie.SIR.varG(tmax, baselineparams, initial_state),
         mc.cores=4) -> out_cont_var_gamma

mclapply(seeds,
         function(s) gillespie.SIR.varA(tmax, baselineparams, initial_state),
         mc.cores=4) -> out_cont_var_alpha

mclapply(seeds,
         function(s) gillespie.SIR.varS(tmax, baselineparams, initial_state),
         mc.cores=4) -> out_cont_var_shed

### SIMS WITH CONTINUOUS UNIFORM VARIATION
mclapply(seeds,
         function(s) gillespie.SIR.varB.uniform(tmax, baselineparams, initial_state),
         mc.cores=4) -> out_uniform_var_beta

mclapply(seeds,
         function(s) gillespie.SIR.varG.uniform(tmax, baselineparams, initial_state),
         mc.cores=4) -> out_uniform_var_gamma

mclapply(seeds,
         function(s) gillespie.SIR.varA.uniform(tmax, baselineparams, initial_state),
         mc.cores=4) -> out_uniform_var_alpha

mclapply(seeds,
         function(s) gillespie.SIR.varS.uniform(tmax, baselineparams, initial_state),
         mc.cores=4) -> out_uniform_var_shed

############################       ALPHA     #################################
############################ STRAT VARIATION #################################

## infected
timeSeq <- 0:150
storeMatrix.alpha.I <- array(NA, dim=c(length(timeSeq),length(out_strat_var_alpha)))
for (j in 1:length(out_strat_var_alpha)) {
  o <- out_strat_var_alpha[[j]]
  storeMatrix.alpha.I[,j] <- o[,3] # num infected
}

## susceptible
storeMatrix.alpha.S <- array(NA, dim=c(length(timeSeq),length(out_strat_var_alpha)))
for (j in 1:length(out_strat_var_alpha)) {
  o <- out_strat_var_alpha[[j]]
  storeMatrix.alpha.S[,j] <- o[,2] # num susceptible
}

## recovered
storeMatrix.alpha.R <- array(NA, dim=c(length(timeSeq),length(out_strat_var_alpha)))
for (j in 1:length(out_strat_var_alpha)) {
  o <- out_strat_var_alpha[[j]]
  storeMatrix.alpha.R[,j] <- o[,4] # num recovered
}

par(mfrow=c(1,3))
plot(0:150, apply(storeMatrix.alpha.I, 1, mean), col="red", lwd=1.75, type="l", ylim=c(0,200), ylab="N", xlab="Time", main="Variation in Alpha")
lines(0:150, apply(storeMatrix.alpha.S, 1, mean), col="blue", lwd=1.75, type="l", ylim=c(0,150), ylab="N", xlab="Time")
lines(0:150, apply(storeMatrix.alpha.R, 1, mean), col="green", lwd=1.75, type="l", ylim=c(0,150), ylab="N", xlab="Time")
legend("topright",legend=c("S","I","R"),fill=c("blue","red","green"), cex=0.5)

############################       BETA     #################################
############################ STRAT VARIATION #################################

## infected
timeSeq <- 0:150
storeMatrix.beta.I <- array(NA, dim=c(length(timeSeq),length(out_strat_var_beta)))
for (j in 1:length(out_strat_var_beta)) {
  o <- out_strat_var_beta[[j]]
  storeMatrix.beta.I[,j] <- o[,3] # num infected
}

## susceptible
storeMatrix.beta.S <- array(NA, dim=c(length(timeSeq),length(out_strat_var_beta)))
for (j in 1:length(out_strat_var_beta)) {
  o <- out_strat_var_beta[[j]]
  storeMatrix.beta.S[,j] <- o[,2] # num susceptible
}

## recovered
storeMatrix.beta.R <- array(NA, dim=c(length(timeSeq),length(out_strat_var_beta)))
for (j in 1:length(out_strat_var_beta)) {
  o <- out_strat_var_beta[[j]]
  storeMatrix.beta.R[,j] <- o[,4] # num recovered
}

plot(0:150, apply(storeMatrix.beta.I, 1, mean), col="red", lwd=1.75, type="l", ylim=c(0,250), ylab="N", xlab="Time", main="Variation in Beta")
lines(0:150, apply(storeMatrix.beta.S, 1, mean), col="blue", lwd=1.75, type="l", ylim=c(0,150), ylab="N", xlab="Time")
lines(0:150, apply(storeMatrix.beta.R, 1, mean), col="green", lwd=1.75, type="l", ylim=c(0,150), ylab="N", xlab="Time")
legend("topright",legend=c("S","I","R"),fill=c("blue","red","green"), cex=0.5)

############################       GAMMA     #################################
############################ STRAT VARIATION #################################

## infected
timeSeq <- 0:150
storeMatrix.gamma.I <- array(NA, dim=c(length(timeSeq),length(out_strat_var_gamma)))
for (j in 1:length(out_strat_var_gamma)) {
  o <- out_strat_var_gamma[[j]]
  storeMatrix.gamma.I[,j] <- o[,3] # num infected
}

## susceptible
storeMatrix.gamma.S <- array(NA, dim=c(length(timeSeq),length(out_strat_var_gamma)))
for (j in 1:length(out_strat_var_gamma)) {
  o <- out_strat_var_gamma[[j]]
  storeMatrix.gamma.S[,j] <- o[,2] # num susceptible
}

## recovered
storeMatrix.gamma.R <- array(NA, dim=c(length(timeSeq),length(out_strat_var_gamma)))
for (j in 1:length(out_strat_var_gamma)) {
  o <- out_strat_var_gamma[[j]]
  storeMatrix.gamma.R[,j] <- o[,4] # num recovered
}

plot(0:150, apply(storeMatrix.gamma.I, 1, mean), col="red", lwd=1.75, type="l", ylim=c(0,250), ylab="N", xlab="Time", main="Variation in Gamma")
lines(0:150, apply(storeMatrix.gamma.S, 1, mean), col="blue", lwd=1.75, type="l", ylim=c(0,150), ylab="N", xlab="Time")
lines(0:150, apply(storeMatrix.gamma.R, 1, mean), col="green", lwd=1.75, type="l", ylim=c(0,150), ylab="N", xlab="Time")
legend("topright",legend=c("S","I","R"),fill=c("blue","red","green"), cex=0.5)

############################  TAU & ALPHA   #################################
############################ NO COVARIATION #################################

## infected
timeSeq <- 0:150
storeMatrix.taualpha.I <- array(NA, dim=c(length(timeSeq),length(out_nocov_taualpha)))
for (j in 1:length(out_nocov_taualpha)) {
  o <- out_nocov_taualpha[[j]]
  storeMatrix.taualpha.I[,j] <- o[,3] # num infected
}

## susceptible
storeMatrix.taualpha.S <- array(NA, dim=c(length(timeSeq),length(out_nocov_taualpha)))
for (j in 1:length(out_nocov_taualpha)) {
  o <- out_nocov_taualpha[[j]]
  storeMatrix.taualpha.S[,j] <- o[,2] # num susceptible
}

## recovered
storeMatrix.taualpha.R <- array(NA, dim=c(length(timeSeq),length(out_nocov_taualpha)))
for (j in 1:length(out_nocov_taualpha)) {
  o <- out_nocov_taualpha[[j]]
  storeMatrix.taualpha.R[,j] <- o[,4] # num susceptible
}

############################ NEG COVARIATION #################################
timeSeq <- 0:150
storeMatrix.taualpha.S2 <- array(NA, dim=c(length(timeSeq),length(out_negcov_taualpha)))
for (j in 1:length(out_negcov_taualpha)) {
  o <- out_negcov_taualpha2[[j]]
  storeMatrix.taualpha.S2[,j] <- o[,2] # num susceptible
}

storeMatrix.taualpha.I2 <- array(NA, dim=c(length(timeSeq),length(out_negcov_taualpha)))
for (j in 1:length(out_negcov_taualpha)) {
  o <- out_negcov_taualpha2[[j]]
  storeMatrix.taualpha.I2[,j] <- o[,3] # num infected
}

storeMatrix.taualpha.R2 <- array(NA, dim=c(length(timeSeq),length(out_negcov_taualpha)))
for (j in 1:length(out_negcov_taualpha)) {
  o <- out_negcov_taualpha[[j]]
  storeMatrix.taualpha.R2[,j] <- o[,4] # num recovered
}


############################ POS COVARIATION #################################
timeSeq <- 0:150
storeMatrix.taualpha.S3 <- array(NA, dim=c(length(timeSeq),length(out_poscov_taualpha)))
for (j in 1:length(out_poscov_taualpha)) {
  o <- out_poscov_taualpha[[j]]
  storeMatrix.taualpha.S3[,j] <- o[,2] # num susceptible
}

storeMatrix.taualpha.I3 <- array(NA, dim=c(length(timeSeq),length(out_poscov_taualpha)))
for (j in 1:length(out_poscov_taualpha)) {
  o <- out_poscov_taualpha[[j]]
  storeMatrix.taualpha.I3[,j] <- o[,3] # num infected
}

storeMatrix.taualpha.R3 <- array(NA, dim=c(length(timeSeq),length(out_poscov_taualpha)))
for (j in 1:length(out_poscov_taualpha)) {
  o <- out_poscov_taualpha[[j]]
  storeMatrix.taualpha.R3[,j] <- o[,4] # num recovered
}

#### PLOTS
par(mfrow=c(1,4))
## mean values 

# NO COV
plot(0:150, apply(storeMatrix.taualpha.I, 1, mean), col="red", lwd=1.75, type="l", ylim=c(0,200), ylab="N", xlab="Time", main="No Cov")
lines(0:150, apply(storeMatrix.taualpha.S, 1, mean), col="blue", lwd=1.75, type="l", ylim=c(0,150), ylab="N", xlab="Time", main="Susceptible")
lines(0:150, apply(storeMatrix.taualpha.R, 1, mean), col="green", lwd=1.75, type="l", ylim=c(0,150), ylab="N", xlab="Time", main="Recovered")
legend("topright",legend=c("S","I","R"),fill=c("blue","red","green"), cex=0.5)

# NEG COV
plot(0:150, apply(storeMatrix.taualpha.I2, 1, mean), col="red", lwd=1.75, type="l", ylim=c(0,200), ylab="N", xlab="Time", main="Neg Cov")
lines(0:150, apply(storeMatrix.taualpha.S2, 1, mean), col="blue", lwd=1.75, type="l", ylim=c(0,150), ylab="N", xlab="Time", main="Susceptible")
lines(0:150, apply(storeMatrix.taualpha.R2, 1, mean), col="green", lwd=1.75, type="l", ylim=c(0,150), ylab="N", xlab="Time", main="Recovered")
legend("topright",legend=c("S","I","R"),fill=c("blue","red","green"), cex=0.5)

# POS COV
plot(0:150, apply(storeMatrix.taualpha.I3, 1, mean), col="red", lwd=1.75, type="l", ylim=c(0,200), ylab="N", xlab="Time", main="Pos Cov")
lines(0:150, apply(storeMatrix.taualpha.S3, 1, mean), col="blue", lwd=1.75, type="l", ylim=c(0,150), ylab="N", xlab="Time", main="Susceptible")
lines(0:150, apply(storeMatrix.taualpha.R3, 1, mean), col="green", lwd=1.75, type="l", ylim=c(0,150), ylab="N", xlab="Time", main="Recovered")
legend("topright",legend=c("S","I","R"),fill=c("blue","red","green"), cex=0.55)

# DETERMINISTIC
out=lsoda(xstart,times, deterministic.SIR, tau.alpha.params2)

S=out[,2]
I=out[,3]
R=out[,4]

plot.ts(S, col="blue", ylim=c(0,200), ylab="N", xlab="Time", main="Deterministic")
lines(I, col="red")
lines(R, col="green")
legend("topright",legend=c("S","I","R"),fill=c("blue","red","green"), cex=.5)


############################  TAU & GAMMA   #################################
############################ NO COVARIATION #################################
storeMatrix.taugamma.S <- array(NA, dim=c(length(timeSeq),length(out_nocov_taugamma)))
for (j in 1:length(out_nocov_taugamma)) {
  o <- out_nocov_taugamma[[j]]
  storeMatrix.taugamma.S[,j] <- o[,2] # num susceptible
}

storeMatrix.taugamma.I <- array(NA, dim=c(length(timeSeq),length(out_nocov_taugamma)))
for (j in 1:length(out_nocov_taugamma)) {
  o <- out_nocov_taugamma[[j]]
  storeMatrix.taugamma.I[,j] <- o[,3] # num infected
}

storeMatrix.taugamma.R <- array(NA, dim=c(length(timeSeq),length(out_nocov_taugamma)))
for (j in 1:length(out_nocov_taugamma)) {
  o <- out_nocov_taugamma[[j]]
  storeMatrix.taugamma.R[,j] <- o[,4] # num recovered
}
############################ NEG COVARIATION #################################
storeMatrix.taugamma.S2 <- array(NA, dim=c(length(timeSeq),length(out_negcov_taugamma)))
for (j in 1:length(out_negcov_taugamma)) {
  o <- out_negcov_taugamma[[j]]
  storeMatrix.taugamma.S2[,j] <- o[,2] # num susceptible
}

storeMatrix.taugamma.I2 <- array(NA, dim=c(length(timeSeq),length(out_negcov_taugamma)))
for (j in 1:length(out_negcov_taugamma)) {
  o <- out_negcov_taugamma[[j]]
  storeMatrix.taugamma.I2[,j] <- o[,3] # num infected
}

storeMatrix.taugamma.R2 <- array(NA, dim=c(length(timeSeq),length(out_negcov_taugamma)))
for (j in 1:length(out_negcov_taugamma)) {
  o <- out_negcov_taugamma[[j]]
  storeMatrix.taugamma.R2[,j] <- o[,4] # num recovered
}

############################ POS COVARIATION #################################
storeMatrix.taugamma.S3 <- array(NA, dim=c(length(timeSeq),length(out_poscov_taugamma)))
for (j in 1:length(out_poscov_taugamma)) {
  o <- out_poscov_taugamma[[j]]
  storeMatrix.taugamma.S3[,j] <- o[,2] # num susceptible
}

storeMatrix.taugamma.I3 <- array(NA, dim=c(length(timeSeq),length(out_poscov_taugamma)))
for (j in 1:length(out_poscov_taugamma)) {
  o <- out_poscov_taugamma[[j]]
  storeMatrix.taugamma.I3[,j] <- o[,3] # num infected
}

storeMatrix.taugamma.R3 <- array(NA, dim=c(length(timeSeq),length(out_poscov_taugamma)))
for (j in 1:length(out_poscov_taugamma)) {
  o <- out_poscov_taugamma[[j]]
  storeMatrix.taugamma.R3[,j] <- o[,4] # num recovered
}

## PLOTS
par(mfrow=c(1,4))
## mean values 

# NO COV
plot(0:150, apply(storeMatrix.taugamma.I, 1, mean), col="red", lwd=1.75, type="l", ylim=c(0,200), ylab="N", xlab="Time", main="No cov")
lines(0:150, apply(storeMatrix.taugamma.S, 1, mean), col="blue", lwd=1.75, type="l", ylim=c(0,150), ylab="N", xlab="Time", main="Susceptible")
lines(0:150, apply(storeMatrix.taugamma.R, 1, mean), col="green", lwd=1.75, type="l", ylim=c(0,150), ylab="N", xlab="Time", main="Recovered")
legend("topright",legend=c("S","I","R"),fill=c("blue","red","green"), cex=0.75)

# NEG COV
plot(0:150, apply(storeMatrix.taugamma.I2, 1, mean), col="red", lwd=1.75, type="l", ylim=c(0,200), ylab="N", xlab="Time", main="Neg cov")
lines(0:150, apply(storeMatrix.taugamma.S2, 1, mean), col="blue", lwd=1.75, type="l", ylim=c(0,150), ylab="N", xlab="Time", main="Susceptible")
lines(0:150, apply(storeMatrix.taugamma.R2, 1, mean), col="green", lwd=1.75, type="l", ylim=c(0,150), ylab="N", xlab="Time", main="Recovered")
legend("topright",legend=c("S","I","R"),fill=c("blue","red","green"), cex=0.75)

# POS COV
plot(0:150, apply(storeMatrix.taugamma.I3, 1, mean), col="red", lwd=1.75, type="l", ylim=c(0,200), ylab="N", xlab="Time", main="Pos cov")
lines(0:150, apply(storeMatrix.taugamma.S3, 1, mean), col="blue", lwd=1.75, type="l", ylim=c(0,150), ylab="N", xlab="Time", main="Susceptible")
lines(0:150, apply(storeMatrix.taugamma.R3, 1, mean), col="green", lwd=1.75, type="l", ylim=c(0,150), ylab="N", xlab="Time", main="Recovered")
legend("topright",legend=c("S","I","R"),fill=c("blue","red","green"), cex=0.75)

# DETERMINISTIC
plot.ts(S, col="blue", ylim=c(0,200), ylab="N", xlab="Time", main="Deterministic")
lines(I, col="red")
lines(R, col="green")
legend("topright",legend=c("S","I","R"),fill=c("blue","red","green"), cex=.7)

############################  CONTACT & GAMMA   #################################
############################   NO COVARIATION   #################################

storeMatrix.cgamma.S <- array(NA, dim=c(length(timeSeq),length(out_nocov_cgamma)))
for (j in 1:length(out_nocov_cgamma)) {
  o <- out_nocov_cgamma[[j]]
  storeMatrix.cgamma.S[,j] <- o[,2] # num susceptible
}

storeMatrix.cgamma.I <- array(NA, dim=c(length(timeSeq),length(out_nocov_cgamma)))
for (j in 1:length(out_nocov_cgamma)) {
  o <- out_nocov_cgamma[[j]]
  storeMatrix.cgamma.I[,j] <- o[,3] # num infected
}

storeMatrix.cgamma.R <- array(NA, dim=c(length(timeSeq),length(out_nocov_cgamma)))
for (j in 1:length(out_nocov_cgamma)) {
  o <- out_nocov_cgamma[[j]]
  storeMatrix.cgamma.R[,j] <- o[,4] # num recovered
}

############################  NEG COVARIATION   #################################
storeMatrix.cgamma.S2 <- array(NA, dim=c(length(timeSeq),length(out_negcov_cgamma)))
for (j in 1:length(out_negcov_cgamma)) {
  o <- out_negcov_cgamma[[j]]
  storeMatrix.cgamma.S2[,j] <- o[,2] # num susceptible
}

storeMatrix.cgamma.I2 <- array(NA, dim=c(length(timeSeq),length(out_negcov_cgamma)))
for (j in 1:length(out_negcov_cgamma)) {
  o <- out_negcov_cgamma[[j]]
  storeMatrix.cgamma.I2[,j] <- o[,3] # num infected
}

storeMatrix.cgamma.R2 <- array(NA, dim=c(length(timeSeq),length(out_negcov_cgamma)))
for (j in 1:length(out_negcov_cgamma)) {
  o <- out_negcov_cgamma[[j]]
  storeMatrix.cgamma.R2[,j] <- o[,4] # num recovered
}

############################  POS COVARIATION   #################################
storeMatrix.cgamma.S3 <- array(NA, dim=c(length(timeSeq),length(out_poscov_cgamma)))
for (j in 1:length(out_poscov_cgamma)) {
  o <- out_poscov_cgamma[[j]]
  storeMatrix.cgamma.S3[,j] <- o[,2] # num susceptible
}

storeMatrix.cgamma.I3 <- array(NA, dim=c(length(timeSeq),length(out_poscov_cgamma)))
for (j in 1:length(out_poscov_cgamma)) {
  o <- out_poscov_cgamma[[j]]
  storeMatrix.cgamma.I3[,j] <- o[,3] # num infected
}

storeMatrix.cgamma.R3 <- array(NA, dim=c(length(timeSeq),length(out_poscov_cgamma)))
for (j in 1:length(out_poscov_cgamma)) {
  o <- out_poscov_cgamma[[j]]
  storeMatrix.cgamma.R3[,j] <- o[,4] # num recovered
}

# PLOTS
par(mfrow=c(1,4))
## mean values

# NO COV
plot(0:150, apply(storeMatrix.cgamma.I, 1, mean), col="red", lwd=1.75, type="l", ylim=c(0,200), ylab="N", xlab="Time", main="No cov")
lines(0:150, apply(storeMatrix.cgamma.S, 1, mean), col="blue", lwd=1.75, type="l", ylim=c(0,150), ylab="N", xlab="Time", main="Susceptible")
lines(0:150, apply(storeMatrix.cgamma.R, 1, mean), col="green", lwd=1.75, type="l", ylim=c(0,150), ylab="N", xlab="Time", main="Recovered")
legend("topright",legend=c("S","I","R"), fill=c("blue","red","green"), cex=.5)

# NEG COV
plot(0:150, apply(storeMatrix.cgamma.I2, 1, mean), col="red", lwd=1.75, type="l", ylim=c(0,200), ylab="N", xlab="Time", main="Neg cov")
lines(0:150, apply(storeMatrix.cgamma.S2, 1, mean), col="blue", lwd=1.75, type="l", ylim=c(0,150), ylab="N", xlab="Time", main="Susceptible")
lines(0:150, apply(storeMatrix.cgamma.R2, 1, mean), col="green", lwd=1.75, type="l", ylim=c(0,150), ylab="N", xlab="Time", main="Recovered")
legend("topright",legend=c("S","I","R"), fill=c("blue","red","green"), cex=.5)

# POS COV
plot(0:150, apply(storeMatrix.cgamma.I3, 1, mean), col="red", lwd=1.75, type="l", ylim=c(0,200), ylab="N", xlab="Time", main="Pos cov")
lines(0:150, apply(storeMatrix.cgamma.S3, 1, mean), col="blue", lwd=1.75, type="l", ylim=c(0,150), ylab="N", xlab="Time", main="Susceptible")
lines(0:150, apply(storeMatrix.cgamma.R3, 1, mean), col="green", lwd=1.75, type="l", ylim=c(0,150), ylab="N", xlab="Time", main="Recovered")
legend("topright",legend=c("S","I","R"), fill=c("blue","red","green"), cex=.5)

# DETERMINISTIC
plot.ts(S, col="blue", ylim=c(0,200), ylab="N", xlab="Time", main= "Deterministic")
lines(I, col="red")
lines(R, col="green")
legend("topright",legend=c("S","I","R"),fill=c("blue","red","green"), cex=.7)

############################  CONTACT & ALPHA   #################################
############################   NO COVARIATION   #################################
storeMatrix.calpha.S <- array(NA, dim=c(length(timeSeq),length(out_nocov_calpha)))
for (j in 1:length(out_nocov_calpha)) {
  o <- out_nocov_calpha[[j]]
  storeMatrix.calpha.S[,j] <- o[,2] # num susceptible
}

storeMatrix.calpha.I <- array(NA, dim=c(length(timeSeq),length(out_nocov_calpha)))
for (j in 1:length(out_nocov_calpha)) {
  o <- out_nocov_calpha[[j]]
  storeMatrix.calpha.I[,j] <- o[,3] # num infected
}

storeMatrix.calpha.R <- array(NA, dim=c(length(timeSeq),length(out_nocov_calpha)))
for (j in 1:length(out_nocov_calpha)) {
  o <- out_nocov_calpha[[j]]
  storeMatrix.calpha.R[,j] <- o[,4] # num recovered
}

############################   NEG COVARIATION   #################################
storeMatrix.calpha.S2 <- array(NA, dim=c(length(timeSeq),length(out_negcov_calpha)))
for (j in 1:length(out_negcov_calpha)) {
  o <- out_negcov_calpha[[j]]
  storeMatrix.calpha.S2[,j] <- o[,2] # num susceptible
}

storeMatrix.calpha.I2 <- array(NA, dim=c(length(timeSeq),length(out_negcov_calpha)))
for (j in 1:length(out_negcov_calpha)) {
  o <- out_negcov_calpha[[j]]
  storeMatrix.calpha.I2[,j] <- o[,3] # num infected
}

storeMatrix.calpha.R2 <- array(NA, dim=c(length(timeSeq),length(out_negcov_calpha)))
for (j in 1:length(out_negcov_calpha)) {
  o <- out_negcov_calpha[[j]]
  storeMatrix.calpha.R2[,j] <- o[,4] # num recovered
}

############################   POS COVARIATION   #################################
storeMatrix.calpha.S3 <- array(NA, dim=c(length(timeSeq),length(out_poscov_calpha)))
for (j in 1:length(out_poscov_calpha)) {
  o <- out_poscov_calpha[[j]]
  storeMatrix.calpha.S3[,j] <- o[,2] # num susceptible
}

storeMatrix.calpha.I3 <- array(NA, dim=c(length(timeSeq),length(out_poscov_calpha)))
for (j in 1:length(out_poscov_calpha)) {
  o <- out_poscov_calpha[[j]]
  storeMatrix.calpha.I3[,j] <- o[,3] # num infected
}

storeMatrix.calpha.R3 <- array(NA, dim=c(length(timeSeq),length(out_poscov_calpha)))
for (j in 1:length(out_poscov_calpha)) {
  o <- out_poscov_calpha[[j]]
  storeMatrix.calpha.R3[,j] <- o[,4] # num recovered
}

# PLOTS
par(mfrow=c(1,4))
## mean values 

# NO COV
plot(0:150, apply(storeMatrix.calpha.I, 1, mean), col="red", lwd=1.75, type="l", ylim=c(0,200), ylab="N", xlab="Time", main="No cov")
lines(0:150, apply(storeMatrix.calpha.S, 1, mean), col="blue", lwd=1.75, type="l", ylim=c(0,150), ylab="N", xlab="Time", main="Susceptible")
lines(0:150, apply(storeMatrix.calpha.R, 1, mean), col="green", lwd=1.75, type="l", ylim=c(0,150), ylab="N", xlab="Time", main="Recovered")
legend("topright",legend=c("S","I","R"),fill=c("blue","red","green"), cex=0.5)

# NEG COV
plot(0:150, apply(storeMatrix.calpha.I2, 1, mean), col="red", lwd=1.75, type="l", ylim=c(0,200), ylab="N", xlab="Time", main="Neg cov")
lines(0:150, apply(storeMatrix.calpha.S2, 1, mean), col="blue", lwd=1.75, type="l", ylim=c(0,150), ylab="N", xlab="Time", main="Susceptible")
lines(0:150, apply(storeMatrix.calpha.R2, 1, mean), col="green", lwd=1.75, type="l", ylim=c(0,150), ylab="N", xlab="Time", main="Recovered")
legend("topright",legend=c("S","I","R"),fill=c("blue","red","green"), cex=0.5)

# POS COV
plot(0:150, apply(storeMatrix.calpha.I3, 1, mean), col="red", lwd=1.75, type="l", ylim=c(0,200), ylab="N", xlab="Time", main="Pos cov")
lines(0:150, apply(storeMatrix.calpha.S3, 1, mean), col="blue", lwd=1.75, type="l", ylim=c(0,150), ylab="N", xlab="Time", main="Susceptible")
lines(0:150, apply(storeMatrix.calpha.R3, 1, mean), col="green", lwd=1.75, type="l", ylim=c(0,150), ylab="N", xlab="Time", main="Recovered")
legend("topright",legend=c("S","I","R"),fill=c("blue","red","green"), cex=0.5)

# DETERMINISTIC
out=lsoda(xstart,times, deterministic.SIR, contact.alpha.params2)
S=out[,2]
I=out[,3]
R=out[,4]

plot.ts(S, col="blue", ylim=c(0,200), ylab="N", xlab="Time", main="Deterministic")
lines(I, col="red")
lines(R, col="green")
legend("topright",legend=c("S","I","R"),fill=c("blue","red","green"), cex=.5)


############################   CONTACT & TAU    #################################
############################   NO COVARIATION   #################################
storeMatrix.ctau.S <- array(NA, dim=c(length(timeSeq),length(out_nocov_ctau)))
for (j in 1:length(out_nocov_ctau)) {
  o <- out_nocov_ctau[[j]]
  storeMatrix.ctau.S[,j] <- o[,2] # num susceptible
}

storeMatrix.ctau.I <- array(NA, dim=c(length(timeSeq),length(out_nocov_ctau)))
for (j in 1:length(out_nocov_ctau)) {
  o <- out_nocov_ctau[[j]]
  storeMatrix.ctau.I[,j] <- o[,3] # num infected
}

storeMatrix.ctau.R <- array(NA, dim=c(length(timeSeq),length(out_nocov_ctau)))
for (j in 1:length(out_nocov_ctau)) {
  o <- out_nocov_ctau[[j]]
  storeMatrix.ctau.R[,j] <- o[,4] # num recovered
}

############################  NEG COVARIATION   #################################
storeMatrix.ctau.S2 <- array(NA, dim=c(length(timeSeq),length(out_negcov_ctau)))
for (j in 1:length(out_negcov_ctau)) {
  o <- out_negcov_ctau[[j]]
  storeMatrix.ctau.S2[,j] <- o[,2] # num susceptible
}

storeMatrix.ctau.I2 <- array(NA, dim=c(length(timeSeq),length(out_negcov_ctau)))
for (j in 1:length(out_negcov_ctau)) {
  o <- out_negcov_ctau[[j]]
  storeMatrix.ctau.I2[,j] <- o[,3] # num infected
}

storeMatrix.ctau.R2 <- array(NA, dim=c(length(timeSeq),length(out_negcov_ctau)))
for (j in 1:length(out_negcov_ctau)) {
  o <- out_negcov_ctau[[j]]
  storeMatrix.ctau.R2[,j] <- o[,4] # num recovered
}

############################  POS COVARIATION   #################################
storeMatrix.ctau.S3 <- array(NA, dim=c(length(timeSeq),length(out_poscov_ctau)))
for (j in 1:length(out_poscov_ctau)) {
  o <- out_poscov_ctau[[j]]
  storeMatrix.ctau.S3[,j] <- o[,2] # num susceptible
}

storeMatrix.ctau.I3 <- array(NA, dim=c(length(timeSeq),length(out_poscov_ctau)))
for (j in 1:length(out_poscov_ctau)) {
  o <- out_poscov_ctau[[j]]
  storeMatrix.ctau.I3[,j] <- o[,3] # num infected
}

storeMatrix.ctau.R3 <- array(NA, dim=c(length(timeSeq),length(out_poscov_ctau)))
for (j in 1:length(out_poscov_ctau)) {
  o <- out_poscov_ctau[[j]]
  storeMatrix.ctau.R3[,j] <- o[,4] # num recovered
}

# PLOTS
par(mfrow=c(1,4))
## mean values 

# NO COV
plot(0:150, apply(storeMatrix.ctau.I, 1, mean), col="red", lwd=1.75, type="l", ylim=c(0,250), ylab="N", xlab="Time", main="No cov")
lines(0:150, apply(storeMatrix.ctau.S, 1, mean), col="blue", lwd=1.75, type="l", ylim=c(0,150), ylab="N", xlab="Time", main="Susceptible")
lines(0:150, apply(storeMatrix.ctau.R, 1, mean), col="green", lwd=1.75, type="l", ylim=c(0,150), ylab="N", xlab="Time", main="Recovered")
legend("topright",legend=c("S","I","R"),fill=c("blue","red","green"), cex=0.5)

# NEG COV
plot(0:150, apply(storeMatrix.ctau.I2, 1, mean), col="red", lwd=1.75, type="l", ylim=c(0,250), ylab="N", xlab="Time", main="Neg cov")
lines(0:150, apply(storeMatrix.ctau.S2, 1, mean), col="blue", lwd=1.75, type="l", ylim=c(0,150), ylab="N", xlab="Time", main="Susceptible")
lines(0:150, apply(storeMatrix.ctau.R2, 1, mean), col="green", lwd=1.75, type="l", ylim=c(0,150), ylab="N", xlab="Time", main="Recovered")
legend("topright",legend=c("S","I","R"),fill=c("blue","red","green"), cex=0.5)

# POS COV
plot(0:150, apply(storeMatrix.ctau.I3, 1, mean), col="red", lwd=1.75, type="l", ylim=c(0,250), ylab="N", xlab="Time", main="Pos cov")
lines(0:150, apply(storeMatrix.ctau.S3, 1, mean), col="blue", lwd=1.75, type="l", ylim=c(0,150), ylab="N", xlab="Time", main="Susceptible")
lines(0:150, apply(storeMatrix.ctau.R3, 1, mean), col="green", lwd=1.75, type="l", ylim=c(0,150), ylab="N", xlab="Time", main="Recovered")
legend("topright",legend=c("S","I","R"),fill=c("blue","red","green"), cex=0.5)

# DETERMINISTIC
out=lsoda(xstart,times, deterministic.SIR, contact.tau.params2)
S=out[,2]
I=out[,3]
R=out[,4]

plot.ts(S, col="blue", ylim=c(0,250), ylab="N", xlab="Time", main="Deterministic")
lines(I, col="red")
lines(R, col="green")
legend("topright",legend=c("S","I","R"),fill=c("blue","red","green"), cex=.5)



