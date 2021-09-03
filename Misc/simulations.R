### simulations and plots for contact and shedding gem sir models

# packages
library(tidyverse)
library(parallel)

# seeds and parameter values
seeds <- floor(runif(20,1,1e5)) # set seeds
x = c(S=70, I=10, R=0)
tmax <- 150

params = c(c=.2, shed=.2, sd_s=.01, sd_a=.01, sd_c=.01, sd_g=.01, 
           h=.1, alpha=.01, gamma=.3, b=2.5, d=.4, bs=.01)
params2 = c(c=.2, shed=.2, sd_s=.2, sd_a=.2, sd_c=.2, sd_g=.2, 
           h=.1, alpha=.01, gamma=.3, b=2.5, d=.4, bs=.01)
nocorr <- matrix(c(1,0,0,1), nrow=2, byrow=T)
negcorr <- matrix(c(1,-.5,-.5,1), nrow=2, byrow=T)
poscorr <- matrix(c(1,.5,.5,1), nrow=2, byrow=T)

#  model simulations
## no covariance
mclapply(seeds,
         function(s) gillespie.SIR.cov_taualpha(tmax, params, nocorr, x),
         mc.cores=4) -> out_nocov_taualpha

mclapply(seeds,
         function(s) gillespie.SIR.cov_taugamma(tmax, params, nocorr, x),
         mc.cores=4) -> out_nocov_taugamma

mclapply(seeds,
         function(s) gillespie.SIR.cov_cgamma(tmax, params, nocorr, x),
         mc.cores=4) -> out_nocov_cgamma

mclapply(seeds,
         function(s) gillespie.SIR.cov_calpha(tmax, params, nocorr, x),
         mc.cores=4) -> out_nocov_calpha

mclapply(seeds,
         function(s) gillespie.SIR.cov_ctau(tmax, params, nocorr, x),
         mc.cores=4) -> out_nocov_ctau

## negative covariance
mclapply(seeds,
         function(s) gillespie.SIR.cov_taualpha(tmax, params2, negcorr, x),
         mc.cores=4) -> out_negcov_taualpha2

mclapply(seeds,
         function(s) gillespie.SIR.cov_taugamma(tmax, params, negcorr, x),
         mc.cores=4) -> out_negcov_taugamma

mclapply(seeds,
         function(s) gillespie.SIR.cov_cgamma(tmax, params, negcorr, x),
         mc.cores=4) -> out_negcov_cgamma

mclapply(seeds,
         function(s) gillespie.SIR.cov_calpha(tmax, params, negcorr, x),
         mc.cores=4) -> out_negcov_calpha

mclapply(seeds,
         function(s) gillespie.SIR.cov_ctau(tmax, params, negcorr, x),
         mc.cores=4) -> out_negcov_ctau

## positive cov
mclapply(seeds,
         function(s) gillespie.SIR.cov_taualpha(tmax, params2, poscorr, x),
         mc.cores=4) -> out_poscov_taualpha2

mclapply(seeds,
         function(s) gillespie.SIR.cov_taugamma(tmax, params, poscorr, x),
         mc.cores=4) -> out_poscov_taugamma

mclapply(seeds,
         function(s) gillespie.SIR.cov_cgamma(tmax, params, poscorr, x),
         mc.cores=4) -> out_poscov_cgamma

mclapply(seeds,
         function(s) gillespie.SIR.cov_calpha(tmax, params, poscorr, x),
         mc.cores=4) -> out_poscov_calpha

mclapply(seeds,
         function(s) gillespie.SIR.cov_ctau(tmax, params, poscorr, x),
         mc.cores=4) -> out_poscov_ctau


## get S, I, & R for each seed
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

# plot
par(mfrow=c(1,3))

##infected
plot(0:150, apply(storeMatrix.taualpha.I, 1, mean), col="red", lwd=1.25, type="l", ylim=c(0,100), ylab="N", xlab="Time", main="Infected")
for (i in 1:ncol(storeMatrix.taualpha.I)) {
  lines(0:150, storeMatrix.taualpha.I[,i], col=grey(.1, alpha=.15))
}
lines(0:150, apply(storeMatrix.taualpha.I, 1, mean), col="red", lwd=1.25, type="l", ylim=c(0,100), ylab="N", xlab="Time", main="Infected")
lines(lowerbnd.I, lwd=1,lty=1, col="black")
lines(upperbnd.I, lwd=1,lty=1,col="black")

lowerbnd.I = apply(storeMatrix.taualpha.I,1,function(x){quantile(x,c(.025,.975))})[1,]
upperbnd.I = apply(storeMatrix.taualpha.I,1,function(x){quantile(x,c(.025,.975))})[2,]

## susceptible
plot(0:150, apply(storeMatrix.taualpha.S, 1, mean), col="red", lwd=1.25, type="l", ylim=c(0,100), ylab="No cov", xlab="Time", main="Susceptible")
for (i in 1:ncol(storeMatrix.taualpha.S)) {
  lines(0:150, storeMatrix.taualpha.S[,i], col=grey(.1, alpha=.15))
}
lines(0:150, apply(storeMatrix.taualpha.S, 1, mean), col="red", lwd=1.25, type="l", ylim=c(0,100), ylab="N", xlab="Time", main="Susceptible")
lines(lowerbnd.S, lwd=1,lty=1, col="black")
lines(upperbnd.S, lwd=1,lty=1,col="black")

## recovered
plot(0:150, apply(storeMatrix.taualpha.R, 1, mean), col="red", lwd=1.25, type="l", ylim=c(0,100), ylab="N", xlab="Time", main="Recovered")
for (i in 1:ncol(storeMatrix.taualpha.R)) {
  lines(0:150, storeMatrix.taualpha.R[,i], col=grey(.1, alpha=.15))
}
lines(0:150, apply(storeMatrix.taualpha.R, 1, mean), col="red", lwd=1.25, type="l", ylim=c(0,100), ylab="N", xlab="Time", main="Recovered")
lines(lowerbnd.R, lwd=1,lty=1, col="black")
lines(upperbnd.R, lwd=1,lty=1,col="black")

lowerbnd.R = apply(storeMatrix.taualpha.R,1,function(x){quantile(x,c(.025,.975))})[1,]
upperbnd.R = apply(storeMatrix.taualpha.R,1,function(x){quantile(x,c(.025,.975))})[2,]

par(mfrow=c(1,1))
## mean values 
plot(0:150, apply(storeMatrix.taualpha.I, 1, mean), col="red", lwd=1.75, type="l", ylim=c(0,150), ylab="N", xlab="Time", main="Simulation means")
lines(0:150, apply(storeMatrix.taualpha.S, 1, mean), col="blue", lwd=1.75, type="l", ylim=c(0,150), ylab="N", xlab="Time", main="Susceptible")
lines(0:150, apply(storeMatrix.taualpha.R, 1, mean), col="green", lwd=1.75, type="l", ylim=c(0,150), ylab="N", xlab="Time", main="Recovered")
legend("topright",legend=c("S","I","R"),fill=c("blue","red","green"), cex=0.75)


#### tau gamma non cov
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

## plot
par(mfrow=c(1,3))

## susceptible
plot(0:150, apply(storeMatrix.taugamma.S, 1, mean), col="red", lwd=1.25, type="l", ylim=c(0,200), ylab="No cov", xlab="Time", main="Susceptible")
for (i in 1:ncol(storeMatrix.taugamma.S)) {
  lines(0:150, storeMatrix.taugamma.S[,i], col=grey(.1, alpha=.15))
}
lines(0:150, apply(storeMatrix.taugamma.S, 1, mean), col="red", lwd=1.25, type="l", ylim=c(0,200), ylab="N", xlab="Time", main="Susceptible")
lines(lowerbnd.Stg, lwd=1,lty=1, col="black")
lines(upperbnd.Stg, lwd=1,lty=1,col="black")

lowerbnd.Stg = apply(storeMatrix.taugamma.S,1,function(x){quantile(x,c(.025,.975))})[1,]
upperbnd.Stg = apply(storeMatrix.taugamma.S,1,function(x){quantile(x,c(.025,.975))})[2,]

## infected
plot(0:150, apply(storeMatrix.taugamma.I, 1, mean), col="red", lwd=1.25, type="l", ylim=c(0,200), ylab="N", xlab="Time", main="Infected")
for (i in 1:ncol(storeMatrix.taugamma.I)) {
  lines(0:150, storeMatrix.taugamma.I[,i], col=grey(.1, alpha=.15))
}
lines(0:150, apply(storeMatrix.taugamma.I, 1, mean), col="red", lwd=1.25, type="l", ylim=c(0,200), ylab="N", xlab="Time", main="Infected")
lines(lowerbnd.Itg, lwd=1,lty=1, col="black")
lines(upperbnd.Itg, lwd=1,lty=1,col="black")

lowerbnd.Itg = apply(storeMatrix.taugamma.I,1,function(x){quantile(x,c(.025,.975))})[1,]
upperbnd.Itg = apply(storeMatrix.taugamma.I,1,function(x){quantile(x,c(.025,.975))})[2,]

## recovered
plot(0:150, apply(storeMatrix.taugamma.R, 1, mean), col="red", lwd=1.25, type="l", ylim=c(0,200), ylab="N", xlab="Time", main="Recovered")
for (i in 1:ncol(storeMatrix.taugamma.R)) {
  lines(0:150, storeMatrix.taugamma.R[,i], col=grey(.1, alpha=.15))
}
lines(0:150, apply(storeMatrix.taugamma.R, 1, mean), col="red", lwd=1.25, type="l", ylim=c(0,200), ylab="N", xlab="Time", main="Recovered")
lines(lowerbnd.Rtg, lwd=1,lty=1, col="black")
lines(upperbnd.Rtg, lwd=1,lty=1,col="black")

lowerbnd.Rtg = apply(storeMatrix.taugamma.R,1,function(x){quantile(x,c(.025,.975))})[1,]
upperbnd.Rtg = apply(storeMatrix.taugamma.R,1,function(x){quantile(x,c(.025,.975))})[2,]

par(mfrow=c(1,1))
## mean values 
plot(0:150, apply(storeMatrix.taugamma.I, 1, mean), col="red", lwd=1.75, type="l", ylim=c(0,200), ylab="N", xlab="Time", main="Simulation means")
lines(0:150, apply(storeMatrix.taugamma.S, 1, mean), col="blue", lwd=1.75, type="l", ylim=c(0,150), ylab="N", xlab="Time", main="Susceptible")
lines(0:150, apply(storeMatrix.taugamma.R, 1, mean), col="green", lwd=1.75, type="l", ylim=c(0,150), ylab="N", xlab="Time", main="Recovered")
legend("topright",legend=c("S","I","R"),fill=c("blue","red","green"), cex=0.75)


### no cov contact and gamma

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

## plot
par(mfrow=c(1,3))

## susceptible
lowerbnd.Scg = apply(storeMatrix.cgamma.S,1,function(x){quantile(x,c(.025,.975))})[1,]
upperbnd.Scg = apply(storeMatrix.cgamma.S,1,function(x){quantile(x,c(.025,.975))})[2,]

plot(0:150, apply(storeMatrix.cgamma.S, 1, mean), col="red", lwd=1.25, type="l", ylim=c(0,200), ylab="No cov", xlab="Time", main="Susceptible")
for (i in 1:ncol(storeMatrix.cgamma.S)) {
  lines(0:150, storeMatrix.cgamma.S[,i], col=grey(.1, alpha=.15))
}
lines(0:150, apply(storeMatrix.cgamma.S, 1, mean), col="red", lwd=1.25, type="l", ylim=c(0,200), ylab="N", xlab="Time", main="Susceptible")
lines(lowerbnd.Scg, lwd=1,lty=1, col="black")
lines(upperbnd.Scg, lwd=1,lty=1,col="black")

## infected
lowerbnd.Icg = apply(storeMatrix.cgamma.I,1,function(x){quantile(x,c(.025,.975))})[1,]
upperbnd.Icg = apply(storeMatrix.cgamma.I,1,function(x){quantile(x,c(.025,.975))})[2,]

plot(0:150, apply(storeMatrix.cgamma.I, 1, mean), col="red", lwd=1.25, type="l", ylim=c(0,200), ylab="N", xlab="Time", main="Infected")
for (i in 1:ncol(storeMatrix.cgamma.I)) {
  lines(0:150, storeMatrix.cgamma.I[,i], col=grey(.1, alpha=.15))
}
lines(0:150, apply(storeMatrix.cgamma.I, 1, mean), col="red", lwd=1.25, type="l", ylim=c(0,200), ylab="N", xlab="Time", main="Infected")
lines(lowerbnd.Icg, lwd=1,lty=1, col="black")
lines(upperbnd.Icg, lwd=1,lty=1,col="black")

## recovered
lowerbnd.Rcg = apply(storeMatrix.cgamma.R,1,function(x){quantile(x,c(.025,.975))})[1,]
upperbnd.Rcg = apply(storeMatrix.cgamma.R,1,function(x){quantile(x,c(.025,.975))})[2,]

plot(0:150, apply(storeMatrix.cgamma.R, 1, mean), col="red", lwd=1.25, type="l", ylim=c(0,200), ylab="N", xlab="Time", main="Recovered")
for (i in 1:ncol(storeMatrix.cgamma.R)) {
  lines(0:150, storeMatrix.cgamma.R[,i], col=grey(.1, alpha=.15))
}
lines(0:150, apply(storeMatrix.cgamma.R, 1, mean), col="red", lwd=1.25, type="l", ylim=c(0,200), ylab="N", xlab="Time", main="Recovered")
lines(lowerbnd.Rcg, lwd=1,lty=1, col="black")
lines(upperbnd.Rcg, lwd=1,lty=1,col="black")

par(mfrow=c(1,1))
## mean values 
plot(0:150, apply(storeMatrix.cgamma.I, 1, mean), col="red", lwd=1.75, type="l", ylim=c(0,200), ylab="N", xlab="Time", main="Simulation means")
lines(0:150, apply(storeMatrix.cgamma.S, 1, mean), col="blue", lwd=1.75, type="l", ylim=c(0,150), ylab="N", xlab="Time", main="Susceptible")
lines(0:150, apply(storeMatrix.cgamma.R, 1, mean), col="green", lwd=1.75, type="l", ylim=c(0,150), ylab="N", xlab="Time", main="Recovered")
legend("topright",legend=c("S","I","R"),fill=c("blue","red","green"), cex=0.75)

### no cov contact and alpha

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

## plot
par(mfrow=c(1,3))

## susceptible
lowerbnd.Sca = apply(storeMatrix.calpha.S,1,function(x){quantile(x,c(.025,.975))})[1,]
upperbnd.Sca = apply(storeMatrix.calpha.S,1,function(x){quantile(x,c(.025,.975))})[2,]

plot(0:150, apply(storeMatrix.calpha.S, 1, mean), col="red", lwd=1.25, type="l", ylim=c(0,200), ylab="No cov", xlab="Time", main="Susceptible")
for (i in 1:ncol(storeMatrix.calpha.S)) {
  lines(0:150, storeMatrix.calpha.S[,i], col=grey(.1, alpha=.15))
}
lines(0:150, apply(storeMatrix.calpha.S, 1, mean), col="red", lwd=1.25, type="l", ylim=c(0,200), ylab="N", xlab="Time", main="Susceptible")
lines(lowerbnd.Sca, lwd=1,lty=1, col="black")
lines(upperbnd.Sca, lwd=1,lty=1,col="black")

## infected
lowerbnd.Ica = apply(storeMatrix.calpha.I,1,function(x){quantile(x,c(.025,.975))})[1,]
upperbnd.Ica = apply(storeMatrix.calpha.I,1,function(x){quantile(x,c(.025,.975))})[2,]

plot(0:150, apply(storeMatrix.calpha.I, 1, mean), col="red", lwd=1.25, type="l", ylim=c(0,200), ylab="N", xlab="Time", main="Infected")
for (i in 1:ncol(storeMatrix.calpha.I)) {
  lines(0:150, storeMatrix.calpha.I[,i], col=grey(.1, alpha=.15))
}
lines(0:150, apply(storeMatrix.calpha.I, 1, mean), col="red", lwd=1.25, type="l", ylim=c(0,200), ylab="N", xlab="Time", main="Infected")
lines(lowerbnd.Ica, lwd=1,lty=1, col="black")
lines(upperbnd.Ica, lwd=1,lty=1,col="black")

## recovered
lowerbnd.Rca = apply(storeMatrix.calpha.R,1,function(x){quantile(x,c(.025,.975))})[1,]
upperbnd.Rca = apply(storeMatrix.calpha.R,1,function(x){quantile(x,c(.025,.975))})[2,]

plot(0:150, apply(storeMatrix.calpha.R, 1, mean), col="red", lwd=1.25, type="l", ylim=c(0,200), ylab="N", xlab="Time", main="Recovered")
for (i in 1:ncol(storeMatrix.calpha.R)) {
  lines(0:150, storeMatrix.calpha.R[,i], col=grey(.1, alpha=.15))
}
lines(0:150, apply(storeMatrix.calpha.R, 1, mean), col="red", lwd=1.25, type="l", ylim=c(0,200), ylab="N", xlab="Time", main="Recovered")
lines(lowerbnd.Rca, lwd=1,lty=1, col="black")
lines(upperbnd.Rca, lwd=1,lty=1,col="black")

par(mfrow=c(1,1))
## mean values 
plot(0:150, apply(storeMatrix.calpha.I, 1, mean), col="red", lwd=1.75, type="l", ylim=c(0,200), ylab="N", xlab="Time", main="Simulation means")
lines(0:150, apply(storeMatrix.calpha.S, 1, mean), col="blue", lwd=1.75, type="l", ylim=c(0,150), ylab="N", xlab="Time", main="Susceptible")
lines(0:150, apply(storeMatrix.calpha.R, 1, mean), col="green", lwd=1.75, type="l", ylim=c(0,150), ylab="N", xlab="Time", main="Recovered")
legend("topright",legend=c("S","I","R"),fill=c("blue","red","green"), cex=0.75)

### no cov contact and transmisibility

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

## plot
par(mfrow=c(1,3))

## susceptible
lowerbnd.Sct = apply(storeMatrix.ctau.S,1,function(x){quantile(x,c(.025,.975))})[1,]
upperbnd.Sct = apply(storeMatrix.ctau.S,1,function(x){quantile(x,c(.025,.975))})[2,]

plot(0:150, apply(storeMatrix.ctau.S, 1, mean), col="red", lwd=1.25, type="l", ylim=c(0,200), ylab="No cov", xlab="Time", main="Susceptible")
for (i in 1:ncol(storeMatrix.ctau.S)) {
  lines(0:150, storeMatrix.ctau.S[,i], col=grey(.1, alpha=.15))
}
lines(0:150, apply(storeMatrix.ctau.S, 1, mean), col="red", lwd=1.25, type="l", ylim=c(0,200), ylab="N", xlab="Time", main="Susceptible")
lines(lowerbnd.Sct, lwd=1,lty=1, col="black")
lines(upperbnd.Sct, lwd=1,lty=1,col="black")

## infected
lowerbnd.Ict = apply(storeMatrix.ctau.I,1,function(x){quantile(x,c(.025,.975))})[1,]
upperbnd.Ict = apply(storeMatrix.ctau.I,1,function(x){quantile(x,c(.025,.975))})[2,]

plot(0:150, apply(storeMatrix.ctau.I, 1, mean), col="red", lwd=1.25, type="l", ylim=c(0,200), ylab="N", xlab="Time", main="Infected")
for (i in 1:ncol(storeMatrix.ctau.I)) {
  lines(0:150, storeMatrix.ctau.I[,i], col=grey(.1, alpha=.15))
}
lines(0:150, apply(storeMatrix.ctau.I, 1, mean), col="red", lwd=1.25, type="l", ylim=c(0,200), ylab="N", xlab="Time", main="Infected")
lines(lowerbnd.Ict, lwd=1,lty=1, col="black")
lines(upperbnd.Ict, lwd=1,lty=1,col="black")

## recovered
lowerbnd.Rct = apply(storeMatrix.ctau.R,1,function(x){quantile(x,c(.025,.975))})[1,]
upperbnd.Rct = apply(storeMatrix.ctau.R,1,function(x){quantile(x,c(.025,.975))})[2,]

plot(0:150, apply(storeMatrix.ctau.R, 1, mean), col="red", lwd=1.25, type="l", ylim=c(0,200), ylab="N", xlab="Time", main="Recovered")
for (i in 1:ncol(storeMatrix.ctau.R)) {
  lines(0:150, storeMatrix.ctau.R[,i], col=grey(.1, alpha=.15))
}
lines(0:150, apply(storeMatrix.ctau.R, 1, mean), col="red", lwd=1.25, type="l", ylim=c(0,200), ylab="N", xlab="Time", main="Recovered")
lines(lowerbnd.Rct, lwd=1,lty=1, col="black")
lines(upperbnd.Rct, lwd=1,lty=1,col="black")

par(mfrow=c(1,1))
## mean values 
plot(0:150, apply(storeMatrix.ctau.I, 1, mean), col="red", lwd=1.75, type="l", ylim=c(0,200), ylab="N", xlab="Time", main="Simulation means")
lines(0:150, apply(storeMatrix.ctau.S, 1, mean), col="blue", lwd=1.75, type="l", ylim=c(0,150), ylab="N", xlab="Time", main="Susceptible")
lines(0:150, apply(storeMatrix.ctau.R, 1, mean), col="green", lwd=1.75, type="l", ylim=c(0,150), ylab="N", xlab="Time", main="Recovered")
legend("topright",legend=c("S","I","R"),fill=c("blue","red","green"), cex=0.75)

########### NEG COVARIANCE PLOTS #####################
timeSeq <- 0:250
storeMatrix.taualpha.I2 <- array(NA, dim=c(length(timeSeq),length(out_negcov_taualpha2)))
for (j in 1:length(out_negcov_taualpha2)) {
  o <- out_negcov_taualpha2[[j]]
  storeMatrix.taualpha.I2[,j] <- o[,3] # num infected
}

## susceptible
storeMatrix.taualpha.S2 <- array(NA, dim=c(length(timeSeq),length(out_negcov_taualpha2)))
for (j in 1:length(out_negcov_taualpha2)) {
  o <- out_negcov_taualpha2[[j]]
  storeMatrix.taualpha.S2[,j] <- o[,2] # num infected
}

## recovered
storeMatrix.taualpha.R2 <- array(NA, dim=c(length(timeSeq),length(out_negcov_taualpha2)))
for (j in 1:length(out_negcov_taualpha2)) {
  o <- out_negcov_taualpha2[[j]]
  storeMatrix.taualpha.R2[,j] <- o[,4] # num recovered
}

# plot
par(mfrow=c(1,3))

##infected
lowerbnd.I2 = apply(storeMatrix.taualpha.I2,1,function(x){quantile(x,c(.025,.975))})[1,]
upperbnd.I2 = apply(storeMatrix.taualpha.I2,1,function(x){quantile(x,c(.025,.975))})[2,]

plot(0:150, apply(storeMatrix.taualpha.I2, 1, mean), col="red", lwd=1.25, type="l", ylim=c(0,100), ylab="N", xlab="Time", main="Infected")
for (i in 1:ncol(storeMatrix.taualpha.I2)) {
  lines(0:150, storeMatrix.taualpha.I2[,i], col=grey(.1, alpha=.15))
}
lines(0:150, apply(storeMatrix.taualpha.I2, 1, mean), col="red", lwd=1.25, type="l", ylim=c(0,100), ylab="N", xlab="Time", main="Infected")
lines(lowerbnd.I2, lwd=1,lty=1, col="black")
lines(upperbnd.I2, lwd=1,lty=1,col="black")

## susceptible
lowerbnd.S2 = apply(storeMatrix.taualpha.S2,1,function(x){quantile(x,c(.025,.975))})[1,]
upperbnd.S2 = apply(storeMatrix.taualpha.S2,1,function(x){quantile(x,c(.025,.975))})[2,]

plot(0:150, apply(storeMatrix.taualpha.S2, 1, mean), col="red", lwd=1.25, type="l", ylim=c(0,100), ylab="Neg cov", xlab="Time", main="Susceptible")
for (i in 1:ncol(storeMatrix.taualpha.S2)) {
  lines(0:150, storeMatrix.taualpha.S2[,i], col=grey(.1, alpha=.15))
}
lines(0:150, apply(storeMatrix.taualpha.S2, 1, mean), col="red", lwd=1.25, type="l", ylim=c(0,100), ylab="N", xlab="Time", main="Susceptible")
lines(lowerbnd.S2, lwd=1,lty=1, col="black")
lines(upperbnd.S2, lwd=1,lty=1,col="black")

## recovered
lowerbnd.R2 = apply(storeMatrix.taualpha.R2,1,function(x){quantile(x,c(.025,.975))})[1,]
upperbnd.R2 = apply(storeMatrix.taualpha.R2,1,function(x){quantile(x,c(.025,.975))})[2,]

plot(0:150, apply(storeMatrix.taualpha.R2, 1, mean), col="red", lwd=1.25, type="l", ylim=c(0,100), ylab="N", xlab="Time", main="Recovered")
for (i in 1:ncol(storeMatrix.taualpha.R2)) {
  lines(0:150, storeMatrix.taualpha.R2[,i], col=grey(.1, alpha=.15))
}
lines(0:150, apply(storeMatrix.taualpha.R2, 1, mean), col="red", lwd=1.25, type="l", ylim=c(0,100), ylab="N", xlab="Time", main="Recovered")
lines(lowerbnd.R2, lwd=1,lty=1, col="black")
lines(upperbnd.R2, lwd=1,lty=1,col="black")

par(mfrow=c(1,1))
## mean values 
plot(0:250, apply(storeMatrix.taualpha.I2, 1, mean), col="red", lwd=1.75, type="l", ylim=c(0,150), ylab="N", xlab="Time", main="Simulation means")
lines(0:250, apply(storeMatrix.taualpha.S2, 1, mean), col="blue", lwd=1.75, type="l", ylim=c(0,150), ylab="N", xlab="Time", main="Susceptible")
lines(0:250, apply(storeMatrix.taualpha.R2, 1, mean), col="green", lwd=1.75, type="l", ylim=c(0,150), ylab="N", xlab="Time", main="Recovered")
legend("topright",legend=c("S","I","R"),fill=c("blue","red","green"), cex=0.75)

#### tau gamma neg cov
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

## plot
par(mfrow=c(1,3))

## susceptible
lowerbnd.Stg2 = apply(storeMatrix.taugamma.S2,1,function(x){quantile(x,c(.025,.975))})[1,]
upperbnd.Stg2 = apply(storeMatrix.taugamma.S2,1,function(x){quantile(x,c(.025,.975))})[2,]

plot(0:150, apply(storeMatrix.taugamma.S2, 1, mean), col="red", lwd=1.25, type="l", ylim=c(0,200), ylab="Neg cov", xlab="Time", main="Susceptible")
for (i in 1:ncol(storeMatrix.taugamma.S2)) {
  lines(0:150, storeMatrix.taugamma.S2[,i], col=grey(.1, alpha=.15))
}
lines(0:150, apply(storeMatrix.taugamma.S2, 1, mean), col="red", lwd=1.25, type="l", ylim=c(0,200), ylab="N", xlab="Time", main="Susceptible")
lines(lowerbnd.Stg2, lwd=1,lty=1, col="black")
lines(upperbnd.Stg2, lwd=1,lty=1,col="black")

## infected
lowerbnd.Itg2 = apply(storeMatrix.taugamma.I2,1,function(x){quantile(x,c(.025,.975))})[1,]
upperbnd.Itg2 = apply(storeMatrix.taugamma.I2,1,function(x){quantile(x,c(.025,.975))})[2,]

plot(0:150, apply(storeMatrix.taugamma.I2, 1, mean), col="red", lwd=1.25, type="l", ylim=c(0,200), ylab="N", xlab="Time", main="Infected")
for (i in 1:ncol(storeMatrix.taugamma.I2)) {
  lines(0:150, storeMatrix.taugamma.I2[,i], col=grey(.1, alpha=.15))
}
lines(0:150, apply(storeMatrix.taugamma.I2, 1, mean), col="red", lwd=1.25, type="l", ylim=c(0,200), ylab="N", xlab="Time", main="Infected")
lines(lowerbnd.Itg2, lwd=1,lty=1, col="black")
lines(upperbnd.Itg2, lwd=1,lty=1,col="black")

## recovered
lowerbnd.Rtg2 = apply(storeMatrix.taugamma.R2,1,function(x){quantile(x,c(.025,.975))})[1,]
upperbnd.Rtg2 = apply(storeMatrix.taugamma.R2,1,function(x){quantile(x,c(.025,.975))})[2,]

plot(0:150, apply(storeMatrix.taugamma.R2, 1, mean), col="red", lwd=1.25, type="l", ylim=c(0,200), ylab="N", xlab="Time", main="Recovered")
for (i in 1:ncol(storeMatrix.taugamma.R2)) {
  lines(0:150, storeMatrix.taugamma.R2[,i], col=grey(.1, alpha=.15))
}
lines(0:150, apply(storeMatrix.taugamma.R2, 1, mean), col="red", lwd=1.25, type="l", ylim=c(0,200), ylab="N", xlab="Time", main="Recovered")
lines(lowerbnd.Rtg, lwd=1,lty=1, col="black")
lines(upperbnd.Rtg, lwd=1,lty=1,col="black")

par(mfrow=c(1,1))
## mean values 
plot(0:150, apply(storeMatrix.taugamma.I2, 1, mean), col="red", lwd=1.75, type="l", ylim=c(0,200), ylab="N", xlab="Time", main="Simulation means")
lines(0:150, apply(storeMatrix.taugamma.S2, 1, mean), col="blue", lwd=1.75, type="l", ylim=c(0,150), ylab="N", xlab="Time", main="Susceptible")
lines(0:150, apply(storeMatrix.taugamma.R2, 1, mean), col="green", lwd=1.75, type="l", ylim=c(0,150), ylab="N", xlab="Time", main="Recovered")
legend("topright",legend=c("S","I","R"),fill=c("blue","red","green"), cex=0.75)

### neg cov contact and gamma

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

## plot
par(mfrow=c(1,3))

## susceptible
lowerbnd.Scg2 = apply(storeMatrix.cgamma.S2,1,function(x){quantile(x,c(.025,.975))})[1,]
upperbnd.Scg2 = apply(storeMatrix.cgamma.S2,1,function(x){quantile(x,c(.025,.975))})[2,]

plot(0:150, apply(storeMatrix.cgamma.S2, 1, mean), col="red", lwd=1.25, type="l", ylim=c(0,200), ylab="Neg cov", xlab="Time", main="Susceptible")
for (i in 1:ncol(storeMatrix.cgamma.S2)) {
  lines(0:150, storeMatrix.cgamma.S2[,i], col=grey(.1, alpha=.15))
}
lines(0:150, apply(storeMatrix.cgamma.S2, 1, mean), col="red", lwd=1.25, type="l", ylim=c(0,200), ylab="N", xlab="Time", main="Susceptible")
lines(lowerbnd.Scg2, lwd=1,lty=1, col="black")
lines(upperbnd.Scg2, lwd=1,lty=1,col="black")

## infected
lowerbnd.Icg2 = apply(storeMatrix.cgamma.I2,1,function(x){quantile(x,c(.025,.975))})[1,]
upperbnd.Icg2 = apply(storeMatrix.cgamma.I2,1,function(x){quantile(x,c(.025,.975))})[2,]

plot(0:150, apply(storeMatrix.cgamma.I2, 1, mean), col="red", lwd=1.25, type="l", ylim=c(0,200), ylab="N", xlab="Time", main="Infected")
for (i in 1:ncol(storeMatrix.cgamma.I2)) {
  lines(0:150, storeMatrix.cgamma.I2[,i], col=grey(.1, alpha=.15))
}
lines(0:150, apply(storeMatrix.cgamma.I2, 1, mean), col="red", lwd=1.25, type="l", ylim=c(0,200), ylab="N", xlab="Time", main="Infected")
lines(lowerbnd.Icg2, lwd=1,lty=1, col="black")
lines(upperbnd.Icg2, lwd=1,lty=1,col="black")

## recovered
lowerbnd.Rcg2 = apply(storeMatrix.cgamma.R2,1,function(x){quantile(x,c(.025,.975))})[1,]
upperbnd.Rcg2 = apply(storeMatrix.cgamma.R2,1,function(x){quantile(x,c(.025,.975))})[2,]

plot(0:150, apply(storeMatrix.cgamma.R2, 1, mean), col="red", lwd=1.25, type="l", ylim=c(0,200), ylab="N", xlab="Time", main="Recovered")
for (i in 1:ncol(storeMatrix.cgamma.R2)) {
  lines(0:150, storeMatrix.cgamma.R2[,i], col=grey(.1, alpha=.15))
}
lines(0:150, apply(storeMatrix.cgamma.R2, 1, mean), col="red", lwd=1.25, type="l", ylim=c(0,200), ylab="N", xlab="Time", main="Recovered")
lines(lowerbnd.Rcg2, lwd=1,lty=1, col="black")
lines(upperbnd.Rcg2, lwd=1,lty=1,col="black")

par(mfrow=c(1,1))
## mean values 
plot(0:150, apply(storeMatrix.cgamma.I2, 1, mean), col="red", lwd=1.75, type="l", ylim=c(0,200), ylab="N", xlab="Time", main="Simulation means")
lines(0:150, apply(storeMatrix.cgamma.S2, 1, mean), col="blue", lwd=1.75, type="l", ylim=c(0,150), ylab="N", xlab="Time", main="Susceptible")
lines(0:150, apply(storeMatrix.cgamma.R2, 1, mean), col="green", lwd=1.75, type="l", ylim=c(0,150), ylab="N", xlab="Time", main="Recovered")
legend("topright",legend=c("S","I","R"),fill=c("blue","red","green"), cex=0.75)

### neg cov contact and alpha

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

## plot
par(mfrow=c(1,3))

## susceptible
lowerbnd.Sca2 = apply(storeMatrix.calpha.S2,1,function(x){quantile(x,c(.025,.975))})[1,]
upperbnd.Sca2 = apply(storeMatrix.calpha.S2,1,function(x){quantile(x,c(.025,.975))})[2,]

plot(0:150, apply(storeMatrix.calpha.S2, 1, mean), col="red", lwd=1.25, type="l", ylim=c(0,200), ylab="Neg cov", xlab="Time", main="Susceptible")
for (i in 1:ncol(storeMatrix.calpha.S2)) {
  lines(0:150, storeMatrix.calpha.S2[,i], col=grey(.1, alpha=.15))
}
lines(0:150, apply(storeMatrix.calpha.S2, 1, mean), col="red", lwd=1.25, type="l", ylim=c(0,200), ylab="N", xlab="Time", main="Susceptible")
lines(lowerbnd.Sca2, lwd=1,lty=1, col="black")
lines(upperbnd.Sca2, lwd=1,lty=1,col="black")

## infected
lowerbnd.Ica2 = apply(storeMatrix.calpha.I2,1,function(x){quantile(x,c(.025,.975))})[1,]
upperbnd.Ica2 = apply(storeMatrix.calpha.I2,1,function(x){quantile(x,c(.025,.975))})[2,]

plot(0:150, apply(storeMatrix.calpha.I2, 1, mean), col="red", lwd=1.25, type="l", ylim=c(0,200), ylab="N", xlab="Time", main="Infected")
for (i in 1:ncol(storeMatrix.calpha.I2)) {
  lines(0:150, storeMatrix.calpha.I2[,i], col=grey(.1, alpha=.15))
}
lines(0:150, apply(storeMatrix.calpha.I2, 1, mean), col="red", lwd=1.25, type="l", ylim=c(0,200), ylab="N", xlab="Time", main="Infected")
lines(lowerbnd.Ica2, lwd=1,lty=1, col="black")
lines(upperbnd.Ica2, lwd=1,lty=1,col="black")

## recovered
lowerbnd.Rca2 = apply(storeMatrix.calpha.R2,1,function(x){quantile(x,c(.025,.975))})[1,]
upperbnd.Rca2 = apply(storeMatrix.calpha.R2,1,function(x){quantile(x,c(.025,.975))})[2,]

plot(0:150, apply(storeMatrix.calpha.R2, 1, mean), col="red", lwd=1.25, type="l", ylim=c(0,200), ylab="N", xlab="Time", main="Recovered")
for (i in 1:ncol(storeMatrix.calpha.R2)) {
  lines(0:150, storeMatrix.calpha.R2[,i], col=grey(.1, alpha=.15))
}
lines(0:150, apply(storeMatrix.calpha.R2, 1, mean), col="red", lwd=1.25, type="l", ylim=c(0,200), ylab="N", xlab="Time", main="Recovered")
lines(lowerbnd.Rca2, lwd=1,lty=1, col="black")
lines(upperbnd.Rca2, lwd=1,lty=1,col="black")

par(mfrow=c(1,1))
## mean values 
plot(0:150, apply(storeMatrix.calpha.I2, 1, mean), col="red", lwd=1.75, type="l", ylim=c(0,200), ylab="N", xlab="Time", main="Simulation means")
lines(0:150, apply(storeMatrix.calpha.S2, 1, mean), col="blue", lwd=1.75, type="l", ylim=c(0,150), ylab="N", xlab="Time", main="Susceptible")
lines(0:150, apply(storeMatrix.calpha.R2, 1, mean), col="green", lwd=1.75, type="l", ylim=c(0,150), ylab="N", xlab="Time", main="Recovered")
legend("topright",legend=c("S","I","R"),fill=c("blue","red","green"), cex=0.75)

### neg cov contact and transmisibility

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

## plot
par(mfrow=c(1,3))

## susceptible
lowerbnd.Sct2 = apply(storeMatrix.ctau.S2,1,function(x){quantile(x,c(.025,.975))})[1,]
upperbnd.Sct2 = apply(storeMatrix.ctau.S2,1,function(x){quantile(x,c(.025,.975))})[2,]

plot(0:150, apply(storeMatrix.ctau.S2, 1, mean), col="red", lwd=1.25, type="l", ylim=c(0,200), ylab="Neg cov", xlab="Time", main="Susceptible")
for (i in 1:ncol(storeMatrix.ctau.S2)) {
  lines(0:150, storeMatrix.ctau.S2[,i], col=grey(.1, alpha=.15))
}
lines(0:150, apply(storeMatrix.ctau.S2, 1, mean), col="red", lwd=1.25, type="l", ylim=c(0,200), ylab="N", xlab="Time", main="Susceptible")
lines(lowerbnd.Sct2, lwd=1,lty=1, col="black")
lines(upperbnd.Sct2, lwd=1,lty=1,col="black")

## infected
lowerbnd.Ict2 = apply(storeMatrix.ctau.I2,1,function(x){quantile(x,c(.025,.975))})[1,]
upperbnd.Ict2 = apply(storeMatrix.ctau.I2,1,function(x){quantile(x,c(.025,.975))})[2,]

plot(0:150, apply(storeMatrix.ctau.I2, 1, mean), col="red", lwd=1.25, type="l", ylim=c(0,200), ylab="N", xlab="Time", main="Infected")
for (i in 1:ncol(storeMatrix.ctau.I2)) {
  lines(0:150, storeMatrix.ctau.I2[,i], col=grey(.1, alpha=.15))
}
lines(0:150, apply(storeMatrix.ctau.I2, 1, mean), col="red", lwd=1.25, type="l", ylim=c(0,200), ylab="N", xlab="Time", main="Infected")
lines(lowerbnd.Ict2, lwd=1,lty=1, col="black")
lines(upperbnd.Ict2, lwd=1,lty=1,col="black")

## recovered
lowerbnd.Rct2 = apply(storeMatrix.ctau.R2,1,function(x){quantile(x,c(.025,.975))})[1,]
upperbnd.Rct2 = apply(storeMatrix.ctau.R2,1,function(x){quantile(x,c(.025,.975))})[2,]

plot(0:150, apply(storeMatrix.ctau.R2, 1, mean), col="red", lwd=1.25, type="l", ylim=c(0,200), ylab="N", xlab="Time", main="Recovered")
for (i in 1:ncol(storeMatrix.ctau.R2)) {
  lines(0:150, storeMatrix.ctau.R2[,i], col=grey(.1, alpha=.15))
}
lines(0:150, apply(storeMatrix.ctau.R2, 1, mean), col="red", lwd=1.25, type="l", ylim=c(0,200), ylab="N", xlab="Time", main="Recovered")
lines(lowerbnd.Rct2, lwd=1,lty=1, col="black")
lines(upperbnd.Rct2, lwd=1,lty=1,col="black")

par(mfrow=c(1,1))
## mean values 
plot(0:150, apply(storeMatrix.ctau.I2, 1, mean), col="red", lwd=1.75, type="l", ylim=c(0,200), ylab="N", xlab="Time", main="Simulation means")
lines(0:150, apply(storeMatrix.ctau.S2, 1, mean), col="blue", lwd=1.75, type="l", ylim=c(0,150), ylab="N", xlab="Time", main="Susceptible")
lines(0:150, apply(storeMatrix.ctau.R2, 1, mean), col="green", lwd=1.75, type="l", ylim=c(0,150), ylab="N", xlab="Time", main="Recovered")
legend("topright",legend=c("S","I","R"),fill=c("blue","red","green"), cex=0.75)

############## POSITIVE COVARIANCE PLOTZZ ######################################
timeSeq <- 0:150
storeMatrix.taualpha.I3 <- array(NA, dim=c(length(timeSeq),length(out_poscov_taualpha2)))
for (j in 1:length(out_poscov_taualpha2)) {
  o <- out_poscov_taualpha2[[j]]
  storeMatrix.taualpha.I3[,j] <- o[,3] # num infected
}

## susceptible
storeMatrix.taualpha.S3 <- array(NA, dim=c(length(timeSeq),length(out_poscov_taualpha2)))
for (j in 1:length(out_poscov_taualpha2)) {
  o <- out_poscov_taualpha2[[j]]
  storeMatrix.taualpha.S3[,j] <- o[,2] # num susceptible
}

## recovered
storeMatrix.taualpha.R3 <- array(NA, dim=c(length(timeSeq),length(out_poscov_taualpha2)))
for (j in 1:length(out_poscov_taualpha2)) {
  o <- out_poscov_taualpha2[[j]]
  storeMatrix.taualpha.R3[,j] <- o[,4] # num recovered
}

# plot
par(mfrow=c(3,3))

##infected
lowerbnd.I3 = apply(storeMatrix.taualpha.I3,1,function(x){quantile(x,c(.025,.975))})[1,]
upperbnd.I = apply(storeMatrix.taualpha.I3,1,function(x){quantile(x,c(.025,.975))})[2,]

plot(0:150, apply(storeMatrix.taualpha.I3, 1, mean), col="red", lwd=1.25, type="l", ylim=c(0,100), ylab="N", xlab="Time", main="Infected")
for (i in 1:ncol(storeMatrix.taualpha.I3)) {
  lines(0:150, storeMatrix.taualpha.I3[,i], col=grey(.1, alpha=.15))
}
lines(0:150, apply(storeMatrix.taualpha.I3, 1, mean), col="red", lwd=1.25, type="l", ylim=c(0,100), ylab="N", xlab="Time", main="Infected")
lines(lowerbnd.I3, lwd=1,lty=1, col="black")
lines(upperbnd.I3, lwd=1,lty=1,col="black")

## susceptible
lowerbnd.S3 = apply(storeMatrix.taualpha.S3,1,function(x){quantile(x,c(.025,.975))})[1,]
upperbnd.S3 = apply(storeMatrix.taualpha.S3,1,function(x){quantile(x,c(.025,.975))})[2,]

plot(0:150, apply(storeMatrix.taualpha.S3, 1, mean), col="red", lwd=1.25, type="l", ylim=c(0,100), ylab="Pos cov", xlab="Time", main="Susceptible")
for (i in 1:ncol(storeMatrix.taualpha.S3)) {
  lines(0:150, storeMatrix.taualpha.S3[,i], col=grey(.1, alpha=.15))
}
lines(0:150, apply(storeMatrix.taualpha.S3, 1, mean), col="red", lwd=1.25, type="l", ylim=c(0,100), ylab="N", xlab="Time", main="Susceptible")
lines(lowerbnd.S3, lwd=1,lty=1, col="black")
lines(upperbnd.S3, lwd=1,lty=1,col="black")

## recovered
lowerbnd.R3 = apply(storeMatrix.taualpha.R3,1,function(x){quantile(x,c(.025,.975))})[1,]
upperbnd.R3 = apply(storeMatrix.taualpha.R3,1,function(x){quantile(x,c(.025,.975))})[2,]

plot(0:150, apply(storeMatrix.taualpha.R3, 1, mean), col="red", lwd=1.25, type="l", ylim=c(0,100), ylab="N", xlab="Time", main="Recovered")
for (i in 1:ncol(storeMatrix.taualpha.R3)) {
  lines(0:150, storeMatrix.taualpha.R3[,i], col=grey(.1, alpha=.15))
}
lines(0:150, apply(storeMatrix.taualpha.R3, 1, mean), col="red", lwd=1.25, type="l", ylim=c(0,100), ylab="N", xlab="Time", main="Recovered")
lines(lowerbnd.R3, lwd=1,lty=1, col="black")
lines(upperbnd.R3, lwd=1,lty=1,col="black")

par(mfrow=c(1,1))
## mean values 
plot(0:150, apply(storeMatrix.taualpha.I3, 1, mean), col="red", lwd=1.75, type="l", ylim=c(0,150), ylab="N", xlab="Time", main="Pos cov")
lines(0:150, apply(storeMatrix.taualpha.S3, 1, mean), col="blue", lwd=1.75, type="l", ylim=c(0,150), ylab="N", xlab="Time", main="Susceptible")
lines(0:150, apply(storeMatrix.taualpha.R3, 1, mean), col="green", lwd=1.75, type="l", ylim=c(0,150), ylab="N", xlab="Time", main="Recovered")
legend("topright",legend=c("S","I","R"),fill=c("blue","red","green"), cex=0.75)

#### tau gamma pos cov
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

## plot
par(mfrow=c(3,3))

## susceptible
lowerbnd.Stg3 = apply(storeMatrix.taugamma.S3,1,function(x){quantile(x,c(.025,.975))})[1,]
upperbnd.Stg3 = apply(storeMatrix.taugamma.S3,1,function(x){quantile(x,c(.025,.975))})[2,]

plot(0:150, apply(storeMatrix.taugamma.S3, 1, mean), col="red", lwd=1.25, type="l", ylim=c(0,200), ylab="Pos cov", xlab="Time", main="Susceptible")
for (i in 1:ncol(storeMatrix.taugamma.S3)) {
  lines(0:150, storeMatrix.taugamma.S3[,i], col=grey(.1, alpha=.15))
}
lines(0:150, apply(storeMatrix.taugamma.S3, 1, mean), col="red", lwd=1.25, type="l", ylim=c(0,200), ylab="N", xlab="Time", main="Susceptible")
lines(lowerbnd.Stg3, lwd=1,lty=1, col="black")
lines(upperbnd.Stg3, lwd=1,lty=1,col="black")

## infected
lowerbnd.Itg3 = apply(storeMatrix.taugamma.I3,1,function(x){quantile(x,c(.025,.975))})[1,]
upperbnd.Itg3 = apply(storeMatrix.taugamma.I3,1,function(x){quantile(x,c(.025,.975))})[2,]

plot(0:150, apply(storeMatrix.taugamma.I3, 1, mean), col="red", lwd=1.25, type="l", ylim=c(0,200), ylab="N", xlab="Time", main="Infected")
for (i in 1:ncol(storeMatrix.taugamma.I3)) {
  lines(0:150, storeMatrix.taugamma.I3[,i], col=grey(.1, alpha=.15))
}
lines(0:150, apply(storeMatrix.taugamma.I3, 1, mean), col="red", lwd=1.25, type="l", ylim=c(0,200), ylab="N", xlab="Time", main="Infected")
lines(lowerbnd.Itg3, lwd=1,lty=1, col="black")
lines(upperbnd.Itg3, lwd=1,lty=1,col="black")

## recovered
lowerbnd.Rtg3 = apply(storeMatrix.taugamma.R3,1,function(x){quantile(x,c(.025,.975))})[1,]
upperbnd.Rtg3 = apply(storeMatrix.taugamma.R3,1,function(x){quantile(x,c(.025,.975))})[2,]

plot(0:150, apply(storeMatrix.taugamma.R3, 1, mean), col="red", lwd=1.25, type="l", ylim=c(0,200), ylab="N", xlab="Time", main="Recovered")
for (i in 1:ncol(storeMatrix.taugamma.R3)) {
  lines(0:150, storeMatrix.taugamma.R3[,i], col=grey(.1, alpha=.15))
}
lines(0:150, apply(storeMatrix.taugamma.R3, 1, mean), col="red", lwd=1.25, type="l", ylim=c(0,200), ylab="N", xlab="Time", main="Recovered")
lines(lowerbnd.Rtg3, lwd=1,lty=1, col="black")
lines(upperbnd.Rtg3, lwd=1,lty=1,col="black")

par(mfrow=c(1,3))
## mean values 
plot(0:150, apply(storeMatrix.taugamma.I3, 1, mean), col="red", lwd=1.75, type="l", ylim=c(0,200), ylab="N", xlab="Time", main="Pos cov")
lines(0:150, apply(storeMatrix.taugamma.S3, 1, mean), col="blue", lwd=1.75, type="l", ylim=c(0,150), ylab="N", xlab="Time", main="Susceptible")
lines(0:150, apply(storeMatrix.taugamma.R3, 1, mean), col="green", lwd=1.75, type="l", ylim=c(0,150), ylab="N", xlab="Time", main="Recovered")
legend("topright",legend=c("S","I","R"),fill=c("blue","red","green"), cex=0.75)

### pos cov contact and gamma

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

## plot
par(mfrow=c(3,3))

## susceptible
lowerbnd.Scg3 = apply(storeMatrix.cgamma.S3,1,function(x){quantile(x,c(.025,.975))})[1,]
upperbnd.Scg3 = apply(storeMatrix.cgamma.S3,1,function(x){quantile(x,c(.025,.975))})[2,]

plot(0:150, apply(storeMatrix.cgamma.S3, 1, mean), col="red", lwd=1.25, type="l", ylim=c(0,200), ylab="Pos cov", xlab="Time", main="Susceptible")
for (i in 1:ncol(storeMatrix.cgamma.S3)) {
  lines(0:150, storeMatrix.cgamma.S3[,i], col=grey(.1, alpha=.15))
}
lines(0:150, apply(storeMatrix.cgamma.S3, 1, mean), col="red", lwd=1.25, type="l", ylim=c(0,200), ylab="N", xlab="Time", main="Susceptible")
lines(lowerbnd.Scg3, lwd=1,lty=1, col="black")
lines(upperbnd.Scg3, lwd=1,lty=1,col="black")

## infected
lowerbnd.Icg3 = apply(storeMatrix.cgamma.I3,1,function(x){quantile(x,c(.025,.975))})[1,]
upperbnd.Icg3 = apply(storeMatrix.cgamma.I3,1,function(x){quantile(x,c(.025,.975))})[2,]

plot(0:150, apply(storeMatrix.cgamma.I3, 1, mean), col="red", lwd=1.25, type="l", ylim=c(0,200), ylab="N", xlab="Time", main="Infected")
for (i in 1:ncol(storeMatrix.cgamma.I3)) {
  lines(0:150, storeMatrix.cgamma.I3[,i], col=grey(.1, alpha=.15))
}
lines(0:150, apply(storeMatrix.cgamma.I3, 1, mean), col="red", lwd=1.25, type="l", ylim=c(0,200), ylab="N", xlab="Time", main="Infected")
lines(lowerbnd.Icg3, lwd=1,lty=1, col="black")
lines(upperbnd.Icg3, lwd=1,lty=1,col="black")

## recovered
lowerbnd.Rcg3 = apply(storeMatrix.cgamma.R3,1,function(x){quantile(x,c(.025,.975))})[1,]
upperbnd.Rcg3 = apply(storeMatrix.cgamma.R3,1,function(x){quantile(x,c(.025,.975))})[2,]

plot(0:150, apply(storeMatrix.cgamma.R3, 1, mean), col="red", lwd=1.25, type="l", ylim=c(0,200), ylab="N", xlab="Time", main="Recovered")
for (i in 1:ncol(storeMatrix.cgamma.R3)) {
  lines(0:150, storeMatrix.cgamma.R3[,i], col=grey(.1, alpha=.15))
}
lines(0:150, apply(storeMatrix.cgamma.R3, 1, mean), col="red", lwd=1.25, type="l", ylim=c(0,200), ylab="N", xlab="Time", main="Recovered")
lines(lowerbnd.Rcg3, lwd=1,lty=1, col="black")
lines(upperbnd.Rcg3, lwd=1,lty=1,col="black")

par(mfrow=c(1,3))
## mean values no
plot(0:150, apply(storeMatrix.cgamma.I3, 1, mean), col="red", lwd=1.75, type="l", ylim=c(0,200), ylab="N", xlab="Time", main="Pos cov")
lines(0:150, apply(storeMatrix.cgamma.S3, 1, mean), col="blue", lwd=1.75, type="l", ylim=c(0,150), ylab="N", xlab="Time", main="Susceptible")
lines(0:150, apply(storeMatrix.cgamma.R3, 1, mean), col="green", lwd=1.75, type="l", ylim=c(0,150), ylab="N", xlab="Time", main="Recovered")
legend("topright",legend=c("S","I","R"), fill=c("blue","red","green"), cex=.5)

### pos cov contact and alpha

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

## plot
par(mfrow=c(3,3))

## susceptible
lowerbnd.Sca3 = apply(storeMatrix.calpha.S3,1,function(x){quantile(x,c(.025,.975))})[1,]
upperbnd.Sca3 = apply(storeMatrix.calpha.S3,1,function(x){quantile(x,c(.025,.975))})[2,]

plot(0:150, apply(storeMatrix.calpha.S3, 1, mean), col="red", lwd=1.25, type="l", ylim=c(0,200), ylab="Pos cov", xlab="Time", main="Susceptible")
for (i in 1:ncol(storeMatrix.calpha.S3)) {
  lines(0:150, storeMatrix.calpha.S3[,i], col=grey(.1, alpha=.15))
}
lines(0:150, apply(storeMatrix.calpha.S3, 1, mean), col="red", lwd=1.25, type="l", ylim=c(0,200), ylab="N", xlab="Time", main="Susceptible")
lines(lowerbnd.Sca3, lwd=1,lty=1, col="black")
lines(upperbnd.Sca3, lwd=1,lty=1,col="black")

## infected
lowerbnd.Ica3 = apply(storeMatrix.calpha.I3,1,function(x){quantile(x,c(.025,.975))})[1,]
upperbnd.Ica3 = apply(storeMatrix.calpha.I3,1,function(x){quantile(x,c(.025,.975))})[2,]

plot(0:150, apply(storeMatrix.calpha.I3, 1, mean), col="red", lwd=1.25, type="l", ylim=c(0,200), ylab="N", xlab="Time", main="Infected")
for (i in 1:ncol(storeMatrix.calpha.I3)) {
  lines(0:150, storeMatrix.calpha.I3[,i], col=grey(.1, alpha=.15))
}
lines(0:150, apply(storeMatrix.calpha.I3, 1, mean), col="red", lwd=1.25, type="l", ylim=c(0,200), ylab="N", xlab="Time", main="Infected")
lines(lowerbnd.Ica3, lwd=1,lty=1, col="black")
lines(upperbnd.Ica3, lwd=1,lty=1,col="black")

## recovered
lowerbnd.Rca3 = apply(storeMatrix.calpha.R3,1,function(x){quantile(x,c(.025,.975))})[1,]
upperbnd.Rca3 = apply(storeMatrix.calpha.R3,1,function(x){quantile(x,c(.025,.975))})[2,]

plot(0:150, apply(storeMatrix.calpha.R3, 1, mean), col="red", lwd=1.25, type="l", ylim=c(0,200), ylab="N", xlab="Time", main="Recovered")
for (i in 1:ncol(storeMatrix.calpha.R3)) {
  lines(0:150, storeMatrix.calpha.R3[,i], col=grey(.1, alpha=.15))
}
lines(0:150, apply(storeMatrix.calpha.R3, 1, mean), col="red", lwd=1.25, type="l", ylim=c(0,200), ylab="N", xlab="Time", main="Recovered")
lines(lowerbnd.Rca3, lwd=1,lty=1, col="black")
lines(upperbnd.Rca3, lwd=1,lty=1,col="black")

par(mfrow=c(1,3))
## mean values 
plot(0:150, apply(storeMatrix.calpha.I3, 1, mean), col="red", lwd=1.75, type="l", ylim=c(0,200), ylab="N", xlab="Time", main="Pos cov")
lines(0:150, apply(storeMatrix.calpha.S3, 1, mean), col="blue", lwd=1.75, type="l", ylim=c(0,150), ylab="N", xlab="Time", main="Susceptible")
lines(0:150, apply(storeMatrix.calpha.R3, 1, mean), col="green", lwd=1.75, type="l", ylim=c(0,150), ylab="N", xlab="Time", main="Recovered")
legend("topright",legend=c("S","I","R"),fill=c("blue","red","green"), cex=0.5)

### pos cov contact and transmisibility

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

## plot
par(mfrow=c(3,3))

## susceptible
lowerbnd.Sct3 = apply(storeMatrix.ctau.S3,1,function(x){quantile(x,c(.025,.975))})[1,]
upperbnd.Sct3 = apply(storeMatrix.ctau.S3,1,function(x){quantile(x,c(.025,.975))})[2,]

plot(0:150, apply(storeMatrix.ctau.S3, 1, mean), col="red", lwd=1.25, type="l", ylim=c(0,200), ylab="N", xlab="Time", main="Susceptible")
for (i in 1:ncol(storeMatrix.ctau.S3)) {
  lines(0:150, storeMatrix.ctau.S3[,i], col=grey(.1, alpha=.15))
}
lines(0:150, apply(storeMatrix.ctau.S3, 1, mean), col="red", lwd=1.25, type="l", ylim=c(0,200), ylab="N", xlab="Time", main="Susceptible")
lines(lowerbnd.Sct3, lwd=1,lty=1, col="black")
lines(upperbnd.Sct3, lwd=1,lty=1,col="black")

## infected
lowerbnd.Ict3 = apply(storeMatrix.ctau.I3,1,function(x){quantile(x,c(.025,.975))})[1,]
upperbnd.Ict3 = apply(storeMatrix.ctau.I3,1,function(x){quantile(x,c(.025,.975))})[2,]

plot(0:150, apply(storeMatrix.ctau.I3, 1, mean), col="red", lwd=1.25, type="l", ylim=c(0,200), ylab="N", xlab="Time", main="Infected")
for (i in 1:ncol(storeMatrix.ctau.I3)) {
  lines(0:150, storeMatrix.ctau.I3[,i], col=grey(.1, alpha=.15))
}
lines(0:150, apply(storeMatrix.ctau.I3, 1, mean), col="red", lwd=1.25, type="l", ylim=c(0,200), ylab="N", xlab="Time", main="Infected")
lines(lowerbnd.Ict3, lwd=1,lty=1, col="black")
lines(upperbnd.Ict3, lwd=1,lty=1,col="black")

## recovered
lowerbnd.Rct3 = apply(storeMatrix.ctau.R3,1,function(x){quantile(x,c(.025,.975))})[1,]
upperbnd.Rct3 = apply(storeMatrix.ctau.R3,1,function(x){quantile(x,c(.025,.975))})[2,]

plot(0:150, apply(storeMatrix.ctau.R3, 1, mean), col="red", lwd=1.25, type="l", ylim=c(0,200), ylab="N", xlab="Time", main="Recovered")
for (i in 1:ncol(storeMatrix.ctau.R3)) {
  lines(0:150, storeMatrix.ctau.R3[,i], col=grey(.1, alpha=.15))
}
lines(0:150, apply(storeMatrix.ctau.R3, 1, mean), col="red", lwd=1.25, type="l", ylim=c(0,200), ylab="N", xlab="Time", main="Recovered")
lines(lowerbnd.Rct3, lwd=1,lty=1, col="black")
lines(upperbnd.Rct3, lwd=1,lty=1,col="black")

par(mfrow=c(1,3))
## mean values 
plot(0:150, apply(storeMatrix.ctau.I3, 1, mean), col="red", lwd=1.75, type="l", ylim=c(0,200), ylab="N", xlab="Time", main="Pos cov")
lines(0:150, apply(storeMatrix.ctau.S3, 1, mean), col="blue", lwd=1.75, type="l", ylim=c(0,150), ylab="N", xlab="Time", main="Susceptible")
lines(0:150, apply(storeMatrix.ctau.R3, 1, mean), col="green", lwd=1.75, type="l", ylim=c(0,150), ylab="N", xlab="Time", main="Recovered")
legend("topright",legend=c("S","I","R"),fill=c("blue","red","green"), cex=0.5)
