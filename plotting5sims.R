#### Just my big plots

### Strat Var A
par(mfrow=c(3,3))
plot(0:150, apply(storeMatrix.alpha.I, 1, mean), col="red", lwd=1.75, type="l", ylim=c(0,200), ylab="N", xlab="Time", main="Variation in Alpha")
lines(0:150, apply(storeMatrix.alpha.S, 1, mean), col="blue", lwd=1.75, type="l", ylim=c(0,150), ylab="N", xlab="Time")
lines(0:150, apply(storeMatrix.alpha.R, 1, mean), col="green", lwd=1.75, type="l", ylim=c(0,150), ylab="N", xlab="Time")
legend("topright",legend=c("S","I","R"),fill=c("blue","red","green"), cex=0.25)

### Strat Var B
plot(0:150, apply(storeMatrix.beta.I, 1, mean), col="red", lwd=1.75, type="l", ylim=c(0,250), ylab="N", xlab="Time", main="Variation in Beta")
lines(0:150, apply(storeMatrix.beta.S, 1, mean), col="blue", lwd=1.75, type="l", ylim=c(0,150), ylab="N", xlab="Time")
lines(0:150, apply(storeMatrix.beta.R, 1, mean), col="green", lwd=1.75, type="l", ylim=c(0,150), ylab="N", xlab="Time")
legend("topright",legend=c("S","I","R"),fill=c("blue","red","green"), cex=0.25)

### Strat Var G
plot(0:150, apply(storeMatrix.gamma.I, 1, mean), col="red", lwd=1.75, type="l", ylim=c(0,250), ylab="N", xlab="Time", main="Variation in Gamma")
lines(0:150, apply(storeMatrix.gamma.S, 1, mean), col="blue", lwd=1.75, type="l", ylim=c(0,150), ylab="N", xlab="Time")
lines(0:150, apply(storeMatrix.gamma.R, 1, mean), col="green", lwd=1.75, type="l", ylim=c(0,150), ylab="N", xlab="Time")
legend("topright",legend=c("S","I","R"),fill=c("blue","red","green"), cex=0.25)


### cleaning sims
## infected
timeSeq <- 0:150
storeMatrix.cont.varB.I <- array(NA, dim=c(length(timeSeq),length(out_cont_var_beta)))
for (j in 1:length(out_cont_var_beta)) {
  o <- out_cont_var_beta[[j]]
  storeMatrix.cont.varB.I[,j] <- o[,3] # num infected
}

## susceptible
storeMatrix.cont.varB.S <- array(NA, dim=c(length(timeSeq),length(out_cont_var_beta)))
for (j in 1:length(out_cont_var_beta)) {
  o <- out_cont_var_beta[[j]]
  storeMatrix.cont.varB.S[,j] <- o[,2] # num susceptible
}

## recovered
storeMatrix.cont.varB.R <- array(NA, dim=c(length(timeSeq),length(out_cont_var_beta)))
for (j in 1:length(out_cont_var_beta)) {
  o <- out_cont_var_beta[[j]]
  storeMatrix.cont.varB.R[,j] <- o[,4] # num recovered
}

## cont var b
plot(0:150, apply(storeMatrix.cont.varB.I, 1, mean), col="red", lwd=1.75, type="l", ylim=c(0,250), ylab="N", xlab="Time", main="Cont. Variation in Beta")
lines(0:150, apply(storeMatrix.cont.varB.S, 1, mean), col="blue", lwd=1.75, type="l", ylim=c(0,150), ylab="N", xlab="Time")
lines(0:150, apply(storeMatrix.cont.varB.R, 1, mean), col="green", lwd=1.75, type="l", ylim=c(0,150), ylab="N", xlab="Time")
legend("topright",legend=c("S","I","R"),fill=c("blue","red","green"), cex=0.25)

## cont var G
plot(0:150, apply(storeMatrix.cont.varG.I, 1, mean), col="red", lwd=1.75, type="l", ylim=c(0,250), ylab="N", xlab="Time", main="Cont. Variation in Gamma")
lines(0:150, apply(storeMatrix.cont.varG.S, 1, mean), col="blue", lwd=1.75, type="l", ylim=c(0,150), ylab="N", xlab="Time")
lines(0:150, apply(storeMatrix.cont.varG.R, 1, mean), col="green", lwd=1.75, type="l", ylim=c(0,150), ylab="N", xlab="Time")
legend("topright",legend=c("S","I","R"),fill=c("blue","red","green"), cex=0.25)

## cont var A
plot(0:150, apply(storeMatrix.cont.varA.I, 1, mean), col="red", lwd=1.75, type="l", ylim=c(0,250), ylab="N", xlab="Time", main="Cont. Variation in Alpha")
lines(0:150, apply(storeMatrix.cont.varA.S, 1, mean), col="blue", lwd=1.75, type="l", ylim=c(0,150), ylab="N", xlab="Time")
lines(0:150, apply(storeMatrix.cont.varA.R, 1, mean), col="green", lwd=1.75, type="l", ylim=c(0,150), ylab="N", xlab="Time")
legend("topright",legend=c("S","I","R"),fill=c("blue","red","green"), cex=0.25)

## no var GEM
plot(0:150, apply(storeMatrix.no.var.I, 1, mean), col="red", lwd=1.75, type="l", ylim=c(0,250), ylab="N", xlab="Time", main="No Var GEM")
lines(0:150, apply(storeMatrix.no.var.S, 1, mean), col="blue", lwd=1.75, type="l", ylim=c(0,150), ylab="N", xlab="Time")
lines(0:150, apply(storeMatrix.no.var.R, 1, mean), col="green", lwd=1.75, type="l", ylim=c(0,150), ylab="N", xlab="Time")
legend("topright",legend=c("S","I","R"),fill=c("blue","red","green"), cex=0.25)


## infected
timeSeq <- 0:150
storeMatrix.cont.varG.I <- array(NA, dim=c(length(timeSeq),length(out_cont_var_gamma)))
for (j in 1:length(out_cont_var_gamma)) {
  o <- out_cont_var_gamma[[j]]
  storeMatrix.cont.varG.I[,j] <- o[,3] # num infected
}

## susceptible
storeMatrix.cont.varG.S <- array(NA, dim=c(length(timeSeq),length(out_cont_var_gamma)))
for (j in 1:length(out_cont_var_gamma)) {
  o <- out_cont_var_gamma[[j]]
  storeMatrix.cont.varG.S[,j] <- o[,2] # num susceptible
}

## recovered
storeMatrix.cont.varG.R <- array(NA, dim=c(length(timeSeq),length(out_cont_var_gamma)))
for (j in 1:length(out_cont_var_gamma)) {
  o <- out_cont_var_gamma[[j]]
  storeMatrix.cont.varG.R[,j] <- o[,4] # num recovered
}

## cont. alpha var
## infected
timeSeq <- 0:150
storeMatrix.cont.varA.I <- array(NA, dim=c(length(timeSeq),length(out_cont_var_alpha)))
for (j in 1:length(out_cont_var_alpha)) {
  o <- out_cont_var_alpha[[j]]
  storeMatrix.cont.varA.I[,j] <- o[,3] # num infected
}

## susceptible
storeMatrix.cont.varA.S <- array(NA, dim=c(length(timeSeq),length(out_cont_var_alpha)))
for (j in 1:length(out_cont_var_alpha)) {
  o <- out_cont_var_alpha[[j]]
  storeMatrix.cont.varA.S[,j] <- o[,2] # num susceptible
}

## recovered
storeMatrix.cont.varA.R <- array(NA, dim=c(length(timeSeq),length(out_cont_var_alpha)))
for (j in 1:length(out_cont_var_alpha)) {
  o <- out_cont_var_alpha[[j]]
  storeMatrix.cont.varA.R[,j] <- o[,4] # num recovered
}

## no var GEM
## infected
timeSeq <- 0:150
storeMatrix.no.var.I <- array(NA, dim=c(length(timeSeq),length(out_no_var)))
for (j in 1:length(out_no_var)) {
  o <- out_no_var[[j]]
  storeMatrix.no.var.I[,j] <- o[,3] # num infected
}

## susceptible
storeMatrix.no.var.S <- array(NA, dim=c(length(timeSeq),length(out_no_var)))
for (j in 1:length(out_no_var)) {
  o <- out_no_var[[j]]
  storeMatrix.no.var.S[,j] <- o[,2] # num susceptible
}

## recovered
storeMatrix.no.var.R <- array(NA, dim=c(length(timeSeq),length(out_no_var)))
for (j in 1:length(out_no_var)) {
  o <- out_no_var[[j]]
  storeMatrix.no.var.R[,j] <- o[,4] # num recovered
}
