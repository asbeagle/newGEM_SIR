#### Just my big plots

### Strat Var A
par(mfrow=c(3,3))
plot(0:150, apply(storeMatrix.alpha.I, 1, mean), col="red", lwd=1.75, type="l", ylim=c(0,200), ylab="N", xlab="Time", main="Variation in Alpha")
lines(0:150, apply(storeMatrix.alpha.S, 1, mean), col="blue", lwd=1.75, type="l", ylim=c(0,150), ylab="N", xlab="Time")
lines(0:150, apply(storeMatrix.alpha.R, 1, mean), col="green", lwd=1.75, type="l", ylim=c(0,150), ylab="N", xlab="Time")
legend("topright",legend=c("S","I","R"),fill=c("blue","red","green"), cex=0.25)

### Strat Var B
plot(0:150, apply(storeMatrix.beta.I, 1, mean), col="red", lwd=1.75, type="l", ylim=c(0,200), ylab="N", xlab="Time", main="Variation in Beta")
lines(0:150, apply(storeMatrix.beta.S, 1, mean), col="blue", lwd=1.75, type="l", ylim=c(0,150), ylab="N", xlab="Time")
lines(0:150, apply(storeMatrix.beta.R, 1, mean), col="green", lwd=1.75, type="l", ylim=c(0,150), ylab="N", xlab="Time")
legend("topright",legend=c("S","I","R"),fill=c("blue","red","green"), cex=0.25)

### Strat Var G
plot(0:150, apply(storeMatrix.gamma.I, 1, mean), col="red", lwd=1.75, type="l", ylim=c(0,200), ylab="N", xlab="Time", main="Variation in Gamma")
lines(0:150, apply(storeMatrix.gamma.S, 1, mean), col="blue", lwd=1.75, type="l", ylim=c(0,150), ylab="N", xlab="Time")
lines(0:150, apply(storeMatrix.gamma.R, 1, mean), col="green", lwd=1.75, type="l", ylim=c(0,150), ylab="N", xlab="Time")
legend("topright",legend=c("S","I","R"),fill=c("blue","red","green"), cex=0.25)


### cleaning sims


## cont var b
plot(0:150, apply(storeMatrix.cont.varB.I, 1, mean), col="red", lwd=1.75, type="l", ylim=c(0,200), ylab="N", xlab="Time", main="Cont. Variation in Beta")
lines(0:150, apply(storeMatrix.cont.varB.S, 1, mean), col="blue", lwd=1.75, type="l", ylim=c(0,150), ylab="N", xlab="Time")
lines(0:150, apply(storeMatrix.cont.varB.R, 1, mean), col="green", lwd=1.75, type="l", ylim=c(0,150), ylab="N", xlab="Time")
legend("topright",legend=c("S","I","R"),fill=c("blue","red","green"), cex=0.25)

## cont var G
plot(0:150, apply(storeMatrix.cont.varG.I, 1, mean), col="red", lwd=1.75, type="l", ylim=c(0,200), ylab="N", xlab="Time", main="Cont. Variation in Gamma")
lines(0:150, apply(storeMatrix.cont.varG.S, 1, mean), col="blue", lwd=1.75, type="l", ylim=c(0,150), ylab="N", xlab="Time")
lines(0:150, apply(storeMatrix.cont.varG.R, 1, mean), col="green", lwd=1.75, type="l", ylim=c(0,150), ylab="N", xlab="Time")
legend("topright",legend=c("S","I","R"),fill=c("blue","red","green"), cex=0.25)

## cont var A
plot(0:150, apply(storeMatrix.cont.varA.I, 1, mean), col="red", lwd=1.75, type="l", ylim=c(0,200), ylab="N", xlab="Time", main="Cont. Variation in Alpha")
lines(0:150, apply(storeMatrix.cont.varA.S, 1, mean), col="blue", lwd=1.75, type="l", ylim=c(0,150), ylab="N", xlab="Time")
lines(0:150, apply(storeMatrix.cont.varA.R, 1, mean), col="green", lwd=1.75, type="l", ylim=c(0,150), ylab="N", xlab="Time")
legend("topright",legend=c("S","I","R"),fill=c("blue","red","green"), cex=0.25)

## no var GEM
plot(0:150, apply(storeMatrix.no.var.I, 1, mean), col="red", lwd=1.75, type="l", ylim=c(0,200), ylab="N", xlab="Time", main="No Var GEM")
lines(0:150, apply(storeMatrix.no.var.S, 1, mean), col="blue", lwd=1.75, type="l", ylim=c(0,150), ylab="N", xlab="Time")
lines(0:150, apply(storeMatrix.no.var.R, 1, mean), col="green", lwd=1.75, type="l", ylim=c(0,150), ylab="N", xlab="Time")
legend("topright",legend=c("S","I","R"),fill=c("blue","red","green"), cex=0.25)

### CONTINUOUS BETA
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

### STRAT BETA
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


### UNIFORM BETA
timeSeq <- 0:150
storeMatrix.uni.varB.I <- array(NA, dim=c(length(timeSeq),length(out_uniform_var_beta)))
for (j in 1:length(out_uniform_var_beta)) {
  o <- out_uniform_var_beta[[j]]
  storeMatrix.uni.varB.I[,j] <- o[,3] # num infected
}

## susceptible
storeMatrix.uni.varB.S <- array(NA, dim=c(length(timeSeq),length(out_uniform_var_beta)))
for (j in 1:length(out_uniform_var_beta)) {
  o <- out_uniform_var_beta[[j]]
  storeMatrix.uni.varB.S[,j] <- o[,2] # num susceptible
}

## recovered
storeMatrix.uni.varB.R <- array(NA, dim=c(length(timeSeq),length(out_uniform_var_beta)))
for (j in 1:length(out_uniform_var_beta)) {
  o <- out_uniform_var_beta[[j]]
  storeMatrix.uni.varB.R[,j] <- o[,4] # num recovered
}

par(mfrow=c(1,1))
plot(0:150, apply(storeMatrix.uni.varB.I, 1, mean), col="red", lwd=1.75, type="l", ylim=c(0,200), ylab="N", xlab="Time", main="Uniform Var")
lines(0:150, apply(storeMatrix.uni.varB.S, 1, mean), col="blue", lwd=1.75, type="l", ylim=c(0,150), ylab="N", xlab="Time")
lines(0:150, apply(storeMatrix.uni.varB.R, 1, mean), col="green", lwd=1.75, type="l", ylim=c(0,150), ylab="N", xlab="Time")
legend("topright",legend=c("S","I","R"),fill=c("blue","red","green"), cex=0.25)

plot(0:150, apply(storeMatrix.beta.I, 1, mean), col="red", lwd=1.75, type="l", ylim=c(0,200), ylab="N", xlab="Time", main="Strat Var")
lines(0:150, apply(storeMatrix.beta.S, 1, mean), col="blue", lwd=1.75, type="l", ylim=c(0,150), ylab="N", xlab="Time")
lines(0:150, apply(storeMatrix.beta.R, 1, mean), col="green", lwd=1.75, type="l", ylim=c(0,150), ylab="N", xlab="Time")
legend("topright",legend=c("S","I","R"),fill=c("blue","red","green"), cex=0.25)

plot(0:150, apply(storeMatrix.cont.varB.I, 1, mean), col="red", lwd=1.75, type="l", ylim=c(0,200), ylab="N", xlab="Time", main="Cont. Var")
lines(0:150, apply(storeMatrix.cont.varB.S, 1, mean), col="blue", lwd=1.75, type="l", ylim=c(0,150), ylab="N", xlab="Time")
lines(0:150, apply(storeMatrix.cont.varB.R, 1, mean), col="green", lwd=1.75, type="l", ylim=c(0,150), ylab="N", xlab="Time")
legend("topright",legend=c("S","I","R"),fill=c("blue","red","green"), cex=0.25)

## no var GEM
plot(0:150, apply(storeMatrix.no.var.I, 1, mean), col="red", lwd=1.75, type="l", ylim=c(0,200), ylab="N", xlab="Time", main="No Var")
lines(0:150, apply(storeMatrix.no.var.S, 1, mean), col="blue", lwd=1.75, type="l", ylim=c(0,150), ylab="N", xlab="Time")
lines(0:150, apply(storeMatrix.no.var.R, 1, mean), col="green", lwd=1.75, type="l", ylim=c(0,150), ylab="N", xlab="Time")
legend("topright",legend=c("S","I","R"),fill=c("blue","red","green"), cex=0.25)


### CONTINUOUS GAMMA
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

### STRAT GAMMA
## infected
timeSeq <- 0:150
storeMatrix.varG.I <- array(NA, dim=c(length(timeSeq),length(out_strat_var_gamma)))
for (j in 1:length(out_strat_var_gamma)) {
  o <- out_strat_var_gamma[[j]]
  storeMatrix.varG.I[,j] <- o[,3] # num infected
}

## susceptible
storeMatrix.varG.S <- array(NA, dim=c(length(timeSeq),length(out_strat_var_gamma)))
for (j in 1:length(out_strat_var_gamma)) {
  o <- out_strat_var_gamma[[j]]
  storeMatrix.varG.S[,j] <- o[,2] # num susceptible
}

## recovered
storeMatrix.varG.R <- array(NA, dim=c(length(timeSeq),length(out_strat_var_gamma)))
for (j in 1:length(out_strat_var_gamma)) {
  o <- out_strat_var_gamma[[j]]
  storeMatrix.varG.R[,j] <- o[,4] # num recovered
}

## UNIFORM GAMMA
timeSeq <- 0:150
storeMatrix.uni.varG.I <- array(NA, dim=c(length(timeSeq),length(out_uniform_var_gamma)))
for (j in 1:length(out_uniform_var_gamma)) {
  o <- out_uniform_var_gamma[[j]]
  storeMatrix.uni.varG.I[,j] <- o[,3] # num infected
}

## susceptible
storeMatrix.uni.varG.S <- array(NA, dim=c(length(timeSeq),length(out_uniform_var_gamma)))
for (j in 1:length(out_uniform_var_gamma)) {
  o <- out_uniform_var_gamma[[j]]
  storeMatrix.uni.varG.S[,j] <- o[,2] # num susceptible
}

## recovered
storeMatrix.uni.varG.R <- array(NA, dim=c(length(timeSeq),length(out_uniform_var_gamma)))
for (j in 1:length(out_uniform_var_gamma)) {
  o <- out_uniform_var_gamma[[j]]
  storeMatrix.uni.varG.R[,j] <- o[,4] # num recovered
}

par(mfrow=c(1,4))
plot(0:150, apply(storeMatrix.uni.varG.I, 1, mean), col="red", lwd=1.75, type="l", ylim=c(0,200), ylab="N", xlab="Time", main="Uniform Var")
lines(0:150, apply(storeMatrix.uni.varG.S, 1, mean), col="blue", lwd=1.75, type="l", ylim=c(0,150), ylab="N", xlab="Time")
lines(0:150, apply(storeMatrix.uni.varG.R, 1, mean), col="green", lwd=1.75, type="l", ylim=c(0,150), ylab="N", xlab="Time")
legend("topright",legend=c("S","I","R"),fill=c("blue","red","green"), cex=0.25)

plot(0:150, apply(storeMatrix.varG.I, 1, mean), col="red", lwd=1.75, type="l", ylim=c(0,200), ylab="N", xlab="Time", main="Strat Var")
lines(0:150, apply(storeMatrix.varG.S, 1, mean), col="blue", lwd=1.75, type="l", ylim=c(0,150), ylab="N", xlab="Time")
lines(0:150, apply(storeMatrix.varG.R, 1, mean), col="green", lwd=1.75, type="l", ylim=c(0,150), ylab="N", xlab="Time")
legend("topright",legend=c("S","I","R"),fill=c("blue","red","green"), cex=0.25)

plot(0:150, apply(storeMatrix.cont.varG.I, 1, mean), col="red", lwd=1.75, type="l", ylim=c(0,200), ylab="N", xlab="Time", main="Cont. Var")
lines(0:150, apply(storeMatrix.cont.varG.S, 1, mean), col="blue", lwd=1.75, type="l", ylim=c(0,150), ylab="N", xlab="Time")
lines(0:150, apply(storeMatrix.cont.varG.R, 1, mean), col="green", lwd=1.75, type="l", ylim=c(0,150), ylab="N", xlab="Time")
legend("topright",legend=c("S","I","R"),fill=c("blue","red","green"), cex=0.25)

## no var GEM
plot(0:150, apply(storeMatrix.no.var.I, 1, mean), col="red", lwd=1.75, type="l", ylim=c(0,200), ylab="N", xlab="Time", main="No Var")
lines(0:150, apply(storeMatrix.no.var.S, 1, mean), col="blue", lwd=1.75, type="l", ylim=c(0,150), ylab="N", xlab="Time")
lines(0:150, apply(storeMatrix.no.var.R, 1, mean), col="green", lwd=1.75, type="l", ylim=c(0,150), ylab="N", xlab="Time")
legend("topright",legend=c("S","I","R"),fill=c("blue","red","green"), cex=0.25)



## CONTINUOUS ALPHA
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

### Strat Var alpha
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

## uniform var alpha
## infected
timeSeq <- 0:150
storeMatrix.uni.alpha.I <- array(NA, dim=c(length(timeSeq),length(out_uniform_var_alpha)))
for (j in 1:length(out_uniform_var_alpha)) {
  o <- out_uniform_var_alpha[[j]]
  storeMatrix.uni.alpha.I[,j] <- o[,3] # num infected
}

## susceptible
storeMatrix.uni.alpha.S <- array(NA, dim=c(length(timeSeq),length(out_uniform_var_alpha)))
for (j in 1:length(out_uniform_var_alpha)) {
  o <- out_uniform_var_alpha[[j]]
  storeMatrix.uni.alpha.S[,j] <- o[,2] # num susceptible
}

## recovered
storeMatrix.uni.alpha.R <- array(NA, dim=c(length(timeSeq),length(out_uniform_var_alpha)))
for (j in 1:length(out_uniform_var_alpha)) {
  o <- out_uniform_var_alpha[[j]]
  storeMatrix.uni.alpha.R[,j] <- o[,4] # num recovered
}

par(mfrow=c(1,4))
plot(0:150, apply(storeMatrix.cont.varA.I, 1, mean), col="red", lwd=1.75, type="l", ylim=c(0,200), ylab="N", xlab="Time", main="Uniform Var")
lines(0:150, apply(storeMatrix.cont.varA.S, 1, mean), col="blue", lwd=1.75, type="l", ylim=c(0,150), ylab="N", xlab="Time")
lines(0:150, apply(storeMatrix.cont.varA.R, 1, mean), col="green", lwd=1.75, type="l", ylim=c(0,150), ylab="N", xlab="Time")
legend("topright",legend=c("S","I","R"),fill=c("blue","red","green"), cex=0.25)

plot(0:150, apply(storeMatrix.alpha.I, 1, mean), col="red", lwd=1.75, type="l", ylim=c(0,200), ylab="N", xlab="Time", main="Strat Var")
lines(0:150, apply(storeMatrix.alpha.S, 1, mean), col="blue", lwd=1.75, type="l", ylim=c(0,150), ylab="N", xlab="Time")
lines(0:150, apply(storeMatrix.alpha.R, 1, mean), col="green", lwd=1.75, type="l", ylim=c(0,150), ylab="N", xlab="Time")
legend("topright",legend=c("S","I","R"),fill=c("blue","red","green"), cex=0.25)

plot(0:150, apply(storeMatrix.uni.alpha.I, 1, mean), col="red", lwd=1.75, type="l", ylim=c(0,200), ylab="N", xlab="Time", main="Cont. Var")
lines(0:150, apply(storeMatrix.uni.alpha.S, 1, mean), col="blue", lwd=1.75, type="l", ylim=c(0,150), ylab="N", xlab="Time")
lines(0:150, apply(storeMatrix.uni.alpha.R, 1, mean), col="green", lwd=1.75, type="l", ylim=c(0,150), ylab="N", xlab="Time")
legend("topright",legend=c("S","I","R"),fill=c("blue","red","green"), cex=0.25)

## no var GEM
plot(0:150, apply(storeMatrix.no.var.I, 1, mean), col="red", lwd=1.75, type="l", ylim=c(0,200), ylab="N", xlab="Time", main="No Var")
lines(0:150, apply(storeMatrix.no.var.S, 1, mean), col="blue", lwd=1.75, type="l", ylim=c(0,150), ylab="N", xlab="Time")
lines(0:150, apply(storeMatrix.no.var.R, 1, mean), col="green", lwd=1.75, type="l", ylim=c(0,150), ylab="N", xlab="Time")
legend("topright",legend=c("S","I","R"),fill=c("blue","red","green"), cex=0.25)


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

### variation in shedding
## infected
timeSeq <- 0:150
storeMatrix.uni.varS.I <- array(NA, dim=c(length(timeSeq),length(out_uniform_var_shed)))
for (j in 1:length(out_uniform_var_shed)) {
  o <- out_uniform_var_shed[[j]]
  storeMatrix.uni.varS.I[,j] <- o[,3] # num infected
}

## susceptible
storeMatrix.uni.varS.S <- array(NA, dim=c(length(timeSeq),length(out_uniform_var_shed)))
for (j in 1:length(out_uniform_var_shed)) {
  o <- out_uniform_var_shed[[j]]
  storeMatrix.uni.varS.S[,j] <- o[,2] # num susceptible
}

## recovered
storeMatrix.uni.varS.R <- array(NA, dim=c(length(timeSeq),length(out_uniform_var_shed)))
for (j in 1:length(out_uniform_var_shed)) {
  o <- out_uniform_var_shed[[j]]
  storeMatrix.uni.varS.R[,j] <- o[,4] # num recovered
}

### continuous var shed
timeSeq <- 0:150
storeMatrix.cont.varS.I <- array(NA, dim=c(length(timeSeq),length(out_cont_var_shed)))
for (j in 1:length(out_cont_var_shed)) {
  o <- out_cont_var_shed[[j]]
  storeMatrix.cont.varS.I[,j] <- o[,3] # num infected
}

## susceptible
storeMatrix.cont.varS.S <- array(NA, dim=c(length(timeSeq),length(out_cont_var_shed)))
for (j in 1:length(out_cont_var_shed)) {
  o <- out_cont_var_shed[[j]]
  storeMatrix.cont.varS.S[,j] <- o[,2] # num susceptible
}

## recovered
storeMatrix.cont.varS.R <- array(NA, dim=c(length(timeSeq),length(out_cont_var_shed)))
for (j in 1:length(out_cont_var_shed)) {
  o <- out_cont_var_shed[[j]]
  storeMatrix.cont.varS.R[,j] <- o[,4] # num recovered
}

### strat var shed
timeSeq <- 0:150
storeMatrix.strat.varS.I <- array(NA, dim=c(length(timeSeq),length(out_strat_var_shed)))
for (j in 1:length(out_strat_var_shed)) {
  o <- out_strat_var_shed[[j]]
  storeMatrix.strat.varS.I[,j] <- o[,3] # num infected
}

## susceptible
storeMatrix.strat.varS.S <- array(NA, dim=c(length(timeSeq),length(out_strat_var_shed)))
for (j in 1:length(out_strat_var_shed)) {
  o <- out_strat_var_shed[[j]]
  storeMatrix.strat.varS.S[,j] <- o[,2] # num susceptible
}

## recovered
storeMatrix.strat.varS.R <- array(NA, dim=c(length(timeSeq),length(out_strat_var_shed)))
for (j in 1:length(out_strat_var_shed)) {
  o <- out_strat_var_shed[[j]]
  storeMatrix.strat.varS.R[,j] <- o[,4] # num recovered
}

par(mfrow=c(1,4))
plot(0:150, apply(storeMatrix.cont.varS.I, 1, mean), col="red", lwd=1.75, type="l", ylim=c(0,200), ylab="N", xlab="Time", main="Uniform Var")
lines(0:150, apply(storeMatrix.cont.varS.S, 1, mean), col="blue", lwd=1.75, type="l", ylim=c(0,150), ylab="N", xlab="Time")
lines(0:150, apply(storeMatrix.cont.varS.R, 1, mean), col="green", lwd=1.75, type="l", ylim=c(0,150), ylab="N", xlab="Time")
legend("topright",legend=c("S","I","R"),fill=c("blue","red","green"), cex=0.25)

plot(0:150, apply(storeMatrix.strat.varS.I, 1, mean), col="red", lwd=1.75, type="l", ylim=c(0,200), ylab="N", xlab="Time", main="Strat Var")
lines(0:150, apply(storeMatrix.strat.varS.S, 1, mean), col="blue", lwd=1.75, type="l", ylim=c(0,150), ylab="N", xlab="Time")
lines(0:150, apply(storeMatrix.strat.varS.R, 1, mean), col="green", lwd=1.75, type="l", ylim=c(0,150), ylab="N", xlab="Time")
legend("topright",legend=c("S","I","R"),fill=c("blue","red","green"), cex=0.25)

plot(0:150, apply(storeMatrix.uni.varS.I, 1, mean), col="red", lwd=1.75, type="l", ylim=c(0,200), ylab="N", xlab="Time", main="Cont. Var")
lines(0:150, apply(storeMatrix.uni.varS.S, 1, mean), col="blue", lwd=1.75, type="l", ylim=c(0,150), ylab="N", xlab="Time")
lines(0:150, apply(storeMatrix.uni.varS.R, 1, mean), col="green", lwd=1.75, type="l", ylim=c(0,150), ylab="N", xlab="Time")
legend("topright",legend=c("S","I","R"),fill=c("blue","red","green"), cex=0.25)

## no var GEM
plot(0:150, apply(storeMatrix.no.var.I, 1, mean), col="red", lwd=1.75, type="l", ylim=c(0,200), ylab="N", xlab="Time", main="No Var")
lines(0:150, apply(storeMatrix.no.var.S, 1, mean), col="blue", lwd=1.75, type="l", ylim=c(0,150), ylab="N", xlab="Time")
lines(0:150, apply(storeMatrix.no.var.R, 1, mean), col="green", lwd=1.75, type="l", ylim=c(0,150), ylab="N", xlab="Time")
legend("topright",legend=c("S","I","R"),fill=c("blue","red","green"), cex=0.25)

