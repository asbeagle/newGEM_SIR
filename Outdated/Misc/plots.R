## 
params4all<- c(beta=.01, alpha=.010, gamma=.2,  varB=.001,varG=.001,varA=.001,b=2,bs=.01,d=.4)

# var G, num susceptible
timeSeq <- 0:100
storeMatrix.varG.S <- array(NA, dim=c(length(timeSeq),length(outvarG)))
for (j in 1:length(outvarG)) {
  o <- outvarG[[j]]
  recordTimes <- lapply(o, function(o2) o2[[1]]) %>% unlist
  storeVector <- rep(0, length=length(timeSeq))
  listInds <- sapply(timeSeq, function(t) which(recordTimes >= t) %>% min)
  if (any(is.infinite(listInds)))
    listInds <- listInds[-which(is.infinite(listInds))] 
  storeVector[1:length(listInds)] <- sapply(listInds, function(i) (o[[i]][[2]])) ## number of susceptible
  storeMatrix.varG.S[,j] <- storeVector
}

## var G num infected
timeSeq <- 0:100
storeMatrix.varG.I <- array(NA, dim=c(length(timeSeq),length(outvarG)))
for (j in 1:length(outvarG)) {
  o <- outvarG[[j]]
  recordTimes <- lapply(o, function(o2) o2[[1]]) %>% unlist
  storeVector <- rep(0, length=length(timeSeq))
  listInds <- sapply(timeSeq, function(t) which(recordTimes >= t) %>% min)
  if (any(is.infinite(listInds)))
    listInds <- listInds[-which(is.infinite(listInds))] 
  storeVector[1:length(listInds)] <- sapply(listInds, function(i) length((o[[i]][[3]]))) ## number of susceptible
  storeMatrix.varG.I[,j] <- storeVector
}

## var G num recovered
timeSeq <- 0:100
storeMatrix.varG.R <- array(NA, dim=c(length(timeSeq),length(outvarG)))
for (j in 1:length(outvarG)) {
  o <- outvarG[[j]]
  recordTimes <- lapply(o, function(o2) o2[[1]]) %>% unlist
  storeVector <- rep(0, length=length(timeSeq))
  listInds <- sapply(timeSeq, function(t) which(recordTimes >= t) %>% min)
  if (any(is.infinite(listInds)))
    listInds <- listInds[-which(is.infinite(listInds))] 
  storeVector[1:length(listInds)] <- sapply(listInds, function(i) (o[[i]][[4]])) ## number of susceptible
  storeMatrix.varG.R[,j] <- storeVector
}

### beta variation
library(parallel)
mclapply(seeds,
         function(s) gillespie.SIR.varB(tmax, params4all, x, s),
         mc.cores=4) -> outvarB

# num infected
timeSeq <- 0:100
storeMatrix.varB.I <- array(NA, dim=c(length(timeSeq),length(outvarB)))
for (j in 1:length(outvarB)) {
  o <- outvarB[[j]] # list of the seeds
  recordTimes <- lapply(o, function(o2) o2[[1]]) %>% unlist
  storeVector <- rep(0, length=length(timeSeq))
  listInds <- sapply(timeSeq, function(t) which(recordTimes >= t) %>% min) #rows=time steps, columns=seeds, value of specific run at this time step
  if (any(is.infinite(listInds)))
    listInds <- listInds[-which(is.infinite(listInds))] 
  storeVector[1:length(listInds)] <- sapply(listInds, function(i) length(o[[i]][[3]]))
  storeMatrix.varB.I[,j] <- storeVector
}

# num susceptible
timeSeq <- 0:100
storeMatrix.varB.S <- array(NA, dim=c(length(timeSeq),length(outvarB)))
for (j in 1:length(outvarB)) {
  o <- outvarB[[j]] # list of the seeds
  recordTimes <- lapply(o, function(o2) o2[[1]]) %>% unlist
  storeVector <- rep(0, length=length(timeSeq))
  listInds <- sapply(timeSeq, function(t) which(recordTimes >= t) %>% min) #rows=time steps, columns=seeds, value of specific run at this time step
  if (any(is.infinite(listInds)))
    listInds <- listInds[-which(is.infinite(listInds))] 
  storeVector[1:length(listInds)] <- sapply(listInds, function(i) (o[[i]][[2]])) ## mean trait value for infected at this time step for this seed
  storeMatrix.varB.S[,j] <- storeVector
}

# num recovered
timeSeq <- 0:100
storeMatrix.varB.R <- array(NA, dim=c(length(timeSeq),length(outvarB)))
for (j in 1:length(outvarB)) {
  o <- outvarB[[j]] # list of the seeds
  recordTimes <- lapply(o, function(o2) o2[[1]]) %>% unlist
  storeVector <- rep(0, length=length(timeSeq))
  listInds <- sapply(timeSeq, function(t) which(recordTimes >= t) %>% min) #rows=time steps, columns=seeds, value of specific run at this time step
  if (any(is.infinite(listInds)))
    listInds <- listInds[-which(is.infinite(listInds))] 
  storeVector[1:length(listInds)] <- sapply(listInds, function(i) (o[[i]][[4]])) ## mean trait value for infected at this time step for this seed
  storeMatrix.varB.R[,j] <- storeVector
}

### alpha variation

## number of S
timeSeq <- 0:100
storeMatrix.varA.S <- array(NA, dim=c(length(timeSeq),length(outvarA)))
for (j in 1:length(outvarA)) {
  o <- outvarA[[j]]
  recordTimes <- lapply(o, function(o2) o2[[1]]) %>% unlist
  storeVector <- rep(0, length=length(timeSeq))
  listInds <- sapply(timeSeq, function(t) which(recordTimes >= t) %>% min)
  if (any(is.infinite(listInds)))
    listInds <- listInds[-which(is.infinite(listInds))] 
  storeVector[1:length(listInds)] <- sapply(listInds, function(i) (o[[i]][[2]])) ## number of infecteds
  storeMatrix.varA.S[,j] <- storeVector
}

## num of I
timeSeq <- 0:100
storeMatrix.varA.I <- array(NA, dim=c(length(timeSeq),length(outvarA)))
for (j in 1:length(outvarA)) {
  o <- outvarA[[j]]
  recordTimes <- lapply(o, function(o2) o2[[1]]) %>% unlist
  storeVector <- rep(0, length=length(timeSeq))
  listInds <- sapply(timeSeq, function(t) which(recordTimes >= t) %>% min)
  if (any(is.infinite(listInds)))
    listInds <- listInds[-which(is.infinite(listInds))] 
  storeVector[1:length(listInds)] <- sapply(listInds, function(i) length((o[[i]][[3]]))) ## number of infected
  storeMatrix.varA.I[,j] <- storeVector
}

# num of R
timeSeq <- 0:100
storeMatrix.varA.R <- array(NA, dim=c(length(timeSeq),length(outvarA)))
for (j in 1:length(outvarA)) {
  o <- outvarA[[j]]
  recordTimes <- lapply(o, function(o2) o2[[1]]) %>% unlist
  storeVector <- rep(0, length=length(timeSeq))
  listInds <- sapply(timeSeq, function(t) which(recordTimes >= t) %>% min)
  if (any(is.infinite(listInds)))
    listInds <- listInds[-which(is.infinite(listInds))] 
  storeVector[1:length(listInds)] <- sapply(listInds, function(i) (o[[i]][[4]])) ## number of infected
  storeMatrix.varA.R[,j] <- storeVector
}

### plots

## mort+alpha+gamma > beta

par(mfrow=c(3,3))

# susceptible gamma
plot(0:100, apply(storeMatrix.varG.S, 1, mean), col="red", lwd=1.75, type="l", ylim=c(0,150), ylab="N", xlab="Time", main="Gamma - Susceptible")
for (i in 1:ncol(storeMatrix.varG.S)) {
  lines(0:100, storeMatrix.varG.S[,i], col=grey(.1, alpha=.37))
}

# infected gamma
plot(0:100, apply(storeMatrix.varG.I, 1, mean), col="red", lwd=1.75, type="l", ylim=c(0,150), ylab="N", xlab="Time", main="Gamma - Infected")
for (i in 1:ncol(storeMatrix.varG.I)) {
  lines(0:100, storeMatrix.varG.I[,i], col=grey(.1, alpha=.37))
}

# recovered gamma
plot(0:100, apply(storeMatrix.varG.R, 1, mean), col="red", lwd=1.75, type="l", ylim=c(0,150), ylab="N", xlab="Time", main="Gamma - Recovered")
for (i in 1:ncol(storeMatrix.varG.R)) {
  lines(0:100, storeMatrix.varG.R[,i], col=grey(.1, alpha=.37))
}

# susceptible beta
plot(0:100, apply(storeMatrix.varB.S, 1, mean), col="red", lwd=1.75, type="l", ylim=c(0,200), ylab="N", xlab="Time", main="Beta - Susceptible")
for (i in 1:ncol(storeMatrix.varB.S)) {
  lines(0:100, storeMatrix.varB.S[,i], col=grey(.1, alpha=.37))
}

# infected beta
plot(0:100, apply(storeMatrix.varB.I, 1, mean), col="red", lwd=1.75, type="l", ylim=c(0,200), ylab="N", xlab="Time", main="Beta - Infected")
for (i in 1:ncol(storeMatrix.varB.I)) {
  lines(0:100, storeMatrix.varB.I[,i], col=grey(.1, alpha=.37))
}

# recovered beta
plot(0:100, apply(storeMatrix.varB.R, 1, mean), col="red", lwd=1.75, type="l", ylim=c(0,150), ylab="N", xlab="Time", main="Beta - Recovered")
for (i in 1:ncol(storeMatrix.varB.R)) {
  lines(0:100, storeMatrix.varB.R[,i], col=grey(.1, alpha=.37))
}

# susceptible alpha
plot(0:100, apply(storeMatrix.varA.S, 1, mean), col="red", lwd=1.75, type="l", ylim=c(0,150), ylab="N", xlab="Time", main="Alpha - Susceptible")
for (i in 1:ncol(storeMatrix.varA.S)) {
  lines(0:100, storeMatrix.varA.S[,i], col=grey(.1, alpha=.37))
}

# infected alpha
plot(0:100, apply(storeMatrix.varA.I, 1, mean), col="red", lwd=1.75, type="l", ylim=c(0,150), ylab="N", xlab="Time", main="Alpha - Infected")
for (i in 1:ncol(storeMatrix.varA.I)) {
  lines(0:100, storeMatrix.varA.I[,i], col=grey(.1, alpha=.37))
}

# recovered alpha
plot(0:100, apply(storeMatrix.varA.R, 1, mean),col="red", lwd=1.75, type="l", ylim=c(0,150), ylab="N", xlab="Time", main="Alpha - Recovered")
for (i in 1:ncol(storeMatrix.varA.R)) {
  lines(0:100, storeMatrix.varA.R[,i], col=grey(.1, alpha=.37))
}



#### change param values m + alpha + gamma < beta
params4all2<-c(beta=.5, alpha=.01, gamma=.2, varB=1e-3, varG=1e-3, varA=1e-3, b=2, bs=.01, d=.1)
x = c(S=70, I=10, R=0)
tmax <- 100

# gamma
# run multiple times
library(tidyverse)
seeds <- floor(runif(20,1,1e5)) # set seeds
library(parallel)
detectCores() # detect how many cores there are

mclapply(seeds,
         function(s) gillespie.SIR.varG(tmax, params4all2, x, s),
         mc.cores=4) -> outvarG2

mclapply(seeds,
         function(s) gillespie.SIR.varB(tmax, params4all2, x, s),
         mc.cores=4) -> outvarB2

mclapply(seeds,
         function(s) gillespie.SIR.varA(tmax, params4all2, x, s),
         mc.cores=4) -> outvarA2

### number S,I,R
## gamma variation

# var G, num susceptible
timeSeq <- 0:100
storeMatrix.varG.S2 <- array(NA, dim=c(length(timeSeq),length(outvarG2)))
for (j in 1:length(outvarG2)) {
  o <- outvarG2[[j]]
  recordTimes <- lapply(o, function(o2) o2[[1]]) %>% unlist
  storeVector <- rep(0, length=length(timeSeq))
  listInds <- sapply(timeSeq, function(t) which(recordTimes >= t) %>% min)
  if (any(is.infinite(listInds)))
    listInds <- listInds[-which(is.infinite(listInds))] 
  storeVector[1:length(listInds)] <- sapply(listInds, function(i) (o[[i]][[2]])) ## number of susceptible
  storeMatrix.varG.S2[,j] <- storeVector
}


## var G num infected
timeSeq <- 0:100
storeMatrix.varG.I2 <- array(NA, dim=c(length(timeSeq),length(outvarG2)))
for (j in 1:length(outvarG2)) {
  o <- outvarG2[[j]]
  recordTimes <- lapply(o, function(o2) o2[[1]]) %>% unlist
  storeVector <- rep(0, length=length(timeSeq))
  listInds <- sapply(timeSeq, function(t) which(recordTimes >= t) %>% min)
  if (any(is.infinite(listInds)))
    listInds <- listInds[-which(is.infinite(listInds))] 
  storeVector[1:length(listInds)] <- sapply(listInds, function(i) length((o[[i]][[3]]))) ## number of susceptible
  storeMatrix.varG.I2[,j] <- storeVector
}

## var G num recovered
timeSeq <- 0:100
storeMatrix.varG.R2 <- array(NA, dim=c(length(timeSeq),length(outvarG2)))
for (j in 1:length(outvarG2)) {
  o <- outvarG2[[j]]
  recordTimes <- lapply(o, function(o2) o2[[1]]) %>% unlist
  storeVector <- rep(0, length=length(timeSeq))
  listInds <- sapply(timeSeq, function(t) which(recordTimes >= t) %>% min)
  if (any(is.infinite(listInds)))
    listInds <- listInds[-which(is.infinite(listInds))] 
  storeVector[1:length(listInds)] <- sapply(listInds, function(i) (o[[i]][[4]])) ## number of susceptible
  storeMatrix.varG.R2[,j] <- storeVector
}

### beta variation

# num infected
timeSeq <- 0:100
storeMatrix.varB.I2 <- array(NA, dim=c(length(timeSeq),length(outvarB2)))
for (j in 1:length(outvarB2)) {
  o <- outvarB2[[j]] # list of the seeds
  recordTimes <- lapply(o, function(o2) o2[[1]]) %>% unlist
  storeVector <- rep(0, length=length(timeSeq))
  listInds <- sapply(timeSeq, function(t) which(recordTimes >= t) %>% min) #rows=time steps, columns=seeds, value of specific run at this time step
  if (any(is.infinite(listInds)))
    listInds <- listInds[-which(is.infinite(listInds))] 
  storeVector[1:length(listInds)] <- sapply(listInds, function(i) length(o[[i]][[3]])) ## mean trait value for infected at this time step for this seed
  storeMatrix.varB.I2[,j] <- storeVector
}

# num susceptible
timeSeq <- 0:100
storeMatrix.varB.S2 <- array(NA, dim=c(length(timeSeq),length(outvarB2)))
for (j in 1:length(outvarB2)) {
  o <- outvarB2[[j]] # list of the seeds
  recordTimes <- lapply(o, function(o2) o2[[1]]) %>% unlist
  storeVector <- rep(0, length=length(timeSeq))
  listInds <- sapply(timeSeq, function(t) which(recordTimes >= t) %>% min) #rows=time steps, columns=seeds, value of specific run at this time step
  if (any(is.infinite(listInds)))
    listInds <- listInds[-which(is.infinite(listInds))] 
  storeVector[1:length(listInds)] <- sapply(listInds, function(i) (o[[i]][[2]])) ## mean trait value for infected at this time step for this seed
  storeMatrix.varB.S2[,j] <- storeVector
}

# num recovered
timeSeq <- 0:100
storeMatrix.varB.R2 <- array(NA, dim=c(length(timeSeq),length(outvarB2)))
for (j in 1:length(outvarB2)) {
  o <- outvarB2[[j]] # list of the seeds
  recordTimes <- lapply(o, function(o2) o2[[1]]) %>% unlist
  storeVector <- rep(0, length=length(timeSeq))
  listInds <- sapply(timeSeq, function(t) which(recordTimes >= t) %>% min) #rows=time steps, columns=seeds, value of specific run at this time step
  if (any(is.infinite(listInds)))
    listInds <- listInds[-which(is.infinite(listInds))] 
  storeVector[1:length(listInds)] <- sapply(listInds, function(i) (o[[i]][[4]])) ## mean trait value for infected at this time step for this seed
  storeMatrix.varB.R2[,j] <- storeVector
}

### alpha variation

## number of S
timeSeq <- 0:100
storeMatrix.varA.S2 <- array(NA, dim=c(length(timeSeq),length(outvarA2)))
for (j in 1:length(outvarA2)) {
  o <- outvarA2[[j]]
  recordTimes <- lapply(o, function(o2) o2[[1]]) %>% unlist
  storeVector <- rep(0, length=length(timeSeq))
  listInds <- sapply(timeSeq, function(t) which(recordTimes >= t) %>% min)
  if (any(is.infinite(listInds)))
    listInds <- listInds[-which(is.infinite(listInds))] 
  storeVector[1:length(listInds)] <- sapply(listInds, function(i) (o[[i]][[2]])) ## number of infecteds
  storeMatrix.varA.S2[,j] <- storeVector
}

## num of I
timeSeq <- 0:100
storeMatrix.varA.I2 <- array(NA, dim=c(length(timeSeq),length(outvarA2)))
for (j in 1:length(outvarA2)) {
  o <- outvarA2[[j]]
  recordTimes <- lapply(o, function(o2) o2[[1]]) %>% unlist
  storeVector <- rep(0, length=length(timeSeq))
  listInds <- sapply(timeSeq, function(t) which(recordTimes >= t) %>% min)
  if (any(is.infinite(listInds)))
    listInds <- listInds[-which(is.infinite(listInds))] 
  storeVector[1:length(listInds)] <- sapply(listInds, function(i) length((o[[i]][[3]]))) ## number of infected
  storeMatrix.varA.I2[,j] <- storeVector
}

# num of R
timeSeq <- 0:100
storeMatrix.varA.R2 <- array(NA, dim=c(length(timeSeq),length(outvarA2)))
for (j in 1:length(outvarA2)) {
  o <- outvarA2[[j]]
  recordTimes <- lapply(o, function(o2) o2[[1]]) %>% unlist
  storeVector <- rep(0, length=length(timeSeq))
  listInds <- sapply(timeSeq, function(t) which(recordTimes >= t) %>% min)
  if (any(is.infinite(listInds)))
    listInds <- listInds[-which(is.infinite(listInds))] 
  storeVector[1:length(listInds)] <- sapply(listInds, function(i) (o[[i]][[4]])) ## number of infected
  storeMatrix.varA.R2[,j] <- storeVector
}

## plots

par(mfrow=c(3,3))

# susceptible gamma
plot(0:100, apply(storeMatrix.varG.S2, 1, mean), col="red", lwd=1.75, type="l", ylim=c(0,150), ylab="N", xlab="Time", main="Gamma - Susceptible")
for (i in 1:ncol(storeMatrix.varG.S2)) {
  lines(0:100, storeMatrix.varG.S2[,i], col=grey(.1, alpha=.37))
}

# infected gamma
plot(0:100, apply(storeMatrix.varG.I2, 1, mean), col="red", lwd=1.75, type="l", ylim=c(0,150), ylab="N", xlab="Time", main="Gamma - Infected")
for (i in 1:ncol(storeMatrix.varG.I2)) {
  lines(0:100, storeMatrix.varG.I2[,i], col=grey(.1, alpha=.37))
}

# recovered gamma
plot(0:100, apply(storeMatrix.varG.R2, 1, mean), col="red", lwd=1.75, type="l", ylim=c(0,150), ylab="N", xlab="Time", main="Gamma - Recovered")
for (i in 1:ncol(storeMatrix.varG.R2)) {
  lines(0:100, storeMatrix.varG.R2[,i], col=grey(.1, alpha=.37))
}

# susceptible beta
plot(0:100, apply(storeMatrix.varB.S2, 1, mean), col="red", lwd=1.75, type="l", ylim=c(0,150), ylab="N", xlab="Time", main="Beta - Susceptible")
for (i in 1:ncol(storeMatrix.varB.S2)) {
  lines(0:100, storeMatrix.varB.S2[,i], col=grey(.1, alpha=.37))
}

# infected beta
plot(0:100, apply(storeMatrix.varB.I2, 1, mean), col="red", lwd=1.75, type="l", ylim=c(0,150), ylab="N", xlab="Time", main="Beta - Infected")
for (i in 1:ncol(storeMatrix.varB.I2)) {
  lines(0:100, storeMatrix.varB.I2[,i], col=grey(.1, alpha=.37))
}

# recovered beta
plot(0:100, apply(storeMatrix.varB.R2, 1, mean), col="red", lwd=1.75, type="l", ylim=c(0,150), ylab="N", xlab="Time", main="Beta - Recovered")
for (i in 1:ncol(storeMatrix.varB.R2)) {
  lines(0:100, storeMatrix.varB.R2[,i], col=grey(.1, alpha=.37))
}

# susceptible alpha
plot(0:100, apply(storeMatrix.varA.S2, 1, mean), col="red", lwd=1.75, type="l", ylim=c(0,150), ylab="N", xlab="Time", main="Alpha - Susceptible")
for (i in 1:ncol(storeMatrix.varA.S2)) {
  lines(0:100, storeMatrix.varA.S2[,i], col=grey(.1, alpha=.37))
}

# infected alpha
plot(0:100, apply(storeMatrix.varA.I2, 1, mean), col="red", lwd=1.75, type="l", ylim=c(0,150), ylab="N", xlab="Time", main="Alpha - Infected")
for (i in 1:ncol(storeMatrix.varA.I2)) {
  lines(0:100, storeMatrix.varA.I2[,i], col=grey(.1, alpha=.37))
}

# recovered alpha
plot(0:100, apply(storeMatrix.varA.R2, 1, mean), col="red", lwd=1.75, type="l", ylim=c(0,150), ylab="N", xlab="Time", main="Alpha - Recovered")
for (i in 1:ncol(storeMatrix.varA.R2)) {
  lines(0:100, storeMatrix.varA.R2[,i], col=grey(.1, alpha=.37))
}

############ make d+alpha+gamma = beta
library(parallel)
params4all3<-c(beta=0.35, alpha=0.10,gamma= 0.100,varB=0.001,varA= 0.001,varG= 0.001,b=2.000,bs= 0.010,d= 0.15)

mclapply(seeds,
         function(s) gillespie.SIR.varG(tmax, params4all3, x, s),
         mc.cores=4) -> outvarG3

mclapply(seeds,
         function(s) gillespie.SIR.varB(tmax, params4all3, x, s),
         mc.cores=4) -> outvarB3

mclapply(seeds,
         function(s) gillespie.SIR.varA(tmax, params4all3, x, s),
         mc.cores=4) -> outvarA3

### variation plot

var_beta1<-c(rep(0, 20))
var_I1<-c(rep(0, 20))

for (i in 1:length(outvarB)){
  var_beta1[i]<-var(outvarB[[i]][[3]][[3]])
}

for (i in 1:length(outvarB)){
    var_I1<-length(outvarB[[i]][[50]][[3]])
  }

plot(var_beta1)
plot(var_I1)

# box plot of mean and variance 

lapply(outvarB, function(o) length(tail(o,1)[[1]][[3]])) %>% unlist # get the number of infected at the last time step

lapply(outvarB, function(o) length(tail(o,1)[[1]][[3]])) %>% unlist %>% var()

lapply(outvarB, function(o) length(tail(o,1)[[1]][[3]])) %>% unlist -> is

sd(is)/mean(is) ## coefficient of variance 


boxplot(is)
