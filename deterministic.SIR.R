####### Deterministic model

library(deSolve)

deterministic.SIR=function(t,x,params){
  beta = params["beta"]
  alpha = params["alpha"]
  gamma = params["gamma"]
  b = params["b"]
  bs = params["bs"]
  d = params["d"]
  ds = params["ds"]
  S=x[1]
  I=x[2]
  R=x[3]
  
  dS=-beta*S*I + ((b-bs*(S+I+R))*(S+I+R)) - (d*(S))
  dI= (beta*S*I)-(I*(alpha+gamma)) - (d*(I))
  dR= gamma*I - (d*(R))
  
  list(c(dS,dI,dR))
  
}

# parameter values
param <-c(beta=.1, alpha=.01, gamma = .3, b=2, d=0.4, bs=0.01, ds=0.1)
times = seq(0,100,by=1) # time steps to output 
xstart = c(70,10,0) # beginning population size

out=lsoda(xstart,times, deterministic.SIR, param)

S=out[,2]
I=out[,3]
R=out[,4]

plot.ts(S, col="red", ylim=c(0,130), ylab="N", xlab="Time")
lines(I, col="blue")
lines(R, col="green")
legend("topright",legend=c("S","I","R"),fill=c("red","blue","green"), cex=.7)
