####### Deterministic model, something dumb

library(deSolve)

deterministic.SIR=function(t,x,params){
  beta = params["beta"]
  alpha = params["alpha"]
  gamma = params["gamma"]
  c = params["c"]
  shed = params["shed"]
  h = params["h"]
  b = params["b"]
  bs = params["bs"]
  d = params["d"]
  ds = params["ds"]
  S=x[1]
  I=x[2]
  R=x[3]
  
  #dS=-c*(shed/(h+shed))*S*I + ((b-bs*(S+I+R))*(S+I+R)) - (d*(S))
  #dI= c*(shed/(h+shed))*S*I-(I*(alpha+gamma)) - (d*(I))
  #dR= gamma*I - (d*(R))
  
  dS  = -S*I*beta + ((b-bs*(S+I+R))*(S+I+R)) - (d*(S))
  dI  = S*I*beta - I(alpha+gamma) - d*(I)
  dR  = gamma*(I) - d*R
  
  list(c(dS,dI,dR))
  
}
param <-c(beta=.25, alpha=.15, gamma = .15, b=2.5, d=0.4, bs=0.01, ds=0.1)
times = seq(0,100,by=1) # time steps to output 
xstart = c(70,10,0) # beginning population size

out.det.novar=lsoda(xstart,times, deterministic.SIR, param)

S.det.novar=out.det.novar[,2]
I.det.novar=out.det.novar[,3]
R.det.novar=out.det.novar[,4]

par(mfrow=c(2,2))
plot.ts(S.det.novar, col="blue", ylim=c(0,200), ylab="N", xlab="Time",main="Det. No Var")
lines(I.det.novar, col="red")
lines(R.det.novar, col="green")
legend("topright",legend=c("S","I","R"),fill=c("blue","red","green"), cex=.25)


deterministic.SIR_stratB=function(t,x,params){
  beta1 = params["beta1"]
  beta2 = params["beta2"]
  beta3 = params["beta3"]
  alpha = params["alpha"]
  gamma = params["gamma"]
  c = params["c"]
  shed = params["shed"]
  h = params["h"]
  b = params["b"]
  bs = params["bs"]
  d = params["d"]
  ds = params["ds"]
  S=x[1]
  I1=x[2]
  I2=x[3]
  I3=x[4]
  R=x[5]
  
  #dS=-c*(shed/(h+shed))*S*I + ((b-bs*(S+I+R))*(S+I+R)) - (d*(S))
  #dI= c*(shed/(h+shed))*S*I-(I*(alpha+gamma)) - (d*(I))
  #dR= gamma*I - (d*(R))
  
  dS= -(S*((beta1*I1)+(beta2*I1)+(beta3*I3)))+ ((b-bs*(S+I1+I2+I3+R))*(S+I1+I2+I3+R)) - (d*(S))
  dI1 = (.333)*(S*((beta1*I1)+(beta2*I1)+(beta3*I3))) - alpha*(I1)-gamma*(I1) - d*(I1)
  dI2 = (.333)*(S*((beta1*I1)+(beta2*I1)+(beta3*I3))) - alpha*(I2)-gamma*(I2)- d*(I2)
  dI3 = (.333)*(S*((beta1*I1)+(beta2*I1)+(beta3*I3))) - alpha*(I3)-gamma*(I3)- d*(I3)
  dR = gamma*(I1+I2+I3) - d*R
  
  list(c(dS,dI1,dI2, dI3,dR))
  
}

# parameter values
param <-c(beta1=.15, beta2=.25,beta3=.35, alpha=.15, gamma = .15, b=2.5, d=0.4, bs=0.01, ds=0.1)

params = c(c=.2, shed=.2, sd_s=.1, sd_a=.1, sd_c=.01, sd_g=.01, 
           h=.1, alpha=.01, gamma=.3, b=2.5, d=.4, bs=.01)
params7 = c(c=.07, shed=.07, sd_s=.25, sd_a=.25, sd_c=.25, sd_g=.25, 
            h=.13, alpha=.1, gamma=.15, b=2.5, d=.4, bs=.01)
params6 = c(c=.05, shed=.05, sd_s=.05, sd_a=.1, sd_c=.05, sd_g=.15, 
            h=.15, alpha=.1, gamma=.15, b=2.5, d=.4, bs=.01)

times = seq(0,100,by=1) # time steps to output 
xstart = c(70,4,4,4,0) # beginning population size

out.strat.b=lsoda(xstart,times, deterministic.SIR, param)

S.strat.b=out.strat.b[,2]
I1.strat.b=out.strat.b[,3]
I2.strat.b=out.strat.b[,4]
I3.strat.b=out.strat.b[,5]
R.strat.b=out.strat.b[,6]

par(mfrow=c(1,1))
plot.ts(S.strat.b, col="blue", ylim=c(0,200), ylab="N", xlab="Time",main="Det. Strat. Beta")
lines(I1.strat.b*3, col="red")
lines(R.strat.b, col="green")
legend("topright",legend=c("S","I","R"),fill=c("blue","red","green"), cex=.25)

deterministic.SIR_stratG=function(t,x,params){
  gamma1 = params["gamma1"]
  gamma2 = params["gamma2"]
  gamma3 = params["gamma3"]
  beta = params["beta"]
  alpha = params["alpha"]
  c = params["c"]
  shed = params["shed"]
  h = params["h"]
  b = params["b"]
  bs = params["bs"]
  d = params["d"]
  ds = params["ds"]
  S=x[1]
  I1=x[2]
  I2=x[3]
  I3=x[4]
  R=x[5]

  
  dS= -(S*(I1+I2+I3)*beta)+ ((b-bs*(S+I1+I2+I3+R))*(S+I1+I2+I3+R)) - (d*(S))
  dI1 = (.333)*((S*beta*I1) - alpha*(I1)-((gamma1*I1)+(gamma2*I2)+(gamma3*I3)) - d*(I1))
  dI2 = (.333)*((S*beta*I2) - alpha*(I2)-((gamma1*I1)+(gamma2*I2)+(gamma3*I3))- d*(I2))
  dI3 = (.333)*((S*beta*I3) - alpha*(I3)-((gamma1*I1)+(gamma2*I2)+(gamma3*I3))- d*(I3))
  dR = ((gamma1*I1)+(gamma2*I2)+(gamma3*I3)) - d*R
  
  list(c(dS,dI1,dI2, dI3,dR))
  
}
params<-c(beta=.25, gamma1=.05, gamma2=.15, gamma3=.35,alpha=.15,b=2.5, d=0.4, bs=0.01, ds=0.1)

times = seq(0,100,by=1) # time steps to output 
xstart = c(70,4,4,4,0) # beginning population size

out.det.strat.g=lsoda(xstart,times, deterministic.SIR_stratG, params)

S.strat.g=out.det.strat.g[,2]
I1.strat.g=out.det.strat.g[,3]
I2.strat.g=out.det.strat.g[,4]
I3.strat.g=out.det.strat.g[,5]
R.strat.g=out.det.strat.g[,6]

plot.ts(S.strat.g, col="blue", ylim=c(0,200), ylab="N", xlab="Time", main="Det. Strat. Gamma")
lines((I1.strat.g*3), col="red")
lines(R.strat.g, col="green")
legend("topright",legend=c("S","I","R"),fill=c("blue","red","green"), cex=.25)

deterministic.SIR_stratA=function(t,x,params){
  alpha1 = params["alpha1"]
  alpha2 = params["alpha2"]
  alpha3 = params["alpha3"]
  beta = params["beta"]
  gamma = params["gamma"]
  c = params["c"]
  shed = params["shed"]
  h = params["h"]
  b = params["b"]
  bs = params["bs"]
  d = params["d"]
  ds = params["ds"]
  S=x[1]
  I1=x[2]
  I2=x[3]
  I3=x[4]
  R=x[5]
  
  
  dS= -(S*(I1+I2+I3)*beta)+ ((b-bs*(S+I1+I2+I3+R))*(S+I1+I2+I3+R)) - (d*(S))
  dI1 = (.333)*((S*beta*I1) - ((alpha1*I1)+(alpha2*I2)+(alpha3*I3))-(gamma*I1) - d*(I1))
  dI2 = (.333)*((S*beta*I2) - ((alpha1*I1)+(alpha2*I2)+(alpha3*I3))-(gamma*I2) - d*(I2))
  dI3 = (.333)*((S*beta*I3) - ((alpha1*I1)+(alpha2*I2)+(alpha3*I3))-(gamma*I3) - d*(I3))
  dR =  (gamma*(I1+I2+I3))- d*R
  
  list(c(dS,dI1,dI2, dI3,dR))
  
}

params<-c(beta=.25, alpha1=.05, alpha2=.15, alpha3=.35,b=2.5, d=0.4, bs=0.01, ds=0.1, gamma=.15)

times = seq(0,100,by=1) # time steps to output 
xstart = c(70,4,4,4,0) # beginning population size

out.det.strat.a=lsoda(xstart,times, deterministic.SIR_stratA, params)

S.strat.a=out.det.strat.a[,2]
I1.strat.a=out.det.strat.a[,3]
I2.strat.a=out.det.strat.a[,4]
I3.strat.a=out.det.strat.a[,5]
R.strat.a=out.det.strat.a[,6]

plot.ts(S.strat.a, col="blue", ylim=c(0,200), ylab="N", xlab="Time", main="Det. Strat. Alpha")
lines((I1.strat.a*3), col="red")
lines(R.strat.a, col="green")
legend("topright",legend=c("S","I","R"),fill=c("blue","red","green"), cex=.25)


