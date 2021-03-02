####### Deterministic model, something dumb

library(deSolve)

deterministic.SIR=function(t,x,params){
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

out=lsoda(xstart,times, deterministic.SIR, param)

S=out[,2]
I1=out[,3]
I2=out[,4]
I3=out[,5]
R=out[,6]

par(mfrow=c(1,1))
plot.ts(S, col="blue", ylim=c(0,200), ylab="N", xlab="Time")
lines(I1, col="red")
lines(I2, col="red")
lines(I3, col="red")
lines(R, col="green")
legend("topright",legend=c("S","I","R"),fill=c("blue","red","green"), cex=.7)
