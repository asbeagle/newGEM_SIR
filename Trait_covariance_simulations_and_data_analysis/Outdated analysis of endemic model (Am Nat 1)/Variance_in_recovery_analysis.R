source("GEM_SIR_cov_storage.R")
library(tidyverse)

x <- 0.24647
hivar = c(   c=(1-x)*0.1,    shed=(1-x)*1/9,    alpha=(1+x)*0.1,    gamma=(1+x)*0.1, 
             sd_c=(1-x)*0.5, sd_shed=(1-x)*5/9, sd_alpha=(1+x)*0.5, sd_gamma=(1+x)*0.5, 
             b=2.5, d=.1, bs=.01)
initial_state <- floor(c(S=unname(((hivar["b"]-hivar["d"])/hivar["bs"]))-5, I=5, R=0))

gillespie.SIR.cov_storage(tmax=100, 
                          params=hivar, 
                          corr=0, 
                          x=initial_state, 
                          covParams=c('c','alpha'),
                          seed=1234123) -> out1

gillespie.SIR.cov_storage(tmax=100, 
                          params=hivar, 
                          corr=0, 
                          x=initial_state, 
                          covParams=c('c','gamma'),
                          seed=1234123) -> out2

par(mfrow=c(1,2))
plot(out1[[1]][,c('t','S')], type='l', ylim=c(0,250))
lines(out1[[1]][,c('t','I')], col=2)
lines(out1[[1]][,c('t','R')], col=4)

plot(out2[[1]][,c('t','S')], type='l', lty=2, ylim=c(0,250))
lines(out2[[1]][,c('t','I')], col=2, lty=2)
lines(out2[[1]][,c('t','R')], col=4, lty=2)



library(deSolve)

diseaseModelVar = function(t,y,parms) {
  b = parms["b"]
  d = parms["d"]
  bs = parms["bs"]
  c = parms["c"]
  shed = parms["shed"]
  gamma = parms["gamma"]
  alpha = parms["alpha"]
  gs = parms["gs"]
  as = parms["as"]
  
  S = y[1]
  I1 = y[2]
  I2 = y[3]
  R = y[4]
  
  dSdt = (b - bs*(S+I1+I2+R))*(S+I1+I2+R) - c*(shed/(1+shed))*S*(I1+I2) - d*S
  dI1dt = 1/2*c*(shed/(1+shed))*S*(I1+I2) - (1+gs)*gamma*I1 - (1+as)*alpha*I1 - d*I1
  dI2dt = 1/2*c*(shed/(1+shed))*S*(I1+I2) - (1-gs)*gamma*I2 - (1-as)*alpha*I2 - d*I2
  dRdt = (1+as)*gamma*I1 + (1-as)*gamma*I2 - d*R
  
  return(list(c(dSdt,dI1dt,dI2dt,dRdt),(b - bs*(S+I1+I2+R))*(S+I1+I2+R),c*(shed/(1+shed))*S*(I1+I2)))
  
}

x <- 0.24647
hivar = c(   c=(1-x)*0.1,    shed=(1-x)*1/9,    alpha=(1+x)*0.1,    gamma=(1+x)*0.1, 
             sd_c=(1-x)*0.5, sd_shed=(1-x)*5/9, sd_alpha=(1+x)*0.5, sd_gamma=(1+x)*0.5, 
             b=2.5, d=.1, bs=.01)

y0 = c(S=244, I1=3, I2=3, R=0)
times = seq(0,100,0.1)

out0 = ode(y0, times, diseaseModelVar, c(hivar, gs=0, as=0)) # no variation
outa1 = ode(y0, times, diseaseModelVar, c(hivar, g=0, as=0.4)) # variation in virulence
outg1 = ode(y0, times, diseaseModelVar, c(hivar, gs=0.4, as=0)) # variation in recovery
outa2 = ode(y0, times, diseaseModelVar, c(hivar, gs=0, as=0.8)) # variation in virulence
outg2 = ode(y0, times, diseaseModelVar, c(hivar, gs=0.8, as=0)) # variation in recovery

## Increasing variance increases peak epidemic size across all parameters, but more when there is 
## variation in virulence than when there is variation in recovery. 
c(no_var_max=max(mutate(as.data.frame(out0),I=I1+I2)$I),
  low_avar_max=max(mutate(as.data.frame(outa1),I=I1+I2)$I),
  low_gvar_max=max(mutate(as.data.frame(outg1),I=I1+I2)$I),
  hi_avar_max=max(mutate(as.data.frame(outa2),I=I1+I2)$I),
  hi_gvar_max=max(mutate(as.data.frame(outg2),I=I1+I2)$I))

## This suggests that more variation increases peak epidemics.
## But this is slightly different than what we're actually doing in the simulations
## The difference is that the traits are lognormally distributed, which means that 
## there are a lot of individuals with low trait values and a few with very 
## high trait values. So what if we consider that the more extreme we make the
## variation the larger the maximums, but the rarer they are. 

## Its hard to come up with a clean way to do this, so here's a clunkier version of it

## In this model 50% of individuals have a virulence of 0.08 and 50% have a virulence of 0.12
## so the average virulence is 0.1
diseaseModelLowAvar = function(t,y,parms) {
  b = parms["b"]
  d = parms["d"]
  bs = parms["bs"]
  c = parms["c"]
  shed = parms["shed"]
  gamma = parms["gamma"]
  alpha = parms["alpha"]

  S = y[1]
  I1 = y[2]
  I2 = y[3]
  R = y[4]
  
  dSdt = (b - bs*(S+I1+I2+R))*(S+I1+I2+R) - c*(shed/(1+shed))*S*(I1+I2) - d*S
  dI1dt = 1/2*c*(shed/(1+shed))*S*(I1+I2) - gamma*I1 - 0.12*I1 - d*I1
  dI2dt = 1/2*c*(shed/(1+shed))*S*(I1+I2) - gamma*I2 - 0.08*I2 - d*I2
  dRdt = gamma*I1 + gamma*I2 - d*R
  
  return(list(c(dSdt,dI1dt,dI2dt,dRdt),(b - bs*(S+I1+I2+R))*(S+I1+I2+R),c*(shed/(1+shed))*S*(I1+I2)))
  
}

## In this model 75% of individuals have a virulence of 0.025 and 25% have a virulence of 0.325
## so the average virulence is still 0.1
diseaseModelHiAvar = function(t,y,parms) {
  b = parms["b"]
  d = parms["d"]
  bs = parms["bs"]
  c = parms["c"]
  shed = parms["shed"]
  gamma = parms["gamma"]
  alpha = parms["alpha"]
  
  S = y[1]
  I1 = y[2]
  I2 = y[3]
  R = y[4]
  
  dSdt = (b - bs*(S+I1+I2+R))*(S+I1+I2+R) - c*(shed/(1+shed))*S*(I1+I2) - d*S
  dI1dt = 0.75*c*(shed/(1+shed))*S*(I1+I2) - gamma*I1 - 0.025*I1 - d*I1
  dI2dt = 0.25*c*(shed/(1+shed))*S*(I1+I2) - gamma*I2 - 0.325*I2 - d*I2
  dRdt = gamma*I1 + gamma*I2 - d*R
  
  return(list(c(dSdt,dI1dt,dI2dt,dRdt),(b - bs*(S+I1+I2+R))*(S+I1+I2+R),c*(shed/(1+shed))*S*(I1+I2)))
  
}


## In this model 50% of individuals have a recovery rate of 0.08 and 50% have a recovery of 0.12
## so the average virulence is 0.1
diseaseModelLowGvar = function(t,y,parms) {
  b = parms["b"]
  d = parms["d"]
  bs = parms["bs"]
  c = parms["c"]
  shed = parms["shed"]
  gamma = parms["gamma"]
  alpha = parms["alpha"]
  
  S = y[1]
  I1 = y[2]
  I2 = y[3]
  R = y[4]
  
  dSdt = (b - bs*(S+I1+I2+R))*(S+I1+I2+R) - c*(shed/(1+shed))*S*(I1+I2) - d*S
  dI1dt = 1/2*c*(shed/(1+shed))*S*(I1+I2) - 0.12*I1 - alpha*I1 - d*I1
  dI2dt = 1/2*c*(shed/(1+shed))*S*(I1+I2) - 0.08*I2 - alpha*I2 - d*I2
  dRdt = 0.12*I1 + 0.08*I2 - d*R
  
  return(list(c(dSdt,dI1dt,dI2dt,dRdt),(b - bs*(S+I1+I2+R))*(S+I1+I2+R),c*(shed/(1+shed))*S*(I1+I2)))
  
}

## In this model 75% of individuals have a recovery of 0.025 and 25% have a recovery of 0.325
## so the average virulence is still 0.1
diseaseModelHiGvar = function(t,y,parms) {
  b = parms["b"]
  d = parms["d"]
  bs = parms["bs"]
  c = parms["c"]
  shed = parms["shed"]
  gamma = parms["gamma"]
  alpha = parms["alpha"]
  
  S = y[1]
  I1 = y[2]
  I2 = y[3]
  R = y[4]
  
  dSdt = (b - bs*(S+I1+I2+R))*(S+I1+I2+R) - c*(shed/(1+shed))*S*(I1+I2) - d*S
  dI1dt = 0.75*c*(shed/(1+shed))*S*(I1+I2) - 0.025*I1 - alpha*I1 - d*I1
  dI2dt = 0.25*c*(shed/(1+shed))*S*(I1+I2) - 0.325*I2 - alpha*I2 - d*I2
  dRdt = 0.025*I1 + 0.325*I2 - d*R
  
  return(list(c(dSdt,dI1dt,dI2dt,dRdt),(b - bs*(S+I1+I2+R))*(S+I1+I2+R),c*(shed/(1+shed))*S*(I1+I2)))
  
}

x <- 0.24647
hivar = c(   c=(1-x)*0.1,    shed=(1-x)*1/9,    alpha=(1+x)*0.1,    gamma=(1+x)*0.1, 
             sd_c=(1-x)*0.5, sd_shed=(1-x)*5/9, sd_alpha=(1+x)*0.5, sd_gamma=(1+x)*0.5, 
             b=2.5, d=.1, bs=.01)

y0 = c(S=244, I1=3, I2=3, R=0)
times = seq(0,100,0.1)

## out0 is the no variation simulation to compare against
outLoAvar = ode(y0, times, diseaseModelLowAvar, hivar) # variation in virulence
outLoGvar = ode(y0, times, diseaseModelLowGvar, hivar) # variation in virulence
outHiAvar = ode(y0, times, diseaseModelHiAvar, hivar) # variation in virulence
outHiGvar = ode(y0, times, diseaseModelHiGvar, hivar) # variation in virulence

