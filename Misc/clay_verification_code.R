#### clays code used to verify that changing alpha removal fixes the problem of not seeing any difference
#### between no variation and variation models

## Parameters that should produce a high level of difference between the no-variation and continuous variation cases
params_novar <- c(beta=0.0025, alpha=0.15, gamma=0.001, varG=1e-3, b=2.5, d=0.001, bs=0.01, ds=0.01, varA=0)
params_contvar <- c(beta=0.0025, alpha=0.15, gamma=0.001, varG=1e-3, b=2.5, d=0.001, bs=0.01, ds=0.01, varA=0.15)
params_discvar <- c(beta=0.0025, alpha=0.15, gamma=0.001, varG=1e-3, b=2.5, d=0.001, bs=0.01, ds=0.01, varA=0.15, epsilon=0.1)

## Choose initial values that are close to the disease-free equilibirum
initial_state <- floor(c(S=unname((params_novar["b"]-params_novar["d"])/params_novar["bs"])-5, I=5, R=0)) 

seeds <- floor(runif(100, 1, 1e5))
mclapply(seeds, 
         function(s) gillespie.SIR.varA(tmax=150, params=params_novar, x=initial_state, seed=s),
         mc.cores=4) -> out_novar

mclapply(seeds, 
         function(s) gillespie.SIR.varA(tmax=150, params=params_contvar, x=initial_state, seed=s),
         mc.cores=4) -> out_contvar

mclapply(seeds, 
         function(s) gillespie.SIR.strat.varA(tmax=150, params=params_discvar, x=initial_state, seed=s),
         mc.cores=4) -> out_discvar


## Plot the dynamics of the susceptible population
lapply(out_novar, function(l) l[,2]) %>%
  do.call("cbind.data.frame",.) %>% 
  apply(., 1, function(x) mean(x,na.rm=TRUE)) -> S_dyn_novar

lapply(out_contvar, function(l) l[,2]) %>%
  do.call("cbind.data.frame",.) %>% 
  apply(., 1, function(x) mean(x,na.rm=TRUE)) -> S_dyn_contvar

lapply(out_discvar, function(l) l[,2]) %>%
  do.call("cbind.data.frame",.) %>% 
  apply(., 1, function(x) mean(x,na.rm=TRUE)) -> S_dyn_discvar

## Simulate the deterministic model as well
library(deSolve)
deterministic.SIR.novar=function(t,x,params){
  beta = params["beta"]
  alpha = params["alpha"]
  gamma = params["gamma"]
  b = params["b"]
  bs = params["bs"]
  d = params["d"]
  eps = params["epsilon"]
  S=x[1]
  I=x[2]
  R=x[3]
  
  dS  = -S*I*beta + ((b-bs*(S+I+R))*(S+I+R)) - (d*(S))
  dI  = S*I*beta - I*(alpha+gamma+d)
  dR  = gamma*(I) - d*R
  
  list(c(dS,dI,dR))
}

deterministic.SIR.discvar=function(t,x,params){
  beta = params["beta"]
  alpha = params["alpha"]
  gamma = params["gamma"]
  b = params["b"]
  bs = params["bs"]
  d = params["d"]
  eps = params["epsilon"]
  S=x[1]
  I1=x[2]
  I2=x[3]
  R=x[4]
  
  dS  = -S*(I1+I2)*beta + ((b-bs*(S+I1+I2+R))*(S+I1+I2+R)) - (d*(S))
  dI1 = S*(I1+I2)*beta/2 - I1*(alpha-eps+gamma+d)
  dI2 = S*(I1+I2)*beta/2 - I2*(alpha+eps+gamma+d)
  dR  = gamma*(I1+I2) - d*R
  
  list(c(dS,dI1,dI2,dR))
}

params_discvar <- c(beta=0.0025, alpha=0.15, gamma=0.001, varG=1e-3, b=2.5, d=0.001, bs=0.01, ds=0.01, varA=0.15, epsilon=0.1)
initial_state_novar <- floor(c(S=unname((params_novar["b"]-params_novar["d"])/params_novar["bs"])-4, I=4, R=0)) 
initial_state_discvar <- floor(c(S=unname((params_novar["b"]-params_novar["d"])/params_novar["bs"])-4, I1=2, I2=2, R=0)) 

det_out_novar <- ode(y=initial_state_novar, times=seq(0,150,1), func=deterministic.SIR.novar, parms=params_discvar)
det_out_discvar <- ode(y=initial_state_discvar, times=seq(0,150,1), func=deterministic.SIR.discvar, parms=params_discvar)

## PLOT EVERYTHING
## Something is not right because the stochastic simulation with discrete variation does not come even remotely close to the deterministic expectation
plot(seq(0,150,1), S_dyn_novar, type='l', lwd=2, ylim=c(10, 250), xlab="Time", ylab="Number susceptible")
lines(seq(0,150,1), S_dyn_contvar, lwd=2, col=3)
lines(seq(0,150,1), S_dyn_discvar, lwd=2, col=4)
abline(h=det_out_novar[151,2], lwd=1, lty=2)
abline(h=det_out_discvar[151,2], lwd=1, lty=2, col=4)
legend(x='topright', c("No variation", "Continuous variation", "Discrete variation"), fill=c(1,3,4), bty='n')
}