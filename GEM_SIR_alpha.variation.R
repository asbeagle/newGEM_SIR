### gillespie model with variation in alpha, no evolution

pick_individuals <- function(N0, traitmean, traitsd) {
  meanlog <- log(traitmean^2/sqrt(traitsd^2+traitmean^2))
  sdlog <- sqrt(log(traitsd^2/traitmean^2+1))
  return(rlnorm(N0, meanlog=meanlog, sdlog=sdlog))
}



gillespie.SIR.varA <- function(tmax, params, x, seed=floor(runif(1,1,1e5))) {
  set.seed(seed)
  beta = params["beta"]
  alpha = params["alpha"]
  gamma = params["gamma"]
  varA = params["varA"]
  b = params["b"]
  bs = params["bs"]
  d = params["d"]
  ds = params["ds"]
  
  S=x["S"]
  I=x["I"]
  R=x["R"]
  
  ## draw the traits of our infected individuals
  alpha_i <- pick_individuals(I, traitmean=alpha, traitsd=sqrt(varA))
  
  # start at time 0
  t <- 0 
  
  results <- vector(mode='list', length=100000)
  results[[1]] <- list(t, S, alpha_i, R)
  
  ## row counter
  i <- 2
  
  #start algorithm
  while(t < tmax & length(alpha_i) > 0){
    irate = length(alpha_i)*S*beta
    rrate = gamma*length(alpha_i)
    brate = (b - bs*(S+length(alpha_i)+R)) * (S+length(alpha_i)+R)
    drateS = S*(d)
    drateI = (alpha_i + d)
    drateR = R*(d)
    
    rates<-c(drateI,irate,rrate,brate,drateS,drateR)  
    ## what time does the event happen?
    dt <- rexp(1, rate=sum(rates))
    
    ## update t 
    t <- t + dt
    
    ## "wheel of fortune"
    wheel <- cumsum(rates)/sum(rates)
    
    ## which event happens? Draw a random uniform to determine
    rand <- runif(1)
    
    event <- 1 + sum(rand > wheel) 
    if (event%in%1:length(alpha_i)){### death of an I
      
      alpha_i <- alpha_i[-sample(1:length(alpha_i),1)]
      
    }
    else if(event==length(alpha_i)+1){ # infection
      S <- S-1
      alpha_i <- c(alpha_i, pick_individuals(1, traitmean=alpha, traitsd=sqrt(varA)))
    }
    else if(event==length(alpha_i)+2){
      alpha_i <- alpha_i[-sample(1:length(alpha_i),1)] # recovery
      R <-R+1
    }
    
    else if (event==length(alpha_i)+3){ ### birth
      S <- S+1
    }
    else if (event==length(alpha_i)+4) { ### death of S
      S <- S-1
    }
    else    {   ### death of R
      R <- R-1
    }
    
    
    results[[i]] <- list(t, S, alpha_i, R)
    i <- i +1
    
  }  
  results <- results[1:(i-1)]
  ## simplify what you're storing so R doesn't crash because the datafiles are too big
  results <- data.frame(time=sapply(results, function(r) r[[1]]),
                        S=sapply(results, function(r) r[[2]]), 
                        I=sapply(results, function(r) length(r[[3]])),
                        R=sapply(results, function(r) r[[4]]))
  ## only keep the state of the system at times 0, 0.1, 0.2, ..., tmax-0.2, tmax-0.1, tmax
  results <- results[sapply(seq(0,tmax,1), function(t) min(which(results$time >= t))),]
  ## these two steps reduce the size of 'results' ~1000-fold
  
  return(results)
}


## output
x = c(S=70, I=10, R=0)
tmax <- 150
params4 = c(beta=.25, alpha=.01, gamma=.15, varG=1e-3, b=2, d=0.4, bs=0.01, ds=0.01, varA=1e-3)

out14 <- gillespie.SIR.var.A(tmax, params4, x)

plot.ts(out14[,2], col="blue", ylim=c(-5, 200))
lines(out14[,3], col="red")
lines(out14[,4], col="green")


###### continuous variation with random uniform
gillespie.SIR.varA.uniform <- function(tmax, params, x, seed=floor(runif(1,1,1e5))) {
  set.seed(seed)
  beta = params["beta"]
  alpha = params["alpha"]
  gamma = params["gamma"]
  varA = params["varA"]
  b = params["b"]
  bs = params["bs"]
  d = params["d"]
  ds = params["ds"]
  
  S=x["S"]
  I=x["I"]
  R=x["R"]
  
  ## draw the traits of our infected individuals
  alpha_i <- runif(I, .15,.35)
  
  # start at time 0
  t <- 0 
  
  results <- vector(mode='list', length=100000)
  results[[1]] <- list(t, S, alpha_i, R)
  
  ## row counter
  i <- 2
  
  #start algorithm
  while(t < tmax & length(alpha_i) > 0){
    irate = length(alpha_i)*S*beta
    rrate = gamma*length(alpha_i)
    brate = (b - bs*(S+length(alpha_i)+R)) * (S+length(alpha_i)+R)
    drateS = S*(d)
    drateI = (alpha_i + d)
    drateR = R*(d)
    
    rates<-c(drateI,irate,rrate,brate,drateS,drateR)  
    
    ## what time does the event happen?
    dt <- rexp(1, rate=sum(rates))
    
    ## update t 
    t <- t + dt
    
    ## "wheel of fortune"
    wheel <- cumsum(rates)/sum(rates)
    
    ## which event happens? Draw a random uniform to determine
    rand <- runif(1)
    
    event <- 1 + sum(rand > wheel) 
    if (event%in%1:length(alpha_i)){### death of an I
      
      alpha_i <- alpha_i[-sample(1:length(alpha_i),1)]
      
    }
    else if(event==length(alpha_i)+1){ # infection
      S <- S-1
      alpha_i <- c(alpha_i, runif(1, .15,.35))
    }
    else if(event==length(alpha_i)+2){
      alpha_i <- alpha_i[-sample(1:length(alpha_i),1)] # recovery
      R <-R+1
    }
    
    else if (event==length(alpha_i)+3){ ### birth
      S <- S+1
    }
    else if (event==length(alpha_i)+4) { ### death of S
      S <- S-1
    }
    else    {   ### death of R
      R <- R-1
    }
    
    
    results[[i]] <- list(t, S, alpha_i, R)
    i <- i +1
    
  }  
  results <- results[1:(i-1)]
  ## simplify what you're storing so R doesn't crash because the datafiles are too big
  results <- data.frame(time=sapply(results, function(r) r[[1]]),
                        S=sapply(results, function(r) r[[2]]), 
                        I=sapply(results, function(r) length(r[[3]])),
                        R=sapply(results, function(r) r[[4]]))
  ## only keep the state of the system at times 0, 0.1, 0.2, ..., tmax-0.2, tmax-0.1, tmax
  results <- results[sapply(seq(0,tmax,1), function(t) min(which(results$time >= t))),]
  ## these two steps reduce the size of 'results' ~1000-fold
  
  return(results)
}

## This is a modification of the code in GEM_SIR_strat.varA.R that allows only two trait values
## at extremes of alpha+eps and alpha-eps
gillespie.SIR.strat.varA <- function(tmax, params, x, seed=floor(runif(1,1,1e5))) {
  set.seed(seed)
  b = params["b"]
  bs = params["bs"]
  d = params["d"]
  beta = params["beta"]
  gamma = params["gamma"]
  alpha = params["alpha"]
  eps = params["epsilon"]
  
  S=x["S"]
  I=x["I"]
  R=x["R"]
  
  ## draw the traits of our infected individuals
  alpha_values<-c(alpha-eps,alpha+eps)
  alpha_i<-sample(alpha_values, I, replace=TRUE)
  
  # start at time 0
  t <- 0 
  
  results <- vector(mode='list', length=100000)
  results[[1]] <- list(t, S, alpha_i, R) 
  
  ## row counter
  i <- 2
  
  #start algorithm
  while(t < tmax & length(alpha_i)>0 ){ 
    irate = beta*S*(length(alpha_i))
    rrate = gamma*(length(alpha_i))
    brate <- (b - bs*(S+R+(length(alpha_i)))) * (S+R+(length(alpha_i)))
    drateS <- S *(d)
    drateI <-(alpha_i + d)
    drateR <- R*(d)
    
    rates<-c(drateI,irate,rrate,brate,drateS,drateR)
    #print(rates)
    ## what time does the event happen?
    dt <- rexp(1, rate=sum(rates))
    
    ## update t 
    t <- t + dt
    
    ## "wheel of fortune"
    wheel <- cumsum(rates)/sum(rates)
    
    ## which event happens? Draw a random uniform to determine
    rand <- runif(1)
    
    #(drateI,irate,rrate,brate,drateS,drateR)
    event <- 1 + sum(rand > wheel) 
    if (event%in%1:length(alpha_i)){### death of an I
      
      alpha_i <- alpha_i[-sample(1:length(alpha_i),1)]
      
    }
    else if(event==length(alpha_i)+1){ # infection
      S <- S-1
      # If you want to let things evolve, use these lines
      #this_gamma_r <- gamma_r[sample(1:length(gamma_r), 1)]
      #new_gamma_r <- pick_individuals(1, traitmean=this_gamma_r,traitsd=sqrt(varG))
      # If you don't want evolution, use this line
      alpha_i <- c(alpha_i, sample(alpha_values, 1, replace=TRUE))
    }
    else if(event==length(alpha_i)+2){
      alpha_i <- alpha_i[-sample(1:length(alpha_i),1)] # recovery
      R <-R+1
    }
    
    else if (event==length(alpha_i)+3){ ### birth
      S <- S+1
    }
    else if (event==length(alpha_i)+4) { ### death of S
      S <- S-1
    }
    else    {   ### death of R
      R <- R-1
    }
    
    results[[i]] <- list(t, S, alpha_i , R)
    i <- i +1
    
  }  
  results <- results[1:(i-1)]
  ## simplify what you're storing so R doesn't crash because the datafiles are too big
  results <- data.frame(time=sapply(results, function(r) r[[1]]),
                        S=sapply(results, function(r) r[[2]]), 
                        I=sapply(results, function(r) length(r[[3]])),
                        R=sapply(results, function(r) r[[4]]))
  ## only keep the state of the system at times 0, 0.1, 0.2, ..., tmax-0.2, tmax-0.1, tmax
  results <- results[sapply(seq(0,tmax,1), function(t) min(which(results$time >= t))),]
  ## these two steps reduce the size of 'results' ~1000-fold
  
  return(results)
}


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
plot(seq(0,150,1), S_dyn_novar, type='l', lwd=2, ylim=c(30, 250))
lines(seq(0,150,1), S_dyn_contvar, lwd=2, col=2)
lines(seq(0,150,1), S_dyn_discvar, lwd=2, col=4)
lines(det_out_novar[,c(1,2)], lwd=2, lty=2, col=gray(0.5))
lines(det_out_discvar[,c(1,2)], lwd=2, lty=2, col=2)

