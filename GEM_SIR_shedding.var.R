### GEM SIR WITH VARIATION IN SHEDDING

pick_individuals <- function(N0, traitmean, traitsd) {
  meanlog <- log(traitmean^2/sqrt(traitsd^2+traitmean^2))
  sdlog <- sqrt(log(traitsd^2/traitmean^2+1))
  return(rlnorm(N0, meanlog=meanlog, sdlog=sdlog))
}


gillespie.SIR.varS <- function(tmax, params, x, seed=floor(runif(1,1,1e5))) {
  set.seed(seed)
  beta = params["beta"]
  alpha = params["alpha"]
  gamma = params["gamma"]
  shed = params["shed"]
  c = params["c"]
  varS = params["varS"]
  b = params["b"]
  bs = params["bs"]
  d = params["d"]
  ds = params["ds"]
  
  S=x["S"]
  I=x["I"]
  R=x["R"]
  
  ## draw the traits of our infected individuals
  new_i<-pick_individuals(I, traitmean=shed, traitsd=sqrt(varS))
  shed_i <- c * (new_i)/(1+new_i)
  #beta_i <- runif(I, 0.15, 0.35) # continuous variation w/ uniform distribution
  # vary shedding with c * s/(1+s) -> beta, c=.5, s= 1, .5, 2
  # varrying beta is equivalent to varrying c with a constant s
  
  # start at time 0
  t <- 0 
  
  results <- vector(mode='list', length=100000)
  results[[1]] <- list(t, S, shed_i, R)
  
  ## row counter
  i <- 2
  
  #start algorithm
  while(t < tmax & length(shed_i) > 0){
    irate = shed_i*S
    rrate = gamma*length(shed_i)
    brate <- (b - bs*(S+length(shed_i)+R)) * (S+length(shed_i)+R)
    drateS <-S*(d)
    drateI <-length(shed_i)*(d+alpha)
    drateR <- R*(d)
    
    rates<-c(irate,rrate,brate,drateS,drateI,drateR)  
    
    ## what time does the event happen?
    dt <- rexp(1, rate=sum(rates))
    
    ## update t 
    t <- t + dt
    
    ## "wheel of fortune"
    wheel <- cumsum(rates)/sum(rates)
    
    ## which event happens? Draw a random uniform to determine
    rand <- runif(1)
    
    event <- 1 + sum(rand > wheel) 
    if (event%in%1:length(shed_i)){### infection event
      S <- S-1
      new_i<-pick_individuals(1, traitmean=shed, traitsd=sqrt(varS))
      new_shed_i <- c * (new_i)/(1+new_i)
      shed_i <- c(shed_i, new_shed_i)
    }
    
    else if (event==length(shed_i)+1){ ### recover, randomly remove a beta_i
      R <- R+1
      shed_i <- shed_i[-sample(1:length(shed_i),1)]
    }
    else if (event==length(shed_i)+2){ ### birth
      S <- S+1
    }
    else if (event==length(shed_i)+3) { ### death of S
      S <- S-1
    }
    else if (event==length(shed_i)+4) {### death of I
      shed_i <- shed_i[-sample(1:length(shed_i),1)]
    }
    else    {   ### death of R
      R <- R-1 
    }
    
    results[[i]] <- list(t, S, shed_i, R)
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

x = c(S=70, I=10, R=0)
tmax <- 150
params4 = c(beta=.25, alpha=.15, gamma=.15, varB=1e-3, b=2, d=0.4, bs=0.01, ds=0.01, c=.5, shed=.15, 
            varS=1e-3)

out15 <- gillespie.SIR.varS(tmax, params4, x)

plot.ts(out15[,2], col="blue", ylim=c(-5, 200))
lines(out15[,3], col="red")
lines(out15[,4], col="green")

##### VARIATION IN SHEDDING WITH CONT UNIFORM DIST
gillespie.SIR.varS.uniform <- function(tmax, params, x, seed=floor(runif(1,1,1e5))) {
  set.seed(seed)
  beta = params["beta"]
  alpha = params["alpha"]
  gamma = params["gamma"]
  shed = params["shed"]
  c = params["c"]
  varS = params["varS"]
  b = params["b"]
  bs = params["bs"]
  d = params["d"]
  ds = params["ds"]
  
  S=x["S"]
  I=x["I"]
  R=x["R"]
  
  ## draw the traits of our infected individuals
  new_i<-runif(I, .2, 2)
  shed_i <- c * (new_i)/(1+new_i)
  #beta_i <- runif(I, 0.15, 0.35) # continuous variation w/ uniform distribution
  # vary shedding with c * s/(1+s) -> beta, c=.5, s= 1, .5, 2
  # varrying beta is equivalent to varrying c with a constant s
  
  # start at time 0
  t <- 0 
  
  results <- vector(mode='list', length=100000)
  results[[1]] <- list(t, S, shed_i, R)
  
  ## row counter
  i <- 2
  
  #start algorithm
  while(t < tmax & length(shed_i) > 0){
    irate = shed_i*S
    rrate = gamma*length(shed_i)
    brate <- (b - bs*(S+length(shed_i)+R)) * (S+length(shed_i)+R)
    drateS <-S*(d)
    drateI <-length(shed_i)*(d+alpha)
    drateR <- R*(d)
    
    rates<-c(irate,rrate,brate,drateS,drateI,drateR)  
    
    ## what time does the event happen?
    dt <- rexp(1, rate=sum(rates))
    
    ## update t 
    t <- t + dt
    
    ## "wheel of fortune"
    wheel <- cumsum(rates)/sum(rates)
    
    ## which event happens? Draw a random uniform to determine
    rand <- runif(1)
    
    event <- 1 + sum(rand > wheel) 
    if (event%in%1:length(shed_i)){### infection event
      S <- S-1
      new_i<-runif(1, 1, 2)
      new_shed_i <- c * (new_i)/(1+new_i)
      shed_i <- c(shed_i, new_shed_i)
    }
    
    else if (event==length(shed_i)+1){ ### recover, randomly remove a beta_i
      R <- R+1
      shed_i <- shed_i[-sample(1:length(shed_i),1)]
    }
    else if (event==length(shed_i)+2){ ### birth
      S <- S+1
    }
    else if (event==length(shed_i)+3) { ### death of S
      S <- S-1
    }
    else if (event==length(shed_i)+4) {### death of I
      shed_i <- shed_i[-sample(1:length(shed_i),1)]
    }
    else    {   ### death of R
      R <- R-1 
    }
    
    results[[i]] <- list(t, S, shed_i, R)
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

x = c(S=70, I=10, R=0)
tmax <- 150
params4 = c(beta=.25, alpha=.15, gamma=.15, varB=1e-3, b=2, d=0.4, bs=0.01, ds=0.01, c=.5, shed=.15, 
            varS=1e-3)

out15 <- gillespie.SIR.varS.uniform(tmax, params4, x)

plot.ts(out15[,2], col="blue", ylim=c(-5, 200))
lines(out15[,3], col="red")
lines(out15[,4], col="green")

##### VARIATION IN SHEDDING WITH STRAT VARIATION
gillespie.SIR.strat.varS <- function(tmax, params, x, seed=floor(runif(1,1,1e5))) {
  set.seed(seed)
  beta = params["beta"]
  alpha = params["alpha"]
  gamma = params["gamma"]
  shed = params["shed"]
  c = params["c"]
  varS = params["varS"]
  b = params["b"]
  bs = params["bs"]
  d = params["d"]
  ds = params["ds"]
  
  S=x["S"]
  I=x["I"]
  R=x["R"]
  
  ## draw the traits of our infected individuals
  shed_value<-c(.2,1,2)
  new_i<-sample(shed_value, I, replace=TRUE)
  shed_i <- c * (new_i)/(1+new_i)
  #beta_i <- runif(I, 0.15, 0.35) # continuous variation w/ uniform distribution
  # vary shedding with c * s/(1+s) -> beta, c=.5, s= 1, .5, 2
  # varrying beta is equivalent to varrying c with a constant s
  
  # start at time 0
  t <- 0 
  
  results <- vector(mode='list', length=100000)
  results[[1]] <- list(t, S, shed_i, R)
  
  ## row counter
  i <- 2
  
  #start algorithm
  while(t < tmax & length(shed_i) > 0){
    irate = shed_i*S
    rrate = gamma*length(shed_i)
    brate <- (b - bs*(S+length(shed_i)+R)) * (S+length(shed_i)+R)
    drateS <-S*(d)
    drateI <-length(shed_i)*(d+alpha)
    drateR <- R*(d)
    
    rates<-c(irate,rrate,brate,drateS,drateI,drateR)  
    
    ## what time does the event happen?
    dt <- rexp(1, rate=sum(rates))
    
    ## update t 
    t <- t + dt
    
    ## "wheel of fortune"
    wheel <- cumsum(rates)/sum(rates)
    
    ## which event happens? Draw a random uniform to determine
    rand <- runif(1)
    
    event <- 1 + sum(rand > wheel) 
    if (event%in%1:length(shed_i)){### infection event
      S <- S-1
      new_i<-sample(shed_value, 1, 2)
      new_shed_i<- c * (new_i)/(1+new_i)
      shed_i <- c(shed_i, new_shed_i)
    }
    
    else if (event==length(shed_i)+1){ ### recover, randomly remove a beta_i
      R <- R+1
      shed_i <- shed_i[-sample(1:length(shed_i),1)]
    }
    else if (event==length(shed_i)+2){ ### birth
      S <- S+1
    }
    else if (event==length(shed_i)+3) { ### death of S
      S <- S-1
    }
    else if (event==length(shed_i)+4) {### death of I
      shed_i <- shed_i[-sample(1:length(shed_i),1)]
    }
    else    {   ### death of R
      R <- R-1 
    }
    
    results[[i]] <- list(t, S, shed_i, R)
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

out15 <- gillespie.SIR.strat.varS(tmax, params4, x)

plot.ts(out15[,2], col="blue", ylim=c(-5, 200))
lines(out15[,3], col="red")
lines(out15[,4], col="green")
