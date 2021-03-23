### gillespie model with variation in beta, no evolution

pick_individuals <- function(N0, traitmean, traitsd) {
  meanlog <- log(traitmean^2/sqrt(traitsd^2+traitmean^2))
  sdlog <- sqrt(log(traitsd^2/traitmean^2+1))
  return(rlnorm(N0, meanlog=meanlog, sdlog=sdlog))
}


gillespie.SIR.varB <- function(tmax, params, x, seed=floor(runif(1,1,1e5))) {
  set.seed(seed)
  beta = params["beta"]
  alpha = params["alpha"]
  gamma = params["gamma"]
  varB = params["varB"]
  b = params["b"]
  bs = params["bs"]
  d = params["d"]
  ds = params["ds"]
  
  S=x["S"]
  I=x["I"]
  R=x["R"]
  
  ## draw the traits of our infected individuals
  beta_i <- pick_individuals(I, traitmean=beta, traitsd=sqrt(varB))
  #beta_i <- runif(I, 0.15, 0.35) # continuous variation w/ uniform distribution
  # vary shedding with c * s/(1+s) -> beta, c=.5, s= 1, .5, 2
  # varrying beta is equivalent to varrying c with a constant s
  
  # start at time 0
  t <- 0 
  
  results <- vector(mode='list', length=100000)
  results[[1]] <- list(t, S, beta_i, R)
  
  ## row counter
  i <- 2
  
  #start algorithm
  while(t < tmax & length(beta_i) > 0){
    irate = beta_i*S
    rrate = gamma*length(beta_i)
    brate <- (b - bs*(S+length(beta_i)+R)) * (S+length(beta_i)+R)
    drateS <-S*(d)
    drateI <-length(beta_i)*(d+alpha)
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
    if (event%in%1:length(beta_i)){### infection event
      S <- S-1
      beta_i <- c(beta_i, pick_individuals(1, traitmean=beta, traitsd=sqrt(varB)))
    }

    else if (event==length(beta_i)+1){ ### recover, remove a beta_i
      R <- R+1
      beta_i <- beta_i[-event]
    }
    else if (event==length(beta_i)+2){ ### birth
      S <- S+1
    }
    else if (event==length(beta_i)+3) { ### death of S
      S <- S-1
    }
    else if (event==length(beta_i)+4) {### death of I
      beta_i <- beta_i[-sample(1:length(beta_i),1)]
    }
    else    {   ### death of R
      R <- R-1 
    }
    
    results[[i]] <- list(t, S, beta_i, R)
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
params4 = c(beta=.25, alpha=.15, gamma=.15, varB=1e-3, b=2, d=0.4, bs=0.01, ds=0.01)

out15 <- gillespie.SIR.varB(tmax, params4, x)

plot.ts(out15[,2], col="blue", ylim=c(-5, 200))
lines(out15[,3], col="red")
lines(out15[,4], col="green")



##### CONTINUOUS VARIATION WITH UNIFORM DISTRIBUTION 
gillespie.SIR.varB.uniform <- function(tmax, params, x, seed=floor(runif(1,1,1e5))) {
  set.seed(seed)
  beta = params["beta"]
  alpha = params["alpha"]
  gamma = params["gamma"]
  varB = params["varB"]
  b = params["b"]
  bs = params["bs"]
  d = params["d"]
  ds = params["ds"]
  
  S=x["S"]
  I=x["I"]
  R=x["R"]
  
  ## draw the traits of our infected individuals
  beta_i <- runif(I, .02,.04) # continuous with uniform distribution

  # vary shedding with c * s/(1+s) -> beta, c=.5, s= 1, .5, 2
  # varrying beta is equivalent to varrying c with a constant s
  
  # start at time 0
  t <- 0 
  
  results <- vector(mode='list', length=100000)
  results[[1]] <- list(t, S, beta_i, R)
  
  ## row counter
  i <- 2
  
  #start algorithm
  while(t < tmax & length(beta_i) > 0){
    irate = beta_i*S
    rrate = gamma*length(beta_i)
    brate <- (b - bs*(S+length(beta_i)+R)) * (S+length(beta_i)+R)
    drateS <-S*(d)
    drateI <-length(beta_i)*(d+alpha)
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
    if (event%in%1:length(beta_i)){### infection event
      S <- S-1
      beta_i <- c(beta_i, runif(1, .02,.04))
    }
    
    else if (event==length(beta_i)+1){ ### recover, randomly remove a beta_i
      R <- R+1
      beta_i <- beta_i[-event]
    }
    else if (event==length(beta_i)+2){ ### birth
      S <- S+1
    }
    else if (event==length(beta_i)+3) { ### death of S
      S <- S-1
    }
    else if (event==length(beta_i)+4) {### death of I
      beta_i <- beta_i[-sample(1:length(beta_i),1)]
    }
    else    {   ### death of R
      R <- R-1 
    }
    
    results[[i]] <- list(t, S, beta_i, R)
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


##### STRATIFIED BETA VARIATION

gillespie.SIR.strat.varB <- function(tmax, params, x, seed=floor(runif(1,1,1e5))) {
  set.seed(seed)
  beta = params["beta"]
  alpha = params["alpha"]
  gamma = params["gamma"]
  varB = params["varB"]
  b = params["b"]
  bs = params["bs"]
  d = params["d"]
  ds = params["ds"]
  eps_b = params["epsilon_b"]
  
  S=x["S"]
  I=x["I"]
  R=x["R"]
  
  ## draw the traits of our infected individuals
  beta_values<-c(beta-eps_b,beta+eps_b)
  beta_i<-sample(beta_values,10, replace=TRUE)
  
  
  # start at time 0
  t <- 0 
  
  results <- vector(mode='list', length=100000)
  results[[1]] <- list(t, S, beta_i, R)
  
  ## row counter
  i <- 2
  
  #start algorithm
  while(t < tmax & length(beta_i) > 0){
    irate = beta_i*S
    rrate = gamma*length(beta_i)
    brate <- (b - bs*(S+length(beta_i)+R)) * (S+length(beta_i)+R)
    drateS <-S*(d)
    drateI <-length(beta_i)*(d+alpha)
    drateR <- R*(d)
    
    rates<-c(irate,rrate,brate,drateS,drateI,drateR)  
    print(rates)
  
    ## what time does the event happen?
    dt <- rexp(1, rate=sum(rates))
    
    ## update t 
    t <- t + dt
    
    ## "wheel of fortune"
    wheel <- cumsum(rates)/sum(rates)
    
    ## which event happens? Draw a random uniform to determine
    rand <- runif(1)
    
    ## if event==1, infection
    ## if event==2, recovery
    ## if event==3, birth
    ## if event==4, death of S
    ## if event==5, death of I
    ## if event==6, death of R
    event <- 1 + sum(rand > wheel) 
    if (event%in%1:length(beta_i)){### infection event
      S <- S-1
      beta_i <- c(beta_i, sample(beta_values, 1, replace=TRUE))
    }
    
    else if (event==length(beta_i)+1){ ### recover,remove a beta_i
      R <- R+1
      beta_i <- beta_i[-event]
    }
    else if (event==length(beta_i)+2){ ### birth
      S <- S+1
    }
    else if (event==length(beta_i)+3) { ### death of S
      S <- S-1
    }
    else if (event==length(beta_i)+4) {### death of I
      beta_i <- beta_i[-sample(1:length(beta_i),1)]
    }
    else    {   ### death of R
      R <- R-1 
    }
    
    results[[i]] <- list(t, S, beta_i, R)
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








