### gillespie model with variation in alpha, no evolution

pick_individuals <- function(N0, traitmean, traitsd) {
  meanlog <- log(traitmean^2/sqrt(traitsd^2+traitmean^2))
  sdlog <- sqrt(log(traitsd^2/traitmean^2+1))
  return(rlnorm(N0, meanlog=meanlog, sdlog=sdlog))
}


#### alpha continuous variation
gillespie.SIR.varA <- function(tmax, params, x, seed=floor(runif(1,1,1e5))) {
  set.seed(seed)
  beta = params["beta"]
  c = params["c"]
  shed = params["shed"]
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
    irate = length(alpha_i)*S*(c*(shed/(1+shed)))
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
      ## This is the mistake - don't remove an individual at random, remove the *particular* individual
      ## alpha_i <- alpha_i[-sample(1:length(alpha_i),1)]
      alpha_i <- alpha_i[-event]
      
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

#### alpha continuous random uniform variation
gillespie.SIR.varA.uniform <- function(tmax, params, x, seed=floor(runif(1,1,1e5))) {
  set.seed(seed)
  beta = params["beta"]
  c = params["c"]
  shed = params["shed"]
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
    irate = length(alpha_i)*S*(c*(shed/(1+shed)))
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
      
      alpha_i <- alpha_i[-event]
      
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
  c = params["c"]
  shed = params["shed"]
  gamma = params["gamma"]
  alpha = params["alpha"]
  eps = params["epsilon"]
  
  S=x["S"]
  I=x["I"]
  R=x["R"]
  
  ## draw the traits of our infected individuals
  alpha_values<-c((alpha-eps),(alpha+eps))
  alpha_i<-sample(alpha_values, I, replace=TRUE)
  
  # start at time 0
  t <- 0 
  
  results <- vector(mode='list', length=100000)
  results[[1]] <- list(t, S, alpha_i, R) 
  
  ## row counter
  i <- 2
  
  #start algorithm
  while(t < tmax & length(alpha_i)>0 ){ 
    irate = S*(length(alpha_i))*(c*(shed/(1+shed)))
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
      ## here is maybe the mistake - don't remove an individual at *random* - remove it according to its traits
      ##alpha_i <- alpha_i[-sample(1:length(alpha_i),1)]
      alpha_i <- alpha_i[-event]
    }
    else if(event==length(alpha_i)+1){ # infection
      S <- S-1
      # If you want to let things evolve, use these lines
      #this_gamma_r <- gamma_r[sample(1:length(gamma_r), 1)]
      #new_gamma_r <- pick_individuals(1, traitmean=this_gamma_r,traitsd=sqrt(varG))
      # If you don't want evolution, use this line
      alpha_i <- c(alpha_i, sample(alpha_values, 1))
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


