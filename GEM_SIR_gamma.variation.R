### gillespie model with variation in gamma, no evolution

pick_individuals <- function(N0, traitmean, traitsd) {
  meanlog <- log(traitmean^2/sqrt(traitsd^2+traitmean^2))
  sdlog <- sqrt(log(traitsd^2/traitmean^2+1))
  return(rlnorm(N0, meanlog=meanlog, sdlog=sdlog))
}



gillespie.SIR.varG <- function(tmax, params, x, seed=floor(runif(1,1,1e5))) {
  set.seed(seed)
  beta = params["beta"]
  alpha = params["alpha"]
  gamma = params["gamma"]
  varG = params["varG"]
  b = params["b"]
  bs = params["bs"]
  d = params["d"]
  ds = params["ds"]
  
  S=x["S"]
  I=x["I"]
  R=x["R"]
  
  ## draw the traits of our infected individuals
  gamma_i <- pick_individuals(I, traitmean=gamma, traitsd=sqrt(varG))
  
  # start at time 0
  t <- 0 
  
  results <- vector(mode='list', length=100000)
  results[[1]] <- list(t, S, gamma_i, R) 
  
  ## row counter
  i <- 2
  
  #start algorithm
  while(t < tmax & length(gamma_i)>0 ){ 
    irate = beta*S*length(gamma_i)
    rrate = gamma_i
    brate <- (b - bs*(S+R+length(gamma_i))) * (S+R+length(gamma_i))
    drateS <-S*(d)
    drateI <-length(gamma_i)*(d+alpha)
    drateR <- R*(d)
    
    rates<-c(rrate,irate,brate,drateS,drateI,drateR)
    #print(rates)
    ## what time does the event happen?
    dt <- rexp(1, rate=sum(rates))
    
    ## update t 
    t <- t + dt
    
    ## "wheel of fortune"
    wheel <- cumsum(rates)/sum(rates)
    
    ## which event happens? Draw a random uniform to determine
    rand <- runif(1)
    
    ## if event==1, recovery
    ## if event==2, infection
    ## if event==3, birth
    ## if event==4, death of S
    ## if event==5, death of I
    ## if event==6, death of R
    
    #(rrate,irate,brate,drateS,drateI,drateR)
    event <- 1 + sum(rand > wheel) 
    if (event%in%1:length(gamma_i)){### recovery
      gamma_i <- gamma_i[-event]
      R <-R+1
    }
    
     else if(event==length(gamma_i)+1){ # infection
      S<- S-1
      # If you want to let things evolve, use these lines
      #this_gamma_r <- gamma_r[sample(1:length(gamma_r), 1)]
      #new_gamma_r <- pick_individuals(1, traitmean=this_gamma_r,traitsd=sqrt(varG))
      # If you don't want evolution, use this line
      new_gamma_i <- pick_individuals(1, traitmean=gamma, traitsd=sqrt(varG))
      gamma_i <- c(gamma_i, new_gamma_i)
    }
    else if (event==length(gamma_i)+2){ ### birth
      S <- S+1
    }
    else if (event==length(gamma_i)+3) { ### death of S
      S <- S-1
    }
    else if (event==length(gamma_i)+4) {### death of I (not caused by infection)
      gamma_i <- gamma_i[-sample(1:length(gamma_i),1)]
    }
    else    {   ### death of R
      R <- R-1
    }
    
    results[[i]] <- list(t, S,gamma_i , R)
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
params4 = c(beta=.25, alpha=.15, gamma=.15, varG=1e-3, b=2, d=0.4, bs=0.01, ds=0.01)

out14 <- gillespie.SIR.varG(tmax, params4, x)

plot.ts(out14[,2], col="blue", ylim=c(-5, 200))
lines(out14[,3], col="red")
lines(out14[,4], col="green")

#### STRATIFIED VARIATION IN GAMMA
gillespie.SIR.strat.varG <- function(tmax, params, x, seed=floor(runif(1,1,1e5))) {
  set.seed(seed)
  beta = params["beta"]
  alpha = params["alpha"]
  gamma = params["gamma"]
  varG = params["varG"]
  b = params["b"]
  bs = params["bs"]
  d = params["d"]
  ds = params["ds"]
  eps = params["epsilon"]
  
  S=x["S"]
  I=x["I"]
  R=x["R"]
  
  ## draw the traits of our infected individuals
  gamma_values<-c(gamma-eps,gamma+eps)
  gamma_i<-sample(gamma_values, I, replace=TRUE)
  
  # start at time 0
  t <- 0 
  
  results <- vector(mode='list', length=100000)
  results[[1]] <- list(t, S, gamma_i, R) 
  
  ## row counter
  i <- 2
  
  #start algorithm
  while(t < tmax & length(gamma_i)>0 ){ 
    irate = beta*S*length(gamma_i)
    rrate = gamma_i
    brate <- (b - bs*(S+R+length(gamma_i))) * (S+R+length(gamma_i))
    drateS <-S*(d)
    drateI <-length(gamma_i)*(d+alpha)
    drateR <- R*(d)
    
    rates<-c(rrate,irate,brate,drateS,drateI,drateR)
    #print(rates)
    ## what time does the event happen?
    dt <- rexp(1, rate=sum(rates))
    
    ## update t 
    t <- t + dt
    
    ## "wheel of fortune"
    wheel <- cumsum(rates)/sum(rates)
    
    ## which event happens? Draw a random uniform to determine
    rand <- runif(1)
    
    ## if event==1, recovery
    ## if event==2, infection
    ## if event==3, birth
    ## if event==4, death of S
    ## if event==5, death of I
    ## if event==6, death of R
    
    #(rrate,irate,brate,drateS,drateI,drateR)
    event <- 1 + sum(rand > wheel) 
    if (event%in%1:length(gamma_i)){### recovery
      gamma_i <- gamma_i[-event]
      R <-R+1
    }
    
    else if(event==length(gamma_i)+1){ # infection
      S<- S-1
      # If you want to let things evolve, use these lines
      #this_gamma_r <- gamma_r[sample(1:length(gamma_r), 1)]
      #new_gamma_r <- pick_individuals(1, traitmean=this_gamma_r,traitsd=sqrt(varG))
      # If you don't want evolution, use this line
      gamma_i <- c(gamma_i, sample(gamma_values, 1, replace=TRUE))
    }
    else if (event==length(gamma_i)+2){ ### birth
      S <- S+1
    }
    else if (event==length(gamma_i)+3) { ### death of S
      S <- S-1
    }
    else if (event==length(gamma_i)+4) {### death of I 
      gamma_i <- gamma_i[-sample(1:length(gamma_i),1)]
    }
    else    {   ### death of R
      R <- R-1
    }
    
    results[[i]] <- list(t, S,gamma_i , R)
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
params4 = c(beta=.25, alpha=.01, gamma=.3, varG=1e-3, b=2, d=0.4, bs=0.01, ds=0.01)

out14 <- gillespie.SIR.strat.varG(tmax, params4, x)

plot.ts(out14[,2], col="blue", ylim=c(-5, 200))
lines(out14[,3], col="red")
lines(out14[,4], col="green")


#### continous variation with random uniform distribution
gillespie.SIR.varG.uniform <- function(tmax, params, x, seed=floor(runif(1,1,1e5))) {
  set.seed(seed)
  beta = params["beta"]
  alpha = params["alpha"]
  gamma = params["gamma"]
  varG = params["varG"]
  b = params["b"]
  bs = params["bs"]
  d = params["d"]
  ds = params["ds"]
  
  S=x["S"]
  I=x["I"]
  R=x["R"]
  
  ## draw the traits of our infected individuals
  gamma_i <- runif(I, .15, .35)
  
  # start at time 0
  t <- 0 
  
  results <- vector(mode='list', length=100000)
  results[[1]] <- list(t, S, gamma_i, R) 
  
  ## row counter
  i <- 2
  
  #start algorithm
  while(t < tmax & length(gamma_i)>0 ){ 
    irate = beta*S*length(gamma_i)
    rrate = gamma_i
    brate <- (b - bs*(S+R+length(gamma_i))) * (S+R+length(gamma_i))
    drateS <-S*(d)
    drateI <-length(gamma_i)*(d+alpha)
    drateR <- R*(d)
    
    rates<-c(rrate,irate,brate,drateS,drateI,drateR)
    #print(rates)
    ## what time does the event happen?
    dt <- rexp(1, rate=sum(rates))
    
    ## update t 
    t <- t + dt
    
    ## "wheel of fortune"
    wheel <- cumsum(rates)/sum(rates)
    
    ## which event happens? Draw a random uniform to determine
    rand <- runif(1)
    
    ## if event==1, recovery
    ## if event==2, infection
    ## if event==3, birth
    ## if event==4, death of S
    ## if event==5, death of I
    ## if event==6, death of R
    
    #(rrate,irate,brate,drateS,drateI,drateR)
    event <- 1 + sum(rand > wheel) 
    if (event%in%1:length(gamma_i)){### recovery
      gamma_i <- gamma_i[-event]
      R <-R+1
    }
    
    else if(event==length(gamma_i)+1){ # infection
      S<- S-1
      # If you want to let things evolve, use these lines
      #this_gamma_r <- gamma_r[sample(1:length(gamma_r), 1)]
      #new_gamma_r <- pick_individuals(1, traitmean=this_gamma_r,traitsd=sqrt(varG))
      # If you don't want evolution, use this line
      new_gamma_i <- runif(1, .15, .35)
      gamma_i <- c(gamma_i, new_gamma_i)
    }
    else if (event==length(gamma_i)+2){ ### birth
      S <- S+1
    }
    else if (event==length(gamma_i)+3) { ### death of S
      S <- S-1
    }
    else if (event==length(gamma_i)+4) {### death of I (not caused by infection)
      gamma_i <- gamma_i[-sample(1:length(gamma_i),1)]
    }
    else    {   ### death of R
      R <- R-1
    }
    
    results[[i]] <- list(t, S,gamma_i , R)
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

