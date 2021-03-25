### gillespie model with no variation, just stochasticity

gillespie.SIR.noVar <- function(tmax, params, x, seed=floor(runif(1,1,1e5))) {
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
  #beta_i <- pick_individuals(I, traitmean=beta, traitsd=sqrt(varB))
  
  # start at time 0
  t <- 0 
  
  results <- vector(mode='list', length=100000)
  results[[1]] <- list(t, S, I, R)
  
  ## row counter
  i <- 2
  
  #start algorithm
  while(t < tmax & I >0){
    irate = beta*S*I
    rrate = gamma*I
    brate <- (b - bs*(S+I+R)) * (S+I+R)
    drateS <-S*(d)
    drateI <-I*(d+alpha)
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
    
    ## if event==1, infection
    ## if event==2, recovery
    ## if event==3, birth
    ## if event==4, death of S
    ## if event==5, death of I
    ## if event==6, death of R
    event <- 1 + sum(rand > wheel) 
    if (event==1){### infection event
      S <- S-1
      I <- I+1
    }
    
    else if (event==2){ ### recover
      R <- R+1
      I <- I-1
    }
    else if (event==3){ ### birth
      S <- S+1
    }
    else if (event==4) { ### death of S
      S <- S-1
    }
    else if (event==5) {### death of I
      I <- I-1
    }
    else    {   ### death of R
      R <- R-1 
    }
    
    results[[i]] <- list(t, S, I, R)
    i <- i +1
    
} 
  
  results <- results[1:(i-1)]
  ## simplify what you're storing so R doesn't crash because the datafiles are too big
  results <- data.frame(time=sapply(results, function(r) r[[1]]),
                        S=sapply(results, function(r) r[[2]]), 
                        I=sapply(results, function(r) r[[3]]),
                        R=sapply(results, function(r) r[[4]]))
  ## only keep the state of the system at times 0, 0.1, 0.2, ..., tmax-0.2, tmax-0.1, tmax
  results <- results[sapply(seq(0,tmax,1), function(t) min(which(results$time >= t))),]
  ## these two steps reduce the size of 'results' ~1000-fold
  
  return(results)
}

