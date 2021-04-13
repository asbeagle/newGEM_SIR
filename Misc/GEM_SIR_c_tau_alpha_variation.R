pick_individuals <- function(N0, traitmean, traitsd) {
  meanlog <- log(traitmean^2/sqrt(traitsd^2+traitmean^2))
  sdlog <- sqrt(log(traitsd^2/traitmean^2+1))
  return(rlnorm(N0, meanlog=meanlog, sdlog=sdlog))
}

gillespie.SIR.var_c <- function(tmax, params, x, seed=floor(runif(1,1,1e5)), family=c("lognormal","normal")) {
  set.seed(seed)
  c = params["c"]
  shed = params["shed"]
  h = params["h"]
  alpha = params["alpha"]
  gamma = params["gamma"]
  sd_c=params["sd_c"]
  b = params["b"]
  bs = params["bs"]
  d = params["d"]
  
  S=x["S"]
  I=x["I"]
  R=x["R"]
  
  ## draw the traits of our infected individuals
  if (family=="lognormal")
    c_i <- pick_individuals(I, traitmean=c, traitsd=sd_c) #infection
  else if (family=="normal") {
    c_i <- rnorm(I, mean=c, sd=sd_c)
    ## truncate to ensure positive trait values
    c_i[c_i < 0] = 1e-5
  }
    
  
  # start at time 0
  t <- 0 
  
  results <- vector(mode='list', length=100000)
  results[[1]] <- list(t, S, c_i, R)
  
  ## row counter
  i <- 2
  
  #start algorithm
  while(t < tmax & length(c_i) >0 ){ 
    irate = c_i*shed/(h+shed)*S
    rrate = gamma*(length(c_i))
    brate <- (b - bs*(S+R+(length(c_i)))) * (S+R+(length(c_i)))
    drateS <-S*(d)
    drateI <-(d+alpha)*(length(c_i))
    drateR <- R*(d)
    
    rates<-c(irate,drateI,rrate,brate,drateS,drateR)
    #print(rates)
    ## what time does the event happen?
    dt <- rexp(1, rate=sum(rates))
    
    ## update t 
    t <- t + dt
    
    ## "wheel of fortune"
    wheel <- cumsum(rates)/sum(rates)
    
    ## which event happens? Draw a random uniform to determine
    rand <- runif(1)
    
    event <- 1 + sum(rand > wheel) 
    if (event%in%1:length(c_i)){### infection
      S <- S-1
      if (family=="lognormal")
        c_i <- c(c_i, pick_individuals(1, traitmean=c, traitsd=sd_c)) # add to list of c i
      else if (family=="normal") ## truncate to ensure positive trait values
        c_i <- c(c_i, max(1e-5, rnorm(1, mean=c, sd=sd_c)))
    }
    else if(event==((length(c_i)+1))){ # death of I
      ind=sample(1:length(c_i),1)
      c_i<-c_i[-(ind)] 
    } 
    else if(event==(length(c_i)+2)){ # recovery
      R <- R+1
      ind=sample(1:length(c_i),1)
      c_i<-c_i[-(ind)] 
    }
    else if (event==(length(c_i)+3)){ ### birth
      S <- S+1
    }
    else if (event==(length(c_i)+4)) { ### death of S
      S <- S-1
    }
    else    {   ### death of R
      R <- R-1
    }
    
    results[[i]] <- list(t, S, c_i, R)
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

gillespie.SIR.var_tau <- function(tmax, params, x, seed=floor(runif(1,1,1e5)), family=c("lognormal","normal")) {
  set.seed(seed)
  c = params["c"]
  shed = params["shed"]
  h = params["h"]
  alpha = params["alpha"]
  gamma = params["gamma"]
  sd_s=params["sd_s"]
  correlation=corr
  b = params["b"]
  bs = params["bs"]
  d = params["d"]
  
  S=x["S"]
  I=x["I"]
  R=x["R"]
  
  ## draw the traits of our infected individuals
  if (family=="lognormal")
    shed_i <- pick_individuals(I, traitmean=shed, traitsd=sd_s)
  else if (family=="normal") {
    shed_i <- rnorm(I, mean=shed, sd=sd_s)
    ## truncate to ensure positive
    shed_i[shed_i < 0] <- 1e-5
  }
  tau_i <- shed_i/(h+shed_i)
  
  # start at time 0
  t <- 0 
  
  results <- vector(mode='list', length=100000)
  results[[1]] <- list(t, S, tau_i, R)
  
  ## row counter
  i <- 2
  
  #start algorithm
  while(t < tmax & length(tau_i) >0 ){ 
    irate = c*tau_i*S
    rrate = gamma*(length(tau_i))
    brate <- (b - bs*(S+R+(length(tau_i)))) * (S+R+(length(tau_i)))
    drateS <-S*(d)
    drateI <-(d+alpha)*(length(tau_i))
    drateR <- R*(d)
    
    rates<-c(irate,drateI,rrate,brate,drateS,drateR)
    ## what time does the event happen?
    dt <- rexp(1, rate=sum(rates))
    
    ## update t 
    t <- t + dt
    
    ## "wheel of fortune"
    wheel <- cumsum(rates)/sum(rates)
    
    ## which event happens? Draw a random uniform to determine
    rand <- runif(1)
    
    event <- 1 + sum(rand > wheel) 
    if (event%in%1:length(tau_i)){### infection
      S <- S-1
      if (family=="lognormal")
        shed_i <- pick_individuals(1, traitmean=shed, traitsd=sd_s)
      else if (family=="normal")
        shed_i <- max(1e-5, rnorm(1, mean=shed, sd=sd_s)) ## truncate to maintain positive
      tau_i <- c(tau_i, shed_i/(h+shed_i)) # add to list of c i
    }
    else if(event==((length(tau_i)+1))){ # death of I
      ind=sample(1:length(tau_i),1)
      tau_i<-tau_i[-(ind)] 
    } 
    else if(event==(length(tau_i)+2)){ # recovery
      R <- R+1
      ind=sample(1:length(tau_i),1)
      tau_i<-tau_i[-(ind)] 
    }
    else if (event==(length(tau_i)+3)){ ### birth
      S <- S+1
    }
    else if (event==(length(tau_i)+4)) { ### death of S
      S <- S-1
    }
    else    {   ### death of R
      R <- R-1
    }
    
    results[[i]] <- list(t, S, tau_i, R)
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

gillespie.SIR.var_alpha <- function(tmax, params, x, seed=floor(runif(1,1,1e5)), family=c("lognormal","normal")) {
  set.seed(seed)
  c = params["c"]
  shed = params["shed"]
  h = params["h"]
  alpha = params["alpha"]
  gamma = params["gamma"]
  sd_a=params["sd_a"]
  correlation=corr
  b = params["b"]
  bs = params["bs"]
  d = params["d"]
  
  S=x["S"]
  I=x["I"]
  R=x["R"]
  
  ## draw the traits of our infected individuals
  if (family=="lognormal")
    alpha_i <- pick_individuals(I, traitmean=alpha, traitsd=sd_a) #infection
  else if (family=="normal") {
    alpha_i <- rnorm(I, mean=alpha, sd=sd_a)
    alpha_i[alpha_i < 0] <- 1e-5 ## truncate to keep positive
  }

  # start at time 0
  t <- 0 
  
  results <- vector(mode='list', length=100000)
  results[[1]] <- list(t, S, alpha_i, R)
  
  ## row counter
  i <- 2
  
  #start algorithm
  while(t < tmax & length(alpha_i) >0 ){ 
    irate = c*S*length(alpha_i)
    rrate = gamma*(length(alpha_i))
    brate <- (b - bs*(S+R+(length(alpha_i)))) * (S+R+(length(alpha_i)))
    drateS <-S*(d)
    drateI <-(d+alpha_i)
    drateR <- R*(d)
    
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
    if (event%in%1:length(alpha_i)){### death of I
      alpha_i <- alpha_i[-event]
     }
    else if(event==((length(alpha_i)+1))){ # infection
      S <- S-1
      if (family=="lognormal")
        alpha_i <- c(alpha_i, pick_individuals(1, traitmean=alpha, traitsd=sd_a))
      else if (family=="normal")
        alpha_i <- c(alpha_i, max(1e-5, rnorm(1, mean=alpha, sd=sd_a))) ## truncated normal
    } 
    else if(event==(length(alpha_i)+2)){ # recovery
      R <- R+1
      ind=sample(1:length(alpha_i),1)
      alpha_i<-alpha_i[-(ind)] 
    }
    else if (event==(length(alpha_i)+3)){ ### birth
      S <- S+1
    }
    else if (event==(length(alpha_i)+4)) { ### death of S
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

