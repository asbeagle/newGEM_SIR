pick_individuals_multivariate <-function(N, traitmeans, traitsds, corr, seed=floor(runif(1,1,1e5))) {
  set.seed(seed)
  
  ## Generate Lognormal random variables with zeros. 
  ## traitmeans - Vector of mean trait values
  ## traitsds - Vector of trait standard deviations
  ## R - Correlation matrix - should have 1s along the diagonal and off-diagonal elements reflecting correlations between different variables
  
  ## Number of traits
  p <- length(traitmeans);
  mu <- matrix(0,p,1);
  sigma <- matrix(0,p,p);
  for (i in 1:p) {
    mu[i] = log(traitmeans[i]^2 / sqrt(traitsds[i]^2+traitmeans[i]^2))
    sigma[i,i] = log(1 + traitsds[i]^2/traitmeans[i]^2);
  }
  
  Cov <- sqrt(sigma)%*%corr%*%sqrt(sigma);
  X <- matrix(rnorm(p*N,0,1), nrow = N, ncol = p);
  for (i in 1:N){
    X[i,] <- exp(mu + chol(Cov)%*%X[i,])
  }
  colnames(X) <- names(traitmeans) ## give the columns names
  return(X);
}

## New version that returns information on the traits and fitness of every infected individual
## This version can be used to simulate an SIR with covariance between *any* two parameters, 
## specified by the new input variable 'covParams' which must be a vector of two of the 
## four parameters ('c', 'shed', 'alpha', 'gamma'). Thus you must also specify standard 
## deviations for all four parameters in the 'params' input variable, and these *must*
## be specified with the full parameter name (e.g., "sd_alpha" rather than "sd_a")
##
## This model differs from the previous model in that the disease does not become endemic: 
## because there is no replenishment of susceptibles, the pathogen eventually goes extinct
gillespie.SIR.cov_storage <- function(tmax, params, corr, x, covParams, seed=floor(runif(1,1,1e5))) {
  set.seed(seed)

  S=x["S"]
  I=x["I"]
  R=x["R"]
  
  ## Setting up storage is actually pretty important for this model because so much information is being retained that it can slow computation a lot
  ## In particular, you need to store data for every infected individual that ever lived, not just the individuals that are currently alive
  ## It's also important to store beginning and end times for each infection - that way you can actually calculate how many individuals are alive at any moment in time.
  ## Because we want to keep track of how many infections are caused by each infected person, we need to give every infected individual a full suite of traits (shedding, contact, virulence, and recovery), even though only two of those will vary at any one time. 
  storage <- data.frame(ID=1:1e5, c=rep(NA,1e5), shed=rep(NA,1e5), alpha=rep(NA,1e5), gamma=rep(NA,1e5), numInf=rep(0,1e5), tInf=rep(NA,1e5), tEnd=rep(Inf,1e5))
  
  ## also set up some storage that just keeps track of the population sizes at every moment in time
  popSizes <- vector(mode='list', 1e6)
  
  ## draw the covarying traits of our initial infected individuals and store them
  new_i <- pick_individuals_multivariate(I, traitmeans=params[covParams], traitsds=params[paste('sd',covParams,sep="_")], corr=corr) 
  storage[1:I,c('c','shed','alpha','gamma','tInf')] <- data.frame(c=ifelse(rep('c' %in% covParams,I), new_i[,'c'], rep(params['c'],I)),
                                                                  shed=ifelse(rep('shed' %in% covParams,I), new_i[,'shed'], rep(params['shed'],I)),
                                                                  alpha=ifelse(rep('alpha' %in% covParams,I), new_i[,'alpha'], rep(params['alpha'],I)),
                                                                  gamma=ifelse(rep('gamma' %in% covParams,I), new_i[,'gamma'], rep(params['gamma'],I)),
                                                                  tInf=rep(0,I))
  
  ## start at time 0
  t <- 0 
  ## ID counter (the ID to assign to the next infection)
  ID <- I+1
  ## popSize list counter
  i <- 1
  
  #start algorithm
  while(t < tmax & I > 0 ) { 
    #print(t)
    ## count how many infected people there are right now - these are individuals with birth times <= the current time and death times that are still Inf
    ## Do this by first collecting all of their IDs
    aliveID <- storage[which(with(storage, tInf<=t & is.infinite(tEnd))),'ID']
    #print(aliveID)
    ## number alive is just the length of that vector
    I <- length(aliveID)
    ## Store the current population size and increment the counter
    popSizes[[i]] <- c(t, S, I, R)
    i <- i + 1
    
    ## Calculate the rates
    irate = storage[aliveID,'c']*storage[aliveID,'shed']/(1 + storage[aliveID,'shed'])*S ## individual infection rates
    rrate = storage[aliveID,'gamma'] ## individual recovery rates
    drateI <- storage[aliveID,'alpha'] ## individual mortality rates
    ## Put all the rates into a single vector
    rates<-c(irate,rrate,drateI)
    ## what time does the event happen?
    dt <- rexp(1, rate=sum(rates))
    ## update t 
    t <- t + dt
    ## "wheel of fortune"
    wheel <- cumsum(rates)/sum(rates)
    ## which event happens? Draw a random uniform to determine
    rand <- runif(1)
    event <- 1 + sum(rand > wheel) 
    #print(event)
    if (event%in%1:I){ ## the first I events are infections 
      S <- S-1
      ## the individual who caused the infection has ID equal to aliveID[event] - increase their numInf by 1
      storage[aliveID[event],'numInf'] <- storage[aliveID[event],'numInf'] + 1
      ## generate traits for the new individual
      new_i <- pick_individuals_multivariate(1, traitmeans=params[covParams], traitsds=params[paste('sd',covParams,sep="_")], corr=corr) 
      ## assign it to the next ID and store it
      storage[ID,c('c','shed','alpha','gamma','tInf')] <- data.frame(c=ifelse('c' %in% covParams, new_i[,'c'], params['c']),
                                                                     shed=ifelse('shed' %in% covParams, new_i[,'shed'], params['shed']),
                                                                     alpha=ifelse('alpha' %in% covParams, new_i[,'alpha'], params['alpha']),
                                                                     gamma=ifelse('gamma' %in% covParams, new_i[,'gamma'], params['gamma']),
                                                                     tInf=t)
      ## increment ID by 1
      ID <- ID+1
    }
    else if(event%in%((I+1):(2*I))) { ## events I+1:2*I are all recoveries
      R <- R+1
      ## the individual who recovered has ID equal to aliveID[event-I] - set their tEnd to the current time
      storage[aliveID[event-I],'tEnd'] <- t
    } 
    
    else { ## events 2*I+1:3*I are all deaths
      ## the individual who recovery has ID equal to aliveID[event-2*I] - set their tEnd to the current time
      storage[aliveID[event-2*I],'tEnd'] <- t
    }
    
    if (ID > nrow(storage)) ## add more rows to storage if necessary
      storage <- rbind(storage, data.frame(ID=(1e5+1):2e5, c=rep(NA,1e5), shed=rep(NA,1e5), alpha=rep(NA,1e5), gamma=rep(NA,1e5), numInf=rep(0,1e5), tInf=rep(NA,1e5), tEnd=rep(Inf,1e5)))
    
    ## Figure out how many individuals are alive now so you can break the while loop if the epidemic has ended
    aliveID <- storage[which(with(storage, tInf<=t & is.infinite(tEnd))),'ID']
    #print(aliveID)
    ## number alive is just the length of that vector
    I <- length(aliveID)
    
    
    #print(storage[1:ID,])
    #readline(prompt="Press [enter] to continue")
  }
  ## Store the final population size
  popSizes[[i]] <- c(t, S, I, R)
  popSizes <- popSizes[1:i]
  popSizes <- do.call('rbind.data.frame',popSizes)
  colnames(popSizes) <- c('t','S','I','R')
  ## Trim this down to store the population sizes at particular points in time (for now, just every timestep)
  ## Find the population sizes just before timepoint 1, 2, 3, etc. and store those
  popSizes <- popSizes[c(1,sapply(1:ceiling(t), function(T) max(which(popSizes$t < T)))),]
  popSizes[,1] <- 0:ceiling(t)
  ## Trim the blank row from storage
  storage <- storage[1:ID,]
  results <- list(popSizes, storage)
  return(results)
}

