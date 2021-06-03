pick_individuals_multivariate <-function(N, traitmeans, traitsds, corr){
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
gillespie.SIR.cov_storage <- function(tmax, params, corr, x, covParams, seed=floor(runif(1,1,1e5))) {
  set.seed(seed)
  ## Set the demographic parameters that do not vary
  b = params["b"]
  bs = params["bs"]
  d = params["d"]
  
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
    drateI <- d + storage[aliveID,'alpha'] ## individual mortality rates
    brate <- (b - bs*(S+I+R)) * (S+I+R) ## population-level birth rate
    drateS <- d*S ## population-level death rates for susceptible and recovered hosts
    drateR <- d*R
    ## Put all the rates into a single vector
    rates<-c(irate,rrate,drateI,brate,drateS,drateR)
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
    
    else if(event%in%((2*I+1):(3*I))) { ## events 2*I+1:3*I are all deaths
      ## the individual who recovery has ID equal to aliveID[event-2*I] - set their tEnd to the current time
      storage[aliveID[event-2*I],'tEnd'] <- t
    }
    else if (event==(3*I+1)) { ## birth of S
      S <- S+1
    }
    else if (event==(3*I+2)) { ## death of S
      S <- S-1
    }
    else {   ### death of R
      R <- R-1
    }
    if (ID > nrow(storage)) ## add more rows to storage if necessary
      storage <- rbind(storage, data.frame(ID=(1e5+1):2e5, c=rep(NA,1e5), shed=rep(NA,1e5), alpha=rep(NA,1e5), gamma=rep(NA,1e5), numInf=rep(0,1e5), tInf=rep(NA,1e5), tEnd=rep(Inf,1e5)))
    
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
  popSizes <- popSizes[c(1,sapply(1:tmax, function(T) max(which(popSizes$t < T)))),]
  ## Trim the blank row from storage
  storage <- storage[1:ID,]
  results <- list(popSizes, storage)
  return(results)
}

params = c(c=.1, shed=.05, alpha=.1, gamma=.12, sd_c=0.05, sd_shed=0.025, sd_alpha=0.05, sd_gamma=0.06, b=2.5, d=.1, bs=.01)
gillespie.SIR.cov_storage(tmax=20, 
                          params=params, 
                          corr=matrix(c(1,0,0,1), nrow=2, byrow=T), 
                          x=c(S=140,I=10,R=0), 
                          covParams=c('c','alpha')) -> example1

gillespie.SIR.cov_storage(tmax=20, 
                          params=params, 
                          corr=matrix(c(1,0,0,1), nrow=2, byrow=T), 
                          x=c(S=140,I=10,R=0), 
                          covParams=c('shed','gamma')) -> example2

library(parallel)
mclapply(1:4, 
         function(i) gillespie.SIR.cov_storage(tmax=20, 
                          params=params, 
                          corr=matrix(c(1,0,0,1), nrow=2, byrow=T), 
                          x=c(S=140,I=10,R=0), 
                          covParams=c('alpha','gamma')),
         mc.cores=4) -> example

## Organizing the population size results for plotting
## Create a dataframe containing time in the first column and the susceptible population sizes in each replicate simulation in the other columns
timeSeq <- 0:tmax
Sdataframe <- array(NA, dim=c(length(timeSeq), length(example)+1))
Sdataframe[,1] <- timeSeq
for (i in 1:length(example)) Sdataframe[,i+1] <- example[[i]][[1]]$S




## Plots you can make
## Easy stuff like population dynamics
plot(example3[[1]]$t, example3[[1]]$S, type='l', xlab='Time', ylab='S(t)')
## More interesting stuff like histograms of fitness (probably should censor any individual that was still alive, since its fitness is incomplete)
with(subset(example3[[2]], !is.infinite(tInf)), hist(numInf, breaks=max(numInf)+1)) ## most individuals have 0 fitness
