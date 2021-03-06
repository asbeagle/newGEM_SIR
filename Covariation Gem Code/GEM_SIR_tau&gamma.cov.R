## covariance in shedding/transmissabitly (tau) and gamma (sigmoidal)
library(tidyverse)

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
  return(X);
}

gillespie.SIR.cov_shedgamma <- function(tmax, params, corr, x, seed=floor(runif(1,1,1e5))) {
  set.seed(seed)
  c = params["c"]
  shed = params["shed"]
  alpha = params["alpha"]
  gamma = params["gamma"]
  sd_s=params["sd_s"]
  sd_g=params["sd_g"]
  correlation=corr
  b = params["b"]
  bs = params["bs"]
  d = params["d"]
  
  S=x["S"]
  I=x["I"]
  R=x["R"]
  
  ## draw the traits of our infected individuals
  new_i <- pick_individuals_multivariate(I, traitmeans=c(shed,gamma), traitsds=c(sd_s,sd_g), corr=correlation) #infection
  #print(new_i)
  shed_i<-new_i[,1]
  gamma_i <- new_i[,2]
  
  # start at time 0
  t <- 0 
  
  results <- vector(mode='list', length=100000)
  results[[1]] <- list(t, S, shed_i, gamma_i, R)
  
  ## row counter
  i <- 2
  
  #start algorithm
  while(t < tmax & length(shed_i) >0 ){ 
    irate = c*(shed_i/(1+shed_i))*S
    rrate = gamma_i
    brate <- (b - bs*(S+R+(length(shed_i)))) * (S+R+(length(shed_i)))
    drateS <-S*(d)
    drateI <-(d+alpha)*(length(shed_i))
    drateR <- R*(d)
    
    rates<-c(irate,rrate,drateI,brate,drateS,drateR)
    #print(rates)
    ## what time does the event happen?
    dt <- rexp(1, rate=sum(rates))
    
    ## update t 
    t <- t + dt
    
    ## "wheel of fortune"
    wheel <- cumsum(rates)/sum(rates)
    
    ## which event happens? Draw a random uniform to determine
    rand <- runif(1)
    
    #(irate,drateI,rrate,brate,drateS,drateR)
    
    event <- 1 + sum(rand > wheel) 
    if (event%in%1:length(shed_i)){### infection
      S <- S-1
      infection <- pick_individuals_multivariate(1,traitmeans=c(shed, gamma), traitsds=c(sd_s, sd_g), corr=correlation)
      shed_i <- c(shed_i, infection[,1]) # add to list of shed i
      gamma_i <- c(gamma_i, infection[,2]) # add to list of gamma i
    }
    else if(event%in%((length(shed_i)+1):(length(shed_i)+length(gamma_i)))){ # recovery
      R <- R+1
      shed_i<-shed_i[-(event-length(shed_i))]
      gamma_i<-gamma_i[-(event-length(gamma_i))] 
    }
    else if(event==(2*(length(shed_i)))+1){ # death of I
      ind=sample(1:length(shed_i),1)
      shed_i<-shed_i[-(ind)] 
      gamma_i<-gamma_i[-(ind)]
    } 
    
    else if (event==(2*(length(shed_i)))+2){ ### birth
      S <- S+1
    }
    else if (event==(2*(length(shed_i)))+3) { ### death of S
      S <- S-1
    }
    else    {   ### death of R
      R <- R-1
    }
    
    results[[i]] <- list(t, S, shed_i, gamma_i , R)
    i <- i +1
    
  }  
  results <- results[1:(i-1)]
  ## simplify what you're storing so R doesn't crash because the datafiles are too big
  results <- data.frame(time=sapply(results, function(r) r[[1]]),
                        S=sapply(results, function(r) r[[2]]), 
                        I=sapply(results, function(r) length(r[[3]])),
                        R=sapply(results, function(r) r[[5]]))
  ## only keep the state of the system at times 0, 0.1, 0.2, ..., tmax-0.2, tmax-0.1, tmax
  results <- results[sapply(seq(0,tmax,1), function(t) min(which(results$time >= t))),]
  ## these two steps reduce the size of 'results' ~1000-fold

  return(results)
}



