## covariance in contacts and gamma (sigmoidal)
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


gillespie.SIR.cov_cgamma <- function(tmax, params, corr, x, seed=floor(runif(1,1,1e5))) {
  set.seed(seed)
  c = params["c"]
  shed = params["shed"]
  alpha = params["alpha"]
  gamma = params["gamma"]
  sd_c=params["sd_c"]
  sd_g=params["sd_g"]
  correlation=corr
  b = params["b"]
  bs = params["bs"]
  d = params["d"]
  
  S=x["S"]
  I=x["I"]
  R=x["R"]
  
  ## draw the traits of our infected individuals
  new_i <- pick_individuals_multivariate(I, traitmeans=c(c,gamma), traitsds=c(sd_c,sd_g), corr=correlation) #infection
  #print(new_i)
  c_i<-new_i[,1]
  gamma_i <- new_i[,2]
  
  # start at time 0
  t <- 0 
  
  results <- vector(mode='list', length=100000)
  results[[1]] <- list(t, S, c_i, gamma_i, R)
  
  ## row counter
  i <- 2
  
  #start algorithm
  while(t < tmax & length(c_i) >0 ){ 
    irate = c_i*((shed)/(1+shed))*S
    rrate = gamma_i
    brate <- (b - bs*(S+R+(length(c_i)))) * (S+R+(length(c_i)))
    drateS <-S*(d)
    drateI <-(d+alpha)*(length(c_i))
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
    if (event%in%1:length(c_i)){### infection
      S <- S-1
      infection <- pick_individuals_multivariate(1,traitmeans=c(c, gamma), traitsds=c(sd_c, sd_g), corr=correlation)
      c_i <- c(c_i, infection[,1]) # add to list of c i
      gamma_i <- c(gamma_i, infection[,2]) # add to list of gamma i
    }
    else if(event%in%((length(c_i)+1):(length(c_i)+length(gamma_i)))){ # recovery
      R <- R+1
      c_i<-c_i[-(event-length(c_i))]
      gamma_i<-gamma_i[-(event-length(gamma_i))] 
    }
    else if(event==(2*(length(c_i)))+1){ # death of I
      ind=sample(1:length(c_i),1)
      c_i<-c_i[-(ind)] 
      gamma_i<-gamma_i[-(ind)]
    } 
    
    else if (event==(2*(length(c_i)))+2){ ### birth
      S <- S+1
    }
    else if (event==(2*(length(c_i)))+3) { ### death of S
      S <- S-1
    }
    else    {   ### death of R
      R <- R-1
    }
    
    results[[i]] <- list(t, S, c_i, gamma_i , R)
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



#x = c(S=70, I=10, R=0)
#tmax <- 150

#cov_parms = c(c=.1, shed=.05, sd_c=.01, sd_g=.1, alpha=.1, gamma=.1, b=2.5, d=.1, bs=.01)
#corr <- matrix(c(1,0,0,1), nrow=2, byrow=T)


#cov_out2 <- gillespie.SIR.cov_cgamma(tmax, cov_parms, corr, initial_state)

#par(mfrow=c(1,1))
#plot.ts(cov_out2[,2], col="blue", ylim=c(-5, 250))
#lines(cov_out2[,3], col="red")
#lines(cov_out2[,4], col="green")

