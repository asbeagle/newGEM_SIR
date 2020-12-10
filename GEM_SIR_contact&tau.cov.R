## covariance in contacts and transmissibility as function of shedding
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


gillespie.SIR.cov_csig <- function(tmax, params, corr, x, seed=floor(runif(1,1,1e5))) {
  set.seed(seed)
  c = params["c"]
  shed = params["shed"]
  h = params["h"]
  alpha = params["alpha"]
  gamma = params["gamma"]
  sd_c=params["sd_c"]
  sd_s=params["sd_s"]
  correlation=corr
  b = params["b"]
  bs = params["bs"]
  d = params["d"]
  
  S=x["S"]
  I=x["I"]
  R=x["R"]
  
  ## draw the traits of our infected individuals
  new_i <- pick_individuals_multivariate(I, traitmeans=c(c,shed), traitsds=c(sd_c,sd_s), corr=correlation) #infection
  #print(new_i)
  c_i<-new_i[,1]
  tau_i <- (new_i[,2]^2)/((h^2)+(new_i[,2]^2))
  
  # start at time 0
  t <- 0 
  
  results <- vector(mode='list', length=100000)
  results[[1]] <- list(t, S, c_i, tau_i, R)
  
  ## row counter
  i <- 2
  
  #start algorithm
  while(t < tmax & length(c_i) >0 ){ 
    irate = c_i*tau_i*S
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
    
    #(irate,drateI,rrate,brate,drateS,drateR)
    
    event <- 1 + sum(rand > wheel) 
    if (event%in%1:length(c_i)){### infection
      S <- S-1
      infection <- pick_individuals_multivariate(1,traitmeans=c(c, shed), traitsds=c(sd_c, sd_s), corr=correlation)
      c_i <- c(c_i, infection[,1]) # add to list of c i
      tau_i <- c(tau_i, (infection[,2]^2)/((h^2)+(infection[,2]^2))) # add to list of tau i
    }
    else if(event==((length(c_i)+1))){ # death of I
      ind=sample(1:length(c_i),1)
      c_i<-c_i[-(ind)] 
      tau_i<-tau_i[-(ind)]
    } 
    
    else if(event==(length(c_i)+2)){ # recovery
      R <- R+1
      ind=sample(1:length(c_i),1)
      c_i<-c_i[-(ind)] 
      tau_i<-tau_i[-(ind)]
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
    
    results[[i]] <- list(t, S, c_i, tau_i , R)
    i <- i +1
    
  }  
  
  results <- results[1:(i-1)]
  return(results)
}

x = c(S=70, I=10, R=0)
tmax <- 150

all.new.params = c(c=.02, shed=.2, sd_c=.01, sd_s=.01, h=.1, alpha=.01, gamma=.3, b=2.5, d=.4, bs=.01)
corr <- matrix(c(1,0,0,1), nrow=2, byrow=T)


new_out <- gillespie.SIR.cov_csig(tmax, all.new.params, corr, x)
length(out17)

lapply(new_out, function(l) l[[2]])%>%unlist%>%plot


