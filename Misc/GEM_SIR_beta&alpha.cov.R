## covariance in beta & alpha SIR GEM covary

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


gillespie.SIR.covB_A <- function(tmax, params, corr, x, seed=floor(runif(1,1,1e5))) {
  set.seed(seed)
  beta = params["beta"]
  alpha = params["alpha"]
  gamma = params["gamma"]
  sd_B=params["sd_B"]
  sd_A=params["sd_A"]
  correlation=corr
  b = params["b"]
  bs = params["bs"]
  d = params["d"]
  
  S=x["S"]
  I=x["I"]
  R=x["R"]
  
  ## draw the traits of our infected individuals
  new_i <- pick_individuals_multivariate(I, traitmeans=c(beta,alpha), traitsds=c(sd_B,sd_A), corr=correlation) #infection
  print(new_i)
  beta_i<-new_i[,1]
  alpha_i <-new_i[,2]
  
  # start at time 0
  t <- 0 
  
  results <- vector(mode='list', length=100000)
  results[[1]] <- list(t, S, beta_i, alpha_i, R)
  
  ## row counter
  i <- 2
  
  #start algorithm
  while(t < tmax & length(beta_i) >0 ){ 
    irate = beta_i*S
    rrate = gamma*(length(beta_i))
    brate <- (b - bs*(S+R+(length(beta_i)))) * (S+R+(length(beta_i)))
    drateS <-S*(d)
    drateI <-(d+alpha_i)
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
    if (event%in%1:length(beta_i)){### infection
      S <- S-1
      infection <- pick_individuals_multivariate(1,traitmeans=c(beta, alpha), traitsds=c(sd_B, sd_A), corr=correlation)
      beta_i <- c(beta_i, infection[,1]) # add to list of beta i
      alpha_i<-c(alpha_i, infection[,2]) # add to list of alpha i
    }
    else if(event%in%((length(beta_i)+1):(length(beta_i)+length(alpha_i)))){ # death of I
      beta_i <- beta_i[-(event-length(alpha_i))]
      alpha_i<- alpha_i[-(event-length(alpha_i))]
    } 
      
    else if(event==(2*length(beta_i)+1)){ # recovery
      R <- R+1
      ind=sample(1:length(beta_i),1)
      beta_i<-beta_i[-(ind)] 
      alpha_i<-alpha_i[-(ind)]
    }

    else if (event==(2*length(beta_i)+2)){ ### birth
      S <- S+1
    }
    else if (event==(2*length(beta_i)+3)) { ### death of S
      S <- S-1
    }
    else    {   ### death of R
      R <- R-1
    }
    
    results[[i]] <- list(t, S,beta_i, alpha_i , R)
    i <- i +1
    
  }  
  
  results <- results[1:(i-1)]
  return(results)
}

x = c(S=70, I=10, R=0)
tmax <- 150
params5 = c(beta=.1, alpha=.01, gamma=.3, sd_B=.01, sd_A=.001, b=2, d=0.4, bs=0.01)
corr <- matrix(c(1,0,0,1), nrow=2, byrow=T)

out17 <- gillespie.SIR.covB_A(tmax, params5, corr, x)
length(out17)


# compare pick individuals functions, look for other package?