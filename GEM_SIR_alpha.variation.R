### gillespie model with variation in alpha, no evolution

pick_individuals <- function(N0, traitmean, traitsd) {
  meanlog <- log(traitmean^2/sqrt(traitsd^2+traitmean^2))
  sdlog <- sqrt(log(traitsd^2/traitmean^2+1))
  return(rlnorm(N0, meanlog=meanlog, sdlog=sdlog))
}



gillespie.SIR.varA <- function(tmax, params, x, seed=floor(runif(1,1,1e5))) {
  set.seed(seed)
  beta = params["beta"]
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
    irate = length(alpha_i)*S*beta
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
      
      alpha_i <- alpha_i[-sample(1:length(alpha_i),1)]
      
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

## output
x = c(S=70, I=10, R=0)
tmax <- 150
params4 = c(beta=.25, alpha=.01, gamma=.15, varG=1e-3, b=2, d=0.4, bs=0.01, ds=0.01, varA=1e-3)

out14 <- gillespie.SIR.var.A(tmax, params4, x)

plot.ts(out14[,2], col="blue", ylim=c(-5, 200))
lines(out14[,3], col="red")
lines(out14[,4], col="green")