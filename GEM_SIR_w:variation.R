### gillespie model with variation, no evolution

pick_individuals <- function(N0, traitmean, traitsd) {
  meanlog <- log(traitmean^2/sqrt(traitsd^2+traitmean^2))
  sdlog <- sqrt(log(traitsd^2/traitmean^2+1))
  return(rlnorm(N0, meanlog=meanlog, sdlog=sdlog))
}


gillespie.SIR4 <- function(tmax, params, x, seed=floor(runif(1,1,1e5))) {
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
  beta_i <- pick_individuals(I, traitmean=beta, traitsd=sqrt(varB))
  
  # start at time 0
  t <- 0 
  
  results <- vector(mode='list', length=100000)
  results[[1]] <- list(t, S, beta_i, R)
  
  ## row counter
  i <- 2
  
  #start algorithm
  while(t < tmax & length(beta_i) > 0){
    irate = beta_i*S
    dmrate = alpha*beta_i^2 # split between beta values
    rrate = gamma*length(beta_i)
    brate <- (b - bs*(S+length(beta_i)+R)) * (S+length(beta_i)+R)
    drateS <-S*(d+ds*(S+length(beta_i)+R))
    drateI <-length(beta_i)*(d+ds*(S+length(beta_i)+R))
    drateR <- R*(d+ds*(S+length(beta_i)+R))
    
    
    rates<-c(irate,dmrate,rrate,brate,drateS,drateI,drateR)   #irate=infection, dmrate=disease mort, rrate=recovery
    
    ## what time does the event happen?
    dt <- rexp(1, rate=sum(rates))
    
    ## update t 
    t <- t + dt
    
    ## "wheel of fortune"
    wheel <- cumsum(rates)/sum(rates)
    
    ## which event happens? Draw a random uniform to determine
    rand <- runif(1)
    
    ## if event==1, infection
    ## if event==2, disease mort
    ## if event==3, recovery
    ## if event==4, birth
    ## if event==5, death
    event <- 1 + sum(rand > wheel) 
    if (event%in%1:length(beta_i)){### infection event
      S <- S-1
      beta_i <- c(beta_i, pick_individuals(1, traitmean=beta_i[event], traitsd=sqrt(varB)))
    }
    
    else if (event%in%((length(beta_i)+1):(2*length(beta_i)))) {
      beta_i <- beta_i[-(event-length(beta_i))]
    }
    else if (event==2*length(beta_i)+1){ ### recover, randomly remove a beta_i
      R <- R+1
      beta_i <- beta_i[-sample(1:length(beta_i),1)]
    }
    else if (event==2*length(beta_i)+2){ ### birth
      S <- S+1
    }
    else if (event==2*length(beta_i)+3) { ### death of S
      S <- S-1
    }
    else if (event==2*length(beta_i)+4) {### death of I (not caused by infection)
      beta_i <- beta_i[-sample(1:length(beta_i),1)]
      
    }
    else    {   ### death of R
      R <- R-1 
      
    }
    
    results[[i]] <- list(t, S, beta_i, R)
    i <- i +1
    
  }  
  
  results <- results[1:(i-1)]
  return(results)
}



## output

x = c(S=70, I=10, R=0)
tmax <- 400
params4 = c(beta=.3, alpha=.3, gamma=.3, varB=1e-3, b=2, d=0.4, bs=0.01, ds=0.01)

out11 <- gillespie.SIR3(tmax, params4, x, 2354324)

seeds <- floor(runif(20,1,1e5)) # set seeds
library(parallel)
detectCores() # detect how many cores there are
mclapply(seeds,
         function(s) gillespie.SIR3(tmax, params4, x, s),
         mc.cores=4) -> out

## Compute number of individuals at times 0, 1, 2, ... 400
timeSeq <- 0:400
storeMatrix <- array(NA, dim=c(length(timeSeq),length(out)))
for (j in 1:length(out)) {
  o <- out[[j]]
  recordTimes <- lapply(o, function(o2) o2[[1]]) %>% unlist
  storeVector <- rep(0, length=length(timeSeq))
  listInds <- sapply(timeSeq, function(t) which(recordTimes >= t) %>% min)
  if (any(is.infinite(listInds)))
    listInds <- listInds[-which(is.infinite(listInds))] 
  storeVector[1:length(listInds)] <- sapply(listInds, function(i) mean(o[[i]][[3]])) ## number of infecteds
  storeMatrix[,j] <- storeVector
}

plot.new()
plot.window(xlim=c(0,400), ylim=c(0,2.3))
axis(1);axis(2);box('plot')
for (i in 1:ncol(storeMatrix)) lines(0:400, storeMatrix[,i], col=gray(0.4))
lines(0:400, apply(storeMatrix, 1, mean), col=2, lwd=2)

lapply(out, function(o) lapply(o, function(o2) o2[[1]]) %>% unlist)


#### output
plot(unlist(lapply(out11, function(y) y[[1]])), unlist(lapply(out11, function(y) length(y[[3]]))), col="red",
     type='l',lwd=1,xlab='Time', ylab="N", ylim=c(-5,100), xlim=c(0,400)) #main="beta=0.1, alpha=0.1, gamma=0.5,varB=1e-3, b=2.5, d=0.4")
# plot number of susceptible
lines(unlist(lapply(out11, function(y) y[[1]])), unlist(lapply(out11, function(y) y[[2]])),col="blue", 
      type='l',lwd=1, xlab="Time", ylab="Number susceptible")
# plot number of recovered
lines(unlist(lapply(out11, function(y) y[[1]])), unlist(lapply(out11, function(y) y[[4]])), col="green",
      type='l',lwd=1, xlab="Time", ylab="Number recovered")
legend("topright",legend=c("S","I","R"),fill=c("blue","red","green"))