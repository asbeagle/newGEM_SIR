### gillespie model with stratified variation in beta

pick_individuals <- function(N0, traitmean, traitsd) {
  meanlog <- log(traitmean^2/sqrt(traitsd^2+traitmean^2))
  sdlog <- sqrt(log(traitsd^2/traitmean^2+1))
  return(rlnorm(N0, meanlog=meanlog, sdlog=sdlog))
}


gillespie.SIR.strat.varB <- function(tmax, params, x, seed=floor(runif(1,1,1e5))) {
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
  beta_values<-c(.15,.25,.35)
  beta_i<-sample(beta_values,10, replace=TRUE)
  
  
  # start at time 0
  t <- 0 
  
  results <- vector(mode='list', length=100000)
  results[[1]] <- list(t, S, beta_i, R)
  
  ## row counter
  i <- 2
  
  #start algorithm
  while(t < tmax & length(beta_i) > 0){
    irate = beta_i*S
    rrate = gamma*length(beta_i)
    brate <- (b - bs*(S+length(beta_i)+R)) * (S+length(beta_i)+R)
    drateS <-S*(d)
    drateI <-length(beta_i)*(d+alpha)
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
    if (event%in%1:length(beta_i)){### infection event
      S <- S-1
      beta_i <- c(beta_i, sample(beta_values, 1, replace=TRUE))
    }
    
    else if (event==length(beta_i)+1){ ### recover, randomly remove a beta_i
      R <- R+1
      beta_i <- beta_i[-sample(1:length(beta_i),1)]
    }
    else if (event==length(beta_i)+2){ ### birth
      S <- S+1
    }
    else if (event==length(beta_i)+3) { ### death of S
      S <- S-1
    }
    else if (event==length(beta_i)+4) {### death of I
      beta_i <- beta_i[-sample(1:length(beta_i),1)]
    }
    else    {   ### death of R
      R <- R-1 
    }
    
    results[[i]] <- list(t, S, beta_i, R)
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
params4 = c(beta=.7, alpha=.01, gamma=.3, varB=1e-3, b=2, d=0.4, bs=0.01, ds=0.01)

out15 <- gillespie.SIR.strat.varB(tmax, params4, x)

plot.ts(out15[,2], col="blue", ylim=c(-5, 200))
lines(out15[,3], col="red")
lines(out15[,4], col="green")

## plot S, I, R
plot(unlist(lapply(out15, function(y) y[[1]])), unlist(lapply(out15, function(y) length(y[[3]]))), col="red",
     type='l',lwd=1,xlab='Time', ylab="N", ylim=c(0,300), xlim=c(0,20)) #main="beta=0.1, alpha=0.1, gamma=0.5,varB=1e-3, b=2.5, d=0.4")
# plot number of susceptible
lines(unlist(lapply(out15, function(y) y[[1]])), unlist(lapply(out15, function(y) y[[2]])),col="blue", 
      type='l',lwd=1, xlab="Time", ylab="Number susceptible")
# plot number of recovered
lines(unlist(lapply(out15, function(y) y[[1]])), unlist(lapply(out15, function(y) y[[4]])), col="green",
      type='l',lwd=1, xlab="Time", ylab="Number recovered")
legend("topright",legend=c("S","I","R"),fill=c("blue","red","green"))


# multiple runs

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

