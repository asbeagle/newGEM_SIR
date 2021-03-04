### gillespie model with variation in gamma, no evolution

pick_individuals <- function(N0, traitmean, traitsd) {
  meanlog <- log(traitmean^2/sqrt(traitsd^2+traitmean^2))
  sdlog <- sqrt(log(traitsd^2/traitmean^2+1))
  return(rlnorm(N0, meanlog=meanlog, sdlog=sdlog))
}



gillespie.SIR.varG <- function(tmax, params, x, seed=floor(runif(1,1,1e5))) {
  set.seed(seed)
  beta = params["beta"]
  alpha = params["alpha"]
  gamma = params["gamma"]
  varG = params["varG"]
  b = params["b"]
  bs = params["bs"]
  d = params["d"]
  ds = params["ds"]
  
  S=x["S"]
  I=x["I"]
  R=x["R"]
  
  ## draw the traits of our infected individuals
  gamma_i <- pick_individuals(I, traitmean=gamma, traitsd=sqrt(varG))
  
  # start at time 0
  t <- 0 
  
  results <- vector(mode='list', length=100000)
  results[[1]] <- list(t, S, gamma_i, R) 
  
  ## row counter
  i <- 2
  
  #start algorithm
  while(t < tmax & length(gamma_i)>0 ){ 
    irate = beta*S*length(gamma_i)
    rrate = gamma_i
    brate <- (b - bs*(S+R+length(gamma_i))) * (S+R+length(gamma_i))
    drateS <-S*(d)
    drateI <-length(gamma_i)*(d+alpha)
    drateR <- R*(d)
    
    rates<-c(rrate,irate,brate,drateS,drateI,drateR)
    #print(rates)
    ## what time does the event happen?
    dt <- rexp(1, rate=sum(rates))
    
    ## update t 
    t <- t + dt
    
    ## "wheel of fortune"
    wheel <- cumsum(rates)/sum(rates)
    
    ## which event happens? Draw a random uniform to determine
    rand <- runif(1)
    
    ## if event==1, recovery
    ## if event==2, infection
    ## if event==3, birth
    ## if event==4, death of S
    ## if event==5, death of I
    ## if event==6, death of R
    
    #(rrate,irate,brate,drateS,drateI,drateR)
    event <- 1 + sum(rand > wheel) 
    if (event%in%1:length(gamma_i)){### recovery
      gamma_i <- gamma_i[-event]
      R <-R+1
    }
    
     else if(event==length(gamma_i)+1){ # infection
      S<- S-1
      # If you want to let things evolve, use these lines
      #this_gamma_r <- gamma_r[sample(1:length(gamma_r), 1)]
      #new_gamma_r <- pick_individuals(1, traitmean=this_gamma_r,traitsd=sqrt(varG))
      # If you don't want evolution, use this line
      new_gamma_i <- pick_individuals(1, traitmean=gamma, traitsd=sqrt(varG))
      gamma_i <- c(gamma_i, new_gamma_i)
    }
    else if (event==length(gamma_i)+2){ ### birth
      S <- S+1
    }
    else if (event==length(gamma_i)+3) { ### death of S
      S <- S-1
    }
    else if (event==length(gamma_i)+4) {### death of I (not caused by infection)
      gamma_i <- gamma_i[-sample(1:length(gamma_i),1)]
    }
    else    {   ### death of R
      R <- R-1
    }
    
    results[[i]] <- list(t, S,gamma_i , R)
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
params4 = c(beta=.25, alpha=.15, gamma=.15, varG=1e-3, b=2, d=0.4, bs=0.01, ds=0.01)

out14 <- gillespie.SIR.varG(tmax, params4, x)

plot.ts(out14[,2], col="blue", ylim=c(-5, 200))
lines(out14[,3], col="red")
lines(out14[,4], col="green")

## plot S,I,R
plot(unlist(lapply(out14, function(y) y[[1]])), unlist(lapply(out14, function(y) length(y[[3]]))), col="red",
     type='l',lwd=1,xlab='Time', ylab="N") 
# plot number of susceptible
lines(unlist(lapply(out14, function(y) y[[1]])), unlist(lapply(out14, function(y) y[[2]])),col="blue", 
      type='l',lwd=1, xlab="Time", ylab="Number susceptible")
# plot number of recovered
lines(unlist(lapply(out14, function(y) y[[1]])), unlist(lapply(out14, function(y) y[[4]])), col="green",
      type='l',lwd=1, xlab="Time", ylab="Number recovered")
legend("topright",legend=c("S","I","R"),fill=c("blue","red","green"))




# run multiple times
seeds <- floor(runif(20,1,1e5)) # set seeds
library(parallel)
detectCores() # detect how many cores there are
mclapply(seeds,
         function(s) gillespie.SIR5(tmax, params4, x, s),
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
plot(unlist(lapply(out14, function(y) y[[1]])), unlist(lapply(out14, function(y) length(y[[3]]))), col="red",
     type='l',lwd=1,xlab='Time', ylab="N") 
# plot number of susceptible
lines(unlist(lapply(out14, function(y) y[[1]])), unlist(lapply(out14, function(y) y[[2]])),col="blue", 
      type='l',lwd=1, xlab="Time", ylab="Number susceptible")
# plot number of recovered
lines(unlist(lapply(out14, function(y) y[[1]])), unlist(lapply(out14, function(y) y[[4]])), col="green",
      type='l',lwd=1, xlab="Time", ylab="Number recovered")
legend("topright",legend=c("S","I","R"),fill=c("blue","red","green"))
