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

pick_individuals <- function(N0, traitmean, traitsd) {
  meanlog <- log(traitmean^2/sqrt(traitsd^2+traitmean^2))
  sdlog <- sqrt(log(traitsd^2/traitmean^2+1))
  return(rlnorm(N0, meanlog=meanlog, sdlog=sdlog))
}


## To make sure everything is working correctly, I want to generate variation in only one variable. This will allow me to compare:
## (1) John's single-variable lognormal pick_individuals function, which converts the mean and standard deviation on the natural scale to the meanlog and sdlog used by the built-in R function 'rlnorm'
traits1 <- pick_individuals(100000, 0.1, 0.1)

## (2) My multivariate lognormal pick_individuals_multivariate function, which converts the means, standard deviations, and correlations on the natural scale to meanlog (the means of the logarithms) and varlog (the covariance between the logarithms)
traitmeans <- c(0.1,0.1)
traitsds <- c(0.1,1e-8) ## variance only in one trait to facilitate the comparison
corr <- matrix(c(1,0,0,1), nrow=2, byrow=TRUE) ## I specify the correlation matrix
traits2 <- pick_individuals_multivariate(100000, traitmeans, traitsds, corr)

## (3) The multivariate lognormal function rlnorm.rplus in the compositions package. I have to specify the covariance matrix of the logarithms and the means of the logarithms. I do this following the same approach used in the pick_individuals function, since there is no covariance between the two variables.
library(compositions)
logcov <- matrix(c(log(traitsds[1]^2/traitmeans[1]^2+1),0,
                     0,log(traitsds[2]^2/traitmeans[2]^2+1)),
                   nrow=2,byrow=TRUE)
logmeans <- c(log(traitmeans[1]^2/sqrt(traitsds[1]^2+traitmeans[1]^2)),
              log(traitmeans[2]^2/sqrt(traitsds[2]^2+traitmeans[2]^2)))
traits3 <- rlnorm.rplus(100000, meanlog=logmeans, varlog=logcov)

## Plot histograms of the trait distributions
par(mfrow=c(1,3))
hist(traits1)
hist(traits2[,1])
hist(traits3[,1])

mean(traits1)     
mean(traits2[,1]) 
mean(traits3[,1]) 
var(traits1)      
var(traits2[,1])  
var(traits3[,1])  





corr <- matrix(c(1,0,0,1), nrow=2, byrow=TRUE) ## I specify the correlation matrix
traits2 <- pick_individuals_multivariate(100000, traitmeans, traitsds, corr)

logcov <- matrix(c(log(traitsds[1]^2/traitmeans[1]^2+1),log(0.5),
                   log(0.5),log(traitsds[2]^2/traitmeans[2]^2+1)),
                 nrow=2,byrow=TRUE)
logmeans <- c(log(traitmeans[1]^2/sqrt(traitsds[1]^2+traitmeans[1]^2)),
              log(traitmeans[2]^2/sqrt(traitsds[2]^2+traitmeans[2]^2)))
traits3 <- rlnorm.rplus(100000, meanlog=logmeans, varlog=logcov)


test <- function(traitmeans, traitsds, corr) {
  a <- rep(-Inf,2)
  while(any(a < 0)) {
    a <- as.numeric(rmvnorm(1, mean=traitmeans, sigma=matrix(c(traitsds[1]^2,rep(corr*traitsds[1]*traitsds[2],2),traitsds[2]^2),nrow=2,byrow=TRUE)))
  }
  return(a)
}

a <- t(sapply(1:1000, function(i) test(traitmeans=c(0.1,0.1), traitsds=c(0.5,0.5), corr=0.5)))

