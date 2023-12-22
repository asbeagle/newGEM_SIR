## Stochastic run in Rcpp
sourceCpp("SIRcov.cpp")
library(deSolve) 
library(MASS)
library(tidyverse)
library(ggplot2)
library(magrittr)

## Compute the deterministic expectation for peak size for different parameter sets
SIRdet <- function(t, y, params) {
  S = y[1];
  I = y[2];
  R = y[3];
  b = params["b"];
  bs = params["bs"];
  d = params["d"]
  c = params["c"]
  s = params["s"]
  a = params["a"]
  g = params["g"]
  
  dSdt = (b-bs*(S+I+R))*(S+I+R) - c*s/(1+s)*S*I - d*S
  dIdt = c*s/(1+s)*S*I - (d+a+g)*I
  dRdt = g*I - d*R
  
  return(list(c(dSdt,dIdt,dRdt)))
}
detPeaks <- c()
for (R in c(8,4,1)) {
  if (R==8) x <- 0
  else if (R==4) x <- 0.24647
  else x <- 0.594804
  
  params = c(b=2.5, d=.1, bs=.01, c=(1-x)*0.1, s=(1-x)*1/9, a=(1+x)*0.1, g=(1+x)*0.1)
  initial_state <- floor(c(S=unname(((params["b"]-params["d"])/params["bs"]))-5, I=5, R=0))
  
  out <- ode(initial_state, seq(0,100,0.1), SIRdet, params)
  detPeaks <- c(detPeaks, max(out[,3]))
}

## Compute a stochastic, no variation expectation for peak size, fadeout probability, fadeout timing, and distribution of secondary infections
Rseq <- seq(1,8,0.1)
stochExp <- vector(mode='list', length=length(Rseq))
for (j in 1:length(Rseq)) {
  R <- Rseq[j]
  print(R)
  ## baseline parameters
  b = 2.5; bs=0.01; d=0.1; c=0.1; s=1/9; a=0.1; g=0.1
  ## scaling coefficient x to scale c, s, a, and g to achieve particular R value
  x = 1/(2*(-b*c*s+c*d*s-a*bs*R*s-bs*g*R*s)) * (-a*bs*R-bs*g*R-2*b*c*s+2*c*d*s+bs*d*R*s+sqrt((a*bs*R+bs*g*R+2*b*c*s-2*c*d*s-bs*d*R*s)^2-4*(-b*c*s+c*d*s-a*bs*R*s-bs*g*R*s)*(a*bs*R+bs*d*R+bs*g*R-b*c*s+c*d*s+a*bs*R*s+bs*d*R*s+bs*g*R*s)))
  ## simulation parameters
  params = c(b=2.5, bs=.01, d=.1, c=(1-x)*0.1, s=(1-x)*1/9, a=(1+x)*0.1, g=(1+x)*0.1)
  ## initial conditions
  initial_state <- floor(c(S=unname(((params["b"]-params["d"])/params["bs"]))-5, I=5, R=0))
  ## simulate stochastic model 1000 times
  lapply(1:1000, function(i) SIRnovar(params, initial_state, 50)) -> out
  ## Calculate the peak infected population size and the timing of fadeouts from each of the 1000 runs
  lapply(out, 
         function(o) 
           data.frame(peak=max(o[[2]][,2]),
                      fadeout=ifelse(any(o[[2]][,2]==0),1,0),
                      fadeoutT=ifelse(any(o[[2]][,2]==0),o[[1]][min(which(o[[2]][,2]==0))],50)) %>%
           mutate(., ## compute the dispersion parameter for each run that did not fade out
                  disp=ifelse(fadeout==0, glm.nb((o[[3]][is.na(o[[3]][,3]),2])~1)$theta, NA))) %>%
    do.call("rbind.data.frame",.) %>% 
    mutate(rep=1:1000, R0=R) -> stochExp[[j]]
}
stochExp %<>% do.call("rbind.data.frame",.)
saveRDS(stochExp, file="Stochastic_no_variation_simulations.RDS")

p1 <- stochExp %>% 
  ggplot(., aes(x=R0, y=peak)) + 
  stat_summary(fun.data = "median_hilow", geom="ribbon", fill=gray(0.75)) + 
  stat_summary(fun="mean", geom="line", col="red") + 
  ylab("Epidemic peak") + 
  xlab("") +
  theme_bw()

p2 <- stochExp %>% 
  ggplot(., aes(x=R0, y=disp)) + 
  stat_summary(fun.data = "median_hilow", geom="ribbon", fill=gray(0.75)) + 
  stat_summary(fun="mean", geom="line", col="red") + 
  ylab("Dispersion (k)") + 
  xlab("") + 
  theme_bw()

p3 <- stochExp %>% 
  group_by(R0) %>% 
  summarise(fadeprob=sum(fadeout)/1000) %>% 
  ggplot(., aes(x=R0, y=fadeprob)) + 
  geom_line(col="red") + 
  ylab("Prob of fadeout") + 
  xlab("R0") + 
  theme_bw() 

p4 <- stochExp %>% 
  ggplot(., aes(x=R0, y=fadeoutT)) + 
  stat_summary(fun.data = "median_hilow", geom="ribbon", fill=gray(0.75)) + 
  stat_summary(fun="mean", geom="line", col="red") + 
  ylab("Time to fadeout") + 
  xlab("R0") + 
  theme_bw()
  
library(patchwork)
png(filename="./Stochastic_no_variation_expectations.png", width=5, height=5, units='in', res=300)
(p1 + p2)/(p3 + p4)
dev.off()

outC <- outS <- outA <- outG <- vector(mode='list', length=9)
j = 1
for (R in c(1,4,8)) {
  for (CV in c(.2, 1, 5)) {
    print(CV)
    print(R)
    ## baseline parameters
    b = 2.5; bs=0.01; d=0.1; c=0.1; s=1/9; a=0.1; g=0.1
    ## scaling coefficient x to scale c, s, a, and g to achieve particular R value
    x = 1/(2*(-b*c*s+c*d*s-a*bs*R*s-bs*g*R*s)) * (-a*bs*R-bs*g*R-2*b*c*s+2*c*d*s+bs*d*R*s+sqrt((a*bs*R+bs*g*R+2*b*c*s-2*c*d*s-bs*d*R*s)^2-4*(-b*c*s+c*d*s-a*bs*R*s-bs*g*R*s)*(a*bs*R+bs*d*R+bs*g*R-b*c*s+c*d*s+a*bs*R*s+bs*d*R*s+bs*g*R*s)))
    ## simulation parameters
    params = c(b=2.5, bs=0.01, d=0.1, 
                c=(1-x)*0.1, s=(1-x)*1/9, a=(1+x)*0.1, g=(1+x)*0.1,
                SD=CV*(1-x)*0.1)
    ## initial conditions
    initial_state <- floor(c(S=unname(((params["b"]-params["d"])/params["bs"]))-5, I=5, R=0))
    ## simulate!
    lapply(1:1000, function(i) SIRvarC(params, initial_state, 50)) -> outC[[j]]
    lapply(1:1000, function(i) SIRvarS(params, initial_state, 50)) -> outS[[j]]
    lapply(1:1000, function(i) SIRvarA(params, initial_state, 50)) -> outA[[j]]
    lapply(1:1000, function(i) SIRvarG(params, initial_state, 50)) -> outG[[j]]
    j = j+1
  }
}
 
lapply(outC, 
       function(O)
         lapply(O, 
                function(o)
                  data.frame(peak=max(o[[2]][,2]),
                             totalI=sum(o[[3]][,3]),
                             fadeout=ifelse(any(o[[2]][,2]==0),1,0),
                             fadeoutT=ifelse(any(o[[2]][,2]==0),o[[1]][min(which(o[[2]][,2]==0))],50)) %>%
                  mutate(., ## compute the dispersion parameter for each run that did not fade out
                         disp=ifelse(fadeout==0, 
                                     ifelse(inherits(try(glm.nb((o[[3]][is.na(o[[3]][,4]),3])~1), silent=TRUE),"try-error"), 
                                            NA, 
                                            glm.nb((o[[3]][is.na(o[[3]][,4]),3])~1)$theta), 
                                     NA))
                ) %>% do.call("rbind.data.frame",.)
       ) -> OutC
         

 
  
  
  var1=c('c','c','c','shed','shed','alpha')
  var2=c('shed','alpha','gamma','alpha','gamma','gamma')
  
  for (j in 1:6) { ## loop over the six different covariance combinations
    for (corr in c("nocorr", "negcorr", "poscorr")) {
      for (varLevel in c("hivar","medvar","lowvar")) {
        print(paste("out",corr,varLevel,var1[j],var2[j],sep="_"))
        mclapply(seeds, 
                 function(i) gillespie.SIR.cov_storage(tmax=100, 
                                                       params=get(varLevel), 
                                                       corr=get(corr), 
                                                       x=initial_state, 
                                                       covParams=c(var1[j],var2[j]),
                                                       seed=i),
                 mc.cores=10) -> z
        saveRDS(z, file=paste0(paste(paste0("out_R=",R),corr,varLevel,var1[j],var2[j],sep="_"),".RDS"))
      }
    }
  }
}

secInf = (o[[3]][is.na(o[[3]][,3]),2])
fit = fitdistr((o[[3]][is.na(o[[3]][,3]),2]), "negative binomial")
hist(secInf, prob=TRUE, main="")
lines(seq(0,15), dnbinom(seq(0,15), size = fit$estimate[1], mu = fit$estimate[2]), lwd=2, col=2)
