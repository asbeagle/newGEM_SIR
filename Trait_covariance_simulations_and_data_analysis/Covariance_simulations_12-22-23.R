library(Rcpp)
library(RcppEigen)
library(deSolve) 
library(MASS)
library(tidyverse)
library(ggplot2)
library(magrittr)
library(parallel)
sourceCpp("SIRcov.cpp")

var1=c('c','c','c','s','s','a')
var2=c('s','a','g','a','g','g')

for (R in c(1,4,8)) {
  for (CV in c(.2, 1, 5)) {
    for (corr in seq(-0.9, 0.9, 0.1)) {
      for (parcomb in seq(1,6)) {
        print(CV)
        print(R)
        print(corr)
        print(var1[parcomb])
        print(var2[parcomb])
        corr = round(corr,1) ## force to be precisely the correct value
        ## baseline parameters
        b = 2.5; bs=0.01; d=0.1; c=0.1; s=1/9; a=0.1; g=0.1
        ## scaling coefficient x to scale c, s, a, and g to achieve particular R value
        x = 1/(2*(-b*c*s+c*d*s-a*bs*R*s-bs*g*R*s)) * (-a*bs*R-bs*g*R-2*b*c*s+2*c*d*s+bs*d*R*s+sqrt((a*bs*R+bs*g*R+2*b*c*s-2*c*d*s-bs*d*R*s)^2-4*(-b*c*s+c*d*s-a*bs*R*s-bs*g*R*s)*(a*bs*R+bs*d*R+bs*g*R-b*c*s+c*d*s+a*bs*R*s+bs*d*R*s+bs*g*R*s)))
        ## simulation parameters depends on the parameter combination
        params = c(b=2.5, bs=0.01, d=0.1, 
                   c=(1-x)*0.1, s=(1-x)*1/9, a=(1+x)*0.1, g=(1+x)*0.1,
                   SD1=0, SD2=0, corr=corr)
        params["SD1"] = unname(CV*params[var1[parcomb]])
        params["SD2"] = unname(CV*params[var2[parcomb]])
        
        ## initial conditions
        initial_state <- floor(c(S=unname(((params["b"]-params["d"])/params["bs"]))-5, I=5, R=0))
        ## simulate!
        if (parcomb==1)
          mclapply(1:1000,
                   function(i) SIRcovCS(params, initial_state, 50),
                   mc.cores=10) -> out
        else if (parcomb==2)
          mclapply(1:1000,
                   function(i) SIRcovCA(params, initial_state, 50),
                   mc.cores=10) -> out
        else if (parcomb==3)
          mclapply(1:1000,
                   function(i) SIRcovCG(params, initial_state, 50),
                   mc.cores=10) -> out
        else if (parcomb==4)
          mclapply(1:1000,
                   function(i) SIRcovSA(params, initial_state, 50),
                   mc.cores=10) -> out
        else if (parcomb==5)
          mclapply(1:1000,
                   function(i) SIRcovSG(params, initial_state, 50),
                   mc.cores=10) -> out
        else 
          mclapply(1:1000,
                   function(i) SIRcovAG(params, initial_state, 50),
                   mc.cores=10) -> out
        
        
        ## Because out is very, very large, it must be analyzed before saving
        lapply(out, 
               function(o)
                 data.frame(peak=max(o[[2]][,2]),
                            totalI=sum(o[[3]][,4]),
                            fadeout=ifelse(any(o[[2]][,2]==0),1,0),
                            fadeoutT=ifelse(any(o[[2]][,2]==0),o[[1]][min(which(o[[2]][,2]==0))],50)) %>%
                 mutate(., ## compute the dispersion parameter for each run that did not fade out
                        disp=ifelse(fadeout==0, 
                                    ifelse(inherits(try(glm.nb((o[[3]][is.na(o[[3]][,5]),4])~1), silent=TRUE),"try-error"), 
                                           NA, 
                                           glm.nb((o[[3]][is.na(o[[3]][,5]),4])~1)$theta), 
                                    NA))
        ) %>% 
          do.call("rbind.data.frame",.) %>%
          mutate(rep=1:1000, 
                 R0=R,
                 CV=CV,
                 corr=corr, 
                 cov=paste0(var1[parcomb],"-",var2[parcomb])) -> outSummary
        
            
        saveRDS(outSummary, file=paste0("out_R0=",R,"_CV=",CV,"_corr=",corr,"_cov=",var1[parcomb],"-",var2[parcomb],".RDS"))
        rm(out) ## prevents accidentally saving the same thing twice and saves some storage space
        rm(outSummary)
      }
    }
  }
}
