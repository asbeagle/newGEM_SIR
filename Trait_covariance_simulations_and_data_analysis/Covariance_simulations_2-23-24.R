library(Rcpp)
library(RcppEigen)
library(MASS)
library(tidyverse)
library(ggplot2)
library(magrittr)
library(parallel)
sourceCpp("SIRcov2.cpp")

var1=c('c','c','c','s','s','a')
var2=c('s','a','g','a','g','g')

for (R in c(1,2,4,8)) {
  for (S in c(100, 250, 500, 1000)) { ## vary iniital number of susceptibles
    for (CV in c(.2, 1, 5)) {
      for (corr in seq(-0.9, 0.9, 0.1)) {
        for (parcomb in seq(1,6)) {
          print(R)
          print(S)
          print(CV)
          print(corr)
          print(var1[parcomb])
          print(var2[parcomb])
          corr = round(corr,1) ## force to be precisely the correct value
          ## baseline parameters
          params = c(c=2*R/S, s=2*R/(S-2*R), a=2*R/S, g=2*R/S,
                     SD1=0, SD2=0, corr=corr)
          params["SD1"] = unname(CV*params[var1[parcomb]])
          params["SD2"] = unname(CV*params[var2[parcomb]])
        
          ## initial conditions
          initial_state <- c(S=S, I=S/50, R=0)
          ## simulate!
          if (parcomb==1)
            mclapply(1:1000,
                     function(i) SIRcovCS(params, initial_state),
                     mc.cores=10) -> out
          else if (parcomb==2)
            mclapply(1:1000,
                     function(i) SIRcovCA(params, initial_state),
                     mc.cores=10) -> out
          else if (parcomb==3)
            mclapply(1:1000,
                     function(i) SIRcovCG(params, initial_state),
                     mc.cores=10) -> out
          else if (parcomb==4)
            mclapply(1:1000,
                     function(i) SIRcovSA(params, initial_state),
                     mc.cores=10) -> out
          else if (parcomb==5)
            mclapply(1:1000,
                     function(i) SIRcovSG(params, initial_state),
                     mc.cores=10) -> out
          else 
            mclapply(1:1000,
                     function(i) SIRcovAG(params, initial_state),
                     mc.cores=10) -> out
        
          ## Because out is very, very large, it must be analyzed before saving
          lapply(out, 
                 function(o)
                   data.frame(peak=max(o[[2]][,2]),
                              totalI=sum(o[[3]][,4]),
                              tEnd=tail(o[[1]],1),
                              disp=ifelse(inherits(try(glm.nb((o[[3]][is.na(o[[3]][,5]),4])~1), silent=TRUE),"try-error"), 
                                          NA, 
                                          glm.nb((o[[3]][is.na(o[[3]][,5]),4])~1)$theta))
          ) %>% 
            do.call("rbind.data.frame",.) %>%
            mutate(rep=1:1000, 
                   R0=R,
                   S=S,
                   CV=CV,
                   corr=corr, 
                   cov=paste0(var1[parcomb],"-",var2[parcomb])) -> outSummary
          
        
        saveRDS(outSummary, file=paste0("out_R0=",R,"_S=",S,"_CV=",CV,"_corr=",corr,"_cov=",var1[parcomb],"-",var2[parcomb],".RDS"))
        rm(out) ## prevents accidentally saving the same thing twice and saves some storage space
        rm(outSummary)
        }
      }
    }
  }
}


# 
# 
# 
# ## baseline parameters to produce R0 assuming an initial S size of 1000
# ## c*s/(1+s)*S/(a+g) = 2; we want all of the parameters to be equal so that none has an undue influence
# ## on the results
# R = 2
# params = c(c=R/500, s=R/(500-R), a=R/500, g=R/500,
#            SD1=0, SD2=0, corr=0)
# ## Set the SDs based on the CV
# CV = 1
# params["SD1"] = unname(CV*params["c"])
# params["SD2"] = unname(CV*params["s"])
# ## Set the initial conditions
# initial_state = c(S=1000, I=10, R=0)
# 
# out = SIRcovCS(params, initial_state)
# glm.nb((out[[3]][is.na(out[[3]][,5]),4])~1)$theta
# 
# params["SD1"] = unname(CV*params["c"])
# params["SD2"] = unname(CV*params["a"])
# out = SIRcovCA(params, initial_state)
# tail(out[[2]])
# glm.nb((out[[3]][is.na(out[[3]][,5]),4])~1)$theta
# 
# params["SD1"] = unname(CV*params["c"])
# params["SD2"] = unname(CV*params["g"])
# out = SIRcovCG(params, initial_state)
# tail(out[[2]])
# glm.nb((out[[3]][is.na(out[[3]][,5]),4])~1)$theta
# 
# params["SD1"] = unname(CV*params["s"])
# params["SD2"] = unname(CV*params["a"])
# out = SIRcovSA(params, initial_state)
# tail(out[[2]])
# glm.nb((out[[3]][is.na(out[[3]][,5]),4])~1)$theta
# 
# params["SD1"] = unname(CV*params["s"])
# params["SD2"] = unname(CV*params["g"])
# out = SIRcovSG(params, initial_state)
# tail(out[[2]])
# glm.nb((out[[3]][is.na(out[[3]][,5]),4])~1)$theta
# 
# params["SD1"] = unname(CV*params["a"])
# params["SD2"] = unname(CV*params["g"])
# out = SIRcovAG(params, initial_state)
# tail(out[[2]])
# glm.nb((out[[3]][is.na(out[[3]][,5]),4])~1)$theta
# 
# 
# params2 = c(c=R/500, s=R/(500-R), a=R/500, g=R/500, SD=0)
# params2["SD"] = unname(CV*params["g"])
# out = SIRvarG(params2, initial_state)
# tail(out[[2]])
# glm.nb((out[[3]][is.na(out[[3]][,4]),3])~1)$theta
# 
# 
# params2 = c(c=R/500, s=R/(500-R), a=R/500, g=R/500, SD=0)
# params2["SD"] = unname(CV*params["a"])
# out = SIRvarA(params2, initial_state)
# tail(out[[2]])
# glm.nb((out[[3]][is.na(out[[3]][,4]),3])~1)$theta
# 
# params2 = c(c=R/500, s=R/(500-R), a=R/500, g=R/500, SD=0)
# params2["SD"] = unname(CV*params["s"])
# out = SIRvarS(params2, initial_state)
# tail(out[[2]])
# glm.nb((out[[3]][is.na(out[[3]][,4]),3])~1)$theta
# 
# params2 = c(c=R/500, s=R/(500-R), a=R/500, g=R/500, SD=0)
# params2["SD"] = unname(CV*params["c"])
# out = SIRvarC(params2, initial_state)
# tail(out[[2]])
# glm.nb((out[[3]][is.na(out[[3]][,4]),3])~1)$theta
# 
# params3 = c(c=R/500, s=R/(500-R), a=R/500, g=R/500)
# out = SIRnovar(params3, initial_state)
# tail(out[[2]])
# glm.nb((out[[3]][is.na(out[[3]][,3]),2])~1)$theta
# 
# 
# 
