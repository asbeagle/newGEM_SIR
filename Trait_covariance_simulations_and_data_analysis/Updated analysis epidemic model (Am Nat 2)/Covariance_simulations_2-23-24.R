library(Rcpp)
library(RcppEigen)
library(MASS)
library(tidyverse)
library(ggplot2)
library(magrittr)
library(parallel)
sourceCpp("SIRcov2.cpp")

# var1=c('c','c','c','s','s','a')
# var2=c('s','a','g','a','g','g')
# 
# for (R in c(1,2,4,8)) {
#   for (S in c(100, 250, 500, 1000)) { ## vary initial number of susceptibles
#     for (CV in c(.2, 1, 5)) {
#       for (corr in seq(-0.9, 0.9, 0.1)) {
#         for (parcomb in seq(1,6)) {
#           print(R)
#           print(S)
#           print(CV)
#           print(corr)
#           print(var1[parcomb])
#           print(var2[parcomb])
#           corr = round(corr,1) ## force to be precisely the correct value
#           ## baseline parameters
#           params = c(c=2*R/S, s=2*R/(S-2*R), a=2*R/S, g=2*R/S,
#                      SD1=0, SD2=0, corr=corr)
#           params["SD1"] = unname(CV*params[var1[parcomb]])
#           params["SD2"] = unname(CV*params[var2[parcomb]])
#         
#           ## initial conditions
#           initial_state <- c(S=S, I=S/50, R=0)
#           ## simulate!
#           if (parcomb==1)
#             mclapply(1:1000,
#                      function(i) SIRcovCS(params, initial_state),
#                      mc.cores=10) -> out
#           else if (parcomb==2)
#             mclapply(1:1000,
#                      function(i) SIRcovCA(params, initial_state),
#                      mc.cores=10) -> out
#           else if (parcomb==3)
#             mclapply(1:1000,
#                      function(i) SIRcovCG(params, initial_state),
#                      mc.cores=10) -> out
#           else if (parcomb==4)
#             mclapply(1:1000,
#                      function(i) SIRcovSA(params, initial_state),
#                      mc.cores=10) -> out
#           else if (parcomb==5)
#             mclapply(1:1000,
#                      function(i) SIRcovSG(params, initial_state),
#                      mc.cores=10) -> out
#           else 
#             mclapply(1:1000,
#                      function(i) SIRcovAG(params, initial_state),
#                      mc.cores=10) -> out
#         
#           ## Because out is very, very large, it must be analyzed before saving
#           lapply(out, 
#                  function(o)
#                    data.frame(peak=max(o[[2]][,2]),
#                               totalI=sum(o[[3]][,4]),
#                               tEnd=tail(o[[1]],1),
#                               disp=ifelse(inherits(try(glm.nb((o[[3]][is.na(o[[3]][,5]),4])~1), silent=TRUE),"try-error"), 
#                                           NA, 
#                                           glm.nb((o[[3]][is.na(o[[3]][,5]),4])~1)$theta))
#           ) %>% 
#             do.call("rbind.data.frame",.) %>%
#             mutate(rep=1:1000, 
#                    R0=R,
#                    S=S,
#                    CV=CV,
#                    corr=corr, 
#                    cov=paste0(var1[parcomb],"-",var2[parcomb])) -> outSummary
#           
#         
#         saveRDS(outSummary, file=paste0("out_R0=",R,"_S=",S,"_CV=",CV,"_corr=",corr,"_cov=",var1[parcomb],"-",var2[parcomb],".RDS"))
#         rm(out) ## prevents accidentally saving the same thing twice and saves some storage space
#         rm(outSummary)
#         }
#       }
#     }
#   }
# }
# 

# for (R in c(1,2,4,8)) {
#   for (S in c(100, 250, 500, 1000)) { ## vary initial number of susceptibles
#     print(R)
#     print(S)
#     initial_state <- c(S=S, I=S/50, R=0)
#     params = c(c=2*R/S, s=2*R/(S-2*R), a=2*R/S, g=2*R/S)
#     mclapply(1:1000,
#              function(i) SIRnovar(params, initial_state),
#              mc.cores=10) -> out
# 
#     lapply(out,
#            function(o)
#              data.frame(peak=max(o[[2]][,2]),
#                         totalI=sum(o[[3]][,2]),
#                         tEnd=tail(o[[1]],1),
#                         disp=ifelse(inherits(try(glm.nb((o[[3]][is.na(o[[3]][,3]),2])~1), silent=TRUE),"try-error"),
#                                     NA,
#                                     glm.nb((o[[3]][is.na(o[[3]][,3]),2])~1)$theta))
#     ) %>%
#       do.call("rbind.data.frame",.) %>%
#       mutate(rep=1:1000,
#              R0=R,
#              S=S) -> outSummary
# 
#     saveRDS(outSummary, file=paste0("out_R0=",R,"_S=",S,"_no_variation.RDS"))
#     rm(out) ## prevents accidentally saving the same thing twice and saves some storage space
#     rm(outSummary)
#   }
# }
# 
# for (R in c(1,2,4,8)) {
#   for (S in c(100, 250, 500, 1000)) { ## vary initial number of susceptibles
#     for (CV in c(.2, 1, 5)) {
#       print(R)
#       print(S)
#       print(CV)
#       
#       initial_state <- c(S=S, I=S/50, R=0)
#       paramsC = c(c=2*R/S, s=2*R/(S-2*R), a=2*R/S, g=2*R/S, SD=CV*2*R/S)
#       paramsS = c(c=2*R/S, s=2*R/(S-2*R), a=2*R/S, g=2*R/S, SD=CV*2*R/(S-2*R))
#       paramsA = c(c=2*R/S, s=2*R/(S-2*R), a=2*R/S, g=2*R/S, SD=CV*2*R/S)
#       paramsG = c(c=2*R/S, s=2*R/(S-2*R), a=2*R/S, g=2*R/S, SD=CV*2*R/S)
# 
#       mclapply(1:1000,
#                function(i) SIRvarC(paramsC, initial_state),
#                mc.cores=10) -> outC
#       
#       mclapply(1:1000,
#                function(i) SIRvarS(paramsS, initial_state),
#                mc.cores=10) -> outS
#       
#       mclapply(1:1000,
#                function(i) SIRvarA(paramsA, initial_state),
#                mc.cores=10) -> outA
#       
#       mclapply(1:1000,
#                function(i) SIRvarG(paramsG, initial_state),
#                mc.cores=10) -> outG
#       
#       lapply(outC, 
#              function(o)
#                data.frame(peak=max(o[[2]][,2]),
#                           totalI=sum(o[[3]][,3]),
#                           tEnd=tail(o[[1]],1),
#                           disp=ifelse(inherits(try(glm.nb((o[[3]][is.na(o[[3]][,4]),3])~1), silent=TRUE),"try-error"), 
#                                       NA, 
#                                       glm.nb((o[[3]][is.na(o[[3]][,4]),3])~1)$theta))
#       ) %>% 
#         do.call("rbind.data.frame",.) %>%
#         mutate(rep=1:1000, 
#                R0=R,
#                S=S,
#                CV=CV,
#                var='c') -> outSummaryC
#       saveRDS(outSummaryC, file=paste0("out_R0=",R,"_S=",S,"_CV=",CV,"_var=c.RDS"))
#       
#       lapply(outS, 
#              function(o)
#                data.frame(peak=max(o[[2]][,2]),
#                           totalI=sum(o[[3]][,3]),
#                           tEnd=tail(o[[1]],1),
#                           disp=ifelse(inherits(try(glm.nb((o[[3]][is.na(o[[3]][,4]),3])~1), silent=TRUE),"try-error"), 
#                                       NA, 
#                                       glm.nb((o[[3]][is.na(o[[3]][,4]),3])~1)$theta))
#       ) %>% 
#         do.call("rbind.data.frame",.) %>%
#         mutate(rep=1:1000, 
#                R0=R,
#                S=S,
#                CV=CV,
#                var='s') -> outSummaryS
#       saveRDS(outSummaryS, file=paste0("out_R0=",R,"_S=",S,"_CV=",CV,"_var=s.RDS"))
#       
#       lapply(outA, 
#              function(o)
#                data.frame(peak=max(o[[2]][,2]),
#                           totalI=sum(o[[3]][,3]),
#                           tEnd=tail(o[[1]],1),
#                           disp=ifelse(inherits(try(glm.nb((o[[3]][is.na(o[[3]][,4]),3])~1), silent=TRUE),"try-error"), 
#                                       NA, 
#                                       glm.nb((o[[3]][is.na(o[[3]][,4]),3])~1)$theta))
#       ) %>% 
#         do.call("rbind.data.frame",.) %>%
#         mutate(rep=1:1000, 
#                R0=R,
#                S=S,
#                CV=CV,
#                var='a') -> outSummaryA
#       saveRDS(outSummaryA, file=paste0("out_R0=",R,"_S=",S,"_CV=",CV,"_var=a.RDS"))
#       
#       lapply(outG, 
#              function(o)
#                data.frame(peak=max(o[[2]][,2]),
#                           totalI=sum(o[[3]][,3]),
#                           tEnd=tail(o[[1]],1),
#                           disp=ifelse(inherits(try(glm.nb((o[[3]][is.na(o[[3]][,4]),3])~1), silent=TRUE),"try-error"), 
#                                       NA, 
#                                       glm.nb((o[[3]][is.na(o[[3]][,4]),3])~1)$theta))
#       ) %>% 
#         do.call("rbind.data.frame",.) %>%
#         mutate(rep=1:1000, 
#                R0=R,
#                S=S,
#                CV=CV,
#                var='g') -> outSummaryG
#       saveRDS(outSummaryG, file=paste0("out_R0=",R,"_S=",S,"_CV=",CV,"_var=g.RDS"))
#       
#       rm(outC) ## prevents accidentally saving the same thing twice and saves some storage space
#       rm(outSummaryC)
#       rm(outS) ## prevents accidentally saving the same thing twice and saves some storage space
#       rm(outSummaryS)
#       rm(outA) ## prevents accidentally saving the same thing twice and saves some storage space
#       rm(outSummaryA)
#       rm(outG) ## prevents accidentally saving the same thing twice and saves some storage space
#       rm(outSummaryG)
#     }
#   }
# }
#       
    


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
var1=c('c','c','c','s','s','a')
var2=c('s','a','g','a','g','g')

for (R in c(1,2,4,8)) {
    for (CV in c(.2, 1, 5)) {
      for (corr in seq(-0.9, 0.9, 0.1)) {
        for (parcomb in seq(1,6)) {
          S = 1000
          print(R)
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
          initial_state <- c(S=S, I=S/100, R=0)
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
                                          glm.nb((o[[3]][is.na(o[[3]][,5]),4])~1)$theta),
                              top20=sum(sort(o[[3]][,4], decreasing=T)[1:floor(0.2*nrow(o[[3]]))])/sum(o[[3]][,4]))
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

for (R in c(1,2,4,8)) {
  for (CV in c(.2, 1, 5)) {
    S = 1000
    print(R)
    print(S)
    print(CV)
    
    initial_state <- c(S=S, I=S/100, R=0)
    paramsC = c(c=2*R/S, s=2*R/(S-2*R), a=2*R/S, g=2*R/S, SD=CV*2*R/S)
    paramsS = c(c=2*R/S, s=2*R/(S-2*R), a=2*R/S, g=2*R/S, SD=CV*2*R/(S-2*R))
    paramsA = c(c=2*R/S, s=2*R/(S-2*R), a=2*R/S, g=2*R/S, SD=CV*2*R/S)
    paramsG = c(c=2*R/S, s=2*R/(S-2*R), a=2*R/S, g=2*R/S, SD=CV*2*R/S)
    
    mclapply(1:1000,
             function(i) SIRvarC(paramsC, initial_state),
             mc.cores=10) -> outC
    
    mclapply(1:1000,
             function(i) SIRvarS(paramsS, initial_state),
             mc.cores=10) -> outS
    
    mclapply(1:1000,
             function(i) SIRvarA(paramsA, initial_state),
             mc.cores=10) -> outA
    
    mclapply(1:1000,
             function(i) SIRvarG(paramsG, initial_state),
             mc.cores=10) -> outG
    
    lapply(outC,
           function(o)
             data.frame(peak=max(o[[2]][,2]),
                        totalI=sum(o[[3]][,3]),
                        tEnd=tail(o[[1]],1),
                        disp=ifelse(inherits(try(glm.nb((o[[3]][is.na(o[[3]][,4]),3])~1), silent=TRUE),"try-error"),
                                    NA,
                                    glm.nb((o[[3]][is.na(o[[3]][,4]),3])~1)$theta),
                        top20=sum(sort(o[[3]][,3], decreasing=T)[1:floor(0.2*nrow(o[[3]]))])/sum(o[[3]][,3]))
    ) %>%
      do.call("rbind.data.frame",.) %>%
      mutate(rep=1:1000,
             R0=R,
             S=S,
             CV=CV,
             var='c') -> outSummaryC
    saveRDS(outSummaryC, file=paste0("out_R0=",R,"_S=",S,"_CV=",CV,"_var=c.RDS"))
    
    lapply(outS,
           function(o)
             data.frame(peak=max(o[[2]][,2]),
                        totalI=sum(o[[3]][,3]),
                        tEnd=tail(o[[1]],1),
                        disp=ifelse(inherits(try(glm.nb((o[[3]][is.na(o[[3]][,4]),3])~1), silent=TRUE),"try-error"),
                                    NA,
                                    glm.nb((o[[3]][is.na(o[[3]][,4]),3])~1)$theta),
                        top20=sum(sort(o[[3]][,3], decreasing=T)[1:floor(0.2*nrow(o[[3]]))])/sum(o[[3]][,3]))
    ) %>%
      do.call("rbind.data.frame",.) %>%
      mutate(rep=1:1000,
             R0=R,
             S=S,
             CV=CV,
             var='s') -> outSummaryS
    saveRDS(outSummaryS, file=paste0("out_R0=",R,"_S=",S,"_CV=",CV,"_var=s.RDS"))
    
    lapply(outA,
           function(o)
             data.frame(peak=max(o[[2]][,2]),
                        totalI=sum(o[[3]][,3]),
                        tEnd=tail(o[[1]],1),
                        disp=ifelse(inherits(try(glm.nb((o[[3]][is.na(o[[3]][,4]),3])~1), silent=TRUE),"try-error"),
                                    NA,
                                    glm.nb((o[[3]][is.na(o[[3]][,4]),3])~1)$theta),
                        top20=sum(sort(o[[3]][,3], decreasing=T)[1:floor(0.2*nrow(o[[3]]))])/sum(o[[3]][,3]))
    ) %>%
      do.call("rbind.data.frame",.) %>%
      mutate(rep=1:1000,
             R0=R,
             S=S,
             CV=CV,
             var='a') -> outSummaryA
    saveRDS(outSummaryA, file=paste0("out_R0=",R,"_S=",S,"_CV=",CV,"_var=a.RDS"))
    
    lapply(outG,
           function(o)
             data.frame(peak=max(o[[2]][,2]),
                        totalI=sum(o[[3]][,3]),
                        tEnd=tail(o[[1]],1),
                        disp=ifelse(inherits(try(glm.nb((o[[3]][is.na(o[[3]][,4]),3])~1), silent=TRUE),"try-error"),
                                    NA,
                                    glm.nb((o[[3]][is.na(o[[3]][,4]),3])~1)$theta),
                        top20=sum(sort(o[[3]][,3], decreasing=T)[1:floor(0.2*nrow(o[[3]]))])/sum(o[[3]][,3]))
    ) %>%
      do.call("rbind.data.frame",.) %>%
      mutate(rep=1:1000,
             R0=R,
             S=S,
             CV=CV,
             var='g') -> outSummaryG
    saveRDS(outSummaryG, file=paste0("out_R0=",R,"_S=",S,"_CV=",CV,"_var=g.RDS"))
    
    rm(outC) ## prevents accidentally saving the same thing twice and saves some storage space
    rm(outSummaryC)
    rm(outS) ## prevents accidentally saving the same thing twice and saves some storage space
    rm(outSummaryS)
    rm(outA) ## prevents accidentally saving the same thing twice and saves some storage space
    rm(outSummaryA)
    rm(outG) ## prevents accidentally saving the same thing twice and saves some storage space
    rm(outSummaryG)
  }
}

## Initial number of susceptible hosts
S = 1000
## Expected R0
R = 4
## Magnitude of variation
CV = 1
## Trait correlation
corr = -0.5
## Covarying traits
var1 = "s"
var2 = "a"
## Set parameters and initial conditions
params = c(c=2*R/S, s=2*R/(S-2*R), a=2*R/S, g=2*R/S,
           SD1=CV*2*R/(S-2*R), SD2=CV*2*R/S, corr=corr)
initial_state <- c(S=S, I=S/100, R=0)
## Simulate the stochastic epidemic model
out = SIRcovSA(params, initial_state)


prev = out[[2]] %>% apply(., 1, function(x) x[2]/sum(x))

png("Conceptual_figure_panel_1.png", height=5, width=5, units='in', res=450)
plot(out[[1]], prev, type='l', xlab="Time", ylab="Infection prevalence")
arrows(x0=(out[[1]][which(prev==max(prev))]+100), 
       x1=out[[1]][which(prev==max(prev))], 
       y0=max(prev), 
       y1=max(prev), 
       length=0.1, col="red")
text(x=(out[[1]][which(prev==max(prev))]+110), 
     y=max(prev), 
     "Peak prevalence = 0.55", 
     adj=0, col="red")
arrows(x0=max(out[[1]]), 
       x1=max(out[[1]]), 
       y0=0.08, 
       y1=-0.02, 
       length=0.1, col="red")
text(x=max(out[[1]]), 
     y=0.11, 
     "Epidemic\nduration = 813", 
     adj=1, col="red")
dev.off()

## Fit negative binomial distribution to observed no of secondary case counts
count = out[[3]][,4]
fit <- glm.nb(count ~ 1)
mu_val <- unname(exp(coef(fit)[1]))
theta_val <- fit$theta
n_obs <- length(count)
bw <- 1 

x_vals <- seq(min(count), max(count), by = 1)
fitted_df <- data.frame(
  count = x_vals,
  pred = n_obs * bw * dnbinom(x_vals, size = theta_val, mu = mu_val)
)

df = data.frame(count=count)
png(file="Conceptual_figure_panel_2.png", height=5, width=5, units='in', res=450)
plot.new()
plot.window(xlim=c(-0.5,19.5), ylim=c(0,600))
axis(1)
axis(2)
box('plot')
rect(xleft=seq(-0.5,18.5,1), 
     ybottom=rep(0,20), 
     xright=seq(0.5,19.5,1), 
     ytop=sapply(seq(0,max(count)), function(n) sum(count==n)),
     col="gray75")
mtext(side=1, "No. of secondary cases caused", line=3)
mtext(side=2, "Frequency", line=3)
lines(fitted_df$count, fitted_df$pred, col="red", lwd=1.5)
points(fitted_df$count, fitted_df$pred, col="red", pch=21, bg="red")
text(x=2, y=590, "Dispersion estimate = 0.45", col="red", adj=0)
dev.off()

png(file="Conceptual_figure_panel_3.png", height=5, width=5, units='in', res=450)
plot(seq(1,length(count))/length(count), cumsum(sort(count, decreasing=TRUE))/max(cumsum(count)), 
     type='l', xlab="Infectiousness (ranked percentile of cases)", ylab="Proportion of total transmission")
lines(x=c(0.2,0.2), y=c(0,0.771), col=2)
arrows(x0=0.2, x1=-0.035, y0=0.771, y1=0.771, length=0.1, col="red")
text(0.22, 0.771, "Proportion of transmission from\ntop 20% most infectious cases = 0.77", col="red", adj=0)
dev.off()
