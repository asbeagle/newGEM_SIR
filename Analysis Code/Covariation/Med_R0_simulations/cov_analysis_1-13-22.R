source("GEM_SIR_cov_storage.R")
library(tidyverse)
library(parallel)

## set up
set.seed(123432)
seeds <- floor(runif(50,1,1e5)) # set seeds

nocorr <- matrix(c(1,0,0,1), nrow=2, byrow=T)
negcorr <- matrix(c(1,-.5,-.5,1), nrow=2, byrow=T)
poscorr <- matrix(c(1,.5,.5,1), nrow=2, byrow=T)

hivar = c(   c=0.1,    shed=0.05,    alpha=0.1,    gamma=0.1, 
          sd_c=0.5, sd_shed=0.25, sd_alpha=0.5, sd_gamma=0.5, 
             b=2.5, d=.1, bs=.01)
medvar = c(   c=0.1,    shed=0.05,    alpha=0.1,    gamma=0.1, 
           sd_c=0.1, sd_shed=0.05, sd_alpha=0.1, sd_gamma=0.1, 
              b=2.5, d=.1, bs=.01)
lowvar = c(   c=0.1,     shed=0.05,    alpha=0.1,     gamma=0.1, 
           sd_c=0.02, sd_shed=0.01, sd_alpha=0.02, sd_gamma=0.02, 
              b=2.5, d=.1, bs=.01)

initial_state <- floor(c(S=unname(((hivar["b"]-hivar["d"])/hivar["bs"]))-5, I=5, R=0))

var1=c('c','c','c','shed','shed','alpha')
var2=c('shed','alpha','gamma','alpha','gamma','gamma')

for (j in 1:6) { ## loop over the six different covariance combinations
  for (covMatrix in c("nocorr", "negcorr", "poscorr")) {
    for (varLevel in c("hivar","medvar","lowvar")) {
      print(paste("out",covMatrix,varLevel,var1[j],var2[j],sep="_"))
      mclapply(seeds, 
               function(i) gillespie.SIR.cov_storage(tmax=100, 
                                                     params=get(varLevel), 
                                                     corr=get(covMatrix), 
                                                     x=initial_state, 
                                                     covParams=c(var1[j],var2[j]),
                                                     seed=i),
               mc.cores=10) -> z
      saveRDS(z, file=paste0(paste("out",covMatrix,varLevel,var1[j],var2[j],sep="_"),".RDS"))
    }
  }
}
  

var1=c('c','c','c','shed','shed','alpha')
var2=c('shed','alpha','gamma','alpha','gamma','gamma')

data <- vector(mode='list', length=54)
i <- 1
for (j in 1:6) { ## loop over the six different covariance combinations
  for (covMatrix in c("nocorr", "negcorr", "poscorr")) {
    for (varLevel in c("hivar","medvar","lowvar")) {
      out <- readRDS(paste0(paste("out",covMatrix,varLevel,var1[j],var2[j],sep="_"),".RDS"))
      data.frame(traits=rep(paste(var1[j],var2[j],sep='-'),50),
                 cov=rep(covMatrix,50),
                 var=rep(varLevel,50),
                 top20=unlist(lapply(out, function(o) sum(sort(o[[2]]$numInf)[round(0.8*nrow(o[[2]])):nrow(o[[2]])])/sum(o[[2]]$numInf))),
                 meanI=unlist(lapply(out, function(o) mean(tail(o[[1]]$I,20)))),
                 maxI=unlist(lapply(out, function(o) max(o[[1]]$I))),
                 medNumInfTop20=unlist(lapply(out, function(o) median(sort(o[[2]]$numInf)[round(0.8*nrow(o[[2]])):nrow(o[[2]])]))),
                 meanNumInfTop20=unlist(lapply(out, function(o) mean(sort(o[[2]]$numInf)[round(0.8*nrow(o[[2]])):nrow(o[[2]])]))),
                 sdNumInfTop20=unlist(lapply(out, function(o) sd(sort(o[[2]]$numInf)[round(0.8*nrow(o[[2]])):nrow(o[[2]])]))),
                 maxNumInf=unlist(lapply(out, function(o) log(max(sort(o[[2]]$numInf)[round(0.8*nrow(o[[2]])):nrow(o[[2]])])))),
                 peakTime=unlist(lapply(out, function(o) which.max(o[[1]]$I)-1))
                 ) -> data[[i]]
      i <- i + 1
    }
  }
}
data <- do.call("rbind.data.frame",data)

## Epidemiological dynamics in general
## This figure is interesting - shows that *having* traits covary affects the peak size of the epidemic but the *direction* of the covariance doesn't actually matter
ggplot(subset(data, meanI>10), aes(x=var,y=maxI,shape=cov,colour=cov)) + 
  geom_boxplot() + 
  facet_wrap(~traits)


## Super-spreading figures 
## This one shows the fraction of total cases caused by the top 20% of spreaders
ggplot(subset(data, meanI>10), aes(x=var,y=top20,shape=cov,colour=cov)) + 
  geom_boxplot() + 
  facet_wrap(~traits)
## This one shows the mean number of infections caused by the individuals in the top 20% of spreaders
ggplot(subset(data, meanI>10), aes(x=var,y=meanNumInfTop20,shape=cov,colour=cov)) + 
  geom_boxplot() + 
  facet_wrap(~traits)


## Effect of super-spreading on epidemiological dynamics
ggplot(subset(data, meanI>10), aes(x=meanNumInfTop20,y=maxI,group=var,shape=var)) + 
  geom_point(aes(colour=cov)) + 
  facet_wrap(~traits)

ggplot(subset(data, meanI>10), aes(x=top20,y=meanI,group=var,shape=var)) + 
  geom_point(aes(colour=cov)) + 
  facet_wrap(~traits)



# ## Plot median trajectories of S (expected equilibrium with no variation is S=63)
# png(filename="~/Desktop/Covariance_results.png", height=10, width=5, units='in', res=600)
# par(mfrow=c(6,3), mar=c(1,1,1,1), oma=c(3,3,0,0))
# for (j in 1:6) { ## loop over the six different covariance combinations
#   for (covMatrix in c("nocorr", "negcorr", "poscorr")) {
#       plot.new()
#       plot.window(xlim=c(0,100), ylim=c(0,235))
#       axis(1); axis(2); box('plot')
#       z1 = get(paste("out",covMatrix,"hivar",var1[j],var2[j],sep="_"))
#       z2 = get(paste("out",covMatrix,"medvar",var1[j],var2[j],sep="_"))
#       z3 = get(paste("out",covMatrix,"lowvar",var1[j],var2[j],sep="_"))
#       abline(h=63, lty=2)
#       lines(0:100, lapply(z1, function(l) l[[1]]$S) %>% do.call("cbind.data.frame",.) %>% apply(., 1, median), lwd=2)
#       lines(0:100, lapply(z2, function(l) l[[1]]$S) %>% do.call("cbind.data.frame",.) %>% apply(., 1, median), lwd=2, col=4)
#       lines(0:100, lapply(z3, function(l) l[[1]]$S) %>% do.call("cbind.data.frame",.) %>% apply(., 1, median), lwd=2, col=2)
#       legend(x='topright', paste(var1[j], var2[j], covMatrix, sep='-'), bty='n')
#   }
# }
# dev.off()
# 
#     
      


