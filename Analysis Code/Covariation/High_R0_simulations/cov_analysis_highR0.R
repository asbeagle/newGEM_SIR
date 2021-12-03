source("GEM_SIR_cov_storage.R")
library(tidyverse)
library(parallel)
library(Hmisc)

## set up
set.seed(123432)
seeds <- floor(runif(50,1,1e5)) # set seeds

nocorr <- matrix(c(1,0,0,1), nrow=2, byrow=T)
negcorr <- matrix(c(1,-.5,-.5,1), nrow=2, byrow=T)
poscorr <- matrix(c(1,.5,.5,1), nrow=2, byrow=T)

hivar = c(   c=0.15,    shed=0.1,    alpha=0.15,    gamma=0.15, 
             sd_c=0.75, sd_shed=0.5, sd_alpha=0.75, sd_gamma=0.75, 
             b=2.5, d=.1, bs=.01)
medvar = c(   c=0.15,    shed=0.1,    alpha=0.15,    gamma=0.15, 
              sd_c=0.15, sd_shed=0.1, sd_alpha=0.15, sd_gamma=0.15, 
              b=2.5, d=.1, bs=.01)
lowvar = c(   c=0.15,     shed=0.1,    alpha=0.15,     gamma=0.15, 
              sd_c=0.03, sd_shed=0.02, sd_alpha=0.03, sd_gamma=0.03, 
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


# ## Plot median trajectories of S (expected equilibrium with no variation is S=63)
# png(filename="~/Desktop/Covariance_results.png", height=10, width=5, units='in', res=600)
par(mfrow=c(6,3), mar=c(1,1,1,1), oma=c(3,3,0,0))
for (j in 1:6) { ## loop over the six different covariance combinations
  for (covMatrix in c("nocorr", "negcorr", "poscorr")) {
    plot.new()
    plot.window(xlim=c(0,100), ylim=c(0,235))
    axis(1); axis(2); box('plot')
    z1 = get(paste("out",covMatrix,"hivar",var1[j],var2[j],sep="_"))
    z2 = get(paste("out",covMatrix,"medvar",var1[j],var2[j],sep="_"))
    z3 = get(paste("out",covMatrix,"lowvar",var1[j],var2[j],sep="_"))
    abline(h=63, lty=2)
    lines(0:100, lapply(z1, function(l) l[[1]]$S) %>% do.call("cbind.data.frame",.) %>% apply(., 1, median), lwd=2)
    lines(0:100, lapply(z2, function(l) l[[1]]$S) %>% do.call("cbind.data.frame",.) %>% apply(., 1, median), lwd=2, col=4)
    lines(0:100, lapply(z3, function(l) l[[1]]$S) %>% do.call("cbind.data.frame",.) %>% apply(., 1, median), lwd=2, col=2)
    legend(x='topright', paste(var1[j], var2[j], covMatrix, sep='-'), bty='n')
  }
}
# dev.off()
# 
#     

## To create an equivalent plot to the one below, just in ggplot, use the following code
var1=c('c','c','c','shed','shed','alpha')
var2=c('shed','alpha','gamma','alpha','gamma','gamma')
i <- 1
data <- vector(mode='list', length=18)
for (j in 1:6) { ## loop over the six different covariance combinations
  for (covMatrix in c("nocorr", "negcorr", "poscorr")) {
    z1 = readRDS(paste0(paste("out",covMatrix,"hivar",var1[j],var2[j],sep="_"),".RDS"))
    z2 = readRDS(paste0(paste("out",covMatrix,"medvar",var1[j],var2[j],sep="_"),".RDS"))
    z3 = readRDS(paste0(paste("out",covMatrix,"lowvar",var1[j],var2[j],sep="_"),".RDS"))
    
    ## Give a nicer name to covMatrix
    covMatrix <- switch(covMatrix,poscorr="(+) Cov",negcorr="(-) Cov",nocorr="(0) Cov")
    
    z <- vector(mode='list',length=3)
    lapply(1:length(z1), function(i) mutate(z1[[i]][[1]], rep=i, t=round(t))) %>% 
      do.call("rbind.data.frame",.) %>%
      mutate(., Variance="high",cov=covMatrix, traits=paste(var1[j],var2[j],sep="-")) -> z[[1]]
    lapply(1:length(z2), function(i) mutate(z2[[i]][[1]], rep=i, t=round(t))) %>% 
      do.call("rbind.data.frame",.) %>%
      mutate(., Variance="med",cov=covMatrix, traits=paste(var1[j],var2[j],sep="-")) -> z[[2]]
    lapply(1:length(z3), function(i) mutate(z3[[i]][[1]], rep=i, t=round(t))) %>% 
      do.call("rbind.data.frame",.) %>%
      mutate(., Variance="low",cov=covMatrix, traits=paste(var1[j],var2[j],sep="-")) -> z[[3]]
    data[[i]] <- do.call("rbind.data.frame",z)
    i <- i+1
  }
}
data <- do.call("rbind.data.frame",data)
saveRDS(data, file="population_dynamics_all_simulations.RDS")    

## To plot the dynamics of S across all trait pairs, covariance structures, and variance levels 
ggplot(data, aes(x=t, y=S, group=rep)) +
  stat_summary(aes(group=Variance), geom="ribbon", fun.data=mean_cl_normal, alpha=.45, fun.args=list(conf.int=0.95)) +
  stat_summary(aes(group=Variance,colour=Variance), geom="line", fun=mean) + 
  facet_grid(traits~cov) + 
  ylim(0,50) + 
  theme_bw()

ggplot(data, aes(x=t, y=I, group=rep)) +
  stat_summary(aes(group=Variance), geom="ribbon", fun.data=mean_cl_normal, alpha=.45, fun.args=list(conf.int=0.95)) +
  stat_summary(aes(group=Variance,colour=Variance), geom="line", fun=mean) + 
  facet_grid(traits~cov) + 
  ylim(0,200) + 
  theme_bw()

## To break out a particular trait pair
ggplot(subset(data, traits="c-gamma"), aes(x=t, y=S, group=rep)) +
  stat_summary(aes(group=Variance), geom="ribbon", fun.data=mean_cl_normal, alpha=.45, fun.args=list(conf.int=0.95)) +
  stat_summary(aes(group=Variance,colour=Variance), geom="line", fun=mean) + 
  facet_grid(.~cov) + 
  ylim(0,50) + 
  theme_bw()

ggplot(subset(data, traits="c-shed"), aes(x=t, y=S, group=rep)) +
  stat_summary(aes(group=Variance), geom="ribbon", fun.data=mean_cl_normal, alpha=.45, fun.args=list(conf.int=0.95)) +
  stat_summary(aes(group=Variance,colour=Variance), geom="line", fun=mean) + 
  facet_grid(.~cov) + 
  ylim(0,50) + 
  theme_bw()

ggplot(subset(data, traits%in%c("c-gamma","c-alpha")), aes(x=t, y=S, group=rep)) +
  stat_summary(aes(group=Variance), geom="ribbon", fun.data=mean_cl_normal, fun.args=list(conf.int=0.95)) +
  stat_summary(aes(group=Variance,colour=Variance), geom="line", fun=mean) + 
  facet_grid(cov~traits) + 
  ylim(0,50) + 
  theme_bw()

ggplot(subset(data, traits%in%c("c-gamma","c-alpha")), aes(x=t, y=I, group=rep)) +
  stat_summary(aes(group=Variance), geom="ribbon", fun.data=mean_cl_normal, fun.args=list(conf.int=0.95), alpha = .45) +
  stat_summary(aes(group=Variance,colour=Variance), geom="line", fun=mean) + 
  facet_grid(cov~traits) + 
  theme_bw()




