source("GEM_SIR_cov_storage.R")
library(tidyverse)
library(parallel)

## set up
set.seed(123432)
seeds <- floor(runif(100,1,1e5)) # set seeds

nocorr <- 0
negcorr <- -0.5
poscorr <- 0.5

## The baseline parameters give R = 8. To scale the parameters for R = 4 and R = 1 use
## Solve[(((1 - p) c ((1 - p) s)/(1 + (1 - p) s) (b - d)/bs)/(d + (1 + p) a + (1 + p) g) /. params) == R, p]

for (R in c(8,4,1)) {
  if (R==8) x <- 0
  else if (R==4) x <- 0.24647
  else x <- 0.594804
  
  hivar = c(   c=(1-x)*0.1,    shed=(1-x)*1/9,    alpha=(1+x)*0.1,    gamma=(1+x)*0.1, 
               sd_c=(1-x)*0.5, sd_shed=(1-x)*5/9, sd_alpha=(1+x)*0.5, sd_gamma=(1+x)*0.5, 
               b=2.5, d=.1, bs=.01)
  medvar = c(  c=(1-x)*0.1,    shed=(1-x)*1/9,    alpha=(1+x)*0.1,    gamma=(1+x)*0.1, 
               sd_c=(1-x)*0.1, sd_shed=(1-x)*1/9, sd_alpha=(1+x)*0.1, sd_gamma=(1+x)*0.1, 
               b=2.5, d=.1, bs=.01)
  lowvar = c(  c=(1-x)*0.1,    shed=(1-x)*1/9,      alpha=(1+x)*0.1,    gamma=(1+x)*0.1, 
               sd_c=(1-x)*0.02, sd_shed=(1-x)*1/45, sd_alpha=(1+x)*0.02, sd_gamma=(1+x)*0.02, 
               b=2.5, d=.1, bs=.01)

  initial_state <- floor(c(S=unname(((hivar["b"]-hivar["d"])/hivar["bs"]))-5, I=5, R=0))
  
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

# var1=c('c','c','c','shed','shed','alpha')
# var2=c('shed','alpha','gamma','alpha','gamma','gamma')
# i <- 1
# data <- vector(mode='list', length=18)
# for (j in 1:6) { ## loop over the six different covariance combinations
#   for (covMatrix in c("nocorr", "negcorr", "poscorr")) {
#     z1 = readRDS(paste0(paste("out",covMatrix,"hivar",var1[j],var2[j],sep="_"),".RDS"))
#     z2 = readRDS(paste0(paste("out",covMatrix,"medvar",var1[j],var2[j],sep="_"),".RDS"))
#     z3 = readRDS(paste0(paste("out",covMatrix,"lowvar",var1[j],var2[j],sep="_"),".RDS"))
#     
#     ## Give a nicer name to covMatrix
#     covMatrix <- switch(covMatrix,poscorr="(+) Cov",negcorr="(-) Cov",nocorr="(0) Cov")
#     
#     z <- vector(mode='list',length=3)
#     lapply(1:length(z1), function(i) mutate(z1[[i]][[1]][1:(ifelse(any(z1[[i]][[1]]$I==0),min(which(z1[[i]][[1]]$I==0)),length(z1[[i]][[1]]$I))),],
#                                             rep=i,t=ceiling(t))) %>% 
#       do.call("rbind.data.frame",.) %>%
#       mutate(., Variance="high",cov=covMatrix, traits=paste(var1[j],var2[j],sep="-")) -> z[[1]]
#     lapply(1:length(z2), function(i) mutate(z2[[i]][[1]][1:(ifelse(any(z2[[i]][[1]]$I==0),min(which(z2[[i]][[1]]$I==0)),length(z2[[i]][[1]]$I))),],
#                                             rep=i,t=ceiling(t))) %>% 
#       do.call("rbind.data.frame",.) %>%
#       mutate(., Variance="med",cov=covMatrix, traits=paste(var1[j],var2[j],sep="-")) -> z[[2]]
#     lapply(1:length(z3), function(i) mutate(z3[[i]][[1]][1:(ifelse(any(z3[[i]][[1]]$I==0),min(which(z3[[i]][[1]]$I==0)),length(z3[[i]][[1]]$I))),],
#                                              rep=i,t=ceiling(t))) %>% 
#       do.call("rbind.data.frame",.) %>%
#       mutate(., Variance="low",cov=covMatrix, traits=paste(var1[j],var2[j],sep="-")) -> z[[3]]
#     data[[i]] <- do.call("rbind.data.frame",z)
#     i <- i+1
#   }
# }
# data <- do.call("rbind.data.frame",data)
# saveRDS(data, file="population_dynamics_all_simulations_R=1.RDS")    
# 
# data2= filter(data, Variance=="high", cov=="(0) Cov", traits=="c-gamma")
# ggplot(data2, aes(x=t, y=S, group=rep)) + geom_line()
# ggplot(data2, aes(x=t, y=S, group=rep)) + stat_summary(aes(group=Variance), geom="line", fun=mean)
# data2 %>% group_by(t) %>% summarise(meanS=mean(S),medS=median(S)) %>% ggplot(.,aes(x=t,y=medS)) + geom_line()
# 
# ## To plot the dynamics of S across all trait pairs, covariance structures, and variance levels
# ggplot(data, aes(x=t, y=S, group=rep)) +
#   stat_summary(aes(group=Variance), geom="ribbon", fun.data=mean_cl_normal, alpha=.45, fun.args=list(conf.int=0.95)) +
#   stat_summary(aes(group=Variance,colour=Variance), geom="line", fun=mean) + 
#   facet_grid(traits~cov) + 
#   ylim(0,250) + 
#   theme_bw()
# 
# ggplot(data, aes(x=t, y=I, group=rep)) +
#   stat_summary(aes(group=Variance), geom="ribbon", fun.data=mean_cl_normal, alpha=.45, fun.args=list(conf.int=0.95)) +
#   stat_summary(aes(group=Variance,colour=Variance), geom="line", fun=mean) + 
#   facet_grid(traits~cov) + 
#   ylim(0,200) + 
#   theme_bw()
# 
# ## To break out a particular trait pair
# ggplot(subset(data, traits="c-gamma"), aes(x=t, y=S, group=rep)) +
#   stat_summary(aes(group=Variance), geom="ribbon", fun.data=mean_cl_normal, alpha=.45, fun.args=list(conf.int=0.95)) +
#   stat_summary(aes(group=Variance,colour=Variance), geom="line", fun=mean) + 
#   facet_grid(.~cov) +
#   theme_bw()
# 
# ggplot(data, aes(x=t, y=S, group=rep)) +
#   stat_summary(aes(group=Variance), geom="ribbon", fun.data=mean_cl_normal, alpha=.45, fun.args=list(conf.int=0.95)) +
#   stat_summary(aes(group=Variance,colour=Variance), geom="line", fun=mean) + 
#   facet_grid(.~cov) +
#   theme_bw()
# 
# 
# 
# low_R0_S<- ggplot(subset(data, traits=="c-shed"), aes(x=t, y=S, group=rep)) +
#   stat_summary(aes(group=Variance), geom="ribbon", fun.data=mean_cl_normal, alpha=.45, fun.args=list(conf.int=0.95)) +
#   stat_summary(aes(group=Variance,colour=Variance), geom="line", fun=mean) + 
#   facet_grid(.~cov) + 
#   scale_color_manual(values=c("darkblue", "darkgreen", "pink"))+
#   theme_bw()+
#   ylim(0,250)
# 
# low_R0_I<-ggplot(subset(data, traits=="c-shed"), aes(x=t, y=I, group=rep)) +
#   stat_summary(aes(group=Variance), geom="ribbon", fun.data=mean_cl_normal, alpha=.45, fun.args=list(conf.int=0.95)) +
#   stat_summary(aes(group=Variance,colour=Variance), geom="line", fun=mean) + 
#   facet_grid(.~cov) + 
#   scale_color_manual(values=c("darkblue", "darkgreen", "pink"))+
#   theme_bw()+
#   ylim(0,200)
# 
# low_R0<-ggarrange(low_R0_S, low_R0_I, nrow=2, ncol=1, common.legend = TRUE, legend = "top",
#           labels = "Low R0",hjust = -.85, vjust=-0.05, font.label=list(size=15, face="plain") )
# 
# low_R0
# 
# med_R0_S<- ggplot(subset(data, traits="c-shed"), aes(x=t, y=S, group=rep)) +
#   stat_summary(aes(group=Variance), geom="ribbon", fun.data=mean_cl_normal, alpha=.45, fun.args=list(conf.int=0.95)) +
#   stat_summary(aes(group=Variance,colour=Variance), geom="line", fun=mean) + 
#   facet_grid(.~cov) + 
#   scale_color_manual(values=c("darkblue", "darkgreen", "pink"))+
#   theme_bw()+
#   ylim(0,250)
# 
# 
# med_R0_I<-ggplot(subset(data, traits="c-shed"), aes(x=t, y=I, group=rep)) +
#   stat_summary(aes(group=Variance), geom="ribbon", fun.data=mean_cl_normal, alpha=.45, fun.args=list(conf.int=0.95)) +
#   stat_summary(aes(group=Variance,colour=Variance), geom="line", fun=mean) + 
#   facet_grid(.~cov) + 
#   scale_color_manual(values=c("darkblue", "darkgreen", "pink"))+
#   theme_bw()+
#   ylim(0,200)
# med_R0_I
# 
# med_R0<-ggarrange(med_R0_S, med_R0_I, nrow=2, ncol=1, common.legend = TRUE, legend = "top",
#                   labels = "Moderate R0",hjust = -.5, vjust=-0.05, font.label=list(size=15, face="plain") )
# med_R0
# 
# ggplot(subset(data, traits="c-shed"), aes(x=t, y=S, group=rep)) +
#   stat_summary(aes(group=Variance), geom="ribbon", fun.data=mean_cl_normal, alpha=.45, fun.args=list(conf.int=0.95)) +
#   stat_summary(aes(group=Variance,colour=Variance), geom="line", fun=mean) + 
#   facet_grid(.~cov) + 
#   scale_color_manual(values=c("darkblue", "darkgreen", "pink"))+
#   theme_bw()
# 
# 
# hi_R0_I<-ggplot(subset(data, traits="c-shed"), aes(x=t, y=I, group=rep)) +
#   stat_summary(aes(group=Variance), geom="ribbon", fun.data=mean_cl_normal, alpha=.45, fun.args=list(conf.int=0.95)) +
#   stat_summary(aes(group=Variance,colour=Variance), geom="line", fun=mean) + 
#   facet_grid(.~cov) + 
#   scale_color_manual(values=c("darkblue", "darkgreen", "pink"))+
#   theme_bw()
# 
# hi_R0<-ggarrange(hi_R0_S, hi_R0_I, nrow=2, ncol=1, common.legend = TRUE, legend = "top",
#                   labels = "High R0",hjust = -.8, vjust=-0.05, font.label=list(size=15, face="plain") )
# hi_R0
# 
# 
# 
# 
# 
# 
# 
# ggplot(subset(data, traits%in%c("c-gamma","c-alpha")), aes(x=t, y=S, group=rep)) +
#   stat_summary(aes(group=Variance), geom="ribbon", fun.data=mean_cl_normal, fun.args=list(conf.int=0.95)) +
#   stat_summary(aes(group=Variance,colour=Variance), geom="line", fun=mean) + 
#   facet_grid(cov~traits) + 
#   scale_color_manual(values=c("darkblue", "darkgreen", "pink"))+
#   ylim(0,230)+
#   theme_bw()
# 
# ggplot(subset(data, traits%in%c("c-gamma","c-alpha")), aes(x=t, y=I, group=rep)) +
#   stat_summary(aes(group=Variance), geom="ribbon", fun.data=mean_cl_normal, fun.args=list(conf.int=0.95), alpha = .45) +
#   stat_summary(aes(group=Variance,colour=Variance), geom="line", fun=mean) + 
#   facet_grid(cov~traits) + 
#   theme_bw()
# 
#       
# 
# 
