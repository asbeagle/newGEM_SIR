### S & I dynamic plots that I want to include in the manuscript

### Hi R0
## To create an equivalent plot to the one below, just in ggplot, use the following code
var1=c('c','c','c','shed','shed','alpha')
var2=c('shed','alpha','gamma','alpha','gamma','gamma')
i <- 1
data_hiR0 <- vector(mode='list', length=18)
for (j in 1:6) { ## loop over the six different covariance combinations
  for (covMatrix in c("nocorr", "negcorr", "poscorr")) {
    z1 = readRDS(paste0(paste("out",covMatrix,"hivar",var1[j],var2[j],sep="_"),".RDS"))
    z2 = readRDS(paste0(paste("out",covMatrix,"medvar",var1[j],var2[j],sep="_"),".RDS"))
    z3 = readRDS(paste0(paste("out",covMatrix,"lowvar",var1[j],var2[j],sep="_"),".RDS"))
    
    ## Give a nicer name to covMatrix
    covMatrix <- switch(covMatrix,poscorr="(+) Cov",negcorr="(-) Cov",nocorr="(0) Cov")
    
    z <- vector(mode='list',length=3)
    lapply(1:length(z1), function(i) mutate(z1[[i]][[1]][1:(ifelse(any(z1[[i]][[1]]$I==0),min(which(z1[[i]][[1]]$I==0)),length(z1[[i]][[1]]$I))),],
                                            rep=i,t=ceiling(t))) %>% 
      do.call("rbind.data.frame",.) %>%
      mutate(., Variance="high",cov=covMatrix, traits=paste(var1[j],var2[j],sep="-")) -> z[[1]]
    lapply(1:length(z2), function(i) mutate(z2[[i]][[1]][1:(ifelse(any(z2[[i]][[1]]$I==0),min(which(z2[[i]][[1]]$I==0)),length(z2[[i]][[1]]$I))),],
                                            rep=i,t=ceiling(t))) %>% 
      do.call("rbind.data.frame",.) %>%
      mutate(., Variance="med",cov=covMatrix, traits=paste(var1[j],var2[j],sep="-")) -> z[[2]]
    lapply(1:length(z3), function(i) mutate(z3[[i]][[1]][1:(ifelse(any(z3[[i]][[1]]$I==0),min(which(z3[[i]][[1]]$I==0)),length(z3[[i]][[1]]$I))),],
                                            rep=i,t=ceiling(t))) %>% 
      do.call("rbind.data.frame",.) %>%
      mutate(., Variance="low",cov=covMatrix, traits=paste(var1[j],var2[j],sep="-")) -> z[[3]]
    data_hiR0[[i]] <- do.call("rbind.data.frame",z)
    i <- i+1
  }
}
data_hiR0 <- do.call("rbind.data.frame",data_hiR0)
saveRDS(data, file="population_dynamics_all_simulations_lowvar.RDS")

### moderate R0
## To create an equivalent plot to the one below, just in ggplot, use the following code
var1=c('c','c','c','shed','shed','alpha')
var2=c('shed','alpha','gamma','alpha','gamma','gamma')
i <- 1
data_medR0 <- vector(mode='list', length=18)
for (j in 1:6) { ## loop over the six different covariance combinations
  for (covMatrix in c("nocorr", "negcorr", "poscorr")) {
    z1 = readRDS(paste0(paste("out",covMatrix,"hivar",var1[j],var2[j],sep="_"),".RDS"))
    z2 = readRDS(paste0(paste("out",covMatrix,"medvar",var1[j],var2[j],sep="_"),".RDS"))
    z3 = readRDS(paste0(paste("out",covMatrix,"lowvar",var1[j],var2[j],sep="_"),".RDS"))
    
    ## Give a nicer name to covMatrix
    covMatrix <- switch(covMatrix,poscorr="(+) Cov",negcorr="(-) Cov",nocorr="(0) Cov")
    
    z <- vector(mode='list',length=3)
    lapply(1:length(z1), function(i) mutate(z1[[i]][[1]][1:(ifelse(any(z1[[i]][[1]]$I==0),min(which(z1[[i]][[1]]$I==0)),length(z1[[i]][[1]]$I))),],
                                            rep=i,t=ceiling(t))) %>% 
      do.call("rbind.data.frame",.) %>%
      mutate(., Variance="high",cov=covMatrix, traits=paste(var1[j],var2[j],sep="-")) -> z[[1]]
    lapply(1:length(z2), function(i) mutate(z2[[i]][[1]][1:(ifelse(any(z2[[i]][[1]]$I==0),min(which(z2[[i]][[1]]$I==0)),length(z2[[i]][[1]]$I))),],
                                            rep=i,t=ceiling(t))) %>% 
      do.call("rbind.data.frame",.) %>%
      mutate(., Variance="med",cov=covMatrix, traits=paste(var1[j],var2[j],sep="-")) -> z[[2]]
    lapply(1:length(z3), function(i) mutate(z3[[i]][[1]][1:(ifelse(any(z3[[i]][[1]]$I==0),min(which(z3[[i]][[1]]$I==0)),length(z3[[i]][[1]]$I))),],
                                            rep=i,t=ceiling(t))) %>% 
      do.call("rbind.data.frame",.) %>%
      mutate(., Variance="low",cov=covMatrix, traits=paste(var1[j],var2[j],sep="-")) -> z[[3]]
    data_medR0[[i]] <- do.call("rbind.data.frame",z)
    i <- i+1
  }
}
data_medR0 <- do.call("rbind.data.frame",data_medR0)
saveRDS(data, file="population_dynamics_all_simulations_lowvar.RDS")

### low R0
## To create an equivalent plot to the one below, just in ggplot, use the following code
var1=c('c','c','c','shed','shed','alpha')
var2=c('shed','alpha','gamma','alpha','gamma','gamma')
i <- 1
data_lowR0 <- vector(mode='list', length=18)
for (j in 1:6) { ## loop over the six different covariance combinations
  for (covMatrix in c("nocorr", "negcorr", "poscorr")) {
    z1 = readRDS(paste0(paste("out",covMatrix,"hivar",var1[j],var2[j],sep="_"),".RDS"))
    z2 = readRDS(paste0(paste("out",covMatrix,"medvar",var1[j],var2[j],sep="_"),".RDS"))
    z3 = readRDS(paste0(paste("out",covMatrix,"lowvar",var1[j],var2[j],sep="_"),".RDS"))
    
    ## Give a nicer name to covMatrix
    covMatrix <- switch(covMatrix,poscorr="(+) Cov",negcorr="(-) Cov",nocorr="(0) Cov")
    
    z <- vector(mode='list',length=3)
    lapply(1:length(z1), function(i) mutate(z1[[i]][[1]][1:(ifelse(any(z1[[i]][[1]]$I==0),min(which(z1[[i]][[1]]$I==0)),length(z1[[i]][[1]]$I))),],
                                            rep=i,t=ceiling(t))) %>% 
      do.call("rbind.data.frame",.) %>%
      mutate(., Variance="high",cov=covMatrix, traits=paste(var1[j],var2[j],sep="-")) -> z[[1]]
    lapply(1:length(z2), function(i) mutate(z2[[i]][[1]][1:(ifelse(any(z2[[i]][[1]]$I==0),min(which(z2[[i]][[1]]$I==0)),length(z2[[i]][[1]]$I))),],
                                            rep=i,t=ceiling(t))) %>% 
      do.call("rbind.data.frame",.) %>%
      mutate(., Variance="med",cov=covMatrix, traits=paste(var1[j],var2[j],sep="-")) -> z[[2]]
    lapply(1:length(z3), function(i) mutate(z3[[i]][[1]][1:(ifelse(any(z3[[i]][[1]]$I==0),min(which(z3[[i]][[1]]$I==0)),length(z3[[i]][[1]]$I))),],
                                            rep=i,t=ceiling(t))) %>% 
      do.call("rbind.data.frame",.) %>%
      mutate(., Variance="low",cov=covMatrix, traits=paste(var1[j],var2[j],sep="-")) -> z[[3]]
    data_lowR0[[i]] <- do.call("rbind.data.frame",z)
    i <- i+1
  }
}
data_lowR0 <- do.call("rbind.data.frame",data_lowR0)
saveRDS(data, file="population_dynamics_all_simulations_lowvar.RDS")

## To plot the dynamics of S across all trait pairs, covariance structures, and variance levels 
ggplot(data, aes(x=t, y=S, group=rep)) +
  stat_summary(aes(group=Variance), geom="ribbon", fun.data=mean_cl_normal, alpha=.45, fun.args=list(conf.int=0.95)) +
  stat_summary(aes(group=Variance,colour=Variance), geom="line", fun=mean) + 
  facet_grid(traits~cov) + 
  ylim(0,250) + 
  theme_bw()

plotz <- ggplot(data, aes(x=t, y=I, group=rep)) +
  stat_summary(aes(group=Variance), geom="ribbon", fun.data=mean_cl_normal, alpha=.45, fun.args=list(conf.int=0.95)) +
  stat_summary(aes(group=Variance,colour=Variance), geom="line", fun=mean) + 
  facet_grid(traits~cov) + 
  ylim(0,200) + 
  theme_bw() 

### plot individual cases
## contact - shed
low_R0_S<- ggplot(subset(data, traits=="c-shed"), aes(x=t, y=S, group=rep)) +
  stat_summary(aes(group=Variance), geom="ribbon", fun.data=mean_cl_normal, alpha=.45, fun.args=list(conf.int=0.95)) +
  stat_summary(aes(group=Variance,colour=Variance), geom="line", fun=mean) + 
  facet_grid(.~cov) + 
  scale_color_manual(values=c("darkblue", "darkgreen", "pink"))+
  theme_bw()+
  ylim(0,250)

low_R0_I<-ggplot(subset(data, traits=="c-shed"), aes(x=t, y=I, group=rep)) +
  stat_summary(aes(group=Variance), geom="ribbon", fun.data=mean_cl_normal, alpha=.45, fun.args=list(conf.int=0.95)) +
  stat_summary(aes(group=Variance,colour=Variance), geom="line", fun=mean) + 
  facet_grid(.~cov) + 
  scale_color_manual(values=c("darkblue", "darkgreen", "pink"))+
  theme_bw()+
  ylim(0,250)

low_R0<-ggarrange(low_R0_S, low_R0_I, nrow=2, ncol=1, common.legend = TRUE, legend = "top",
                  labels = "Low R0",hjust = -.85, vjust=-0.05, font.label=list(size=15, face="plain") )

low_R0

med_R0_S<- ggplot(subset(data, traits=="c-shed"), aes(x=t, y=S, group=rep)) +
  stat_summary(aes(group=Variance), geom="ribbon", fun.data=mean_cl_normal, alpha=.45, fun.args=list(conf.int=0.95)) +
  stat_summary(aes(group=Variance,colour=Variance), geom="line", fun=mean) + 
  facet_grid(.~cov) + 
  scale_color_manual(values=c("darkblue", "darkgreen", "pink"))+
  theme_bw()


med_R0_I<-ggplot(subset(data, traits=="c-shed"), aes(x=t, y=I, group=rep)) +
  stat_summary(aes(group=Variance), geom="ribbon", fun.data=mean_cl_normal, alpha=.45, fun.args=list(conf.int=0.95)) +
  stat_summary(aes(group=Variance,colour=Variance), geom="line", fun=mean) + 
  facet_grid(.~cov) + 
  scale_color_manual(values=c("darkblue", "darkgreen", "pink"))+
  theme_bw()
med_R0_I

med_R0<-ggarrange(med_R0_S, med_R0_I, nrow=2, ncol=1, common.legend = TRUE, legend = "top",
                  labels = "Moderate R0",hjust = -.5, vjust=-0.05, font.label=list(size=15, face="plain") )
med_R0

hi_R0_S<-ggplot(subset(data, traits=="c-shed"), aes(x=t, y=S, group=rep)) +
  stat_summary(aes(group=Variance), geom="ribbon", fun.data=mean_cl_normal, alpha=.45, fun.args=list(conf.int=0.95)) +
  stat_summary(aes(group=Variance,colour=Variance), geom="line", fun=mean) + 
  facet_grid(.~cov) + 
  scale_color_manual(values=c("darkblue", "darkgreen", "pink"))+
  theme_bw()


hi_R0_I<-ggplot(subset(data, traits="c-shed"), aes(x=t, y=I, group=rep)) +
  stat_summary(aes(group=Variance), geom="ribbon", fun.data=mean_cl_normal, alpha=.45, fun.args=list(conf.int=0.95)) +
  stat_summary(aes(group=Variance,colour=Variance), geom="line", fun=mean) + 
  facet_grid(.~cov) + 
  scale_color_manual(values=c("darkblue", "darkgreen", "pink"))+
  theme_bw()

hi_R0<-ggarrange(hi_R0_S, hi_R0_I, nrow=2, ncol=1, common.legend = TRUE, legend = "top",
                 labels = "High R0",hjust = -.8, vjust=-0.05, font.label=list(size=15, face="plain") )
hi_R0

ggarrange(low_R0, med_R0, hi_R0, nrow=3, ncol=1)


#### MAKE PLOTS WITHOUT COVARIATION PANELS :D

low_R0_I_nocov<-ggplot(subset(data_lowR0, cov=="(0) Cov"), aes(x=t, y=I, group=rep)) +
  stat_summary(aes(group=Variance), geom="ribbon", fun.data=mean_cl_normal, alpha=.45, fun.args=list(conf.int=0.95)) +
  stat_summary(aes(group=Variance,colour=Variance), geom="line", fun=mean) + 
  facet_grid(traits~cov) + 
  ylim(0,200) + 
  scale_color_manual(values=c("darkblue", "darkgreen", "pink"))+
  theme_bw()+
  ggtitle("Low R0")
low_R0_I_nocov


med_R0_I_nocov<-ggplot(subset(data_medR0, cov=="(0) Cov"), aes(x=t, y=I, group=rep)) +
  stat_summary(aes(group=Variance), geom="ribbon", fun.data=mean_cl_normal, alpha=.45, fun.args=list(conf.int=0.95)) +
  stat_summary(aes(group=Variance,colour=Variance), geom="line", fun=mean) + 
  facet_grid(traits~cov) + 
  ylim(0,200) + 
  scale_color_manual(values=c("darkblue", "darkgreen", "pink"))+
  theme_bw()+
  ggtitle("Moderate R0")
med_R0_I_nocov


hi_R0_I_nocov<-ggplot(subset(data_hiR0, cov=="(0) Cov"), aes(x=t, y=I, group=rep)) +
  stat_summary(aes(group=Variance), geom="ribbon", fun.data=mean_cl_normal, alpha=.45, fun.args=list(conf.int=0.95)) +
  stat_summary(aes(group=Variance,colour=Variance), geom="line", fun=mean) + 
  facet_grid(traits~cov) + 
  #ylim(0,200) + 
  scale_color_manual(values=c("darkblue", "darkgreen", "pink"))+
  theme_bw()+
  ggtitle("High R0")
hi_R0_I_nocov


ggarrange(low_R0_I_nocov, med_R0_I_nocov,hi_R0_I_nocov, common.legend = TRUE, nrow=1, ncol=3)


