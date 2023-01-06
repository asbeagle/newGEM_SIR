## frac 0 ALL PLOTS

## shed gamma
i <- 1
sumdata1 <- vector(mode='list')
for (var in c("low","med","hi")) {
  for (corr in c("no","neg","pos")) {
    data <- readRDS(paste0("out_",corr,"corr_",var,"var_shed_gamma.RDS"))
    lapply(1:length(data), function(d)
      data.frame(R=subset(data[[d]][[2]], tInf <= data[[d]][[1]]$t[max(which(data[[d]][[1]]$S > 0.75*235))])$numInf,
                 epiSpeed=data[[d]][[1]]$t[max(which(data[[d]][[1]]$S > 0.75*235))],
                 Corr = as.factor(corr),
                 Var = as.factor(var),
                 Rep = d)) %>% do.call("rbind.data.frame",.) -> sumdata1[[i]]
    i <- i+1
  }
}
sumdata_use1 <-do.call("rbind.data.frame", sumdata1)


sumdata_use1 %>% group_by(Corr, Var) %>% summarise(meanR = mean(R), varR=var(R), fracSS=sum(R>(mean(R)+2*sd(R)))/length(R), 
                                                   fracSDS=sum(R>(mean(R)+4*sd(R)))/length(R), frac0=sum(R==0)/length(R),
                                                   meanSpeed=mean(epiSpeed)) -> sumtable_shed_gamma

shed_gamma<-ggplot(sumtable_shed_gamma, aes(x=frac0, y=meanR))+
  geom_point(data=sumtable_shed_gamma, aes(color=Var, shape=Corr), alpha = .85, size = 3.5)+
  scale_color_manual(values=c("darkblue", "darkgreen", "pink"))+
  labs(title="shed-gamma")+
  labs(x = "% 0 Reff", y = "mean Reff")+
  ylim(0.2, 3)+
  xlim(0.25, 1)+
  theme_bw()
shed_gamma

shed_gamma<-ggplot(sumtable_shed_gamma, aes(x=fracSDS, y=frac0))+
  geom_point(data=sumtable_shed_gamma, aes(color=Var, shape=Corr), alpha = .85, size = 3.5)+
  scale_color_manual(values=c("darkblue", "darkgreen", "pink"))+
  labs(title="shed-gamma")+
  labs(x = "% SDS", y = "% 0 Reff")+
  theme_bw()

shed_gamma

shed_gamma<-ggplot(sumtable_shed_gamma, aes(x=fracSS, y=frac0))+
  geom_point(data=sumtable_shed_gamma, aes(color=Var, shape=Corr), alpha = .85, size = 3.5)+
  scale_color_manual(values=c("darkblue", "darkgreen", "pink"))+
  labs(title="shed-gamma")+
  labs(x = "% SS", y = "% 0 Reff")+
  theme_bw()
shed_gamma

## shed alpha
i <- 1
sumdata <- vector(mode='list')
for (var in c("low","med","hi")) {
  for (corr in c("no","neg","pos")) {
    data <- readRDS(paste0("out_",corr,"corr_",var,"var_shed_alpha.RDS"))
    lapply(1:length(data), function(d)
      data.frame(R=subset(data[[d]][[2]], tInf <= data[[d]][[1]]$t[max(which(data[[d]][[1]]$S > 0.75*235))])$numInf,
                 epiSpeed=data[[d]][[1]]$t[max(which(data[[d]][[1]]$S > 0.75*235))],
                 Corr = as.factor(corr),
                 Var = as.factor(var),
                 Rep = d)) %>% do.call("rbind.data.frame",.) -> sumdata[[i]]
    i <- i+1
  }
}
sumdata_use <-do.call("rbind.data.frame", sumdata)


sumdata_use %>% group_by(Corr, Var) %>% summarise(meanR = mean(R), varR=var(R), fracSS=sum(R>(mean(R)+2*sd(R)))/length(R), 
                                                   fracSDS=sum(R>(mean(R)+4*sd(R)))/length(R), frac0=sum(R==0)/length(R),
                                                   meanSpeed=mean(epiSpeed)) -> sumtable_shed_alpha

shed_alpha<-ggplot(sumtable_shed_alpha, aes(x=frac0, y=meanR))+
  geom_point(data=sumtable_shed_alpha, aes(color=Var, shape=Corr), alpha = .85, size = 3.5)+
  scale_color_manual(values=c("darkblue", "darkgreen", "pink"))+
  labs(title="shed-alpha")+
  labs(x = "% 0 Reff", y = "mean Reff")+
  ylim(0.2, 3)+
  xlim(0.25, 1)+
  theme_bw()
shed_alpha

shed_alpha<-ggplot(sumtable_shed_alpha, aes(x=fracSS, y=frac0))+
  geom_point(data=sumtable_shed_alpha, aes(color=Var, shape=Corr), alpha = .85, size = 3.5)+
  scale_color_manual(values=c("darkblue", "darkgreen", "pink"))+
  labs(title="shed-alpha")+
  labs(x = "% SS", y = "% 0 Reff")+
  theme_bw()
shed_alpha

shed_alpha<-ggplot(sumtable_shed_alpha, aes(x=fracSDS, y=frac0))+
  geom_point(data=sumtable_shed_alpha, aes(color=Var, shape=Corr), alpha = .85, size = 3.5)+
  scale_color_manual(values=c("darkblue", "darkgreen", "pink"))+
  labs(title="shed-alpha")+
  labs(x = "% SDS", y = "% 0 Reff")+
  theme_bw()
shed_alpha

## contact shed
i <- 1
sumdata <- vector(mode='list')
for (var in c("low","med","hi")) {
  for (corr in c("no","neg","pos")) {
    data <- readRDS(paste0("out_",corr,"corr_",var,"var_c_shed.RDS"))
    lapply(1:length(data), function(d)
      data.frame(R=subset(data[[d]][[2]], tInf <= data[[d]][[1]]$t[max(which(data[[d]][[1]]$S > 0.75*235))])$numInf,
                 epiSpeed=data[[d]][[1]]$t[max(which(data[[d]][[1]]$S > 0.75*235))],
                 Corr = as.factor(corr),
                 Var = as.factor(var),
                 Rep = d)) %>% do.call("rbind.data.frame",.) -> sumdata[[i]]
    i <- i+1
  }
}
sumdata_use <-do.call("rbind.data.frame", sumdata)


sumdata_use %>% group_by(Corr, Var) %>% summarise(meanR = mean(R), varR=var(R), fracSS=sum(R>(mean(R)+2*sd(R)))/length(R), 
                                                  fracSDS=sum(R>(mean(R)+4*sd(R)))/length(R), frac0=sum(R==0)/length(R),
                                                  meanSpeed=mean(epiSpeed)) -> sumtable_c_shed

c_shed<-ggplot(sumtable_c_shed, aes(x=frac0, y=meanR))+
  geom_point(data=sumtable_c_shed, aes(color=Var, shape=Corr), alpha = .85, size = 3.5)+
  scale_color_manual(values=c("darkblue", "darkgreen", "pink"))+
  labs(title="contact-shed")+
  labs(x = "% 0 Reff", y = "mean Reff")+
  ylim(0.2, 3)+
  xlim(0.25, 1)+
  theme_bw()
c_shed

c_shed<-ggplot(sumtable_c_shed, aes(x=fracSDS, y=frac0))+
  geom_point(data=sumtable_c_shed, aes(color=Var, shape=Corr), alpha = .85, size = 3.5)+
  scale_color_manual(values=c("darkblue", "darkgreen", "pink"))+
  labs(title="contact-shed")+
  labs(x = "% SDS", y = "% 0 Reff")+
  theme_bw()
c_shed

c_shed<-ggplot(sumtable_c_shed, aes(x=fracSS, y=frac0))+
  geom_point(data=sumtable_c_shed, aes(color=Var, shape=Corr), alpha = .85, size = 3.5)+
  scale_color_manual(values=c("darkblue", "darkgreen", "pink"))+
  labs(title="contact-shed")+
  labs(x = "% SS", y = "% 0 Reff")+
  theme_bw()
c_shed

## contact gamma
i <- 1
sumdata <- vector(mode='list')
for (var in c("low","med","hi")) {
  for (corr in c("no","neg","pos")) {
    data <- readRDS(paste0("out_",corr,"corr_",var,"var_c_gamma.RDS"))
    lapply(1:length(data), function(d)
      data.frame(R=subset(data[[d]][[2]], tInf <= data[[d]][[1]]$t[max(which(data[[d]][[1]]$S > 0.75*235))])$numInf,
                 epiSpeed=data[[d]][[1]]$t[max(which(data[[d]][[1]]$S > 0.75*235))],
                 Corr = as.factor(corr),
                 Var = as.factor(var),
                 Rep = d)) %>% do.call("rbind.data.frame",.) -> sumdata[[i]]
    i <- i+1
  }
}
sumdata_use <-do.call("rbind.data.frame", sumdata)


sumdata_use %>% group_by(Corr, Var) %>% summarise(meanR = mean(R), varR=var(R), fracSS=sum(R>(mean(R)+2*sd(R)))/length(R), 
                                                  fracSDS=sum(R>(mean(R)+4*sd(R)))/length(R), frac0=sum(R==0)/length(R),
                                                  meanSpeed=mean(epiSpeed)) -> sumtable_c_gamma

c_gamma<-ggplot(sumtable_c_gamma, aes(x=frac0, y=meanR))+
  geom_point(data=sumtable_c_gamma, aes(color=Var, shape=Corr), alpha = .85, size = 3.5)+
  scale_color_manual(values=c("darkblue", "darkgreen", "pink"))+
  labs(title="contact-gamma")+
  labs(x = "% 0 Reff", y = "mean Reff")+
  ylim(0.2, 3)+
  xlim(0.25, 1)+
  theme_bw()
c_gamma

c_gamma<- ggplot(sumtable_c_gamma, aes(x=fracSDS, y=frac0))+
  geom_point(data=sumtable_c_gamma, aes(color=Var, shape=Corr), alpha = .85, size = 3.5)+
  scale_color_manual(values=c("darkblue", "darkgreen", "pink"))+
  labs(title="contact-gamma")+
  labs(x = "% SDS", y = "% 0 Reff")+
  theme_bw()
c_gamma

c_gamma<- ggplot(sumtable_c_gamma, aes(x=fracSS, y=frac0))+
  geom_point(data=sumtable_c_gamma, aes(color=Var, shape=Corr), alpha = .85, size = 3.5)+
  scale_color_manual(values=c("darkblue", "darkgreen", "pink"))+
  labs(title="contact-gamma")+
  labs(x = "% SS", y = "% 0 Reff")+
  #xlim(0, 0.02)+
  theme_bw()
c_gamma

## contact alpha
i <- 1
sumdata <- vector(mode='list')
for (var in c("low","med","hi")) {
  for (corr in c("no","neg","pos")) {
    data <- readRDS(paste0("out_",corr,"corr_",var,"var_c_alpha.RDS"))
    lapply(1:length(data), function(d)
      data.frame(R=subset(data[[d]][[2]], tInf <= data[[d]][[1]]$t[max(which(data[[d]][[1]]$S > 0.75*235))])$numInf,
                 epiSpeed=data[[d]][[1]]$t[max(which(data[[d]][[1]]$S > 0.75*235))],
                 Corr = as.factor(corr),
                 Var = as.factor(var),
                 Rep = d)) %>% do.call("rbind.data.frame",.) -> sumdata[[i]]
    i <- i+1
  }
}
sumdata_use <-do.call("rbind.data.frame", sumdata)


sumdata_use %>% group_by(Corr, Var) %>% summarise(meanR = mean(R), varR=var(R), fracSS=sum(R>(mean(R)+2*sd(R)))/length(R), 
                                                  fracSDS=sum(R>(mean(R)+4*sd(R)))/length(R), frac0=sum(R==0)/length(R),
                                                  meanSpeed=mean(epiSpeed)) -> sumtable_c_alpha

c_alpha<- ggplot(sumtable_c_alpha, aes(x=frac0, y=meanR))+
  geom_point(data=sumtable_c_alpha, aes(color=Var, shape=Corr), alpha = .85, size = 3.5)+
  scale_color_manual(values=c("darkblue", "darkgreen", "pink"))+
  labs(title="contact-alpha")+
  labs(x = "% 0 Reff", y = "mean Reff")+
  ylim(0.2, 3)+
  xlim(0.25, 1)+
  theme_bw()
c_alpha

c_alpha<-ggplot(sumtable_c_alpha, aes(x=fracSDS, y=frac0))+
  geom_point(data=sumtable_c_alpha, aes(color=Var, shape=Corr), alpha = .85, size = 3.5)+
  scale_color_manual(values=c("darkblue", "darkgreen", "pink"))+
  labs(title="contact-alpha")+
  labs(x = "% SDS", y = "% 0 Reff")+
  theme_bw()
c_alpha

c_alpha<-ggplot(sumtable_c_alpha, aes(x=fracSS, y=frac0))+
  geom_point(data=sumtable_c_alpha, aes(color=Var, shape=Corr), alpha = .85, size = 3.5)+
  scale_color_manual(values=c("darkblue", "darkgreen", "pink"))+
  labs(title="contact-alpha")+
  labs(x = "% SS", y = "% 0 Reff")+
  #xlim(0, 0.02)+
  theme_bw()
c_alpha

## alpha gamma
i <- 1
sumdata <- vector(mode='list')
for (var in c("low","med","hi")) {
  for (corr in c("no","neg","pos")) {
    data <- readRDS(paste0("out_",corr,"corr_",var,"var_alpha_gamma.RDS"))
    lapply(1:length(data), function(d)
      data.frame(R=subset(data[[d]][[2]], tInf <= data[[d]][[1]]$t[max(which(data[[d]][[1]]$S > 0.75*235))])$numInf,
                 epiSpeed=data[[d]][[1]]$t[max(which(data[[d]][[1]]$S > 0.75*235))],
                 Corr = as.factor(corr),
                 Var = as.factor(var),
                 Rep = d)) %>% do.call("rbind.data.frame",.) -> sumdata[[i]]
    i <- i+1
  }
}
sumdata_use <-do.call("rbind.data.frame", sumdata)


sumdata_use %>% group_by(Corr, Var) %>% summarise(meanR = mean(R), varR=var(R), fracSS=sum(R>(mean(R)+2*sd(R)))/length(R), 
                                                  fracSDS=sum(R>(mean(R)+4*sd(R)))/length(R), frac0=sum(R==0)/length(R),
                                                  meanSpeed=mean(epiSpeed)) -> sumtable_alpha_gamma

alpha_gamma<- ggplot(sumtable_alpha_gamma, aes(x=frac0, y=meanR))+
  geom_point(data=sumtable_alpha_gamma, aes(color=Var, shape=Corr), alpha = .85, size = 3.5)+
  scale_color_manual(values=c("darkblue", "darkgreen", "pink"))+
  labs(title="alpha-gamma")+
  labs(x = "% 0 Reff", y = "mean Reff")+
  ylim(0.2, 3)+
  xlim(0.25, 1)+
  theme_bw()
alpha_gamma

alpha_gamma<-ggplot(sumtable_alpha_gamma, aes(x=fracSDS, y=frac0))+
  geom_point(data=sumtable_alpha_gamma, aes(color=Var, shape=Corr), alpha = .85, size = 3.5)+
  scale_color_manual(values=c("darkblue", "darkgreen", "pink"))+
  labs(title="alpha-gamma")+
  labs(x = "% SDS", y = "% 0 Reff")+
  theme_bw()
alpha_gamma

alpha_gamma<-ggplot(sumtable_alpha_gamma, aes(x=fracSS, y=frac0))+
  geom_point(data=sumtable_alpha_gamma, aes(color=Var, shape=Corr), alpha = .85, size = 3.5)+
  scale_color_manual(values=c("darkblue", "darkgreen", "pink"))+
  labs(title="alpha-gamma")+
  labs(x = "% SS", y = "% 0 Reff")+
  #xlim(0, 0.02)+
  theme_bw()
alpha_gamma


library(tidyverse)
library(cowplot)
library(ggplot2)
library(ggpubr)

ggarrange(c_shed, c_alpha, c_gamma,  alpha_gamma, shed_alpha, shed_gamma, nrow=2, ncol=3, common.legend = TRUE, legend="right")
