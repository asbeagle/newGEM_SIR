### plots with all levels of R0

library(tidyverse)

## low R0
## hi var
## shed gamma
i <- 1
sumdata1 <- vector(mode='list')
for (var in c("hi")) {
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
                                                   meanSpeed=mean(epiSpeed)) %>% mutate(., R0 = "low")-> sumtable_shed_gamma_lowR0

## med R0
i <- 1
sumdata1 <- vector(mode='list')
for (var in c("hi")) {
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
                                                   meanSpeed=mean(epiSpeed)) %>% mutate(., R0 = "med")-> sumtable_shed_gamma_medR0

## hi R0
i <- 1
sumdata1 <- vector(mode='list')
for (var in c("hi")) {
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
                                                   meanSpeed=mean(epiSpeed)) %>% mutate(., R0 = "hi")-> sumtable_shed_gamma_hiR0


hivar_allR0_shed_gamma<- rbind(sumtable_shed_gamma_hiR0,sumtable_shed_gamma_lowR0,sumtable_shed_gamma_medR0)

shed_gamma_hivar<-ggplot(hivar_allR0_shed_gamma, aes(x=fracSDS, y=frac0))+
  geom_point(data=hivar_allR0_shed_gamma, aes(color=R0, shape=Corr), gamma = .85, size = 3.5)+
  scale_color_manual(values=c("darkblue", "darkgreen", "pink"))+
  labs(title="shed-gamma-hivar")+
  labs(x = "% SDS", y = "% 0 Reff")+
  theme_bw()+
  xlim(0, 0.02)+
  ylim(0, 0.8)
shed_gamma_hivar

## MED VAR
## low R0
## shed gamma
i <- 1
sumdata1 <- vector(mode='list')
for (var in c("med")) {
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
                                                   meanSpeed=mean(epiSpeed)) %>% mutate(., R0 = "low")-> sumtable_shed_gamma_lowR0

## med R0
i <- 1
sumdata1 <- vector(mode='list')
for (var in c("med")) {
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
                                                   meanSpeed=mean(epiSpeed)) %>% mutate(., R0 = "med")-> sumtable_shed_gamma_medR0

## hi R0
i <- 1
sumdata1 <- vector(mode='list')
for (var in c("med")) {
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
                                                   meanSpeed=mean(epiSpeed)) %>% mutate(., R0 = "hi")-> sumtable_shed_gamma_hiR0


medvar_allR0_shed_gamma<- rbind(sumtable_shed_gamma_hiR0,sumtable_shed_gamma_lowR0,sumtable_shed_gamma_medR0)

shed_gamma_medvar<-ggplot(medvar_allR0_shed_gamma, aes(x=fracSDS, y=frac0))+
  geom_point(data=medvar_allR0_shed_gamma, aes(color=R0, shape=Corr), gamma = .85, size = 3.5)+
  scale_color_manual(values=c("darkblue", "darkgreen", "pink"))+
  labs(title="shed-gamma-medvar")+
  labs(x = "% SDS", y = "% 0 Reff")+
  theme_bw() +
  xlim(0, 0.02)+
  ylim(0, 0.8)
shed_gamma_medvar

## LOW VAR
## low R0
## shed gamma
i <- 1
sumdata1 <- vector(mode='list')
for (var in c("low")) {
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
                                                   meanSpeed=mean(epiSpeed)) %>% mutate(., R0 = "low")-> sumtable_shed_gamma_lowR0

## med R0
i <- 1
sumdata1 <- vector(mode='list')
for (var in c("low")) {
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
                                                   meanSpeed=mean(epiSpeed)) %>% mutate(., R0 = "med")-> sumtable_shed_gamma_medR0

## hi R0
i <- 1
sumdata1 <- vector(mode='list')
for (var in c("low")) {
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
                                                   meanSpeed=mean(epiSpeed)) %>% mutate(., R0 = "hi")-> sumtable_shed_gamma_hiR0


lowvar_allR0_shed_gamma<- rbind(sumtable_shed_gamma_hiR0,sumtable_shed_gamma_lowR0,sumtable_shed_gamma_medR0)

shed_gamma_lowvar<-ggplot(lowvar_allR0_shed_gamma, aes(x=fracSDS, y=frac0))+
  geom_point(data=lowvar_allR0_shed_gamma, aes(color=R0, shape=Corr), gamma = .85, size = 3.5)+
  scale_color_manual(values=c("darkblue", "darkgreen", "pink"))+
  labs(title="shed-gamma-lowvar")+
  labs(x = "% SDS", y = "% 0 Reff")+
  theme_bw()+
  xlim(0, 0.02)+
  ylim(0, 0.8)
shed_gamma_lowvar

ggarrange(shed_gamma_lowvar,shed_gamma_medvar,shed_gamma_hivar, nrow = 1, ncol=3)


### Plots :) 
## super duper spreaders
shed_gamma_lowvar<-ggplot(lowvar_allR0_shed_gamma, aes(x=fracSDS, y=frac0))+
  geom_point(data=lowvar_allR0_shed_gamma, aes(color=R0, shape=Corr), alpha = .85, size = 3.5)+
  scale_color_manual(values=c("darkblue", "darkgreen", "pink"))+
  labs(title="shed-gamma-lowvar")+
  labs(x = "% SDS", y = "% 0 Reff")+
  theme_bw()+
  xlim(0, 0.02)+
  ylim(0, 0.8)
shed_gamma_lowvar

shed_gamma_medvar<-ggplot(medvar_allR0_shed_gamma, aes(x=fracSDS, y=frac0))+
  geom_point(data=medvar_allR0_shed_gamma, aes(color=R0, shape=Corr), alpha = .85, size = 3.5)+
  scale_color_manual(values=c("darkblue", "darkgreen", "pink"))+
  labs(title="shed-gamma-medvar")+
  labs(x = "% SDS", y = "% 0 Reff")+
  theme_bw() +
  xlim(0, 0.02)+
  ylim(0, 0.8)
shed_gamma_medvar

shed_gamma_hivar<-ggplot(hivar_allR0_shed_gamma, aes(x=fracSDS, y=frac0))+
  geom_point(data=hivar_allR0_shed_gamma, aes(color=R0, shape=Corr), alpha = .85, size = 3.5)+
  scale_color_manual(values=c("darkblue", "darkgreen", "pink"))+
  labs(title="shed-gamma-hivar")+
  labs(x = "% SDS", y = "% 0 Reff")+
  theme_bw()+
  xlim(0, 0.02)+
  ylim(0, 0.8)
shed_gamma_hivar

ggarrange(shed_gamma_lowvar,shed_gamma_medvar,shed_gamma_hivar, nrow = 1, ncol=3)

## super spreaders
shed_gamma_lowvar_SS<-ggplot(lowvar_allR0_shed_gamma, aes(x=fracSS, y=frac0))+
  geom_point(data=lowvar_allR0_shed_gamma, aes(color=R0, shape=Corr), alpha = .85, size = 3.5)+
  scale_color_manual(values=c("darkblue", "darkgreen", "pink"))+
  labs(title="shed-gamma-lowvar")+
  labs(x = "% SS", y = "% 0 Reff")+
  theme_bw()+
  ylim(0, 0.85)+
  xlim(0, 0.065)
shed_gamma_lowvar_SS

shed_gamma_medvar_SS<-ggplot(medvar_allR0_shed_gamma, aes(x=fracSS, y=frac0))+
  geom_point(data=medvar_allR0_shed_gamma, aes(color=R0, shape=Corr), alpha = .85, size = 3.5)+
  scale_color_manual(values=c("darkblue", "darkgreen", "pink"))+
  labs(title="shed-gamma-medvar")+
  labs(x = "% SS", y = "% 0 Reff")+
  theme_bw() +
  ylim(0, 0.85)+
  xlim(0, 0.065)
shed_gamma_medvar

shed_gamma_hivar_SS<-ggplot(hivar_allR0_shed_gamma, aes(x=fracSS, y=frac0))+
  geom_point(data=hivar_allR0_shed_gamma, aes(color=R0, shape=Corr), alpha = .85, size = 3.5)+
  scale_color_manual(values=c("darkblue", "darkgreen", "pink"))+
  labs(title="shed-gamma-hivar")+
  labs(x = "% SS", y = "% 0 Reff")+
  theme_bw()+
  ylim(0, 0.85)+
  xlim(0, 0.065)
shed_gamma_hivar

ggarrange(shed_gamma_lowvar_SS,shed_gamma_medvar_SS,shed_gamma_hivar_SS, nrow = 1, ncol=3)

rbind(rbind(lowvar_allR0_alpha_gamma, medvar_allR0_alpha_gamma),hivar_allR0_alpha_gamma) -> newdata_alpha_gamma
alpha_gamma_new<-ggplot(newdata, aes(x=Var, y=frac0))+
  geom_point(data=newdata_alpha_gamma, aes(color=R0, shape=Corr), alpha = .85, size = 3.5)+
  scale_color_manual(values=c("darkblue", "darkgreen", "pink"))+
  labs(title="alpha-gamma")+
  labs(x = "Var", y = "% 0 Reff")+
  theme_bw()+
  ylim(0, 1)

rbind(rbind(lowvar_allR0_c_gamma, medvar_allR0_c_gamma),hivar_allR0_c_gamma) -> newdata_c_gamma
c_gamma_new<-ggplot(newdata, aes(x=Var, y=frac0))+
  geom_point(data=newdata_c_gamma, aes(color=R0, shape=Corr), alpha = .85, size = 3.5)+
  scale_color_manual(values=c("darkblue", "darkgreen", "pink"))+
  labs(title="contact-gamma")+
  labs(x = "Var", y = "% 0 Reff")+
  theme_bw()+
  ylim(0, 1)

rbind(rbind(lowvar_allR0_c_alpha, medvar_allR0_c_alpha),hivar_allR0_c_alpha) -> newdata_c_alpha
c_alpha_new<-ggplot(newdata_c_alpha, aes(x=Var, y=frac0))+
  geom_point(data=newdata_c_gamma, aes(color=R0, shape=Corr), alpha = .85, size = 3.5)+
  scale_color_manual(values=c("darkblue", "darkgreen", "pink"))+
  labs(title="contact-alpha")+
  labs(x = "Var", y = "% 0 Reff")+
  theme_bw()+
  ylim(0, 1)

rbind(rbind(lowvar_allR0_shed_gamma, medvar_allR0_shed_gamma),hivar_allR0_shed_gamma) -> newdata_shed_gamma
shed_gamma_new<-ggplot(newdata, aes(x=Var, y=frac0))+
  geom_point(data=newdata_shed_gamma, aes(color=R0, shape=Corr), alpha = .85, size = 3.5)+
  scale_color_manual(values=c("darkblue", "darkgreen", "pink"))+
  labs(title="shed-gamma")+
  labs(x = "Var", y = "% 0 Reff")+
  theme_bw()+
  ylim(0, 1)

rbind(rbind(lowvar_allR0_shed_alpha, medvar_allR0_shed_alpha),hivar_allR0_shed_alpha) -> newdata_shed_alpha
shed_alpha_new<-ggplot(newdata, aes(x=Var, y=frac0))+
  geom_point(data=newdata_shed_alpha, aes(color=R0, shape=Corr), alpha = .85, size = 3.5)+
  scale_color_manual(values=c("darkblue", "darkgreen", "pink"))+
  labs(title="shed-alpha")+
  labs(x = "Var", y = "% 0 Reff")+
  theme_bw()+
  ylim(0, 1)

rbind(rbind(lowvar_allR0_c_shed, medvar_allR0_c_shed),hivar_allR0_c_shed) -> newdata_c_shed
c_shed_new<-ggplot(newdata, aes(x=Var, y=frac0))+
  geom_point(data=newdata_c_shed, aes(color=R0, shape=Corr), alpha = .85, size = 3.5)+
  scale_color_manual(values=c("darkblue", "darkgreen", "pink"))+
  labs(title="contact-shed")+
  labs(x = "Var", y = "% 0 Reff")+
  theme_bw()+
  ylim(0, 1)

ggarrange(c_shed_new, alpha_gamma_new, c_alpha_new, c_gamma_new, shed_alpha_new, shed_gamma_new, nrow=3, ncol=2,common.legend = TRUE, legend="right" )

