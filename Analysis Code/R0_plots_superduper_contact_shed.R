### contact shed

## low R0
## hi var
i <- 1
sumdata1 <- vector(mode='list')
for (var in c("hi")) {
  for (corr in c("no","neg","pos")) {
    data <- readRDS(paste0("out_",corr,"corr_",var,"var_c_shed.RDS"))
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
                                                   meanSpeed=mean(epiSpeed)) %>% mutate(., R0 = "low")-> sumtable_c_shed_lowR0

## med R0
i <- 1
sumdata1 <- vector(mode='list')
for (var in c("hi")) {
  for (corr in c("no","neg","pos")) {
    data <- readRDS(paste0("out_",corr,"corr_",var,"var_c_shed.RDS"))
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
                                                   meanSpeed=mean(epiSpeed)) %>% mutate(., R0 = "med")-> sumtable_c_shed_medR0

## hi R0
i <- 1
sumdata1 <- vector(mode='list')
for (var in c("hi")) {
  for (corr in c("no","neg","pos")) {
    data <- readRDS(paste0("out_",corr,"corr_",var,"var_c_shed.RDS"))
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
                                                   meanSpeed=mean(epiSpeed)) %>% mutate(., R0 = "hi")-> sumtable_c_shed_hiR0


hivar_allR0_c_shed<- rbind(sumtable_c_shed_hiR0,sumtable_c_shed_lowR0,sumtable_c_shed_medR0)

c_shed_hivar<-ggplot(hivar_allR0_c_shed, aes(x=fracSDS, y=frac0))+
  geom_point(data=hivar_allR0_c_shed, aes(color=R0, shape=Corr), alpha = .85, size = 3.5)+
  scale_color_manual(values=c("darkblue", "darkgreen", "pink"))+
  labs(title="c-shed-hivar")+
  labs(x = "% SDS", y = "% 0 Reff")+
  theme_bw()+
  xlim(0, 0.02)+
  ylim(0, 0.9)
c_shed_hivar

c_shed_hivar_SS<-ggplot(hivar_allR0_c_shed, aes(x=fracSS, y=frac0))+
  geom_point(data=hivar_allR0_c_shed, aes(color=R0, shape=Corr), alpha = .85, size = 3.5)+
  scale_color_manual(values=c("darkblue", "darkgreen", "pink"))+
  labs(title="c-shed-hivar")+
  labs(x = "% SS", y = "% 0 Reff")+
  theme_bw()+
  ylim(0, 0.9)+
  xlim(0, 0.07)
c_shed_hivar_SS

## MED VAR
## low R0
## c shed
i <- 1
sumdata1 <- vector(mode='list')
for (var in c("med")) {
  for (corr in c("no","neg","pos")) {
    data <- readRDS(paste0("out_",corr,"corr_",var,"var_c_shed.RDS"))
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
                                                   meanSpeed=mean(epiSpeed)) %>% mutate(., R0 = "low")-> sumtable_c_shed_lowR0

## med R0
i <- 1
sumdata1 <- vector(mode='list')
for (var in c("med")) {
  for (corr in c("no","neg","pos")) {
    data <- readRDS(paste0("out_",corr,"corr_",var,"var_c_shed.RDS"))
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
                                                   meanSpeed=mean(epiSpeed)) %>% mutate(., R0 = "med")-> sumtable_c_shed_medR0

## hi R0
i <- 1
sumdata1 <- vector(mode='list')
for (var in c("med")) {
  for (corr in c("no","neg","pos")) {
    data <- readRDS(paste0("out_",corr,"corr_",var,"var_c_shed.RDS"))
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
                                                   meanSpeed=mean(epiSpeed)) %>% mutate(., R0 = "hi")-> sumtable_c_shed_hiR0


medvar_allR0_c_shed<- rbind(sumtable_c_shed_hiR0,sumtable_c_shed_lowR0,sumtable_c_shed_medR0)

c_shed_medvar<-ggplot(medvar_allR0_c_shed, aes(x=fracSDS, y=frac0))+
  geom_point(data=medvar_allR0_c_shed, aes(color=R0, shape=Corr), alpha = .85, size = 3.5)+
  scale_color_manual(values=c("darkblue", "darkgreen", "pink"))+
  labs(title="c-shed-medvar")+
  labs(x = "% SDS", y = "% 0 Reff")+
  theme_bw() +
  xlim(0, 0.02)+
  ylim(0, 0.9)
c_shed_medvar

c_shed_medvar_SS<-ggplot(medvar_allR0_c_shed, aes(x=fracSS, y=frac0))+
  geom_point(data=medvar_allR0_c_shed, aes(color=R0, shape=Corr), alpha = .85, size = 3.5)+
  scale_color_manual(values=c("darkblue", "darkgreen", "pink"))+
  labs(title="c-shed-medvar")+
  labs(x = "% SS", y = "% 0 Reff")+
  theme_bw()+
  ylim(0, 0.9)+
  xlim(0, 0.07)
c_shed_medvar_SS

## LOW VAR
## low R0
## c shed
i <- 1
sumdata1 <- vector(mode='list')
for (var in c("low")) {
  for (corr in c("no","neg","pos")) {
    data <- readRDS(paste0("out_",corr,"corr_",var,"var_c_shed.RDS"))
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
                                                   meanSpeed=mean(epiSpeed)) %>% mutate(., R0 = "low")-> sumtable_c_shed_lowR0

## med R0
i <- 1
sumdata1 <- vector(mode='list')
for (var in c("low")) {
  for (corr in c("no","neg","pos")) {
    data <- readRDS(paste0("out_",corr,"corr_",var,"var_c_shed.RDS"))
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
                                                   meanSpeed=mean(epiSpeed)) %>% mutate(., R0 = "med")-> sumtable_c_shed_medR0

## hi R0
i <- 1
sumdata1 <- vector(mode='list')
for (var in c("low")) {
  for (corr in c("no","neg","pos")) {
    data <- readRDS(paste0("out_",corr,"corr_",var,"var_c_shed.RDS"))
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
                                                   meanSpeed=mean(epiSpeed)) %>% mutate(., R0 = "hi")-> sumtable_c_shed_hiR0


lowvar_allR0_c_shed<- rbind(sumtable_c_shed_hiR0,sumtable_c_shed_lowR0,sumtable_c_shed_medR0)

c_shed_lowvar<-ggplot(lowvar_allR0_c_shed, aes(x=fracSDS, y=frac0))+
  geom_point(data=lowvar_allR0_c_shed, aes(color=R0, shape=Corr), alpha = .85, size = 3.5)+
  scale_color_manual(values=c("darkblue", "darkgreen", "pink"))+
  labs(title="c-shed-lowvar")+
  labs(x = "% SDS", y = "% 0 Reff")+
  theme_bw()+
  xlim(0, 0.02)+
  ylim(0, 0.9)
c_shed_lowvar

c_shed_lowvar_SS<-ggplot(lowvar_allR0_c_shed, aes(x=fracSS, y=frac0))+
  geom_point(data=lowvar_allR0_c_shed, aes(color=R0, shape=Corr), alpha = .85, size = 3.5)+
  scale_color_manual(values=c("darkblue", "darkgreen", "pink"))+
  labs(title="c-shed-lowvar")+
  labs(x = "% SS", y = "% 0 Reff")+
  theme_bw()+
  ylim(0, 0.9)+
  xlim(0, 0.07)
c_shed_lowvar_SS

ggarrange(c_shed_lowvar,c_shed_medvar,c_shed_hivar, nrow = 1, ncol=3)

ggarrange(c_shed_lowvar_SS,c_shed_medvar_SS,c_shed_hivar_SS, nrow = 1, ncol=3)
