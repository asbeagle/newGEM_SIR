library(tidyverse)
library(magrittr)
library(MASS)

## Empirical cumulative density function plots
## (Technically, 1-ECDF!)
## R0 = 1 plots
var1=c('c',    'c',    'shed', 'shed', 'c',   'alpha')
var2=c('gamma','alpha','gamma','alpha','shed','gamma')
i <- 1
data <- vector(mode='list', length=18)
for (j in 1:6) { ## loop over the six different covariance combinations
  for (covMatrix in c("nocorr", "negcorr", "poscorr")) {
    z1 = readRDS(paste0(paste("out_R=1",covMatrix,"hivar",var1[j],var2[j],sep="_"),".RDS"))
    z2 = readRDS(paste0(paste("out_R=1",covMatrix,"medvar",var1[j],var2[j],sep="_"),".RDS"))
    z3 = readRDS(paste0(paste("out_R=1",covMatrix,"lowvar",var1[j],var2[j],sep="_"),".RDS"))

    ## Provide a label for each plot panel
    panel <- paste0(LETTERS[j],switch(covMatrix,poscorr="(3)",negcorr="(1)",nocorr="(2)"))
    
    ## Give nicer names
    covMatrix <- switch(covMatrix,poscorr="(+) Cov",negcorr="(-) Cov",nocorr="(0) Cov")
    trait1 <- switch(var1[j],c="Contact",shed="Infectiousness",alpha="Virulence")
    trait2 <- switch(var2[j],shed="Infectiousness",alpha="Virulence",gamma="Recovery")
    
    z <- vector(mode='list',length=3)
    lapply(1:length(z1), function(i) data.frame(peakSize=max(z1[[i]][[1]]$I,na.rm=TRUE), ## peak epidemic size
                                       peakPrev=max(z1[[i]][[1]]$I/(z1[[i]][[1]]$S+z1[[i]][[1]]$I+z1[[i]][[1]]$R),na.rm=TRUE),
                                       disp=ifelse(max(z1[[i]][[1]]$t)>99, glm.nb(z1[[i]][[2]]$numInf~1)$theta, NA), ## dispersion parameter of negative binomial distribution fit to numInf
                                       fadeout=ifelse(length(z1[[i]][[1]]$t)<99, 1, 0), # binary: did this replicate go extinct?
                                       rep=i)) %>% 
      do.call("rbind.data.frame",.) %>%
      mutate(., Variance="high",cov=covMatrix, traits=paste(trait1,trait2,sep="-"), panelLab=panel) -> z[[1]]
    lapply(1:length(z2), function(i) data.frame(peakSize=max(z2[[i]][[1]]$I,na.rm=TRUE), ## peak epidemic size
                                       peakPrev=max(z2[[i]][[1]]$I/(z2[[i]][[1]]$S+z2[[i]][[1]]$I+z2[[i]][[1]]$R),na.rm=TRUE),
                                       disp=ifelse(max(z2[[i]][[1]]$t)>99, glm.nb(z2[[i]][[2]]$numInf~1)$theta, NA), ## dispersion parameter of negative binomial distribution fit to numInf
                                       fadeout=ifelse(length(z2[[i]][[1]]$t)<99, 1, 0),
                                       rep=i)) %>% 
      do.call("rbind.data.frame",.) %>%
      mutate(., Variance="med",cov=covMatrix, traits=paste(trait1,trait2,sep="-"), panelLab=panel) -> z[[2]]
    lapply(1:length(z3), function(i) data.frame(peakSize=max(z3[[i]][[1]]$I,na.rm=TRUE), ## peak epidemic size
                                       peakPrev=max(z3[[i]][[1]]$I/(z3[[i]][[1]]$S+z3[[i]][[1]]$I+z3[[i]][[1]]$R),na.rm=TRUE),
                                       disp=ifelse(max(z3[[i]][[1]]$t)>99, glm.nb(z3[[i]][[2]]$numInf~1)$theta, NA), ## dispersion parameter of negative binomial distribution fit to numInf
                                       fadeout=ifelse(length(z3[[i]][[1]]$t)<99, 1, 0),
                                       rep=i)) %>% 
      do.call("rbind.data.frame",.) %>%
      mutate(., Variance="low",cov=covMatrix, traits=paste(trait1,trait2,sep="-"), panelLab=panel) -> z[[3]]
    data[[i]] <- do.call("rbind.data.frame",z)
    i <- i+1
  }
}
data <- do.call("rbind.data.frame",data)
data$cov <- factor(data$cov, levels=c("(-) Cov","(0) Cov","(+) Cov"))
data$Variance <- factor(data$Variance, levels=c("low","med","high"))
data$traits <- factor(data$traits, levels=c("Contact-Recovery","Contact-Virulence","Infectiousness-Recovery","Infectiousness-Virulence","Contact-Infectiousness","Virulence-Recovery"))
# compute the fraction of replicates that make it to each epidemic size
data %>% group_by(Variance, cov, traits) %>% summarise(size=seq(0,max(data$peakSize)+2), ECDF=sapply(seq(0,max(data$peakSize)+2), function(i) sum(peakSize>=i)/100)) -> data2
# compute the fraction of replicates that make it to each epidemic prevalence
data %>% group_by(Variance, cov, traits) %>% summarise(size=seq(0,max(data$peakPrev)+0.02,0.01), ECDF=sapply(seq(0,max(data$peakPrev)+0.02,0.01), function(i) sum(peakPrev>=i)/100)) -> data3

png(filename="./All_ECDF_R0=1.png", height=6.5, width=9, units='in', res=300)
ggplot(data2, aes(x=size, y=ECDF, group=Variance, colour=Variance)) + 
  geom_line(size=0.8) + 
  facet_grid(cov~traits) +
  ylim(0,1.05) + 
  geom_text(data=subset(data %>% group_by(Variance, cov, traits) %>% summarise(pan=unique(panelLab)),Variance=='high'),
            mapping=aes(x=2, y=1.03, label=pan),color='black', size=2.5, hjust=0) +  
  geom_text(data=subset(data %>% 
                           group_by(Variance, cov, traits) %>% 
                           summarise(fadeoutFrac=sum(fadeout==1)/100) %>% 
                           mutate(lab=paste0("Pr","=",fadeoutFrac)),
                         Variance=="low"),
             mapping=aes(x=120, y=1.03, label=lab, color=Variance), size=2.5, hjust=0) +
  geom_text(data=subset(data %>% 
                          group_by(Variance, cov, traits) %>% 
                          summarise(fadeoutFrac=sum(fadeout==1)/100) %>% 
                          mutate(lab=paste0("Pr","=",fadeoutFrac)),
                        Variance=="med"),
            mapping=aes(x=120, y=0.93, label=lab, color=Variance), size=2.5, hjust=0) +
  geom_text(data=subset(data %>% 
                          group_by(Variance, cov, traits) %>% 
                          summarise(fadeoutFrac=sum(fadeout==1)/100) %>% 
                          mutate(lab=paste0("Pr","=",fadeoutFrac)),
                        Variance=="high"),
            mapping=aes(x=120, y=0.83, label=lab, color=Variance), size=2.5, hjust=0) +
  xlab("Epidemic size") + 
  ylab(expression("%"~of~simulations~(italic(R)[0]==1))) +           
  ggtitle("Intergroup trait pairings                                                                                                                             Intragroup trait pairings") + 
  theme_bw() + 
  theme(axis.text=element_text(size=6), axis.title=element_text(size=8), legend.text=element_text(size=7), legend.title=element_text(size=8), plot.title = element_text(size=8), strip.text.x = element_text(size = 7), strip.text.y = element_text(size = 7)) +
  scale_color_manual(values=c("cornflowerblue", "pink", "sienna")) +
  scale_fill_manual(values=c("cornflowerblue", "pink", "sienna")) 
dev.off()

png(filename="./All_ECDF_peak_prev_R0=1.png", height=6.5, width=9, units='in', res=300)
ggplot(data3, aes(x=size, y=ECDF, group=Variance, colour=Variance)) + 
  geom_line(size=0.8) + 
  facet_grid(cov~traits) +
  ylim(0,1.05) + 
  geom_text(data=subset(data %>% group_by(Variance, cov, traits) %>% summarise(pan=unique(panelLab)),Variance=='high'),
            mapping=aes(x=0.01, y=1.03, label=pan),color='black', size=2.5, hjust=0) +  
  geom_text(data=subset(data %>% 
                          group_by(Variance, cov, traits) %>% 
                          summarise(fadeoutFrac=sum(fadeout==1)/100) %>% 
                          mutate(lab=paste0("Pr","=",fadeoutFrac)),
                        Variance=="low"),
            mapping=aes(x=0.52, y=1.03, label=lab, color=Variance), size=2.5, hjust=0) +
  geom_text(data=subset(data %>% 
                          group_by(Variance, cov, traits) %>% 
                          summarise(fadeoutFrac=sum(fadeout==1)/100) %>% 
                          mutate(lab=paste0("Pr","=",fadeoutFrac)),
                        Variance=="med"),
            mapping=aes(x=0.52, y=0.93, label=lab, color=Variance), size=2.5, hjust=0) +
  geom_text(data=subset(data %>% 
                          group_by(Variance, cov, traits) %>% 
                          summarise(fadeoutFrac=sum(fadeout==1)/100) %>% 
                          mutate(lab=paste0("Pr","=",fadeoutFrac)),
                        Variance=="high"),
            mapping=aes(x=0.52, y=0.83, label=lab, color=Variance), size=2.5, hjust=0) +
  xlab("Epidemic prevalence") + 
  ylab(expression("%"~of~simulations~(italic(R)[0]==1))) + 
  ggtitle("Intergroup trait pairings                                                                                                                             Intragroup trait pairings") + 
  theme_bw() + 
  theme(axis.text=element_text(size=6), axis.title=element_text(size=8), legend.text=element_text(size=7), legend.title=element_text(size=8), plot.title = element_text(size=8), strip.text.x = element_text(size = 7), strip.text.y = element_text(size = 7)) +
  scale_color_manual(values=c("cornflowerblue", "pink", "sienna")) +
  scale_fill_manual(values=c("cornflowerblue", "pink", "sienna")) 
dev.off()


png(filename="./All_fadeout_peaks_R0=1.png", height=6.5, width=9, units='in', res=300)
ggplot(subset(data, fadeout==1), aes(x=peakSize)) + 
  geom_histogram(aes(group=Variance,color=Variance,fill=Variance),bins=10,alpha=0.2,position='identity')+
  facet_grid(cov~traits) +
  ylim(0,105) + 
  geom_text(data=subset(data %>% group_by(Variance, cov, traits) %>% summarise(pan=unique(panelLab)),Variance=='high'),
            mapping=aes(x=-10, y=103, label=pan),color='black', size=2.5, hjust=0) +  
  xlab("Peak") + 
  ylab(expression("#"~of~simulations~(italic(R)[0]==1))) + 
  ylab("# of simulations") +
  ggtitle("Intergroup trait pairings                                                                                                                               Intragroup trait pairings") + 
  theme_bw() + 
  theme(axis.text=element_text(size=6), axis.title=element_text(size=8), legend.text=element_text(size=7), legend.title=element_text(size=8), plot.title = element_text(size=8), strip.text.x = element_text(size = 7), strip.text.y = element_text(size = 7)) +
  scale_color_manual(values=c("cornflowerblue", "pink", "sienna")) +
  scale_fill_manual(values=c("cornflowerblue", "pink", "sienna")) 
dev.off()

## Very low fraction of fadeouts that get very large - can actually compute this
## What fraction of simulations that fadeout every get to a prevalence above
data %>% filter(fadeout==1) %>% # look only at epidemics that faded out
  group_by(Variance, cov, traits) %>% # summarize across replicates
  summarize(prop1=sum(peakPrev > 0.2)/100,
            prop2=sum(peakSize > 50)/100) -> fadeout_data # compute the fraction of replicates with a prevalence above any particular percentage

fadeout_data %>% print(n=54) # show all the results



## R0 = 4 plots
var1=c('c',    'c',    'shed', 'shed', 'c',   'alpha')
var2=c('gamma','alpha','gamma','alpha','shed','gamma')
i <- 1
data <- vector(mode='list', length=18)
for (j in 1:6) { ## loop over the six different covariance combinations
  for (covMatrix in c("nocorr", "negcorr", "poscorr")) {
    z1 = readRDS(paste0(paste("out_R=4",covMatrix,"hivar",var1[j],var2[j],sep="_"),".RDS"))
    z2 = readRDS(paste0(paste("out_R=4",covMatrix,"medvar",var1[j],var2[j],sep="_"),".RDS"))
    z3 = readRDS(paste0(paste("out_R=4",covMatrix,"lowvar",var1[j],var2[j],sep="_"),".RDS"))
    
    ## Provide a label for each plot panel
    panel <- paste0(LETTERS[j],switch(covMatrix,poscorr="(3)",negcorr="(1)",nocorr="(2)"))
    
    ## Give nicer names
    covMatrix <- switch(covMatrix,poscorr="(+) Cov",negcorr="(-) Cov",nocorr="(0) Cov")
    trait1 <- switch(var1[j],c="Contact",shed="Infectiousness",alpha="Virulence")
    trait2 <- switch(var2[j],shed="Infectiousness",alpha="Virulence",gamma="Recovery")
    
    z <- vector(mode='list',length=3)
    lapply(1:length(z1), function(i) data.frame(peakSize=max(z1[[i]][[1]]$I,na.rm=TRUE), ## peak epidemic size
                                                peakPrev=max(z1[[i]][[1]]$I/(z1[[i]][[1]]$S+z1[[i]][[1]]$I+z1[[i]][[1]]$R),na.rm=TRUE),
                                                fadeout=ifelse(length(z1[[i]][[1]]$t)<99, 1, 0)) %>% # binary: did this replicate go extinct?
                                                mutate(.,
                                                       disp=ifelse(fadeout==0, ifelse(inherits(try(glm.nb(z1[[i]][[2]]$numInf~1)$theta),'try-error'), NA, glm.nb(z1[[i]][[2]]$numInf~1)$theta), NA), ## dispersion parameter of negative binomial distribution fit to numInf
                                                       dispAlt=ifelse(peakSize > 10, ifelse(inherits(try(glm.nb(z1[[i]][[2]]$numInf~1)$theta),'try-error'), NA, glm.nb(z1[[i]][[2]]$numInf~1)$theta), NA), ## dispersion parameter of negative binomial distribution fit to numInf
                                                       rep=i)) %>% 
      do.call("rbind.data.frame",.) %>%
      mutate(., Variance="high",cov=covMatrix, traits=paste(trait1,trait2,sep="-"), panelLab=panel) -> z[[1]]
    lapply(1:length(z2), function(i) data.frame(peakSize=max(z2[[i]][[1]]$I,na.rm=TRUE), ## peak epidemic size
                                                peakPrev=max(z2[[i]][[1]]$I/(z2[[i]][[1]]$S+z2[[i]][[1]]$I+z2[[i]][[1]]$R),na.rm=TRUE),
                                                fadeout=ifelse(length(z2[[i]][[1]]$t)<99, 1, 0)) %>% # binary: did this replicate go extinct?
                                                mutate(., 
                                                       disp=ifelse(fadeout==0, ifelse(inherits(try(glm.nb(z2[[i]][[2]]$numInf~1)$theta),'try-error'), NA, glm.nb(z2[[i]][[2]]$numInf~1)$theta), NA), ## dispersion parameter of negative binomial distribution fit to numInf
                                                       dispAlt=ifelse(peakSize > 10, ifelse(inherits(try(glm.nb(z2[[i]][[2]]$numInf~1)$theta),'try-error'), NA, glm.nb(z2[[i]][[2]]$numInf~1)$theta), NA), ## dispersion parameter of negative binomial distribution fit to numInf
                                                       rep=i)) %>% 
      do.call("rbind.data.frame",.) %>%
      mutate(., Variance="med",cov=covMatrix, traits=paste(trait1,trait2,sep="-"), panelLab=panel) -> z[[2]]
    lapply(1:length(z3), function(i) data.frame(peakSize=max(z3[[i]][[1]]$I,na.rm=TRUE), ## peak epidemic size
                                                peakPrev=max(z3[[i]][[1]]$I/(z3[[i]][[1]]$S+z3[[i]][[1]]$I+z3[[i]][[1]]$R),na.rm=TRUE),
                                                fadeout=ifelse(length(z3[[i]][[1]]$t)<99, 1, 0)) %>% # binary: did this replicate go extinct?
                                                mutate(., 
                                                       disp=ifelse(fadeout==0, ifelse(inherits(try(glm.nb(z3[[i]][[2]]$numInf~1)$theta),'try-error'), NA, glm.nb(z3[[i]][[2]]$numInf~1)$theta), NA), ## dispersion parameter of negative binomial distribution fit to numInf
                                                       dispAlt=ifelse(peakSize > 10, ifelse(inherits(try(glm.nb(z3[[i]][[2]]$numInf~1)$theta),'try-error'), NA, glm.nb(z3[[i]][[2]]$numInf~1)$theta), NA), ## dispersion parameter of negative binomial distribution fit to numInf
                                                       rep=i)) %>% 
      do.call("rbind.data.frame",.) %>%
      mutate(., Variance="low",cov=covMatrix, traits=paste(trait1,trait2,sep="-"), panelLab=panel) -> z[[3]]
    data[[i]] <- do.call("rbind.data.frame",z)
    i <- i+1
  }
}
data <- do.call("rbind.data.frame",data)
data$cov <- factor(data$cov, levels=c("(-) Cov","(0) Cov","(+) Cov"))
data$Variance <- factor(data$Variance, levels=c("low","med","high"))
data$traits <- factor(data$traits, levels=c("Contact-Recovery","Contact-Virulence","Infectiousness-Recovery","Infectiousness-Virulence","Contact-Infectiousness","Virulence-Recovery"))
# compute the fraction of replicates that make it to each epidemic size
data %>% group_by(Variance, cov, traits) %>% summarise(size=seq(0,max(data$peakSize)+2), ECDF=sapply(seq(0,max(data$peakSize)+2), function(i) sum(peakSize>=i)/100)) -> data2
# compute the fraction of replicates that make it to each epidemic prevalence
data %>% group_by(Variance, cov, traits) %>% summarise(size=seq(0,max(data$peakPrev)+0.02,0.01), ECDF=sapply(seq(0,max(data$peakPrev)+0.02,0.01), function(i) sum(peakPrev>=i)/100)) -> data3


subset(data, traits=="Contact-Recovery" & cov=="(-) Cov" & Variance=="high")$disp
subset(data, traits=="Contact-Recovery" & cov=="(-) Cov" & Variance=="med")$disp

data %>% group_by(traits, cov, Variance) %>% summarise(meanDisp=mean(dispAlt,na.rm=T), meanPeak=mean(peakSize), fadeoutFrac=sum(fadeout==1)/100) %>% 
  ggplot(., aes(x=meanDisp, y=fadeoutFrac, color=Variance, fill=Variance)) + 
  geom_point() + 
  facet_grid(cov~traits)

png(filename="./All_ECDF_R0=4.png", height=6.5, width=9, units='in', res=300)
ggplot(data2, aes(x=size, y=ECDF, group=Variance, colour=Variance)) + 
  geom_line(size=0.8) + 
  facet_grid(cov~traits) +
  ylim(0,1.05) + 
  geom_text(data=subset(data %>% group_by(Variance, cov, traits) %>% summarise(pan=unique(panelLab)),Variance=='high'),
            mapping=aes(x=2, y=1.03, label=pan),color='black', size=2.5, hjust=0) +  
  geom_text(data=subset(data %>% 
                          group_by(Variance, cov, traits) %>% 
                          summarise(fadeoutFrac=sum(fadeout==1)/100) %>% 
                          mutate(lab=paste0("Pr","=",fadeoutFrac)),
                        Variance=="low"),
            mapping=aes(x=140, y=1.03, label=lab, color=Variance), size=2.5, hjust=0) +
  geom_text(data=subset(data %>% 
                          group_by(Variance, cov, traits) %>% 
                          summarise(fadeoutFrac=sum(fadeout==1)/100) %>% 
                          mutate(lab=paste0("Pr","=",fadeoutFrac)),
                        Variance=="med"),
            mapping=aes(x=140, y=0.93, label=lab, color=Variance), size=2.5, hjust=0) +
  geom_text(data=subset(data %>% 
                          group_by(Variance, cov, traits) %>% 
                          summarise(fadeoutFrac=sum(fadeout==1)/100) %>% 
                          mutate(lab=paste0("Pr","=",fadeoutFrac)),
                        Variance=="high"),
            mapping=aes(x=140, y=0.83, label=lab, color=Variance), size=2.5, hjust=0) +
  xlab("Epidemic size") + 
  ylab(expression("%"~of~simulations~(italic(R)[0]==4))) + 
  ggtitle("Intergroup trait pairings                                                                                                                             Intragroup trait pairings") + 
  theme_bw() + 
  theme(axis.text=element_text(size=6), axis.title=element_text(size=8), legend.text=element_text(size=7), legend.title=element_text(size=8), plot.title = element_text(size=8), strip.text.x = element_text(size = 7), strip.text.y = element_text(size = 7)) +
  scale_color_manual(values=c("cornflowerblue", "pink", "sienna")) +
  scale_fill_manual(values=c("cornflowerblue", "pink", "sienna")) 
dev.off()

png(filename="./All_fadeout_peaks_R0=4.png", height=6.5, width=9, units='in', res=300)
ggplot(subset(data, fadeout==1), aes(x=peakSize)) + 
  geom_histogram(aes(group=Variance,color=Variance,fill=Variance),bins=10,alpha=0.3,position='identity')+
  facet_grid(cov~traits) +
  ylim(0,95) + 
  geom_text(data=subset(data %>% group_by(Variance, cov, traits) %>% summarise(pan=unique(panelLab)),Variance=='high'),
            mapping=aes(x=-10, y=93, label=pan),color='black', size=2.5, hjust=0) +  
  xlab("Peak") + 
  ylab(expression("#"~of~simulations~(italic(R)[0]==4))) + 
  ggtitle("Intergroup trait pairings                                                                                                                               Intragroup trait pairings") + 
  theme_bw() + 
  theme(axis.text=element_text(size=6), axis.title=element_text(size=8), legend.text=element_text(size=7), legend.title=element_text(size=8), plot.title = element_text(size=8), strip.text.x = element_text(size = 7), strip.text.y = element_text(size = 7)) +
  scale_color_manual(values=c("pink", "sienna")) +
  scale_fill_manual(values=c("pink", "sienna")) 
dev.off()

png(filename="./All_dispersion_values_R0=4.png", height=6.5, width=9, units='in', res=300)
ggplot(data, aes(x=dispAlt)) + 
  geom_density(aes(group=Variance,color=Variance,fill=Variance),alpha=0.3,position='identity')+
  geom_vline(xintercept = 0.16, linetype="dashed") +
  facet_grid(cov~traits) +
  ylim(0,63) + xlim(0,1.2) +
  geom_text(data=subset(data %>% group_by(Variance, cov, traits) %>% summarise(pan=unique(panelLab)),Variance=='high'),
            mapping=aes(x=0.06, y=61, label=pan),color='black', size=2.5, hjust=0) +  
  xlab(expression(Dispersion~(italic(k))))+
  ylab(expression(Density~(italic(R)[0]==4))) + 
  ggtitle("Intergroup trait pairings                                                                                                                               Intragroup trait pairings") + 
  theme_bw() + 
  theme(axis.text=element_text(size=6), axis.title=element_text(size=8), legend.text=element_text(size=7), legend.title=element_text(size=8), plot.title = element_text(size=8), strip.text.x = element_text(size = 7), strip.text.y = element_text(size = 7)) +
  scale_color_manual(values=c("cornflowerblue", "pink", "sienna")) +
  scale_fill_manual(values=c("cornflowerblue", "pink", "sienna")) 
dev.off()

png(filename="./All_peak_against_dispersion_R0=4.png", height=4, width=9, units='in', res=300)
ggplot(mutate(data, Covariance=cov), aes(x=dispAlt, y=peakSize, shape=Variance, color=Covariance)) + 
  geom_point(size=0.8) + 
  geom_vline(xintercept = 0.16, linetype="dashed") +
  facet_grid(~traits) +
  #ylim(80,205) + 
  geom_text(data=subset(mutate(data, Covariance=cov) %>% group_by(Variance, Covariance, traits) %>% summarise(pan=unique(strsplit(panelLab,"\\(")[[1]][1])),Variance=='high' & Covariance=="(0) Cov"),
            mapping=aes(x=0.075, y=203, label=pan),color='black', size=2.5, hjust=0) +  
  xlab(expression(Dispersion~(italic(k)))) +
  ylab(expression(Peak~size~(italic(R)[0]==4))) + 
  ggtitle("Intergroup trait pairings                                                                                                                               Intragroup trait pairings") + 
  theme_bw() + 
  theme(axis.text=element_text(size=6), axis.title=element_text(size=8), legend.text=element_text(size=7), legend.title=element_text(size=8), plot.title = element_text(size=8), strip.text.x = element_text(size = 7), strip.text.y = element_text(size = 7)) +
  scale_color_manual(values=c("cornflowerblue", "pink", "sienna")) +
  scale_fill_manual(values=c("cornflowerblue", "pink", "sienna")) 
dev.off()



## R0 = 8 plots
var1=c('c',    'c',    'shed', 'shed', 'c',   'alpha')
var2=c('gamma','alpha','gamma','alpha','shed','gamma')
i <- 1
data <- vector(mode='list', length=18)
for (j in 1:6) { ## loop over the six different covariance combinations
  for (covMatrix in c("nocorr", "negcorr", "poscorr")) {
    z1 = readRDS(paste0(paste("out_R=8",covMatrix,"hivar",var1[j],var2[j],sep="_"),".RDS"))
    z2 = readRDS(paste0(paste("out_R=8",covMatrix,"medvar",var1[j],var2[j],sep="_"),".RDS"))
    z3 = readRDS(paste0(paste("out_R=8",covMatrix,"lowvar",var1[j],var2[j],sep="_"),".RDS"))
    
    ## Provide a label for each plot panel
    panel <- paste0(LETTERS[j],switch(covMatrix,poscorr="(3)",negcorr="(1)",nocorr="(2)"))
    
    ## Give nicer names
    covMatrix <- switch(covMatrix,poscorr="(+) Cov",negcorr="(-) Cov",nocorr="(0) Cov")
    trait1 <- switch(var1[j],c="Contact",shed="Infectiousness",alpha="Virulence")
    trait2 <- switch(var2[j],shed="Infectiousness",alpha="Virulence",gamma="Recovery")
    
    z <- vector(mode='list',length=3)
    lapply(1:length(z1), function(i) data.frame(peakSize=max(z1[[i]][[1]]$I,na.rm=TRUE), ## peak epidemic size
                                                peakPrev=max(z1[[i]][[1]]$I/(z1[[i]][[1]]$S+z1[[i]][[1]]$I+z1[[i]][[1]]$R),na.rm=TRUE),
                                                disp=ifelse(max(z1[[i]][[1]]$t)>99, glm.nb(z1[[i]][[2]]$numInf~1)$theta, NA), ## dispersion parameter of negative binomial distribution fit to numInf
                                                fadeout=ifelse(length(z1[[i]][[1]]$t)<99, 1, 0), # binary: did this replicate go extinct?
                                                rep=i)) %>% 
      do.call("rbind.data.frame",.) %>%
      mutate(., Variance="high",cov=covMatrix, traits=paste(trait1,trait2,sep="-"), panelLab=panel) -> z[[1]]
    lapply(1:length(z2), function(i) data.frame(peakSize=max(z2[[i]][[1]]$I,na.rm=TRUE), ## peak epidemic size
                                                peakPrev=max(z2[[i]][[1]]$I/(z2[[i]][[1]]$S+z2[[i]][[1]]$I+z2[[i]][[1]]$R),na.rm=TRUE),
                                                disp=ifelse(max(z2[[i]][[1]]$t)>99, glm.nb(z2[[i]][[2]]$numInf~1)$theta, NA), ## dispersion parameter of negative binomial distribution fit to numInf
                                                fadeout=ifelse(length(z2[[i]][[1]]$t)<99, 1, 0),
                                                rep=i)) %>% 
      do.call("rbind.data.frame",.) %>%
      mutate(., Variance="med",cov=covMatrix, traits=paste(trait1,trait2,sep="-"), panelLab=panel) -> z[[2]]
    lapply(1:length(z3), function(i) data.frame(peakSize=max(z3[[i]][[1]]$I,na.rm=TRUE), ## peak epidemic size
                                                peakPrev=max(z3[[i]][[1]]$I/(z3[[i]][[1]]$S+z3[[i]][[1]]$I+z3[[i]][[1]]$R),na.rm=TRUE),
                                                disp=ifelse(max(z3[[i]][[1]]$t)>99, glm.nb(z3[[i]][[2]]$numInf~1)$theta, NA), ## dispersion parameter of negative binomial distribution fit to numInf
                                                fadeout=ifelse(length(z3[[i]][[1]]$t)<99, 1, 0),
                                                rep=i)) %>% 
      do.call("rbind.data.frame",.) %>%
      mutate(., Variance="low",cov=covMatrix, traits=paste(trait1,trait2,sep="-"), panelLab=panel) -> z[[3]]
    data[[i]] <- do.call("rbind.data.frame",z)
    i <- i+1
  }
}
data <- do.call("rbind.data.frame",data)
data$cov <- factor(data$cov, levels=c("(-) Cov","(0) Cov","(+) Cov"))
data$Variance <- factor(data$Variance, levels=c("low","med","high"))
data$traits <- factor(data$traits, levels=c("Contact-Recovery","Contact-Virulence","Infectiousness-Recovery","Infectiousness-Virulence","Contact-Infectiousness","Virulence-Recovery"))
# compute the fraction of replicates that make it to each epidemic size
data %>% group_by(Variance, cov, traits) %>% summarise(size=seq(0,max(data$peakSize)+2), ECDF=sapply(seq(0,max(data$peakSize)+2), function(i) sum(peakSize>=i)/100)) -> data2
# compute the fraction of replicates that make it to each epidemic prevalence
data %>% group_by(Variance, cov, traits) %>% summarise(size=seq(0,max(data$peakPrev)+0.02,0.01), ECDF=sapply(seq(0,max(data$peakPrev)+0.02,0.01), function(i) sum(peakPrev>=i)/100)) -> data3

png(filename="./All_ECDF_R0=8.png", height=6.5, width=9, units='in', res=300)
ggplot(data2, aes(x=size, y=ECDF, group=Variance, colour=Variance)) + 
  geom_line(size=0.8) + 
  facet_grid(cov~traits) +
  ylim(0,1.05) + 
  geom_text(data=subset(data %>% group_by(Variance, cov, traits) %>% summarise(pan=unique(panelLab)),Variance=='high'),
            mapping=aes(x=2, y=1.03, label=pan),color='black', size=2.5, hjust=0) +  
  geom_text(data=subset(data %>% 
                          group_by(Variance, cov, traits) %>% 
                          summarise(fadeoutFrac=sum(fadeout==1)/100) %>% 
                          mutate(lab=paste0("Pr","=",fadeoutFrac)),
                        Variance=="low"),
            mapping=aes(x=90, y=0.9, label=lab, color=Variance), size=2.5, hjust=0) +
  geom_text(data=subset(data %>% 
                          group_by(Variance, cov, traits) %>% 
                          summarise(fadeoutFrac=sum(fadeout==1)/100) %>% 
                          mutate(lab=paste0("Pr","=",fadeoutFrac)),
                        Variance=="med"),
            mapping=aes(x=90, y=0.8, label=lab, color=Variance), size=2.5, hjust=0) +
  geom_text(data=subset(data %>% 
                          group_by(Variance, cov, traits) %>% 
                          summarise(fadeoutFrac=sum(fadeout==1)/100) %>% 
                          mutate(lab=paste0("Pr","=",fadeoutFrac)),
                        Variance=="high"),
            mapping=aes(x=90, y=0.7, label=lab, color=Variance), size=2.5, hjust=0) +
  xlab("Epidemic size") + 
  ylab(expression("%"~of~simulations~(italic(R)[0]==8))) + 
  ggtitle("Intergroup trait pairings                                                                                                                               Intragroup trait pairings") + 
  theme_bw() + 
  theme(axis.text=element_text(size=6), axis.title=element_text(size=8), legend.text=element_text(size=7), legend.title=element_text(size=8), plot.title = element_text(size=8), strip.text.x = element_text(size = 7), strip.text.y = element_text(size = 7)) +
  scale_color_manual(values=c("cornflowerblue", "pink", "sienna")) +
  scale_fill_manual(values=c("cornflowerblue", "pink", "sienna")) 
dev.off()

png(filename="./All_fadeout_peaks_R0=8.png", height=6.5, width=9, units='in', res=300)
ggplot(subset(data, fadeout==1), aes(x=peakSize)) + 
  geom_histogram(aes(group=Variance,color=Variance,fill=Variance),bins=10,alpha=0.3,position='identity')+
  facet_grid(cov~traits) +
  ylim(0,66) + 
  geom_text(data=subset(data %>% group_by(Variance, cov, traits) %>% summarise(pan=unique(panelLab)),Variance=='high'),
            mapping=aes(x=-5, y=64, label=pan),color='black', size=2.5, hjust=0) +  
  xlab("Peak") + 
  ylab(expression("#"~of~simulations~(italic(R)[0]==8))) + 
  ggtitle("Intergroup trait pairings                                                                                                                               Intragroup trait pairings") + 
  theme_bw() + 
  theme(axis.text=element_text(size=6), axis.title=element_text(size=8), legend.text=element_text(size=7), legend.title=element_text(size=8), plot.title = element_text(size=8), strip.text.x = element_text(size = 7), strip.text.y = element_text(size = 7)) +
  scale_color_manual(values=c("sienna")) +
  scale_fill_manual(values=c("sienna")) 
dev.off()

png(filename="./All_dispersion_values_R0=8.png", height=6.5, width=9, units='in', res=300)
ggplot(filter(data, fadeout==0, disp < 2), aes(x=disp)) + 
  geom_density(aes(group=Variance,color=Variance,fill=Variance),alpha=0.3,position='identity')+
  geom_vline(xintercept = 0.16, linetype="dashed") +
  facet_grid(cov~traits) +
  ylim(0,63) + 
  geom_text(data=subset(data %>% group_by(Variance, cov, traits) %>% summarise(pan=unique(panelLab)),Variance=='high'),
            mapping=aes(x=0.02, y=61, label=pan),color='black', size=2.5, hjust=0) +  
  xlab(expression(Dispersion~(italic(k)))) +
  ylab(expression(Density~(italic(R)[0]==8))) + 
  ggtitle("Intergroup trait pairings                                                                                                                               Intragroup trait pairings") + 
  theme_bw() + 
  theme(axis.text=element_text(size=6), axis.title=element_text(size=8), legend.text=element_text(size=7), legend.title=element_text(size=8), plot.title = element_text(size=8), strip.text.x = element_text(size = 7), strip.text.y = element_text(size = 7)) +
  scale_color_manual(values=c("cornflowerblue", "pink", "sienna")) +
  scale_fill_manual(values=c("cornflowerblue", "pink", "sienna")) 
dev.off()


png(filename="./All_peak_against_dispersion_R0=8.png", height=2.5, width=9, units='in', res=300)
ggplot(filter(mutate(data, Covariance=cov), fadeout==0, disp < 2), aes(x=disp, y=peakSize, shape=Variance, color=Covariance)) + 
  geom_point(size=0.8) + 
  geom_vline(xintercept = 0.16, linetype="dashed") +
  facet_grid(~traits) +
  ylim(115,220) + 
  geom_text(data=subset(mutate(data, Covariance=cov) %>% group_by(Variance, Covariance, traits) %>% summarise(pan=unique(strsplit(panelLab,"\\(")[[1]][1])),Variance=='high' & Covariance=="(0) Cov"),
            mapping=aes(x=0.075, y=219, label=pan),color='black', size=2.5, hjust=0) +  
  xlab(expression(Dispersion~(italic(k)))) +
  ylab(expression(Peak~size~(italic(R)[0]==8))) + 
  ggtitle("Intergroup trait pairings                                                                                                                               Intragroup trait pairings") + 
  theme_bw() + 
  theme(axis.text=element_text(size=6), axis.title=element_text(size=8), legend.text=element_text(size=7), legend.title=element_text(size=8), plot.title = element_text(size=8), strip.text.x = element_text(size = 7), strip.text.y = element_text(size = 7)) +
  scale_color_manual(values=c("cornflowerblue", "pink", "sienna")) +
  scale_fill_manual(values=c("cornflowerblue", "pink", "sienna")) 
dev.off()





## Supplementary Figure showing equilibrium R
var1=c('c',    'c',    'shed', 'shed')
var2=c('gamma','alpha','gamma','alpha')
i <- 1
data <- vector(mode='list', length=12)
for (j in 1:4) { ## loop over the six different covariance combinations
  for (covMatrix in c("nocorr", "negcorr", "poscorr")) {
    z1 = readRDS(paste0(paste("out_R=4",covMatrix,"hivar",var1[j],var2[j],sep="_"),".RDS"))
    z2 = readRDS(paste0(paste("out_R=4",covMatrix,"medvar",var1[j],var2[j],sep="_"),".RDS"))
    z3 = readRDS(paste0(paste("out_R=4",covMatrix,"lowvar",var1[j],var2[j],sep="_"),".RDS"))
    
    ## Give nicer names
    covMatrix <- switch(covMatrix,poscorr="(+) Cov",negcorr="(-) Cov",nocorr="(0) Cov")
    trait1 <- switch(var1[j],c="Contact",shed="Infectiousness")
    trait2 <- switch(var2[j],alpha="Virulence",gamma="Recovery")
    
    z <- vector(mode='list',length=3)
    lapply(1:length(z1), function(i) data.frame(peakI=max(z1[[i]][[1]]$I,na.rm=T), 
                                                minS=min(z1[[i]][[1]]$S,na.rm=T), 
                                                equilR=mean(tail(z1[[i]][[1]]$R,10),na.rm=T), 
                                                equilI=mean(tail(z1[[i]][[1]]$I,10),na.rm=T),
                                                equilPrev=mean(tail(z3[[i]][[1]]$I/(z3[[i]][[1]]$S+z3[[i]][[1]]$I+z3[[i]][[1]]$R),10),na.rm=T),
                                                disp=ifelse(max(z1[[i]][[1]]$t)>99, glm.nb(z1[[i]][[2]]$numInf~1)$theta, NA), ## dispersion parameter of negative binomial distribution fit to numInf
                                                rep=i)) %>% 
      do.call("rbind.data.frame",.) %>%
      mutate(., Variance="high",cov=covMatrix, traits=paste(trait1,trait2,sep="-")) -> z[[1]]
    
    lapply(1:length(z2), function(i) data.frame(peakI=max(z2[[i]][[1]]$I,na.rm=T), 
                                                minS=min(z2[[i]][[1]]$S,na.rm=T), 
                                                equilR=mean(tail(z2[[i]][[1]]$R,10),na.rm=T), 
                                                equilI=mean(tail(z2[[i]][[1]]$I,10),na.rm=T),
                                                equilPrev=mean(tail(z3[[i]][[1]]$I/(z3[[i]][[1]]$S+z3[[i]][[1]]$I+z3[[i]][[1]]$R),10),na.rm=T),
                                                disp=ifelse(max(z2[[i]][[1]]$t)>99, glm.nb(z2[[i]][[2]]$numInf~1)$theta, NA), ## dispersion parameter of negative binomial distribution fit to numInf
                                                rep=i)) %>% 
      do.call("rbind.data.frame",.) %>%
      mutate(., Variance="med",cov=covMatrix, traits=paste(trait1,trait2,sep="-")) -> z[[2]]
    
    lapply(1:length(z3), function(i) data.frame(peakI=max(z3[[i]][[1]]$I,na.rm=T), 
                                                minS=min(z3[[i]][[1]]$S,na.rm=T), 
                                                equilR=mean(tail(z3[[i]][[1]]$R,10),na.rm=T), 
                                                equilI=mean(tail(z3[[i]][[1]]$I,10),na.rm=T),
                                                equilPrev=mean(tail(z3[[i]][[1]]$I/(z3[[i]][[1]]$S+z3[[i]][[1]]$I+z3[[i]][[1]]$R),10),na.rm=T),
                                                disp=ifelse(max(z3[[i]][[1]]$t)>99, glm.nb(z3[[i]][[2]]$numInf~1)$theta, NA), ## dispersion parameter of negative binomial distribution fit to numInf
                                                rep=i)) %>% 
      do.call("rbind.data.frame",.) %>%
      mutate(., Variance="low",cov=covMatrix, traits=paste(trait1,trait2,sep="-")) -> z[[3]]
    data[[i]] <- do.call("rbind.data.frame",z)
    i <- i+1
  }
}
data <- do.call("rbind.data.frame",data)
data$cov <- factor(data$cov, levels=c("(-) Cov","(0) Cov","(+) Cov"))
data$Variance <- factor(data$Variance, levels=c("low","med","high"))
data$traits <- factor(data$traits, levels=c("Contact-Recovery","Contact-Virulence","Infectiousness-Recovery","Infectiousness-Virulence"))

png(filename="./Intergroup_equilibrium_R_R0=4.png", height=2.5, width=9, units='in', res=300)
ggplot(subset(data,peakI>50), aes(x=equilR)) + 
  geom_density(aes(group=Variance,color=Variance,fill=Variance),alpha=0.3,position='identity') +
  facet_grid(cov~traits) +
  xlab("# of recovered hosts") +
  ylab(expression(Density~(italic(R)[0]==4))) + 
  ggtitle("Intergroup trait pairings") + 
  theme_bw() + 
  theme(axis.text=element_text(size=6), axis.title=element_text(size=8), legend.text=element_text(size=7), legend.title=element_text(size=8), plot.title = element_text(size=8), strip.text.x = element_text(size = 7), strip.text.y = element_text(size = 7)) +
  scale_color_manual(values=c("cornflowerblue", "pink", "sienna")) +
  scale_fill_manual(values=c("cornflowerblue", "pink", "sienna")) 
dev.off()

ggplot(subset(data,peakI>50), aes(x=disp, y=equilPrev, group=Variance,color=Variance,fill=Variance)) +
  geom_point() + 
  facet_grid(cov~traits)
  
  
  
  geom_density(aes(group=Variance,color=Variance,fill=Variance),alpha=0.3,position='identity') +
  facet_grid(cov~traits) +
  xlab("# of recovered hosts") +
  ylab(expression(Density~(italic(R)[0]==4))) + 
  ggtitle("Intergroup trait pairings") + 
  theme_bw() + 
  theme(axis.text=element_text(size=6), axis.title=element_text(size=8), legend.text=element_text(size=7), legend.title=element_text(size=8), plot.title = element_text(size=8), strip.text.x = element_text(size = 7), strip.text.y = element_text(size = 7)) +
  scale_color_manual(values=c("cornflowerblue", "pink", "sienna")) +
  scale_fill_manual(values=c("cornflowerblue", "pink", "sienna")) 
