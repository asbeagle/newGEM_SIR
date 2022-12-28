
## R0 = 1 plots
var1=c('c','c','c','shed','shed','alpha')
var2=c('shed','alpha','gamma','alpha','gamma','gamma')
i <- 1
data <- vector(mode='list', length=18)
for (j in 1:6) { ## loop over the six different covariance combinations
  for (covMatrix in c("nocorr", "negcorr", "poscorr")) {
    z1 = readRDS(paste0(paste("out_R=1",covMatrix,"hivar",var1[j],var2[j],sep="_"),".RDS"))
    z2 = readRDS(paste0(paste("out_R=1",covMatrix,"medvar",var1[j],var2[j],sep="_"),".RDS"))
    z3 = readRDS(paste0(paste("out_R=1",covMatrix,"lowvar",var1[j],var2[j],sep="_"),".RDS"))

    ## Give nicer names
    covMatrix <- switch(covMatrix,poscorr="(+)",negcorr="(-)",nocorr="(0)")
    trait1 <- switch(var1[j],c="Contact",shed="Infectiousness",alpha="Virulence")
    trait2 <- switch(var2[j],shed="Infectiousness",alpha="Virulence",gamma="Recovery")

    z <- vector(mode='list',length=3)
    lapply(1:length(z1), function(i) mutate(z1[[i]][[1]][1:(ifelse(any(z1[[i]][[1]]$I==0),min(which(z1[[i]][[1]]$I==0)),length(z1[[i]][[1]]$I))),],
                                            rep=i,t=ceiling(t))) %>%
      do.call("rbind.data.frame",.) %>%
      mutate(., Variance="high",cov=covMatrix, traits=paste(trait1,trait2,sep="-")) -> z[[1]]
    lapply(1:length(z2), function(i) mutate(z2[[i]][[1]][1:(ifelse(any(z2[[i]][[1]]$I==0),min(which(z2[[i]][[1]]$I==0)),length(z2[[i]][[1]]$I))),],
                                            rep=i,t=ceiling(t))) %>%
      do.call("rbind.data.frame",.) %>%
      mutate(., Variance="med",cov=covMatrix, traits=paste(trait1,trait2,sep="-")) -> z[[2]]
    lapply(1:length(z3), function(i) mutate(z3[[i]][[1]][1:(ifelse(any(z3[[i]][[1]]$I==0),min(which(z3[[i]][[1]]$I==0)),length(z3[[i]][[1]]$I))),],
                                             rep=i,t=ceiling(t))) %>%
      do.call("rbind.data.frame",.) %>%
      mutate(., Variance="low",cov=covMatrix, traits=paste(trait1,trait2,sep="-")) -> z[[3]]
    data[[i]] <- do.call("rbind.data.frame",z)
    i <- i+1
  }
}
data <- do.call("rbind.data.frame",data)
data$cov <- factor(data$cov, levels=c("(-)","(0)","(+)"))
data$Variance <- factor(data$Variance, levels=c("low","med","high"))
data$traits <- factor(data$traits, levels=c("Contact-Infectiousness","Virulence-Recovery","Contact-Recovery","Contact-Virulence","Infectiousness-Recovery","Infectiousness-Virulence"))

data %>% group_by(Variance, cov, traits, rep) %>% summarise(maxI=max(I), maxT=max(t)) -> data2

data2 %>% group_by(Variance, cov, traits) %>% summarise(peak=5:200, ECDF=sapply(5:200, function(i) 1-sum(maxI>=i)/100)) -> data3

png(filename="./Intragroup_ECDF_R0=1.png", height=3.9, width=3.9, units='in', res=300)
ggplot(subset(data3, traits%in%c("Contact-Infectiousness","Virulence-Recovery")), aes(x=peak, y=ECDF, group=Variance, colour=Variance)) + 
  geom_line(size=0.8) + 
  facet_grid(cov~traits) +
  xlab("Peak") + ylab("% of simulations") + ggtitle("Intragroup trait pairings") + theme(axis.text=element_text(size=6), axis.title=element_text(size=6), legend.text=element_text(size=6), legend.title=element_text(size=6), plot.title = element_text(hjust = 0.5, size=8), strip.text.x = element_text(size = 6), strip.text.y = element_text(size = 6)) +
  geom_vline(data=(data2 %>% group_by(Variance, cov, traits) %>% summarise(fadeoutPeak=ifelse(any(maxT<100),max(maxI[maxT<100]),NA)) %>% filter(traits%in%c("Contact-Infectiousness","Virulence-Recovery"))), aes(xintercept=fadeoutPeak, colour=Variance), linetype='dotted', size=0.6) +  
  geom_hline(data=(data2 %>% group_by(Variance, cov, traits) %>% summarise(fadeoutFrac=ifelse(any(maxT<100),sum(maxT<100)/100,0)) %>% filter(traits%in%c("Contact-Infectiousness","Virulence-Recovery"))), aes(yintercept=fadeoutFrac, colour=Variance), linetype='dashed', size=0.3)
dev.off()

png(filename="./Intergroup_ECDF_R0=1.png", height=3.9, width=6.5, units='in', res=300)
ggplot(subset(data3, traits%in%c("Contact-Virulence","Contact-Recovery","Infectiousness-Virulence","Infectiousness-Recovery")), aes(x=peak, y=ECDF, group=Variance, colour=Variance)) + 
  geom_line(size=0.8) + 
  facet_grid(cov~traits) +
  xlab("Peak") + ylab("% of simulations") + ggtitle("Intergroup trait pairings") + theme(axis.text=element_text(size=6), axis.title=element_text(size=6), legend.text=element_text(size=6), legend.title=element_text(size=6), plot.title = element_text(hjust = 0.5, size=8), strip.text.x = element_text(size = 6), strip.text.y = element_text(size = 6)) +
  geom_vline(data=(data2 %>% group_by(Variance, cov, traits) %>% summarise(fadeoutPeak=ifelse(any(maxT<100),max(maxI[maxT<100]),NA)) %>% filter(traits%in%c("Contact-Virulence","Contact-Recovery","Infectiousness-Virulence","Infectiousness-Recovery"))), aes(xintercept=fadeoutPeak, colour=Variance), linetype='dotted', size=0.6) +  
  geom_hline(data=(data2 %>% group_by(Variance, cov, traits) %>% summarise(fadeoutFrac=ifelse(any(maxT<100),sum(maxT<100)/100,0)) %>% filter(traits%in%c("Contact-Virulence","Contact-Recovery","Infectiousness-Virulence","Infectiousness-Recovery"))), aes(yintercept=fadeoutFrac, colour=Variance), linetype='dashed', size=0.3)
dev.off()

png(filename="./All_ECDF_R0=1.png", height=6.5, width=9, units='in', res=300)
ggplot(data3, aes(x=peak, y=ECDF, group=Variance, colour=Variance)) + 
  geom_line(size=0.8) + 
  facet_grid(cov~traits) +
  xlab("Peak") + ylab("% of simulations") + ggtitle("Intragroup trait pairings                                              Intergroup trait pairings") + theme(axis.text=element_text(size=6), axis.title=element_text(size=6), legend.text=element_text(size=6), legend.title=element_text(size=6), plot.title = element_text(size=8), strip.text.x = element_text(size = 6), strip.text.y = element_text(size = 6)) +
  geom_vline(data=(data2 %>% group_by(Variance, cov, traits) %>% summarise(fadeoutPeak=ifelse(any(maxT<100),max(maxI[maxT<100]),NA))), aes(xintercept=fadeoutPeak, colour=Variance), linetype='dotted', size=0.6) +  
  geom_hline(data=(data2 %>% group_by(Variance, cov, traits) %>% summarise(fadeoutFrac=ifelse(any(maxT<100),sum(maxT<100)/100,0))), aes(yintercept=fadeoutFrac, colour=Variance), linetype='dashed', size=0.3)
dev.off()

  

## R0 = 4 plots
var1=c('c','c','c','shed','shed','alpha')
var2=c('shed','alpha','gamma','alpha','gamma','gamma')
i <- 1
data <- vector(mode='list', length=18)
for (j in 1:6) { ## loop over the six different covariance combinations
  for (covMatrix in c("nocorr", "negcorr", "poscorr")) {
    z1 = readRDS(paste0(paste("out_R=4",covMatrix,"hivar",var1[j],var2[j],sep="_"),".RDS"))
    z2 = readRDS(paste0(paste("out_R=4",covMatrix,"medvar",var1[j],var2[j],sep="_"),".RDS"))
    z3 = readRDS(paste0(paste("out_R=4",covMatrix,"lowvar",var1[j],var2[j],sep="_"),".RDS"))
    
    ## Give nicer names
    covMatrix <- switch(covMatrix,poscorr="(+)",negcorr="(-)",nocorr="(0)")
    trait1 <- switch(var1[j],c="Contact",shed="Infectiousness",alpha="Virulence")
    trait2 <- switch(var2[j],shed="Infectiousness",alpha="Virulence",gamma="Recovery")
    
    z <- vector(mode='list',length=3)
    lapply(1:length(z1), function(i) mutate(z1[[i]][[1]][1:(ifelse(any(z1[[i]][[1]]$I==0),min(which(z1[[i]][[1]]$I==0)),length(z1[[i]][[1]]$I))),],
                                            rep=i,t=ceiling(t))) %>%
      do.call("rbind.data.frame",.) %>%
      mutate(., Variance="high",cov=covMatrix, traits=paste(trait1,trait2,sep="-")) -> z[[1]]
    lapply(1:length(z2), function(i) mutate(z2[[i]][[1]][1:(ifelse(any(z2[[i]][[1]]$I==0),min(which(z2[[i]][[1]]$I==0)),length(z2[[i]][[1]]$I))),],
                                            rep=i,t=ceiling(t))) %>%
      do.call("rbind.data.frame",.) %>%
      mutate(., Variance="med",cov=covMatrix, traits=paste(trait1,trait2,sep="-")) -> z[[2]]
    lapply(1:length(z3), function(i) mutate(z3[[i]][[1]][1:(ifelse(any(z3[[i]][[1]]$I==0),min(which(z3[[i]][[1]]$I==0)),length(z3[[i]][[1]]$I))),],
                                            rep=i,t=ceiling(t))) %>%
      do.call("rbind.data.frame",.) %>%
      mutate(., Variance="low",cov=covMatrix, traits=paste(trait1,trait2,sep="-")) -> z[[3]]
    data[[i]] <- do.call("rbind.data.frame",z)
    i <- i+1
  }
}
data <- do.call("rbind.data.frame",data)
data$cov <- factor(data$cov, levels=c("(-)","(0)","(+)"))
data$Variance <- factor(data$Variance, levels=c("low","med","high"))
data$traits <- factor(data$traits, levels=c("Contact-Infectiousness","Virulence-Recovery","Contact-Recovery","Contact-Virulence","Infectiousness-Recovery","Infectiousness-Virulence"))

data %>% group_by(Variance, cov, traits, rep) %>% summarise(maxI=max(I), maxT=max(t)) -> data2

data2 %>% group_by(Variance, cov, traits) %>% summarise(peak=5:200, ECDF=sapply(5:200, function(i) 1-sum(maxI>=i)/100)) -> data3

png(filename="./Intragroup_ECDF_R0=4.png", height=3.9, width=3.9, units='in', res=300)
ggplot(subset(data3, traits%in%c("Contact-Infectiousness","Virulence-Recovery")), aes(x=peak, y=ECDF, group=Variance, colour=Variance)) + 
  geom_line(size=0.8) + 
  facet_grid(cov~traits) +
  xlab("Peak") + ylab("% of simulations") + ggtitle("Intragroup trait pairings") + theme(axis.text=element_text(size=6), axis.title=element_text(size=6), legend.text=element_text(size=6), legend.title=element_text(size=6), plot.title = element_text(hjust = 0.5, size=8), strip.text.x = element_text(size = 6), strip.text.y = element_text(size = 6)) +
  geom_vline(data=(data2 %>% group_by(Variance, cov, traits) %>% summarise(fadeoutPeak=ifelse(any(maxT<100),max(maxI[maxT<100]),NA)) %>% filter(traits%in%c("Contact-Infectiousness","Virulence-Recovery"))), aes(xintercept=fadeoutPeak, colour=Variance), linetype='dotted', size=0.6) +  
  geom_hline(data=(data2 %>% group_by(Variance, cov, traits) %>% summarise(fadeoutFrac=ifelse(any(maxT<100),sum(maxT<100)/100,0)) %>% filter(traits%in%c("Contact-Infectiousness","Virulence-Recovery"))), aes(yintercept=fadeoutFrac, colour=Variance), linetype='dashed', size=0.3)
dev.off()

png(filename="./Intergroup_ECDF_R0=4.png", height=3.9, width=6.5, units='in', res=300)
ggplot(subset(data3, traits%in%c("Contact-Virulence","Contact-Recovery","Infectiousness-Virulence","Infectiousness-Recovery")), aes(x=peak, y=ECDF, group=Variance, colour=Variance)) + 
  geom_line(size=0.8) + 
  facet_grid(cov~traits) +
  xlab("Peak") + ylab("% of simulations") + ggtitle("Intergroup trait pairings") + theme(axis.text=element_text(size=6), axis.title=element_text(size=6), legend.text=element_text(size=6), legend.title=element_text(size=6), plot.title = element_text(hjust = 0.5, size=8), strip.text.x = element_text(size = 6), strip.text.y = element_text(size = 6)) +
  geom_vline(data=(data2 %>% group_by(Variance, cov, traits) %>% summarise(fadeoutPeak=ifelse(any(maxT<100),max(maxI[maxT<100]),NA)) %>% filter(traits%in%c("Contact-Virulence","Contact-Recovery","Infectiousness-Virulence","Infectiousness-Recovery"))), aes(xintercept=fadeoutPeak, colour=Variance), linetype='dotted', size=0.6) +  
  geom_hline(data=(data2 %>% group_by(Variance, cov, traits) %>% summarise(fadeoutFrac=ifelse(any(maxT<100),sum(maxT<100)/100,0)) %>% filter(traits%in%c("Contact-Virulence","Contact-Recovery","Infectiousness-Virulence","Infectiousness-Recovery"))), aes(yintercept=fadeoutFrac, colour=Variance), linetype='dashed', size=0.3)
dev.off()

png(filename="./All_ECDF_R0=4.png", height=6.5, width=9, units='in', res=300)
ggplot(data3, aes(x=peak, y=ECDF, group=Variance, colour=Variance)) + 
  geom_line(size=0.8) + 
  facet_grid(cov~traits) +
  xlab("Peak") + ylab("% of simulations") + ggtitle("Intragroup trait pairings                                              Intergroup trait pairings") + theme(axis.text=element_text(size=6), axis.title=element_text(size=6), legend.text=element_text(size=6), legend.title=element_text(size=6), plot.title = element_text(size=8), strip.text.x = element_text(size = 6), strip.text.y = element_text(size = 6)) +
  geom_vline(data=(data2 %>% group_by(Variance, cov, traits) %>% summarise(fadeoutPeak=ifelse(any(maxT<100),max(maxI[maxT<100]),NA))), aes(xintercept=fadeoutPeak, colour=Variance), linetype='dotted', size=0.6) +  
  geom_hline(data=(data2 %>% group_by(Variance, cov, traits) %>% summarise(fadeoutFrac=ifelse(any(maxT<100),sum(maxT<100)/100,0))), aes(yintercept=fadeoutFrac, colour=Variance), linetype='dashed', size=0.3)
dev.off()


## R0 = 8 plots
var1=c('c','c','c','shed','shed','alpha')
var2=c('shed','alpha','gamma','alpha','gamma','gamma')
i <- 1
data <- vector(mode='list', length=18)
for (j in 1:6) { ## loop over the six different covariance combinations
  for (covMatrix in c("nocorr", "negcorr", "poscorr")) {
    z1 = readRDS(paste0(paste("out_R=8",covMatrix,"hivar",var1[j],var2[j],sep="_"),".RDS"))
    z2 = readRDS(paste0(paste("out_R=8",covMatrix,"medvar",var1[j],var2[j],sep="_"),".RDS"))
    z3 = readRDS(paste0(paste("out_R=8",covMatrix,"lowvar",var1[j],var2[j],sep="_"),".RDS"))
    
    ## Give nicer names
    covMatrix <- switch(covMatrix,poscorr="(+)",negcorr="(-)",nocorr="(0)")
    trait1 <- switch(var1[j],c="Contact",shed="Infectiousness",alpha="Virulence")
    trait2 <- switch(var2[j],shed="Infectiousness",alpha="Virulence",gamma="Recovery")
    
    z <- vector(mode='list',length=3)
    lapply(1:length(z1), function(i) mutate(z1[[i]][[1]][1:(ifelse(any(z1[[i]][[1]]$I==0),min(which(z1[[i]][[1]]$I==0)),length(z1[[i]][[1]]$I))),],
                                            rep=i,t=ceiling(t))) %>%
      do.call("rbind.data.frame",.) %>%
      mutate(., Variance="high",cov=covMatrix, traits=paste(trait1,trait2,sep="-")) -> z[[1]]
    lapply(1:length(z2), function(i) mutate(z2[[i]][[1]][1:(ifelse(any(z2[[i]][[1]]$I==0),min(which(z2[[i]][[1]]$I==0)),length(z2[[i]][[1]]$I))),],
                                            rep=i,t=ceiling(t))) %>%
      do.call("rbind.data.frame",.) %>%
      mutate(., Variance="med",cov=covMatrix, traits=paste(trait1,trait2,sep="-")) -> z[[2]]
    lapply(1:length(z3), function(i) mutate(z3[[i]][[1]][1:(ifelse(any(z3[[i]][[1]]$I==0),min(which(z3[[i]][[1]]$I==0)),length(z3[[i]][[1]]$I))),],
                                            rep=i,t=ceiling(t))) %>%
      do.call("rbind.data.frame",.) %>%
      mutate(., Variance="low",cov=covMatrix, traits=paste(trait1,trait2,sep="-")) -> z[[3]]
    data[[i]] <- do.call("rbind.data.frame",z)
    i <- i+1
  }
}
data <- do.call("rbind.data.frame",data)
data$cov <- factor(data$cov, levels=c("(-)","(0)","(+)"))
data$Variance <- factor(data$Variance, levels=c("low","med","high"))
data$traits <- factor(data$traits, levels=c("Contact-Infectiousness","Virulence-Recovery","Contact-Recovery","Contact-Virulence","Infectiousness-Recovery","Infectiousness-Virulence"))

data %>% group_by(Variance, cov, traits, rep) %>% summarise(maxI=max(I), maxT=max(t)) -> data2

data2 %>% group_by(Variance, cov, traits) %>% summarise(peak=5:200, ECDF=sapply(5:200, function(i) 1-sum(maxI>=i)/100)) -> data3

png(filename="./Intragroup_ECDF_R0=8.png", height=3.9, width=3.9, units='in', res=300)
ggplot(subset(data3, traits%in%c("Contact-Infectiousness","Virulence-Recovery")), aes(x=peak, y=ECDF, group=Variance, colour=Variance)) + 
  geom_line(size=0.8) + 
  facet_grid(cov~traits) +
  xlab("Peak") + ylab("% of simulations") + ggtitle("Intragroup trait pairings") + theme(axis.text=element_text(size=6), axis.title=element_text(size=6), legend.text=element_text(size=6), legend.title=element_text(size=6), plot.title = element_text(hjust = 0.5, size=8), strip.text.x = element_text(size = 6), strip.text.y = element_text(size = 6)) +
  geom_vline(data=(data2 %>% group_by(Variance, cov, traits) %>% summarise(fadeoutPeak=ifelse(any(maxT<100),max(maxI[maxT<100]),NA)) %>% filter(traits%in%c("Contact-Infectiousness","Virulence-Recovery"))), aes(xintercept=fadeoutPeak, colour=Variance), linetype='dotted', size=0.6) +  
  geom_hline(data=(data2 %>% group_by(Variance, cov, traits) %>% summarise(fadeoutFrac=ifelse(any(maxT<100),sum(maxT<100)/100,0)) %>% filter(traits%in%c("Contact-Infectiousness","Virulence-Recovery"))), aes(yintercept=fadeoutFrac, colour=Variance), linetype='dashed', size=0.3)
dev.off()

png(filename="./Intergroup_ECDF_R0=8.png", height=3.9, width=6.5, units='in', res=300)
ggplot(subset(data3, traits%in%c("Contact-Virulence","Contact-Recovery","Infectiousness-Virulence","Infectiousness-Recovery")), aes(x=peak, y=ECDF, group=Variance, colour=Variance)) + 
  geom_line(size=0.8) + 
  facet_grid(cov~traits) +
  xlab("Peak") + ylab("% of simulations") + ggtitle("Intergroup trait pairings") + theme(axis.text=element_text(size=6), axis.title=element_text(size=6), legend.text=element_text(size=6), legend.title=element_text(size=6), plot.title = element_text(hjust = 0.5, size=8), strip.text.x = element_text(size = 6), strip.text.y = element_text(size = 6)) +
  geom_vline(data=(data2 %>% group_by(Variance, cov, traits) %>% summarise(fadeoutPeak=ifelse(any(maxT<100),max(maxI[maxT<100]),NA)) %>% filter(traits%in%c("Contact-Virulence","Contact-Recovery","Infectiousness-Virulence","Infectiousness-Recovery"))), aes(xintercept=fadeoutPeak, colour=Variance), linetype='dotted', size=0.6) +  
  geom_hline(data=(data2 %>% group_by(Variance, cov, traits) %>% summarise(fadeoutFrac=ifelse(any(maxT<100),sum(maxT<100)/100,0)) %>% filter(traits%in%c("Contact-Virulence","Contact-Recovery","Infectiousness-Virulence","Infectiousness-Recovery"))), aes(yintercept=fadeoutFrac, colour=Variance), linetype='dashed', size=0.3)
dev.off()

png(filename="./All_ECDF_R0=8.png", height=6.5, width=9, units='in', res=300)
ggplot(data3, aes(x=peak, y=ECDF, group=Variance, colour=Variance)) + 
  geom_line(size=0.8) + 
  facet_grid(cov~traits) +
  xlab("Peak") + ylab("% of simulations") + ggtitle("Intragroup trait pairings                                              Intergroup trait pairings") + theme(axis.text=element_text(size=6), axis.title=element_text(size=6), legend.text=element_text(size=6), legend.title=element_text(size=6), plot.title = element_text(size=8), strip.text.x = element_text(size = 6), strip.text.y = element_text(size = 6)) +
  geom_vline(data=(data2 %>% group_by(Variance, cov, traits) %>% summarise(fadeoutPeak=ifelse(any(maxT<100),max(maxI[maxT<100]),NA))), aes(xintercept=fadeoutPeak, colour=Variance), linetype='dotted', size=0.6) +  
  geom_hline(data=(data2 %>% group_by(Variance, cov, traits) %>% summarise(fadeoutFrac=ifelse(any(maxT<100),sum(maxT<100)/100,0))), aes(yintercept=fadeoutFrac, colour=Variance), linetype='dashed', size=0.3)
dev.off()
