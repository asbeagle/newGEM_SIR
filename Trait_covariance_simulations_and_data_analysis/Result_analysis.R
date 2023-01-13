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
    lapply(1:length(z1), function(i) z1[[i]][[1]] %>% 
             mutate(N=S+I+R) %>% 
             summarise(peakI=max(I,na.rm=T), 
                       minS=min(S,na.rm=T), 
                       minN=min(N,na.rm=T), 
                       peakR=max(R,na.rm=T),
                       peakIt=t[which(I==peakI)[1]], 
                       minSt=t[which(S==minS)[1]], 
                       minNt=t[which(N==minN)[1]], 
                       rep=i)) %>% 
      do.call("rbind.data.frame",.) %>%
      mutate(., Variance="high",cov=covMatrix, traits=paste(trait1,trait2,sep="-")) -> z[[1]]
    
    lapply(1:length(z2), function(i) z2[[i]][[1]] %>% 
             mutate(N=S+I+R) %>% 
             summarise(peakI=max(I,na.rm=T), 
                       minS=min(S,na.rm=T), 
                       minN=min(N,na.rm=T), 
                       peakR=max(R,na.rm=T),
                       peakIt=t[which(I==peakI)[1]], 
                       minSt=t[which(S==minS)[1]], 
                       minNt=t[which(N==minN)[1]], 
                       rep=i)) %>% 
      do.call("rbind.data.frame",.) %>%
      mutate(., Variance="med",cov=covMatrix, traits=paste(trait1,trait2,sep="-")) -> z[[2]]
    
    lapply(1:length(z3), function(i) z3[[i]][[1]] %>% 
             mutate(N=S+I+R) %>% 
             summarise(peakI=max(I,na.rm=T), 
                       minS=min(S,na.rm=T), 
                       minN=min(N,na.rm=T), 
                       peakR=max(R,na.rm=T),
                       equilR=mean(tail(R,100),na.rm=T), 
                       rep=i)) %>% 
      do.call("rbind.data.frame",.) %>%
      mutate(., Variance="low",cov=covMatrix, traits=paste(trait1,trait2,sep="-")) -> z[[3]]
    data[[i]] <- do.call("rbind.data.frame",z)
    i <- i+1
  }
}
data <- do.call("rbind.data.frame",data)

ggplot(subset(data,peakI>50), aes(x=minS)) + 
  geom_density(aes(group=Variance,color=Variance,fill=Variance),alpha=0.3,position='identity') +
  xlim(0,70) + 
  facet_grid(traits~cov) 
  
ggplot(subset(data,peakI>50), aes(x=peakI)) + 
  geom_density(aes(group=Variance,color=Variance,fill=Variance),alpha=0.3,position='identity') +
  facet_grid(traits~cov) 

ggplot(subset(data,peakI>50), aes(x=minN)) + 
  geom_density(aes(group=Variance,color=Variance,fill=Variance),alpha=0.3,position='identity') +
  facet_grid(traits~cov) 

ggplot(subset(data,peakI>50), aes(x=peakR)) + 
  geom_density(aes(group=Variance,color=Variance,fill=Variance),alpha=0.3,position='identity') +
  facet_grid(traits~cov) 



data$cov <- factor(data$cov, levels=c("(-) Cov","(0) Cov","(+) Cov"))
data$Variance <- factor(data$Variance, levels=c("low","med","high"))
data$traits <- factor(data$traits, levels=c("Contact-Recovery","Contact-Virulence","Infectiousness-Recovery","Infectiousness-Virulence","Contact-Infectiousness","Virulence-Recovery"))
# compute the fraction of replicates that make it to each epidemic size
data %>% group_by(Variance, cov, traits) %>% summarise(size=seq(0,max(data$peakSize)+2), ECDF=sapply(seq(0,max(data$peakSize)+2), function(i) sum(peakSize>=i)/100)) -> data2
# compute the fraction of replicates that make it to each epidemic prevalence
data %>% group_by(Variance, cov, traits) %>% summarise(size=seq(0,max(data$peakPrev)+0.02,0.01), ECDF=sapply(seq(0,max(data$peakPrev)+0.02,0.01), function(i) sum(peakPrev>=i)/100)) -> data3
