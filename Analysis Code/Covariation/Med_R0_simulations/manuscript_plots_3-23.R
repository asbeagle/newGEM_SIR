###plot code for manuscript

library(tidyverse)
library(MASS)
library(ggplot2)
library(cowplot)
library(ggpubr)
library(ggtext)
library(magrittr)

## Stats on super spreading (k) and peak epidemic size
var1=c('c','c','c','shed','shed','alpha')
var2=c('shed','alpha','gamma','alpha','gamma','gamma')
i <- 1
data <- vector(mode='list', length=18)
for (j in 1:6) { ## loop over the six different covariance combinations
  for (covMatrix in c("nocorr", "negcorr", "poscorr")) {
    z1 = readRDS(paste0(paste("out","R=1",covMatrix,"hivar",var1[j],var2[j],sep="_"),".RDS"))
    z2 = readRDS(paste0(paste("out","R=1",covMatrix,"medvar",var1[j],var2[j],sep="_"),".RDS"))
    z3 = readRDS(paste0(paste("out","R=1",covMatrix,"lowvar",var1[j],var2[j],sep="_"),".RDS"))
    z <- vector(mode='list',length=3)
    lapply(1:length(z1), function(i) c(max(z1[[i]][[1]]$I), ## peak epidemic size
                                       sum(sort(z1[[i]][[2]]$numInf,decreasing=TRUE)[1:round(0.2*nrow(z1[[i]][[2]]))])/sum(z1[[i]][[2]]$numInf), ## fraction of cases caused by top 20% of spreaders
                                       ifelse(max(z1[[i]][[1]]$t)>99, glm.nb(z1[[i]][[2]]$numInf~1)$theta, NA), ## dispersion parameter of negative binomial distribution fit to numInf
                                       ifelse(length(z1[[i]][[1]]$t)<99, 1, 0),
                                       i)) %>% 
      do.call("rbind.data.frame",.) %>%
      mutate(., Variance="high",cov=covMatrix, traits=paste(var1[j],var2[j],sep="-")) -> z[[1]]
    colnames(z[[1]])[1:5] <- c("peak","top20SS","disp","fadeout", "rep")
    
    lapply(1:length(z2), function(i) c(max(z2[[i]][[1]]$I), ## peak epidemic size
                                       sum(sort(z2[[i]][[2]]$numInf,decreasing=TRUE)[1:round(0.2*nrow(z2[[i]][[2]]))])/sum(z2[[i]][[2]]$numInf), ## fraction of cases caused by top 20% of spreaders
                                       ifelse(max(z2[[i]][[1]]$t)>99, glm.nb(z2[[i]][[2]]$numInf~1)$theta, NA), ## dispersion parameter of negative binomial distribution fit to numInf
                                       ifelse(length(z2[[i]][[1]]$t)<99, 1, 0),
                                       i)) %>% 
      do.call("rbind.data.frame",.) %>%
      mutate(., Variance="med",cov=covMatrix, traits=paste(var1[j],var2[j],sep="-")) -> z[[2]]
    colnames(z[[2]])[1:5] <- c("peak","top20SS","disp","fadeout","rep")
    
    lapply(1:length(z3), function(i) c(max(z3[[i]][[1]]$I), ## peak epidemic size
                                       sum(sort(z3[[i]][[2]]$numInf,decreasing=TRUE)[1:round(0.2*nrow(z3[[i]][[2]]))])/sum(z3[[i]][[2]]$numInf), ## fraction of cases caused by top 20% of spreaders
                                       ifelse(max(z3[[i]][[1]]$t)>99, glm.nb(z3[[i]][[2]]$numInf~1)$theta, NA), ## dispersion parameter of negative binomial distribution fit to numInf
                                       ifelse(length(z3[[i]][[1]]$t)<99, 1, 0),
                                       i)) %>% 
      do.call("rbind.data.frame",.) %>%
      mutate(., Variance="low",cov=covMatrix, traits=paste(var1[j],var2[j],sep="-")) -> z[[3]]
    colnames(z[[3]])[1:5] <- c("peak","top20SS","disp","fadeout","rep")
    
    data[[i]] <- do.call("rbind.data.frame",z)
    i <- i+1
  }
}
data %<>% do.call("rbind.data.frame",.) 

data %>% group_by(Variance, cov, traits) %>% summarise(nFade=sum(fadeout))


# order categories 
data$traits<-factor(data$traits, levels=c("alpha-gamma", "c-shed", "c-alpha","c-gamma", "shed-alpha", "shed-gamma"))

##create a column with names rather param symbols
mutate(data, trait_name=ifelse(traits=="alpha-gamma","virulence-recovery",
                               ifelse(traits=="c-shed","contact-infectiousness", 
                                ifelse(traits=="c-alpha","contact-virulence",
                                  ifelse(traits=="c-gamma","contact-recovery",
                                    ifelse(traits=="shed-alpha","infectiousness-virulence","infectiousness-recovery")))))) %>%
  mutate(., cov=ifelse(cov=="nocorr","none",
                              ifelse(cov=="poscorr","positive","negative"))) %>%
  rename(Covariance=cov) %>%
  mutate(., Variance=ifelse(Variance=="low", "low",
                            ifelse(Variance=="med", "moderate", "high"))) -> data

data$Variance<-factor(data$Variance, levels = c("low", "moderate", "high"))
data$trait_name<-factor(data$trait_name, levels=c("virulence-recovery", "contact-infectiousness", 
                                                  "contact-virulence","contact-recovery", "infectiousness-virulence", 
                                                  "infectiousness-recovery"))


# create column with covariation signs rather than words                                           
mutate(data, cov_sign=ifelse(Covariance=="negative","(-)",
                             ifelse(Covariance=="none", "(0)", "(+)")))->data
data$cov_sign<-factor(data$cov_sign,levels=c("(-)", "(0)","(+)"))

as.factor(data$fadeout)

ggplot(data, aes(x=fadeout, fill=Variance))+
  geom_bar( position=position_dodge())+
  facet_grid(cov_sign~traits)+
  theme_bw()+
  scale_x_discrete(limits=c(0, 1))


####### SUPER SPREADING PLOTS
## need to bin extreme disp values
2 -> data$disp[data$disp>= 2]
disp_data<-subset(data, disp!="NA")

# max disp R=1: 10707.47
# max disp R=4: 34713.93
# max disp R=8: 1.202398
# super spreading via dispersion plots, all trait pairs
ggplot(data, aes(x=disp))+
  geom_density(aes(group=Variance,color=Variance,fill=Variance),lwd=.5, alpha=0.9)+
  facet_grid(cov~traits)+
  theme_bw()+
  scale_color_manual(values=c("cornflowerblue", "pink", "sienna"))+
  scale_fill_manual(values=c("cornflowerblue", "pink", "sienna"))+
  labs(y="Density", x="Dispersion (k)")

# super spreading via dispersion plots, some pairs
intra_k<-ggplot(subset(disp_data, trait_name%in%c("contact-infectiousness","virulence-recovery")), aes(x=disp))+
  geom_density(aes(group=Variance,color=Variance,fill=Variance),lwd=.5, alpha=.9)+
  facet_grid(cov_sign~trait_name)+
  theme_bw()+
  scale_color_manual(values=c("cornflowerblue", "pink", "sienna3"))+
  scale_fill_manual(values=c("cornflowerblue", "pink", "sienna3"))+
  labs(y="Density", x="Dispersion (*k*)")+
  ggtitle("Intragroup Trait Pairings")+
  theme(text = element_text(size=15),
        plot.title = element_text(hjust = 0.5),
        plot.margin = margin(t = 10,  # Top margin
                             r = 20,  # Right margin
                             b = 10,  # Bottom margin
                             l = 40))+
  geom_vline(xintercept = 0.16, linetype="dashed")+
  theme(axis.title.x = ggtext::element_markdown())

pdf("intra_k_plots.pdf", width = 10, height = 10)
print(intra_k) 
dev.off() 

inter_k_a<-ggplot(subset(disp_data, trait_name%in%c("contact-virulence", "contact-recovery","infectiousness-virulence", "infectiousness-recovery")), aes(x=disp))+
  geom_density(aes(group=Variance,color=Variance,fill=Variance),lwd=.5, alpha=0.9)+
  facet_grid(cov_sign~trait_name)+
  theme_bw()+
  scale_color_manual(values=c( "cornflowerblue", "pink", "sienna3"))+
  scale_fill_manual(values=c( "cornflowerblue", "pink", "sienna3"))+
  labs(y="Density", x="Dispersion (*k*)")+
  ggtitle("Intergroup Trait Pairings")+
  theme(text = element_text(size=15),
        plot.title = element_text(hjust = 0.5),
        plot.margin = margin(t = 10,  # Top margin
                             r = 20,  # Right margin
                             b = 10,  # Bottom margin
                             l = 40))+
  geom_vline(xintercept = 0.16, linetype="dashed")+
  theme(axis.title.x = ggtext::element_markdown())

inter_k_b<-ggplot(subset(data, trait_name%in%c("infectiousness-virulence")), aes(x=disp))+
  geom_density(aes(group=Variance,color=Variance,fill=Variance),lwd=.5, alpha=0.9)+
  facet_grid(cov_sign~trait_name)+
  theme_bw()+
  scale_color_manual(values=c( "cornflowerblue", "pink", "sienna3"))+
  scale_fill_manual(values=c( "cornflowerblue", "pink", "sienna3"))+
  labs(y="Density", x="Dispersion (*k*)")+
  ggtitle("Intergroup Trait Pairings (b)")+
  theme(text = element_text(size=15),
        plot.title = element_text(hjust = 0.5),
        plot.margin = margin(t = 10,  # Top margin
                             r = 20,  # Right margin
                             b = 10,  # Bottom margin
                             l = 40))+
  geom_vline(xintercept = 0.16, linetype="dashed")+
  theme(axis.title.x = ggtext::element_markdown())


pdf("inter_k_plots.pdf", width = 10, height = 10)
print(inter_k) 
dev.off() 

####### PEAK EPIDEMIC SIZE PLOTS
data1<-subset(data, peak>5) # remove peaks smaller than 5
data2<-subset(data, peak>10)

# making colors var and x-cov
inter_peak<-ggplot(subset(data, trait_name%in%c("contact-virulence","contact-recovery", "infectiousness-virulence", "infectiousness-recovery")),aes(x=cov_sign, y=peak,shape=Variance,fill=Variance))+
  geom_violin()  +
  facet_grid(~trait_name)+
  theme_bw()+
  labs(y="Peak", x="Covariance")+
  scale_fill_manual(values=c("cornflowerblue", "mistyrose", "sienna3"))+
  ggtitle("Intergroup Trait Pairings")+
  theme(text = element_text(size=15),
        plot.title = element_text(hjust = 0.5),
        plot.margin = margin(t = 10,  # Top margin
                             r = 20,  # Right margin
                             b = 10,  # Bottom margin
                             l = 10))

pdf("inter_peak_plots.pdf", width = 10, height = 5)
print(inter_peak) 
dev.off() 

# intragroup trait pairings
intra_peak<- ggplot(subset(data, trait_name%in%c("contact-infectiousness", "virulence-recovery")), aes(x=cov_sign, y=peak,shape=Variance, fill=Variance))+
  geom_violin()  +
  facet_grid(~trait_name)+
  theme_bw()+
  labs(y="Peak", x="Covariation")+
  scale_fill_manual(values=c("cornflowerblue", "mistyrose", "sienna3"))+
  ggtitle("Intragroup Trait Pairings")+
  theme(text = element_text(size=15),
        plot.title = element_text(hjust = 0.5),
        plot.margin = margin(t = 10,  # Top margin
                             r = 10,  # Right margin
                             b = 10,  # Bottom margin
                             l = 10))
  

pdf("intra_peak_plots.pdf", width = 10, height = 5)
print(intra_peak) 
dev.off() 

#ggarrange(intra, inter, common.legend = TRUE, ncol=2, nrow=1, legend="right")

###### EFFECT OF SUPER SPREADING ON PEAK EPIDEMIC SIZE
data3<-subset(data, peak>50) 
ggplot(subset(data3, trait_name%in%c("virulence-recovery","contact-infectiousness")), aes(x=disp, y=peak, shape=Variance, color=Covariance))+
  geom_point()+
  facet_grid(~trait_name)+
  theme_bw()+
  labs(y="Peak", x="Dispersion (k)")+
  scale_color_manual(values=c("cornflowerblue", "pink", "sienna3"))+
  theme(text = element_text(size=15),
        plot.title = element_text(hjust = 0.5),
        plot.margin = margin(t = 10,  # Top margin
                             r = 10,  # Right margin
                             b = 10,  # Bottom margin
                             l = 10))+
  geom_vline(xintercept = 0.16, linetype="dashed")

ggplot(subset(data, trait_name%in%c("contact-virulence","contact-recovery", "infectiousness-virulence", "infectiousness-recovery")), aes(x=disp, y=peak, shape=Variance, color=cov_sign))+
  geom_point()+
  facet_grid(~traits)+
  theme_bw()+
  labs(y="Peak", x="Dispersion (k)")+
  scale_color_manual(values=c("cornflowerblue", "pink", "sienna3"))+
  theme(text = element_text(size=15),
        plot.title = element_text(hjust = 0.5),
        plot.margin = margin(t = 10,  # Top margin
                             r = 10,  # Right margin
                             b = 10,  # Bottom margin
                             l = 10))+
  geom_vline(xintercept = 0.16, linetype="dashed")

dispXpeak_plot<-ggplot(subset(data, trait_name%in%c("virulence-recovery","contact-infectiousness","contact-recovery","infectiousness-recovery")), 
                       aes(x=disp, y=peak, shape=Variance, color=Covariance))+
  geom_point()+
  facet_grid(~trait_name)+
  theme_bw()+
  labs(y="Peak", x="Dispersion (*k*)")+
  scale_color_manual(values=c("cornflowerblue", "pink", "sienna3"))+
  theme(text = element_text(size=15),
        plot.title = element_text(hjust = 0.5),
        plot.margin = margin(t = 10,  # Top margin
                             r = 10,  # Right margin
                             b = 10,  # Bottom margin
                             l = 10),
        panel.spacing.x=unit(5,"mm"),
        axis.title.x = ggtext::element_markdown())+
  geom_vline(xintercept = 0.16, linetype="dashed")

pdf("dispXpeak_plot.pdf", width =10, height = 7)
print(dispXpeak_plot) 
dev.off() 

### Calculate epidemic fade out



###### INFECTED DYNAMMIC PLOTS
# format output for plotting in ggplot
## To create an equivalent plot to the one below, just in ggplot, use the following code
var1=c('c','c','c','shed','shed','alpha')
var2=c('shed','alpha','gamma','alpha','gamma','gamma')
i <- 1
data2 <- vector(mode='list', length=18)
for (j in 1:6) { ## loop over the six different covariance combinations
  for (covMatrix in c("nocorr", "negcorr", "poscorr")) {
    z1 = readRDS(paste0(paste("out","R=1",covMatrix,"hivar",var1[j],var2[j],sep="_"),".RDS"))
    z2 = readRDS(paste0(paste("out","R=1",covMatrix,"medvar",var1[j],var2[j],sep="_"),".RDS"))
    z3 = readRDS(paste0(paste("out","R=1",covMatrix,"lowvar",var1[j],var2[j],sep="_"),".RDS"))
    
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
    data2[[i]] <- do.call("rbind.data.frame",z)
    i <- i+1
  }
}
data2 <- do.call("rbind.data.frame",data2)

##create a column with names rather param symbols
mutate(data2, trait_name=ifelse(traits=="alpha-gamma","virulence-recovery",
                               ifelse(traits=="c-shed","contact-infectiousness", 
                                      ifelse(traits=="c-alpha","contact-virulence",
                                             ifelse(traits=="c-gamma","contact-recovery",
                                                    ifelse(traits=="shed-alpha","infectiousness-virulence","infectiousness-recovery")))))) %>%
  mutate(., cov=ifelse(cov=="(0) Cov","none",
                       ifelse(cov=="(+) Cov","positive","negative"))) %>%
  rename(Covariance=cov) %>%
  mutate(., Variance=ifelse(Variance=="low", "low",
                            ifelse(Variance=="med", "moderate", "high"))) -> data2

data2$Variance<-factor(data2$Variance, levels = c("low", "moderate", "high"))
data2$trait_name<-factor(data2$trait_name, levels=c("virulence-recovery", "contact-infectiousness", 
                                                  "contact-virulence","contact-recovery", "infectiousness-virulence", 
                                                  "infectiousness-recovery"))


# create column with covariation signs rather than words                                           
mutate(data2, cov_sign=ifelse(Covariance=="negative","(-)",
                             ifelse(Covariance=="none", "(0)", "(+)")))->data2
data2$cov_sign<-factor(data2$cov_sign,levels=c("(-)", "(0)","(+)"))


## To plot the dynamics of I across all trait pairs, covariance structures, and variance levels
intergroup_I<-ggplot(subset(data2, trait_name%in%c("contact-recovery","contact-virulence","infectiousness-recovery","infectiousness-virulence")), aes(x=t, y=I, group=rep)) +
  stat_summary(aes(group=Variance), geom="ribbon", fun.data=mean_cl_normal, alpha=.45, fun.args=list(conf.int=0.95)) +
  stat_summary(aes(group=Variance,colour=Variance), geom="line", fun=mean) + 
  facet_grid(trait_name~cov_sign) + 
  ylim(0,175) + 
  xlim(0,75)+
  scale_color_manual(values=c("cornflowerblue", "pink", "sienna3"))+
  labs(x = "time", y = "average # infected individuals")+
  theme_bw()+
  theme(text = element_text(size=15),
        plot.title = element_text(hjust = 0.5),
        plot.margin = margin(t = 10,  # Top margin
                             r = 20,  # Right margin
                             b = 10,  # Bottom margin
                             l = 40))+
  ggtitle("Intergroup Trait Pairings")
#theme(strip.background.y=element_rect(fill="lightgreen"))

pdf("intergroup_I.pdf", width =8, height = 9)
print(intergroup_I) 
dev.off() 

intragroup_I<-ggplot(subset(data2, trait_name%in%c("virulence-recovery", "contact-infectiousness")), aes(x=t, y=I, group=rep)) +
  stat_summary(aes(group=Variance), geom="ribbon", fun.data=mean_cl_normal, alpha=.45, fun.args=list(conf.int=0.95)) +
  stat_summary(aes(group=Variance,colour=Variance), geom="line", fun=mean) + 
  facet_grid(trait_name~cov_sign) + 
  ylim(0,175) + 
 # xlim(0,75)+
  scale_color_manual(values=c("cornflowerblue", "pink", "sienna3"))+
  labs(x = "time", y = "average # infected individuals")+
  theme_bw()+
  theme(text = element_text(size=15),
        plot.title = element_text(hjust = 0.5),
        plot.margin = margin(t = 10,  # Top margin
                             r = 20,  # Right margin
                             b = 10,  # Bottom margin
                             l = 40))+
  ggtitle("Intragroup Trait Pairings")
#theme(strip.background.y=element_rect(fill="lightblue"))

## pdf of all plots
pdf("R=1_plots.pdf", width=9, height = 7)
print(intragroup_I)
print(intergroup_I)
print(inter_k_a)
print(intra_k)
print(inter_peak)
print(intra_peak)
print(dispXpeak_plot)
dev.off() 

max(data$disp)

