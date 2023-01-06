nocorr <- matrix(c(1,0,0,1), nrow=2, byrow=T)
negcorr <- matrix(c(1,-.5,-.5,1), nrow=2, byrow=T)
poscorr <- matrix(c(1,.5,.5,1), nrow=2, byrow=T)

hivar = c(   c=0.1,    shed=0.05,    alpha=0.1,    gamma=0.1, 
             sd_c=0.5, sd_shed=0.25, sd_alpha=0.5, sd_gamma=0.5, 
             b=2.5, d=.1, bs=.01)
medvar = c(   c=0.1,    shed=0.05,    alpha=0.1,    gamma=0.1, 
              sd_c=0.1, sd_shed=0.05, sd_alpha=0.1, sd_gamma=0.1, 
              b=2.5, d=.1, bs=.01)
lowvar = c(   c=0.1,     shed=0.05,    alpha=0.1,     gamma=0.1, 
              sd_c=0.02, sd_shed=0.01, sd_alpha=0.02, sd_gamma=0.02, 
              b=2.5, d=.1, bs=.01)

## Statistical analysis
data <- readRDS("population_dynamics_all_simulations_medvar.RDS")

library(emmeans)
library(tidyverse)
library(MASS)
var1=c('c','c','c','shed','shed','alpha')
var2=c('shed','alpha','gamma','alpha','gamma','gamma')
contrasts1 <- vector(mode='list', length=18)
contrasts2 <- vector(mode='list', length=18)
i <- 1
for (j in 1:6) {
  for (covMatrix in c("(0) Cov", "(-) Cov", "(+) Cov")) {
    filter(data, traits==paste(var1[j],var2[j],sep='-') & cov==covMatrix) %>% group_by(Variance, rep) %>% summarise(peak=ifelse(max(I)>50,max(I),NA), Iequil=ifelse(min(I)>1,mean(tail(I),50),NA)) %>%
      as.data.frame() -> d 
    with(d, lm(log(peak)~Variance)) -> model1
    with(d, lm(Iequil~Variance)) -> model2
    normal.test1 <- shapiro.test(model1$residuals)$p.value
    normal.test2 <- shapiro.test(model2$residuals)$p.value
    model1 %>% emmeans(pairwise~Variance) -> comparison1 ## see https://timmastny.com/blog/tests-pairwise-categorical-mean-emmeans-contrast/
    model2 %>% emmeans(pairwise~Variance) -> comparison2 ## see https://timmastny.com/blog/tests-pairwise-categorical-mean-emmeans-contrast/
    comparison1$contrasts %>% as.data.frame %>% mutate(., estimate=exp(estimate), normalTest=normal.test1, Covariance=covMatrix, traits=paste(var1[j],var2[j],sep='-')) -> contrast1
    comparison2$contrasts %>% as.data.frame %>% mutate(., estimate=exp(estimate), normalTest=normal.test2, Covariance=covMatrix, traits=paste(var1[j],var2[j],sep='-')) -> contrast2
    contrasts1[[i]] <- contrast1
    contrasts2[[i]] <- contrast2
    i <- i + 1
    
  }
}
contrasts1 %>% do.call("rbind.data.frame", .) -> contrasts1
contrasts2 %>% do.call("rbind.data.frame", .) -> contrasts2
write.csv(contrasts1, file="Var_effects_on_peak_stats.csv")
write.csv(contrasts2, file="Var_effects_on_equil_stats.csv")


contrasts <- vector(mode='list', length=18)
i <- 1
for (j in 1:6) {
  for (var in c("high", "med", "low")) {
    filter(data, traits==paste(var1[j],var2[j],sep='-') & Variance==var) %>% group_by(cov, rep) %>% summarise(peak=ifelse(max(I)>50,max(I),NA)) %>%
      as.data.frame() -> summ
    summ$cov <- factor(sapply(summ$cov, function(c) switch(c, "(-) Cov"="neg", "(0) Cov"="zero", "(+) Cov"="pos")) %>% unname, levels=c("pos","zero","neg"))
    with(summ, lm(log(peak)~cov)) -> model
    normal.test <- shapiro.test(model$residuals)$p.value
    model %>% emmeans(pairwise~cov) -> comparison ## see https://timmastny.com/blog/tests-pairwise-categorical-mean-emmeans-contrast/
    comparison$contrasts %>% as.data.frame %>% mutate(., normalTest=normal.test, Variance=var, traits=paste(var1[j],var2[j],sep='-')) -> contrast
    contrasts[[i]] <- contrast
    i <- i + 1
  }
}
contrasts %>% do.call("rbind.data.frame", .) -> contrasts2
write.csv(contrasts2, file="Cov_effects_on_peak_stats.csv")

## Stats on super spreading via dispersion (K)
var1=c('c','c','c','shed','shed','alpha')
var2=c('shed','alpha','gamma','alpha','gamma','gamma')
i <- 1
data <- vector(mode='list', length=18)
for (j in 1:6) { ## loop over the six different covariance combinations
  for (covMatrix in c("nocorr", "negcorr", "poscorr")) {
    z1 = readRDS(paste0(paste("out",covMatrix,"hivar",var1[j],var2[j],sep="_"),".RDS"))
    z2 = readRDS(paste0(paste("out",covMatrix,"medvar",var1[j],var2[j],sep="_"),".RDS"))
    z3 = readRDS(paste0(paste("out",covMatrix,"lowvar",var1[j],var2[j],sep="_"),".RDS"))
    z <- vector(mode='list',length=3)
    lapply(1:length(z1), function(i) c(max(z1[[i]][[1]]$I), ## peak epidemic size
                                       sum(sort(z1[[i]][[2]]$numInf,decreasing=TRUE)[1:round(0.2*nrow(z1[[i]][[2]]))])/sum(z1[[i]][[2]]$numInf), ## fraction of cases caused by top 20% of spreaders
                                       ifelse(max(z1[[i]][[1]]$t)>99, glm.nb(z1[[i]][[2]]$numInf~1)$theta, NA), ## dispersion parameter of negative binomial distribution fit to numInf
                                       i)) %>% 
      do.call("rbind.data.frame",.) %>%
      mutate(., Variance="high",cov=covMatrix, traits=paste(var1[j],var2[j],sep="-")) -> z[[1]]
    colnames(z[[1]])[1:4] <- c("peak","top20SS","disp","rep")
    
    lapply(1:length(z2), function(i) c(max(z2[[i]][[1]]$I), ## peak epidemic size
                                       sum(sort(z2[[i]][[2]]$numInf,decreasing=TRUE)[1:round(0.2*nrow(z2[[i]][[2]]))])/sum(z2[[i]][[2]]$numInf), ## fraction of cases caused by top 20% of spreaders
                                       ifelse(max(z2[[i]][[1]]$t)>99, glm.nb(z2[[i]][[2]]$numInf~1)$theta, NA), ## dispersion parameter of negative binomial distribution fit to numInf
                                       i)) %>% 
      do.call("rbind.data.frame",.) %>%
      mutate(., Variance="med",cov=covMatrix, traits=paste(var1[j],var2[j],sep="-")) -> z[[2]]
    colnames(z[[2]])[1:4] <- c("peak","top20SS","disp","rep")
    
    lapply(1:length(z3), function(i) c(max(z3[[i]][[1]]$I), ## peak epidemic size
                                       sum(sort(z3[[i]][[2]]$numInf,decreasing=TRUE)[1:round(0.2*nrow(z3[[i]][[2]]))])/sum(z3[[i]][[2]]$numInf), ## fraction of cases caused by top 20% of spreaders
                                       ifelse(max(z3[[i]][[1]]$t)>99, glm.nb(z3[[i]][[2]]$numInf~1)$theta, NA), ## dispersion parameter of negative binomial distribution fit to numInf
                                       i)) %>% 
      do.call("rbind.data.frame",.) %>%
      mutate(., Variance="low",cov=covMatrix, traits=paste(var1[j],var2[j],sep="-")) -> z[[3]]
    colnames(z[[3]])[1:4] <- c("peak","top20SS","disp","rep")
    
    data[[i]] <- do.call("rbind.data.frame",z)
    i <- i+1
  }
}
data %>% do.call("rbind.data.frame",.) -> data

data$Variance<-factor(data$Variance, levels = c("low", "med", "high"))
data$traits<-factor(data$traits, levels=c("alpha-gamma", "c-shed", "c-alpha","c-gamma", "shed-alpha", "shed-gamma"))

ggplot(data, aes(x=disp))+
  geom_density(aes(group=Variance,color=Variance,fill=Variance),lwd=.5, alpha=0.45)+
  facet_grid(cov~traits)+
  theme_bw()+
  scale_color_manual(values=c("pink", "darkblue", "darkgreen"))+
  scale_fill_manual(values=c("pink", "darkblue", "darkgreen"))+
  labs(y="Density", x="Dispersion (K)")


## plot only the moderate variation dispersion plots 
var_mod_data<-subset(data, Variance%in%c("med"))
ggplot(subset(var_mod_data, traits%in%c("alpha-gamma","c-shed", "c-alpha", "c-gamma")), aes(x=disp))+
  geom_density(aes(group=Variance,color=Variance,fill=Variance),lwd=.5, alpha=0.45)+
  facet_grid(cov~traits)+
  theme_bw()+
  scale_color_manual(values=("darkblue"))+
  scale_fill_manual(values=("darkblue"))+
  labs(y="Density", x="Dispersion (K)")

ggplot(subset(data, traits%in%c("alpha-gamma","c-shed", "shed-alpha", "shed-gamma")), aes(x=disp))+
  geom_density(aes(group=Variance,color=Variance,fill=Variance),lwd=.5, alpha=0.45)+
  facet_grid(cov~traits)+
  theme_bw()+
  scale_color_manual(values=c("pink", "darkblue", "darkgreen"))+
  scale_fill_manual(values=c("pink", "darkblue", "darkgreen"))+
  labs(y="Density", x="Dispersion (K)")

ggplot(subset(data, traits%in%c("alpha-gamma","c-shed", "shed-alpha", "shed-gamma")), aes(x=peak))+
  geom_density(aes(group=Variance,color=Variance,fill=Variance),lwd=.5, alpha=0.45)+
  facet_grid(cov~traits)+
  theme_bw()+
  scale_color_manual(values=c("pink", "darkblue", "darkgreen"))+
  scale_fill_manual(values=c("pink", "darkblue", "darkgreen"))+
  labs(y="Density", x="Peak")



ggplot(subset(data, traits%in%c("alpha-gamma","c-shed", "shed-alpha", "shed-gamma")), aes(x=disp, y=peak, shape=Variance, color=cov))+
  geom_point()+
  facet_grid(~traits)+
  theme_bw()+
  scale_color_manual(values=c("pink", "darkblue", "darkgreen"))

ggplot(subset(data, traits%in%c("alpha-gamma","c-shed", "c-alpha", "c-gamma")), aes(x=Variance, y=disp,shape=cov,colour=cov))+
  geom_boxplot()  +
  facet_grid(~traits)+
  theme_bw()+
  scale_color_manual(values=c("pink", "darkblue", "darkgreen"))

data<-subset(data, peak>5)

ggplot(subset(data, traits%in%c("alpha-gamma","c-shed", "c-alpha", "c-gamma")), aes(x=Variance, y=peak,shape=cov,colour=cov, fill=cov))+
  geom_violin()  +
  facet_grid(~traits)+
  theme_bw()+
  scale_color_manual(values=c("pink", "darkblue", "darkgreen"))
  

ggplot(subset(data, traits%in%c("c-gamma", "shed-gamma")), aes(x=Variance, y=peak,shape=cov,fill=cov))+
  geom_violin()  +
  facet_grid(~traits)+
  theme_bw()+
  scale_fill_manual(values=c("pink", "darkblue", "darkgreen"))

ggplot(subset(data, traits%in%c("alpha-gamma","c-shed")), aes(x=Variance, y=peak,shape=cov,colour=cov))+
  geom_violin()  +
  facet_grid(~traits)+
  theme_bw()+
  scale_color_manual(values=c("pink", "darkblue", "darkgreen"))


contrasts <- vector(mode='list', length=18)
i <- 1
for (j in 1:6) {
  for (covMatrix in c("nocorr", "negcorr", "poscorr")) {
    filter(data, traits==paste(var1[j],var2[j],sep='-') & cov==covMatrix) %>% 
      as.data.frame() %>%
      with(., lm(top20SS~Variance)) %>% emmeans(pairwise~Variance)-> comparison ## see https://timmastny.com/blog/tests-pairwise-categorical-mean-emmeans-contrast/
    comparison$contrasts %>% as.data.frame %>% mutate(., Covariance=covMatrix, traits=paste(var1[j],var2[j],sep='-')) -> contrast
    contrasts[[i]] <- contrast
    i <- i + 1
  }
}
contrasts %>% do.call("rbind.data.frame", .) -> contrasts2
write.csv(contrasts2, file="Var_effects_on_SS_stats.csv")

contrasts <- vector(mode='list', length=18)
i <- 1
for (j in 1:6) {
  for (var in c("high", "med", "low")) {
    filter(data, traits==paste(var1[j],var2[j],sep='-') & Variance==var) %>% 
      as.data.frame() -> summ
    summ$cov <- factor(sapply(summ$cov, function(c) switch(c, "negcorr"="neg", "nocorr"="zero", "poscorr"="pos")) %>% unname, levels=c("pos","zero","neg"))
    with(summ, lm(top20SS~cov)) %>% emmeans(pairwise~cov)-> comparison ## see https://timmastny.com/blog/tests-pairwise-categorical-mean-emmeans-contrast/
    comparison$contrasts %>% as.data.frame %>% mutate(., Variance=var, traits=paste(var1[j],var2[j],sep='-')) -> contrast
    contrasts[[i]] <- contrast
    i <- i + 1
  }
}
contrasts %>% do.call("rbind.data.frame", .) -> contrasts2
write.csv(contrasts2, file="Cov_effects_on_SS_stats.csv")


## Is there any effect of proportion of super-spreaders on peak epidemic size?
## For each trait pairing, you could just regress peak epidemic size against proportion of super-spreaders
glm(peak~top20SS,data=subset(data,traits=="c-shed" & peak > 50)) %>% summary ## significant
glm(peak~top20SS,data=subset(data,traits=="c-alpha" & peak > 50)) %>% summary ## significant
glm(peak~top20SS,data=subset(data,traits=="c-gamma" & peak > 50)) %>% summary ## significant
glm(peak~top20SS,data=subset(data,traits=="shed-alpha" & peak > 50)) %>% summary ## NOT significant
glm(peak~top20SS,data=subset(data,traits=="shed-gamma" & peak > 50)) %>% summary ## significant
glm(peak~top20SS,data=subset(data,traits=="alpha-gamma" & peak > 50)) %>% summary ## significant

## Of course, there is ALWAYS a really strong correlation between top20SS and variance, so it could be that once you've accounted for the effect of variance on top20SS, there is no effect of top20SS on the peak
aov(top20SS~Variance,data=subset(data,traits=="c-shed" & peak > 50)) %>% summary ## significant
aov(top20SS~Variance,data=subset(data,traits=="c-alpha" & peak > 50)) %>% summary ## significant
aov(top20SS~Variance,data=subset(data,traits=="c-gamma" & peak > 50)) %>% summary ## significant
aov(top20SS~Variance,data=subset(data,traits=="shed-alpha" & peak > 50)) %>% summary ## significant
aov(top20SS~Variance,data=subset(data,traits=="shed-gamma" & peak > 50)) %>% summary ## significant
aov(top20SS~Variance,data=subset(data,traits=="alpha-gamma" & peak > 50)) %>% summary ## significant

## c-shed influence of top20SS on peak after accounting for variance
glm(peak~top20SS,data=subset(data,traits=="c-shed" & peak > 50 & Variance=="low")) %>% summary ## significant
glm(peak~top20SS,data=subset(data,traits=="c-shed" & peak > 50 & Variance=="med")) %>% summary ## significant
glm(peak~top20SS,data=subset(data,traits=="c-shed" & peak > 50 & Variance=="high")) %>% summary ## NOT SIGNIFICANT
glm(peak~top20SS+Variance,data=subset(data,traits=="c-shed" & peak > 50)) %>% summary ## significant effect of SS after accounting for variance

aov(top20SS~cov,data=subset(data,traits=="c-shed" & peak > 50 & Variance=="low")) %>% summary
aov(top20SS~cov,data=subset(data,traits=="c-shed" & peak > 50 & Variance=="med")) %>% summary
aov(top20SS~cov,data=subset(data,traits=="c-shed" & peak > 50 & Variance=="high")) %>% summary
glm(peak~top20SS+Variance+cov,data=subset(data,traits=="c-shed" & peak > 50)) %>% summary ## significant effect of SS after accounting for variance and covariance

## c-alpha influence of top20SS on peak after accounting for variance
glm(peak~top20SS,data=subset(data,traits=="c-alpha" & peak > 50 & Variance=="low")) %>% summary ## significant
glm(peak~top20SS,data=subset(data,traits=="c-alpha" & peak > 50 & Variance=="med")) %>% summary ## significant
glm(peak~top20SS,data=subset(data,traits=="c-alpha" & peak > 50 & Variance=="high")) %>% summary ## NOT SIGNIFICANT
#glm(peak~top20SS+Variance-1,data=subset(data,traits=="c-alpha" & peak > 50)) %>% summary ## significant effect of SS after accounting for variance

## of course, this could just be that there is a strong relationship between SS and covariance! Need to check this too.
aov(top20SS~cov,data=subset(data,traits=="c-alpha" & peak > 50 & Variance=="low")) %>% summary ## significant
aov(top20SS~cov,data=subset(data,traits=="c-alpha" & peak > 50 & Variance=="med")) %>% summary ## significant 
aov(top20SS~cov,data=subset(data,traits=="c-alpha" & peak > 50 & Variance=="high")) %>% summary ## significant
#glm(peak~top20SS+Variance+cov,data=subset(data,traits=="c-alpha" & peak > 50)) %>% summary ## significant effect of SS after accounting for variance and covariance

## Must perform separate regressions after accounting for both variance and covariance, since these have a significant effect on top20SS
peakSSreg <- expand.grid(cov=c("nocorr","negcorr","poscorr"),
                         Variance=c("low","med","high"),
                         traits=c("c-shed","c-alpha","c-gamma","shed-alpha","shed-gamma","alpha-gamma"),
                         slope=NA,
                         signif=NA,
                         normal=NA
                         ) %>% as.data.frame
for (i in 1:nrow(peakSSreg)) {
  mod <- try(lm(peak~top20SS,data=subset(data,traits==peakSSreg$traits[i] & peak > 50 & Variance==peakSSreg$Variance[i] & cov==peakSSreg$cov[i])))
  if(!inherits(mod,'try-error')) { ## if all the simulations went extinct, don't try to do a regression
    peakSSreg$slope[i] <- mod$coefficients[2]
    if (!is.na(mod$coefficients[2])) { ## if there was no variation in top20SS then don't try to do a regression
      peakSSreg$signif[i] <- signif(summary(mod)$coefficients[2,4],3)
      peakSSreg$normal[i] <- ifelse(signif(shapiro.test(mod$residuals)$p.value,3) > 0.05, TRUE, FALSE)
    }
  }
}

