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
var1=c('c','c','c','shed','shed','alpha')
var2=c('shed','alpha','gamma','alpha','gamma','gamma')
contrasts <- vector(mode='list', length=18)
i <- 1
for (j in 1:6) {
  for (covMatrix in c("(0) Cov", "(-) Cov", "(+) Cov")) {
    filter(data, traits==paste(var1[j],var2[j],sep='-') & cov==covMatrix) %>% group_by(Variance, rep) %>% summarise(peak=ifelse(max(I)>50,max(I),NA)) %>%
      as.data.frame() %>%
      with(., lm(peak~Variance)) %>% emmeans(pairwise~Variance)-> comparison ## see https://timmastny.com/blog/tests-pairwise-categorical-mean-emmeans-contrast/
    comparison$contrasts %>% as.data.frame %>% mutate(., Covariance=covMatrix, traits=paste(var1[j],var2[j],sep='-')) -> contrast
    contrasts[[i]] <- contrast
    i <- i + 1
  }
}
contrasts %>% do.call("rbind.data.frame", .) -> contrasts2
write.csv(contrasts2, file="Var_effects_on_peak_stats.csv")


contrasts <- vector(mode='list', length=18)
i <- 1
for (j in 1:6) {
  for (var in c("high", "med", "low")) {
    filter(data, traits==paste(var1[j],var2[j],sep='-') & Variance==var) %>% group_by(cov, rep) %>% summarise(peak=ifelse(max(I)>50,max(I),NA)) %>%
      as.data.frame() -> summ
    summ$cov <- factor(sapply(summ$cov, function(c) switch(c, "(-) Cov"="neg", "(0) Cov"="zero", "(+) Cov"="pos")) %>% unname, levels=c("pos","zero","neg"))
    with(summ, lm(peak~cov)) %>% emmeans(pairwise~cov)-> comparison ## see https://timmastny.com/blog/tests-pairwise-categorical-mean-emmeans-contrast/
    comparison$contrasts %>% as.data.frame %>% mutate(., Variance=var, traits=paste(var1[j],var2[j],sep='-')) -> contrast
    contrasts[[i]] <- contrast
    i <- i + 1
  }
}
contrasts %>% do.call("rbind.data.frame", .) -> contrasts2
write.csv(contrasts2, file="Cov_effects_on_peak_stats.csv")

## Stats on the fraction of cases caused by the 20% of spreaders
var1=c('c','c','c','shed','shed','alpha')
var2=c('shed','alpha','gamma','alpha','gamma','gamma')
i <- 1
data <- vector(mode='list', length=18)
for (j in 1:6) { ## loop over the six different covariance combinations
  for (covMatrix in c("nocorr", "negcorr", "poscorr")) {
    z1 = readRDS(paste0(paste("out",covMatrix,"hivar",var1[j],var2[j],sep="_"),".RDS"))
    z2 = readRDS(paste0(paste("out",covMatrix,"medvar",var1[j],var2[j],sep="_"),".RDS"))
    z3 = readRDS(paste0(paste("out",covMatrix,"lowvar",var1[j],var2[j],sep="_"),".RDS"))
    
