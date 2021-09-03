### create distributions of each level of variation and type of covariation
library(tidyverse)
library(ggplot2)
library(gridExtra)
library(cowplot)

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
### distribution of traits
d=data.frame(picks_...)
mutate(d, r0=(shed/1+shed))

## contact-shed covariation
## high var
#no cor
data.frame(pick_individuals_multivariate(1000, traitmeans = c(hivar["c"],hivar["shed"]), 
                                        traitsds = c(hivar["sd_c"], hivar["sd_shed"]),corr=nocorr))%>% 
  mutate(., R0 = (shed/(1+shed)/(hivar["alpha"]+hivar["gamma"]+hivar["d"]))) ->c_shed_nocor_hivar

#pos cor
data.frame(pick_individuals_multivariate(1000, traitmeans = c(hivar["c"],hivar["shed"]), 
                                         traitsds = c(hivar["sd_c"], hivar["sd_shed"]),corr=poscorr))%>% 
  mutate(., R0 = (shed/(1+shed)/(hivar["alpha"]+hivar["gamma"]+hivar["d"]))) ->c_shed_poscor_hivar

#pos cor
data.frame(pick_individuals_multivariate(1000, traitmeans = c(hivar["c"],hivar["shed"]), 
                                         traitsds = c(hivar["sd_c"], hivar["sd_shed"]),corr=negcorr))%>% 
  mutate(., R0 = (shed/(1+shed)/(hivar["alpha"]+hivar["gamma"]+hivar["d"]))) ->c_shed_negcor_hivar

## med var
#no cor
data.frame(pick_individuals_multivariate(1000, traitmeans = c(medvar["c"],medvar["shed"]), 
                                         traitsds = c(medvar["sd_c"], medvar["sd_shed"]),corr=nocorr))%>% 
  mutate(., R0 = (shed/(1+shed)/(medvar["alpha"]+medvar["gamma"]+medvar["d"]))) ->c_shed_nocor_medvar

#pos cor
data.frame(pick_individuals_multivariate(1000, traitmeans = c(medvar["c"],medvar["shed"]), 
                                         traitsds = c(medvar["sd_c"], medvar["sd_shed"]),corr=poscorr))%>% 
  mutate(., R0 = (shed/(1+shed)/(medvar["alpha"]+medvar["gamma"]+medvar["d"]))) ->c_shed_poscor_medvar

#pos cor
data.frame(pick_individuals_multivariate(1000, traitmeans = c(medvar["c"],medvar["shed"]), 
                                         traitsds = c(medvar["sd_c"], medvar["sd_shed"]),corr=negcorr))%>% 
  mutate(., R0 = (shed/(1+shed)/(medvar["alpha"]+medvar["gamma"]+medvar["d"]))) ->c_shed_negcor_medvar

## low var
#no cor
data.frame(pick_individuals_multivariate(1000, traitmeans = c(lowvar["c"],lowvar["shed"]), 
                                         traitsds = c(lowvar["sd_c"], lowvar["sd_shed"]),corr=nocorr))%>% 
  mutate(., R0 = (shed/(1+shed)/(lowvar["alpha"]+lowvar["gamma"]+lowvar["d"]))) ->c_shed_nocor_lowvar

#pos cor
data.frame(pick_individuals_multivariate(1000, traitmeans = c(lowvar["c"],lowvar["shed"]), 
                                         traitsds = c(lowvar["sd_c"], lowvar["sd_shed"]),corr=poscorr))%>% 
  mutate(., R0 = (shed/(1+shed)/(lowvar["alpha"]+lowvar["gamma"]+lowvar["d"]))) ->c_shed_poscor_lowvar

#pos cor
data.frame(pick_individuals_multivariate(1000, traitmeans = c(lowvar["c"],lowvar["shed"]), 
                                         traitsds = c(lowvar["sd_c"], lowvar["sd_shed"]),corr=negcorr))%>% 
  mutate(., R0 = (shed/(1+shed)/(lowvar["alpha"]+lowvar["gamma"]+lowvar["d"]))) ->c_shed_negcor_lowvar

### plots
## histogram of each covariation with increasing variation
# no covariation, contact values
c_nocor<- ggplot(c_shed_nocor_lowvar, aes(x=c))+
  geom_histogram(color="pink", fill="pink", binwidth=.1, alpha=.4)+
  geom_histogram(data=c_shed_nocor_medvar, aes(x=c), binwidth =.1, alpha = .4, color="darkgreen",fill="darkgreen")+
  geom_histogram(data=c_shed_nocor_hivar, aes(x=c),  alpha = .2, binwidth = .1, color="darkblue",fill="darkblue")+
  labs(title="contact-shed-nocor")
c_nocor

shed_nocor<- ggplot(c_shed_nocor_lowvar, aes(x=shed))+
  geom_histogram(color="pink", fill="pink", binwidth=.1, alpha=.4)+
  geom_histogram(data=c_shed_nocor_medvar, aes(x=shed), binwidth =.1, alpha = .4, color="darkgreen",fill="darkgreen")+
  geom_histogram(data=c_shed_nocor_hivar, aes(x=shed),  alpha = .4, binwidth = .1, color="darkblue",fill="darkblue")+
  labs(title="contact-shed-nocor")

R0_nocor<-ggplot(c_shed_nocor_lowvar, aes(x=R0))+
  geom_histogram(color="pink", fill="pink", alpha= .4, binwidth = 0.1)+
  geom_histogram(data=c_shed_nocor_medvar, aes(x=R0), binwidth = 0.1, alpha = .4, color="darkgreen", fill="darkgreen")+
  geom_histogram(data=c_shed_nocor_hivar, aes(x=R0), binwidth = 0.1, alpha = .4, color="darkblue", fill="darkblue")
R0_nocor

plot_grid(c_nocor, shed_nocor,R0_nocor, ncol=3, nrow=1)

c_poscor<- ggplot(c_shed_poscor_lowvar, aes(x=c))+
  geom_histogram(color="pink", fill="pink", binwidth=.1, alpha=1)+
  geom_histogram(data=c_shed_poscor_medvar, aes(x=c), binwidth =.1, alpha = .4, color="darkgreen",fill="darkgreen")+
  geom_histogram(data=c_shed_poscor_hivar, aes(x=c),  alpha = .2, binwidth = .1, color="darkblue",fill="darkblue")+
  labs(title="contact-shed-poscor")
c_poscor

shed_poscor<- ggplot(c_shed_poscor_lowvar, aes(x=shed))+
  geom_histogram(color="pink", fill="pink", binwidth=.1, alpha=1)+
  geom_histogram(data=c_shed_poscor_medvar, aes(x=shed), binwidth =.1, alpha = .4, color="darkgreen",fill="darkgreen")+
  geom_histogram(data=c_shed_poscor_hivar, aes(x=shed),  alpha = .2, binwidth = .1, color="darkblue",fill="darkblue")+
  labs(title="contact-shed-poscor")
shed_poscor

R0_poscor<-ggplot(c_shed_poscor_lowvar, aes(x=R0))+
  geom_histogram(color="pink", fill="pink", alpha= .4, binwidth = 0.1)+
  geom_histogram(data=c_shed_poscor_medvar, aes(x=R0), binwidth = 0.1, alpha = .4, color="darkgreen", fill="darkgreen")+
  geom_histogram(data=c_shed_poscor_hivar, aes(x=R0), binwidth = 0.1, alpha = .4, color="darkblue", fill="darkblue")
R0_poscor

plot_grid(c_poscor, shed_poscor, R0_poscor, ncol=3, nrow=1)

c_negcor<- ggplot(c_shed_negcor_lowvar, aes(x=c))+
  geom_histogram(color="pink", fill="pink", binwidth=.1, alpha=1)+
  geom_histogram(data=c_shed_negcor_medvar, aes(x=c), binwidth =.1, alpha = .4, color="darkgreen",fill="darkgreen")+
  geom_histogram(data=c_shed_negcor_hivar, aes(x=c),  alpha = .2, binwidth = .1, color="darkblue",fill="darkblue")+
  labs(title="contact-shed-negcor")
c_poscor

shed_negcor<- ggplot(c_shed_negcor_lowvar, aes(x=shed))+
  geom_histogram(color="pink", fill="pink", binwidth=.1, alpha=1)+
  geom_histogram(data=c_shed_negcor_medvar, aes(x=shed), binwidth =.1, alpha = .4, color="darkgreen",fill="darkgreen")+
  geom_histogram(data=c_shed_negcor_hivar, aes(x=shed),  alpha = .2, binwidth = .1, color="darkblue",fill="darkblue")+
  labs(title="contact-shed-negcor")
shed_negcor

R0_negcor<-ggplot(c_shed_negcor_lowvar, aes(x=R0))+
  geom_histogram(color="pink", fill="pink", alpha= .4, binwidth = 0.1)+
  geom_histogram(data=c_shed_negcor_medvar, aes(x=R0), binwidth = 0.1, alpha = .4, color="darkgreen", fill="darkgreen")+
  geom_histogram(data=c_shed_negcor_hivar, aes(x=R0), binwidth = 0.1, alpha = .4, color="darkblue", fill="darkblue")
R0_negcor

plot_grid(c_nocor, shed_nocor, R0_nocor, c_poscor, shed_poscor, R0_poscor, 
          c_negcor, shed_negcor, R0_negcor, ncol=3, nrow=3)
