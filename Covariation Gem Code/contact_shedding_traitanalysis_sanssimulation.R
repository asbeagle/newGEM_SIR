## contact and shedding
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

R0 = 235*((((hivar["shed"]/(1+hivar["shed"]))*hivar["c"])/(hivar["alpha"]+hivar["gamma"]+hivar["d"])))

### distribution of traits
d=data.frame(picks_...)
mutate(d, r0=(shed/1+shed))

## contact-shed covariation
## high var
#no cor
data.frame(pick_individuals_multivariate(1000, traitmeans = c(hivar["c"],hivar["shed"]), 
                                        traitsds = c(hivar["sd_c"], hivar["sd_shed"]),corr=nocorr))%>% 
  mutate(., R0 = (235*(((shed/(1+shed))*c)/(hivar["alpha"]+hivar["gamma"]+hivar["d"])))) ->c_shed_nocor_hivar

c_shed_nocor_hivar_mean_R0<-mean(c_shed_nocor_hivar[,3])

#pos cor
data.frame(pick_individuals_multivariate(1000, traitmeans = c(hivar["c"],hivar["shed"]), 
                                         traitsds = c(hivar["sd_c"], hivar["sd_shed"]),corr=poscorr))%>% 
  mutate(., R0 = (235*(((shed/(1+shed))*c)/(hivar["alpha"]+hivar["gamma"]+hivar["d"])))) ->c_shed_poscor_hivar

c_shed_poscor_hivar_mean_R0<-mean(c_shed_poscor_hivar[,3])

#pos cor
data.frame(pick_individuals_multivariate(1000, traitmeans = c(hivar["c"],hivar["shed"]), 
                                         traitsds = c(hivar["sd_c"], hivar["sd_shed"]),corr=negcorr))%>% 
  mutate(., R0 = (235*(((shed/(1+shed))*c)/(hivar["alpha"]+hivar["gamma"]+hivar["d"])))) ->c_shed_negcor_hivar

c_shed_negcor_hivar_mean_R0<-mean(c_shed_negcor_hivar[,3])

## med var
#no cor
data.frame(pick_individuals_multivariate(1000, traitmeans = c(medvar["c"],medvar["shed"]), 
                                         traitsds = c(medvar["sd_c"], medvar["sd_shed"]),corr=nocorr))%>% 
  mutate(., R0 = (235*(((shed/(1+shed))*c)/(hivar["alpha"]+hivar["gamma"]+hivar["d"])))) ->c_shed_nocor_medvar

c_shed_nocor_medvar_mean_R0<-mean(c_shed_nocor_medvar[,3])

#pos cor
data.frame(pick_individuals_multivariate(1000, traitmeans = c(medvar["c"],medvar["shed"]), 
                                         traitsds = c(medvar["sd_c"], medvar["sd_shed"]),corr=poscorr))%>% 
  mutate(., R0 = (235*(((shed/(1+shed))*c)/(hivar["alpha"]+hivar["gamma"]+hivar["d"])))) ->c_shed_poscor_medvar

c_shed_poscor_medvar_mean_R0<-mean(c_shed_poscor_medvar[,3])

# neg cor
data.frame(pick_individuals_multivariate(1000, traitmeans = c(medvar["c"],medvar["shed"]), 
                                         traitsds = c(medvar["sd_c"], medvar["sd_shed"]),corr=negcorr))%>% 
  mutate(., R0 = (235*(((shed/(1+shed))*c)/(hivar["alpha"]+hivar["gamma"]+hivar["d"])))) ->c_shed_negcor_medvar

c_shed_negcor_medvar_mean_R0<-mean(c_shed_negcor_medvar[,3])

## low var
#no cor
data.frame(pick_individuals_multivariate(1000, traitmeans = c(lowvar["c"],lowvar["shed"]), 
                                         traitsds = c(lowvar["sd_c"], lowvar["sd_shed"]),corr=nocorr))%>% 
  mutate(., R0 = (235*(((shed/(1+shed))*c)/(hivar["alpha"]+hivar["gamma"]+hivar["d"])))) ->c_shed_nocor_lowvar

c_shed_nocor_lowvar_mean_R0<-mean(c_shed_nocor_lowvar[,3])

#pos cor
data.frame(pick_individuals_multivariate(1000, traitmeans = c(lowvar["c"],lowvar["shed"]), 
                                         traitsds = c(lowvar["sd_c"], lowvar["sd_shed"]),corr=poscorr))%>% 
  mutate(., R0 = (235*(((shed/(1+shed))*c)/(hivar["alpha"]+hivar["gamma"]+hivar["d"])))) ->c_shed_poscor_lowvar

c_shed_poscor_lowvar_mean_R0<-mean(c_shed_poscor_lowvar[,3])

#pos cor
data.frame(pick_individuals_multivariate(1000, traitmeans = c(lowvar["c"],lowvar["shed"]), 
                                         traitsds = c(lowvar["sd_c"], lowvar["sd_shed"]),corr=negcorr))%>% 
  mutate(., R0 = (235*(((shed/(1+shed))*c)/(hivar["alpha"]+hivar["gamma"]+hivar["d"])))) ->c_shed_negcor_lowvar

c_shed_negcor_lowvar_mean_R0<-mean(c_shed_negcor_lowvar[,3])

### plots
library(ggpubr)
## bin extreme values
# no cov
c_shed_nocor_hivar$R0[c_shed_nocor_hivar$R0 >= 50] <-50
c_shed_nocor_medvar$R0[c_shed_nocor_medvar$R0 >= 50] <-50

# pos cov
c_shed_poscor_hivar$R0[c_shed_poscor_hivar$R0 >= 50] <- 50
c_shed_poscor_medvar$R0[c_shed_poscor_medvar$R0 >= 50] <-50

# neg cov
c_shed_poscor_hivar$R0[c_shed_poscor_hivar$R0 >= 50] <- 50
c_shed_poscor_medvar$R0[c_shed_poscor_medvar$R0 >= 50] <-50

## R0 distrubutions
colors<- c("Hi Var" = "darkblue", "Med Var"="darkgreen", "Low Var" = "pink")
# no cov
R0_nocor<-ggplot(c_shed_nocor_lowvar, aes(x=R0))+
  geom_histogram(aes(color="Low Var"), fill = "pink", alpha= .4, binwidth = 2)+
  geom_histogram(data=c_shed_nocor_medvar, aes(x=R0, color = "Med Var"), fill="darkgreen", binwidth = 2, alpha = .4)+
  geom_histogram(data=c_shed_nocor_hivar, aes(x=R0, color = "Hi Var"),fill = "darkblue", binwidth = 2, alpha = .4)+
  geom_vline(xintercept = R0, linetype="dotdash", color = "red", size=.5)+
  geom_vline(xintercept = c_shed_nocor_hivar_mean_R0, linetype="dotdash", color = "darkblue", size=.5)+
  geom_vline(xintercept = c_shed_nocor_lowvar_mean_R0, linetype="dotdash", color = "pink", size=.5)+
  geom_vline(xintercept = c_shed_nocor_medvar_mean_R0, linetype="dotdash", color = "darkgreen", size=.5)+
  labs(title="c-shed-nocor",
       color = "Legend")+
  scale_x_continuous(breaks=c(seq(0, 50, by=5)), labels=c(seq(0,45, by=5),"50+"))+
  ylim(0, 800)+
  scale_color_manual(values = c("darkblue", "pink", "darkgreen"))+
  guides(color = guide_legend(override.aes = list(fill = c("darkblue", "pink", "darkgreen"))))
R0_nocor

# pos cov
R0_poscor<-ggplot(c_shed_poscor_lowvar, aes(x=R0))+
  geom_histogram(color="pink", fill="pink", alpha= .4, binwidth = 2)+
  geom_histogram(data=c_shed_poscor_medvar, aes(x=R0),binwidth = 2, alpha = .4, color="darkgreen", fill="darkgreen")+
  geom_histogram(data=c_shed_poscor_hivar, aes(x=R0), binwidth = 2, alpha = .4, color="darkblue", fill="darkblue")+
  geom_vline(xintercept = R0, linetype="dotdash", color = "red", size=.5)+
  geom_vline(xintercept = c_shed_poscor_hivar_mean_R0, linetype="dotdash", color = "darkblue", size=.5)+
  geom_vline(xintercept = c_shed_poscor_lowvar_mean_R0, linetype="dotdash", color = "pink", size=.5)+
  geom_vline(xintercept = c_shed_poscor_medvar_mean_R0, linetype="dotdash", color = "darkgreen", size=.5)+
  labs(title="c-shed-poscor")+
  scale_x_continuous(breaks=c(seq(0, 50, by=5)), labels=c(seq(0,45, by=5),"50+"))+
  ylim(0, 800)

# neg cov
R0_negcor<-ggplot(c_shed_negcor_lowvar, aes(x=R0))+
  geom_histogram(color="pink", fill="pink", alpha= .4, binwidth = 2)+
  geom_histogram(data=c_shed_negcor_medvar, aes(x=R0), binwidth = 2, alpha = .4, color="darkgreen", fill="darkgreen")+
  geom_histogram(data=c_shed_negcor_hivar, aes(x=R0), binwidth = 2, alpha = .4, color="darkblue", fill="darkblue")+
  geom_vline(xintercept = R0, linetype="dotdash", color = "red", size=.5)+
  geom_vline(xintercept = c_shed_negcor_hivar_mean_R0, linetype="dotdash", color = "darkblue", size=.5)+
  geom_vline(xintercept = c_shed_negcor_lowvar_mean_R0, linetype="dotdash", color = "pink", size=.5)+
  geom_vline(xintercept = c_shed_negcor_medvar_mean_R0, linetype="dotdash", color = "darkgreen", size=.5)+
  labs(title="c-shed-negcor")+
  scale_x_continuous(breaks=c(seq(0, 50, by=5)), labels=c(seq(0,45, by=5),"50+"))+
  ylim(0, 800)

ggarrange(R0_nocor, R0_poscor, R0_negcor, ncol=3, nrow=1, common.legend = TRUE, legend="right")

### distrubution of all trait values
## histogram of each covariation with increasing variation
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

c_shed_nocor_hivar$R0[c_shed_nocor_hivar$R0 >= 50] <-50
c_shed_nocor_medvar$R0[c_shed_nocor_medvar$R0 >= 50] <-50

R0_nocor<-ggplot(c_shed_nocor_lowvar, aes(x=R0))+
  geom_histogram(color="pink", fill="pink", alpha= .4, binwidth = 2)+
  geom_histogram(data=c_shed_nocor_medvar, aes(x=R0), binwidth = 2, alpha = .4, color="darkgreen", fill="darkgreen")+
  geom_histogram(data=c_shed_nocor_hivar, aes(x=R0), binwidth = 2, alpha = .4, color="darkblue", fill="darkblue")+
  geom_vline(xintercept = R0, linetype="dotdash", color = "red", size=.5)+
  geom_vline(xintercept = c_shed_nocor_hivar_mean_R0, linetype="dotdash", color = "darkblue", size=.5)+
  geom_vline(xintercept = c_shed_nocor_lowvar_mean_R0, linetype="dotdash", color = "pink", size=.5)+
  geom_vline(xintercept = c_shed_nocor_medvar_mean_R0, linetype="dotdash", color = "darkgreen", size=.5)+
  labs(title="c-shed-nocor")+
  scale_x_continuous(breaks=c(seq(0, 50, by=5)), labels=c(seq(0,45, by=5),"50+"))+
  ylim(0, 800)
R0_nocor

plot_grid(c_nocor, shed_nocor,R0_nocor, ncol=3, nrow=1)

c_poscor<- ggplot(c_shed_poscor_lowvar, aes(x=c))+
  geom_histogram(color="pink", fill="pink", binwidth=.1, alpha=.4)+
  geom_histogram(data=c_shed_poscor_medvar, aes(x=c), binwidth =.1, alpha = .4, color="darkgreen",fill="darkgreen")+
  geom_histogram(data=c_shed_poscor_hivar, aes(x=c),  alpha = .2, binwidth = .1, color="darkblue",fill="darkblue")+
  labs(title="contact-shed-poscor")
c_poscor

shed_poscor<- ggplot(c_shed_poscor_lowvar, aes(x=shed))+
  geom_histogram(color="pink", fill="pink", binwidth=.1, alpha=.4)+
  geom_histogram(data=c_shed_poscor_medvar, aes(x=shed), binwidth =.1, alpha = .4, color="darkgreen",fill="darkgreen")+
  geom_histogram(data=c_shed_poscor_hivar, aes(x=shed),  alpha = .4, binwidth = .1, color="darkblue",fill="darkblue")+
  labs(title="contact-shed-poscor")
shed_poscor


c_shed_poscor_hivar$R0[c_shed_poscor_hivar$R0 >= 50] <- 50
c_shed_poscor_medvar$R0[c_shed_poscor_medvar$R0 >= 50] <-50

R0_poscor<-ggplot(c_shed_poscor_lowvar, aes(x=R0))+
  geom_histogram(color="pink", fill="pink", alpha= .4, binwidth = 2)+
  geom_histogram(data=c_shed_poscor_medvar, aes(x=R0),binwidth = 2, alpha = .4, color="darkgreen", fill="darkgreen")+
  geom_histogram(data=c_shed_poscor_hivar, aes(x=R0), binwidth = 2, alpha = .4, color="darkblue", fill="darkblue")+
  geom_vline(xintercept = R0, linetype="dotdash", color = "red", size=.5)+
  geom_vline(xintercept = c_shed_poscor_hivar_mean_R0, linetype="dotdash", color = "darkblue", size=.5)+
  geom_vline(xintercept = c_shed_poscor_lowvar_mean_R0, linetype="dotdash", color = "pink", size=.5)+
  geom_vline(xintercept = c_shed_poscor_medvar_mean_R0, linetype="dotdash", color = "darkgreen", size=.5)+
  labs(title="c-shed-poscor")+
  scale_x_continuous(breaks=c(seq(0, 50, by=5)), labels=c(seq(0,45, by=5),"50+"))+
  ylim(0, 800)
R0_poscor

plot_grid(c_poscor, shed_poscor, R0_poscor, ncol=3, nrow=1)

c_negcor<- ggplot(c_shed_negcor_lowvar, aes(x=c))+
  geom_histogram(color="pink", fill="pink", binwidth=.1, alpha=.4)+
  geom_histogram(data=c_shed_negcor_medvar, aes(x=c), binwidth =.1, alpha = .4, color="darkgreen",fill="darkgreen")+
  geom_histogram(data=c_shed_negcor_hivar, aes(x=c),  alpha = .2, binwidth = .1, color="darkblue",fill="darkblue")+
  labs(title="contact-shed-negcor")


shed_negcor<- ggplot(c_shed_negcor_lowvar, aes(x=shed))+
  geom_histogram(color="pink", fill="pink", binwidth=.1, alpha=.4)+
  geom_histogram(data=c_shed_negcor_medvar, aes(x=shed), binwidth =.1, alpha = .4, color="darkgreen",fill="darkgreen")+
  geom_histogram(data=c_shed_negcor_hivar, aes(x=shed),  alpha = .2, binwidth = .1, color="darkblue",fill="darkblue")+
  labs(title="contact-shed-negcor")
shed_negcor

c_shed_negcor_hivar$R0[c_shed_negcor_hivar$R0 >= 50] <-50
c_shed_negcor_medvar$R0[c_shed_negcor_medvar$R0 >= 50] <-50


R0_negcor<-ggplot(c_shed_negcor_lowvar, aes(x=R0))+
  geom_histogram(color="pink", fill="pink", alpha= .4, binwidth = 2)+
  geom_histogram(data=c_shed_negcor_medvar, aes(x=R0), binwidth = 2, alpha = .4, color="darkgreen", fill="darkgreen")+
  geom_histogram(data=c_shed_negcor_hivar, aes(x=R0), binwidth = 2, alpha = .4, color="darkblue", fill="darkblue")+
  geom_vline(xintercept = R0, linetype="dotdash", color = "red", size=.5)+
  geom_vline(xintercept = c_shed_negcor_hivar_mean_R0, linetype="dotdash", color = "darkblue", size=.5)+
  geom_vline(xintercept = c_shed_negcor_lowvar_mean_R0, linetype="dotdash", color = "pink", size=.5)+
  geom_vline(xintercept = c_shed_negcor_medvar_mean_R0, linetype="dotdash", color = "darkgreen", size=.5)+
  labs(title="c-shed-negcor")+
  scale_x_continuous(breaks=c(seq(0, 50, by=5)), labels=c(seq(0,45, by=5),"50+"))+
  ylim(0, 800)
R0_negcor

plot_grid(c_nocor, shed_nocor, R0_nocor, c_poscor, shed_poscor, R0_poscor, 
          c_negcor, shed_negcor, R0_negcor, ncol=3, nrow=3)

plot_grid(R0_nocor,R0_poscor,R0_negcor, ncol=3, nrow=1)


