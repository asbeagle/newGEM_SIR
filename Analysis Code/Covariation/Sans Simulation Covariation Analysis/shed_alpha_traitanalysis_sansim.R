## shed alpha

library(tidyverse)
library(cowplot)
library(ggplot2)

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


## high var
#no cor
data.frame(pick_individuals_multivariate(1000, traitmeans = c(hivar["shed"],hivar["alpha"]), 
                                         traitsds = c(hivar["sd_shed"], hivar["sd_alpha"]),corr=nocorr))%>% 
  mutate(., R0 = (235*(((shed/(1+shed))*(hivar["c"]))/(hivar["gamma"] + alpha + hivar["d"])))) ->shed_alpha_nocor_hivar

shed_alpha_nocor_hivar_mean_R0<-mean(shed_alpha_nocor_hivar[,3])

#pos cor
data.frame(pick_individuals_multivariate(1000, traitmeans = c(hivar["shed"],hivar["alpha"]), 
                                         traitsds = c(hivar["sd_shed"], hivar["sd_alpha"]),corr=poscorr))%>% 
  mutate(., R0 = (235*(((shed/(1+shed))*(hivar["c"]))/(hivar["gamma"] + alpha + hivar["d"])))) ->shed_alpha_poscor_hivar

shed_alpha_poscor_hivar_mean_R0<-mean(shed_alpha_poscor_hivar[,3])

#neg cor
data.frame(pick_individuals_multivariate(1000, traitmeans = c(hivar["shed"],hivar["alpha"]), 
                                         traitsds = c(hivar["sd_shed"], hivar["sd_alpha"]),corr=negcorr))%>% 
  mutate(., R0 = (235*(((shed/(1+shed))*(hivar["c"]))/(hivar["gamma"] + alpha + hivar["d"])))) ->shed_alpha_negcor_hivar

shed_alpha_negcor_hivar_mean_R0<-mean(shed_alpha_negcor_hivar[,3])

## med var
#no cor
data.frame(pick_individuals_multivariate(1000, traitmeans = c(medvar["shed"],medvar["alpha"]), 
                                         traitsds = c(medvar["sd_shed"], medvar["sd_alpha"]),corr=nocorr))%>% 
  mutate(., R0 = (235*(((shed/(1+shed))*(hivar["c"]))/(hivar["gamma"] + alpha + hivar["d"])))) ->shed_alpha_nocor_medvar

shed_alpha_nocor_medvar_mean_R0<-mean(shed_alpha_nocor_medvar[,3])

#pos cor
data.frame(pick_individuals_multivariate(1000, traitmeans = c(medvar["shed"],medvar["alpha"]), 
                                         traitsds = c(medvar["sd_shed"], medvar["sd_alpha"]),corr=poscorr))%>% 
  mutate(., R0 = (235*(((shed/(1+shed))*(hivar["c"]))/(hivar["gamma"] + alpha + hivar["d"]))))->shed_alpha_poscor_medvar

shed_alpha_poscor_medvar_mean_R0<-mean(shed_alpha_poscor_medvar[,3])

#neg cor
data.frame(pick_individuals_multivariate(1000, traitmeans = c(medvar["shed"],medvar["alpha"]), 
                                         traitsds = c(medvar["sd_shed"], medvar["sd_alpha"]),corr=negcorr))%>% 
  mutate(., R0 = (235*(((shed/(1+shed))*(hivar["c"]))/(hivar["gamma"] + alpha + hivar["d"])))) ->shed_alpha_negcor_medvar

shed_alpha_negcor_medvar_mean_R0<-mean(shed_alpha_negcor_medvar[,3])

## low var
#no cor
data.frame(pick_individuals_multivariate(1000, traitmeans = c(lowvar["shed"],lowvar["alpha"]), 
                                         traitsds = c(lowvar["sd_shed"], lowvar["sd_alpha"]),corr=nocorr))%>% 
  mutate(., R0 = (235*(((shed/(1+shed))*(hivar["c"]))/(hivar["gamma"] + alpha + hivar["d"])))) ->shed_alpha_nocor_lowvar

shed_alpha_nocor_lowvar_mean_R0<-mean(shed_alpha_nocor_lowvar[,3])

#pos cor
data.frame(pick_individuals_multivariate(1000, traitmeans = c(lowvar["shed"],lowvar["alpha"]), 
                                         traitsds = c(lowvar["sd_shed"], lowvar["sd_alpha"]),corr=poscorr))%>% 
  mutate(., R0 = (235*(((shed/(1+shed))*(hivar["c"]))/(hivar["gamma"] + alpha + hivar["d"])))) ->shed_alpha_poscor_lowvar

shed_alpha_poscor_lowvar_mean_R0<-mean(shed_alpha_poscor_lowvar[,3])

#neg cor
data.frame(pick_individuals_multivariate(1000, traitmeans = c(lowvar["shed"],lowvar["alpha"]), 
                                         traitsds = c(lowvar["sd_shed"], lowvar["sd_alpha"]),corr=negcorr))%>% 
  mutate(., R0 = (235*(((shed/(1+shed))*(hivar["c"]))/(hivar["gamma"] + alpha + hivar["d"])))) ->shed_alpha_negcor_lowvar

shed_alpha_negcor_lowvar_mean_R0<-mean(shed_alpha_negcor_lowvar[,3])

### plots
colors<- c("Hi Var" = "darkblue", "Med Var"="darkgreen", "Low Var" = "pink")
## bin extreme values
# no cov
shed_alpha_nocor_hivar$R0[shed_alpha_nocor_hivar$R0 >= 50] <-50
shed_alpha_nocor_medvar$R0[shed_alpha_nocor_medvar$R0 >= 50] <-50

# pos cov
shed_alpha_poscor_hivar$R0[shed_alpha_poscor_hivar$R0 >= 50] <- 50
shed_alpha_poscor_medvar$R0[shed_alpha_poscor_medvar$R0 >= 50] <-50

# neg cov
shed_alpha_negcor_hivar$R0[shed_alpha_negcor_hivar$R0 >= 50] <- 50
shed_alpha_negcor_medvar$R0[shed_alpha_negcor_medvar$R0 >= 50] <-50

## R0 distrubution plots
# no cov
R0_nocor<-ggplot(shed_alpha_nocor_lowvar, aes(x=R0))+
  geom_histogram(aes(color="Low Var"), fill = "pink", alpha= .4, binwidth = 2)+
  geom_histogram(data=shed_alpha_nocor_medvar, aes(x=R0,color="Med Var"), binwidth = 2, alpha = .4, fill="darkgreen")+
  geom_histogram(data=shed_alpha_nocor_hivar, aes(x=R0, color="Hi Var"), binwidth = 2, alpha = .4, fill="darkblue")+
  geom_vline(xintercept = R0, linetype="dotdash", color = "red", size=.5)+
  geom_vline(xintercept = shed_alpha_nocor_hivar_mean_R0, linetype="dotdash", color = "darkblue", size=.5)+
  geom_vline(xintercept = shed_alpha_nocor_lowvar_mean_R0, linetype="dotdash", color = "pink", size=.5)+
  geom_vline(xintercept = shed_alpha_nocor_medvar_mean_R0, linetype="dotdash", color = "darkgreen", size=.5)+
  labs(title="shed-alpha-nocor", 
       color = "Legend")+
  scale_x_continuous(breaks=c(seq(0, 50, by=5)), labels=c(seq(0,45, by=5),"50+"))+
  ylim(0, 800)+
  scale_color_manual(values = c("darkblue", "pink", "darkgreen"))+
  guides(color = guide_legend(override.aes = list(fill = c("darkblue", "pink", "darkgreen"))))

# pos cov
R0_poscor<-ggplot(shed_alpha_poscor_lowvar, aes(x=R0))+
  geom_histogram(color="pink", fill="pink", alpha= .4, binwidth = 2)+
  geom_histogram(data=shed_alpha_poscor_medvar, aes(x=R0), binwidth = 2, alpha = .4, color="darkgreen", fill="darkgreen")+
  geom_histogram(data=shed_alpha_poscor_hivar, aes(x=R0), binwidth = 2, alpha = .4, color="darkblue", fill="darkblue")+
  geom_vline(xintercept = R0, linetype="dotdash", color = "red", size=.5)+
  geom_vline(xintercept = shed_alpha_poscor_hivar_mean_R0, linetype="dotdash", color = "darkblue", size=.5)+
  geom_vline(xintercept = shed_alpha_poscor_lowvar_mean_R0, linetype="dotdash", color = "pink", size=.5)+
  geom_vline(xintercept = shed_alpha_poscor_medvar_mean_R0, linetype="dotdash", color = "darkgreen", size=.5)+
  labs(title="shed-alpha-poscor")+
  scale_x_continuous(breaks=c(seq(0, 50, by=5)), labels=c(seq(0,45, by=5),"50+"))+
  ylim(0, 800)

#neg cov
R0_negcor<-ggplot(shed_alpha_negcor_lowvar, aes(x=R0))+
  geom_histogram(color="pink", fill="pink", alpha= .4, binwidth = 2)+
  geom_histogram(data=shed_alpha_negcor_medvar, aes(x=R0), binwidth = 2, alpha = .4, color="darkgreen", fill="darkgreen")+
  geom_histogram(data=shed_alpha_negcor_hivar, aes(x=R0), binwidth = 2, alpha = .4, color="darkblue", fill="darkblue")+
  geom_vline(xintercept = R0, linetype="dotdash", color = "red", size=.5)+
  geom_vline(xintercept = shed_alpha_negcor_hivar_mean_R0, linetype="dotdash", color = "darkblue", size=.5)+
  geom_vline(xintercept = shed_alpha_negcor_lowvar_mean_R0, linetype="dotdash", color = "pink", size=.5)+
  geom_vline(xintercept = shed_alpha_negcor_medvar_mean_R0, linetype="dotdash", color = "darkgreen", size=.5)+
  labs(title="shed-alpha-negcor")+
  scale_x_continuous(breaks=c(seq(0, 50, by=5)), labels=c(seq(0,45, by=5),"50+"))+
  ylim(0, 800)

ggarrange(R0_nocor, R0_poscor, R0_negcor, ncol=3, nrow=1, common.legend = TRUE, legend="right")

## histogram of each covariation with increasing variation
shed_nocor<- ggplot(shed_alpha_nocor_lowvar, aes(x=shed))+
  geom_histogram(color="pink", fill="pink", binwidth=.1, alpha=.4)+
  geom_histogram(data=shed_alpha_nocor_medvar, aes(x=shed), binwidth =.1, alpha = .4, color="darkgreen",fill="darkgreen")+
  geom_histogram(data=shed_alpha_nocor_hivar, aes(x=shed),  alpha = .2, binwidth = .1, color="darkblue",fill="darkblue")+
  labs(title="shed-alpha-nocor")


alpha_nocor<- ggplot(c_alpha_nocor_lowvar, aes(x=alpha))+
  geom_histogram(color="pink", fill="pink", binwidth=.1, alpha=.4)+
  geom_histogram(data=c_alpha_nocor_medvar, aes(x=alpha), binwidth =.1, alpha = .4, color="darkgreen",fill="darkgreen")+
  geom_histogram(data=c_alpha_nocor_hivar, aes(x=alpha),  alpha = .4, binwidth = .1, color="darkblue",fill="darkblue")+
  labs(title="shed-alpha-nocor")


R0_nocor<-ggplot(shed_alpha_nocor_lowvar, aes(x=R0))+
  geom_histogram(color="pink", fill="pink", alpha= .4, binwidth = 0.01)+
  geom_histogram(data=shed_alpha_nocor_medvar, aes(x=R0), binwidth = 0.01, alpha = .4, color="darkgreen", fill="darkgreen")+
  geom_histogram(data=shed_alpha_nocor_hivar, aes(x=R0), binwidth = 0.01, alpha = .4, color="darkblue", fill="darkblue")+
  geom_vline(xintercept = R0, linetype="solid", color = "green", size=.4)

shed_poscor<- ggplot(shed_alpha_poscor_lowvar, aes(x=shed))+
  geom_histogram(color="pink", fill="pink", binwidth=.1, alpha=.4)+
  geom_histogram(data=shed_alpha_poscor_medvar, aes(x=shed), binwidth =.1, alpha = .4, color="darkgreen",fill="darkgreen")+
  geom_histogram(data=shed_alpha_poscor_hivar, aes(x=shed),  alpha = .2, binwidth = .1, color="darkblue",fill="darkblue")+
  labs(title="shed-alpha-poscor")


alpha_poscor<- ggplot(shed_alpha_poscor_lowvar, aes(x=alpha))+
  geom_histogram(color="pink", fill="pink", binwidth=.1, alpha= .4)+
  geom_histogram(data=shed_alpha_poscor_medvar, aes(x=alpha), binwidth =.1, alpha = .4, color="darkgreen",fill="darkgreen")+
  geom_histogram(data=shed_alpha_poscor_hivar, aes(x=alpha),  alpha = .4, binwidth = .1, color="darkblue",fill="darkblue")+
  labs(title="shed-alpha-poscor")


R0_poscor<-ggplot(shed_alpha_poscor_lowvar, aes(x=R0))+
  geom_histogram(color="pink", fill="pink", alpha= .4, binwidth = 0.01)+
  geom_histogram(data=shed_alpha_poscor_medvar, aes(x=R0), binwidth = 0.01, alpha = .4, color="darkgreen", fill="darkgreen")+
  geom_histogram(data=shed_alpha_poscor_hivar, aes(x=R0), binwidth = 0.01, alpha = .4, color="darkblue", fill="darkblue")+
  geom_vline(xintercept = R0, linetype="solid", color = "green", size=.4)

shed_negcor<- ggplot(shed_alpha_negcor_lowvar, aes(x=shed))+
  geom_histogram(color="pink", fill="pink", binwidth=.1, alpha=.4)+
  geom_histogram(data=shed_alpha_negcor_medvar, aes(x=shed), binwidth =.1, alpha = .4, color="darkgreen",fill="darkgreen")+
  geom_histogram(data=shed_alpha_negcor_hivar, aes(x=shed),  alpha = .2, binwidth = .1, color="darkblue",fill="darkblue")+
  labs(title="shed-alpha-negcor")

alpha_negcor<- ggplot(shed_alpha_negcor_lowvar, aes(x=alpha))+
  geom_histogram(color="pink", fill="pink", binwidth=.1, alpha=.4)+
  geom_histogram(data=shed_alpha_negcor_medvar, aes(x=alpha), binwidth =.1, alpha = .4, color="darkgreen",fill="darkgreen")+
  geom_histogram(data=shed_alpha_negcor_hivar, aes(x=alpha),  alpha = .2, binwidth = .1, color="darkblue",fill="darkblue")+
  labs(title="shed-alpha-negcor")

R0_negcor<-ggplot(shed_alpha_negcor_lowvar, aes(x=R0))+
  geom_histogram(color="pink", fill="pink", alpha= .4, binwidth = 0.01)+
  geom_histogram(data=shed_alpha_negcor_medvar, aes(x=R0), binwidth = 0.01, alpha = .4, color="darkgreen", fill="darkgreen")+
  geom_histogram(data=shed_alpha_negcor_hivar, aes(x=R0), binwidth = 0.01, alpha = .4, color="darkblue", fill="darkblue")+
  geom_vline(xintercept = R0, linetype="solid", color = "green", size=.4)

plot_grid(shed_nocor, alpha_nocor, R0_nocor, shed_poscor, alpha_poscor, R0_poscor, 
          shed_negcor, alpha_negcor, R0_negcor, ncol=3, nrow=3)

