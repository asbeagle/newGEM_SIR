## alpha gamma
library(tidyverse)
library(cowplot)
library(ggplot2)

nocorr <- matrix(c(1,0,0,1), nrow=2, byrow=T)
negcorr <- matrix(c(1,-.5,-.5,1), nrow=2, byrow=T)
poscorr <- matrix(c(1,.5,.5,1), nrow=2, byrow=T)

## med R0 parms
hivar = c(   c=0.1,    shed=0.05,    alpha=0.1,    gamma=0.1, 
             sd_c=0.5, sd_shed=0.25, sd_alpha=0.5, sd_gamma=0.5, 
             b=2.5, d=.1, bs=.01)
medvar = c(   c=0.1,    shed=0.05,    alpha=0.1,    gamma=0.1, 
              sd_c=0.1, sd_shed=0.05, sd_alpha=0.1, sd_gamma=0.1, 
              b=2.5, d=.1, bs=.01)
lowvar = c(   c=0.1,     shed=0.05,    alpha=0.1,     gamma=0.1, 
              sd_c=0.02, sd_shed=0.01, sd_alpha=0.02, sd_gamma=0.02, 
              b=2.5, d=.1, bs=.01)

# high R0 parms
hivar = c(   c=0.15,    shed=0.1,    alpha=0.15,    gamma=0.15, 
             sd_c=0.75, sd_shed=0.5, sd_alpha=0.75, sd_gamma=0.75, 
             b=2.5, d=.1, bs=.01)
medvar = c(   c=0.15,    shed=0.1,    alpha=0.15,    gamma=0.15, 
              sd_c=0.15, sd_shed=0.1, sd_alpha=0.15, sd_gamma=0.15, 
              b=2.5, d=.1, bs=.01)
lowvar = c(   c=0.15,     shed=0.1,    alpha=0.15,     gamma=0.15, 
              sd_c=0.03, sd_shed=0.02, sd_alpha=0.03, sd_gamma=0.03, 
              b=2.5, d=.1, bs=.01)
R0 = 8


# low R0 parms
hivar = c(   c=0.1,    shed=0.02,    alpha=0.1,    gamma=0.1, 
             sd_c=0.5, sd_shed=0.1, sd_alpha=0.5, sd_gamma=0.5, 
             b=2.5, d=.1, bs=.01)
medvar = c(   c=0.1,    shed=0.02,    alpha=0.1,    gamma=0.1, 
              sd_c=0.1, sd_shed=0.02, sd_alpha=0.1, sd_gamma=0.1, 
              b=2.5, d=.1, bs=.01)
lowvar = c(   c=0.1,     shed=0.02,    alpha=0.1,     gamma=0.1, 
              sd_c=0.02, sd_shed=0.004, sd_alpha=0.02, sd_gamma=0.02, 
              b=2.5, d=.1, bs=.01)

R0 = 1.5

## high var
#no cor
data.frame(pick_individuals_multivariate(1000, traitmeans = c(hivar["alpha"],hivar["gamma"]), 
                                         traitsds = c(hivar["sd_alpha"], hivar["sd_gamma"]),corr=nocorr))%>% 
  mutate(., R0 = (235*(((hivar["shed"]/(1+hivar["shed"]))*(hivar["c"]))/(gamma + alpha + hivar["d"])))) -> alpha_gamma_nocor_hivar

alpha_gamma_nocor_hivar_mean_R0<-mean(alpha_gamma_nocor_hivar[,3])

#pos cor
data.frame(pick_individuals_multivariate(1000, traitmeans = c(hivar["alpha"],hivar["gamma"]), 
                                         traitsds = c(hivar["sd_alpha"], hivar["sd_gamma"]),corr=poscorr))%>% 
  mutate(., R0 = (235*(((hivar["shed"]/(1+hivar["shed"]))*(hivar["c"]))/(gamma + alpha + hivar["d"])))) -> alpha_gamma_poscor_hivar

alpha_gamma_poscor_hivar_mean_R0<-mean(alpha_gamma_poscor_hivar[,3])

#neg cor
data.frame(pick_individuals_multivariate(1000, traitmeans = c(hivar["alpha"],hivar["gamma"]), 
                                         traitsds = c(hivar["sd_alpha"], hivar["sd_gamma"]),corr=negcorr))%>% 
  mutate(., R0 = (235*(((hivar["shed"]/(1+hivar["shed"]))*(hivar["c"]))/(gamma + alpha + hivar["d"])))) -> alpha_gamma_negcor_hivar

alpha_gamma_negcor_hivar_mean_R0<-mean(alpha_gamma_negcor_hivar[,3])

## med var
#no cor
data.frame(pick_individuals_multivariate(1000, traitmeans = c(medvar["alpha"],medvar["gamma"]), 
                                         traitsds = c(medvar["sd_alpha"], medvar["sd_gamma"]),corr=nocorr))%>% 
  mutate(., R0 = (235*(((hivar["shed"]/(1+hivar["shed"]))*(hivar["c"]))/(gamma + alpha + hivar["d"])))) -> alpha_gamma_nocor_medvar

alpha_gamma_nocor_medvar_mean_R0<-mean(alpha_gamma_nocor_medvar[,3])

#pos cor
data.frame(pick_individuals_multivariate(1000, traitmeans = c(medvar["alpha"],medvar["gamma"]), 
                                         traitsds = c(medvar["sd_alpha"], medvar["sd_gamma"]),corr=poscorr))%>% 
  mutate(., R0 = (235*(((hivar["shed"]/(1+hivar["shed"]))*(hivar["c"]))/(gamma + alpha + hivar["d"])))) -> alpha_gamma_poscor_medvar

alpha_gamma_poscor_medvar_mean_R0<-mean(alpha_gamma_poscor_medvar[,3])

#neg cor
data.frame(pick_individuals_multivariate(1000, traitmeans = c(medvar["alpha"],medvar["gamma"]), 
                                         traitsds = c(medvar["sd_alpha"], medvar["sd_gamma"]),corr=negcorr))%>% 
  mutate(., R0 = (235*(((hivar["shed"]/(1+hivar["shed"]))*(hivar["c"]))/(gamma + alpha + hivar["d"])))) -> alpha_gamma_negcor_medvar

alpha_gamma_negcor_medvar_mean_R0<-mean(alpha_gamma_negcor_medvar[,3])

## low var
#no cor
data.frame(pick_individuals_multivariate(1000, traitmeans = c(lowvar["alpha"],lowvar["gamma"]), 
                                         traitsds = c(lowvar["sd_alpha"], lowvar["sd_gamma"]),corr=nocorr))%>% 
  mutate(., R0 = (235*(((hivar["shed"]/(1+hivar["shed"]))*(hivar["c"]))/(gamma + alpha + hivar["d"])))) -> alpha_gamma_nocor_lowvar

alpha_gamma_nocor_lowvar_mean_R0<-mean(alpha_gamma_nocor_lowvar[,3])

#pos cor
data.frame(pick_individuals_multivariate(1000, traitmeans = c(lowvar["alpha"],lowvar["gamma"]), 
                                         traitsds = c(lowvar["sd_alpha"], lowvar["sd_gamma"]),corr=poscorr))%>% 
  mutate(., R0 = (235*(((hivar["shed"]/(1+hivar["shed"]))*(hivar["c"]))/(gamma + alpha + hivar["d"])))) -> alpha_gamma_poscor_lowvar

alpha_gamma_poscor_lowvar_mean_R0<-mean(alpha_gamma_poscor_lowvar[,3])

#neg cor
data.frame(pick_individuals_multivariate(1000, traitmeans = c(lowvar["alpha"],lowvar["gamma"]), 
                                         traitsds = c(lowvar["sd_alpha"], lowvar["sd_gamma"]),corr=negcorr))%>% 
  mutate(., R0 = (235*(((hivar["shed"]/(1+hivar["shed"]))*(hivar["c"]))/(gamma + alpha + hivar["d"]))))-> alpha_gamma_negcor_lowvar

alpha_gamma_negcor_lowvar_mean_R0<-mean(alpha_gamma_negcor_lowvar[,3])

### plots
colors<- c("Hi Var" = "darkblue", "Med Var"="darkgreen", "Low Var" = "pink")
## bin extreme values
# no cov
alpha_gamma_nocor_hivar$R0[alpha_gamma_nocor_hivar$R0 >= 50] <-50
alpha_gamma_nocor_medvar$R0[alpha_gamma_nocor_medvar$R0 >= 50] <-50

# pos cov
alpha_gamma_poscor_hivar$R0[alpha_gamma_poscor_hivar$R0 >= 50] <- 50
alpha_gamma_poscor_medvar$R0[alpha_gamma_poscor_medvar$R0 >= 50] <-50

# neg cov
alpha_gamma_negcor_hivar$R0[alpha_gamma_negcor_hivar$R0 >= 50] <- 50
alpha_gamma_negcor_medvar$R0[alpha_gamma_negcor_medvar$R0 >= 50] <-50

## R0 distrubution plots

# no cov
R0_nocor<-ggplot(alpha_gamma_nocor_lowvar, aes(x=R0))+
  geom_histogram(aes(color="Low Var"), fill = "pink", alpha= .4, binwidth = .5)+
  geom_histogram(data=alpha_gamma_nocor_medvar, aes(x=R0,color="Med Var"), binwidth = .5, alpha = .4, fill="darkgreen")+
  geom_histogram(data=alpha_gamma_nocor_hivar, aes(x=R0, color="Hi Var"), binwidth = .5, alpha = .4, fill="darkblue")+
  geom_vline(xintercept = R0, linetype="dotdash", color = "red", size=.5)+
  geom_vline(xintercept = alpha_gamma_nocor_hivar_mean_R0, linetype="dotdash", color = "darkblue", size=.5)+
  geom_vline(xintercept = alpha_gamma_nocor_lowvar_mean_R0, linetype="dotdash", color = "pink", size=.5)+
  geom_vline(xintercept = alpha_gamma_nocor_medvar_mean_R0, linetype="dotdash", color = "darkgreen", size=.5)+
  labs(title="alpha-gamma-nocor", 
       color = "Legend")+
  scale_x_continuous(breaks=c(seq(0, 50, by=5)), labels=c(seq(0,45, by=5),"50+"))+
  ylim(0, 800)+
  scale_color_manual(values = c("darkblue", "pink", "darkgreen"))+
  guides(color = guide_legend(override.aes = list(fill = c("darkblue", "pink", "darkgreen"))))
R0_nocor

# pos cov
R0_poscor<-ggplot(alpha_gamma_poscor_lowvar, aes(x=R0))+
  geom_histogram(color="pink", fill="pink", alpha= .4, binwidth = .5)+
  geom_histogram(data=alpha_gamma_poscor_medvar, aes(x=R0), binwidth = .5, alpha = .4, color="darkgreen", fill="darkgreen")+
  geom_histogram(data=alpha_gamma_poscor_hivar, aes(x=R0), binwidth = .5, alpha = .4, color="darkblue", fill="darkblue")+
  geom_vline(xintercept = R0, linetype="dotdash", color = "red", size=.5)+
  geom_vline(xintercept = alpha_gamma_poscor_hivar_mean_R0, linetype="dotdash", color = "darkblue", size=.5)+
  geom_vline(xintercept = alpha_gamma_poscor_lowvar_mean_R0, linetype="dotdash", color = "pink", size=.5)+
  geom_vline(xintercept = alpha_gamma_poscor_medvar_mean_R0, linetype="dotdash", color = "darkgreen", size=.5)+
  labs(title="alpha-gamma-poscor")+
  scale_x_continuous(breaks=c(seq(0, 50, by=5)), labels=c(seq(0,45, by=5),"50+"))+
  ylim(0, 800)

# neg cov
R0_negcor<-ggplot(alpha_gamma_negcor_lowvar, aes(x=R0))+
  geom_histogram(color="pink", fill="pink", alpha= .4, binwidth = .5)+
  geom_histogram(data=alpha_gamma_negcor_medvar, aes(x=R0), binwidth = .5, alpha = .4, color="darkgreen", fill="darkgreen")+
  geom_histogram(data=alpha_gamma_negcor_hivar, aes(x=R0), binwidth = .5, alpha = .4, color="darkblue", fill="darkblue")+
  geom_vline(xintercept = R0, linetype="dotdash", color = "red", size=.5)+
  geom_vline(xintercept = alpha_gamma_negcor_hivar_mean_R0, linetype="dotdash", color = "darkblue", size=.5)+
  geom_vline(xintercept = alpha_gamma_negcor_lowvar_mean_R0, linetype="dotdash", color = "pink", size=.5)+
  geom_vline(xintercept = alpha_gamma_negcor_medvar_mean_R0, linetype="dotdash", color = "darkgreen", size=.5)+
  labs(title="alpha-gamma-negcor")+
  scale_x_continuous(breaks=c(seq(0, 50, by=5)), labels=c(seq(0,45, by=5),"50+"))+
  ylim(0, 800)
R0_negcor

ggarrange(R0_nocor, R0_poscor, R0_negcor, ncol=3, nrow=1, common.legend = TRUE, legend="right")


## histogram of each covariation with increasing variation
alpha_nocor<- ggplot(alpha_gamma_nocor_lowvar, aes(x=alpha))+
  geom_histogram(color="pink", fill="pink", binwidth=.1, alpha=.4)+
  geom_histogram(data=alpha_gamma_nocor_medvar, aes(x=alpha), binwidth =.1, alpha = .4, color="darkgreen",fill="darkgreen")+
  geom_histogram(data=alpha_gamma_nocor_hivar, aes(x=alpha),  alpha = .4, binwidth = .1, color="darkblue",fill="darkblue")+
  labs(title="alpha-gamma-nocor")

gamma_nocor<- ggplot(alpha_gamma_nocor_lowvar, aes(x=gamma))+
  geom_histogram(color="pink", fill="pink", binwidth=.1, alpha=.4)+
  geom_histogram(data=alpha_gamma_nocor_medvar, aes(x=gamma), binwidth =.1, alpha = .4, color="darkgreen",fill="darkgreen")+
  geom_histogram(data=alpha_gamma_nocor_hivar, aes(x=gamma),  alpha = .4, binwidth = .1, color="darkblue",fill="darkblue")+
  labs(title="alpha-gamma-nocor")

R0_nocor<-ggplot(alpha_gamma_nocor_lowvar, aes(x=R0))+
  geom_histogram(color="pink", fill="pink", alpha= .4, binwidth = 0.01)+
  geom_histogram(data=alpha_gamma_nocor_medvar, aes(x=R0), binwidth = 0.01, alpha = .4, color="darkgreen", fill="darkgreen")+
  geom_histogram(data=alpha_gamma_nocor_hivar, aes(x=R0), binwidth = 0.01, alpha = .4, color="darkblue", fill="darkblue")+
  geom_vline(xintercept = R0, linetype="dotdash", color = "red", size=.5)+
  geom_vline(xintercept = alpha_gamma_nocor_hivar_mean_R0, linetype="dotdash", color = "darkblue", size=.5)+
  geom_vline(xintercept = alpha_gamma_nocor_lowvar_mean_R0, linetype="dotdash", color = "pink", size=.5)+
  geom_vline(xintercept = alpha_gamma_nocor_medvar_mean_R0, linetype="dotdash", color = "darkgreen", size=.5)+
  labs(title="alpha-gamma-nocor")+
  ylim(0, 800)

R0_nocor

alpha_poscor<- ggplot(alpha_gamma_poscor_lowvar, aes(x=alpha))+
  geom_histogram(color="pink", fill="pink", binwidth=.1, alpha=.4)+
  geom_histogram(data=alpha_gamma_poscor_medvar, aes(x=alpha), binwidth =.1, alpha = .4, color="darkgreen",fill="darkgreen")+
  geom_histogram(data=alpha_gamma_poscor_hivar, aes(x=alpha),  alpha = .4, binwidth = .1, color="darkblue",fill="darkblue")+
  labs(title="alpha-gamma-poscor")

gamma_poscor<- ggplot(alpha_gamma_poscor_lowvar, aes(x=gamma))+
  geom_histogram(color="pink", fill="pink", binwidth=.1, alpha=.4)+
  geom_histogram(data=alpha_gamma_poscor_medvar, aes(x=gamma), binwidth =.1, alpha = .4, color="darkgreen",fill="darkgreen")+
  geom_histogram(data=alpha_gamma_poscor_hivar, aes(x=gamma),  alpha = .4, binwidth = .1, color="darkblue",fill="darkblue")+
  labs(title="alpha-gamma-poscor")

R0_poscor<-ggplot(alpha_gamma_poscor_lowvar, aes(x=R0))+
  geom_histogram(color="pink", fill="pink", alpha= .4, binwidth = 0.01)+
  geom_histogram(data=alpha_gamma_poscor_medvar, aes(x=R0), binwidth = 0.01, alpha = .4, color="darkgreen", fill="darkgreen")+
  geom_histogram(data=alpha_gamma_poscor_hivar, aes(x=R0), binwidth = 0.01, alpha = .4, color="darkblue", fill="darkblue")+
  geom_vline(xintercept = R0, linetype="dotdash", color = "red", size=.5)+
  geom_vline(xintercept = alpha_gamma_poscor_hivar_mean_R0, linetype="dotdash", color = "darkblue", size=.5)+
  geom_vline(xintercept = alpha_gamma_poscor_lowvar_mean_R0, linetype="dotdash", color = "pink", size=.5)+
  geom_vline(xintercept = alpha_gamma_poscor_medvar_mean_R0, linetype="dotdash", color = "darkgreen", size=.5)+
  labs(title="alpha-gamma-negcor")+
  ylim(0, 800)
R0_poscor

alpha_negcor<- ggplot(alpha_gamma_negcor_lowvar, aes(x=alpha))+
  geom_histogram(color="pink", fill="pink", binwidth=.1, alpha=.4)+
  geom_histogram(data=alpha_gamma_negcor_medvar, aes(x=alpha), binwidth =.1, alpha = .4, color="darkgreen",fill="darkgreen")+
  geom_histogram(data=alpha_gamma_negcor_hivar, aes(x=alpha),  alpha = .2, binwidth = .1, color="darkblue",fill="darkblue")+
  labs(title="alpha-gamma-negcor")

gamma_negcor<- ggplot(alpha_gamma_negcor_lowvar, aes(x=gamma))+
  geom_histogram(color="pink", fill="pink", binwidth=.1, alpha=.4)+
  geom_histogram(data=alpha_gamma_negcor_medvar, aes(x=gamma), binwidth =.1, alpha = .4, color="darkgreen",fill="darkgreen")+
  geom_histogram(data=alpha_gamma_negcor_hivar, aes(x=gamma),  alpha = .2, binwidth = .1, color="darkblue",fill="darkblue")+
  labs(title="alpha-gamma-negcor")

R0_negcor<-ggplot(alpha_gamma_negcor_lowvar, aes(x=R0))+
  geom_histogram(color="pink", fill="pink", alpha= .4, binwidth = 0.01)+
  geom_histogram(data=alpha_gamma_negcor_medvar, aes(x=R0), binwidth = 0.01, alpha = .4, color="darkgreen", fill="darkgreen")+
  geom_histogram(data=alpha_gamma_negcor_hivar, aes(x=R0), binwidth = 0.01, alpha = .4, color="darkblue", fill="darkblue")+
  geom_vline(xintercept = R0, linetype="dotdash", color = "red", size=.5)+
  geom_vline(xintercept = alpha_gamma_negcor_hivar_mean_R0, linetype="dotdash", color = "darkblue", size=.5)+
  geom_vline(xintercept = alpha_gamma_negcor_lowvar_mean_R0, linetype="dotdash", color = "pink", size=.5)+
  geom_vline(xintercept = alpha_gamma_negcor_medvar_mean_R0, linetype="dotdash", color = "darkgreen", size=.5)+
  labs(title="alpha-gamma-negcor")+
  ylim(0, 800)
R0_negcor

plot_grid(alpha_nocor, gamma_nocor, R0_nocor, alpha_poscor, gamma_poscor, R0_poscor, 
          alpha_negcor, gamma_negcor, R0_negcor, ncol=3, nrow=3)

plot_grid(R0_nocor, R0_poscor,R0_negcor, ncol=3, nrow=1)
