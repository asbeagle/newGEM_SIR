## Load dataframes with the results from 1,000 stochastic replicate simulations
## Dataframes contain the following variables:
## peak = max. number of infected individuals
## fadeout = 0/1 did fadeout occur?
## fadeoutT = when did fadeout occur (if it did)?
## disp = dispersion parameter of a negative binomial distribution fit to the observed number of secondary infections caused by each infected individual
## R0 = deterministic R0 of the parameter set

## All individuals have identical trait values
noVar = readRDS(file="Stochastic_no_variation_simulations.RDS")

## These dataframes have additional information
## totalI = total number of infected individuals (not just peak)
## CV = coefficient of variation in trait that varies among individuals
## Individuals vary in contact rate
## Remove unnecessary text to turn R0 and CV into numerical rather than strings
cVar = readRDS(file="Variance_in_contact_results_12-28-23.RDS")
#cVar$R0 = (strsplit(cVar$R0,"=") %>% lapply(., function(a) as.numeric(a[2])) %>% unlist)
#cVar$CV = (strsplit(cVar$CV,"=") %>% lapply(., function(a) as.numeric(a[2])) %>% unlist)
## Individuals vary in shedding/infectiousness
sVar = readRDS(file="Variance_in_shedding_results_12-28-23.RDS")
#sVar$R0 = (strsplit(sVar$R0,"=") %>% lapply(., function(a) as.numeric(a[2])) %>% unlist)
#sVar$CV = (strsplit(sVar$CV,"=") %>% lapply(., function(a) as.numeric(a[2])) %>% unlist)
## Individuals vary in virulence rate
aVar = readRDS(file="Variance_in_virulence_results_12-28-23.RDS")
#aVar$R0 = (strsplit(aVar$R0,"=") %>% lapply(., function(a) as.numeric(a[2])) %>% unlist)
#aVar$CV = (strsplit(aVar$CV,"=") %>% lapply(., function(a) as.numeric(a[2])) %>% unlist)
## Individuals vary in recovery rate
gVar = readRDS(file="Variance_in_recovery_results_12-28-23.RDS")
#gVar$R0 = (strsplit(gVar$R0,"=") %>% lapply(., function(a) as.numeric(a[2])) %>% unlist)
#gVar$CV = (strsplit(gVar$CV,"=") %>% lapply(., function(a) as.numeric(a[2])) %>% unlist)

## These dataframes have even more information
## These simulations vary two traits simultaneously, with a correlation that varies between -0.9 and 0.9.
## corr = correlation value
## cov = character vector telling which traits are covarying
# R = 1, 4, 8; corr = -0.9 to 0.9; CV = 0.2, 1, 5, vars = c, s, a, g
csCov = caCov = cgCov = saCov = sgCov = agCov = vector('list', length=18*9)
i = 1;
for (corr in c(seq(-0.9, -0.1, 0.1), seq(0.1, 0.9, 0.1))) {
  for (R in c(1,4,8)) {
    for (CV in c(0.2,1,5)) {
      csCov[[i]] = readRDS(file=paste0("out_R0=",R,"_CV=",CV,"_corr=",corr,"_cov=c-s.RDS"))
      caCov[[i]] = readRDS(file=paste0("out_R0=",R,"_CV=",CV,"_corr=",corr,"_cov=c-a.RDS"))
      cgCov[[i]] = readRDS(file=paste0("out_R0=",R,"_CV=",CV,"_corr=",corr,"_cov=c-g.RDS"))
      saCov[[i]] = readRDS(file=paste0("out_R0=",R,"_CV=",CV,"_corr=",corr,"_cov=s-a.RDS"))
      sgCov[[i]] = readRDS(file=paste0("out_R0=",R,"_CV=",CV,"_corr=",corr,"_cov=s-g.RDS"))
      agCov[[i]] = readRDS(file=paste0("out_R0=",R,"_CV=",CV,"_corr=",corr,"_cov=a-g.RDS"))
      i = i+1
    }
  }
}
csCov = do.call("rbind.data.frame",csCov)
caCov = do.call("rbind.data.frame",caCov)
cgCov = do.call("rbind.data.frame",cgCov)
saCov = do.call("rbind.data.frame",saCov)
sgCov = do.call("rbind.data.frame",sgCov)
agCov = do.call("rbind.data.frame",agCov)

## Round correlation so that everything plots correctly
csCov$corr = round(csCov$corr,1)
caCov$corr = round(csCov$corr,1)
cgCov$corr = round(csCov$corr,1)
saCov$corr = round(csCov$corr,1)
sgCov$corr = round(csCov$corr,1)
agCov$corr = round(csCov$corr,1)

## Turn R0 and CV into strings for nice plotting
csCov$R0 = paste0("R0=",csCov$R0)
csCov$CV = paste0("CV=",csCov$CV)
caCov$R0 = paste0("R0=",caCov$R0)
caCov$CV = paste0("CV=",caCov$CV)
cgCov$R0 = paste0("R0=",cgCov$R0)
cgCov$CV = paste0("CV=",cgCov$CV)
saCov$R0 = paste0("R0=",saCov$R0)
saCov$CV = paste0("CV=",saCov$CV)
sgCov$R0 = paste0("R0=",sgCov$R0)
sgCov$CV = paste0("CV=",sgCov$CV)
agCov$R0 = paste0("R0=",agCov$R0)
agCov$CV = paste0("CV=",agCov$CV)


## Plots that I want to make
## 1. Probability of fadeout
## 2. Peak infections when R0 = 1
## 3. Peak infections for R0 = 4 and 8 (ignoring fadeouts)
## 4. Dispersion parameters (superspreading)

## 1.
png(file="Fadeout_probability_c-s.png", width=5, height=5, units='in', res=450)
csCov %>% 
  group_by(corr, R0, CV) %>% 
  summarize(fadeoutProb=sum(fadeout)/1000) %>%
  ggplot(., aes(x=corr, y=fadeoutProb, color="Trait covariation")) + 
  geom_hline(data=(cVar %>% group_by(R0,CV) %>% summarize(fadeoutProb=sum(fadeout)/1000)), mapping=aes(yintercept=fadeoutProb, color="Contact variation"), linewidth=1, alpha=0.5) +
  geom_hline(data=(sVar %>% group_by(R0,CV) %>% summarize(fadeoutProb=sum(fadeout)/1000)), mapping=aes(yintercept=fadeoutProb, color="Shedding variation"), linewidth=1, alpha=0.5) +
  geom_line(linewidth=1) + 
  scale_color_manual(values = c(
    "Contact variation" = "red",
    "Shedding variation" = "blue",
    "Trait covariation" = "black")) + 
  facet_grid(R0~CV) +  
  xlab("Trait correlation") + 
  ylab("Probability of fadeout") + 
  theme_bw() + 
  theme(legend.position="bottom", legend.title=element_blank())
dev.off()

png(file="Fadeout_probability_a-g.png", width=5, height=5, units='in', res=450)
agCov %>% 
  group_by(corr, R0, CV) %>% 
  summarize(fadeoutProb=sum(fadeout)/1000) %>%
  ggplot(., aes(x=corr, y=fadeoutProb, color="Trait covariation")) + 
  geom_hline(data=(aVar %>% group_by(R0,CV) %>% summarize(fadeoutProb=sum(fadeout)/1000)), mapping=aes(yintercept=fadeoutProb, color="Mortality variation"), linewidth=1, alpha=0.5) +
  geom_hline(data=(gVar %>% group_by(R0,CV) %>% summarize(fadeoutProb=sum(fadeout)/1000)), mapping=aes(yintercept=fadeoutProb, color="Recovery variation"), linewidth=1, alpha=0.5) +
  geom_line(linewidth=1) + 
  scale_color_manual(values = c(
    "Mortality variation" = "red",
    "Recovery variation" = "blue",
    "Trait covariation" = "black")) + facet_grid(R0~CV) +  
  xlab("Trait correlation") + 
  ylab("Probability of fadeout") + 
  theme_bw() + 
  theme(legend.position="bottom", legend.title=element_blank())
dev.off()

png(file="Fadeout_probability_s-a.png", width=5, height=5, units='in', res=450)
saCov %>% 
  group_by(corr, R0, CV) %>% 
  summarize(fadeoutProb=sum(fadeout)/1000) %>%
  ggplot(., aes(x=corr, y=fadeoutProb, color="Trait covariation")) + 
  geom_hline(data=(sVar %>% group_by(R0,CV) %>% summarize(fadeoutProb=sum(fadeout)/1000)), mapping=aes(yintercept=fadeoutProb, color="Shedding variation"), linewidth=1, alpha=0.5) +
  geom_hline(data=(aVar %>% group_by(R0,CV) %>% summarize(fadeoutProb=sum(fadeout)/1000)), mapping=aes(yintercept=fadeoutProb, color="Mortality variation"), linewidth=1, alpha=0.5) +
  geom_line(linewidth=1) + 
  scale_color_manual(values = c(
    "Shedding variation" = "blue",
    "Mortality variation" = "red",
    "Trait covariation" = "black")) + facet_grid(R0~CV) +  
  xlab("Trait correlation") + 
  ylab("Probability of fadeout") + 
  theme_bw() + 
  theme(legend.position="bottom", legend.title=element_blank())
dev.off()

png(file="Fadeout_probability_c-a.png", width=5, height=5, units='in', res=450)
caCov %>% 
  group_by(corr, R0, CV) %>% 
  summarize(fadeoutProb=sum(fadeout)/1000) %>%
  ggplot(., aes(x=corr, y=fadeoutProb, color="Trait covariation")) + 
  geom_hline(data=(cVar %>% group_by(R0,CV) %>% summarize(fadeoutProb=sum(fadeout)/1000)), mapping=aes(yintercept=fadeoutProb, color="Contact variation"), linewidth=1, alpha=0.5) +
  geom_hline(data=(aVar %>% group_by(R0,CV) %>% summarize(fadeoutProb=sum(fadeout)/1000)), mapping=aes(yintercept=fadeoutProb, color="Mortality variation"), linewidth=1, alpha=0.5) +
  geom_line(linewidth=1) + 
  scale_color_manual(values = c(
    "Contact variation" = "blue",
    "Mortality variation" = "red",
    "Trait covariation" = "black")) + facet_grid(R0~CV) +  
  xlab("Trait correlation") + 
  ylab("Probability of fadeout") + 
  theme_bw() + 
  theme(legend.position="bottom", legend.title=element_blank())
dev.off()

png(file="Fadeout_probability_s-g.png", width=5, height=5, units='in', res=450)
saCov %>% 
  group_by(corr, R0, CV) %>% 
  summarize(fadeoutProb=sum(fadeout)/1000) %>%
  ggplot(., aes(x=corr, y=fadeoutProb, color="Trait covariation")) + 
  geom_hline(data=(sVar %>% group_by(R0,CV) %>% summarize(fadeoutProb=sum(fadeout)/1000)), mapping=aes(yintercept=fadeoutProb, color="Shedding variation"), linewidth=1, alpha=0.5) +
  geom_hline(data=(gVar %>% group_by(R0,CV) %>% summarize(fadeoutProb=sum(fadeout)/1000)), mapping=aes(yintercept=fadeoutProb, color="Recovery variation"), linewidth=1, alpha=0.5) +
  geom_line(linewidth=1) + 
  scale_color_manual(values = c(
    "Shedding variation" = "blue",
    "Recovery variation" = "red",
    "Trait covariation" = "black")) + facet_grid(R0~CV) +  
  xlab("Trait correlation") + 
  ylab("Probability of fadeout") + 
  theme_bw() + 
  theme(legend.position="bottom", legend.title=element_blank())
dev.off()

png(file="Fadeout_probability_c-g.png", width=5, height=5, units='in', res=450)
caCov %>% 
  group_by(corr, R0, CV) %>% 
  summarize(fadeoutProb=sum(fadeout)/1000) %>%
  ggplot(., aes(x=corr, y=fadeoutProb, color="Trait covariation")) + 
  geom_hline(data=(cVar %>% group_by(R0,CV) %>% summarize(fadeoutProb=sum(fadeout)/1000)), mapping=aes(yintercept=fadeoutProb, color="Contact variation"), linewidth=1, alpha=0.5) +
  geom_hline(data=(gVar %>% group_by(R0,CV) %>% summarize(fadeoutProb=sum(fadeout)/1000)), mapping=aes(yintercept=fadeoutProb, color="Recovery variation"), linewidth=1, alpha=0.5) +
  geom_line(linewidth=1) + 
  scale_color_manual(values = c(
    "Contact variation" = "blue",
    "Recovery variation" = "red",
    "Trait covariation" = "black")) + facet_grid(R0~CV) +  
  xlab("Trait correlation") + 
  ylab("Probability of fadeout") + 
  theme_bw() + 
  theme(legend.position="bottom", legend.title=element_blank())
dev.off()

## Fig. 1
rbind(filter(csCov, R0=="R0=1"),
      filter(agCov, R0=="R0=1"),
      filter(saCov, R0=="R0=1")) -> df1
df1$cov = factor(df1$cov, 
                 levels=c("c-s", "a-g", "s-a"),
                 labels=c("Shedding/Contact",
                          "Recovery/Mortality",
                          "Shedding/Mortality"))
df1 %>% 
  group_by(CV, corr, cov) %>%
  summarize(fadeoutProb=sum(fadeout)/1000) %>%
  mutate(., line="Trait covariation") -> df1

filter(sVar, R0=="R0=1") %>% 
  group_by(CV) %>% 
  summarize(fadeoutProb=sum(fadeout)/1000) %>%
  merge(., data.frame(corr=c(seq(-0.9,-0.1,0.1),seq(0.1,0.9,0.1)))) %>% 
  mutate(., cov="Shedding/Contact", line="Shedding variation") -> df2

filter(cVar, R0=="R0=1") %>% 
  group_by(CV) %>% 
  summarize(fadeoutProb=sum(fadeout)/1000) %>%
  merge(., data.frame(corr=c(seq(-0.9,-0.1,0.1),seq(0.1,0.9,0.1)))) %>% 
  mutate(., cov="Shedding/Contact", line="Contact variation") -> df3

filter(gVar, R0=="R0=1") %>% 
  group_by(CV) %>% 
  summarize(fadeoutProb=sum(fadeout)/1000) %>%
  merge(., data.frame(corr=c(seq(-0.9,-0.1,0.1),seq(0.1,0.9,0.1)))) %>% 
  mutate(., cov="Recovery/Mortality", line="Recovery variation") -> df4

filter(aVar, R0=="R0=1") %>% 
  group_by(CV) %>% 
  summarize(fadeoutProb=sum(fadeout)/1000) %>%
  merge(., data.frame(corr=c(seq(-0.9,-0.1,0.1),seq(0.1,0.9,0.1)))) %>% 
  mutate(., cov="Recovery/Mortality", line="Mortality variation") -> df5

filter(sVar, R0=="R0=1") %>% 
  group_by(CV) %>% 
  summarize(fadeoutProb=sum(fadeout)/1000) %>%
  merge(., data.frame(corr=c(seq(-0.9,-0.1,0.1),seq(0.1,0.9,0.1)))) %>% 
  mutate(., cov="Shedding/Mortality", line="Shedding variation") -> df6

filter(aVar, R0=="R0=1") %>% 
  group_by(CV) %>% 
  summarize(fadeoutProb=sum(fadeout)/1000) %>%
  merge(., data.frame(corr=c(seq(-0.9,-0.1,0.1),seq(0.1,0.9,0.1)))) %>% 
  mutate(., cov="Shedding/Mortality", line="Mortality variation") -> df7

rbind(df1,
      df2[,colnames(df1)],
      df3[,colnames(df1)],
      df4[,colnames(df1)],
      df5[,colnames(df1)],
      df6[,colnames(df1)],
      df7[,colnames(df1)]) -> df
df$cov = factor(df$cov, levels=c("Shedding/Contact",
                                 "Recovery/Mortality",
                                 "Shedding/Mortality"))

png("Fadeout_prob_R0=1.png", height=5, width=5, units='in', res=450)          
df %>%
  ggplot(., aes(x=corr, y=fadeoutProb, group=line, color=line)) + 
  geom_line(alpha=0.7, linewidth=1) +
  scale_color_manual(values = c(
    "Shedding variation" = "#0072B2",
    "Recovery variation" = "#009E73",
    "Mortality variation" = "#D55E00",
    "Contact variation" = "#E69F00",
    "Trait covariation" = "black")) + 
  facet_grid(cov~CV) + 
  xlab("Trait correlation") + 
  ylab("Probability of fadeout") + 
  theme_bw() + 
  theme(legend.position="bottom", legend.title=element_blank()) + 
  guides(color = guide_legend(nrow=2))
dev.off()

## Fig. 2
rbind(filter(csCov, CV=="CV=5"),
      filter(agCov, CV=="CV=5"),
      filter(saCov, CV=="CV=5")) -> df1
df1$cov = factor(df1$cov, 
                 levels=c("c-s", "a-g", "s-a"),
                 labels=c("Shedding/Contact",
                          "Recovery/Mortality",
                          "Shedding/Mortality"))
df1 %>% 
  group_by(R0, corr, cov) %>%
  summarize(fadeoutProb=sum(fadeout)/1000) %>%
  mutate(., line="Trait covariation") -> df1

filter(sVar, CV=="CV=5") %>% 
  group_by(R0) %>% 
  summarize(fadeoutProb=sum(fadeout)/1000) %>%
  merge(., data.frame(corr=c(seq(-0.9,-0.1,0.1),seq(0.1,0.9,0.1)))) %>% 
  mutate(., cov="Shedding/Contact", line="Shedding variation") -> df2

filter(cVar, CV=="CV=5") %>% 
  group_by(R0) %>% 
  summarize(fadeoutProb=sum(fadeout)/1000) %>%
  merge(., data.frame(corr=c(seq(-0.9,-0.1,0.1),seq(0.1,0.9,0.1)))) %>% 
  mutate(., cov="Shedding/Contact", line="Contact variation") -> df3

filter(gVar, CV=="CV=5") %>% 
  group_by(R0) %>% 
  summarize(fadeoutProb=sum(fadeout)/1000) %>%
  merge(., data.frame(corr=c(seq(-0.9,-0.1,0.1),seq(0.1,0.9,0.1)))) %>% 
  mutate(., cov="Recovery/Mortality", line="Recovery variation") -> df4

filter(aVar, CV=="CV=5") %>% 
  group_by(R0) %>% 
  summarize(fadeoutProb=sum(fadeout)/1000) %>%
  merge(., data.frame(corr=c(seq(-0.9,-0.1,0.1),seq(0.1,0.9,0.1)))) %>% 
  mutate(., cov="Recovery/Mortality", line="Mortality variation") -> df5

filter(sVar, CV=="CV=5") %>% 
  group_by(R0) %>% 
  summarize(fadeoutProb=sum(fadeout)/1000) %>%
  merge(., data.frame(corr=c(seq(-0.9,-0.1,0.1),seq(0.1,0.9,0.1)))) %>% 
  mutate(., cov="Shedding/Mortality", line="Shedding variation") -> df6

filter(aVar, CV=="CV=5") %>% 
  group_by(R0) %>% 
  summarize(fadeoutProb=sum(fadeout)/1000) %>%
  merge(., data.frame(corr=c(seq(-0.9,-0.1,0.1),seq(0.1,0.9,0.1)))) %>% 
  mutate(., cov="Shedding/Mortality", line="Mortality variation") -> df7

rbind(df1,
      df2[,colnames(df1)],
      df3[,colnames(df1)],
      df4[,colnames(df1)],
      df5[,colnames(df1)],
      df6[,colnames(df1)],
      df7[,colnames(df1)]) -> df
df$cov = factor(df$cov, levels=c("Shedding/Contact",
                                 "Recovery/Mortality",
                                 "Shedding/Mortality"))
      
png("Fadeout_prob_CV=5.png", height=5, width=5, units='in', res=450)          
df %>%
  ggplot(., aes(x=corr, y=fadeoutProb, group=line, color=line)) + 
  geom_line(alpha=0.7, linewidth=1) +
  scale_color_manual(values = c(
    "Shedding variation" = "#0072B2",
    "Recovery variation" = "#009E73",
    "Mortality variation" = "#D55E00",
    "Contact variation" = "#E69F00",
    "Trait covariation" = "black")) + 
  facet_grid(cov~R0) + 
  xlab("Trait correlation") + 
  ylab("Probability of fadeout") + 
  theme_bw() + 
  theme(legend.position="bottom", legend.title=element_blank()) + 
  guides(color = guide_legend(nrow=2))
dev.off()



CVfig2 %>%
  group_by(corr, R0, cov) %>% 
  summarize(fadoutProb=sum(fadeout)/1000) %>% 
  ggplot(., aes(x=corr, y=fadoutProb, color="Trait covariation")) + 
  geom_line(linewidth=1) + 
  facet_grid(R0~cov) + 
  theme_bw() 




## Fig. 2 medians and interquartile range 
png(file="Peak_epidemic_size_c-s.png", width=5, height=5, units='in', res=450)
ggplot() + 
  geom_rect(data=(cVar %>% 
                    group_by(R0, CV) %>% 
                    summarize(minPeak=sort(peak)[250], maxPeak=sort(peak)[750])),
            mapping=aes(ymin=minPeak, ymax=maxPeak, xmin=-0.9, xmax=0.9),
            fill="pink",
            alpha=0.6) +
  geom_rect(data=(cVar %>% 
                    group_by(R0, CV) %>% 
                    summarize(minPeak=sort(peak)[25], maxPeak=sort(peak)[975])),
            mapping=aes(ymin=minPeak, ymax=maxPeak, xmin=-0.9, xmax=0.9),
            fill="pink",
            alpha=0.2) +
  geom_hline(data=(cVar %>% 
                     group_by(R0, CV) %>% 
                     summarize(medianPeak=median(peak))),
             mapping=aes(yintercept=medianPeak,color="Contact variation")) +
  geom_rect(data=(sVar %>% 
                    group_by(R0, CV) %>% 
                    summarize(minPeak=sort(peak)[250], maxPeak=sort(peak)[750])),
            mapping=aes(ymin=minPeak, ymax=maxPeak, xmin=-0.9, xmax=0.9),
            fill="lightblue",
            alpha=0.6) +
  geom_rect(data=(sVar %>% 
                    group_by(R0, CV) %>% 
                    summarize(minPeak=sort(peak)[25], maxPeak=sort(peak)[975])),
            mapping=aes(ymin=minPeak, ymax=maxPeak, xmin=-0.9, xmax=0.9),
            fill="lightblue",
            alpha=0.2) +
  geom_hline(data=(sVar %>% 
                     group_by(R0, CV) %>% 
                     summarize(medianPeak=median(peak))),
             mapping=aes(yintercept=medianPeak,color="Shedding variation")) +
  geom_ribbon(data=csCov %>% 
                group_by(corr, R0, CV) %>% 
                summarize(minPeak=sort(peak)[250], maxPeak=sort(peak)[750]),
              aes(x=corr, ymin=minPeak, ymax=maxPeak),
              fill="grey",
              alpha=0.6) +
  geom_ribbon(data=csCov %>% 
                group_by(corr, R0, CV) %>% 
                summarize(minPeak=sort(peak)[25], maxPeak=sort(peak)[975]),
              aes(x=corr, ymin=minPeak, ymax=maxPeak),
              fill="grey",
              alpha=0.2) +
  geom_line(data=csCov %>% 
              group_by(corr, R0, CV) %>% 
              summarize(medianPeak=median(peak)),
            aes(x=corr, y=medianPeak, color="Trait covariation")) + 
  scale_color_manual(values = c(
    "Contact variation" = "red",
    "Shedding variation" = "blue",
    "Trait covariation" = "black")) + 
  facet_grid(R0~CV) + 
  xlab("Trait correlation") + 
  ylab("Peak no. of infecteds") + 
  theme_bw() + 
  theme(legend.position="bottom", legend.title=element_blank())
dev.off()

png(file="Peak_epidemic_size_a-g.png", width=5, height=5, units='in', res=450)
ggplot() + 
  geom_rect(data=(aVar %>% 
                    group_by(R0, CV) %>% 
                    summarize(minPeak=sort(peak)[250], maxPeak=sort(peak)[750])),
            mapping=aes(ymin=minPeak, ymax=maxPeak, xmin=-0.9, xmax=0.9),
            fill="pink",
            alpha=0.6) +
  geom_rect(data=(aVar %>% 
                    group_by(R0, CV) %>% 
                    summarize(minPeak=sort(peak)[25], maxPeak=sort(peak)[975])),
            mapping=aes(ymin=minPeak, ymax=maxPeak, xmin=-0.9, xmax=0.9),
            fill="pink",
            alpha=0.2) +
  geom_hline(data=(aVar %>% 
                     group_by(R0, CV) %>% 
                     summarize(medianPeak=median(peak))),
             mapping=aes(yintercept=medianPeak,color="Mortality variation")) +
  geom_rect(data=(gVar %>% 
                    group_by(R0, CV) %>% 
                    summarize(minPeak=sort(peak)[250], maxPeak=sort(peak)[750])),
            mapping=aes(ymin=minPeak, ymax=maxPeak, xmin=-0.9, xmax=0.9),
            fill="lightblue",
            alpha=0.6) +
  geom_rect(data=(gVar %>% 
                    group_by(R0, CV) %>% 
                    summarize(minPeak=sort(peak)[25], maxPeak=sort(peak)[975])),
            mapping=aes(ymin=minPeak, ymax=maxPeak, xmin=-0.9, xmax=0.9),
            fill="lightblue",
            alpha=0.2) +
  geom_hline(data=(gVar %>% 
                     group_by(R0, CV) %>% 
                     summarize(medianPeak=median(peak))),
             mapping=aes(yintercept=medianPeak,color="Recovery variation")) +
  geom_ribbon(data=agCov %>% 
                group_by(corr, R0, CV) %>% 
                summarize(minPeak=sort(peak)[250], maxPeak=sort(peak)[750]),
              aes(x=corr, ymin=minPeak, ymax=maxPeak),
              fill="grey",
              alpha=0.6) +
  geom_ribbon(data=agCov %>% 
                group_by(corr, R0, CV) %>% 
                summarize(minPeak=sort(peak)[25], maxPeak=sort(peak)[975]),
              aes(x=corr, ymin=minPeak, ymax=maxPeak),
              fill="grey",
              alpha=0.2) +
  geom_line(data=agCov %>% 
              group_by(corr, R0, CV) %>% 
              summarize(medianPeak=median(peak)),
            aes(x=corr, y=medianPeak, color="Trait covariation")) + 
  scale_color_manual(values = c(
    "Mortality variation" = "red",
    "Recovery variation" = "blue",
    "Trait covariation" = "black")) + 
  facet_grid(R0~CV) + 
  xlab("Trait correlation") + 
  ylab("Peak no. of infecteds") + 
  theme_bw() + 
  theme(legend.position="bottom", legend.title=element_blank())
dev.off()

png(file="Peak_epidemic_size_s-a.png", width=5, height=5, units='in', res=450)
ggplot() + 
  geom_rect(data=(sVar %>% 
                    group_by(R0, CV) %>% 
                    summarize(minPeak=sort(peak)[250], maxPeak=sort(peak)[750])),
            mapping=aes(ymin=minPeak, ymax=maxPeak, xmin=-0.9, xmax=0.9),
            fill="lightblue",
            alpha=0.6) +
  geom_rect(data=(sVar %>% 
                    group_by(R0, CV) %>% 
                    summarize(minPeak=sort(peak)[25], maxPeak=sort(peak)[975])),
            mapping=aes(ymin=minPeak, ymax=maxPeak, xmin=-0.9, xmax=0.9),
            fill="lightblue",
            alpha=0.2) +
  geom_hline(data=(sVar %>% 
                     group_by(R0, CV) %>% 
                     summarize(medianPeak=median(peak))),
             mapping=aes(yintercept=medianPeak,color="Shedding variation")) +
  geom_rect(data=(aVar %>% 
                    group_by(R0, CV) %>% 
                    summarize(minPeak=sort(peak)[250], maxPeak=sort(peak)[750])),
            mapping=aes(ymin=minPeak, ymax=maxPeak, xmin=-0.9, xmax=0.9),
            fill="pink",
            alpha=0.6) +
  geom_rect(data=(aVar %>% 
                    group_by(R0, CV) %>% 
                    summarize(minPeak=sort(peak)[25], maxPeak=sort(peak)[975])),
            mapping=aes(ymin=minPeak, ymax=maxPeak, xmin=-0.9, xmax=0.9),
            fill="pink",
            alpha=0.2) +
  geom_hline(data=(aVar %>% 
                     group_by(R0, CV) %>% 
                     summarize(medianPeak=median(peak))),
             mapping=aes(yintercept=medianPeak,color="Mortality variation")) +
  geom_ribbon(data=saCov %>% 
                group_by(corr, R0, CV) %>% 
                summarize(minPeak=sort(peak)[250], maxPeak=sort(peak)[750]),
              aes(x=corr, ymin=minPeak, ymax=maxPeak),
              fill="grey",
              alpha=0.6) +
  geom_ribbon(data=saCov %>% 
                group_by(corr, R0, CV) %>% 
                summarize(minPeak=sort(peak)[25], maxPeak=sort(peak)[975]),
              aes(x=corr, ymin=minPeak, ymax=maxPeak),
              fill="grey",
              alpha=0.2) +
  geom_line(data=saCov %>% 
              group_by(corr, R0, CV) %>% 
              summarize(medianPeak=median(peak)),
            aes(x=corr, y=medianPeak, color="Trait covariation")) + 
  scale_color_manual(values = c(
    "Shedding variation" = "blue",
    "Mortality variation" = "red",
    "Trait covariation" = "black")) + 
  facet_grid(R0~CV) + 
  xlab("Trait correlation") + 
  ylab("Peak no. of infecteds") + 
  theme_bw() + 
  theme(legend.position="bottom", legend.title=element_blank())
dev.off()



## Do stats
statsy = c()
for (r in c(1, 4, 8)) {
  for (cv in c(0.2, 1, 5)) {
    for (cor in c(seq(-0.9, -0.1, 0.1), seq(0.1, 0.9, 0.1))) {
      cor = round(cor,1)
      peakTTest = t.test(filter(cVar, R0==r, CV==cv)$peak, filter(csCov, R0==r, CV==cv, corr==cor)$peak)
      dispTTest = try(t.test(filter(cVar, R0==r, CV==cv)$disp, filter(csCov, R0==r, CV==cv, corr==cor)$disp))
      if(inherits(dispTTest,"try-error"))
        statsy = rbind(statsy, c(r, cv, cor, peakTTest$estimate, peakTTest$p.value, mean(filter(cVar, R0==r, CV==cv)$disp, na.rm=T), mean(filter(csCov, R0==r, CV==cv, corr==cor)$disp,na.rm=T), NA))
      else
        statsy = rbind(statsy, c(r, cv, cor, peakTTest$estimate, peakTTest$p.value, dispTTest$estimate, dispTTest$p.value))
    }
  }  
}




