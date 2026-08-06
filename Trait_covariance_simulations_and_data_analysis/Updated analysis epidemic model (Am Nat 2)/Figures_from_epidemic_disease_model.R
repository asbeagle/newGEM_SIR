library(tidyverse)
library(ggplot2)
library(magrittr)

out = list()
i = 1
for (pmtr in c("c","s","a","g")) {
  for (R0 in c(1,2,4,8)) {
    for (CV in c(0.2, 1, 5)) { 
      out[[i]] = readRDS(file=paste0("out_R0=",R0,"_S=1000_CV=",CV,"_var=",pmtr,".RDS"))
      i = i + 1
    }
  }
}
out = do.call("rbind.data.frame", out)
out$CV = paste0("CV==",out$CV)
out$R0 = paste0("R[0]==",out$R0)

out2 = list()
i = 1
for (pmtr in c("c-s","c-a","c-g","s-a","s-g","a-g")) {
  for (R0 in c(1,2,4,8)) {
    for (CV in c(0.2, 1, 5)) { 
      for (corr in seq(-0.9,0.9,0.1)) { 
        out2[[i]] = readRDS(file=paste0("out_R0=",R0,"_S=1000_CV=",CV,"_corr=",corr,"_cov=",pmtr,".RDS"))
        i = i + 1
      }
    }
  }
}
out2 = do.call("rbind.data.frame", out2)
out2$CV = paste0("CV==",out2$CV)
out2$R0 = paste0("R[0]==",out2$R0)

## There are two ways to quantify the severity of an epidemic. 
## One is the peak number of infected individuals.
## The other is how long the epidemic lasts, which gives a sense of how "explosive" the epidemic was.
var1 = c("c","c","c","s","s","a")
var2 = c("s","a","g","a","g","g")
for (i in 1:6) {
  name1 = ifelse(var1[i]=="c","Contact variation", ifelse(var1[i]=="s","Shedding variation","Mortality variation"))
  name2 = ifelse(var2[i]=="s","Shedding variation", ifelse(var2[i]=="a","Mortality variation","Recovery variation"))
  
  mycols <- c(
    "Contact variation" = "#0072B2", # blue
    "Shedding variation" = "#D55E00", # vermillion
    "Mortality variation" = "#009E73", # bluish green
    "Recovery variation" = "#CC79A7", # reddish purple
    "Trait covariation" = "#4D4D4D" # dark grey
  )
  
  png(file=paste0("Epidemic_size_",var1[i],"-",var2[i],".png"), width=5, height=6, units='in', res=450)
  p = ggplot() + 
    geom_rect(data=(filter(out, var==var1[i]) %>% 
                      group_by(R0, CV) %>% 
                      summarize(minPeak=sort(peak)[250], maxPeak=sort(peak)[750])),
              mapping=aes(ymin=minPeak, ymax=maxPeak, xmin=-0.9, xmax=0.9),
              fill=mycols[name1],
              alpha=0.4) +
    geom_hline(data=(filter(out, var==var1[i]) %>% 
                       group_by(R0, CV) %>% 
                       summarize(medianPeak=median(peak))),
               mapping=aes(yintercept=medianPeak,color=name1)) +
    geom_rect(data=(filter(out, var==var2[i]) %>% 
                      group_by(R0, CV) %>% 
                      summarize(minPeak=sort(peak)[250], maxPeak=sort(peak)[750])),
              mapping=aes(ymin=minPeak, ymax=maxPeak, xmin=-0.9, xmax=0.9),
              fill=mycols[name2],
              alpha=0.4) +
    geom_hline(data=(filter(out, var==var2[i]) %>% 
                       group_by(R0, CV) %>% 
                       summarize(medianPeak=median(peak))),
               mapping=aes(yintercept=medianPeak,color=name2)) +
    geom_ribbon(data=filter(out2, cov==paste0(var1[i],"-",var2[i])) %>% 
                  group_by(corr, R0, CV) %>% 
                  summarize(minPeak=sort(peak)[250], maxPeak=sort(peak)[750]),
                aes(x=corr, ymin=minPeak, ymax=maxPeak),
                fill=mycols["Trait covariation"],
                alpha=0.4) +
    geom_line(data=filter(out2, cov==paste0(var1[i],"-",var2[i])) %>% 
                group_by(corr, R0, CV) %>% 
                summarize(medianPeak=median(peak)),
              aes(x=corr, y=medianPeak, color="Trait covariation")) + 
    scale_color_manual(values = mycols) + 
    facet_grid(R0~CV, labeller=label_parsed) + 
    xlab("Trait correlation") + 
    ylab("Peak no. of infecteds") + 
    theme_bw() + 
    theme(legend.position="bottom", legend.title=element_blank())
  print(p)
  dev.off()
}  
  
for (i in 1:6) {
  name1 = ifelse(var1[i]=="c","Contact variation", ifelse(var1[i]=="s","Shedding variation","Mortality variation"))
  name2 = ifelse(var2[i]=="s","Shedding variation", ifelse(var2[i]=="a","Mortality variation","Recovery variation"))
  
  mycols <- c(
    "Contact variation" = "#0072B2", # blue
    "Shedding variation" = "#D55E00", # vermillion
    "Mortality variation" = "#009E73", # bluish green
    "Recovery variation" = "#CC79A7", # reddish purple
    "Trait covariation" = "#4D4D4D" # dark grey
  )
  
  png(file=paste0("Epidemic_duration_",var1[i],"-",var2[i],".png"), width=5, height=6, units='in', res=450)
  p = ggplot() + 
    geom_rect(data=(filter(out, var==var1[i]) %>% 
                      group_by(R0, CV) %>% 
                      summarize(minDuration=sort(tEnd)[250], maxDuration=sort(tEnd)[750])),
              mapping=aes(ymin=minDuration, ymax=maxDuration, xmin=-0.9, xmax=0.9),
              fill=mycols[name1],
              alpha=0.4) +
    geom_hline(data=(filter(out, var==var1[i]) %>% 
                       group_by(R0, CV) %>% 
                       summarize(medianDuration=median(tEnd))),
               mapping=aes(yintercept=medianDuration,color=name1)) +
    geom_rect(data=(filter(out, var==var2[i]) %>% 
                      group_by(R0, CV) %>% 
                      summarize(minDuration=sort(tEnd)[250], maxDuration=sort(tEnd)[750])),
              mapping=aes(ymin=minDuration, ymax=maxDuration, xmin=-0.9, xmax=0.9),
              fill=mycols[name2],
              alpha=0.4) +
    geom_hline(data=(filter(out, var==var2[i]) %>% 
                       group_by(R0, CV) %>% 
                       summarize(medianDuration=median(tEnd))),
               mapping=aes(yintercept=medianDuration,color=name2)) +
    geom_ribbon(data=filter(out2, cov==paste0(var1[i],"-",var2[i])) %>% 
                  group_by(corr, R0, CV) %>% 
                  summarize(minDuration=sort(tEnd)[250], maxDuration=sort(tEnd)[750]),
                aes(x=corr, ymin=minDuration, ymax=maxDuration),
                fill=mycols["Trait covariation"],
                alpha=0.4) +
    geom_line(data=filter(out2, cov==paste0(var1[i],"-",var2[i])) %>% 
                group_by(corr, R0, CV) %>% 
                summarize(medianDuration=median(tEnd)),
              aes(x=corr, y=medianDuration, color="Trait covariation")) + 
    scale_color_manual(values = mycols) + 
    facet_grid(R0~CV, labeller=label_parsed) + 
    xlab("Trait correlation") + 
    ylab("Epidemic duration") + 
    theme_bw() + 
    theme(legend.position="bottom", legend.title=element_blank())
  print(p)
  dev.off()
}  

## Superspreading can be quantified using the estimated dispersion parameter
## Dispersion parameter estimates
for (i in 1:6) {
  name1 = ifelse(var1[i]=="c","Contact variation", ifelse(var1[i]=="s","Shedding variation","Mortality variation"))
  name2 = ifelse(var2[i]=="s","Shedding variation", ifelse(var2[i]=="a","Mortality variation","Recovery variation"))
  
  mycols <- c(
    "Contact variation" = "#0072B2", # blue
    "Shedding variation" = "#D55E00", # vermillion
    "Mortality variation" = "#009E73", # bluish green
    "Recovery variation" = "#CC79A7", # reddish purple
    "Trait covariation" = "#4D4D4D" # dark grey
  )
  
  png(file=paste0("Dispersion_parameter_",var1[i],"-",var2[i],".png"), width=5, height=6, units='in', res=450)
  p = ggplot() + 
    geom_rect(data=(filter(out, var==var1[i], !is.na(disp)) %>% 
                      group_by(R0, CV) %>% 
                      summarize(minDispersion=quantile(disp,0.25), maxDispersion=quantile(disp,0.75))),
              mapping=aes(ymin=minDispersion, ymax=maxDispersion, xmin=-0.9, xmax=0.9),
              fill=mycols[name1],
              alpha=0.4) +
    geom_hline(data=(filter(out, var==var1[i], !is.na(disp)) %>% 
                       group_by(R0, CV) %>% 
                       summarize(medianDispersion=median(disp))),
               mapping=aes(yintercept=medianDispersion,color=name1)) +
    geom_rect(data=(filter(out, var==var2[i], !is.na(disp)) %>% 
                      group_by(R0, CV) %>% 
                      summarize(minDispersion=quantile(disp,0.25), maxDispersion=quantile(disp,0.75))),
              mapping=aes(ymin=minDispersion, ymax=maxDispersion, xmin=-0.9, xmax=0.9),
              fill=mycols[name2],
              alpha=0.4) +
    geom_hline(data=(filter(out, var==var2[i], !is.na(disp)) %>% 
                       group_by(R0, CV) %>% 
                       summarize(medianDispersion=median(disp))),
               mapping=aes(yintercept=medianDispersion,color=name2)) +
    geom_ribbon(data=(filter(out2, cov==paste0(var1[i],"-",var2[i]), !is.na(disp)) %>% 
                  group_by(corr, R0, CV) %>% 
                  summarize(minDispersion=quantile(disp,0.25), maxDispersion=quantile(disp,0.75))),
                aes(x=corr, ymin=minDispersion, ymax=maxDispersion),
                fill=mycols["Trait covariation"],
                alpha=0.4) +
    geom_line(data=filter(out2, cov==paste0(var1[i],"-",var2[i]), !is.na(disp)) %>% 
                group_by(corr, R0, CV) %>% 
                summarize(medianDispersion=median(disp)),
              aes(x=corr, y=medianDispersion, color="Trait covariation")) + 
    scale_color_manual(values = mycols) + 
    facet_grid(R0~CV, labeller=label_parsed) + 
    ylim(0,2) + 
    xlab("Trait correlation") + 
    ylab("Estimated dispersion parameter") + 
    theme_bw() + 
    theme(legend.position="bottom", legend.title=element_blank())
  print(p)
  dev.off()
}  

## A different way of defining superspreading is via the proportion of transmission due to the most infectious 20% of cases (a la fig. 1c in Lloyd-Smith et al)
for (i in 1:6) {
  name1 = ifelse(var1[i]=="c","Contact variation", ifelse(var1[i]=="s","Shedding variation","Mortality variation"))
  name2 = ifelse(var2[i]=="s","Shedding variation", ifelse(var2[i]=="a","Mortality variation","Recovery variation"))
  
  mycols <- c(
    "Contact variation" = "#0072B2", # blue
    "Shedding variation" = "#D55E00", # vermillion
    "Mortality variation" = "#009E73", # bluish green
    "Recovery variation" = "#CC79A7", # reddish purple
    "Trait covariation" = "#4D4D4D" # dark grey
  )
  
  png(file=paste0("Epidemic_proportion_top_20_",var1[i],"-",var2[i],".png"), width=5, height=6, units='in', res=450)
  p = ggplot() + 
    geom_rect(data=(filter(out, var==var1[i], !is.na(top20)) %>% 
                      group_by(R0, CV) %>% 
                      summarize(minProportion=quantile(top20,0.25), maxProportion=quantile(top20,0.75))),
              mapping=aes(ymin=minProportion, ymax=maxProportion, xmin=-0.9, xmax=0.9),
              fill=mycols[name1],
              alpha=0.4) +
    geom_hline(data=(filter(out, var==var1[i], !is.na(top20)) %>% 
                       group_by(R0, CV) %>% 
                       summarize(medianProportion=median(top20))),
               mapping=aes(yintercept=medianProportion,color=name1)) +
    geom_rect(data=(filter(out, var==var2[i], !is.na(top20)) %>% 
                      group_by(R0, CV) %>% 
                      summarize(minProportion=quantile(top20,0.25), maxProportion=quantile(top20,0.75))),
              mapping=aes(ymin=minProportion, ymax=maxProportion, xmin=-0.9, xmax=0.9),
              fill=mycols[name2],
              alpha=0.4) +
    geom_hline(data=(filter(out, var==var2[i], !is.na(top20)) %>% 
                       group_by(R0, CV) %>% 
                       summarize(medianProportion=median(top20))),
               mapping=aes(yintercept=medianProportion,color=name2)) +
    geom_ribbon(data=(filter(out2, cov==paste0(var1[i],"-",var2[i]), !is.na(top20)) %>% 
                        group_by(corr, R0, CV) %>% 
                        summarize(minProportion=quantile(top20,0.25), maxProportion=quantile(top20,0.75))),
                aes(x=corr, ymin=minProportion, ymax=maxProportion),
                fill=mycols["Trait covariation"],
                alpha=0.4) +
    geom_line(data=filter(out2, cov==paste0(var1[i],"-",var2[i]), !is.na(top20)) %>% 
                group_by(corr, R0, CV) %>% 
                summarize(medianProportion=median(top20)),
              aes(x=corr, y=medianProportion, color="Trait covariation")) + 
    scale_color_manual(values = mycols) + 
    facet_grid(R0~CV, labeller=label_parsed) + 
    xlab("Trait correlation") + 
    ylab("Proportion of transmission due to\nmost infectious 20% of cases") + 
    theme_bw() + 
    theme(legend.position="bottom", legend.title=element_blank())
  print(p)
  dev.off()
}  

## Epidemic sizes
list(
  filter(out, var=="c") %>% 
    group_by(CV, R0) %>% 
    summarize(medianPeak=median(peak), 
              minPeak=quantile(peak,0.25), 
              maxPeak=quantile(peak,0.75),
              medianTop20=median(top20[peak>=50]),
              minTop20=quantile(top20[peak>=50],0.25),
              maxTop20=quantile(top20[peak>=50],0.75),
              medianDisp=median(disp[peak>=50]),
              minDisp=quantile(disp[peak>=50],0.25),
              maxDisp=quantile(disp[peak>=50],0.75),
              fadeoutProb=sum(peak<50)/1000) %>% 
    merge(., data.frame(corr=seq(-0.9,0.9,0.1))) %>% 
    mutate(., 
           covName="Shedding/Contact",
           varName="Contact Variation"),
  filter(out, var=="s") %>% 
    group_by(CV, R0) %>% 
    summarize(medianPeak=median(peak), 
              minPeak=quantile(peak,0.25), 
              maxPeak=quantile(peak,0.75),
              medianTop20=median(top20[peak>=50]),
              minTop20=quantile(top20[peak>=50],0.25),
              maxTop20=quantile(top20[peak>=50],0.75),
              medianDisp=median(disp[peak>=50]),
              minDisp=quantile(disp[peak>=50],0.25),
              maxDisp=quantile(disp[peak>=50],0.75),
              fadeoutProb=sum(peak<50)/1000) %>% 
    merge(., data.frame(corr=seq(-0.9,0.9,0.1))) %>% 
    mutate(., 
           covName="Shedding/Contact",
           varName="Shedding Variation"),
  filter(out, var=="s") %>% 
    group_by(CV, R0) %>% 
    summarize(medianPeak=median(peak), 
              minPeak=quantile(peak,0.25), 
              maxPeak=quantile(peak,0.75),
              medianTop20=median(top20[peak>=50]),
              minTop20=quantile(top20[peak>=50],0.25),
              maxTop20=quantile(top20[peak>=50],0.75),
              medianDisp=median(disp[peak>=50]),
              minDisp=quantile(disp[peak>=50],0.25),
              maxDisp=quantile(disp[peak>=50],0.75),
              fadeoutProb=sum(peak<50)/1000) %>% 
    merge(., data.frame(corr=seq(-0.9,0.9,0.1))) %>% 
    mutate(., 
           covName="Shedding/Mortality",
           varName="Shedding Variation"),
  filter(out, var=="a") %>% 
    group_by(CV, R0) %>% 
    summarize(medianPeak=median(peak), 
              minPeak=quantile(peak,0.25), 
              maxPeak=quantile(peak,0.75),
              medianTop20=median(top20[peak>=50]),
              minTop20=quantile(top20[peak>=50],0.25),
              maxTop20=quantile(top20[peak>=50],0.75),
              medianDisp=median(disp[peak>=50]),
              minDisp=quantile(disp[peak>=50],0.25),
              maxDisp=quantile(disp[peak>=50],0.75),
              fadeoutProb=sum(peak<50)/1000) %>% 
    merge(., data.frame(corr=seq(-0.9,0.9,0.1))) %>% 
    mutate(., 
           covName="Shedding/Mortality",
           varName="Mortality Variation"),
  filter(out, var=="a") %>% 
    group_by(CV, R0) %>% 
    summarize(medianPeak=median(peak), 
              minPeak=quantile(peak,0.25), 
              maxPeak=quantile(peak,0.75),
              medianTop20=median(top20[peak>=50]),
              minTop20=quantile(top20[peak>=50],0.25),
              maxTop20=quantile(top20[peak>=50],0.75),
              medianDisp=median(disp[peak>=50]),
              minDisp=quantile(disp[peak>=50],0.25),
              maxDisp=quantile(disp[peak>=50],0.75),
              fadeoutProb=sum(peak<50)/1000) %>% 
    merge(., data.frame(corr=seq(-0.9,0.9,0.1))) %>% 
    mutate(., 
           covName="Recovery/Mortality",
           varName="Mortality Variation"),
  filter(out, var=="g") %>% 
    group_by(CV, R0) %>% 
    summarize(medianPeak=median(peak), 
              minPeak=quantile(peak,0.25), 
              maxPeak=quantile(peak,0.75),
              medianTop20=median(top20[peak>=50]),
              minTop20=quantile(top20[peak>=50],0.25),
              maxTop20=quantile(top20[peak>=50],0.75),
              medianDisp=median(disp[peak>=50]),
              minDisp=quantile(disp[peak>=50],0.25),
              maxDisp=quantile(disp[peak>=50],0.75),
              fadeoutProb=sum(peak<50)/1000) %>% 
    merge(., data.frame(corr=seq(-0.9,0.9,0.1))) %>% 
    mutate(., 
           covName="Recovery/Mortality",
           varName="Recovery Variation")) %>% 
  do.call("rbind.data.frame",.) -> df1

list(
  filter(out2, cov=="a-g") %>% 
    group_by(CV, R0, corr) %>% 
    summarize(medianPeak=median(peak), 
              minPeak=quantile(peak,0.25), 
              maxPeak=quantile(peak,0.75),
              medianTop20=median(top20[peak>=50]),
              minTop20=quantile(top20[peak>=50],0.25),
              maxTop20=quantile(top20[peak>=50],0.75),
              medianDisp=median(disp[peak>=50]),
              minDisp=quantile(disp[peak>=50],0.25),
              maxDisp=quantile(disp[peak>=50],0.75),
              fadeoutProb=sum(peak<50)/1000) %>% 
    mutate(.,
           covName="Recovery/Mortality",
           varName="Trait Covariation"),
  filter(out2, cov=="c-s") %>% 
    group_by(CV, R0, corr) %>% 
    summarize(medianPeak=median(peak), 
              minPeak=quantile(peak,0.25), 
              maxPeak=quantile(peak,0.75),
              medianTop20=median(top20[peak>=50]),
              minTop20=quantile(top20[peak>=50],0.25),
              maxTop20=quantile(top20[peak>=50],0.75),
              medianDisp=median(disp[peak>=50]),
              minDisp=quantile(disp[peak>=50],0.25),
              maxDisp=quantile(disp[peak>=50],0.75),
              fadeoutProb=sum(peak<50)/1000) %>% 
    mutate(.,
           covName="Shedding/Contact",
           varName="Trait Covariation"),
  filter(out2, cov=="s-a") %>% 
    group_by(CV, R0, corr) %>% 
    summarize(medianPeak=median(peak), 
              minPeak=quantile(peak,0.25), 
              maxPeak=quantile(peak,0.75),
              medianTop20=median(top20[peak>=50]),
              minTop20=quantile(top20[peak>=50],0.25),
              maxTop20=quantile(top20[peak>=50],0.75),
              medianDisp=median(disp[peak>=50]),
              minDisp=quantile(disp[peak>=50],0.25),
              maxDisp=quantile(disp[peak>=50],0.75),
              fadeoutProb=sum(peak<50)/1000) %>% 
    mutate(.,
           covName="Shedding/Mortality",
           varName="Trait Covariation")
) %>% do.call("rbind.data.frame", .) -> df2

combo = rbind(df1[colnames(df2)], df2)
combo$covName = factor(combo$covName, levels=c("Shedding/Contact", "Recovery/Mortality", "Shedding/Mortality"))
  
png(file="Fig1_peak_prevalence_R0=1.png", height=6, width=5.5, units='in', res=450)
combo %>% 
  filter(., R0=="R[0]==1") %>%
  ggplot(., aes(x=corr, y=medianPeak/1000, color=varName)) +
  geom_line() + 
  geom_ribbon(aes(ymin=minPeak/1000, ymax=maxPeak/1000, fill=varName), linetype="blank", alpha=0.2) + 
  facet_grid(covName~CV, labeller=label_parsed) + 
  xlab("Trait Correlation") + 
  ylab("Peak Prevalence") + 
  theme_bw() + 
  theme(legend.position="bottom", legend.title=element_blank()) + 
  guides(color = guide_legend(nrow=2))
dev.off()

png(file="Fig2_superspreading_R0=1.png", height=6, width=5.5, units='in', res=450)
combo %>% 
  filter(., R0=="R[0]==1") %>%
  ggplot(., aes(x=corr, y=medianTop20, color=varName)) +
  geom_line(aes(linetype="Transmission From Top 20% of Cases")) + 
  geom_ribbon(aes(ymin=minTop20, ymax=maxTop20, fill=varName), alpha=0.2, linetype='blank') + 
  geom_line(aes(x=corr, y=fadeoutProb, color=varName, linetype="Epidemic Fadeout")) + 
  scale_linetype_manual(values = c("Transmission From Top 20% of Cases" = 1,
                                   "Epidemic Fadeout" = 2)) + 
  facet_grid(covName~CV, labeller=label_parsed) + 
  xlab("Trait Correlation") + 
  ylab("Proportion") + 
  theme_bw() + 
  theme(legend.position="bottom", legend.title=element_blank(), legend.box="vertical") + 
  guides(color = guide_legend(nrow=2), 
         linetype = guide_legend(nrow=1))
dev.off()

png(file="Fig3_peak_prevalence_R0=2.png", height=6, width=5.5, units='in', res=450)
combo %>% 
  filter(., R0=="R[0]==2") %>%
  ggplot(., aes(x=corr, y=medianPeak/1000, color=varName)) +
  geom_line() + 
  geom_ribbon(aes(ymin=minPeak/1000, ymax=maxPeak/1000, fill=varName), linetype="blank", alpha=0.2) + 
  facet_grid(covName~CV, labeller=label_parsed) + 
  xlab("Trait Correlation") + 
  ylab("Peak Prevalence") + 
  theme_bw() + 
  theme(legend.position="bottom", legend.title=element_blank()) + 
  guides(color = guide_legend(nrow=2))
dev.off()

png(file="Fig3_peak_prevalence_R0=4.png", height=6, width=5.5, units='in', res=450)
combo %>% 
  filter(., R0=="R[0]==4") %>%
  ggplot(., aes(x=corr, y=medianPeak/1000, color=varName)) +
  geom_line() + 
  geom_ribbon(aes(ymin=minPeak/1000, ymax=maxPeak/1000, fill=varName), linetype="blank", alpha=0.2) + 
  facet_grid(covName~CV, labeller=label_parsed) + 
  xlab("Trait Correlation") + 
  ylab("Peak Prevalence") + 
  theme_bw() + 
  theme(legend.position="bottom", legend.title=element_blank()) + 
  guides(color = guide_legend(nrow=2))
dev.off()

png(file="Fig3_peak_prevalence_R0=8.png", height=6, width=5.5, units='in', res=450)
combo %>% 
  filter(., R0=="R[0]==8") %>%
  ggplot(., aes(x=corr, y=medianPeak/1000, color=varName)) +
  geom_line() + 
  geom_ribbon(aes(ymin=minPeak/1000, ymax=maxPeak/1000, fill=varName), linetype="blank", alpha=0.2) + 
  facet_grid(covName~CV, labeller=label_parsed) + 
  xlab("Trait Correlation") + 
  ylab("Peak Prevalence") + 
  theme_bw() + 
  theme(legend.position="bottom", legend.title=element_blank()) + 
  guides(color = guide_legend(nrow=2))
dev.off()

png(file="Fig3_peak_prevalence_CV=1.png", height=6, width=5.5, units='in', res=450)
combo %>% 
  filter(., CV=="CV==1") %>%
  ggplot(., aes(x=corr, y=medianPeak/1000, color=varName)) +
  geom_line() + 
  geom_ribbon(aes(ymin=minPeak/1000, ymax=maxPeak/1000, fill=varName), linetype="blank", alpha=0.2) + 
  facet_grid(covName~R0, labeller=label_parsed) + 
  xlab("Trait Correlation") + 
  ylab("Peak Prevalence") + 
  theme_bw() + 
  theme(legend.position="bottom", legend.title=element_blank()) + 
  guides(color = guide_legend(nrow=2))
dev.off()

png(file="Fig3_peak_prevalence_CV=5.png", height=6, width=5.5, units='in', res=450)
combo %>% 
  filter(., CV=="CV==5") %>%
  ggplot(., aes(x=corr, y=medianPeak/1000, color=varName)) +
  geom_line() + 
  geom_ribbon(aes(ymin=minPeak/1000, ymax=maxPeak/1000, fill=varName), linetype="blank", alpha=0.2) + 
  facet_grid(covName~R0, labeller=label_parsed) + 
  xlab("Trait Correlation") + 
  ylab("Peak Prevalence") + 
  theme_bw() + 
  theme(legend.position="bottom", legend.title=element_blank()) + 
  guides(color = guide_legend(nrow=2))
dev.off()

png(file="Fig4_superspreading_R0=4.png", height=6, width=5.5, units='in', res=450)
combo %>% 
  filter(., R0=="R[0]==4") %>%
  ggplot(., aes(x=corr, y=medianTop20, color=varName)) +
  geom_line(aes(linetype="Transmission From Top 20% of Cases")) + 
  geom_ribbon(aes(ymin=minTop20, ymax=maxTop20, fill=varName), alpha=0.2, linetype='blank') + 
  geom_line(aes(x=corr, y=fadeoutProb, color=varName, linetype="Epidemic Fadeout")) + 
  scale_linetype_manual(values = c("Transmission From Top 20% of Cases" = 1,
                                   "Epidemic Fadeout" = 2)) + 
  facet_grid(covName~CV, labeller=label_parsed) + 
  xlab("Trait Correlation") + 
  ylab("Proportion") + 
  theme_bw() + 
  theme(legend.position="bottom", legend.title=element_blank(), legend.box="vertical") + 
  guides(color = guide_legend(nrow=2), 
         linetype = guide_legend(nrow=1))
dev.off()

png(file="Fig4_superspreading_R0=2.png", height=6, width=5.5, units='in', res=450)
combo %>% 
  filter(., R0=="R[0]==2") %>%
  ggplot(., aes(x=corr, y=medianTop20, color=varName)) +
  geom_line(aes(linetype="Transmission From Top 20% of Cases")) + 
  geom_ribbon(aes(ymin=minTop20, ymax=maxTop20, fill=varName), alpha=0.2, linetype='blank') + 
  geom_line(aes(x=corr, y=fadeoutProb, color=varName, linetype="Epidemic Fadeout")) + 
  scale_linetype_manual(values = c("Transmission From Top 20% of Cases" = 1,
                                   "Epidemic Fadeout" = 2)) + 
  facet_grid(covName~CV, labeller=label_parsed) + 
  xlab("Trait Correlation") + 
  ylab("Proportion") + 
  theme_bw() + 
  theme(legend.position="bottom", legend.title=element_blank(), legend.box="vertical") + 
  guides(color = guide_legend(nrow=2), 
         linetype = guide_legend(nrow=1))
dev.off()

png(file="Fig4_superspreading_R0=8.png", height=6, width=5.5, units='in', res=450)
combo %>% 
  filter(., R0=="R[0]==8") %>%
  ggplot(., aes(x=corr, y=medianTop20, color=varName)) +
  geom_line(aes(linetype="Transmission From Top 20% of Cases")) + 
  geom_ribbon(aes(ymin=minTop20, ymax=maxTop20, fill=varName), alpha=0.2, linetype='blank') + 
  geom_line(aes(x=corr, y=fadeoutProb, color=varName, linetype="Epidemic Fadeout")) + 
  scale_linetype_manual(values = c("Transmission From Top 20% of Cases" = 1,
                                   "Epidemic Fadeout" = 2)) + 
  facet_grid(covName~CV, labeller=label_parsed) + 
  xlab("Trait Correlation") + 
  ylab("Proportion") + 
  theme_bw() + 
  theme(legend.position="bottom", legend.title=element_blank(), legend.box="vertical") + 
  guides(color = guide_legend(nrow=2), 
         linetype = guide_legend(nrow=1))
dev.off()

png(file="Fig4_superspreading_CV=1.png", height=6, width=5.5, units='in', res=450)
combo %>% 
  filter(., CV=="CV==1") %>%
  ggplot(., aes(x=corr, y=medianTop20, color=varName)) +
  geom_line(aes(linetype="Transmission From Top 20% of Cases")) + 
  geom_ribbon(aes(ymin=minTop20, ymax=maxTop20, fill=varName), alpha=0.2, linetype='blank') + 
  geom_line(aes(x=corr, y=fadeoutProb, color=varName, linetype="Epidemic Fadeout")) + 
  scale_linetype_manual(values = c("Transmission From Top 20% of Cases" = 1,
                                   "Epidemic Fadeout" = 2)) + 
  facet_grid(covName~R0, labeller=label_parsed) + 
  xlab("Trait Correlation") + 
  ylab("Proportion") + 
  theme_bw() + 
  theme(legend.position="bottom", legend.title=element_blank(), legend.box="vertical") + 
  guides(color = guide_legend(nrow=2), 
         linetype = guide_legend(nrow=1))
dev.off()

png(file="Fig4_superspreading_CV=5.png", height=6, width=5.5, units='in', res=450)
combo %>% 
  filter(., CV=="CV==5") %>%
  ggplot(., aes(x=corr, y=medianTop20, color=varName)) +
  geom_line(aes(linetype="Transmission From Top 20% of Cases")) + 
  geom_ribbon(aes(ymin=minTop20, ymax=maxTop20, fill=varName), alpha=0.2, linetype='blank') + 
  geom_line(aes(x=corr, y=fadeoutProb, color=varName, linetype="Epidemic Fadeout")) + 
  scale_linetype_manual(values = c("Transmission From Top 20% of Cases" = 1,
                                   "Epidemic Fadeout" = 2)) + 
  facet_grid(covName~R0, labeller=label_parsed) + 
  xlab("Trait Correlation") + 
  ylab("Proportion") + 
  theme_bw() + 
  theme(legend.position="bottom", legend.title=element_blank(), legend.box="vertical") + 
  guides(color = guide_legend(nrow=2), 
         linetype = guide_legend(nrow=1))
dev.off()

png(file="Fig5_superspreading_metrics.png", height=6, width=5.5, units='in', res=450)
combo %>% 
  ggplot(., aes(x=medianDisp, y=medianTop20, color=corr, shape=CV)) + 
  geom_point(size=2) + 
  scale_color_gradient2(low="#2166AC", mid="gray75", high="#B2182B", midpoint=0) + 
  facet_grid(covName~.) + 
  xlab("Dispersion parameter") + 
  ylab("Proportion of transmission from\ntop 20% most infectious cases") + 
  theme_bw() + 
  theme(legend.title=element_blank()) 
dev.off()  


