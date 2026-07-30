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

filter(out2, CV=="CV==1", totalI > 500, cov%in%c("a-g", "c-s", "s-a")) %>% 
  ggplot(., aes(x=disp, y=top20, color=corr)) + 
  geom_point() + 
  facet_grid(cov~R0, labeller=label_parsed) + 
  xlab("Dispersion parameter (k)") + 
  ylab("Proportion of transmission due to\nmost infectious 20% of cases") + 
  theme_bw()

filter(out2, R0=="R[0]==4", totalI > 500, cov%in%c("a-g", "c-s", "s-a")) %>% 
  ggplot(., aes(x=disp, y=top20, color=corr)) + 
  geom_point() + 
  facet_grid(cov~CV, labeller=label_parsed) + 
  xlab("Dispersion parameter (k)") + 
  ylab("Proportion of transmission due to\nmost infectious 20% of cases") + 
  theme_bw()

