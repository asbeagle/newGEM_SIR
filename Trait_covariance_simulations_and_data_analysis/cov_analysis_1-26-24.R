noVar = readRDS(file="Stochastic_no_variation_simulations.RDS")
cVar = readRDS(file="Variance_in_contact_results_12-28-23.RDS")
sVar = readRDS(file="Variance_in_shedding_results_12-28-23.RDS")
aVar = readRDS(file="Variance_in_virulence_results_12-28-23.RDS")
gVar = readRDS(file="Variance_in_recovery_results_12-28-23.RDS")
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

csCov$corr = round(csCov$corr,1)

cVar$R0 = (strsplit(cVar$R0,"=") %>% lapply(., function(a) as.numeric(a[2])) %>% unlist)
cVar$CV = (strsplit(cVar$CV,"=") %>% lapply(., function(a) as.numeric(a[2])) %>% unlist)
cVarPeaks = (cVar %>% group_by(R0,CV) %>% summarize(meanPeak=mean(peak)))

sVar$R0 = (strsplit(sVar$R0,"=") %>% lapply(., function(a) as.numeric(a[2])) %>% unlist)
sVar$CV = (strsplit(sVar$CV,"=") %>% lapply(., function(a) as.numeric(a[2])) %>% unlist)
sVarPeaks = (sVar %>% group_by(R0,CV) %>% summarize(meanPeak=mean(peak)))

ggplot(csCov, aes(x=corr, y=peak)) + 
  geom_point() + 
  geom_hline(data=cVarPeaks, mapping=aes(yintercept=meanPeak, color="red")) +
  geom_hline(data=sVarPeaks, mapping=aes(yintercept=meanPeak, color="blue")) +
  facet_grid(R0~CV) +
  theme(legend.position="none")
  
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




