var1=c('c','c','c','s','s','a')
var2=c('s','a','g','a','g','g')

for (R in c(1,2,4,8)) {
  for (S in c(100, 250, 500, 1000)) { ## vary iniital number of susceptibles
    for (CV in c(.2, 1, 5)) {
      for (parcomb in seq(1,6)) {
        ## Pull the results across all correlations
        lapply(seq(-0.9, 0.9, 0.1), function(corr) readRDS(file=paste0("out_R0=",R,"_S=",S,"_CV=",CV,"_corr=",corr,"_cov=",var1[parcomb],"-",var2[parcomb],".RDS"))) %>%
          do.call("rbind.data.frame",.) %>% 
          mutate(trait=paste(gsub("-",",",cov), "covary")) -> outCov
        ## Pull the results for single variation with the same CV
        readRDS(file=paste0("out_R0=",R,"_S=",S,"_CV=",CV,"_var=",var1[parcomb],".RDS")) %>% 
          replicate(19, ., simplify=FALSE) %>%
          do.call("rbind",.) %>%
          mutate(corr=rep(seq(-0.9,0.9,0.1), each=1000),
                 trait=paste(var, "varies")) -> outVar1
        readRDS(file=paste0("out_R0=",R,"_S=",S,"_CV=",CV,"_var=",var2[parcomb],".RDS")) %>% 
          replicate(19, ., simplify=FALSE) %>%
          do.call("rbind",.) %>%
          mutate(corr=rep(seq(-0.9,0.9,0.1), each=1000),
                 trait=paste(var, "varies")) -> outVar2
        ## Pull the results for no variation
        readRDS(file=paste0("out_R0=",R,"_S=",S,"_no_variation.RDS")) %>% 
          replicate(19, ., simplify=FALSE) %>%
          do.call("rbind",.) %>%
          mutate(corr=rep(seq(-0.9,0.9,0.1), each=1000),
                 trait="no variation") -> outNoVar
        
        bind_rows(outCov[,colnames(outNoVar)], 
                  outVar1[,colnames(outNoVar)], 
                  outVar2[,colnames(outNoVar)], 
                  outNoVar) -> out
        
        out %>% group_by(R0, S, corr, trait) %>% summarise(medPeak=median(peak), loPeak=quantile(peak,0.025), hiPeak=quantile(peak,0.975)) %>%
          ggplot(., aes(x=corr, y=medPeak, group=trait, color=trait)) + 
          geom_line() + 
          geom_ribbon(aes(ymin=loPeak, ymax=hiPeak), alpha=0.3)
        
        
        
        ggplot(out, aes(x=corr, y=peak, group=trait, colour=trait)) + 
          geom_ribbon()
        
        
        
        
        
          