## shed gamma
library(tidyverse)

for(corr in c("poscorr","nocorr","negcorr")) {
  for (var in c("hi","med","low")) {
    data <- readRDS(paste0("out_",corr,"_",var,"var_c_alpha.RDS"))
    lapply(data, function(d)
      data.frame(frac0Reff= sum(d[[2]]$numInf==0)/nrow(d[[2]]),
                 Corr = corr,
                 Var = var)) -> sumdata
    }
}


sumdata_use<-do.call("rbind.data.frame", sumdata)

for (var in c("low","med","hi")) {
  for (corr in c("no","neg","pos")) {
    data <- readRDS(paste0("out_",corr,"corr_",var,"var_c_shed.RDS"))
    lapply(data, function(d)
      data.frame(frac0Reff= sum(d[[2]]$numInf==0)/nrow(d[[2]]),
                 Corr = as.factor(corr),
                 Var = as.factor(var)))  -> sumdata1
  }
}

sumdata_use1 <-do.call("rbind.data.frame", sumdata1)

shed_gamma_frac0<-rbind(sumdata_use,sumdata_use1)

sumdata_use %>%
  group_by(Corr,Var)

sumdata$Corr <-factor(sumdata)


ggplot(data = sumdata, mapping = aes(x=Var, y=frac0Reff))+
  geom_point()


## no covariation
# hi var
shed_gamma_nocorr_hi<-vector(mode='list', 50)
for (i in 1:length(out_nocorr_hivar_shed_gamma)) {
  shed_gamma_nocorr_hi[[i]]<-out_nocorr_hivar_shed_gamma[[i]][[2]]$numInf
}
all_inf<-lapply(shed_gamma_nocorr_hi, function(x) x[x != 0 & !is.na(x)])

length(unlist(shed_gamma_nocorr_hi))
length(unlist(all_inf))

shed_gamma_secinf_nocorr_hi <- 1-(length(unlist(all_inf))/length(unlist(shed_gamma_nocorr_hi)))

mean(unlist(all_inf))
median(unlist(all_inf))

# med var
shed_gamma_nocorr_med<-vector(mode='list', 50)
for (i in 1:length(out_nocorr_medvar_shed_gamma)) {
  shed_gamma_nocorr_med[[i]]<-out_nocorr_medvar_shed_gamma[[i]][[2]]$numInf
}
all_inf<-lapply(shed_gamma_nocorr_med, function(x) x[x != 0 & !is.na(x)])

length(unlist(shed_gamma_nocorr_med))
length(unlist(all_inf))
1-(length(unlist(all_inf))/length(unlist(shed_gamma_nocorr_med)))

mean(unlist(all_inf))
median(unlist(all_inf))

# low var
shed_gamma_nocorr_low<-vector(mode='list', 50)
for (i in 1:length(out_nocorr_lowvar_shed_gamma)) {
  shed_gamma_nocorr_low[[i]]<-out_nocorr_lowvar_shed_gamma[[i]][[2]]$numInf
}
all_inf<-lapply(shed_gamma_nocorr_low, function(x) x[x != 0 & !is.na(x)])

length(unlist(shed_gamma_nocorr_low))
length(unlist(all_inf))
1-(length(unlist(all_inf))/length(unlist(shed_gamma_nocorr_low)))

mean(unlist(all_inf))
median(unlist(all_inf))

## positive covariation
# hi var
shed_gamma_poscorr_hi<-vector(mode='list', 50)
for (i in 1:length(out_poscorr_hivar_shed_gamma)) {
  shed_gamma_poscorr_hi[[i]]<-out_poscorr_hivar_shed_gamma[[i]][[2]]$numInf
}
all_inf<-lapply(shed_gamma_poscorr_hi, function(x) x[x != 0 & !is.na(x)])

length(unlist(shed_gamma_poscorr_hi))
length(unlist(all_inf))
1-(length(unlist(all_inf))/length(unlist(shed_gamma_poscorr_hi)))

mean(unlist(all_inf))
median(unlist(all_inf))

# med var
shed_gamma_poscorr_med<-vector(mode='list', 50)
for (i in 1:length(out_poscorr_medvar_shed_gamma)) {
  shed_gamma_poscorr_med[[i]]<-out_poscorr_medvar_shed_gamma[[i]][[2]]$numInf
}
all_inf<-lapply(shed_gamma_poscorr_med, function(x) x[x != 0 & !is.na(x)])

length(unlist(shed_gamma_poscorr_med))
length(unlist(all_inf))

1-(length(unlist(all_inf))/length(unlist(shed_gamma_poscorr_med)))

mean(unlist(all_inf))
median(unlist(all_inf))

# low var
shed_gamma_poscorr_low<-vector(mode='list', 50)
for (i in 1:length(out_poscorr_lowvar_shed_gamma)) {
  shed_gamma_poscorr_low[[i]]<-out_poscorr_lowvar_shed_gamma[[i]][[2]]$numInf
}
all_inf<-lapply(shed_gamma_poscorr_low, function(x) x[x != 0 & !is.na(x)])

length(unlist(shed_gamma_poscorr_low))
length(unlist(all_inf))

1-(length(unlist(all_inf))/length(unlist(shed_gamma_poscorr_low)))

mean(unlist(all_inf))
median(unlist(all_inf))

## neg cov
# hi var
shed_gamma_negcorr_hi<-vector(mode='list', 50)
for (i in 1:length(out_negcorr_hivar_shed_gamma)) {
  shed_gamma_negcorr_hi[[i]]<-out_negcorr_hivar_shed_gamma[[i]][[2]]$numInf
}
all_inf<-lapply(shed_gamma_negcorr_hi, function(x) x[x != 0 & !is.na(x)])

length(unlist(shed_gamma_negcorr_hi))
length(unlist(all_inf))

1-(length(unlist(all_inf))/length(unlist(shed_gamma_negcorr_hi)))

mean(unlist(all_inf))
median(unlist(all_inf))

# med var
shed_gamma_negcorr_med<-vector(mode='list', 50)
for (i in 1:length(out_negcorr_medvar_shed_gamma)) {
  shed_gamma_negcorr_med[[i]]<-out_negcorr_medvar_shed_gamma[[i]][[2]]$numInf
}
all_inf<-lapply(shed_gamma_negcorr_med, function(x) x[x != 0 & !is.na(x)])

length(unlist(shed_gamma_negcorr_med))
length(unlist(all_inf))

1-(length(unlist(all_inf))/length(unlist(shed_gamma_negcorr_med)))

mean(unlist(all_inf))
median(unlist(all_inf))

# low var
shed_gamma_negcorr_low<-vector(mode='list', 50)
for (i in 1:length(out_negcorr_lowvar_shed_gamma)) {
  shed_gamma_negcorr_low[[i]]<-out_negcorr_lowvar_shed_gamma[[i]][[2]]$numInf
}
all_inf<-lapply(shed_gamma_negcorr_low, function(x) x[x != 0 & !is.na(x)])

length(unlist(shed_gamma_negcorr_low))
length(unlist(all_inf))

1-(length(unlist(all_inf))/length(unlist(shed_gamma_negcorr_low)))

mean(unlist(all_inf))
median(unlist(all_inf))

