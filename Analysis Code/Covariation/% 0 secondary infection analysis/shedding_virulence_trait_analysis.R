#### shedding and virulence

## all trait values in poss cor
## hi var
# alpha values for each simulation
store_alpha_hi<-vector(mode='list', 50)
for (i in 1:length(out_poscorr_hivar_shed_alpha)) {
  store_alpha_hi[[i]]<-out_poscorr_hivar_shed_alpha[[i]][[2]]$alpha
}

# contact rates for each simulation
store_shed_hi<-vector(mode='list', 50)
for (i in 1:length(out_poscorr_hivar_shed_alpha)) {
  store_shed_hi[[i]]<-out_poscorr_hivar_shed_alpha[[i]][[2]]$shed
}

## med var
# alpha values for each simulation
store_alpha_med<-vector(mode='list', 50)
for (i in 1:length(out_poscorr_medvar_shed_alpha)) {
  store_alpha_med[[i]]<-out_poscorr_medvar_shed_alpha[[i]][[2]]$alpha
}

# contact rates for each simulation
store_shed_med<-vector(mode='list', 50)
for (i in 1:length(out_poscorr_medvar_shed_alpha)) {
  store_shed_med[[i]]<-out_poscorr_medvar_shed_alpha[[i]][[2]]$shed
}

## low var
# alpha values for each simulation
store_alpha_low<-vector(mode='list', 50)
for (i in 1:length(out_poscorr_lowvar_shed_alpha)) {
  store_alpha_low[[i]]<-out_poscorr_lowvar_shed_alpha[[i]][[2]]$alpha
}

# contact rates for each simulation
store_shed_low<-vector(mode='list', 50)
for (i in 1:length(out_poscorr_lowvar_shed_alpha)) {
  store_shed_low[[i]]<-out_poscorr_lowvar_shed_alpha[[i]][[2]]$shed
}

par(mfrow=c(3,2))
hist(unlist(store_alpha_low))
hist(unlist(store_shed_low))
hist(unlist(store_alpha_med))
hist(unlist(store_shed_med))
hist(unlist(store_alpha_hi))
hist(unlist(store_shed_hi))

## all trait values in no cor
## hi var
# alpha values for each simulation
store_alpha_hi<-vector(mode='list', 50)
for (i in 1:length(out_nocorr_hivar_shed_alpha)) {
  store_alpha_hi[[i]]<-out_nocorr_hivar_shed_alpha[[i]][[2]]$alpha
}

# contact rates for each simulation
store_shed_hi<-vector(mode='list', 50)
for (i in 1:length(out_nocorr_hivar_shed_alpha)) {
  store_shed_hi[[i]]<-out_nocorr_hivar_shed_alpha[[i]][[2]]$shed
}

## med var
# alpha values for each simulation
store_alpha_med<-vector(mode='list', 50)
for (i in 1:length(out_nocorr_medvar_shed_alpha)) {
  store_alpha_med[[i]]<-out_nocorr_medvar_shed_alpha[[i]][[2]]$alpha
}

# contact rates for each simulation
store_shed_med<-vector(mode='list', 50)
for (i in 1:length(out_nocorr_medvar_shed_alpha)) {
  store_shed_med[[i]]<-out_nocorr_medvar_shed_alpha[[i]][[2]]$shed
}

## low var
# alpha values for each simulation
store_alpha_low<-vector(mode='list', 50)
for (i in 1:length(out_nocorr_lowvar_shed_alpha)) {
  store_alpha_low[[i]]<-out_nocorr_lowvar_shed_alpha[[i]][[2]]$alpha
}

# contact rates for each simulation
store_shed_low<-vector(mode='list', 50)
for (i in 1:length(out_nocorr_lowvar_shed_alpha)) {
  store_shed_low[[i]]<-out_nocorr_lowvar_shed_alpha[[i]][[2]]$shed
}

par(mfrow=c(3,2))
hist(unlist(store_alpha_low))
hist(unlist(store_shed_low))
hist(unlist(store_alpha_med))
hist(unlist(store_shed_med))
hist(unlist(store_alpha_hi))
hist(unlist(store_shed_hi))

## max trait values in poss cor
## hi var
# alpha values for each simulation
store_maxalpha_hi<-vector(mode='list', 50)
for (i in 1:length(out_poscorr_hivar_shed_alpha)) {
  store_maxalpha_hi[[i]]<-max(out_poscorr_hivar_shed_alpha[[i]][[2]]$alpha, na.rm=T)
}

# contact rates for each simulation
store_maxshed_hi<-vector(mode='list', 50)
for (i in 1:length(out_poscorr_hivar_shed_alpha)) {
  store_maxshed_hi[[i]]<-max(out_poscorr_hivar_shed_alpha[[i]][[2]]$shed,na.rm=T)
}

## med var
# alpha values for each simulation
store_maxalpha_med<-vector(mode='list', 50)
for (i in 1:length(out_poscorr_medvar_shed_alpha)) {
  store_maxalpha_med[[i]]<-max(out_poscorr_medvar_shed_alpha[[i]][[2]]$alpha,na.rm=T)
}

# contact rates for each simulation
store_maxshed_med<-vector(mode='list', 50)
for (i in 1:length(out_poscorr_medvar_shed_alpha)) {
  store_maxshed_med[[i]]<-max(out_poscorr_medvar_shed_alpha[[i]][[2]]$shed,na.rm=T)
}

## low var
# alpha values for each simulation
store_maxalpha_low<-vector(mode='list', 50)
for (i in 1:length(out_poscorr_lowvar_shed_alpha)) {
  store_maxalpha_low[[i]]<-max(out_poscorr_lowvar_shed_alpha[[i]][[2]]$alpha,na.rm=T)
}

# contact rates for each simulation
store_maxshed_low<-vector(mode='list', 50)
for (i in 1:length(out_poscorr_lowvar_shed_alpha)) {
  store_maxshed_low[[i]]<-max(out_poscorr_lowvar_shed_alpha[[i]][[2]]$shed,na.rm=T)
}

par(mfrow=c(3,2))
hist(unlist(store_maxalpha_low))
hist(unlist(store_maxshed_low))
hist(unlist(store_maxalpha_med))
hist(unlist(store_maxshed_med))
hist(unlist(store_maxalpha_hi))
hist(unlist(store_maxshed_hi))


### distribution of traits
d=data.frame(picks_...)
mutate(d, r0=(shed/1+shed))
       
data.frame(pick_indiv_multi(...))%>% 
  mutate(., R0 = ,,,) %>%
  with(., hist(R0))



## all trait values in no cor
## hi var
# alpha values for each simulation
store_maxalpha_hi<-vector(mode='list', 50)
for (i in 1:length(out_nocorr_hivar_shed_alpha)) {
  store_maxalpha_hi[[i]]<-max(out_nocorr_hivar_shed_alpha[[i]][[2]]$alpha, na.rm=T)
}

# contact rates for each simulation
store_maxshed_hi<-vector(mode='list', 50)
for (i in 1:length(out_nocorr_hivar_shed_alpha)) {
  store_maxshed_hi[[i]]<-max(out_nocorr_hivar_shed_alpha[[i]][[2]]$shed,na.rm=T)
}

## med var
# alpha values for each simulation
store_maxalpha_med<-vector(mode='list', 50)
for (i in 1:length(out_nocorr_medvar_shed_alpha)) {
  store_maxalpha_med[[i]]<-max(out_nocorr_medvar_shed_alpha[[i]][[2]]$alpha,na.rm=T)
}

# contact rates for each simulation
store_maxshed_med<-vector(mode='list', 50)
for (i in 1:length(out_nocorr_medvar_shed_alpha)) {
  store_maxshed_med[[i]]<-max(out_nocorr_medvar_shed_alpha[[i]][[2]]$shed,na.rm=T)
}

## low var
# alpha values for each simulation
store_maxalpha_low<-vector(mode='list', 50)
for (i in 1:length(out_nocorr_lowvar_shed_alpha)) {
  store_maxalpha_low[[i]]<-max(out_nocorr_lowvar_shed_alpha[[i]][[2]]$alpha,na.rm=T)
}

# contact rates for each simulation
store_maxshed_low<-vector(mode='list', 50)
for (i in 1:length(out_nocorr_lowvar_shed_alpha)) {
  store_maxshed_low[[i]]<-max(out_nocorr_lowvar_shed_alpha[[i]][[2]]$shed,na.rm=T)
}

par(mfrow=c(1,2))
hist(unlist(store_maxalpha_low))
hist(unlist(store_maxshed_low))
hist(unlist(store_maxalpha_med))
hist(unlist(store_maxshed_med))
hist(unlist(store_maxalpha_hi))
hist(unlist(store_maxshed_hi))


### max secondary infections
# pos cor hi var
store_maxnumInf_hi<-vector(mode='list', 50)
for (i in 1:length(out_poscorr_hivar_c_alpha)) {
  store_maxnumInf_hi[[i]]<-max(out_poscorr_hivar_shed_alpha[[i]][[2]]$numInf,na.rm=T)
}

store_maxnumInf_med<-vector(mode='list', 50)
for (i in 1:length(out_poscorr_hivar_c_alpha)) {
  store_maxnumInf_med[[i]]<-max(out_poscorr_medvar_shed_alpha[[i]][[2]]$numInf,na.rm=T)
}

store_maxnumInf_low<-vector(mode='list', 50)
for (i in 1:length(out_poscorr_hivar_c_alpha)) {
  store_maxnumInf_low[[i]]<-max(out_poscorr_lowvar_shed_alpha[[i]][[2]]$numInf,na.rm=T)
}

par(mfrow=c(1,3))
hist(unlist(store_maxnumInf_low))
hist(unlist(store_maxnumInf_med))
hist(unlist(store_maxnumInf_hi))

# no cor
store_maxnumInf_hi<-vector(mode='list', 50)
for (i in 1:length(out_nocorr_hivar_c_alpha)) {
  store_maxnumInf_hi[[i]]<-max(out_nocorr_hivar_shed_alpha[[i]][[2]]$numInf,na.rm=T)
}

store_maxnumInf_med<-vector(mode='list', 50)
for (i in 1:length(out_nocorr_medvar_c_alpha)) {
  store_maxnumInf_med[[i]]<-max(out_nocorr_medvar_shed_alpha[[i]][[2]]$numInf,na.rm=T)
}

store_maxnumInf_low<-vector(mode='list', 50)
for (i in 1:length(out_nocorr_lowvar_c_alpha)) {
  store_maxnumInf_low[[i]]<-max(out_nocorr_lowvar_shed_alpha[[i]][[2]]$numInf,na.rm=T)
}

hist(unlist(store_maxnumInf_low))
hist(unlist(store_maxnumInf_med))
hist(unlist(store_maxnumInf_hi))

