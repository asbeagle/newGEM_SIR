### contact-alpha-poscor-hivar
## all trait values
# shed values for each simulation
store_alpha_hi<-vector(mode='list', 50)
for (i in 1:length(out_poscorr_hivar_c_alpha)) {
  store_alpha_hi[[i]]<-out_poscorr_hivar_c_alpha[[i]][[2]]$alpha
}

# contact rates for each simulation
store_c_hi<-vector(mode='list', 50)
for (i in 1:length(out_poscorr_hivar_c_alpha)) {
  store_c_hi[[i]]<-out_poscorr_hivar_c_alpha[[i]][[2]]$c
}

hist(unlist(store_alpha_hi))
hist(unlist(store_c_hi))

# epidemic duration of all simulations
store_epidemic_duration_hi <- vector(mode='list', 50)
for (i in 1:length(out_poscorr_hivar_c_alpha)) {
  store_epidemic_duration_hi[[i]]<-max(out_poscorr_hivar_c_alpha[[i]][[1]]$t)
}

# secondary infection in each simulation
store_numInf_hi<-vector(mode='list', 50)
for (i in 1:length(out_poscorr_hivar_c_alpha)) {
  store_numInf_hi[[i]]<-out_poscorr_hivar_c_alpha[[i]][[2]]$numInf
}

# max secondary infection in each simulation
store_maxnumInf_hi<-vector(mode='list', 50)
for (i in 1:length(out_poscorr_hivar_c_alpha)) {
  store_maxnumInf_hi[[i]]<-max(out_poscorr_hivar_c_alpha[[i]][[2]]$numInf,na.rm=T)
}

### contact-alpha-poscor-lowvar
## all trait values
# shed values for each simulation
store_alpha_low<-vector(mode='list', 50)
for (i in 1:length(out_poscorr_lowvar_c_alpha)) {
  store_alpha_low[[i]]<-out_poscorr_lowvar_c_alpha[[i]][[2]]$alpha
}

# contact rates for each simulation
store_c_low<-vector(mode='list', 50)
for (i in 1:length(out_poscorr_lowvar_c_alpha)) {
  store_c_low[[i]]<-out_poscorr_lowvar_c_alpha[[i]][[2]]$c
}

hist(unlist(store_alpha_low))
hist(unlist(store_c_low))

store_epidemic_duration_low <- vector(mode='list', 50)
for (i in 1:length(out_poscorr_lowvar_c_alpha)) {
  store_epidemic_duration_low[[i]]<-max(out_poscorr_lowvar_c_alpha[[i]][[1]]$t)
}

# secondary infection in each simulation
store_numInf_low<-vector(mode='list', 50)
for (i in 1:length(out_poscorr_lowvar_c_alpha)) {
  store_numInf_low[[i]]<-out_poscorr_lowvar_c_alpha[[i]][[2]]$numInf
}

# max secondary infection in each simulation
store_maxnumInf_low<-vector(mode='list', 50)
for (i in 1:length(out_poscorr_lowvar_c_alpha)) {
  store_maxnumInf_low[[i]]<-max(out_poscorr_lowvar_c_alpha[[i]][[2]]$numInf,na.rm=T)
}

### contact-alpha-poscor-medvar
## all trait values
# shed values for each simulation
store_alpha_med<-vector(mode='list', 50)
for (i in 1:length(out_poscorr_medvar_c_alpha)) {
  store_alpha_med[[i]]<-out_poscorr_medvar_c_alpha[[i]][[2]]$alpha
}

# contact rates for each simulation
store_c_med<-vector(mode='list', 50)
for (i in 1:length(out_poscorr_medvar_c_alpha)) {
  store_c_med[[i]]<-out_poscorr_medvar_c_alpha[[i]][[2]]$c
}

par(mfrow=c(1,3))

hist(unlist(store_alpha_med))
hist(unlist(store_c_med))

store_epidemic_duration_med <- vector(mode='list', 50)
for (i in 1:length(out_poscorr_medvar_c_alpha)) {
  store_epidemic_duration_med[[i]]<-max(out_poscorr_medvar_c_alpha[[i]][[1]]$t)
}

# secondary infection in each simulation
store_numInf_med<-vector(mode='list', 50)
for (i in 1:length(out_poscorr_medvar_c_alpha)) {
  store_numInf_med[[i]]<-out_poscorr_medvar_c_alpha[[i]][[2]]$numInf
}

# max secondary infection in each simulation
store_maxnumInf_med<-vector(mode='list', 50)
for (i in 1:length(out_poscorr_medvar_c_alpha)) {
  store_maxnumInf_med[[i]]<-max(out_poscorr_medvar_c_alpha[[i]][[2]]$numInf,na.rm=T)
}

######### cov-alpha-contact-hivar pos corr
# maximum shed values for each simulation
store_maxalpha_hi<-vector(mode='list', 50)
for (i in 1:length(out_poscorr_hivar_c_alpha)) {
  store_maxalpha_hi[[i]]<-max(out_poscorr_hivar_c_alpha[[i]][[2]]$alpha,na.rm=T)
}

# maximum contact rates for each simulation
store_maxc_hi<-vector(mode='list', 50)
for (i in 1:length(out_poscorr_hivar_c_alpha)) {
  store_maxc_hi[[i]]<-max(out_poscorr_hivar_c_alpha[[i]][[2]]$c,na.rm=T)
}

par(mfrow=c(1,2))
hist(unlist(store_maxalpha_hi))
hist(unlist(store_maxc_hi))

######### cov-alpha-contact-medvar pos corr
# maximum shed values for each simulation
store_maxalpha_med<-vector(mode='list', 50)
for (i in 1:length(out_poscorr_medvar_c_alpha)) {
  store_maxalpha_med[[i]]<-max(out_poscorr_medvar_c_alpha[[i]][[2]]$alpha,na.rm=T)
}

# maximum contact rates for each simulation
store_maxc_med<-vector(mode='list', 50)
for (i in 1:length(out_poscorr_medvar_c_alpha)) {
  store_maxc_med[[i]]<-max(out_poscorr_medvar_c_alpha[[i]][[2]]$c,na.rm=T)
}

par(mfrow=c(1,2))
hist(unlist(store_maxalpha_med))
hist(unlist(store_maxc_med))

######### cov-alpha-contact-lowvar pos corr
# maximum shed values for each simulation
store_maxalpha_low<-vector(mode='list', 50)
for (i in 1:length(out_poscorr_lowvar_c_alpha)) {
  store_maxalpha_low[[i]]<-max(out_poscorr_lowvar_c_alpha[[i]][[2]]$alpha,na.rm=T)
}

# maximum contact rates for each simulation
store_maxc_low<-vector(mode='list', 50)
for (i in 1:length(out_poscorr_lowvar_c_alpha)) {
  store_maxc_low[[i]]<-max(out_poscorr_lowvar_c_alpha[[i]][[2]]$c,na.rm=T)
}

par(mfrow=c(3,2))
hist(unlist(store_maxalpha_low))
hist(unlist(store_maxc_low))
