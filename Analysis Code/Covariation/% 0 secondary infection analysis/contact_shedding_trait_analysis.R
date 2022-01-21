### Covariation between contact rate and shedding rate

## high var neg corr contact shed
# get all contact rates for all simulations
store_c<-vector(mode='list', 50)
for (i in 1:length(out_negcorr_hivar_c_shed)) {
  store_c[[i]]<-out_negcorr_hivar_c_shed[[i]][[2]]$c
}

# get all shedding rates for all simulations
store_shed<-vector(mode='list', 50)
for (i in 1:length(out_negcorr_hivar_c_shed)) {
  store_shed[[i]]<-out_negcorr_hivar_c_shed[[i]][[2]]$shed
}

# epidemic duration of all simulations
store_epidemic_duration <- vector(mode='list', 50)
for (i in 1:length(out_negcorr_hivar_c_shed)) {
  store_epidemic_duration[[i]]<-max(out_negcorr_hivar_c_shed[[i]][[1]]$t)
}

# maximum shed values for each simulation
store_maxshed<-vector(mode='list', 50)
for (i in 1:length(out_negcorr_hivar_c_shed)) {
  store_maxshed[[i]]<-max(out_negcorr_hivar_c_shed[[i]][[2]]$shed,na.rm=T)
}

# maximum contact rates for each simulation
store_maxc<-vector(mode='list', 50)
for (i in 1:length(out_negcorr_hivar_c_shed)) {
  store_maxc[[i]]<-max(out_negcorr_hivar_c_shed[[i]][[2]]$c,na.rm=T)
}

store_c<-vector(mode='list', 50)
for (i in 1:length(out_negcorr_hivar_c_shed)) {
  store_c[[i]]<-out_negcorr_hivar_c_shed[[i]][[2]]$c
}

store_shed<-vector(mode='list', 50)
for (i in 1:length(out_negcorr_hivar_c_shed)) {
  store_shed[[i]]<-out_negcorr_hivar_c_shed[[i]][[2]]$c
}

# secondary infection in each simulation
store_numInf<-vector(mode='list', 50)
for (i in 1:length(out_negcorr_hivar_c_shed)) {
  store_numInf[[i]]<-out_negcorr_hivar_c_shed[[i]][[2]]$numInf
}

par(mfrow=c(1,2))
plot(unlist(store_epidemic_duration), unlist(store_maxshed), ylab="max shed", xlab="epidemic duration", main="contact-shed-hivar-negcorr")
plot(unlist(store_epidemic_duration), unlist(store_maxc), ylab="max c", xlab="epidemic duration", main="contact-shed-hivar-negcorr")

cor.test(unlist(store_epidemic_duration),unlist(store_maxshed))
cor.test(unlist(store_epidemic_duration),unlist(store_maxc))

hist(unlist(store_maxshed))
hist(unlist(store_maxc), xlim=c(0,15))

# plot secondary infection against beta
plot(unlist(store_c)*unlist(store_shed)/(1+unlist(store_shed)), unlist(store_numInf), xlim=c(0,0.1))

# linear model 
# max shed value against epidemic duration
lm(unlist(store_epidemic_duration)~unlist(store_maxshed))

# max contact values against epidemic duration
lm(unlist(store_epidemic_duration)~unlist(store_maxc))

#### contact-shedding-posscorr-hivar

# epidemic duration of all simulations
store_epidemic_duration_poscor <- vector(mode='list', 50)
for (i in 1:length(out_poscorr_hivar_c_shed)) {
  store_epidemic_duration_poscor[[i]]<-max(out_poscorr_hivar_c_shed[[i]][[1]]$t)
}

# maximum shed values for each simulation
store_maxshed_poscor<-vector(mode='list', 50)
for (i in 1:length(out_poscorr_hivar_c_shed)) {
  store_maxshed_poscor[[i]]<-max(out_poscorr_hivar_c_shed[[i]][[2]]$shed,na.rm=T)
}

# maximum contact rates for each simulation
store_maxc_poscor<-vector(mode='list', 50)
for (i in 1:length(out_poscorr_hivar_c_shed)) {
  store_maxc_poscor[[i]]<-max(out_poscorr_hivar_c_shed[[i]][[2]]$c,na.rm=T)
}

plot(unlist(store_epidemic_duration_poscor),unlist(store_maxshed_poscor),ylab="max shed", xlab="epidemic duration", main="contact-shed-hivar-poscorr")
plot(unlist(store_epidemic_duration_poscor),unlist(store_maxc_poscor),ylab="max c", xlab="epidemic duration", main="contact-shed-hivar-poscorr")

hist(unlist(store_maxc_poscor))
hist(unlist(store_maxshed_poscor))

# linear model 
# max shed value against epidemic duration
lm(unlist(store_epidemic_duration_poscor)~unlist(store_maxshed_poscor))

# max contact values against epidemic duration
lm(unlist(store_epidemic_duration_poscor)~unlist(store_maxc_poscor))

cor.test(unlist(store_epidemic_duration_poscor),unlist(store_maxc_poscor))
cor.test(unlist(store_epidemic_duration_poscor),unlist(store_maxshed_poscor))

#### contact-shedding-nocorr-hivar

# epidemic duration of all simulations
store_epidemic_duration_nocor <- vector(mode='list', 50)
for (i in 1:length(out_nocorr_hivar_c_shed)) {
  store_epidemic_duration_nocor[[i]]<-max(out_nocorr_hivar_c_shed[[i]][[1]]$t)
}

# maximum shed values for each simulation
store_maxshed_nocor<-vector(mode='list', 50)
for (i in 1:length(out_nocorr_hivar_c_shed)) {
  store_maxshed_nocor[[i]]<-max(out_nocorr_hivar_c_shed[[i]][[2]]$shed,na.rm=T)
}

# maximum contact rates for each simulation
store_maxc_nocor<-vector(mode='list', 50)
for (i in 1:length(out_nocorr_hivar_c_shed)) {
  store_maxc_nocor[[i]]<-max(out_nocorr_hivar_c_shed[[i]][[2]]$c,na.rm=T)
}

plot(unlist(store_epidemic_duration_nocor),unlist(store_maxshed_nocor),ylab="max shed", xlab="epidemic duration", main="contact-shed-hivar-nocorr")
plot(unlist(store_epidemic_duration_nocor),unlist(store_maxc_nocor),ylab="max c", xlab="epidemic duration", main="contact-shed-hivar-nocorr")

hist(unlist(store_maxshed_nocor))
hist(unlist(store_maxc_nocor))

# linear model 
# max shed value against epidemic duration
lm(unlist(store_epidemic_duration_nocor)~unlist(store_maxshed_nocor))

# max contact values against epidemic duration
lm(unlist(store_epidemic_duration_nocor)~unlist(store_maxc_nocor))
cor.test(unlist(store_epidemic_duration_nocor), unlist(store_maxc_nocor))
cor.test(unlist(store_epidemic_duration_nocor), unlist(store_maxshed_nocor))

### pos cor has the least amount of extinctions
### no cor has second least
### neg cor none make it to tmax

#### contact-virulence
# epidemic duration of all simulations
store_epidemic_duration_poscor_hi <- vector(mode='list', 50)
for (i in 1:length(out_poscorr_hivar_c_alpha)) {
  store_epidemic_duration_poscor_hi[[i]]<-max(out_poscorr_hivar_c_alpha[[i]][[1]]$t)
}

# maximum alpha values for each simulation
store_maxalpha_poscor_hi<-vector(mode='list', 50)
for (i in 1:length(out_poscorr_hivar_c_alpha)) {
  store_maxalpha_poscor_hi[[i]]<-max(out_poscorr_hivar_c_alpha[[i]][[2]]$alpha,na.rm=T)
}

# maximum contact rates for each simulation
store_maxc_poscor_hi<-vector(mode='list', 50)
for (i in 1:length(out_poscorr_hivar_c_alpha)) {
  store_maxc_poscor_hi[[i]]<-max(out_poscorr_hivar_c_alpha[[i]][[2]]$c,na.rm=T)
}

## all alpha values
store_alpha_poscor_hi<-vector(mode='list', 50)
for (i in 1:length(out_poscorr_hivar_c_alpha)) {
  store_alpha_poscor_hi[[i]]<-out_poscorr_hivar_c_alpha[[i]][[2]]$alpha
}

## all contact values
store_c_poscor_hi<-vector(mode='list', 50)
for (i in 1:length(out_poscorr_hivar_c_alpha)) {
  store_c_poscor_hi[[i]]<-out_poscorr_hivar_c_alpha[[i]][[2]]$c
}

plot(unlist(store_epidemic_duration_poscor_hi),unlist(store_maxc_poscor_hi))
plot(unlist(store_epidemic_duration_poscor_hi),unlist(store_maxalpha_poscor_hi))

plot(unlist(store_maxc_poscor_hi),unlist(store_maxalpha_poscor_hi))

hist(unlist(store_alpha_poscor_hi))
hist(unlist(store_c_poscor_hi))

## secondary infection
store_inf_poscor_hi<-vector(mode='list', 50)
for (i in 1:length(out_poscorr_hivar_c_alpha)) {
  store_inf_poscor_hi[[i]]<-out_poscorr_hivar_c_alpha[[i]][[2]]$numInf
}

hist(unlist(store_inf_poscor_hi))

plot((unlist(store_c_poscor_hi))*(.05/(1+.05)), unlist(store_inf_poscor_hi))


######### cov-shed-contact-lowvar neg corr
# epidemic duration of all simulations
store_epidemic_duration_low <- vector(mode='list', 50)
for (i in 1:length(out_negcorr_lowvar_c_shed)) {
  store_epidemic_duration_low[[i]]<-max(out_negcorr_lowvar_c_shed[[i]][[1]]$t)
}

# maximum shed values for each simulation
store_maxshed_low<-vector(mode='list', 50)
for (i in 1:length(out_negcorr_lowvar_c_shed)) {
  store_maxshed_low[[i]]<-max(out_negcorr_lowvar_c_shed[[i]][[2]]$shed,na.rm=T)
}

# maximum contact rates for each simulation
store_maxc_low<-vector(mode='list', 50)
for (i in 1:length(out_negcorr_lowvar_c_shed)) {
  store_maxc_low[[i]]<-max(out_negcorr_lowvar_c_shed[[i]][[2]]$c,na.rm=T)
}

######### cov-shed-contact-medvar neg corr
# epidemic duration of all simulations
store_epidemic_duration_med <- vector(mode='list', 50)
for (i in 1:length(out_negcorr_medvar_c_shed)) {
  store_epidemic_duration_med[[i]]<-max(out_negcorr_medvar_c_shed[[i]][[1]]$t)
}

# maximum shed values for each simulation
store_maxshed_med<-vector(mode='list', 50)
for (i in 1:length(out_negcorr_medvar_c_shed)) {
  store_maxshed_med[[i]]<-max(out_negcorr_medvar_c_shed[[i]][[2]]$shed,na.rm=T)
}

# maximum contact rates for each simulation
store_maxc_med<-vector(mode='list', 50)
for (i in 1:length(out_negcorr_medvar_c_shed)) {
  store_maxc_med[[i]]<-max(out_negcorr_medvar_c_shed[[i]][[2]]$c,na.rm=T)
}

######### cov-shed-contact-lowvar pos corr
# epidemic duration of all simulations
store_epidemic_duration_low <- vector(mode='list', 50)
for (i in 1:length(out_poscorr_lowvar_c_shed)) {
  store_epidemic_duration_low[[i]]<-max(out_poscorr_lowvar_c_shed[[i]][[1]]$t)
}

# maximum shed values for each simulation
store_maxshed_low<-vector(mode='list', 50)
for (i in 1:length(out_poscorr_lowvar_c_shed)) {
  store_maxshed_low[[i]]<-max(out_poscorr_lowvar_c_shed[[i]][[2]]$shed,na.rm=T)
}

# maximum contact rates for each simulation
store_maxc_low<-vector(mode='list', 50)
for (i in 1:length(out_poscorr_lowvar_c_shed)) {
  store_maxc_low[[i]]<-max(out_poscorr_lowvar_c_shed[[i]][[2]]$c,na.rm=T)
}

######### cov-shed-contact-medvar neg corr
# epidemic duration of all simulations
store_epidemic_duration_med <- vector(mode='list', 50)
for (i in 1:length(out_poscorr_medvar_c_shed)) {
  store_epidemic_duration_med[[i]]<-max(out_poscorr_medvar_c_shed[[i]][[1]]$t)
}

# maximum shed values for each simulation
store_maxshed_med<-vector(mode='list', 50)
for (i in 1:length(out_poscorr_medvar_c_shed)) {
  store_maxshed_med[[i]]<-max(out_poscorr_medvar_c_shed[[i]][[2]]$shed,na.rm=T)
}

# maximum contact rates for each simulation
store_maxc_med<-vector(mode='list', 50)
for (i in 1:length(out_poscorr_medvar_c_shed)) {
  store_maxc_med[[i]]<-max(out_poscorr_medvar_c_shed[[i]][[2]]$c,na.rm=T)
}

######### cov-shed-contact-lowvar pos corr
# epidemic duration of all simulations
store_epidemic_duration_low <- vector(mode='list', 50)
for (i in 1:length(out_poscorr_lowvar_c_shed)) {
  store_epidemic_duration_low[[i]]<-max(out_poscorr_lowvar_c_shed[[i]][[1]]$t)
}

# maximum shed values for each simulation
store_maxshed_low<-vector(mode='list', 50)
for (i in 1:length(out_poscorr_lowvar_c_shed)) {
  store_maxshed_low[[i]]<-max(out_poscorr_lowvar_c_shed[[i]][[2]]$shed,na.rm=T)
}

# maximum contact rates for each simulation
store_maxc_low<-vector(mode='list', 50)
for (i in 1:length(out_poscorr_lowvar_c_shed)) {
  store_maxc_low[[i]]<-max(out_poscorr_lowvar_c_shed[[i]][[2]]$c,na.rm=T)
}

######### cov-shed-contact-medvar neg corr
# epidemic duration of all simulations
store_epidemic_duration_med <- vector(mode='list', 50)
for (i in 1:length(out_poscorr_medvar_c_shed)) {
  store_epidemic_duration_med[[i]]<-max(out_poscorr_medvar_c_shed[[i]][[1]]$t)
}

# maximum shed values for each simulation
store_maxshed_med<-vector(mode='list', 50)
for (i in 1:length(out_poscorr_medvar_c_shed)) {
  store_maxshed_med[[i]]<-max(out_poscorr_medvar_c_shed[[i]][[2]]$shed,na.rm=T)
}

# maximum contact rates for each simulation
store_maxc_med<-vector(mode='list', 50)
for (i in 1:length(out_poscorr_medvar_c_shed)) {
  store_maxc_med[[i]]<-max(out_poscorr_medvar_c_shed[[i]][[2]]$c,na.rm=T)
}

######### cov-shed-contact-lowvar no corr
# epidemic duration of all simulations
store_epidemic_duration_low <- vector(mode='list', 50)
for (i in 1:length(out_nocorr_lowvar_c_shed)) {
  store_epidemic_duration_low[[i]]<-max(out_nocorr_lowvar_c_shed[[i]][[1]]$t)
}

# maximum shed values for each simulation
store_maxshed_low<-vector(mode='list', 50)
for (i in 1:length(out_nocorr_lowvar_c_shed)) {
  store_maxshed_low[[i]]<-max(out_nocorr_lowvar_c_shed[[i]][[2]]$shed,na.rm=T)
}

# maximum contact rates for each simulation
store_maxc_low<-vector(mode='list', 50)
for (i in 1:length(out_nocorr_lowvar_c_shed)) {
  store_maxc_low[[i]]<-max(out_nocorr_lowvar_c_shed[[i]][[2]]$c,na.rm=T)
}

######### cov-shed-contact-medvar neg corr
# epidemic duration of all simulations
store_epidemic_duration_med <- vector(mode='list', 50)
for (i in 1:length(out_nocorr_medvar_c_shed)) {
  store_epidemic_duration_med[[i]]<-max(out_nocorr_medvar_c_shed[[i]][[1]]$t)
}

# maximum shed values for each simulation
store_maxshed_med<-vector(mode='list', 50)
for (i in 1:length(out_nocorr_medvar_c_shed)) {
  store_maxshed_med[[i]]<-max(out_nocorr_medvar_c_shed[[i]][[2]]$shed,na.rm=T)
}

# maximum contact rates for each simulation
store_maxc_med<-vector(mode='list', 50)
for (i in 1:length(out_nocorr_medvar_c_shed)) {
  store_maxc_med[[i]]<-max(out_nocorr_medvar_c_shed[[i]][[2]]$c,na.rm=T)
}


hist(unlist(store_maxshed_low))
hist(unlist(store_maxc_low))
hist(unlist(store_maxshed_med))
hist(unlist(store_maxc_med))
