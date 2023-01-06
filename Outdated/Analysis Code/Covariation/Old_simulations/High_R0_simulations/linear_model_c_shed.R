### Covariation between contact rate and shedding rate

#### HIGH R0, LM OF CONTACT/EPIDEMIC DURATION, SHEDDING/EPIDEMIC DURATION

# epidemic duration of all simulations
store_epidemic_duration <- vector(mode='list', 50)
for (i in 1:length(out_nocorr_hivar_c_shed)) {
  store_epidemic_duration[[i]]<-max(out_nocorr_hivar_c_shed[[i]][[1]]$t)
}

# maximum shed values for each simulation
store_maxshed<-vector(mode='list', 50)
for (i in 1:length(out_nocorr_hivar_c_shed)) {
  store_maxshed[[i]]<-max(out_nocorr_hivar_c_shed[[i]][[2]]$shed,na.rm=T)
}

# maximum contact rates for each simulation
store_maxc<-vector(mode='list', 50)
for (i in 1:length(out_nocorr_hivar_c_shed)) {
  store_maxc[[i]]<-max(out_nocorr_hivar_c_shed[[i]][[2]]$c,na.rm=T)
}

# max secondary infection in each simulation
store_maxnumInf<-vector(mode='list', 50)
for (i in 1:length(out_nocorr_hivar_c_shed)) {
  store_maxnumInf[[i]]<-max(out_nocorr_hivar_c_shed[[i]][[2]]$numInf, na.rm=T)
}

par(mfrow=c(1,2))
plot(unlist(store_epidemic_duration), unlist(store_maxshed), ylab="max shed", xlab="epidemic duration")
plot(unlist(store_epidemic_duration), unlist(store_maxc), ylab="max c", xlab="epidemic duration")

par(mfrow=c(1,1))
plot(unlist(store_epidemic_duration), unlist(store_maxnumInf), ylab="max num inf", xlab="epidemic duration")

cor.test(unlist(store_epidemic_duration),unlist(store_maxnumInf))
summary(lm(unlist(store_epidemic_duration)~unlist(store_maxnumInf)))

cor.test(unlist(store_epidemic_duration),unlist(store_maxshed))
cor.test(unlist(store_epidemic_duration),unlist(store_maxc))


# linear model 
# max shed value against epidemic duration
maxshed_lm_lowR0<-summary(lm(unlist(store_epidemic_duration)~unlist(store_maxshed)))

# max contact values against epidemic duration
maxc_lm_lowR0<-summary(lm(unlist(store_epidemic_duration)~unlist(store_maxc)))


hist(unlist(store_maxshed))
hist(unlist(store_maxc), xlim=c(0,15))