### shed alpha

# epidemic duration of all simulations
store_epidemic_duration <- vector(mode='list', 50)
for (i in 1:length(out_nocorr_hivar_shed_alpha)) {
  store_epidemic_duration[[i]]<-max(out_nocorr_hivar_shed_alpha[[i]][[1]]$t)
}

# maximum shed values for each simulation
store_maxshed<-vector(mode='list', 50)
for (i in 1:length(out_nocorr_hivar_shed_alpha)) {
  store_maxshed[[i]]<-max(out_nocorr_hivar_shed_alpha[[i]][[2]]$shed,na.rm=T)
}

# maximum contact rates for each simulation
store_maxc<-vector(mode='list', 50)
for (i in 1:length(out_nocorr_hivar_shed_alpha)) {
  store_maxc[[i]]<-max(out_nocorr_hivar_shed_alpha[[i]][[2]]$c,na.rm=T)
}

# max secondary infection in each simulation
store_maxnumInf<-vector(mode='list', 50)
for (i in 1:length(out_nocorr_hivar_shed_alpha)) {
  store_maxnumInf[[i]]<-max(out_nocorr_hivar_shed_alpha[[i]][[2]]$numInf, na.rm=T)
}


par(mfrow=c(1,2))
plot(unlist(store_epidemic_duration), unlist(store_maxshed), ylab="max shed", xlab="epidemic duration")
plot(unlist(store_epidemic_duration), unlist(store_maxc), ylab="max c", xlab="epidemic duration")

shed_alpha_cor<-cor.test(unlist(store_epidemic_duration),unlist(store_maxshed))
cor.test(unlist(store_epidemic_duration),unlist(store_maxc))



# epidemic duration of all simulations
store_epidemic_duration <- vector(mode='list', 50)
for (i in 1:length(out_nocorr_hivar_c_alpha)) {
  store_epidemic_duration[[i]]<-max(out_nocorr_hivar_c_alpha[[i]][[1]]$t)
}

# maximum contact rates for each simulation
store_maxc<-vector(mode='list', 50)
for (i in 1:length(out_nocorr_hivar_c_alpha)) {
  store_maxc[[i]]<-max(out_nocorr_hivar_c_alpha[[i]][[2]]$c,na.rm=T)
}

contact_alpha_cor<-cor.test(unlist(store_epidemic_duration),unlist(store_maxc))

plot(unlist(store_epidemic_duration), unlist(store_maxc), ylab="max c", xlab="epidemic duration")

