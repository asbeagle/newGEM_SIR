library(tidyverse
        )

SC_hivar_c <- array(NA, dim=c(length(timeSeq), length(out_negcorr_hivar_c_shed)+1))
SC_hivar_c[,1] <- timeSeq
for (i in 1:length(SC_hivar_c)) SC_hivar_c[,i+1] <- out_negcorr_hivar_c_shed[[i]][[1]]$c


## Plots you can make
## Easy stuff like population dynamics
plot(example3[[1]]$t, example3[[1]]$S, type='l', xlab='Time', ylab='S(t)')
## More interesting stuff like histograms of fitness (probably should censor any individual that was still alive, since its fitness is incomplete)
with(subset(example3[[2]], !is.infinite(tInf)), hist(numInf, breaks=max(numInf)+1)) ## most individuals have 0 fitness

plot(0:100, lapply(out_negcorr_hivar_c_shed, function(l) l[[1]]$S) %>% do.call("cbind.data.frame",.) %>% apply(., 1, median), lwd=2, type="l")

#### high var neg corr contact shed
store_c<-vector(mode='list', 50)
for (i in 1:length(out_negcorr_hivar_c_shed)) {
  store_c[[i]]<-out_negcorr_hivar_c_shed[[i]][[2]]$c
}

store_shed<-vector(mode='list', 50)
for (i in 1:length(out_negcorr_hivar_c_shed)) {
  store_shed[[i]]<-out_negcorr_hivar_c_shed[[i]][[2]]$shed
}

store_epidemic_duration <- vector(mode='list', 50)
for (i in 1:length(out_negcorr_hivar_c_shed)) {
  store_epidemic_duration[[i]]<-max(out_negcorr_hivar_c_shed[[i]][[1]]$t)
}

store_maxshed<-vector(mode='list', 50)
for (i in 1:length(out_negcorr_hivar_c_shed)) {
  store_maxshed[[i]]<-max(out_negcorr_hivar_c_shed[[i]][[2]]$shed,na.rm=T)
}

store_maxc<-vector(mode='list', 50)
for (i in 1:length(out_negcorr_hivar_c_shed)) {
  store_maxc[[i]]<-max(out_negcorr_hivar_c_shed[[i]][[2]]$c,na.rm=T)
}


plot(unlist(store_c)*unlist(store_shed)/(1+unlist(store_shed)), unlist(store_numInf), xlim=c(0,0.1))


store_numInf<-vector(mode='list', 50)
for (i in 1:length(out_negcorr_hivar_c_shed)) {
  store_numInf[[i]]<-out_negcorr_hivar_c_shed[[i]][[2]]$numInf
}


breaks<-seq(from=0, to=2.5, by=.05)
hist(unlist(store_c), breaks=breaks, main="highvar-c-shed-negcorr", xlab="c",ylim=c(0,50))

store_s<-vector(mode='list', 50)
for (i in 1:length(out_negcorr_hivar_c_shed)) {
  store_s[i]<-out_negcorr_hivar_c_shed[[i]][[2]]$shed
}
hist(unlist(store_s), breaks=breaks, main="highvar-c-shed-negcorr", xlab="shed",ylim=c(0,50))

par(mfrow=c(3,2))

### low var neg corr contact shed
store_c_lowvar_negcor<-vector(mode='list', 50)
for (i in 1:length(out_negcorr_lowvar_c_shed)) {
  store_c_lowvar_negcor[i]<-out_negcorr_lowvar_c_shed[[i]][[2]]$c
}

hist(unlist(store_c_lowvar_negcor), breaks=breaks, main="lowvar-c-shed-negcorr", xlab="c",xlim=c(0,.5))

store_s_lowvar_negcor<-vector(mode='list', 50)
for (i in 1:length(out_negcorr_lowvar_c_shed)) {
  store_s_lowvar_negcor[i]<-out_negcorr_lowvar_c_shed[[i]][[2]]$shed
}
hist(unlist(store_s_lowvar_negcor), breaks=breaks, main="lowvar-c-shed-negcorr", xlab="shed",xlim=c(0,.5), ylim=c(0,30))

### med var neg corr contact shed
store_c_medvar_negcor<-vector(mode='list', 50)
for (i in 1:length(out_negcorr_medvar_c_shed)) {
  store_c_medvar_negcor[i]<-out_negcorr_medvar_c_shed[[i]][[2]]$c
}

hist(unlist(store_c_medvar_negcor), breaks=breaks, main="medvar-c-shed-negcorr", xlab="c", xlim=c(0,1), ylim=c(0,40))

store_s_medvar_negcor<-vector(mode='list', 50)
for (i in 1:length(out_negcorr_medvar_c_shed)) {
  store_s_medvar_negcor[i]<-out_negcorr_medvar_c_shed[[i]][[2]]$shed
}
hist(unlist(store_s_medvar_negcor), breaks=breaks, main="medvar-c-shed-negcorr", xlab="shed",xlim=c(0,1),ylim=c(0,40))

## all together plots
hist(unlist(store_c_lowvar_negcor), breaks=breaks, main="lowvar-c-shed-negcorr", xlab="c",xlim=c(0,2.5),ylim=c(0,40))
hist(unlist(store_s_lowvar_negcor), breaks=breaks, main="lowvar-c-shed-negcorr", xlab="shed",xlim=c(0,2.5), ylim=c(0,40))

hist(unlist(store_c_medvar_negcor), breaks=breaks, main="medvar-c-shed-negcorr", xlab="c", xlim=c(0,2.5), ylim=c(0,40))
hist(unlist(store_s_medvar_negcor), breaks=breaks, main="medvar-c-shed-negcorr", xlab="shed",xlim=c(0,2.5),ylim=c(0,40))

hist(unlist(store_c), breaks=breaks, main="highvar-c-shed-negcorr", xlab="c",xlim=c(0,2.5), ylim=c(0,40))
hist(unlist(store_s), breaks=breaks, main="highvar-c-shed-negcorr", xlab="shed",xlim=c(0,2.5), ylim=c(0,40))

#### high var no corr contact shed
store_c_hivar_nocor<-vector(mode='list', 50)
for (i in 1:length(out_nocorr_hivar_c_shed)) {
  store_c_hivar_nocor[i]<-out_nocorr_hivar_c_shed[[i]][[2]]$c
}

store_s_hivar_nocor<-vector(mode='list', 50)
for (i in 1:length(out_nocorr_hivar_c_shed)) {
  store_s_hivar_nocor[i]<-out_nocorr_hivar_c_shed[[i]][[2]]$shed
}
par(mfrow=c(3,2))

### low var no corr contact shed
store_c_lowvar_nocor<-vector(mode='list', 50)
for (i in 1:length(out_nocorr_lowvar_c_shed)) {
  store_c_lowvar_nocor[i]<-out_nocorr_lowvar_c_shed[[i]][[2]]$c
}

store_s_lowvar_nocor<-vector(mode='list', 50)
for (i in 1:length(out_nocorr_lowvar_c_shed)) {
  store_s_lowvar_nocor[i]<-out_nocorr_lowvar_c_shed[[i]][[2]]$shed
}

### med var neg corr contact shed
store_c_medvar_nocor<-vector(mode='list', 50)
for (i in 1:length(out_nocorr_medvar_c_shed)) {
  store_c_medvar_nocor[i]<-out_nocorr_medvar_c_shed[[i]][[2]]$c
}

store_s_medvar_nocor<-vector(mode='list', 50)
for (i in 1:length(out_nocorr_medvar_c_shed)) {
  store_s_medvar_nocor[i]<-out_nocorr_medvar_c_shed[[i]][[2]]$shed
}

## all together plots
hist(unlist(store_c_lowvar_nocor), breaks=breaks, main="lowvar-c-shed-nocorr", xlab="c",xlim=c(0,2.5),ylim=c(0,45))
hist(unlist(store_s_lowvar_nocor), breaks=breaks, main="lowvar-c-shed-nocorr", xlab="shed",xlim=c(0,2.5), ylim=c(0,45))

hist(unlist(store_c_medvar_nocor), breaks=breaks, main="medvar-c-shed-nocorr", xlab="c", xlim=c(0,2.5), ylim=c(0,45))
hist(unlist(store_s_medvar_nocor), breaks=breaks, main="medvar-c-shed-nocorr", xlab="shed",xlim=c(0,2.5),ylim=c(0,45))

hist(unlist(store_c_hivar_nocor), breaks=breaks, main="highvar-c-shed-nocorr", xlab="c",xlim=c(0,2.5), ylim=c(0,45))
hist(unlist(store_s_hivar_nocor), breaks=breaks, main="highvar-c-shed-nocorr", xlab="shed",xlim=c(0,2.5), ylim=c(0,45))

#### high var pos corr contact shed
store_c_hivar_poscor<-vector(mode='list', 50)
for (i in 1:length(out_poscorr_hivar_c_shed)) {
  store_c_hivar_poscor[i]<-out_poscorr_hivar_c_shed[[i]][[2]]$c
}

store_s_hivar_poscor<-vector(mode='list', 50)
for (i in 1:length(out_poscorr_hivar_c_shed)) {
  store_s_hivar_poscor[i]<-out_poscorr_hivar_c_shed[[i]][[2]]$shed
}
par(mfrow=c(3,2))

### low var pos corr contact shed
store_c_lowvar_poscor<-vector(mode='list', 50)
for (i in 1:length(out_poscorr_lowvar_c_shed)) {
  store_c_lowvar_poscor[i]<-out_poscorr_lowvar_c_shed[[i]][[2]]$c
}

store_s_lowvar_poscor<-vector(mode='list', 2000)
store_s_lowvar_poscor[1]<-
for (i in 1:length(out_poscorr_lowvar_c_shed)) {
  store_s_lowvar_poscor[i+1]<-out_poscorr_lowvar_c_shed[[i+1]][[2]]$shed
}

### med var pos corr contact shed
store_c_medvar_nocor<-vector(mode='list', 50)
for (i in 1:length(out_nocorr_medvar_c_shed)) {
  store_c_medvar_nocor[i]<-out_nocorr_medvar_c_shed[[i]][[2]]$c
}

store_s_medvar_nocor<-vector(mode='list', 50)
for (i in 1:length(out_nocorr_medvar_c_shed)) {
  store_s_medvar_nocor[i]<-out_nocorr_medvar_c_shed[[i]][[2]]$shed
}



## all together plots
hist(unlist(store_c_lowvar_nocor), breaks=breaks, main="lowvar-c-shed-nocorr", xlab="c",xlim=c(0,2.5),ylim=c(0,45))
hist(unlist(store_s_lowvar_nocor), breaks=breaks, main="lowvar-c-shed-nocorr", xlab="shed",xlim=c(0,2.5), ylim=c(0,45))

hist(unlist(store_c_medvar_nocor), breaks=breaks, main="medvar-c-shed-nocorr", xlab="c", xlim=c(0,2.5), ylim=c(0,45))
hist(unlist(store_s_medvar_nocor), breaks=breaks, main="medvar-c-shed-nocorr", xlab="shed",xlim=c(0,2.5),ylim=c(0,45))

hist(unlist(store_c_hivar_nocor), breaks=breaks, main="highvar-c-shed-nocorr", xlab="c",xlim=c(0,2.5), ylim=c(0,45))
hist(unlist(store_s_hivar_nocor), breaks=breaks, main="highvar-c-shed-nocorr", xlab="shed",xlim=c(0,2.5), ylim=c(0,45))


par(mfrow=c(1,1))

hist(out_negcorr_hivar_c_shed[[20]][[2]]$shed, xlim=c(0,1), ylim=c(0,40), breaks=breaks)
hist(out_negcorr_hivar_c_shed[[20]][[2]]$c, xlim=c(0,1), ylim=c(0,40),breaks=breaks)

hist(out_negcorr_lowvar_c_shed[[20]][[2]]$shed, xlim=c(0,1), breaks=breaks)
hist(out_negcorr_lowvar_c_shed[[20]][[2]]$c, xlim=c(0,1),breaks=breaks)

hist(out_negcorr_medvar_c_shed[[20]][[2]]$shed, xlim=c(0,1), breaks=breaks)
hist(out_negcorr_medvar_c_shed[[20]][[2]]$c, xlim=c(0,1),breaks=breaks)

##pos corr
hist(out_poscorr_hivar_c_shed[[20]][[2]]$shed,breaks=seq(from=0,to=7, by=0.1))
hist(out_poscorr_hivar_c_shed[[20]][[2]]$c)

hist(out_poscorr_medvar_c_shed[[20]][[2]]$shed)
hist(out_poscorr_medvar_c_shed[[20]][[2]]$c)

hist(out_poscorr_lowvar_c_shed[[20]][[2]]$shed)
hist(out_poscorr_lowvar_c_shed[[20]][[2]]$c)

##no cor
hist(out_nocorr_hivar_c_shed[[10]][[2]]$c)
hist(out_nocorr_hivar_c_shed[[11]][[2]]$c)
hist(out_nocorr_hivar_c_shed[[12]][[2]]$c)
hist(out_nocorr_hivar_c_shed[[13]][[2]]$c)
hist(out_nocorr_hivar_c_shed[[14]][[2]]$c)
hist(out_nocorr_hivar_c_shed[[15]][[2]]$c)
hist(out_nocorr_hivar_c_shed[[16]][[2]]$c)
hist(out_nocorr_hivar_c_shed[[17]][[2]]$c)
hist(out_nocorr_hivar_c_shed[[18]][[2]]$c)
hist(out_nocorr_hivar_c_shed[[19]][[2]]$c)
hist(out_nocorr_hivar_c_shed[[20]][[2]]$c)
hist(out_nocorr_hivar_c_shed[[20]][[2]]$shed)

hist(out_negcorr_hivar_c_shed[[20]][[2]]$c)
hist(out_negcorr_hivar_c_shed[[20]][[2]]$shed)

hist(out_poscorr_hivar_c_shed[[20]][[2]]$c)
hist(out_poscorr_hivar_c_shed[[20]][[2]]$shed)

hist(out_nocorr_lowvar_c_shed[[10]][[2]]$c)
hist(out_nocorr_lowvar_c_shed[[11]][[2]]$c)
hist(out_nocorr_lowvar_c_shed[[12]][[2]]$c)
hist(out_nocorr_lowvar_c_shed[[13]][[2]]$c)
hist(out_nocorr_lowvar_c_shed[[14]][[2]]$c)
hist(out_nocorr_lowvar_c_shed[[15]][[2]]$c)
hist(out_nocorr_lowvar_c_shed[[16]][[2]]$c)
hist(out_nocorr_lowvar_c_shed[[17]][[2]]$c)
hist(out_nocorr_lowvar_c_shed[[18]][[2]]$c)
hist(out_nocorr_lowvar_c_shed[[19]][[2]]$c)
hist(out_nocorr_lowvar_c_shed[[20]][[2]]$c)

hist(out_poscorr_hivar_c_shed[[10]][[2]]$c)
hist(out_poscorr_hivar_c_shed[[11]][[2]]$c)
hist(out_poscorr_hivar_c_shed[[12]][[2]]$c)
hist(out_poscorr_hivar_c_shed[[13]][[2]]$c)
hist(out_poscorr_hivar_c_shed[[14]][[2]]$c)
hist(out_poscorr_hivar_c_shed[[15]][[2]]$c)
hist(out_poscorr_hivar_c_shed[[16]][[2]]$c)
hist(out_poscorr_hivar_c_shed[[17]][[2]]$c)
hist(out_poscorr_hivar_c_shed[[18]][[2]]$c)
hist(out_poscorr_hivar_c_shed[[19]][[2]]$c)
hist(out_poscorr_hivar_c_shed[[20]][[2]]$c)

#### covariation contact and virulence alpha
hist(out_poscorr_hivar_c_alpha[[20]][[2]]$c)
hist(out_poscorr_hivar_c_alpha[[19]][[2]]$c)
hist(out_poscorr_hivar_c_alpha[[18]][[2]]$c)
hist(out_poscorr_hivar_c_alpha[[17]][[2]]$c)
hist(out_poscorr_hivar_c_alpha[[16]][[2]]$c)
hist(out_poscorr_hivar_c_alpha[[16]][[2]]$c)



