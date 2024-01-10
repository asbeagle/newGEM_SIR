## example plots for SS MS

# r04, nocorr, medvar, contact/shed
dat<-readRDS("out_R=4_nocorr_medvar_c_shed.RDS")

lapply(1:length(z1), function(i) data.frame(peakSize=max(z1[[i]][[1]]$I,na.rm=TRUE), ## peak epidemic size
                                            peakPrev=max(z1[[i]][[1]]$I/(z1[[i]][[1]]$S+z1[[i]][[1]]$I+z1[[i]][[1]]$R),na.rm=TRUE),
                                            disp=ifelse(max(z1[[i]][[1]]$t)>99, glm.nb(z1[[i]][[2]]$numInf~1)$theta, NA), ## dispersion parameter of negative binomial distribution fit to numInf
                                            fadeout=ifelse(length(z1[[i]][[1]]$t)<99, 1, 0), # binary: did this replicate go extinct?
                                            rep=i))

data %>% group_by(Variance, cov, traits) %>% reframe(size=seq(0,max(data$peakSize)+2), ECDF=sapply(seq(0,max(data$peakSize)+2), function(i) sum(peakSize>=i)/100)) -> data2

# plot 20 sims of I over time
subset_dat<-lapply(1:20, function(i) data.frame(numI = dat[[i]][[1]]$I,
                                            time=round(dat[[i]][[1]]$t),
                                            rep=i))


ggplot(bind_rows(subset_dat, .id="rep"), 
                aes(time, numI, color=rep))+
  geom_line(show.legend = FALSE)+
  theme_bw()+
  labs(y="Number Infected", x= "Time")


