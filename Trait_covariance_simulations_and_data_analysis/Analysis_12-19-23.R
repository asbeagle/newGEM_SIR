## Stochastic run in Rcpp
sourceCpp("SIRcov.cpp")
library(deSolve) 
library(MASS)
library(tidyverse)
library(ggplot2)
library(magrittr)

## Compute the deterministic expectation for peak size for different parameter sets
SIRdet <- function(t, y, params) {
  S = y[1];
  I = y[2];
  R = y[3];
  b = params["b"];
  bs = params["bs"];
  d = params["d"]
  c = params["c"]
  s = params["s"]
  a = params["a"]
  g = params["g"]
  
  dSdt = (b-bs*(S+I+R))*(S+I+R) - c*s/(1+s)*S*I - d*S
  dIdt = c*s/(1+s)*S*I - (d+a+g)*I
  dRdt = g*I - d*R
  
  return(list(c(dSdt,dIdt,dRdt)))
}
detPeaks <- c()
for (R in c(8,4,1)) {
  if (R==8) x <- 0
  else if (R==4) x <- 0.24647
  else x <- 0.594804
  
  params = c(b=2.5, d=.1, bs=.01, c=(1-x)*0.1, s=(1-x)*1/9, a=(1+x)*0.1, g=(1+x)*0.1)
  initial_state <- floor(c(S=unname(((params["b"]-params["d"])/params["bs"]))-5, I=5, R=0))
  
  out <- ode(initial_state, seq(0,100,0.1), SIRdet, params)
  detPeaks <- c(detPeaks, max(out[,3]))
}

## Compute a stochastic, no variation expectation for peak size, fadeout probability, fadeout timing, and distribution of secondary infections
Rseq <- seq(1,8,0.1)
stochExp <- vector(mode='list', length=length(Rseq))
for (j in 1:length(Rseq)) {
  R <- Rseq[j]
  print(R)
  ## baseline parameters
  b = 2.5; bs=0.01; d=0.1; c=0.1; s=1/9; a=0.1; g=0.1
  ## scaling coefficient x to scale c, s, a, and g to achieve particular R value
  x = 1/(2*(-b*c*s+c*d*s-a*bs*R*s-bs*g*R*s)) * (-a*bs*R-bs*g*R-2*b*c*s+2*c*d*s+bs*d*R*s+sqrt((a*bs*R+bs*g*R+2*b*c*s-2*c*d*s-bs*d*R*s)^2-4*(-b*c*s+c*d*s-a*bs*R*s-bs*g*R*s)*(a*bs*R+bs*d*R+bs*g*R-b*c*s+c*d*s+a*bs*R*s+bs*d*R*s+bs*g*R*s)))
  ## simulation parameters
  params = c(b=2.5, bs=.01, d=.1, c=(1-x)*0.1, s=(1-x)*1/9, a=(1+x)*0.1, g=(1+x)*0.1)
  ## initial conditions
  initial_state <- floor(c(S=unname(((params["b"]-params["d"])/params["bs"]))-5, I=5, R=0))
  ## simulate stochastic model 1000 times
  lapply(1:1000, function(i) SIRnovar(params, initial_state, 50)) -> out
  ## Calculate the peak infected population size and the timing of fadeouts from each of the 1000 runs
  lapply(out, 
         function(o) 
           data.frame(peak=max(o[[2]][,2]),
                      fadeout=ifelse(any(o[[2]][,2]==0),1,0),
                      fadeoutT=ifelse(any(o[[2]][,2]==0),o[[1]][min(which(o[[2]][,2]==0))],50)) %>%
           mutate(., ## compute the dispersion parameter for each run that did not fade out
                  disp=ifelse(fadeout==0, glm.nb((o[[3]][is.na(o[[3]][,3]),2])~1)$theta, NA))) %>%
    do.call("rbind.data.frame",.) %>% 
    mutate(rep=1:1000, R0=R) -> stochExp[[j]]
}
stochExp %<>% do.call("rbind.data.frame",.)
saveRDS(stochExp, file="Stochastic_no_variation_simulations.RDS")

p1 <- stochExp %>% 
  ggplot(., aes(x=R0, y=peak)) + 
  stat_summary(fun.data = "median_hilow", geom="ribbon", fill=gray(0.75)) + 
  stat_summary(fun="mean", geom="line", col="red") + 
  ylab("Epidemic peak") + 
  xlab("") +
  theme_bw()

p2 <- stochExp %>% 
  ggplot(., aes(x=R0, y=disp)) + 
  stat_summary(fun.data = "median_hilow", geom="ribbon", fill=gray(0.75)) + 
  stat_summary(fun="mean", geom="line", col="red") + 
  ylab("Dispersion (k)") + 
  xlab("") + 
  theme_bw()

p3 <- stochExp %>% 
  group_by(R0) %>% 
  summarise(fadeprob=sum(fadeout)/1000) %>% 
  ggplot(., aes(x=R0, y=fadeprob)) + 
  geom_line(col="red") + 
  ylab("Prob of fadeout") + 
  xlab("R0") + 
  theme_bw() 

p4 <- stochExp %>% 
  ggplot(., aes(x=R0, y=fadeoutT)) + 
  stat_summary(fun.data = "median_hilow", geom="ribbon", fill=gray(0.75)) + 
  stat_summary(fun="mean", geom="line", col="red") + 
  ylab("Time to fadeout") + 
  xlab("R0") + 
  theme_bw()
  
library(patchwork)
png(filename="./Stochastic_no_variation_expectations.png", width=5, height=5, units='in', res=300)
(p1 + p2)/(p3 + p4)
dev.off()

outC <- outS <- outA <- outG <- vector(mode='list', length=9)
j = 1
for (R in c(1,4,8)) {
  for (CV in c(.2, 1, 5)) {
    print(CV)
    print(R)
    ## baseline parameters
    b = 2.5; bs=0.01; d=0.1; c=0.1; s=1/9; a=0.1; g=0.1
    ## scaling coefficient x to scale c, s, a, and g to achieve particular R value
    x = 1/(2*(-b*c*s+c*d*s-a*bs*R*s-bs*g*R*s)) * (-a*bs*R-bs*g*R-2*b*c*s+2*c*d*s+bs*d*R*s+sqrt((a*bs*R+bs*g*R+2*b*c*s-2*c*d*s-bs*d*R*s)^2-4*(-b*c*s+c*d*s-a*bs*R*s-bs*g*R*s)*(a*bs*R+bs*d*R+bs*g*R-b*c*s+c*d*s+a*bs*R*s+bs*d*R*s+bs*g*R*s)))
    ## simulation parameters
    params = c(b=2.5, bs=0.01, d=0.1, 
                c=(1-x)*0.1, s=(1-x)*1/9, a=(1+x)*0.1, g=(1+x)*0.1,
                SD=CV*(1-x)*0.1)
    ## initial conditions
    initial_state <- floor(c(S=unname(((params["b"]-params["d"])/params["bs"]))-5, I=5, R=0))
    ## simulate!
    lapply(1:1000, function(i) SIRvarC(params, initial_state, 50)) -> outC[[j]]
    lapply(1:1000, function(i) SIRvarS(params, initial_state, 50)) -> outS[[j]]
    lapply(1:1000, function(i) SIRvarA(params, initial_state, 50)) -> outA[[j]]
    lapply(1:1000, function(i) SIRvarG(params, initial_state, 50)) -> outG[[j]]
    j = j+1
  }
}
 
lapply(outC, 
       function(O)
         lapply(O, 
                function(o)
                  data.frame(peak=max(o[[2]][,2]),
                             totalI=sum(o[[3]][,3]),
                             fadeout=ifelse(any(o[[2]][,2]==0),1,0),
                             fadeoutT=ifelse(any(o[[2]][,2]==0),o[[1]][min(which(o[[2]][,2]==0))],50)) %>%
                  mutate(., ## compute the dispersion parameter for each run that did not fade out
                         disp=ifelse(fadeout==0, 
                                     ifelse(inherits(try(glm.nb((o[[3]][is.na(o[[3]][,4]),3])~1), silent=TRUE),"try-error"), 
                                            NA, 
                                            glm.nb((o[[3]][is.na(o[[3]][,4]),3])~1)$theta), 
                                     NA))
                ) %>% do.call("rbind.data.frame",.)
       ) -> OutC
OutC %<>% do.call("rbind.data.frame",.)
OutC$R0 <- rep(paste0("R0=",c(1,4,8)), each=3000)
OutC$CV <- rep(rep(paste0("CV=",c(0.2,1,5)), each=1000),3)
saveRDS(OutC, file="Variance_in_contact_results_12-28-23.RDS")

cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

png(file="Variance_in_contact_peaks.png", height=3, width=5, units='in', res=300)
ggplot(OutC, aes(x=peak, group=CV, fill=CV)) + 
  geom_histogram() + 
  facet_grid(~R0) + 
  scale_fill_manual(values=cbPalette) + 
  theme_bw() + 
  theme(legend.position="bottom")
dev.off()

png(file="Variance_in_contact_dispersion.png", height=3, width=5, units='in', res=300)
ggplot(OutC, aes(x=disp, group=CV, fill=CV)) + 
  geom_histogram() + 
  facet_grid(~R0) + 
  scale_fill_manual(values=cbPalette) +
  theme_bw() +
  theme(legend.position="bottom") +
  xlim(0,1.5)  
dev.off()


lapply(outS, 
       function(O)
         lapply(O, 
                function(o)
                  data.frame(peak=max(o[[2]][,2]),
                             totalI=sum(o[[3]][,3]),
                             fadeout=ifelse(any(o[[2]][,2]==0),1,0),
                             fadeoutT=ifelse(any(o[[2]][,2]==0),o[[1]][min(which(o[[2]][,2]==0))],50)) %>%
                  mutate(., ## compute the dispersion parameter for each run that did not fade out
                         disp=ifelse(fadeout==0, 
                                     ifelse(inherits(try(glm.nb((o[[3]][is.na(o[[3]][,4]),3])~1), silent=TRUE),"try-error"), 
                                            NA, 
                                            glm.nb((o[[3]][is.na(o[[3]][,4]),3])~1)$theta), 
                                     NA))
         ) %>% do.call("rbind.data.frame",.)
) -> OutS
OutS %<>% do.call("rbind.data.frame",.)
OutS$R0 <- rep(paste0("R0=",c(1,4,8)), each=3000)
OutS$CV <- rep(rep(paste0("CV=",c(0.2,1,5)), each=1000),3)
saveRDS(OutS, file="Variance_in_shedding_results_12-28-23.RDS")

png(file="Variance_in_shedding_peaks.png", height=3, width=5, units='in', res=300)
ggplot(OutS, aes(x=peak, group=CV, fill=CV)) + 
  geom_histogram() + 
  facet_grid(~R0) + 
  scale_fill_manual(values=cbPalette) + 
  theme_bw() + 
  theme(legend.position="bottom")
dev.off()

png(file="Variance_in_shedding_dispersion.png", height=3, width=5, units='in', res=300)
ggplot(OutS, aes(x=disp, group=CV, fill=CV)) + 
  geom_histogram() + 
  facet_grid(~R0) + 
  scale_fill_manual(values=cbPalette) +
  theme_bw() +
  theme(legend.position="bottom") +
  xlim(0,1.5)  
dev.off()



lapply(outA, 
       function(O)
         lapply(O, 
                function(o)
                  data.frame(peak=max(o[[2]][,2]),
                             totalI=sum(o[[3]][,3]),
                             fadeout=ifelse(any(o[[2]][,2]==0),1,0),
                             fadeoutT=ifelse(any(o[[2]][,2]==0),o[[1]][min(which(o[[2]][,2]==0))],50)) %>%
                  mutate(., ## compute the dispersion parameter for each run that did not fade out
                         disp=ifelse(fadeout==0, 
                                     ifelse(inherits(try(glm.nb((o[[3]][is.na(o[[3]][,4]),3])~1), silent=TRUE),"try-error"), 
                                            NA, 
                                            glm.nb((o[[3]][is.na(o[[3]][,4]),3])~1)$theta), 
                                     NA))
         ) %>% do.call("rbind.data.frame",.)
) -> OutA
OutA %<>% do.call("rbind.data.frame",.)
OutA$R0 <- rep(paste0("R0=",c(1,4,8)), each=3000)
OutA$CV <- rep(rep(paste0("CV=",c(0.2,1,5)), each=1000),3)
saveRDS(OutA, file="Variance_in_virulence_results_12-28-23.RDS")

png(file="Variance_in_virulence_peaks.png", height=3, width=5, units='in', res=300)
ggplot(OutA, aes(x=peak, group=CV, fill=CV)) + 
  geom_histogram() + 
  facet_grid(~R0) + 
  scale_fill_manual(values=cbPalette) + 
  theme_bw() + 
  theme(legend.position="bottom")
dev.off()

png(file="Variance_in_virulence_dispersion.png", height=3, width=5, units='in', res=300)
ggplot(OutA, aes(x=disp, group=CV, fill=CV)) + 
  geom_histogram() + 
  facet_grid(~R0) + 
  scale_fill_manual(values=cbPalette) +
  theme_bw() +
  theme(legend.position="bottom") +
  xlim(0,1.5)  
dev.off()



lapply(outG, 
       function(O)
         lapply(O, 
                function(o)
                  data.frame(peak=max(o[[2]][,2]),
                             totalI=sum(o[[3]][,3]),
                             fadeout=ifelse(any(o[[2]][,2]==0),1,0),
                             fadeoutT=ifelse(any(o[[2]][,2]==0),o[[1]][min(which(o[[2]][,2]==0))],50)) %>%
                  mutate(., ## compute the dispersion parameter for each run that did not fade out
                         disp=ifelse(fadeout==0, 
                                     ifelse(inherits(try(glm.nb((o[[3]][is.na(o[[3]][,4]),3])~1), silent=TRUE),"try-error"), 
                                            NA, 
                                            glm.nb((o[[3]][is.na(o[[3]][,4]),3])~1)$theta), 
                                     NA))
         ) %>% do.call("rbind.data.frame",.)
) -> OutG
OutG %<>% do.call("rbind.data.frame",.)
OutG$R0 <- rep(paste0("R0=",c(1,4,8)), each=3000)
OutG$CV <- rep(rep(paste0("CV=",c(0.2,1,5)), each=1000),3)
saveRDS(OutG, file="Variance_in_recovery_results_12-28-23.RDS")

png(file="Variance_in_recovery_peaks.png", height=3, width=5, units='in', res=300)
ggplot(OutG, aes(x=peak, group=CV, fill=CV)) + 
  geom_histogram() + 
  facet_grid(~R0) + 
  scale_fill_manual(values=cbPalette) + 
  theme_bw() + 
  theme(legend.position="bottom")
dev.off()

png(file="Variance_in_recovery_dispersion.png", height=3, width=5, units='in', res=300)
ggplot(OutG, aes(x=disp, group=CV, fill=CV)) + 
  geom_histogram() + 
  facet_grid(~R0) + 
  scale_fill_manual(values=cbPalette) +
  theme_bw() +
  theme(legend.position="bottom") +
  xlim(0,1.5)  
dev.off()




         

