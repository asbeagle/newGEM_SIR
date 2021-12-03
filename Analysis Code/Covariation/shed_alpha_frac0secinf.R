# frac 0 secondary infection

## shed alpha
# hi var
# pos
for (var in c("low","med","hi")) {
  for (corr in c("no","neg","pos")) {
    data <- readRDS(paste0("out_",corr,"corr_",var,"var_shed_alpha.RDS"))
    lapply(data, function(d)
      data.frame(R=d[[2]]$numInf,
                 frac0Reff= sum(d[[2]]$numInf==0)/nrow(d[[2]]),
                 Corr = as.factor(corr),
                 Var = as.factor(var)))  -> sumdata1
  }
}
sumdata_use1 <-do.call("rbind.data.frame", sumdata1)

sumdata_use1 %>% 
  mutate(R = mean(R), frac0Reff = mean(frac0Reff))%>% 
  distinct(.)-> sumdata_use1


# neg
for (var in c("low","med","hi")) {
  for (corr in c("no","pos","neg")) {
    data <- readRDS(paste0("out_",corr,"corr_",var,"var_shed_alpha.RDS"))
    lapply(data, function(d)
      data.frame(R=d[[2]]$numInf,
                 frac0Reff= sum(d[[2]]$numInf==0)/nrow(d[[2]]),
                 Corr = as.factor(corr),
                 Var = as.factor(var)))  -> sumdata2
  }
}
sumdata_use2 <-do.call("rbind.data.frame", sumdata2)

sumdata_use2 %>% 
  mutate(R = mean(R), frac0Reff = mean(frac0Reff))%>%
  distinct(.)-> sumdata_use2

# no
for (var in c("low","med","hi")) {
  for (corr in c("pos","neg","no")) {
    data <- readRDS(paste0("out_",corr,"corr_",var,"var_shed_alpha.RDS"))
    lapply(data, function(d)
      data.frame(R=d[[2]]$numInf,
                 frac0Reff= sum(d[[2]]$numInf==0)/nrow(d[[2]]),
                 Corr = as.factor(corr),
                 Var = as.factor(var)))  -> sumdata3
  }
}
sumdata_use3 <-do.call("rbind.data.frame", sumdata3)

sumdata_use3 %>% 
  mutate(R = mean(R), frac0Reff = mean(frac0Reff))%>%
  distinct(.)-> sumdata_use3

## med var
# pos
for (var in c("low","hi","med")) {
  for (corr in c("no","neg","pos")) {
    data <- readRDS(paste0("out_",corr,"corr_",var,"var_shed_alpha.RDS"))
    lapply(data, function(d)
      data.frame(R=d[[2]]$numInf,
                 frac0Reff= sum(d[[2]]$numInf==0)/nrow(d[[2]]),
                 Corr = as.factor(corr),
                 Var = as.factor(var)))  -> sumdata4
  }
}
sumdata_use4 <-do.call("rbind.data.frame", sumdata4)

sumdata_use4 %>% 
  mutate(R = mean(R), frac0Reff = mean(frac0Reff))%>%
  distinct(.)-> sumdata_use4

# neg
for (var in c("low","hi","med")) {
  for (corr in c("no","pos","neg")) {
    data <- readRDS(paste0("out_",corr,"corr_",var,"var_shed_alpha.RDS"))
    lapply(data, function(d)
      data.frame(R=d[[2]]$numInf,
                 frac0Reff= sum(d[[2]]$numInf==0)/nrow(d[[2]]),
                 Corr = as.factor(corr),
                 Var = as.factor(var)))  -> sumdata5
  }
}
sumdata_use5 <-do.call("rbind.data.frame", sumdata5)

sumdata_use5 %>% 
  mutate(R = mean(R), frac0Reff = mean(frac0Reff))%>%
  distinct(.)-> sumdata_use5

# no
for (var in c("low","hi","med")) {
  for (corr in c("pos","neg","no")) {
    data <- readRDS(paste0("out_",corr,"corr_",var,"var_shed_alpha.RDS"))
    lapply(data, function(d)
      data.frame(R=d[[2]]$numInf,
                 frac0Reff= sum(d[[2]]$numInf==0)/nrow(d[[2]]),
                 Corr = as.factor(corr),
                 Var = as.factor(var)))  -> sumdata6
  }
}
sumdata_use6 <-do.call("rbind.data.frame", sumdata6)

sumdata_use6 %>% 
  mutate(R = mean(R), frac0Reff = mean(frac0Reff))%>%
  distinct(.)-> sumdata_use6


## low var
# pos
for (var in c("hi","med","low")) {
  for (corr in c("no","neg","pos")) {
    data <- readRDS(paste0("out_",corr,"corr_",var,"var_shed_alpha.RDS"))
    lapply(data, function(d)
      data.frame(R=d[[2]]$numInf,
                 frac0Reff= sum(d[[2]]$numInf==0)/nrow(d[[2]]),
                 Corr = as.factor(corr),
                 Var = as.factor(var)))  -> sumdata7
  }
}
sumdata_use7 <-do.call("rbind.data.frame", sumdata7)

sumdata_use7 %>% 
  mutate(R = mean(R), frac0Reff = mean(frac0Reff))%>%
  distinct(.)-> sumdata_use7

# neg
for (var in c("hi","med","low")) {
  for (corr in c("no","pos","neg")) {
    data <- readRDS(paste0("out_",corr,"corr_",var,"var_shed_alpha.RDS"))
    lapply(data, function(d)
      data.frame(R=d[[2]]$numInf,
                 frac0Reff= sum(d[[2]]$numInf==0)/nrow(d[[2]]),
                 Corr = as.factor(corr),
                 Var = as.factor(var)))  -> sumdata8
  }
}
sumdata_use8 <-do.call("rbind.data.frame", sumdata8)

sumdata_use8 %>% 
  mutate(R = mean(R), frac0Reff = mean(frac0Reff))%>%
  distinct(.)-> sumdata_use8

# no
for (var in c("hi","med","low")) {
  for (corr in c("pos","neg","no")) {
    data <- readRDS(paste0("out_",corr,"corr_",var,"var_shed_alpha.RDS"))
    lapply(data, function(d)
      data.frame(R=d[[2]]$numInf,
                 frac0Reff= sum(d[[2]]$numInf==0)/nrow(d[[2]]),
                 Corr = as.factor(corr),
                 Var = as.factor(var)))  -> sumdata9
  }
}
sumdata_use9 <-do.call("rbind.data.frame", sumdata9)

sumdata_use9 %>% 
  mutate(R = mean(R), frac0Reff = mean(frac0Reff))%>%
  distinct(.)-> sumdata_use9

shed_alpha_data <-rbind(sumdata_use1,sumdata_use2,sumdata_use3,
                           sumdata_use4,sumdata_use5,sumdata_use6,
                           sumdata_use7,sumdata_use8,sumdata_use9)

frac0_shed_alpha <- rbind(sumdata_use1,sumdata_use2,sumdata_use3,sumdata_use4,sumdata_use5,sumdata_use6,sumdata_use7,sumdata_use8,sumdata_use9)

# plot
shed_alpha<-ggplot(data=frac0_shed_alpha, mapping = aes(x=Corr, y = frac0Reff))+
  geom_point(data = frac0_shed_alpha, aes(color=Var, shape = Var), alpha = .7, size = 3.5)+
  scale_color_manual(values=c("darkblue", "darkgreen", "pink"))+
  labs(title = "Shed-Alpha", y = "% no secondary infections")+
  scale_x_discrete(limits=c("no", "neg","pos"))

shed_alpha<- ggplot(shed_alpha_data, mapping = aes(x=frac0Reff, y=R))+
  geom_point(aes(color=Var, shape = Corr), size=4, alpha = .85)+
  scale_color_manual(values=c("darkblue", "darkgreen", "pink"))+
  labs(title="shed-alpha")+
  labs(x = "mean % 0 Reff", y = "mean Reff")+
  theme_bw()

shed_alpha
