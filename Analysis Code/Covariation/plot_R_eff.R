### contact and shedding
## get data frame with all variance level and corr
library(tidyverse)
library(ggplot2)


oo <- vector(mode='list', length=9)
i <- 1
for (var in c("low","med","hi")) {
  for (corr in c("no","neg","pos")) {
    dataName <- paste0("out_",corr,"corr_",var,"var_c_shed")
    oo[[i]] <- data.frame(R=lapply(get(dataName), function(l) l[[2]]$numInf==0) %>% unlist,
                          Var=var,
                          Corr=corr)
    i <- i + 1
  }
}
ooo <- do.call("rbind.data.frame", oo)

ooo$Var <- factor(ooo$Var, levels = c("low", "med", "hi"))

ggplot(data=ooo, mapping = aes(x=Var, y= (length(R==0))/(length(R))))+
  geom_point(data=ooo, color=Corr)

q<-ggplot(data=ooo, mapping=aes(x=Var, y=log(R+1)))+
  geom_boxplot(data=ooo, aes(color=Corr))+
  scale_x_discrete(limits=c("low", "med","hi"))+
  scale_color_manual(values=c("pink", "darkgreen", "darkblue"))+
  labs(title="contact-shed")+
  geom_hline(yintercept = log(2), linetype = "dotdash", color = "black", alpha = .8)
q

p<-ggplot(data=ooo, mapping=aes(x=Corr, y=log(R+1)))+
  geom_boxplot(data=ooo, aes(color=Var))+
  scale_x_discrete(limits=c("no", "neg","pos"))+
  scale_color_manual(values=c("pink", "darkgreen", "darkblue"))+
  labs(title="contact-shed")+
  geom_hline(yintercept = log(2), linetype = "dotdash", color = "black", alpha = .8)
p

ooo %>%
  group_by(Corr,Var)%>%
  summarize(mean=mean(R),var=var(R))-> ooo2

ooo %>%
  group_by(Corr,Var)%>%
  summarize(frac0= (length(which(ooo$R==0)))/(length(ooo$R)))-> ooo2

# contact alpha
## get data frame with all variance level and corr
oo <- vector(mode='list', length=9)
i <- 1
for (var in c("low","med","hi")) {
  for (corr in c("no","neg","pos")) {
    dataName <- paste0("out_",corr,"corr_",var,"var_c_alpha")
    oo[[i]] <- data.frame(R=lapply(get(dataName), function(l) l[[2]]$numInf) %>% unlist,
                          Var=var,
                          Corr=corr,
                          R_0 = lapply(get(dataName), function(l) length(l[[2]]$numInf==0) %>% unlist))
    i <- i + 1
  }
}
ooo <- do.call("rbind.data.frame", oo)


    

ooo$Var <- factor(ooo$Var, levels = c("low", "med", "hi"))

q<-ggplot(data=ooo, mapping=aes(x=Var, y=log(R+1)))+
  geom_boxplot(data=ooo, aes(color=Corr))+
  scale_x_discrete(limits=c("low", "med","hi"))+
  scale_color_manual(values=c("pink", "darkgreen", "darkblue"))+
  labs(title="contact-alpha")+
  geom_hline(yintercept = log(2), linetype = "dotdash", color = "black", alpha = .8)
q

p<-ggplot(data=ooo, mapping=aes(x=Corr, y=log(R+1)))+
  geom_boxplot(data=ooo, aes(color=Var))+
  scale_x_discrete(limits=c("no", "neg","pos"))+
  scale_color_manual(values=c("pink", "darkgreen", "darkblue"))+
  labs(title="contact-alpha")+
  geom_hline(yintercept = log(2), linetype = "dotdash", color = "black", alpha = .8)
p

ooo %>%
  group_by(Corr,Var)%>%
  summarize(mean=mean(R),var=var(R))-> ooo2

ooo %>%
  group_by(Corr,Var) -> ooo2


mutate(ooo, frac0 = (length(which(R==0)))/(length(R))) 

# contact gamma
## get data frame with all variance level and corr
oo <- vector(mode='list', length=9)
i <- 1
for (var in c("low","med","hi")) {
  for (corr in c("no","neg","pos")) {
    dataName <- paste0("out_",corr,"corr_",var,"var_c_gamma")
    oo[[i]] <- data.frame(R=lapply(get(dataName), function(l) l[[2]]$numInf) %>% unlist,
                          Var=var,
                          Corr=corr)
    i <- i + 1
  }
}
ooo <- do.call("rbind.data.frame", oo)

ooo$Var <- factor(ooo$Var, levels = c("low", "med", "hi"))

q<-ggplot(data=ooo, mapping=aes(x=Var, y=log(R+1)))+
  geom_boxplot(data=ooo, aes(color=Corr))+
  scale_x_discrete(limits=c("low", "med","hi"))+
  scale_color_manual(values=c("pink", "darkgreen", "darkblue"))+
  labs(title="contact-gamma")+
  geom_hline(yintercept = log(2), linetype = "dotdash", color = "black", alpha = .8)
q

p<-ggplot(data=ooo, mapping=aes(x=Corr, y=log(R+1)))+
  geom_boxplot(data=ooo, aes(color=Var))+
  scale_x_discrete(limits=c("no", "neg","pos"))+
  scale_color_manual(values=c("pink", "darkgreen", "darkblue"))+
  labs(title="contact-gamma")+
  geom_hline(yintercept = log(2), linetype = "dotdash", color = "black", alpha = .8)
p

ooo %>%
  group_by(Corr,Var)%>%
  summarize(mean=mean(R),var=var(R))-> ooo2

# shed gamma
## get data frame with all variance level and corr
oo <- vector(mode='list', length=9)
i <- 1
for (var in c("low","med","hi")) {
  for (corr in c("no","neg","pos")) {
    dataName <- paste0("out_",corr,"corr_",var,"var_shed_gamma")
    oo[[i]] <- data.frame(R=lapply(get(dataName), function(l) l[[2]]$numInf) %>% unlist,
                          Var=var,
                          Corr=corr)
    i <- i + 1
  }
}
ooo <- do.call("rbind.data.frame", oo)

ooo$Var <- factor(ooo$Var, levels = c("low", "med", "hi"))

q<-ggplot(data=ooo, mapping=aes(x=Var, y=log(R+1)))+
  geom_boxplot(data=ooo, aes(color=Corr))+
  scale_x_discrete(limits=c("low", "med","hi"))+
  scale_color_manual(values=c("pink", "darkgreen", "darkblue"))+
  labs(title="shed-gamma")+
  geom_hline(yintercept = log(2), linetype = "dotdash", color = "black", alpha = .8)
q

p<-ggplot(data=ooo, mapping=aes(x=Corr, y=log(R+1)))+
  geom_boxplot(data=ooo, aes(color=Var))+
  scale_x_discrete(limits=c("no", "neg","pos"))+
  scale_color_manual(values=c("pink", "darkgreen", "darkblue"))+
  labs(title="shed-gamma")+
  geom_hline(yintercept = log(2), linetype = "dotdash", color = "black", alpha = .8)
p

ooo %>%
  group_by(Corr,Var)%>%
  summarize(mean=mean(R),var=var(R))-> ooo2

# shed alpha
## get data frame with all variance level and corr
oo <- vector(mode='list', length=9)
i <- 1
for (var in c("low","med","hi")) {
  for (corr in c("no","neg","pos")) {
    dataName <- paste0("out_",corr,"corr_",var,"var_shed_alpha")
    oo[[i]] <- data.frame(R=lapply(get(dataName), function(l) l[[2]]$numInf) %>% unlist,
                          Var=var,
                          Corr=corr)
    i <- i + 1
  }
}
ooo <- do.call("rbind.data.frame", oo)

ooo$Var <- factor(ooo$Var, levels = c("low", "med", "hi"))

q<-ggplot(data=ooo, mapping=aes(x=Var, y=log(R+1)))+
  geom_boxplot(data=ooo, aes(color=Corr))+
  scale_x_discrete(limits=c("low", "med","hi"))+
  scale_color_manual(values=c("pink", "darkgreen", "darkblue"))+
  labs(title="shed-alpha")+
  geom_hline(yintercept = log(2), linetype = "dotdash", color = "black", alpha = .8)
q

p<-ggplot(data=ooo, mapping=aes(x=Corr, y=log(R+1)))+
  geom_boxplot(data=ooo, aes(color=Var))+
  scale_x_discrete(limits=c("no", "neg","pos"))+
  scale_color_manual(values=c("pink", "darkgreen", "darkblue"))+
  labs(title="shed-alpha")+
  geom_hline(yintercept = log(2), linetype = "dotdash", color = "black", alpha = .8)
p

ooo %>%
  group_by(Corr,Var)%>%
  summarize(mean=mean(R),var=var(R))-> ooo2

# alpha gamma
## get data frame with all variance level and corr
oo <- vector(mode='list', length=9)
i <- 1
for (var in c("low","med","hi")) {
  for (corr in c("no","neg","pos")) {
    dataName <- paste0("out_",corr,"corr_",var,"var_alpha_gamma")
    oo[[i]] <- data.frame(R=lapply(get(dataName), function(l) l[[2]]$numInf) %>% unlist,
                          Var=var,
                          Corr=corr)
    i <- i + 1
  }
}
ooo <- do.call("rbind.data.frame", oo)

ooo$Var <- factor(ooo$Var, levels = c("low", "med", "hi"))

q<-ggplot(data=ooo, mapping=aes(x=Var, y=log(R+1)))+
  geom_boxplot(data=ooo, aes(color=Corr))+
  scale_x_discrete(limits=c("low", "med","hi"))+
  scale_color_manual(values=c("pink", "darkgreen", "darkblue"))+
  labs(title="alpha-gamma")+
  geom_hline(yintercept = log(2), linetype = "dotdash", color = "black", alpha = .8)
q

p<-ggplot(data=ooo, mapping=aes(x=Corr, y=log(R+1)))+
  geom_boxplot(data=ooo, aes(color=Var))+
  scale_x_discrete(limits=c("no", "neg","pos"))+
  scale_color_manual(values=c("pink", "darkgreen", "darkblue"))+
  labs(title="alpha-gamma")+
  geom_hline(yintercept = log(2), linetype = "dotdash", color = "black", alpha = .8)
p

ooo$Var <- factor(ooo$Var, levels = c("low", "med", "hi"))

ooo %>%
  group_by(Corr,Var)%>%
  summarize(mean=mean(R),var=var(R))-> ooo2


