## super spreading question

## Med R0
## SHED GAMMA
i <- 1
sumdata1 <- vector(mode='list')
for (var in c("low","med","hi")) {
  for (corr in c("no","neg","pos")) {
    data <- readRDS(paste0("out_",corr,"corr_",var,"var_shed_gamma.RDS"))
    lapply(1:length(data), function(d)
      data.frame( Final = mean(tail(data[[d]][[1]]$I, n=50)), # average number of I from last 50 time steps
                  maxInf = max(data[[d]][[2]]$numInf),
                  Corr = as.factor(corr),
                  Var = as.factor(var),
                  Rep = d)) %>% do.call("rbind.data.frame",.) -> sumdata1[[i]]
    i <- i+1
  }
}
sumdata_use1 <-do.call("rbind.data.frame", sumdata1)

shed_gamma<- ggplot(sumdata_use1, aes(x=Final, y=maxInf,color= Var, group=Var))+
  geom_point(data=sumdata_use1, aes(shape=Corr, color = Var), alpha = 0.75, size = 2.5)+
  theme_bw()+
  labs(title="med-R0-shed-gamma")+
  labs(x = "Final I", y = "Max Reff")+
  scale_color_manual(values=c("darkblue", "darkgreen", "pink"))+
  geom_smooth(method='lm', alpha=0.35)+
  expand_limits(x=0, y=0)
shed_gamma

ggplot(sumdata_use1, aes(x=meantop20, y=maxInf,color= Var, group=Var))+
  geom_point(data=sumdata_use1, aes(shape=Corr, color = Var), alpha = 0.75, size = 2.5)+
  theme_bw()+
  labs(title="med-R0-shed-gamma")+
  labs(x = "Final I", y = "Mean # Inf in top 20%")+
  scale_color_manual(values=c("darkblue", "darkgreen", "pink"))


### SHED ALPHA
i <- 1
sumdata1 <- vector(mode='list')
for (var in c("low","med","hi")) {
  for (corr in c("no","neg","pos")) {
    data <- readRDS(paste0("out_",corr,"corr_",var,"var_shed_alpha.RDS"))
    lapply(1:length(data), function(d)
      data.frame(  Final = mean(tail(data[[d]][[1]]$I, n=50)),
                  maxInf = max(data[[d]][[2]]$numInf),
                  Corr = as.factor(corr),
                  Var = as.factor(var),
                  Rep = d)) %>% do.call("rbind.data.frame",.) -> sumdata1[[i]]
    i <- i+1
  }
}
sumdata_use1 <-do.call("rbind.data.frame", sumdata1)

shed_alpha<- ggplot(sumdata_use1, aes(x=Final, y=maxInf,color= Var, group=Var))+
  geom_point(data=sumdata_use1, aes(shape=Corr, color = Var), alpha = 0.75, size = 2.5)+
  theme_bw()+
  labs(title="med-R0-shed-alpha")+
  labs(x = "Final I", y = "Max Reff")+
  scale_color_manual(values=c("darkblue", "darkgreen", "pink"))+
  geom_smooth(method='lm',   alpha=0.35)+
  expand_limits(x=0, y=0)
shed_alpha


### CONTACT SHED
i <- 1
sumdata1 <- vector(mode='list')
for (var in c("low","med","hi")) {
  for (corr in c("no","neg","pos")) {
    data <- readRDS(paste0("out_",corr,"corr_",var,"var_c_shed.RDS"))
    lapply(1:length(data), function(d)
      data.frame(  Final = mean(tail(data[[d]][[1]]$I, n=50)),
                  maxInf = max(data[[d]][[2]]$numInf),
                  Corr = as.factor(corr),
                  Var = as.factor(var),
                  Rep = d)) %>% do.call("rbind.data.frame",.) -> sumdata1[[i]]
    i <- i+1
  }
}
sumdata_use1 <-do.call("rbind.data.frame", sumdata1)

c_shed<- ggplot(sumdata_use1, aes(x=Final, y=maxInf,color= Var, group=Var))+
  geom_point(data=sumdata_use1, aes(shape=Corr, color = Var), alpha = 0.75, size = 2.5)+
  theme_bw()+
  labs(title="med-R0-c-shed")+
  labs(x = "Final I", y = "Max Reff")+
  scale_color_manual(values=c("darkblue", "darkgreen", "pink"))+
  geom_smooth(method='lm',   alpha=0.35)+
  expand_limits(x=0, y=0)
c_shed

### CONTACT GAMMA
i <- 1
sumdata1 <- vector(mode='list')
for (var in c("low","med","hi")) {
  for (corr in c("no","neg","pos")) {
    data <- readRDS(paste0("out_",corr,"corr_",var,"var_c_gamma.RDS"))
    lapply(1:length(data), function(d)
      data.frame(  Final = mean(tail(data[[d]][[1]]$I, n=50)),
                  maxInf = max(data[[d]][[2]]$numInf),
                  Corr = as.factor(corr),
                  Var = as.factor(var),
                  Rep = d)) %>% do.call("rbind.data.frame",.) -> sumdata1[[i]]
    i <- i+1
  }
}
sumdata_use1 <-do.call("rbind.data.frame", sumdata1)

c_gamma<- ggplot(sumdata_use1, aes(x=Final, y=maxInf,color= Var, group=Var))+
  geom_point(data=sumdata_use1, aes(shape=Corr, color = Var), alpha = 0.75, size = 2.5)+
  theme_bw()+
  labs(title="med-R0-c-gamma")+
  labs(x = "Final I", y = "Max Reff")+
  scale_color_manual(values=c("darkblue", "darkgreen", "pink"))+
  geom_smooth(method='lm',   alpha=0.35)+
  expand_limits(x=0, y=0)
c_gamma

### CONTACT ALPHA
i <- 1
sumdata1 <- vector(mode='list')
for (var in c("low","med","hi")) {
  for (corr in c("no","neg","pos")) {
    data <- readRDS(paste0("out_",corr,"corr_",var,"var_c_alpha.RDS"))
    lapply(1:length(data), function(d)
      data.frame(  Final = mean(tail(data[[d]][[1]]$I, n=50)),
                  maxInf = max(data[[d]][[2]]$numInf),
                  Corr = as.factor(corr),
                  Var = as.factor(var),
                  Rep = d)) %>% do.call("rbind.data.frame",.) -> sumdata1[[i]]
    i <- i+1
  }
}
sumdata_use1 <-do.call("rbind.data.frame", sumdata1)

c_alpha<- ggplot(sumdata_use1, aes(x=Final, y=maxInf,color= Var, group=Var))+
  geom_point(data=sumdata_use1, aes(shape=Corr, color = Var), alpha = 0.75, size = 2.5)+
  theme_bw()+
  labs(title="med-R0-c-alpha")+
  labs(x = "Final I", y = "Max Reff")+
  scale_color_manual(values=c("darkblue", "darkgreen", "pink"))+
  geom_smooth(method='lm',   alpha=0.35)+
  expand_limits(x=0, y=0)
c_alpha

## ALPHA GAMMA
i <- 1
sumdata1 <- vector(mode='list')
for (var in c("low","med","hi")) {
  for (corr in c("no","neg","pos")) {
    data <- readRDS(paste0("out_",corr,"corr_",var,"var_alpha_gamma.RDS"))
    lapply(1:length(data), function(d)
      data.frame(  Final = mean(tail(data[[d]][[1]]$I, n=50)),
                  maxInf = max(data[[d]][[2]]$numInf),
                  Corr = as.factor(corr),
                  Var = as.factor(var),
                  Rep = d)) %>% do.call("rbind.data.frame",.) -> sumdata1[[i]]
    i <- i+1
  }
}
sumdata_use1 <-do.call("rbind.data.frame", sumdata1)

alpha_gamma<- ggplot(sumdata_use1, aes(x=Final, y=maxInf,color= Var, group=Var))+
  geom_point(data=sumdata_use1, aes(shape=Corr, color = Var), alpha = 0.75, size = 2.5)+
  theme_bw()+
  labs(title="med-R0-alpha-gamma")+
  labs(x = "Final I", y = "Max Reff")+
  scale_color_manual(values=c("darkblue", "darkgreen", "pink"))+
  geom_smooth(method='lm',   alpha=0.35)+
  expand_limits(x=0, y=0)
alpha_gamma

library(tidyverse)
library(cowplot)
library(ggplot2)
library(ggpubr)

ggarrange(shed_gamma, shed_alpha, c_gamma, c_alpha, c_shed, alpha_gamma, common.legend = TRUE, 
          ncol=2, nrow=3, legend = "right")

## Low R0
## SHED GAMMA
i <- 1
sumdata1 <- vector(mode='list')
for (var in c("low","med","hi")) {
  for (corr in c("no","neg","pos")) {
    data <- readRDS(paste0("out_",corr,"corr_",var,"var_shed_gamma.RDS"))
    lapply(1:length(data), function(d)
      data.frame(  Final = mean(tail(data[[d]][[1]]$I, n=50)),
                  maxInf = max(data[[d]][[2]]$numInf),
                  Corr = as.factor(corr),
                  Var = as.factor(var),
                  Rep = d)) %>% do.call("rbind.data.frame",.) -> sumdata1[[i]]
    i <- i+1
  }
}
sumdata_use1 <-do.call("rbind.data.frame", sumdata1)

shed_gamma<- ggplot(sumdata_use1, aes(x=maxInf, y=Final,color= Var, group=Var))+
  geom_point(data=sumdata_use1, aes(shape=Corr, color = Var), alpha = 0.75, size = 2.5)+
  theme_bw()+
  labs(title="Med-R0-shed-gamma")+
  labs(x = "Max Reff", y = "Final I")+
  scale_color_manual(values=c("darkblue", "darkgreen", "pink"))+
  ylim(0,175)+
  xlim(0,425)+
  theme(text = element_text(size=15),
        plot.title = element_text(hjust = 0.5),
        plot.margin = margin(t = 30,  # Top margin
                             r = 20,  # Right margin
                             b = 5,  # Bottom margin
                             l = 40))
shed_gamma


### SHED ALPHA
i <- 1
sumdata1 <- vector(mode='list')
for (var in c("low","med","hi")) {
  for (corr in c("no","neg","pos")) {
    data <- readRDS(paste0("out_",corr,"corr_",var,"var_shed_alpha.RDS"))
    lapply(1:length(data), function(d)
      data.frame(  Final = mean(tail(data[[d]][[1]]$I, n=50)),
                  maxInf = max(data[[d]][[2]]$numInf),
                  Corr = as.factor(corr),
                  Var = as.factor(var),
                  Rep = d)) %>% do.call("rbind.data.frame",.) -> sumdata1[[i]]
    i <- i+1
  }
}
sumdata_use1 <-do.call("rbind.data.frame", sumdata1)

shed_alpha<- ggplot(sumdata_use1, aes(x=maxInf, y=Final,color= Var, group=Var))+
  geom_point(data=sumdata_use1, aes(shape=Corr, color = Var), alpha = 0.75, size = 2.5)+
  theme_bw()+
  labs(title="Med-R0-shed-alpha")+
  labs(x = "Max Reff", y = "Final I")+
  scale_color_manual(values=c("darkblue", "darkgreen", "pink"))+
  ylim(0,175)+
  xlim(0,425)+
  theme(text = element_text(size=15),
        plot.title = element_text(hjust = 0.5),
        plot.margin = margin(t = 20,  # Top margin
                             r = 20,  # Right margin
                             b = 5,  # Bottom margin
                             l = 40))
shed_alpha


### CONTACT SHED
i <- 1
sumdata1 <- vector(mode='list')
for (var in c("low","med","hi")) {
  for (corr in c("no","neg","pos")) {
    data <- readRDS(paste0("out_",corr,"corr_",var,"var_c_shed.RDS"))
    lapply(1:length(data), function(d)
      data.frame(  Final = mean(tail(data[[d]][[1]]$I, n=50)),
                  maxInf = max(data[[d]][[2]]$numInf),
                  Corr = as.factor(corr),
                  Var = as.factor(var),
                  Rep = d)) %>% do.call("rbind.data.frame",.) -> sumdata1[[i]]
    i <- i+1
  }
}
sumdata_use1 <-do.call("rbind.data.frame", sumdata1)

c_shed<- ggplot(sumdata_use1, aes(x=maxInf, y=Final,color= Var, group=Var))+
  geom_point(data=sumdata_use1, aes(shape=Corr, color = Var), alpha = 0.75, size = 2.5)+
  theme_bw()+
  labs(title="Med-R0-c-shed")+
  labs(x = "Max Reff", y = "Final I")+
  scale_color_manual(values=c("darkblue", "darkgreen", "pink"))+
  ylim(0,175)+
  xlim(0,425)+
  theme(text = element_text(size=15),
        plot.title = element_text(hjust = 0.5),
        plot.margin = margin(t = 20,  # Top margin
                             r = 20,  # Right margin
                             b = 5,  # Bottom margin
                             l = 40))
c_shed

### CONTACT GAMMA
i <- 1
sumdata1 <- vector(mode='list')
for (var in c("low","med","hi")) {
  for (corr in c("no","neg","pos")) {
    data <- readRDS(paste0("out_",corr,"corr_",var,"var_c_gamma.RDS"))
    lapply(1:length(data), function(d)
      data.frame(  Final = mean(tail(data[[d]][[1]]$I, n=50)),
                  maxInf = max(data[[d]][[2]]$numInf),
                  Corr = as.factor(corr),
                  Var = as.factor(var),
                  Rep = d)) %>% do.call("rbind.data.frame",.) -> sumdata1[[i]]
    i <- i+1
  }
}
sumdata_use1 <-do.call("rbind.data.frame", sumdata1)

c_gamma<- ggplot(sumdata_use1, aes(x=maxInf, y=Final,color= Var, group=Var))+
  geom_point(data=sumdata_use1, aes(shape=Corr, color = Var), alpha = 0.75, size = 2.5)+
  theme_bw()+
  labs(title="Med-R0-c-gamma")+
  labs(x = "Max Reff", y = "Final I")+
  scale_color_manual(values=c("darkblue", "darkgreen", "pink"))+
  ylim(0,175)+
  xlim(0,425)+
  theme(text = element_text(size=15),
        plot.title = element_text(hjust = 0.5),
        plot.margin = margin(t = 20,  # Top margin
                             r = 20,  # Right margin
                             b = 5,  # Bottom margin
                             l = 40))
c_gamma

### CONTACT ALPHA
i <- 1
sumdata1 <- vector(mode='list')
for (var in c("low","med","hi")) {
  for (corr in c("no","neg","pos")) {
    data <- readRDS(paste0("out_",corr,"corr_",var,"var_c_alpha.RDS"))
    lapply(1:length(data), function(d)
      data.frame(  Final = mean(tail(data[[d]][[1]]$I, n=50)),
                  maxInf = max(data[[d]][[2]]$numInf),
                  Corr = as.factor(corr),
                  Var = as.factor(var),
                  Rep = d)) %>% do.call("rbind.data.frame",.) -> sumdata1[[i]]
    i <- i+1
  }
}
sumdata_use1 <-do.call("rbind.data.frame", sumdata1)

c_alpha<- ggplot(sumdata_use1, aes(x=maxInf, y=Final,color= Var, group=Var))+
  geom_point(data=sumdata_use1, aes(shape=Corr, color = Var), alpha = 0.75, size = 2.5)+
  theme_bw()+
  labs(title="Med-R0-c-alpha")+
  labs(x = "Max Reff", y = "Final I")+
  scale_color_manual(values=c("darkblue", "darkgreen", "pink"))+
  ylim(0,175)+
  xlim(0,425)+
  theme(text = element_text(size=15),
        plot.title = element_text(hjust = 0.5),
        plot.margin = margin(t = 20,  # Top margin
                             r = 20,  # Right margin
                             b = 5,  # Bottom margin
                             l = 40))
c_alpha

## ALPHA GAMMA
i <- 1
sumdata1 <- vector(mode='list')
for (var in c("low","med","hi")) {
  for (corr in c("no","neg","pos")) {
    data <- readRDS(paste0("out_",corr,"corr_",var,"var_alpha_gamma.RDS"))
    lapply(1:length(data), function(d)
      data.frame(  Final = mean(tail(data[[d]][[1]]$I, n=50)),
                  maxInf = max(data[[d]][[2]]$numInf),
                  Corr = as.factor(corr),
                  Var = as.factor(var),
                  Rep = d)) %>% do.call("rbind.data.frame",.) -> sumdata1[[i]]
    i <- i+1
  }
}
sumdata_use1 <-do.call("rbind.data.frame", sumdata1)

alpha_gamma<- ggplot(sumdata_use1, aes(x=maxInf, y=Final,color= Var, group=Var))+
  geom_point(data=sumdata_use1, aes(shape=Corr, color = Var), alpha = 0.75, size = 2.5)+
  theme_bw()+
  labs(title="Med-R0-alpha-gamma")+
  labs(x = "Max Reff", y = "Final I")+
  scale_color_manual(values=c("darkblue", "darkgreen", "pink"))+
  ylim(0,175)+
  xlim(0,425)+
  theme(text = element_text(size=15),
        plot.title = element_text(hjust = 0.5),
        plot.margin = margin(t = 20,  # Top margin
                             r = 20,  # Right margin
                             b = 5,  # Bottom margin
                             l = 40))
alpha_gamma

library(tidyverse)
library(cowplot)
library(ggplot2)
library(ggpubr)

ggarrange(shed_gamma, shed_alpha, c_gamma, c_alpha, c_shed, alpha_gamma, common.legend = TRUE, 
          ncol=2, nrow=3, legend = "right")

### HIHGHHGHHHGHGHGHGHG R NAUGHT

## SHED GAMMA
i <- 1
sumdata1 <- vector(mode='list')
for (var in c("low","med","hi")) {
  for (corr in c("no","neg","pos")) {
    data <- readRDS(paste0("out_",corr,"corr_",var,"var_shed_gamma.RDS"))
    lapply(1:length(data), function(d)
      data.frame( Final = mean(tail(data[[d]][[1]]$I, n=50)),
                  maxInf = max(data[[d]][[2]]$numInf),
                  Corr = as.factor(corr),
                  Var = as.factor(var),
                  Rep = d)) %>% do.call("rbind.data.frame",.) -> sumdata1[[i]]
    i <- i+1
  }
}
sumdata_use1 <-do.call("rbind.data.frame", sumdata1)

shed_gamma<- ggplot(sumdata_use1, aes(x=maxInf, y=Final,color= Var, group=Var))+
  geom_point(data=sumdata_use1, aes(shape=Corr, color = Var), alpha = 0.75, size = 2.5)+
  theme_bw()+
  labs(title="hi-R0-shed-gamma")+
  labs(x = "Max Reff", y = "Final I")+
  scale_color_manual(values=c("darkblue", "darkgreen", "pink"))+
  ylim(0,175)+
  xlim(0,425)
shed_gamma

shed_gamma<- ggplot(sumdata_use1, aes(x=maxInf, y=Final,color= Var, group=Var))+
  geom_boxplot(data=sumdata_use1, aes(shape=Corr, color = Var))+
  theme_bw()+
  labs(title="hi-R0-shed-gamma")+
  labs(x = "Max Reff", y = "Final I")+
  scale_color_manual(values=c("darkblue", "darkgreen", "pink"))+
  ylim(0,150)+
  xlim(0,350)


### SHED ALPHA
i <- 1
sumdata1 <- vector(mode='list')
for (var in c("low","med","hi")) {
  for (corr in c("no","neg","pos")) {
    data <- readRDS(paste0("out_",corr,"corr_",var,"var_shed_alpha.RDS"))
    lapply(1:length(data), function(d)
      data.frame(  Final = mean(tail(data[[d]][[1]]$I, n=50)),
                  maxInf = max(data[[d]][[2]]$numInf),
                  Corr = as.factor(corr),
                  Var = as.factor(var),
                  Rep = d)) %>% do.call("rbind.data.frame",.) -> sumdata1[[i]]
    i <- i+1
  }
}
sumdata_use1 <-do.call("rbind.data.frame", sumdata1)

shed_alpha<- ggplot(sumdata_use1, aes(x=maxInf, y=Final,color= Var, group=Var))+
  geom_point(data=sumdata_use1, aes(shape=Corr, color = Var), alpha = 0.75, size = 2.5)+
  theme_bw()+
  labs(title="hi-R0-shed-alpha")+
  labs(x = "Max Reff", y = "Final I")+
  scale_color_manual(values=c("darkblue", "darkgreen", "pink"))+
  ylim(0,175)+
  xlim(0,425)
shed_alpha


### CONTACT SHED
i <- 1
sumdata1 <- vector(mode='list')
for (var in c("low","med","hi")) {
  for (corr in c("no","neg","pos")) {
    data <- readRDS(paste0("out_",corr,"corr_",var,"var_c_shed.RDS"))
    lapply(1:length(data), function(d)
      data.frame(  Final = mean(tail(data[[d]][[1]]$I, n=50)),
                  maxInf = max(data[[d]][[2]]$numInf),
                  Corr = as.factor(corr),
                  Var = as.factor(var),
                  Rep = d)) %>% do.call("rbind.data.frame",.) -> sumdata1[[i]]
    i <- i+1
  }
}
sumdata_use1 <-do.call("rbind.data.frame", sumdata1)

c_shed<- ggplot(sumdata_use1, aes(x=log10(maxInf), y=log10(Final),color= Var, group=Var))+
  geom_point(data=sumdata_use1, aes(shape=Corr, color = Var), alpha = 0.75, size = 2.5)+
  theme_bw()+
  labs(title="hi-R0-c-shed")+
  labs(x = "Max Reff", y = "Final I")+
  scale_color_manual(values=c("darkblue", "darkgreen", "pink"))

c_shed<- ggplot(sumdata_use1, aes(x=maxInf, y=Final,color= Var, group=Var))+
  geom_boxplot(data=sumdata_use1, aes(color = Var))+
  theme_bw()+
  labs(title="hi-R0-c-shed")+
  labs(x = "Max Reff", y = "Final I")+
  scale_color_manual(values=c("darkblue", "darkgreen", "pink"))+
  ylim(0,150)+
  xlim(0,350)

ggarrange(shed_gamma, c_shed, common.legend = TRUE, nrow=1, ncol=2)

ggplot(sumdata_use1, aes(x=maxInf, y=Final,color= Var, group=Var))+
  geom_violin(data=sumdata_use1, aes(color = Var))+
  theme_bw()+
  labs(title="hi-R0-c-shed")+
  labs(x = "Max Reff", y = "Final I")+
  scale_color_manual(values=c("darkblue", "darkgreen", "pink"))+
  ylim(0,150)+
  xlim(0,350)

c_shed

### CONTACT GAMMA
i <- 1
sumdata1 <- vector(mode='list')
for (var in c("low","med","hi")) {
  for (corr in c("no","neg","pos")) {
    data <- readRDS(paste0("out_",corr,"corr_",var,"var_c_gamma.RDS"))
    lapply(1:length(data), function(d)
      data.frame(  Final = mean(tail(data[[d]][[1]]$I, n=50)),
                  maxInf = max(data[[d]][[2]]$numInf),
                  Corr = as.factor(corr),
                  Var = as.factor(var),
                  Rep = d)) %>% do.call("rbind.data.frame",.) -> sumdata1[[i]]
    i <- i+1
  }
}
sumdata_use1 <-do.call("rbind.data.frame", sumdata1)

c_gamma<- ggplot(sumdata_use1, aes(x=maxInf, y=Final,color= Var, group=Var))+
  geom_point(data=sumdata_use1, aes(shape=Corr, color = Var), alpha = 0.75, size = 2.5)+
  theme_bw()+
  labs(title="hi-R0-c-gamma")+
  labs(x = "Max Reff", y = "Final I")+
  scale_color_manual(values=c("darkblue", "darkgreen", "pink"))+
  ylim(0,175)+
  xlim(0,425)
c_gamma

### CONTACT ALPHA
i <- 1
sumdata1 <- vector(mode='list')
for (var in c("low","med","hi")) {
  for (corr in c("no","neg","pos")) {
    data <- readRDS(paste0("out_",corr,"corr_",var,"var_c_alpha.RDS"))
    lapply(1:length(data), function(d)
      data.frame(  Final = mean(tail(data[[d]][[1]]$I, n=50)),
                  maxInf = max(data[[d]][[2]]$numInf),
                  Corr = as.factor(corr),
                  Var = as.factor(var),
                  Rep = d)) %>% do.call("rbind.data.frame",.) -> sumdata1[[i]]
    i <- i+1
  }
}
sumdata_use1 <-do.call("rbind.data.frame", sumdata1)

c_alpha<- ggplot(sumdata_use1, aes(x=maxInf, y=Final,color= Var, group=Var))+
  geom_point(data=sumdata_use1, aes(shape=Corr, color = Var), alpha = 0.75, size = 2.5)+
  theme_bw()+
  labs(title="hi-R0-c-alpha")+
  labs(x = "Max Reff", y = "Final I")+
  scale_color_manual(values=c("darkblue", "darkgreen", "pink"))+
  ylim(0,175)+
  xlim(0,425)
c_alpha

## ALPHA GAMMA
i <- 1
sumdata1 <- vector(mode='list')
for (var in c("low","med","hi")) {
  for (corr in c("no","neg","pos")) {
    data <- readRDS(paste0("out_",corr,"corr_",var,"var_alpha_gamma.RDS"))
    lapply(1:length(data), function(d)
      data.frame(  Final = mean(tail(data[[d]][[1]]$I, n=50)),
                  maxInf = max(data[[d]][[2]]$numInf),
                  Corr = as.factor(corr),
                  Var = as.factor(var),
                  Rep = d)) %>% do.call("rbind.data.frame",.) -> sumdata1[[i]]
    i <- i+1
  }
}
sumdata_use1 <-do.call("rbind.data.frame", sumdata1)

alpha_gamma<- ggplot(sumdata_use1, aes(x=maxInf, y=Final,color= Var, group=Var))+
  geom_point(data=sumdata_use1, aes(shape=Corr, color = Var), alpha = 0.75, size = 2.5)+
  theme_bw()+
  labs(title="hi-R0-alpha-gamma")+
  labs(x = "Max Reff", y = "Final I")+
  scale_color_manual(values=c("darkblue", "darkgreen", "pink"))+
  ylim(0,175)+
  xlim(0,425)
alpha_gamma

library(tidyverse)
library(cowplot)
library(ggplot2)
library(ggpubr)

ggarrange(shed_gamma, shed_alpha, c_gamma, c_alpha, c_shed, alpha_gamma, common.legend = TRUE, 
          ncol=2, nrow=3, legend = "right")
