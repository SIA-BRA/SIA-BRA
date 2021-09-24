## SIA-BRA_biomes
# Diniz-Reis et al 2021
# GEB
# GLMM 

rm(list=ls())
graphics.off()
set.seed(1)

### Packages #####
#install.packages("multcomp")
#install.packages("MuMIn")
#install.packages("tidytext")
#install.packages("lmerTest")
#install.packages("evaluate")
#install.packages("ggplot2")
#install.packages("effectsize")
#install.packages("huxtable")
#install.packages('flextable')
#install.packages("officer")
#install.packages("effsize")
library(lmerTest)
library(MuMIn) 
library(multcomp)
library(tidyverse)  
library(broom)
library(tibble)
library(lme4)
library(tidytext)
library(evaluate)
library(broom.mixed)
library(ggplot2)
library(effectsize)
library(huxtable)
library(flextable)
library(officer)
library(effsize)
library(ggdist)
library(cowplot)
library(readr)
library(Rmisc)
library(ggpubr)
library(tidyverse)   
library(systemfonts)   
library(scico)       
library(ggtext)      
library(ggforce)     
library(ggdist)      
library(magick)      
library(patchwork)   
library(kableExtra)  
library(rcartocolor)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(Hmisc)

#### Data upload and exploration ######
choose.dir()
getwd()
data<-read.csv("SIA_BRA_biome.csv", header=T, stringsAsFactors = T) %>%
      filter(tissue !="bone")

options(digits = 3) 

id<-data %>% 
  group_by(class,biome)

stat.c<-summarise(id,
                  avg=mean(d13C_pd, na.rm=TRUE),
                  sd=sd(d13C_pd,na.rm=TRUE),
                  med=median(d13C_pd, na.rm=TRUE),
                  IQR=IQR(d13C_pd, na.rm=TRUE),
                  min=min(d13C_pd, na.rm=TRUE),
                  max=max(d13C_pd, na.rm=TRUE),
                  count=n())
stat.c

stat.n<-summarise(id,
                  avg=mean(d15N, na.rm=TRUE),
                  Sd=sd(d15N,na.rm=TRUE),
                  med=median(d15N, na.rm=TRUE),
                  IQR=IQR(d15N, na.rm=TRUE),
                  min=min(d15N, na.rm=TRUE),
                  max=max(d15N, na.rm=TRUE),
                  count=n())
stat.n

#### Modeling GLMM #####

#Filtering of Fish isotopic data
data_fish<-filter(data, class == "Actinopteri")

#Modelos
data_fish_d13C<-filter(data_fish,d13C_pd!="NA")
modeld13C_fish<-lmer(d13C_pd~ biome + (1|author) + (1|family/genera/species), data=data_fish_d13C)
summary(modeld13C_fish)
r.squaredGLMM(modeld13C_fish)

#Tukey pos-hoc teste para compara??o entre biomas
summary(glht(modeld13C_fish, linfct=mcp(biome="Tukey")))

data_fish_d15N<-filter(data_fish,d15N_pd!="NA")
modeld15N_fish<-lmer(d15N_pd ~ biome + (1|author) + (1|family/genera/species), data_fish_d15N)
r.squaredGLMM(modeld15N_fish)
#Tukey pos-hoc teste para compara??o entre biomas
summary(glht(modeld15N_fish, linfct=mcp(biome="Tukey")))


#Filtering of Aves isotopic data
data_aves<-filter(data, class == "Aves")

#Modelos
data_aves_d13C<-filter(data_aves,d13C_pd!="NA")
modeld13C_aves<-lmer(d13C_pd~ biome + (1|author) + (1|family/genera/species), data=data_aves_d13C)
summary(modeld13C_aves)
r.squaredGLMM(modeld13C_aves)
#Tukey pos-hoc teste para compara??o entre biomas
summary(glht(modeld13C_aves, linfct=mcp(biome="Tukey")))

data_aves_d15N<-filter(data_aves,d15N_pd!="NA")
modeld15N_aves<-lmer(d15N_pd ~ biome + (1|author) + (1|family/genera/species), data=data_aves_d15N)
r.squaredGLMM(modeld15N_aves)
#Tukey pos-hoc teste para compara??o entre biomas
summary(glht(modeld15N_aves, linfct=mcp(biome="Tukey")))

#Filtering of Mammal's isotopic data, the Caatinga biome that have only 6 samples
data_mammals<-filter(data, class == "Mammalia")
data_mammals<-filter(data_mammals, biome != "Caatinga")

data_mammals_d13C<-filter(data_mammals,d13C_pd!="NA")
modeld13C_mammals<-lmer(d13C_pd~ biome + (1|author) + (1|family/genera/species), data=data_mammals_d13C)
summary(modeld13C_mammals)
r.squaredGLMM(modeld13C_mammals)
#Tukey pos-hoc teste para compara??o entre biomas
summary(glht(modeld13C_mammals, linfct=mcp(biome="Tukey")))


data_mammals_d15N<-filter(data_mammals,d15N_pd!="NA")
modeld15N_mammals<-lmer(d15N_pd ~ biome + (1|author) + (1|family/genera/species), data=data_mammals_d15N)
r.squaredGLMM(modeld15N_mammals)
#Tukey pos-hoc teste para compara??o entre biomas
summary(glht(modeld15N_mammals, linfct=mcp(biome="Tukey")))


# Graficos
minhas_cores<- c("#000001","#009E73", "#F0E442", "#D55E00", "#CC79A7", "#56B4E9")

##### Graphics #####
# Raincloud plots with mean and confidence interval
#Aves
summary_data<-summarySE(data_aves, "d13C_pd",  "biome", na.rm = TRUE, conf.interval = 0.95, .drop = TRUE)
head(summary_data)

grafico1<-ggplot(data_aves, aes(biome, d13C_pd, fill = biome)) + 
  stat_halfeye(alpha=0.4,adjust =1,justification = -.1, .width = 0, width = 1,point_colour = NA)+
  geom_point(data = summary_data, aes(x = biome, y = d13C_pd), position = position_nudge(.25), size=2, colour = "black") +
  geom_errorbar(data = summary_data, aes(x = biome, y = d13C_pd, ymin = d13C_pd-ci, ymax = d13C_pd +ci), position = position_nudge(.25), colour = "black", width = 0.1, size = 0.5)+
  geom_boxplot(width = .15, outlier.color = NA, alpha = 0.5, color="black") + 
  annotate("text", x=5, y=-40, label= "a", size = 5)+
  scale_colour_manual(values=c("#000001","#009E73", "#D55E00", "#56B4E9")) +
  scale_fill_manual(values=c("#000001","#009E73", "#D55E00", "#56B4E9"))+
  scale_y_continuous(name=expression(paste(delta^{13}, "Cd (\u2030)")),
                                           breaks = seq(-40,-10,5),
                                           limits=c(-40,-10))+
  theme(axis.line.x = element_line(size = 0.5, colour = "black"),
        axis.line.y = element_line(size = 0.5, colour = "black"),
        axis.line = element_line(size=1, colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.text.x=element_text(size=12, color="black"),
        axis.text.y=element_text(size=12, color="black"),
        axis.title.x=element_text(size=18, face="bold"),
        axis.title.y=element_text(size=16, face="bold"),
        legend.title = element_blank(),
        legend.text = element_blank(),
        legend.position = "none")
grafico1<-grafico1+ coord_flip()
grafico1


summary_data<-summarySE(data_aves, "d15N_pd","biome", na.rm = TRUE, conf.interval = 0.95, .drop = TRUE)
head(summary_data)
grafico2<-ggplot(data_aves, aes(biome, d15N_pd, fill = biome)) + 
  stat_halfeye(alpha=0.4,adjust = 1,justification = -.1, .width = 0, width = 1, point_colour = NA)+
  geom_point(data = summary_data, aes(x = biome, y = d15N_pd), position = position_nudge(.25), size=2, colour = "black") +
  geom_errorbar(data = summary_data, aes(x = biome, y = d15N_pd, ymin = d15N_pd-ci, ymax = d15N_pd +ci), position = position_nudge(.25), colour = "black", width = 0.1, size = 0.5)+
  geom_boxplot(width = .15, outlier.color = NA, alpha = 0.5, color="black") + 
  annotate("text", x=5, y=-5, label= "b", size = 5)+
  scale_colour_manual(values=c("#000001","#009E73", "#D55E00", "#56B4E9")) +
  scale_fill_manual(values=c("#000001","#009E73", "#D55E00", "#56B4E9"))+
  scale_y_continuous(name=expression(paste(delta^{15}, "Nd (\u2030)")),
                                     breaks = seq(-5,20,5),
                                     limits=c(-5, 20))+
  theme(axis.line.x = element_line(size = 0.5, colour = "black"),
        axis.line.y = element_line(size = 0.5, colour = "black"),
        axis.line = element_line(size=1, colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.text.x=element_text(size=12, color="black"),
        axis.text.y=element_text(size=12, color="black"),
        axis.title.x=element_text(size=18, face="bold"),
        axis.title.y=element_text(size=16, face="bold"),
        legend.title = element_blank(),
        legend.text = element_blank(),
        legend.position = "none")

grafico2<-grafico2+ coord_flip()
grafico2

#Fishes
summary_data<-summarySE(data_fish, "d13C_pd",  "biome", na.rm = TRUE, conf.interval = 0.95, .drop = TRUE)
head(summary_data)

grafico3<-ggplot(data_fish, aes(biome, d13C_pd, fill = biome)) + 
  stat_halfeye(alpha=0.4,adjust =1,justification = -.1, .width = 0, width = 1,point_colour = NA)+
  geom_point(data = summary_data, aes(x = biome, y = d13C_pd), position = position_nudge(.25), size=2, colour = "black") +
  geom_errorbar(data = summary_data, aes(x = biome, y = d13C_pd, ymin = d13C_pd-ci, ymax = d13C_pd +ci), position = position_nudge(.25), colour = "black", width = 0.1, size = 0.5)+
  geom_boxplot(width = .15, outlier.color = NA, alpha = 0.5, color="black") + 
  annotate("text", x=6, y=-40, label= "c", size = 5)+
  scale_colour_manual(values=c("#000001","#009E73", "#F0E442", "#D55E00", "#56B4E9")) +
  scale_fill_manual(values=c("#000001","#009E73", "#F0E442", "#D55E00", "#56B4E9"))+
  scale_y_continuous(name=expression(paste(delta^{13}, "Cd (\u2030)")),
                                           breaks = seq(-40,-10,5),
                                           limits=c(-40,-10))+
  theme(axis.line.x = element_line(size = 0.5, colour = "black"),
        axis.line.y = element_line(size = 0.5, colour = "black"),
        axis.line = element_line(size=1, colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.text.x=element_text(size=12, color="black"),
        axis.text.y=element_text(size=12, color="black"),
        axis.title.x=element_text(size=18, face="bold"),
        axis.title.y=element_text(size=16, face="bold"),
        legend.title = element_blank(),
        legend.text = element_blank(),
        legend.position = "none")
grafico3<-grafico3+ coord_flip()
grafico3


summary_data<-summarySE(data_fish, "d15N_pd","biome", na.rm = TRUE, conf.interval = 0.95, .drop = TRUE)
head(summary_data)
grafico4<-ggplot(data_fish, aes(biome, d15N_pd, fill = biome)) + 
  stat_halfeye(alpha=0.4,adjust = 1,justification = -.1, .width = 0, width = 1,point_colour = NA)+
  geom_point(data = summary_data, aes(x = biome, y = d15N_pd), position = position_nudge(.25), size=2, colour = "black") +
  geom_errorbar(data = summary_data, aes(x = biome, y = d15N_pd, ymin = d15N_pd-ci, ymax = d15N_pd +ci), position = position_nudge(.25), colour = "black", width = 0.1, size = 0.5)+
  geom_boxplot(width = .15, outlier.color = NA, alpha = 0.5, color="black") + 
  annotate("text", x=6, y=-5, label= "d", size = 5)+
  scale_colour_manual(values=c("#000001","#009E73", "#F0E442", "#D55E00", "#56B4E9")) +
  scale_fill_manual(values=c("#000001","#009E73", "#F0E442", "#D55E00", "#56B4E9"))+
  scale_y_continuous(name=expression(paste(delta^{15}, "Nd (\u2030)")),
                                     breaks = seq(-5,20,5),
                                     limits=c(-5, 20))+
  theme(axis.line.x = element_line(size = 0.5, colour = "black"),
        axis.line.y = element_line(size = 0.5, colour = "black"),
        axis.line = element_line(size=1, colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.text.x=element_text(size=12, color="black"),
        axis.text.y=element_text(size=12, color="black"),
        axis.title.x=element_text(size=18, face="bold"),
        axis.title.y=element_text(size=16, face="bold"),
        legend.title = element_blank(),
        legend.text = element_blank(),
        legend.position = "none")

grafico4<-grafico4+ coord_flip()
grafico4

#Mammals
summary_data<-summarySE(data_mammals, "d13C_pd",  "biome", na.rm = TRUE, conf.interval = 0.95, .drop = TRUE)
head(summary_data)

#Rainclouds with mean and confidence interval
grafico5<-ggplot(data_mammals, aes(biome, d13C_pd, fill = biome)) + 
  stat_halfeye(alpha=0.4,adjust =1,justification = -.1, .width = 0, width = 1,point_colour = NA)+
  geom_point(data = summary_data, aes(x = biome, y = d13C_pd), position = position_nudge(.25), size=2, colour = "black") +
  geom_errorbar(data = summary_data, aes(x = biome, y = d13C_pd, ymin = d13C_pd-ci, ymax = d13C_pd +ci), position = position_nudge(.25), colour = "black", width = 0.1, size = 0.5)+
  geom_boxplot(width = .15, outlier.color = NA, alpha = 0.5, color="black") + 
  annotate("text", x=5, y=-40, label= "e", size = 5)+ 
  scale_colour_manual(values=c("#000001","#009E73", "#D55E00", "#56B4E9")) +
  scale_fill_manual(values=c("#000001","#009E73", "#D55E00", "#56B4E9"))+
  scale_y_continuous(name=expression(paste(delta^{13}, "Cd (\u2030)")),
                                           breaks = seq(-40,-10,5),
                                           limits=c(-40,-10))+
  theme(axis.line.x = element_line(size = 0.5, colour = "black"),
        axis.line.y = element_line(size = 0.5, colour = "black"),
        axis.line = element_line(size=1, colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.text.x=element_text(size=12, color="black"),
        axis.text.y=element_text(size=12, color="black"),
        axis.title.x=element_text(size=18, face="bold"),
        axis.title.y=element_text(size=16, face="bold"),
        legend.title = element_blank(),
        legend.text = element_blank(),
        legend.position = "none")
grafico5<-grafico5+ coord_flip()
grafico5


summary_data<-summarySE(data_mammals, "d15N_pd","biome", na.rm = TRUE, conf.interval = 0.95, .drop = TRUE)
head(summary_data)

grafico6<-ggplot(data_mammals, aes(biome, d15N_pd, fill = biome)) + 
  stat_halfeye(alpha=0.4,adjust = 1,justification = -.1, .width = 0, width = 1,point_colour = NA)+
  geom_point(data = summary_data, aes(x = biome, y = d15N_pd), position = position_nudge(.25), size=2, colour = "black") +
  geom_errorbar(data = summary_data, aes(x = biome, y = d15N_pd, ymin = d15N_pd-ci, ymax = d15N_pd +ci), position = position_nudge(.25), colour = "black", width = 0.1, size = 0.5)+
  geom_boxplot(width = .15, outlier.color = NA, alpha = 0.5, color="black") + 
  annotate("text", x=5, y=-5, label= "f", size = 5)+
  scale_colour_manual(values=c("#000001","#009E73", "#D55E00", "#56B4E9")) +
  scale_fill_manual(values=c("#000001","#009E73", "#D55E00", "#56B4E9"))+
  scale_y_continuous(name=expression(paste(delta^{15}, "Nd (\u2030)")),
                                     breaks = seq(-5,20,5),
                                     limits=c(-5, 20))+
  theme(axis.line.x = element_line(size = 0.5, colour = "black"),
        axis.line.y = element_line(size = 0.5, colour = "black"),
        axis.line = element_line(size=1, colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.text.x=element_text(size=12, color="black"),
        axis.text.y=element_text(size=12, color="black"),
        axis.title.x=element_text(size=18, face="bold"),
        axis.title.y=element_text(size=16, face="bold"),
        legend.title = element_blank(),
        legend.text = element_blank(),
        legend.position = "none")

grafico6<-grafico6+ coord_flip()
grafico6

graficos<-ggarrange(grafico1+rremove("ylab")+rremove("xlab"),grafico2+rremove("ylab")+rremove("xlab"), 
                    grafico3+rremove("ylab")+rremove("xlab"),grafico4+rremove("ylab")+rremove("xlab"),
                    grafico5+rremove("ylab"),grafico6+rremove("ylab"),
                    common.legend = F, legend = F,  ncol = 2, nrow = 3)
graficos

ggsave("graficos_halfeye.jpeg", units="cm", width=20, height=18, dpi=300)
dev.off()