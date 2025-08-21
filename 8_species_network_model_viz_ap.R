rm(list=ls())
setwd('~/Dropbox (University of Oregon)/')
setwd("parasite_networks")

## Script for plotting all of the important explanatory variables.
library(ggpubr)
library(gridExtra)
library(tidyverse)
library(brms)
library(tidybayes)
library(viridis)

source("src/misc.R")
source("src/ggplotThemes.R")

## ***********************************************************************
## Apicystis spp.
## ***********************************************************************

## ***********************************************************************
## 1. Bombus
## ***********************************************************************

## https://www.rensvandeschoot.com/tutorials/generalised-linear-models-with-brms/
## Generate newdata draws

load(file="saved/bombus_ApicystisSpp.Rdata")
bombus.ApicystisSpp.cond.effects <-
  conditional_effects(bombus.ApicystisSpp)


sub.bombus <- list(
  SpApicystisSpp=bombus[bombus$ProjectSubProject != "SF" &
                        bombus$ProjectSubProject != "PN-CA-FIRE",], 
  SpApicystisSpp=bombus[bombus$ProjectSubProject != "SF" &
                            bombus$ProjectSubProject != "PN-CA-FIRE",],
  SpNosemaBombi= bombus[bombus$ProjectSubProject != "SI" &
                        bombus$ProjectSubProject != "PN-CA-FIRE",],
  SpNosemaCeranae=bombus[bombus$ProjectSubProject != "SI" &
                         bombus$ProjectSubProject != "PN-CA-FIRE",]
)


apicystis.bombus.box <- bombus %>%
  ggplot(aes(x = ProjectSubProject, y = SpApicystisSpp / SpScreened)) +
  geom_boxplot() +
  ggtitle("Bombus") + 
  theme_ms() +
  geom_jitter(aes(color = GenusSpecies), linewidth = 2, alpha = 0.9) +
  scale_color_viridis(discrete = TRUE) +
  xlab("Project") + 
  ylab("Apicystis prevalence") +
  theme(
    legend.position = "bottom",
    plot.title = element_text(color = "black")  # Ensure title is black
  ) +
  labs(color = "")


## ***************************************************************************
## Apicystis ~ all network metrics

## betweenness 
bombus.ApicystisSpp_zweighted.betweenness <-
  bombus.ApicystisSpp.cond.effects[["zweighted.betweenness"]]

p1.parasite <- ggplot(bombus.ApicystisSpp_zweighted.betweenness,
                      aes(x = zweighted.betweenness, y= estimate__)) +
  geom_line(aes(x = zweighted.betweenness, y= estimate__), linewidth = 1.5) +
  geom_ribbon(aes(ymin = lower__, ymax = upper__, alpha=0.4))+
  labs(x = "Betweenness centrality", y = "Apicystis prevalence",
       fill = "Credible interval") +
  theme_ms() +
  theme(legend.position = "none") +
  geom_point(data=sub.bombus$SpApicystisSpp,
             aes(y= SpApicystisSpp/SpScreened,
                 x=zweighted.betweenness, color=ProjectSubProject
                 )) +
  scale_color_viridis(discrete=TRUE)

## closeness
bombus.ApicystisSpp_zweighted.closeness <-
  bombus.ApicystisSpp.cond.effects[["zweighted.closeness"]]

p2.parasite <- ggplot(bombus.ApicystisSpp_zweighted.closeness,
                      aes(x = zweighted.closeness, y= estimate__)) +
  geom_line(aes(x = zweighted.closeness, y= estimate__), linewidth = 1.5,
            linetype="solid") +
  geom_ribbon(aes(ymin = lower__, ymax = upper__, alpha=0.4))+
  labs(x = "Closeness centrality", y = "Apicystis prevalence",
       fill = "Credible interval") +
  theme_ms() +
  geom_point(data=sub.bombus$SpApicystisSpp,
             aes(y= SpApicystisSpp/SpScreened,
                 x=zweighted.closeness, color=ProjectSubProject
                 )) +
  scale_color_viridis(discrete=TRUE)


## degree
bombus.ApicystisSpp_zdegree <-
  bombus.ApicystisSpp.cond.effects[["zdegree"]]

p3.parasite <- ggplot(bombus.ApicystisSpp_zdegree,
                      aes(x = zdegree, y= estimate__)) +
  geom_line(aes(x = zdegree, y= estimate__ ), linewidth = 1.5) +
  geom_ribbon(aes(ymin = lower__, ymax = upper__, alpha=0.4))+
  labs(x = "Degree centrality", y = "Apicystis prevalence",
       fill = "Credible interval") +
  theme_ms() +
  geom_point(data=sub.bombus$SpApicystisSpp,
             aes(y= SpApicystisSpp/SpScreened,
                 x=zdegree, color=ProjectSubProject
                 )) +
  scale_color_viridis(discrete=TRUE)

## d
bombus.ApicystisSpp_zd <-
  bombus.ApicystisSpp.cond.effects[["zd"]]

p4.parasite <- ggplot(bombus.ApicystisSpp_zd,
                      aes(x = zd, y= estimate__)) +
  geom_line(aes(x = zd, y= estimate__ ), linewidth = 1.5) +
  geom_ribbon(aes(ymin = lower__, ymax = upper__, alpha=0.4))+
  labs(x = "d", y = "Apicystis prevalence",
       fill = "Credible interval") +
  theme_ms() +
  geom_point(data=sub.bombus$SpApicystisSpp,
             aes(y= SpApicystisSpp/SpScreened,
                 x=zd, color=ProjectSubProject
                 )) +
  scale_color_viridis(discrete=TRUE)


## overlap
bombus.ApicystisSpp_zHBOverlap <-
  bombus.ApicystisSpp.cond.effects[["zHBOverlap"]]

p5.parasite <- ggplot(bombus.ApicystisSpp_zHBOverlap,
                      aes(x = zHBOverlap, y= estimate__)) +
  geom_line(aes(x = zHBOverlap, y= estimate__ ), linewidth = 1.5) +
  geom_ribbon(aes(ymin = lower__, ymax = upper__, alpha=0.4))+
  labs(x = "d", y = "Apicystis prevalence",
       fill = "Credible interval") +
  theme_ms() +
  geom_point(data=sub.bombus$SpApicystisSpp,
             aes(y= SpApicystisSpp/SpScreened,
                 x=zHBOverlap, color=ProjectSubProject
                 )) +
  scale_color_viridis(discrete=TRUE)


leg <- p5.parasite + theme(legend.position = "bottom") +
  guides(alpha = "none") + labs(shape = "Project", color="")


bombus.ApicystisSpp <- ggarrange(p1.parasite, p2.parasite,
                                      p3.parasite,
                                      p4.parasite,
                                      p5.parasite,
                                      labels = c("A", "B", "C","D", "E"), 
                                      ncol = 2, nrow = 3,
                                      common.legend = TRUE, legend="bottom",
                                      legend.grob=get_legend(leg))

ggsave(bombus.ApicystisSpp,
       file="figures/bombus.ApicystisSpp.pdf",
       height=10, width=12)

## ***********************************************************************
## 2. Apis
## ***********************************************************************

## https://www.rensvandeschoot.com/tutorials/generalised-linear-models-with-brms/
## Generate newdata draws

load(file="saved/apis_ApicystisSpp.Rdata")
apis.ApicystisSpp.cond.effects <-
  conditional_effects(apis.ApicystisSpp)

apicystis.apis.box <- apis %>%
  ggplot(aes(x = ProjectSubProject, y = SpApicystisSpp / SpScreened)) +
  geom_boxplot() +
  ggtitle("Apis") + 
  theme_ms() +
  geom_jitter(aes(color = GenusSpecies), linewidth = 2, alpha = 0.9) +
  scale_color_viridis(discrete = TRUE) +
  xlab("Project") + 
  ylab("Apicystis prevalence") +
  theme(
    legend.position = "bottom",
    plot.title = element_text(color = "black")  # Ensure title is black
  ) +
  labs(color = "")


## ***************************************************************************
## Apicystis ~ all network metrics

## betweenness 
apis.ApicystisSpp_zweighted.betweenness <-
  apis.ApicystisSpp.cond.effects[["zweighted.betweenness"]]

p1.parasite <- ggplot(apis.ApicystisSpp_zweighted.betweenness,
                      aes(x = zweighted.betweenness, y= estimate__)) +
  geom_line(aes(x = zweighted.betweenness, y= estimate__),
            linewidth = 1.5, linetype="dashed") +
  geom_ribbon(aes(ymin = lower__, ymax = upper__, alpha=0.4))+
  labs(x = "Betweenness centrality", y = "Apicystis prevalence",
       fill = "Credible interval") +
  theme_ms() +
  theme(legend.position = "none") +
  geom_point(data=sub.apis$SpApicystisSpp,
             aes(y= SpApicystisSpp/SpScreened,
                 x=zweighted.betweenness, color=ProjectSubProject
                 )) +
  scale_color_viridis(discrete=TRUE)

## closeness
apis.ApicystisSpp_zweighted.closeness <-
  apis.ApicystisSpp.cond.effects[["zweighted.closeness"]]

p2.parasite <- ggplot(apis.ApicystisSpp_zweighted.closeness,
                      aes(x = zweighted.closeness, y= estimate__)) +
  geom_line(aes(x = zweighted.closeness, y= estimate__), linewidth = 1.5,
            linetype="dashed") +
  geom_ribbon(aes(ymin = lower__, ymax = upper__, alpha=0.4))+
  labs(x = "Closeness centrality", y = "Apicystis prevalence",
       fill = "Credible interval") +
  theme_ms() +
  geom_point(data=sub.apis$SpApicystisSpp,
             aes(y= SpApicystisSpp/SpScreened,
                 x=zweighted.closeness, color=ProjectSubProject
                 )) +
  scale_color_viridis(discrete=TRUE)


## degree
apis.ApicystisSpp_zdegree <-
  apis.ApicystisSpp.cond.effects[["zdegree"]]

p3.parasite <- ggplot(apis.ApicystisSpp_zdegree,
                      aes(x = zdegree, y= estimate__)) +
  geom_line(aes(x = zdegree, y= estimate__ ), linewidth = 1.5,
            linetype="dashed") +
  geom_ribbon(aes(ymin = lower__, ymax = upper__, alpha=0.4))+
  labs(x = "Degree centrality", y = "Apicystis prevalence",
       fill = "Credible interval") +
  theme_ms() +
  geom_point(data=sub.apis$SpApicystisSpp,
             aes(y= SpApicystisSpp/SpScreened,
                 x=zdegree, color=ProjectSubProject
                 )) +
  scale_color_viridis(discrete=TRUE)

## d
apis.ApicystisSpp_zd <-
  apis.ApicystisSpp.cond.effects[["zd"]]

p4.parasite <- ggplot(apis.ApicystisSpp_zd,
                      aes(x = zd, y= estimate__)) +
  geom_line(aes(x = zd, y= estimate__ ), linewidth = 1.5,
            linetype="dashed") +
  geom_ribbon(aes(ymin = lower__, ymax = upper__, alpha=0.4))+
  labs(x = "d", y = "Apicystis prevalence",
       fill = "Credible interval") +
  theme_ms() +
  geom_point(data=sub.apis$SpApicystisSpp,
             aes(y= SpApicystisSpp/SpScreened,
                 x=zd, color=ProjectSubProject
                 )) +
  scale_color_viridis(discrete=TRUE)



leg <- p4.parasite + theme(legend.position = "bottom") +
  guides(alpha = "none") + labs(shape = "Project", color="")


apis.ApicystisSpp <- ggarrange(p1.parasite, p2.parasite,
                                      p3.parasite,
                                      p4.parasite,
                                      labels = c("A", "B", "C","D"), 
                                      ncol = 2, nrow = 2,
                                      common.legend = TRUE, legend="bottom",
                                      legend.grob=get_legend(leg))

ggsave(apis.ApicystisSpp,
       file="figures/apis.ApicystisSpp.pdf",
       height=8, width=12)

## ***********************************************************************
## 2. Melissodes
## ***********************************************************************

## https://www.rensvandeschoot.com/tutorials/generalised-linear-models-with-brms/
## Generate newdata draws

load(file="saved/melissodes_ApicystisSpp.Rdata")
melissodes.ApicystisSpp.cond.effects <-
  conditional_effects(melissodes.ApicystisSpp)

load(file="saved/melissodes_ApicystisSpp.Rdata")
melissodes.ApicystisSpp.cond.effects <-
  conditional_effects(melissodes.ApicystisSpp)


apicystis.melissodes.box <- melissodes %>%
  ggplot(aes(x = ProjectSubProject, y = SpApicystisSpp / SpScreened)) +
  geom_boxplot() +
  ggtitle("Melissodes") + 
  theme_ms() +
  geom_jitter(aes(color = GenusSpecies), linewidth = 2, alpha = 0.9) +
  scale_color_viridis(discrete = TRUE) +
  xlab("Project") + 
  ylab("Apicystis prevalence") +
  theme(
    legend.position = "bottom",
    plot.title = element_text(color = "black")  # Ensure title is black
  ) +
  labs(color = "")

apicystis.melissodes.box <- melissodes %>%
  ggplot(aes(x = ProjectSubProject, y = SpApicystisSpp / SpScreened)) +
  geom_boxplot() +
  ggtitle("Melissodes") + 
  theme_ms() +
  geom_jitter(aes(color = GenusSpecies), linewidth = 2, alpha = 0.9) +
  scale_color_viridis(discrete = TRUE) +
  xlab("Project") + 
  ylab("Apicystis prevalence") +
  theme(
    legend.position = "bottom",
    plot.title = element_text(color = "black")  # Ensure title is black
  ) +
  labs(color = "")


## ***************************************************************************
## Apicystis ~ all network metrics

## betweenness 
melissodes.ApicystisSpp_zweighted.betweenness <-
  melissodes.ApicystisSpp.cond.effects[["zweighted.betweenness"]]

p1.parasite <- ggplot(melissodes.ApicystisSpp_zweighted.betweenness,
                      aes(x = zweighted.betweenness, y= estimate__)) +
  geom_line(aes(x = zweighted.betweenness, y= estimate__), linewidth = 1.5) +
  geom_ribbon(aes(ymin = lower__, ymax = upper__, alpha=0.4))+
  labs(x = "Betweenness centrality", y = "Apicystis prevalence",
       fill = "Credible interval") +
  theme_ms() +
  theme(legend.position = "none") +
  geom_point(data=melissodes,
             aes(y= SpApicystisSpp/SpScreened,
                 x=zweighted.betweenness, color=ProjectSubProject
                 )) +
  scale_color_viridis(discrete=TRUE)

## closeness
melissodes.ApicystisSpp_zweighted.closeness <-
  melissodes.ApicystisSpp.cond.effects[["zweighted.closeness"]]

p2.parasite <- ggplot(melissodes.ApicystisSpp_zweighted.closeness,
                      aes(x = zweighted.closeness, y= estimate__)) +
  geom_line(aes(x = zweighted.closeness, y= estimate__), linewidth = 1.5) +
  geom_ribbon(aes(ymin = lower__, ymax = upper__, alpha=0.4))+
  labs(x = "Closeness centrality", y = "Apicystis prevalence",
       fill = "Credible interval") +
  theme_ms() +
  geom_point(data=melissodes,
             aes(y= SpApicystisSpp/SpScreened,
                 x=zweighted.betweenness, color=ProjectSubProject
                 )) +
  scale_color_viridis(discrete=TRUE)


## degree
melissodes.ApicystisSpp_zdegree <-
  melissodes.ApicystisSpp.cond.effects[["zdegree"]]

p3.parasite <- ggplot(melissodes.ApicystisSpp_zdegree,
                      aes(x = zdegree, y= estimate__)) +
  geom_line(aes(x = zdegree, y= estimate__ ), linewidth = 1.5) +
  geom_ribbon(aes(ymin = lower__, ymax = upper__, alpha=0.4))+
  labs(x = "Degree centrality", y = "Apicystis prevalence",
       fill = "Credible interval") +
  theme_ms() +
  geom_point(data=melissodes,
             aes(y= SpApicystisSpp/SpScreened,
                 x=zweighted.betweenness, color=ProjectSubProject
                 )) +
  scale_color_viridis(discrete=TRUE)

## d
melissodes.ApicystisSpp_zd <-
  melissodes.ApicystisSpp.cond.effects[["zd"]]

p4.parasite <- ggplot(melissodes.ApicystisSpp_zd,
                      aes(x = zd, y= estimate__)) +
  geom_line(aes(x = zd, y= estimate__ ), linewidth = 1.5) +
  geom_ribbon(aes(ymin = lower__, ymax = upper__, alpha=0.4))+
  labs(x = "d", y = "Apicystis prevalence",
       fill = "Credible interval") +
  theme_ms() +
  geom_point(data=melissodes,
             aes(y= SpApicystisSpp/SpScreened,
                 x=zweighted.betweenness, color=ProjectSubProject
                 )) +
  scale_color_viridis(discrete=TRUE)

leg <- p4.parasite + theme(legend.position = "bottom") +
  guides(alpha = "none") +
  labs(shape = "Project", color="") + labs(shape = "Project", color="")


melissodes.ApicystisSpp <- ggarrange(p1.parasite, p2.parasite,
                                    p3.parasite, p4.parasite,
                                    labels = c("A", "B", "C","D"), 
                                    ncol = 2, nrow = 2,
                                    common.legend = TRUE, legend="bottom",
                                    legend.grob=get_legend(leg))


ggsave(melissodes.ApicystisSpp,
       file="figures/melissodes.ApicystisSpp.pdf",
       height=8, width=12)

## all boxes

apicystis.box <- ggarrange(apicystis.bombus.box, apicystis.apis.box,
                                    apicystis.melissodes.box,
                                    labels = c("A", "B", "C"), 
                                    ncol = 1, nrow = 3)

ggsave(apicystis.box,
       file="figures/box_ApicystisSpp.pdf",
       height=20, width=10)




