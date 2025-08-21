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


load(file="saved/network_bombus_ApicystisSpp.Rdata")
bombus.ApicystisSpp.cond.effects <-
  conditional_effects(bombus.ApicystisSpp)

sub.bombus <- list(
  GenusApicystisSpp=bombus[bombus$ProjectSubProject != "SF",], 
  GenusApicystisSpp=bombus[bombus$ProjectSubProject != "SF",],
  GenusNosemaBombi= bombus[bombus$ProjectSubProject != "SI" &
                             bombus$ProjectSubProject != "PN-CA-FIRE",],
  GenusNosemaCeranae=bombus[bombus$ProjectSubProject != "SI" &
                              bombus$ProjectSubProject != "PN-CA-FIRE",]
)

apicystis.bombus.box <- bombus %>%
  ggplot(aes(x = ProjectSubProject,
             y = GenusApicystisSpp / GenusScreened)) +
  geom_boxplot() +
  ggtitle("Bombus") + 
  theme_ms() +
  geom_jitter(aes(color = Genus), linewidth = 2, alpha = 0.9) +
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

## Connectance and nestedness have strong effects

## connectance 
bombus.ApicystisSpp_connectance <-
  bombus.ApicystisSpp.cond.effects[["connectance"]]

p1.parasite <- ggplot(bombus.ApicystisSpp_connectance,
                      aes(x = connectance, y= estimate__)) +
  geom_line(aes(x = connectance, y= estimate__), linewidth = 1.5,
            linetype="solid") +
  geom_ribbon(aes(ymin = lower__, ymax = upper__, alpha=0.4))+
  labs(x = "Connectance", y = "Apicystis prevalence",
       fill = "Credible interval") +
  theme_ms() +
  theme(legend.position = "none") +
  geom_point(data=sub.bombus$GenusApicystisSpp,
             aes(y= PropGenusApicystisSpp,
                 x=connectance, color=ProjectSubProject
                 )) +
  scale_color_viridis(discrete=TRUE)

## Nestedness
bombus.ApicystisSpp_zweighted.NODF <-
  bombus.ApicystisSpp.cond.effects[["zweighted.NODF"]]

p2.parasite <- ggplot(bombus.ApicystisSpp_zweighted.NODF,
                      aes(x = zweighted.NODF, y= estimate__)) +
  geom_line(aes(x = zweighted.NODF, y= estimate__), linewidth = 1.5,
            linetype="solid") +
  geom_ribbon(aes(ymin = lower__, ymax = upper__, alpha=0.4))+
  labs(x = "Nestedness (NODF)", y = "Apicystis prevalence",
       fill = "Credible interval") +
  theme_ms() +
  geom_point(data=sub.bombus$GenusApicystisSpp[
               sub.bombus$GenusApicystisSpp$zweighted.NODF<0.96,],
             aes(y= PropGenusApicystisSpp,
                 x=zweighted.NODF, color=ProjectSubProject
                 )) +
  scale_color_viridis(discrete=TRUE)


## modularity
bombus.ApicystisSpp_zweighted.cluster.coefficient.HL <-
  bombus.ApicystisSpp.cond.effects[["zweighted.cluster.coefficient.HL"]]

p3.parasite <- ggplot(bombus.ApicystisSpp_zweighted.cluster.coefficient.HL,
                      aes(x = zweighted.cluster.coefficient.HL, y= estimate__)) +
  geom_line(aes(x = zweighted.cluster.coefficient.HL, y= estimate__ ),
                      linewidth = 1.5,
            linetype="dashed") +
  geom_ribbon(aes(ymin = lower__, ymax = upper__, alpha=0.4))+
  labs(x = "Modularity (Clustering coeff.)", y = "Apicystis prevalence",
       fill = "Credible interval") +
  theme_ms() +
  geom_point(data=sub.bombus$GenusApicystisSpp[
    sub.bombus$GenusApicystisSpp$zweighted.cluster.coefficient.HL<7.5,]
             aes(y= PropGenusApicystisSpp,
                 x=zweighted.cluster.coefficient.HL,
                 color=ProjectSubProject
                 )) +
  scale_color_viridis(discrete=TRUE)

## H2
bombus.ApicystisSpp_zH2 <-
  bombus.ApicystisSpp.cond.effects[["zH2"]]

p4.parasite <- ggplot(bombus.ApicystisSpp_zH2,
                      aes(x = zH2, y= estimate__)) +
  geom_line(aes(x = zH2, y= estimate__ ), linewidth = 1.5,
            linetype="dashed") +
  geom_ribbon(aes(ymin = lower__, ymax = upper__, alpha=0.4))+
  labs(x = "Specialization (H2)", y = "Apicystis prevalence",
       fill = "Credible interval") +
  theme_ms() +
  geom_point(data=sub.bombus$GenusApicystisSpp,
             aes(y= PropGenusApicystisSpp,
                 x=zH2, color=ProjectSubProject
                 )) +
  scale_color_viridis(discrete=TRUE)


leg <- p4.parasite + theme(legend.position = "bottom") +
  guides(alpha = "none") + labs(shape = "Project", color="")


bombus.ApicystisSpp.plot <- ggarrange(p1.parasite,
                                      p2.parasite, p3.parasite, p4.parasite,
                                      labels = c("A", "B", "C","D"), 
                                      ncol = 2, nrow = 2,
                                      common.legend = TRUE, legend="bottom",
                                      legend.grob=get_legend(leg))

ggsave(bombus.ApicystisSpp.plot,
       file="figures/network_bombus.ApicystisSpp.pdf",
       height=8, width=12)

## ***********************************************************************
## 2. Apis
## ***********************************************************************

## https://www.rensvandeschoot.com/tutorials/generalised-linear-models-with-brms/
## Generate newdata draws

load(file="saved/network_apis_ApicystisSpp.Rdata")
apis.ApicystisSpp.cond.effects <-
  conditional_effects(apis.ApicystisSpp)

apis.sub <- apis[apis$ProjectSubProject != "PN-CA-FIRE" &
                   apis$ProjectSubProject != "PN-COAST",]

apicystis.apis.box <- apis %>%
  ggplot(aes(x = ProjectSubProject, y = GenusApicystisSpp / GenusScreened)) +
  geom_boxplot() +
  ggtitle("Apis") + 
  theme_ms() +
  geom_jitter(aes(color = Genus), linewidth = 2, alpha = 0.9) +
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

## all metrics but clustering coeff, though that is marginal (0.94)

## connectance 
apis.ApicystisSpp_connectance <-
  apis.ApicystisSpp.cond.effects[["connectance"]]

p1.parasite <- ggplot(apis.ApicystisSpp_connectance,
                      aes(x = connectance, y= estimate__)) +
  geom_line(aes(x = connectance, y= estimate__), linewidth = 1.5,
            linetype="solid") +
  geom_ribbon(aes(ymin = lower__, ymax = upper__, alpha=0.4))+
  labs(x = "Connectance", y = "Apicystis prevalence",
       fill = "Credible interval") +
  theme_ms() +
  theme(legend.position = "none") +
  geom_point(data=apis.sub,
             aes(y= PropGenusApicystisSpp,
                 x=connectance, color=ProjectSubProject)) +
  scale_color_viridis(discrete=TRUE)

## Nestedness
apis.ApicystisSpp_zweighted.NODF <-
  apis.ApicystisSpp.cond.effects[["zweighted.NODF"]]

p2.parasite <- ggplot(apis.ApicystisSpp_zweighted.NODF,
                      aes(x = zweighted.NODF, y= estimate__)) +
  geom_line(aes(x = zweighted.NODF, y= estimate__), linewidth = 1.5,
            linetype="solid") +
  geom_ribbon(aes(ymin = lower__, ymax = upper__, alpha=0.4))+
  labs(x = "Nestedness (NODF)", y = "Apicystis prevalence",
       fill = "Credible interval") +
  theme_ms() +
  geom_point(data=apis,
             aes(y= PropGenusApicystisSpp,
                 x=zweighted.NODF, color=ProjectSubProject)) +
  scale_color_viridis(discrete=TRUE)


## modularity
apis.ApicystisSpp_zweighted.cluster.coefficient.HL <-
  apis.ApicystisSpp.cond.effects[["zweighted.cluster.coefficient.HL"]]

p3.parasite <- ggplot(apis.ApicystisSpp_zweighted.cluster.coefficient.HL,
                      aes(x = zweighted.cluster.coefficient.HL, y= estimate__)) +
  geom_line(aes(x = zweighted.cluster.coefficient.HL, y= estimate__ ), linewidth = 1.5,
            linetype="dashed") +
  geom_ribbon(aes(ymin = lower__, ymax = upper__, alpha=0.4))+
  labs(x = "Modularity (Clustering coeff.)", y = "Apicystis prevalence",
       fill = "Credible interval") +
  theme_ms() +
  geom_point(data=apis,
             aes(y= PropGenusApicystisSpp,
                 x=zweighted.cluster.coefficient.HL, color=ProjectSubProject)) +
  scale_color_viridis(discrete=TRUE)

## H2
apis.ApicystisSpp_zH2 <-
  apis.ApicystisSpp.cond.effects[["zH2"]]

p4.parasite <- ggplot(apis.ApicystisSpp_zH2,
                      aes(x = zH2, y= estimate__)) +
  geom_line(aes(x = zH2, y= estimate__ ), linewidth = 1.5,
            linetype="solid") +
  geom_ribbon(aes(ymin = lower__, ymax = upper__, alpha=0.4))+
  labs(x = "Specialization (H2)", y = "Apicystis prevalence",
       fill = "Credible interval") +
  theme_ms() +
  geom_point(data=apis,
             aes(y= PropGenusApicystisSpp,
                 x=zH2, color=ProjectSubProject)) +
  scale_color_viridis(discrete=TRUE)


leg <- p4.parasite + theme(legend.position = "bottom") +
  guides(alpha = "none") + labs(shape = "Project", color="")


apis.ApicystisSpp <- ggarrange(p1.parasite, p2.parasite,
                                    p3.parasite, p4.parasite,
                                    labels = c("A", "B", "C","D"), 
                                    ncol = 2, nrow = 2,
                                    common.legend = TRUE, legend="bottom",
                                    legend.grob=get_legend(leg))


ggsave(apis.ApicystisSpp,
       file="figures/network_apis.ApicystisSpp.pdf",
       height=8, width=12)

## ***********************************************************************
## 2. Melissodes
## ***********************************************************************

## https://www.rensvandeschoot.com/tutorials/generalised-linear-models-with-brms/
## Generate newdata draws

load(file="saved/network_melissodes_ApicystisSpp.Rdata")
melissodes.ApicystisSpp.cond.effects <-
  conditional_effects(melissodes.ApicystisSpp)

apicystis.melissodes.box <- melissodes %>%
  ggplot(aes(x = ProjectSubProject,
             y = GenusApicystisSpp / GenusScreened)) +
  geom_boxplot() +
  ggtitle("Melissodes") + 
  theme_ms() +
  geom_jitter(aes(color = Genus), linewidth = 2, alpha = 0.9) +
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

## no strong effects

## connectance 
melissodes.ApicystisSpp_connectance <-
  melissodes.ApicystisSpp.cond.effects[["connectance"]]

p1.parasite <- ggplot(melissodes.ApicystisSpp_connectance,
                      aes(x = connectance, y= estimate__)) +
  geom_line(aes(x = connectance, y= estimate__), linewidth = 1.5) +
  geom_ribbon(aes(ymin = lower__, ymax = upper__, alpha=0.4))+
  labs(x = "Connectance", y = "Apicystis prevalence",
       fill = "Credible interval") +
  theme_ms() +
  theme(legend.position = "none") +
  geom_point(data=melissodes,
             aes(y= PropGenusApicystisSpp,
                 x=connectance, color=ProjectSubProject
                 )) +
  scale_color_viridis(discrete=TRUE)

## Nestedness
melissodes.ApicystisSpp_zweighted.NODF <-
  melissodes.ApicystisSpp.cond.effects[["zweighted.NODF"]]

p2.parasite <- ggplot(melissodes.ApicystisSpp_zweighted.NODF,
                      aes(x = zweighted.NODF, y= estimate__)) +
  geom_line(aes(x = zweighted.NODF, y= estimate__), linewidth = 1.5) +
  geom_ribbon(aes(ymin = lower__, ymax = upper__, alpha=0.4))+
  labs(x = "Nestedness (NODF)", y = "Apicystis prevalence",
       fill = "Credible interval") +
  theme_ms() +
  geom_point(data=melissodes,
             aes(y= PropGenusApicystisSpp,
                 x=connectance, color=ProjectSubProject
                 )) +
  scale_color_viridis(discrete=TRUE)


## modularity
melissodes.ApicystisSpp_zweighted.cluster.coefficient.HL <-
  melissodes.ApicystisSpp.cond.effects[["zweighted.cluster.coefficient.HL"]]

p3.parasite <- ggplot(melissodes.ApicystisSpp_zweighted.cluster.coefficient.HL,
                      aes(x = zweighted.cluster.coefficient.HL, y= estimate__)) +
  geom_line(aes(x = zweighted.cluster.coefficient.HL, y= estimate__ ), linewidth = 1.5) +
  geom_ribbon(aes(ymin = lower__, ymax = upper__, alpha=0.4))+
  labs(x = "Modularity (Clustering coeff.)", y = "Apicystis prevalence",
       fill = "Credible interval") +
  theme_ms() +
  geom_point(data=melissodes,
             aes(y= PropGenusApicystisSpp,
                 x=connectance, color=ProjectSubProject
                 )) +
  scale_color_viridis(discrete=TRUE)

## H2
melissodes.ApicystisSpp_zH2 <-
  melissodes.ApicystisSpp.cond.effects[["zH2"]]

p4.parasite <- ggplot(melissodes.ApicystisSpp_zH2,
                      aes(x = zH2, y= estimate__)) +
  geom_line(aes(x = zH2, y= estimate__ ), linewidth = 1.5) +
  geom_ribbon(aes(ymin = lower__, ymax = upper__, alpha=0.4))+
  labs(x = "Specialization (H2)", y = "Apicystis prevalence",
       fill = "Credible interval") +
  theme_ms() +
  geom_point(data=melissodes,
             aes(y= PropGenusApicystisSpp,
                 x=connectance, color=ProjectSubProject
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
       file="figures/network_melissodes.ApicystisSpp.pdf",
       height=8, width=12)



## all boxes
apicystis.box <- ggarrange(apicystis.bombus.box, apicystis.apis.box,
                                    apicystis.melissodes.box,
                                    labels = c("A", "B", "C"), 
                                    ncol = 1, nrow = 3)

ggsave(apicystis.box,
       file="figures/box_ApicystisPresence.pdf",
       height=20, width=10)



