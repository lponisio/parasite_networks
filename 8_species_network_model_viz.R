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
## Crithidia spp.
## ***********************************************************************

## ***********************************************************************
## 1. Bombus
## ***********************************************************************

## https://www.rensvandeschoot.com/tutorials/generalised-linear-models-with-brms/
## Generate newdata draws

load(file="saved/bombus_CrithidiaPresence.Rdata")
bombus.CrithidiaPresence.cond.effects <-
  conditional_effects(bombus.CrithidiaPresence)

load(file="saved/bombus_ApicystisSpp.Rdata")
bombus.ApicystisSpp.cond.effects <-
  conditional_effects(bombus.ApicystisSpp)

load(file="saved/bombus_SpNosemaBombi.Rdata")
bombus.SpNosemaBombi.cond.effects <-
  conditional_effects(bombus.SpNosemaBombi)

crithidia.bombus.box <- bombus %>%
  ggplot(aes(x = ProjectSubProject, y = SpCrithidiaPresence / SpScreened)) +
  geom_boxplot() +
  ggtitle("Bombus") + 
  theme_ms() +
  geom_jitter(aes(color = GenusSpecies), size = 2, alpha = 0.9) +
  scale_color_viridis(discrete = TRUE) +
  xlab("Project") + 
  ylab("Crithidia prevalence") +
  theme(
    legend.position = "bottom",
    plot.title = element_text(color = "black")  # Ensure title is black
  ) +
  labs(color = "")

apicystis.bombus.box <- bombus %>%
  ggplot(aes(x = ProjectSubProject, y = SpApicystisSpp / SpScreened)) +
  geom_boxplot() +
  ggtitle("Bombus") + 
  theme_ms() +
  geom_jitter(aes(color = GenusSpecies), size = 2, alpha = 0.9) +
  scale_color_viridis(discrete = TRUE) +
  xlab("Project") + 
  ylab("Apicystis prevalence") +
  theme(
    legend.position = "bottom",
    plot.title = element_text(color = "black")  # Ensure title is black
  ) +
  labs(color = "")


## ***************************************************************************
## Crithidia ~ all network metrics

## betweenness 
bombus.CrithidiaPresence_zweighted.betweenness <-
  bombus.CrithidiaPresence.cond.effects[["zweighted.betweenness"]]

p1.parasite <- ggplot(bombus.CrithidiaPresence_zweighted.betweenness,
                      aes(x = zweighted.betweenness, y= estimate__)) +
  geom_line(aes(x = zweighted.betweenness, y= estimate__), size = 1.5) +
  geom_ribbon(aes(ymin = lower__, ymax = upper__, alpha=0.4))+
  labs(x = "Betweenness centrality", y = "Crithidia prevalence",
       fill = "Credible interval") +
  theme_ms() +
  theme(legend.position = "none") +
  geom_point(data=bombus,
             aes(y= SpCrithidiaPresence/SpScreened,
                 x=zweighted.betweenness, shape=ProjectSubProject,
                 color=GenusSpecies)) +
  scale_color_viridis(discrete=TRUE)

## closeness
bombus.CrithidiaPresence_zweighted.closeness <-
  bombus.CrithidiaPresence.cond.effects[["zweighted.closeness"]]

p2.parasite <- ggplot(bombus.CrithidiaPresence_zweighted.closeness,
                      aes(x = zweighted.closeness, y= estimate__)) +
  geom_line(aes(x = zweighted.closeness, y= estimate__), size = 1.5,
            linetype="dashed") +
  geom_ribbon(aes(ymin = lower__, ymax = upper__, alpha=0.4))+
  labs(x = "Closeness centrality", y = "Crithidia prevalence",
       fill = "Credible interval") +
  theme_ms() +
  geom_point(data=bombus,
             aes(y= SpCrithidiaPresence/SpScreened,
                 x=zweighted.betweenness, shape=ProjectSubProject,
                 color=GenusSpecies)) +
  scale_color_viridis(discrete=TRUE)


## degree
bombus.CrithidiaPresence_zdegree <-
  bombus.CrithidiaPresence.cond.effects[["zdegree"]]

p3.parasite <- ggplot(bombus.CrithidiaPresence_zdegree,
                      aes(x = zdegree, y= estimate__)) +
  geom_line(aes(x = zdegree, y= estimate__ ), size = 1.5) +
  geom_ribbon(aes(ymin = lower__, ymax = upper__, alpha=0.4))+
  labs(x = "Degree centrality", y = "Crithidia prevalence",
       fill = "Credible interval") +
  theme_ms() +
  geom_point(data=bombus,
             aes(y= SpCrithidiaPresence/SpScreened,
                 x=zweighted.betweenness, shape=ProjectSubProject,
                 color=GenusSpecies)) +
  scale_color_viridis(discrete=TRUE)

## d
bombus.CrithidiaPresence_zd <-
  bombus.CrithidiaPresence.cond.effects[["zd"]]

p4.parasite <- ggplot(bombus.CrithidiaPresence_zd,
                      aes(x = zd, y= estimate__)) +
  geom_line(aes(x = zd, y= estimate__ ), size = 1.5) +
  geom_ribbon(aes(ymin = lower__, ymax = upper__, alpha=0.4))+
  labs(x = "d", y = "Crithidia prevalence",
       fill = "Credible interval") +
  theme_ms() +
  geom_point(data=bombus,
             aes(y= SpCrithidiaPresence/SpScreened,
                 x=zweighted.betweenness, shape=ProjectSubProject,
                 color=GenusSpecies)) +
  scale_color_viridis(discrete=TRUE)


leg <- p4.parasite + theme(legend.position = "bottom") +
  guides(alpha = "none") + labs(shape = "Project", color="")


bombus.CrithidiaPresence <- ggarrange(p1.parasite, p2.parasite, p3.parasite, p4.parasite,
                                      labels = c("A", "B", "C","D"), 
                                      ncol = 2, nrow = 2,
                                      common.legend = TRUE, legend="bottom",
                                      legend.grob=get_legend(leg))

ggsave(bombus.CrithidiaPresence,
       file="figures/bombus.CrithidiaPresence.pdf",
       height=8, width=12)

## ***********************************************************************
## 2. Apis
## ***********************************************************************

## https://www.rensvandeschoot.com/tutorials/generalised-linear-models-with-brms/
## Generate newdata draws

load(file="saved/apis_CrithidiaPresence.Rdata")
apis.CrithidiaPresence.cond.effects <-
  conditional_effects(apis.CrithidiaPresence)

load(file="saved/apis_ApicystisSpp.Rdata")
apis.ApicystisSpp.cond.effects <-
  conditional_effects(apis.ApicystisSpp)

crithidia.apis.box <- apis %>%
  ggplot(aes(x = ProjectSubProject, y = SpCrithidiaPresence / SpScreened)) +
  geom_boxplot() +
  ggtitle("Apis") + 
  theme_ms() +
  geom_jitter(aes(color = GenusSpecies), size = 2, alpha = 0.9) +
  scale_color_viridis(discrete = TRUE) +
  xlab("Project") + 
  ylab("Crithidia prevalence") +
  theme(
    legend.position = "bottom",
    plot.title = element_text(color = "black")  # Ensure title is black
  ) +
  labs(color = "")


apicystis.apis.box <- apis %>%
  ggplot(aes(x = ProjectSubProject, y = SpApicystisSpp / SpScreened)) +
  geom_boxplot() +
  ggtitle("Apis") + 
  theme_ms() +
  geom_jitter(aes(color = GenusSpecies), size = 2, alpha = 0.9) +
  scale_color_viridis(discrete = TRUE) +
  xlab("Project") + 
  ylab("Apicystis prevalence") +
  theme(
    legend.position = "bottom",
    plot.title = element_text(color = "black")  # Ensure title is black
  ) +
  labs(color = "")


## ***************************************************************************
## Crithidia ~ all network metrics

## betweenness 
apis.CrithidiaPresence_zweighted.betweenness <-
  apis.CrithidiaPresence.cond.effects[["zweighted.betweenness"]]

p1.parasite <- ggplot(apis.CrithidiaPresence_zweighted.betweenness,
                      aes(x = zweighted.betweenness, y= estimate__)) +
  geom_line(aes(x = zweighted.betweenness, y= estimate__), size = 1.5) +
  geom_ribbon(aes(ymin = lower__, ymax = upper__, alpha=0.4))+
  labs(x = "Betweenness centrality", y = "Crithidia prevalence",
       fill = "Credible interval") +
  theme_ms() +
  theme(legend.position = "none") +
  geom_point(data=apis,
             aes(y= GenusCrithidiaPresence/GenusScreened,
                 x=zweighted.betweenness, shape=ProjectSubProject)) +
  scale_color_viridis(discrete=TRUE)

## closeness
apis.CrithidiaPresence_zweighted.closeness <-
  apis.CrithidiaPresence.cond.effects[["zweighted.closeness"]]

p2.parasite <- ggplot(apis.CrithidiaPresence_zweighted.closeness,
                      aes(x = zweighted.closeness, y= estimate__)) +
  geom_line(aes(x = zweighted.closeness, y= estimate__), size = 1.5) +
  geom_ribbon(aes(ymin = lower__, ymax = upper__, alpha=0.4))+
  labs(x = "Closeness centrality", y = "Crithidia prevalence",
       fill = "Credible interval") +
  theme_ms() +
  geom_point(data=apis,
             aes(y= GenusCrithidiaPresence/GenusScreened,
                 x=zweighted.closeness, shape=ProjectSubProject)) +
  scale_color_viridis(discrete=TRUE)


## degree
apis.CrithidiaPresence_zdegree <-
  apis.CrithidiaPresence.cond.effects[["zdegree"]]

p3.parasite <- ggplot(apis.CrithidiaPresence_zdegree,
                      aes(x = zdegree, y= estimate__)) +
  geom_line(aes(x = zdegree, y= estimate__ ), size = 1.5) +
  geom_ribbon(aes(ymin = lower__, ymax = upper__, alpha=0.4))+
  labs(x = "Degree centrality", y = "Crithidia prevalence",
       fill = "Credible interval") +
  theme_ms() +
  geom_point(data=apis,
             aes(y= GenusCrithidiaPresence/GenusScreened,
                 x=zdegree, shape=ProjectSubProject)) +
  scale_color_viridis(discrete=TRUE)

## d
apis.CrithidiaPresence_zd <-
  apis.CrithidiaPresence.cond.effects[["zd"]]

p4.parasite <- ggplot(apis.CrithidiaPresence_zd,
                      aes(x = zd, y= estimate__)) +
  geom_line(aes(x = zd, y= estimate__ ), size = 1.5) +
  geom_ribbon(aes(ymin = lower__, ymax = upper__, alpha=0.4))+
  labs(x = "d", y = "Crithidia prevalence",
       fill = "Credible interval") +
  theme_ms() +
  geom_point(data=apis,
             aes(y= GenusCrithidiaPresence/GenusScreened,
                 x=zd, shape=ProjectSubProject)) +
  scale_color_viridis(discrete=TRUE)


leg <- p4.parasite + theme(legend.position = "bottom") +
  guides(alpha = "none") + labs(shape = "Project", color="")


apis.CrithidiaPresence <- ggarrange(p1.parasite, p2.parasite,
                                    p3.parasite, p4.parasite,
                                    labels = c("A", "B", "C","D"), 
                                    ncol = 2, nrow = 2,
                                    common.legend = TRUE, legend="bottom",
                                    legend.grob=get_legend(leg))


ggsave(apis.CrithidiaPresence,
       file="figures/apis.CrithidiaPresence.pdf",
       height=8, width=12)

## ***********************************************************************
## 2. Melissodes
## ***********************************************************************

## https://www.rensvandeschoot.com/tutorials/generalised-linear-models-with-brms/
## Generate newdata draws

load(file="saved/melissodes_CrithidiaPresence.Rdata")
melissodes.CrithidiaPresence.cond.effects <-
  conditional_effects(melissodes.CrithidiaPresence)

load(file="saved/melissodes_ApicystisSpp.Rdata")
melissodes.ApicystisSpp.cond.effects <-
  conditional_effects(melissodes.ApicystisSpp)


crithidia.melissodes.box <- melissodes %>%
  ggplot(aes(x = ProjectSubProject, y = SpCrithidiaPresence / SpScreened)) +
  geom_boxplot() +
  ggtitle("Melissodes") + 
  theme_ms() +
  geom_jitter(aes(color = GenusSpecies), size = 2, alpha = 0.9) +
  scale_color_viridis(discrete = TRUE) +
  xlab("Project") + 
  ylab("Crithidia prevalence") +
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
  geom_jitter(aes(color = GenusSpecies), size = 2, alpha = 0.9) +
  scale_color_viridis(discrete = TRUE) +
  xlab("Project") + 
  ylab("Apicystis prevalence") +
  theme(
    legend.position = "bottom",
    plot.title = element_text(color = "black")  # Ensure title is black
  ) +
  labs(color = "")


## ***************************************************************************
## Crithidia ~ all network metrics

## betweenness 
melissodes.CrithidiaPresence_zweighted.betweenness <-
  melissodes.CrithidiaPresence.cond.effects[["zweighted.betweenness"]]

p1.parasite <- ggplot(melissodes.CrithidiaPresence_zweighted.betweenness,
                      aes(x = zweighted.betweenness, y= estimate__)) +
  geom_line(aes(x = zweighted.betweenness, y= estimate__), size = 1.5) +
  geom_ribbon(aes(ymin = lower__, ymax = upper__, alpha=0.4))+
  labs(x = "Betweenness centrality", y = "Crithidia prevalence",
       fill = "Credible interval") +
  theme_ms() +
  theme(legend.position = "none") +
  geom_point(data=melissodes,
             aes(y= SpCrithidiaPresence/SpScreened,
                 x=zweighted.betweenness, shape=ProjectSubProject,
                 color=GenusSpecies)) +
  scale_color_viridis(discrete=TRUE)

## closeness
melissodes.CrithidiaPresence_zweighted.closeness <-
  melissodes.CrithidiaPresence.cond.effects[["zweighted.closeness"]]

p2.parasite <- ggplot(melissodes.CrithidiaPresence_zweighted.closeness,
                      aes(x = zweighted.closeness, y= estimate__)) +
  geom_line(aes(x = zweighted.closeness, y= estimate__), size = 1.5) +
  geom_ribbon(aes(ymin = lower__, ymax = upper__, alpha=0.4))+
  labs(x = "Closeness centrality", y = "Crithidia prevalence",
       fill = "Credible interval") +
  theme_ms() +
  geom_point(data=melissodes,
             aes(y= SpCrithidiaPresence/SpScreened,
                 x=zweighted.betweenness, shape=ProjectSubProject,
                 color=GenusSpecies)) +
  scale_color_viridis(discrete=TRUE)


## degree
melissodes.CrithidiaPresence_zdegree <-
  melissodes.CrithidiaPresence.cond.effects[["zdegree"]]

p3.parasite <- ggplot(melissodes.CrithidiaPresence_zdegree,
                      aes(x = zdegree, y= estimate__)) +
  geom_line(aes(x = zdegree, y= estimate__ ), size = 1.5) +
  geom_ribbon(aes(ymin = lower__, ymax = upper__, alpha=0.4))+
  labs(x = "Degree centrality", y = "Crithidia prevalence",
       fill = "Credible interval") +
  theme_ms() +
  geom_point(data=melissodes,
             aes(y= SpCrithidiaPresence/SpScreened,
                 x=zweighted.betweenness, shape=ProjectSubProject,
                 color=GenusSpecies)) +
  scale_color_viridis(discrete=TRUE)

## d
melissodes.CrithidiaPresence_zd <-
  melissodes.CrithidiaPresence.cond.effects[["zd"]]

p4.parasite <- ggplot(melissodes.CrithidiaPresence_zd,
                      aes(x = zd, y= estimate__)) +
  geom_line(aes(x = zd, y= estimate__ ), size = 1.5) +
  geom_ribbon(aes(ymin = lower__, ymax = upper__, alpha=0.4))+
  labs(x = "d", y = "Crithidia prevalence",
       fill = "Credible interval") +
  theme_ms() +
  geom_point(data=melissodes,
             aes(y= SpCrithidiaPresence/SpScreened,
                 x=zweighted.betweenness, shape=ProjectSubProject,
                 color=GenusSpecies)) +
  scale_color_viridis(discrete=TRUE)

leg <- p4.parasite + theme(legend.position = "bottom") +
  guides(alpha = "none") +
  labs(shape = "Project", color="") + labs(shape = "Project", color="")


melissodes.CrithidiaPresence <- ggarrange(p1.parasite, p2.parasite,
                                    p3.parasite, p4.parasite,
                                    labels = c("A", "B", "C","D"), 
                                    ncol = 2, nrow = 2,
                                    common.legend = TRUE, legend="bottom",
                                    legend.grob=get_legend(leg))


ggsave(melissodes.CrithidiaPresence,
       file="figures/melissodes.CrithidiaPresence.pdf",
       height=8, width=12)


crithidia.box <- ggarrange(crithidia.bombus.box, crithidia.apis.box,
                                    crithidia.melissodes.box,
                                    labels = c("A", "B", "C"), 
                                    ncol = 1, nrow = 3)

ggsave(crithidia.box,
       file="figures/box_CrithidiaPresence.pdf",
       height=20, width=10)



apicystis.box <- ggarrange(apicystis.bombus.box, apicystis.apis.box,
                                    apicystis.melissodes.box,
                                    labels = c("A", "B", "C"), 
                                    ncol = 1, nrow = 3)

ggsave(apicystis.box,
       file="figures/box_ApicystisPresence.pdf",
       height=20, width=10)



