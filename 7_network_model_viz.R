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

load(file="saved/network_bombus_CrithidiaPresence.Rdata")
bombus.CrithidiaPresence.cond.effects <-
  conditional_effects(bombus.CrithidiaPresence)

load(file="saved/network_bombus_ApicystisSpp.Rdata")
bombus.ApicystisSpp.cond.effects <-
  conditional_effects(bombus.ApicystisSpp)

load(file="saved/network_bombus_NosemaBombi.Rdata")
bombus.SpNosemaBombi.cond.effects <-
  conditional_effects(bombus.NosemaBombi)

crithidia.bombus.box <- bombus %>%
  ggplot(aes(x = ProjectSubProject,
             y = GenusCrithidiaPresence / GenusScreened)) +
  geom_boxplot() +
  ggtitle("Bombus") + 
  theme_ms() +
  geom_jitter(aes(color = Genus), size = 2, alpha = 0.9) +
  scale_color_viridis(discrete = TRUE) +
  xlab("Project") + 
  ylab("Crithidia prevalence") +
  theme(
    legend.position = "bottom",
    plot.title = element_text(color = "black")  # Ensure title is black
  ) +
  labs(color = "")

apicystis.bombus.box <- bombus %>%
  ggplot(aes(x = ProjectSubProject,
             y = GenusApicystisSpp / GenusScreened)) +
  geom_boxplot() +
  ggtitle("Bombus") + 
  theme_ms() +
  geom_jitter(aes(color = Genus), size = 2, alpha = 0.9) +
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

## connectance 
bombus.CrithidiaPresence_connectance <-
  bombus.CrithidiaPresence.cond.effects[["connectance"]]

p1.parasite <- ggplot(bombus.CrithidiaPresence_connectance,
                      aes(x = connectance, y= estimate__)) +
  geom_line(aes(x = connectance, y= estimate__), size = 1.5,
            linetype="dashed") +
  geom_ribbon(aes(ymin = lower__, ymax = upper__, alpha=0.4))+
  labs(x = "Connectance", y = "Crithidia prevalence",
       fill = "Credible interval") +
  theme_ms() +
  theme(legend.position = "none") +
  geom_point(data=bombus,
             aes(y= GenusCrithidiaPresence/GenusScreened,
                 x=connectance, shape=ProjectSubProject,
                 color=Genus)) +
  scale_color_viridis(discrete=TRUE)

## Nestedness
bombus.CrithidiaPresence_zweighted.NODF <-
  bombus.CrithidiaPresence.cond.effects[["zweighted.NODF"]]

p2.parasite <- ggplot(bombus.CrithidiaPresence_zweighted.NODF,
                      aes(x = zweighted.NODF, y= estimate__)) +
  geom_line(aes(x = zweighted.NODF, y= estimate__), size = 1.5,
            linetype="solid") +
  geom_ribbon(aes(ymin = lower__, ymax = upper__, alpha=0.4))+
  labs(x = "Nestedness", y = "Crithidia prevalence",
       fill = "Credible interval") +
  theme_ms() +
  geom_point(data=bombus,
             aes(y= GenusCrithidiaPresence/GenusScreened,
                 x=zweighted.NODF, shape=ProjectSubProject,
                 color=Genus)) +
  scale_color_viridis(discrete=TRUE)


## redundancy
bombus.CrithidiaPresence_zFunRedundancy.Pols <-
  bombus.CrithidiaPresence.cond.effects[["zFunRedundancy.Pols"]]

p3.parasite <- ggplot(bombus.CrithidiaPresence_zFunRedundancy.Pols,
                      aes(x = zFunRedundancy.Pols, y= estimate__)) +
  geom_line(aes(x = zFunRedundancy.Pols, y= estimate__ ), size = 1.5) +
  geom_ribbon(aes(ymin = lower__, ymax = upper__, alpha=0.4))+
  labs(x = "Redundancy", y = "Crithidia prevalence",
       fill = "Credible interval") +
  theme_ms() +
  geom_point(data=bombus,
             aes(y= GenusCrithidiaPresence/GenusScreened,
                 x=zFunRedundancy.Pols,
                 shape=ProjectSubProject,
                 color=Genus)) +
  scale_color_viridis(discrete=TRUE)

## H2
bombus.CrithidiaPresence_zH2 <-
  bombus.CrithidiaPresence.cond.effects[["zH2"]]

p4.parasite <- ggplot(bombus.CrithidiaPresence_zH2,
                      aes(x = zH2, y= estimate__)) +
  geom_line(aes(x = zH2, y= estimate__ ), size = 1.5) +
  geom_ribbon(aes(ymin = lower__, ymax = upper__, alpha=0.4))+
  labs(x = "H2", y = "Crithidia prevalence",
       fill = "Credible interval") +
  theme_ms() +
  geom_point(data=bombus,
             aes(y= GenusCrithidiaPresence/GenusScreened,
                 x=zH2, shape=ProjectSubProject,
                 color=Genus)) +
  scale_color_viridis(discrete=TRUE)


leg <- p4.parasite + theme(legend.position = "bottom") +
  guides(alpha = "none") + labs(shape = "Project", color="")


bombus.CrithidiaPresence <- ggarrange(p1.parasite,
                                      p2.parasite, p3.parasite, p4.parasite,
                                      labels = c("A", "B", "C","D"), 
                                      ncol = 2, nrow = 2,
                                      common.legend = TRUE, legend="bottom",
                                      legend.grob=get_legend(leg))

ggsave(bombus.CrithidiaPresence,
       file="figures/network_bombus.CrithidiaPresence.pdf",
       height=8, width=12)

## ***********************************************************************
## 2. Apis
## ***********************************************************************

## https://www.rensvandeschoot.com/tutorials/generalised-linear-models-with-brms/
## Generate newdata draws

load(file="saved/network_apis_CrithidiaPresence.Rdata")
apis.CrithidiaPresence.cond.effects <-
  conditional_effects(apis.CrithidiaPresence)

load(file="saved/network_apis_ApicystisSpp.Rdata")
apis.ApicystisSpp.cond.effects <-
  conditional_effects(apis.ApicystisSpp)

crithidia.apis.box <- apis %>%
  ggplot(aes(x = ProjectSubProject,
             y = GenusCrithidiaPresence / GenusScreened)) +
  geom_boxplot() +
  ggtitle("Apis") + 
  theme_ms() +
  geom_jitter(aes(color = Genus), size = 2, alpha = 0.9) +
  scale_color_viridis(discrete = TRUE) +
  xlab("Project") + 
  ylab("Crithidia prevalence") +
  theme(
    legend.position = "bottom",
    plot.title = element_text(color = "black")  # Ensure title is black
  ) +
  labs(color = "")


apicystis.apis.box <- apis %>%
  ggplot(aes(x = ProjectSubProject, y = GenusApicystisSpp / GenusScreened)) +
  geom_boxplot() +
  ggtitle("Apis") + 
  theme_ms() +
  geom_jitter(aes(color = Genus), size = 2, alpha = 0.9) +
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

## connectance 
apis.CrithidiaPresence_connectance <-
  apis.CrithidiaPresence.cond.effects[["connectance"]]

p1.parasite <- ggplot(apis.CrithidiaPresence_connectance,
                      aes(x = connectance, y= estimate__)) +
  geom_line(aes(x = connectance, y= estimate__), size = 1.5,
            linetype="dashed") +
  geom_ribbon(aes(ymin = lower__, ymax = upper__, alpha=0.4))+
  labs(x = "Connectance centrality", y = "Crithidia prevalence",
       fill = "Credible interval") +
  theme_ms() +
  theme(legend.position = "none") +
  geom_point(data=apis,
             aes(y= GenusCrithidiaPresence/GenusScreened,
                 x=connectance, shape=ProjectSubProject)) +
  scale_color_viridis(discrete=TRUE)

## Nestedness
apis.CrithidiaPresence_zweighted.NODF <-
  apis.CrithidiaPresence.cond.effects[["zweighted.NODF"]]

p2.parasite <- ggplot(apis.CrithidiaPresence_zweighted.NODF,
                      aes(x = zweighted.NODF, y= estimate__)) +
  geom_line(aes(x = zweighted.NODF, y= estimate__), size = 1.5,
            linetype="solid") +
  geom_ribbon(aes(ymin = lower__, ymax = upper__, alpha=0.4))+
  labs(x = "Nestedness", y = "Crithidia prevalence",
       fill = "Credible interval") +
  theme_ms() +
  geom_point(data=apis,
             aes(y= GenusCrithidiaPresence/GenusScreened,
                 x=zweighted.NODF, shape=ProjectSubProject)) +
  scale_color_viridis(discrete=TRUE)


## redundancy
apis.CrithidiaPresence_zFunRedundancy.Pols <-
  apis.CrithidiaPresence.cond.effects[["zFunRedundancy.Pols"]]

p3.parasite <- ggplot(apis.CrithidiaPresence_zFunRedundancy.Pols,
                      aes(x = zFunRedundancy.Pols, y= estimate__)) +
  geom_line(aes(x = zFunRedundancy.Pols, y= estimate__ ), size = 1.5,
            linetype="dashed") +
  geom_ribbon(aes(ymin = lower__, ymax = upper__, alpha=0.4))+
  labs(x = "Degree centrality", y = "Crithidia prevalence",
       fill = "Credible interval") +
  theme_ms() +
  geom_point(data=apis,
             aes(y= GenusCrithidiaPresence/GenusScreened,
                 x=zFunRedundancy.Pols, shape=ProjectSubProject)) +
  scale_color_viridis(discrete=TRUE)

## H2
apis.CrithidiaPresence_zH2 <-
  apis.CrithidiaPresence.cond.effects[["zH2"]]

p4.parasite <- ggplot(apis.CrithidiaPresence_zH2,
                      aes(x = zH2, y= estimate__)) +
  geom_line(aes(x = zH2, y= estimate__ ), size = 1.5,
            linetype="dashed") +
  geom_ribbon(aes(ymin = lower__, ymax = upper__, alpha=0.4))+
  labs(x = "H2", y = "Crithidia prevalence",
       fill = "Credible interval") +
  theme_ms() +
  geom_point(data=apis,
             aes(y= GenusCrithidiaPresence/GenusScreened,
                 x=zH2, shape=ProjectSubProject)) +
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
       file="figures/network_apis.CrithidiaPresence.pdf",
       height=8, width=12)

## ***********************************************************************
## 2. Melissodes
## ***********************************************************************

## https://www.rensvandeschoot.com/tutorials/generalised-linear-models-with-brms/
## Generate newdata draws

load(file="saved/network_melissodes_CrithidiaPresence.Rdata")
melissodes.CrithidiaPresence.cond.effects <-
  conditional_effects(melissodes.CrithidiaPresence)

load(file="saved/network_melissodes_ApicystisSpp.Rdata")
melissodes.ApicystisSpp.cond.effects <-
  conditional_effects(melissodes.ApicystisSpp)


crithidia.melissodes.box <- melissodes %>%
  ggplot(aes(x = ProjectSubProject,
             y = GenusCrithidiaPresence / GenusScreened)) +
  geom_boxplot() +
  ggtitle("Melissodes") + 
  theme_ms() +
  geom_jitter(aes(color = Genus), size = 2, alpha = 0.9) +
  scale_color_viridis(discrete = TRUE) +
  xlab("Project") + 
  ylab("Crithidia prevalence") +
  theme(
    legend.position = "bottom",
    plot.title = element_text(color = "black")  # Ensure title is black
  ) +
  labs(color = "")

apicystis.melissodes.box <- melissodes %>%
  ggplot(aes(x = ProjectSubProject,
             y = GenusApicystisSpp / GenusScreened)) +
  geom_boxplot() +
  ggtitle("Melissodes") + 
  theme_ms() +
  geom_jitter(aes(color = Genus), size = 2, alpha = 0.9) +
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

## connectance 
melissodes.CrithidiaPresence_connectance <-
  melissodes.CrithidiaPresence.cond.effects[["connectance"]]

p1.parasite <- ggplot(melissodes.CrithidiaPresence_connectance,
                      aes(x = connectance, y= estimate__)) +
  geom_line(aes(x = connectance, y= estimate__), size = 1.5) +
  geom_ribbon(aes(ymin = lower__, ymax = upper__, alpha=0.4))+
  labs(x = "Connectance", y = "Crithidia prevalence",
       fill = "Credible interval") +
  theme_ms() +
  theme(legend.position = "none") +
  geom_point(data=melissodes,
             aes(y= GenusCrithidiaPresence/GenusScreened,
                 x=connectance, shape=ProjectSubProject,
                 color=Genus)) +
  scale_color_viridis(discrete=TRUE)

## Nestedness
melissodes.CrithidiaPresence_zweighted.NODF <-
  melissodes.CrithidiaPresence.cond.effects[["zweighted.NODF"]]

p2.parasite <- ggplot(melissodes.CrithidiaPresence_zweighted.NODF,
                      aes(x = zweighted.NODF, y= estimate__)) +
  geom_line(aes(x = zweighted.NODF, y= estimate__), size = 1.5) +
  geom_ribbon(aes(ymin = lower__, ymax = upper__, alpha=0.4))+
  labs(x = "Nestedness", y = "Crithidia prevalence",
       fill = "Credible interval") +
  theme_ms() +
  geom_point(data=melissodes,
             aes(y= GenusCrithidiaPresence/GenusScreened,
                 x=connectance, shape=ProjectSubProject,
                 color=Genus)) +
  scale_color_viridis(discrete=TRUE)


## redundancy
melissodes.CrithidiaPresence_zFunRedundancy.Pols <-
  melissodes.CrithidiaPresence.cond.effects[["zFunRedundancy.Pols"]]

p3.parasite <- ggplot(melissodes.CrithidiaPresence_zFunRedundancy.Pols,
                      aes(x = zFunRedundancy.Pols, y= estimate__)) +
  geom_line(aes(x = zFunRedundancy.Pols, y= estimate__ ), size = 1.5) +
  geom_ribbon(aes(ymin = lower__, ymax = upper__, alpha=0.4))+
  labs(x = "Degree centrality", y = "Crithidia prevalence",
       fill = "Credible interval") +
  theme_ms() +
  geom_point(data=melissodes,
             aes(y= GenusCrithidiaPresence/GenusScreened,
                 x=connectance, shape=ProjectSubProject,
                 color=Genus)) +
  scale_color_viridis(discrete=TRUE)

## H2
melissodes.CrithidiaPresence_zH2 <-
  melissodes.CrithidiaPresence.cond.effects[["zH2"]]

p4.parasite <- ggplot(melissodes.CrithidiaPresence_zH2,
                      aes(x = zH2, y= estimate__)) +
  geom_line(aes(x = zH2, y= estimate__ ), size = 1.5) +
  geom_ribbon(aes(ymin = lower__, ymax = upper__, alpha=0.4))+
  labs(x = "H2", y = "Crithidia prevalence",
       fill = "Credible interval") +
  theme_ms() +
  geom_point(data=melissodes,
             aes(y= GenusCrithidiaPresence/GenusScreened,
                 x=connectance, shape=ProjectSubProject,
                 color=Genus)) +
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



