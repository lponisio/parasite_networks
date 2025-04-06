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

load(file="saved/popgen_bombus_CrithidiaPresence.Rdata")
bombus.CrithidiaPresence.cond.effects <- conditional_effects(bombus.CrithidiaPresence)

load(file="saved/popgen_bombus_ApicystisSpp.Rdata")
bombus.ApicystisSpp.cond.effects <- conditional_effects(bombus.ApicystisSpp)

## ***************************************************************************
## Crithidia ~ all network metrics
##  Ho
bombus.CrithidiaPresence_Ho <-
  bombus.CrithidiaPresence.cond.effects[["Ho"]]

p1.parasite <- ggplot(bombus.CrithidiaPresence_Ho,
                      aes(x = Ho, y= estimate__)) +
  geom_line(aes(x = Ho, y= estimate__), size = 1.5) +
  geom_ribbon(aes(ymin = lower__, ymax = upper__, alpha=0.4))+
  labs(x = "Observed heterozygocity (Ho)", y = "Crithidia prevalence",
       fill = "Credible interval") +
  theme_ms() +
  theme(legend.position = "none") +

  geom_point(data=bombus,
              aes(y= SpCrithidiaPresence/SpScreened,
                  x=Ho, color=GenusSpecies)) +
  scale_color_viridis(discrete=TRUE, option="turbo")

## Fis
bombus.CrithidiaPresence_Fis <-
  bombus.CrithidiaPresence.cond.effects[["Fis"]]

p2.parasite <- ggplot(bombus.CrithidiaPresence_Fis,
                      aes(x = Fis, y= estimate__)) +
  geom_line(aes(x = Fis, y= estimate__), size = 1.5, linetype = "dashed") +
  geom_ribbon(aes(ymin = lower__, ymax = upper__, alpha=0.4))+
  labs(x = "Inbreeding (Fis)", y = "Crithidia prevalence",
       fill = "Credible interval") +
  theme_ms() +

  geom_point(data=bombus,
              aes(y= SpCrithidiaPresence/SpScreened,
                  x=Fis, color=GenusSpecies)) +
  scale_color_viridis(discrete=TRUE, option="turbo")



## Crithidia ~ all network metrics
##  Ho
bombus.ApicystisSpp_Ho <-
  bombus.ApicystisSpp.cond.effects[["Ho"]]

p3.parasite <- ggplot(bombus.ApicystisSpp_Ho,
                      aes(x = Ho, y= estimate__)) +
  geom_line(aes(x = Ho, y= estimate__), size = 1.5, linetype = "dashed") +
  geom_ribbon(aes(ymin = lower__, ymax = upper__, alpha=0.4))+
  labs(x = "Observed heterozygocity (Ho)", y = "Crithidia prevalence",
       fill = "Credible interval") +
  theme_ms() +
  theme(legend.position = "none") +

  geom_point(data=bombus,
              aes(y= SpApicystisSpp/SpScreened,
                  x=Ho, color=GenusSpecies)) +
  scale_color_viridis(discrete=TRUE, option="turbo")

## Fis
bombus.ApicystisSpp_Fis <-
  bombus.ApicystisSpp.cond.effects[["Fis"]]

p4.parasite <- ggplot(bombus.ApicystisSpp_Fis,
                      aes(x = Fis, y= estimate__)) +
  geom_line(aes(x = Fis, y= estimate__), size = 1.5, linetype = "dashed") +
  geom_ribbon(aes(ymin = lower__, ymax = upper__, alpha=0.4))+
  labs(x = "Inbreeding (Fis)", y = "Crithidia prevalence",
       fill = "Credible interval") +
  theme_ms() +

  geom_point(data=bombus,
              aes(y= SpApicystisSpp/SpScreened,
                  x=Fis, color=GenusSpecies)) +
  scale_color_viridis(discrete=TRUE, option="turbo")


leg <- p4.parasite + theme(legend.position = "bottom") + guides(alpha = "none")

    
bombus.popgen <- ggarrange(p1.parasite, p2.parasite, p3.parasite, p4.parasite,
                            labels = c("A", "B", "C","D"), 
                            ncol = 2, nrow = 2,
                            common.legend = TRUE, legend="bottom",
                            legend.grob=get_legend(leg))

ggsave(bombus.popgen,
       file="figures/bombus_popgen.pdf",
       height=8, width=12)
