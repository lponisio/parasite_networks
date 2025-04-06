library(stringr)
library(tidyverse)
rm(list=ls())
setwd("~/University of Oregon Dropbox/Lauren Ponisio/")
setwd("parasite_networks_saved")
source("../parasite_networks/src/CalcMetrics.R")
source("../parasite_networks/src/SpCalcMetrics.R")
source("../parasite_networks/src/vaznull2.R")
source("../parasite_networks/src/misc.R")

pnw.subprojects <- read.csv("stands.csv")
## ********************************************************
## Parasite prevalence
## ********************************************************
par.files <- file.path("screenings", list.files("screenings"))
par.files
s
load(par.files[[1]])
hja.par <- all.sums

load(par.files[[2]])
pnw.par <- all.sums

load(par.files[[3]])
sf.par <- all.sums

load(par.files[[4]])
si.par <- all.sums

colnames(pnw.par)[colnames(pnw.par)=="Stand"] <- "Site"
colnames(pnw.par)[colnames(pnw.par)=="Round"] <- "SampleRound"

si.par$Year <- as.numeric(si.par$Year)
sf.par$Year <- as.numeric(sf.par$Year)
pnw.par$Year <- as.numeric(pnw.par$Year)
hja.par$Year <- as.numeric(hja.par$Year)

si.par$SampleRound <- as.numeric(si.par$SampleRound)
sf.par$SampleRound <- as.numeric(sf.par$SampleRound)
pnw.par$SampleRound <- as.numeric(pnw.par$SampleRound)
hja.par$SampleRound <- as.numeric(hja.par$SampleRound)

par.mets <- rbind(pnw.par, sf.par, si.par, hja.par)

par.mets <- par.mets[!is.na(par.mets$GenusSpecies),]
par.mets <- par.mets[par.mets$GenusSpecies != "",]

syrphid.genera <- c("Eristalis", "Eristalinus",
                    "Copestylum", "Fazia", "Helophilus",
                    "Paragus", "Syritta", "Syrphus", "Toxomerus",
                    "Eupeodes", "Allograpta", "Melanostoma",
                    "Sphaerophoria", "Erynnis")

par.mets <- par.mets[!par.mets$Genus %in% syrphid.genera,]

par.mets <- par.mets[par.mets$SpScreened != 0,]

sort(table(par.mets$GenusSpecies))
table(par.mets$GenusSpecies, par.mets$Project)

par.mets$PropGenusCrithidiaPresence <-
  par.mets$GenusCrithidiaPresence/par.mets$GenusScreened

par.mets$PropGenusApicystisSpp <-
  par.mets$GenusApicystisSpp/par.mets$GenusScreened


par.mets$PropSpCrithidiaPresence <-
  par.mets$SpCrithidiaPresence/par.mets$SpScreened

par.mets$PropSpApicystisSpp <-
  par.mets$SpApicystisSpp/par.mets$SpScreened

sampling <- data.frame(
  Project = c("PN", "SI", "SF", "HJ"),
  SurveyMin =c(70, 600, 90, 100)
)

dim(par.mets)
par.mets <- merge(par.mets, sampling, key="Project")
dim(par.mets)

par.mets$SubProject <- ""
par.mets$SubProject <- pnw.subprojects$SubProject[match(par.mets$Site,
                                                        pnw.subprojects$StandCode)]
par.mets$SubProject[is.na(par.mets$SubProject)] <- ""

## Remove bombus nest data
par.mets <- par.mets[par.mets$SubProject != "OR-NEST",]

par.mets$ProjectSubProject <- fix.white.space(paste(par.mets$Project,
                                    par.mets$SubProject, sep="-"))

par.mets$ProjectSubProject[par.mets$SubProject == ""] <-
  par.mets$Project[par.mets$SubProject == ""]

save(par.mets,
     file="../parasite_networks/data/sp_genus_site_mets.RData")

## ********************************************************
## Networks
## ********************************************************

net.files <- file.path("networks", list.files("networks"))
net.files

load(net.files[[1]])
nets.hja <- nets

load(net.files[[2]])
nets.pnw <- nets

load(net.files[[3]])
nets.sf <- nets

load(net.files[[4]])
nets.si <- nets

N <-  99
all.nets <- c(nets.pnw, nets.sf, nets.si, nets.hja)
all.nets <- all.nets[apply(sapply(all.nets, dim) > 2, 2, all)]

## mets <- lapply(all.nets, calcNetworkMetrics,
##                N=N)

## save(mets,
##      file="../parasite_networks/data/raw_mets.RData")

load(file="../parasite_networks/data/raw_mets.RData")

cols.to.keep <- c("Project", "ProjectSubProject", "Year", "SampleRound", "SurveyMin",
                  colnames(par.mets)[grep("Site",
                                          colnames(par.mets))],
                  colnames(par.mets)[grep("Genus",
                                          colnames(par.mets))]
                  )

cols.to.keep <- cols.to.keep[cols.to.keep != "GenusSpecies" &
                             cols.to.keep != "SpSiteYear" ]

network.metrics <- prepDat(mets,  par.mets,
                    cols.to.keep=cols.to.keep,
                    net.type="YearSR")

network.metrics$GenusRelativeAbundance <-
 log(network.metrics$GenusAbundance/network.metrics$SurveyMin)

network.metrics$SiteRelativeBeeAbundance <-
 log(network.metrics$SiteBeeAbundance/network.metrics$SurveyMin)

network.metrics$SiteRelativeBeeDiversity <-
 log(network.metrics$SiteBeeDiversity/network.metrics$SurveyMin)

save(network.metrics,
     file="../parasite_networks/data/network_mets.RData")

## *****************************************************
## Species level metrics
## *****************************************************

N <- 99
## sp.mets <- lapply(all.nets, SpCalcNetworkMetrics,
##                N=N, index=c("closeness", "betweenness",
##                             "degree", "d",
##                             "nestedrank",
##                             "proportional generality"))
## save(sp.mets,
##      file="../parasite_networks/data/raw_sp_mets.RData")

load(file="../parasite_networks/data/raw_sp_mets.RData")

cols.to.keep <- unique(c("Project", "ProjectSubProject", "Year",
                         "SampleRound",
                         "Site",
                         "SiteBeeAbundance",
                         "SiteBeeRichness",
                         "SiteBeeDiversity",
                         "SurveyMin",
                  colnames(par.mets)[grep("Sp",
                                          colnames(par.mets))],
                  colnames(par.mets)[grep("Genus",
                                          colnames(par.mets))]
                  ))


sp.network.metrics <- SpPrepDat(sp.mets,  par.mets,
                    cols.to.keep=cols.to.keep,
                    net.type="YearSR")

sp.network.metrics$GenusRelativeAbundance <-
 log(sp.network.metrics$GenusAbundance/sp.network.metrics$SurveyMin)

sp.network.metrics$SpRelativeAbundance <-
 log(sp.network.metrics$SpAbundance/sp.network.metrics$SurveyMin)

sp.network.metrics$SiteRelativeBeeAbundance <-
 log(sp.network.metrics$SiteBeeAbundance/sp.network.metrics$SurveyMin)

sp.network.metrics$SiteRelativeBeeDiversity <-
 log(sp.network.metrics$SiteBeeDiversity/sp.network.metrics$SurveyMin)



save(sp.network.metrics,
     file="../parasite_networks/data/sp_mets.RData")


## *****************************************************
## Pop gen
## *****************************************************

SI.popgen <- read.csv("popgen/SI.csv")

SI.popgen <- SI.popgen %>%
  pivot_longer(
    cols = starts_with("Site_"),
    names_to = "Site",
    names_prefix = "Site_",
    values_to = "Value",
    values_drop_na = TRUE
  )

SI.popgen

SI.popgen <- SI.popgen %>%
  pivot_wider(names_from = Metric, values_from = Value)
SI.popgen

save(SI.popgen,
     file="../parasite_networks/data/SI_popgen.RData")
