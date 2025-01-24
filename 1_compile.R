library(stringr)
library(dplyr)
library(ggplot2)
rm(list=ls())
setwd("~/University of Oregon Dropbox/Lauren Ponisio/")
setwd("parasite_networks_saved")
source("../parasite_networks/src/CalcMetrics.R")
source("../parasite_networks/src/SpCalcMetrics.R")
source("../parasite_networks/src/vaznull2.R")

## # Regular expression
## pattern <- "networks/([A-Z]{2})"

## # Use gregexpr to find matches and regmatches to extract them
## matches <- regmatches(net.files, gregexpr(pattern, net.files))

## # Extract the captured groups (the two capital letters)
## extracted <- str_match_all(net.files, pattern)
## extracted.names <- sapply(extracted, function(x) x[,2])

## names(nets) <- extracted.names

## ********************************************************
## Species metrics
## ********************************************************

sp.files <- file.path("sp_metrics", list.files("sp_metrics"))

load(sp.files[[1]])
pnw.sp <- sp.lev
pnw.sp$Project <- "PN"

load(sp.files[[2]])
sf.sp <- sp.lev
sf.sp$Project <- "SF"

load(sp.files[[3]])
si.sp <- sp.lev
si.sp$Project <- "SI"

sp.mets <- rbind(pnw.sp, sf.sp, si.sp)

## ********************************************************
## Parasite prevalence
## ********************************************************

par.files <- file.path("screenings", list.files("screenings"))

load(par.files[[1]])
pnw.par <- all.sums
pnw.par$Project <- "PN"

load(par.files[[3]])
sf.par <- all.sums
sf.par$Project <- "SF"

load(par.files[[2]])
si.par <- all.sums
si.par$Project <- "SI"

colnames(pnw.par)[colnames(pnw.par)=="Stand"] <- "Site"
colnames(pnw.par)[colnames(pnw.par)=="Round"] <- "SampleRound"

si.par$Year <- as.numeric(si.par$Year)
sf.par$Year <- as.numeric(sf.par$Year)
pnw.par$Year <- as.numeric(pnw.par$Year)

si.par$SampleRound <- as.numeric(si.par$SampleRound)
sf.par$SampleRound <- as.numeric(sf.par$SampleRound)
pnw.par$SampleRound <- as.numeric(pnw.par$SampleRound)

par.mets <- rbind(pnw.par, sf.par, si.par)

par.mets <- par.mets[!is.na(par.mets$GenusSpecies),]
par.mets <- par.mets[par.mets$GenusSpecies != "",]

syrphid.genera <- c("Eristalis", "Eristalinus",
                    "Copestylum", "Fazia", "Helophilus",
                    "Paragus", "Syritta", "Syrphus", "Toxomerus",
                    "Eupeodes", "Allograpta", "Melanostoma",
                    "Sphaerophoria")

par.mets <- par.mets[!par.mets$Genus %in% syrphid.genera,]

par.mets <- par.mets[par.mets$SpScreened != 0,]

par.mets$PropGenusCrithidiaPresence <-
  par.mets$GenusCrithidiaPresence/par.mets$GenusScreened

par.mets$PropGenusApicystisSpp <-
  par.mets$GenusApicystisSpp/par.mets$GenusScreened


par.mets$PropSpCrithidiaPresence <-
  par.mets$SpCrithidiaPresence/par.mets$SpScreened

par.mets$PropSpApicystisSpp <-
  par.mets$SpApicystisSpp/par.mets$SpScreened


save(par.mets,
     file="../parasite_networks/data/sp_genus_site_mets.RData")

## ********************************************************
## Networks
## ********************************************************

net.files <- file.path("networks", list.files("networks"))
net.files

load(net.files[[1]])
nets.pnw <- nets

load(net.files[[2]])
nets.sf <- nets

load(net.files[[3]])
nets.si <- nets

N <-  99
all.nets <- c(nets.pnw, nets.sf, nets.si)
all.nets <- all.nets[apply(sapply(all.nets, dim) > 2, 2, all)]

## mets <- lapply(all.nets, calcNetworkMetrics,
##                N=N)

## save(mets,
##      file="../parasite_networks/data/raw_mets.RData")
load(file="../parasite_networks/data/raw_mets.RData")

cols.to.keep <- c("Project", "Year", "SampleRound",
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

save(network.metrics,
     file="../parasite_networks/data/network_mets.RData")

## *****************************************************
## 
## *****************************************************

## N <- 99
## sp.mets <- lapply(all.nets, SpCalcNetworkMetrics,
##                N=N, index=c("closeness", "betweenness",
##                             "degree", "d",
##                             "nestedrank"))
## save(sp.mets,
##      file="../parasite_networks/data/raw_sp_mets.RData")

load(file="../parasite_networks/data/raw_sp_mets.RData")

cols.to.keep <- unique(c("Project", "Year",
                         "SampleRound",
                         "Site",
                  colnames(par.mets)[grep("Sp",
                                          colnames(par.mets))],
                  colnames(par.mets)[grep("Genus",
                                          colnames(par.mets))]
                  ))


sp.network.metrics <- SpPrepDat(sp.mets,  par.mets,
                    cols.to.keep=cols.to.keep,
                    net.type="YearSR")

save(sp.mets,
     file="../parasite_networks/data/sp_mets.RData")
