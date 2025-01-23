library(stringr)
library(dplyr)
library(ggplot2)
rm(list=ls())
setwd("~/University of Oregon Dropbox/Lauren Ponisio/")
setwd("parasite_networks_saved")
source("../parasite_networks/src/CalcMetrics.R")
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
sp.mets$Year <- as.numeric(sp.mets$Year)

sp.mets$SampleRound <- as.numeric(sp.mets$SampleRound)

poll.mets <- sp.mets[sp.mets$speciesType == "higher.level",]

poll.mets$Genus <- sapply(strsplit(poll.mets$GenusSpecies,
                          "\\s"),
                          function(x) x[1])

poll.mets$weighted.closeness[poll.mets$Genus == "Bombus" &
                             poll.mets$Project == "SI"]

poll.mets$weighted.closeness[poll.mets$Genus == "Bombus" &
                             poll.mets$Project == "PN"]

key <- colnames(par.mets)[colnames(par.mets) %in% colnames(poll.mets)]

par.mets$key <- apply(par.mets[ , key ], 1, paste,
                      collapse = "-" )
poll.mets$key <- apply(poll.mets[ , key ], 1, paste,
                      collapse = "-" )

si.keys <- poll.mets$key[poll.mets$Project == "SI"]
si.par.keys <- par.mets$key[par.mets$Project == "SI"]

si.keys[si.keys %in% si.par.keys]

dim(par.mets)
all.mets <- left_join(par.mets, poll.mets, by=key)
dim(all.mets)



all.mets$PropGenusCrithidiaPresence <-
  all.mets$GenusCrithidiaPresence/all.mets$GenusScreened

## all.mets$PropSpCrithidiaPresence <-
##   all.mets$SpCrithidiaPresence/all.mets$SpScreened


## ggplot(all.mets, aes(x = weighted.closeness,
##                    y = PropSpCrithidiaPresence,
##                     color = Project)) +
##   geom_point(size=2)


## bombus <- all.mets[all.mets$Genus == "Bombus",]

## ggplot(bombus, aes(x = weighted.closeness,
##                    y = PropSpCrithidiaPresence,
##                     color = Project)) +
##   geom_point(size=2)



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

mets <- lapply(all.nets, calcNetworkMetrics,
               N=N)

cols.to.keep <- c("Project", "Year", "SampleRound",
                  colnames(all.mets)[grep("Site", colnames(all.mets))],
                  colnames(all.mets)[grep("Genus", colnames(all.mets))]
                  )

cols.to.keep <- cols.to.keep[cols.to.keep != "GenusSpecies" &
                             cols.to.keep != "SpSiteYear" ]

network.metrics <- prepDat(mets,  all.mets,
                    cols.to.keep=cols.to.keep,
                    net.type="YearSR")

save(network.metrics,
     file="../parasite_networks/data/network_mets.RData")

bombus <- cor.dats[cor.dats$Genus == "Bombus",]

nodf <- ggplot(bombus, aes(x = zweighted.NODF,
                   y = PropGenusCrithidiaPresence,
                    color = Project)) +
  geom_point(size=2)


niche.overlap.hl <- ggplot(bombus,
                        aes(x = zniche.overlap.HL,
                            y = PropGenusCrithidiaPresence,
                            color = Project)) +
  geom_point(size=2)


niche.overlap.ll <- ggplot(bombus,
                        aes(x = zniche.overlap.LL,
                            y = PropGenusCrithidiaPresence,
                            color = Project)) +
  geom_point(size=2)




complementarity.hl <- ggplot(bombus,
                        aes(x = zfunctional.complementarity.HL,
                            y = PropGenusCrithidiaPresence,
                            color = Project)) +
  geom_point(size=2)


complementarity.ll <- ggplot(bombus,
                        aes(x = zfunctional.complementarity.LL,
                            y = PropGenusCrithidiaPresence,
                            color = Project)) +
  geom_point(size=2)

