
rm(list=ls())
setwd("~/University of Oregon Dropbox/Lauren Ponisio/")
setwd("parasite_networks")
library(brms)
library(performance)
library(lme4)
library(car)
ncores <- 3


source('src/writeResultsTable.R')
source("src/init.R")
source("src/runPlotFreqModelDiagnostics.R")

load(file="data/sp_mets.RData")
load(file="data/SI_popgen.RData")

## Subset to bombus since that is the only popgen data we have
bombus <- sp.network.metrics[sp.network.metrics$Genus == "Bombus",]
SI.popgen$Project <- "SI"

par.cols <- c("SpApicystisSpp", "SpCrithidiaPresence")
bombus <- bombus[, c("GenusSpecies", "Project", "Year", "Site",
                     par.cols, "SpScreened")]

## Merge popgen data
dim(bombus)
bombus <- merge(bombus, SI.popgen)
dim(bombus)  
## Drops all projects except SI and all species without popgen data


## *******************************************************
## Build and assess models
## *******************************************************

net.cols <- c("Ho",
             ## "Hs", # sooooo colinear
              "Fis")


net.cols.scale <- paste0("scale(", net.cols, ")")

par.formulas <- vector(mode="list", length=length(par.cols))
freq.par.formulas <- vector(mode="list", length=length(par.cols))
names(freq.par.formulas) <- par.cols
names(par.formulas) <- par.cols

for(par in par.cols){
  par.formulas[[par]] <- bf(formula(paste0(par, "| trials(SpScreened)~", 
                         paste(c(net.cols.scale,
                                "Year", "(1|Site)"),
                               collapse="+"))), family="binomial")

  freq.par.formulas[[par]] <-
    formula(paste0("cbind(", par, ", SpScreened)~", 
                   paste(c(net.cols.scale, "Year",
                           "(1|Site)"),
                         collapse="+")))
}


for(form in freq.par.formulas){
  mod <- glmer(form, data=bombus,
                             family="binomial",
                             glmerControl(optimizer ='optimx',
                                          optCtrl=list(method='nlminb')))
  print(summary(mod))
  print(vif(mod))
}

## VIF look good
## GenusSpecies is super colinear in the SpApicystisSpp model, but not
## the SpCrithidiaPresence. In either model, does come out as sig. 

## *******************************************************
## Run Baysian models using brms
## *******************************************************

bombus.CrithidiaPresence <- brm(par.formulas$SpCrithidiaPresence,
                        bombus,
                        cores=ncores,
                        iter = 10^4,
                        chains =3,
                        thin=1,
                        init=0,
                        open_progress = FALSE,
                        control = list(adapt_delta = 0.999,
                                       stepsize = 0.001,
                                       max_treedepth = 20))
                  
write.ms.table(bombus.CrithidiaPresence,
               "popgen_bombus_CrithidiaPresence")
save(bombus, bombus.CrithidiaPresence,
     file="saved/popgen_bombus_CrithidiaPresence.Rdata")

load(file="saved/popgen_bombus_CrithidiaPresence.Rdata")

plot.res(bombus.CrithidiaPresence, "popgen_bombus_CrithidiaPresence")

summary(bombus.CrithidiaPresence)

bayes_R2(bombus.CrithidiaPresence)

plot(pp_check(bombus.CrithidiaPresence,
              resp="SpCrithidiaPresence", ndraws=10^3))

## *******************************************************

bombus.ApicystisSpp <- brm(par.formulas$SpApicystisSpp,
                        bombus,
                        cores=ncores,
                        iter = 10^4,
                        chains =3,
                        thin=1,
                        init=0,
                        open_progress = FALSE,
                        control = list(adapt_delta = 0.999,
                                       stepsize = 0.001,
                                       max_treedepth = 20))
                  
write.ms.table(bombus.ApicystisSpp, "popgen_bombus_ApicystisSpp")
save(bombus, bombus.ApicystisSpp,
     file="saved/popgen_bombus_ApicystisSpp.Rdata")

load(file="saved/popgen_bombus_ApicystisSpp.Rdata")

plot.res(bombus.ApicystisSpp, "popgen_bombus_ApicystisSpp")

summary(bombus.ApicystisSpp)

bayes_R2(bombus.ApicystisSpp)

plot(pp_check(bombus.ApicystisSpp, resp="SpApicystisSpp", ndraws=10^3))

## Not enough positives  of other parasites to fit models
