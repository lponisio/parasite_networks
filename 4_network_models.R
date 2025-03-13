
rm(list=ls())
setwd("~/University of Oregon Dropbox/Lauren Ponisio/")
setwd("parasite_networks")
source('src/writeResultsTable.R')
source("src/init.R")
library(brms)

load(file="data/sp_mets.RData")
ncores <- 3

net.cols <- c("zweighted.betweenness",
              "zweighted.closeness",
              "zd", "znestedrank",
              "zdegree")

net.cols.scale <- paste0("scale(", net.cols, ")")

par.cols <- c("SpApicystisSpp", "SpCrithidiaPresence",
              "SpNosemaBombi", "SpNosemaCeranae")

par.formulas <- vector(mode="list", length=length(par.cols))
names(par.formulas) <- par.cols
for(par in par.cols){
  par.formulas[[par]] <- bf(formula(paste0(par, "| trials(SpScreened)~", 
                         paste(c(net.cols.scale, "(1|Site)",
                                 "(1|Year)", "Project",
                                 "(1|gr(GenusSpecies, cov = phylo_matrix))"),
                               collapse="+"))), family="binomial")
}

## *******************************************************
## Bombus
## *******************************************************

bombus <- sp.network.metrics[sp.network.metrics$Genus == "Bombus",]
load("data/Bombus_phylogeny.Rdata")
not.in.phylo <- unique(bombus$GenusSpecies[
  !bombus$GenusSpecies
  %in%
  phylo$tip.label])
not.in.phylo

bombus <- bombus[bombus$GenusSpecies != not.in.phylo,] 

## *******************************************************

bombus.CrithidiaPresence <- brm(par.formulas$SpCrithidiaPresence,
                        bombus,
                        data2=list(phylo_matrix=phylo_matrix),
                        cores=ncores,
                        iter = 10^4,
                        chains =3,
                        thin=1,
                        init=0,
                        open_progress = FALSE,
                        control = list(adapt_delta = 0.999,
                                       stepsize = 0.001,
                                       max_treedepth = 20))
                  
write.ms.table(bombus.CrithidiaPresence, "bombus_CrithidiaPresence")
save(bombus,bombus.CrithidiaPresence,
     file="saved/bombus_CrithidiaPresence.Rdata")

load(file="saved/bombus_CrithidiaPresence.Rdata")

plot.res(bombus.CrithidiaPresence, "bombus_CrithidiaPresence")

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
                  
write.ms.table(bombus.ApicystisSpp, "bombus_ApicystisSpp")
save(bombus, bombus.ApicystisSpp,
     file="saved/bombus_ApicystisSpp.Rdata")

load(file="saved/bombus_ApicystisSpp.Rdata")

plot.res(bombus.ApicystisSpp, "bombus_ApicystisSpp")

summary(bombus.ApicystisSpp)

bayes_R2(bombus.ApicystisSpp)

plot(pp_check(bombus.ApicystisSpp, resp="SpApicystisSpp", ndraws=10^3))

## *******************************************************

bombus.SpNosemaBombi <- brm(par.formulas$SpNosemaBombi,
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
                  
write.ms.table(bombus.SpNosemaBombi, "bombus_SpNosemaBombi")
save(bombus,bombus.SpNosemaBombi,
     file="saved/bombus_SpNosemaBombi.Rdata")

load(file="saved/bombus_SpNosemaBombi.Rdata")

plot.res(bombus.SpNosemaBombi, "bombus_SpNosemaBombi")

summary(bombus.SpNosemaBombi)

bayes_R2(bombus.SpNosemaBombi)

plot(pp_check(bombus.SpNosemaBombi, resp="SpSpNosemaBombi", ndraws=10^3))

## *******************************************************
## Melissodes
## *******************************************************

melissodes <- sp.network.metrics[sp.network.metrics$Genus ==
                                 "Melissodes",]
load("data/Melissodes_phylogeny.Rdata")
not.in.phylo <- unique(melissodes$GenusSpecies[
  !melissodes$GenusSpecies
  %in%
  phylo$tip.label])
not.in.phylo

melissodes <- melissodes[melissodes$GenusSpecies != not.in.phylo,]

## *******************************************************

melissodes.CrithidiaPresence <- brm(par.formulas$SpCrithidiaPresence,
                        melissodes,
                        cores=ncores,
                        iter = 10^4,
                        chains =3,
                        thin=1,
                        init=0,
                        open_progress = FALSE,
                        control = list(adapt_delta = 0.999,
                                       stepsize = 0.001,
                                       max_treedepth = 20))
                  
write.ms.table(melissodes.CrithidiaPresence, "melissodes_CrithidiaPresence")
save(melissodes,melissodes.CrithidiaPresence,
     file="saved/melissodes_CrithidiaPresence.Rdata")

load(file="saved/melissodes_CrithidiaPresence.Rdata")

plot.res(melissodes.CrithidiaPresence, "melissodes_CrithidiaPresence")

summary(melissodes.CrithidiaPresence)

bayes_R2(melissodes.CrithidiaPresence)

plot(pp_check(melissodes.CrithidiaPresence,
              resp="SpCrithidiaPresence", ndraws=10^3))

## *******************************************************

melissodes.ApicystisSpp <- brm(par.formulas$SpApicystisSpp,
                        melissodes,
                        cores=ncores,
                        iter = 10^4,
                        chains =3,
                        thin=1,
                        init=0,
                        open_progress = FALSE,
                        control = list(adapt_delta = 0.999,
                                       stepsize = 0.001,
                                       max_treedepth = 20))
                  
write.ms.table(melissodes.ApicystisSpp, "melissodes_ApicystisSpp")
save(melissodes,melissodes.ApicystisSpp,
     file="saved/melissodes_ApicystisSpp.Rdata")

load(file="saved/melissodes_ApicystisSpp.Rdata")

plot.res(melissodes.ApicystisSpp, "melissodes_ApicystisSpp")

summary(melissodes.ApicystisSpp)

bayes_R2(melissodes.ApicystisSpp)

plot(pp_check(melissodes.ApicystisSpp, resp="SpApicystisSpp", ndraws=10^3))

## *******************************************************
## Apis
## *******************************************************

apis <- sp.network.metrics[sp.network.metrics$Genus == "Apis",]

## Drop phylogeny from models
par.formulas <- vector(mode="list", length=length(par.cols))
names(par.formulas) <- par.cols
for(par in par.cols){
  par.formulas[[par]] <- bf(formula(paste0(par, "| trials(SpScreened)~", 
                         paste(c(net.cols.scale, "(1|Site)",
                                 "(1|Year)", "Project",
                               collapse="+"))), family="binomial")
}


## *******************************************************

apis.CrithidiaPresence <- brm(par.formulas$SpCrithidiaPresence,
                        apis,
                        cores=ncores,
                        iter = 10^4,
                        chains =3,
                        thin=1,
                        init=0,
                        open_progress = FALSE,
                        control = list(adapt_delta = 0.999,
                                       stepsize = 0.001,
                                       max_treedepth = 20))
                  
write.ms.table(apis.CrithidiaPresence, "apis_CrithidiaPresence")
save(apis,apis.CrithidiaPresence,
     file="saved/apis_CrithidiaPresence.Rdata")

load(file="saved/apis_CrithidiaPresence.Rdata")

plot.res(apis.CrithidiaPresence, "apis_CrithidiaPresence")

summary(apis.CrithidiaPresence)

bayes_R2(apis.CrithidiaPresence)

plot(pp_check(apis.CrithidiaPresence,
              resp="SpCrithidiaPresence", ndraws=10^3))

## *******************************************************

apis.ApicystisSpp <- brm(par.formulas$SpApicystisSpp,
                        apis,
                        cores=ncores,
                        iter = 10^4,
                        chains =3,
                        thin=1,
                        init=0,
                        open_progress = FALSE,
                        control = list(adapt_delta = 0.999,
                                       stepsize = 0.001,
                                       max_treedepth = 20))
                  
write.ms.table(apis.ApicystisSpp, "apis_ApicystisSpp")
save(apis,apis.ApicystisSpp,
     file="saved/apis_ApicystisSpp.Rdata")

load(file="saved/apis_ApicystisSpp.Rdata")

plot.res(apis.ApicystisSpp, "apis_ApicystisSpp")

summary(apis.ApicystisSpp)

bayes_R2(apis.ApicystisSpp)

plot(pp_check(apis.ApicystisSpp, resp="SpApicystisSpp", ndraws=10^3))

## *******************************************************

apis.SpNosemaCeranae <- brm(par.formulas$SpNosemaCeranae,
                        apis,
                        cores=ncores,
                        iter = 10^4,
                        chains =3,
                        thin=1,
                        init=0,
                        open_progress = FALSE,
                        control = list(adapt_delta = 0.999,
                                       stepsize = 0.001,
                                       max_treedepth = 20))
                  
write.ms.table(apis.SpNosemaCeranae, "apis_SpNosemaCeranae")
save(apis,apis.SpNosemaCeranae,
     file="saved/apis_SpNosemaCeranae.Rdata")

load(file="saved/apis_SpNosemaCeranae.Rdata")

plot.res(apis.SpNosemaCeranae, "apis_SpNosemaCeranae")

summary(apis.SpNosemaCeranae)

bayes_R2(apis.SpNosemaCeranae)

plot(pp_check(apis.SpNosemaCeranae, resp="SpSpNosemaCeranae", ndraws=10^3))
