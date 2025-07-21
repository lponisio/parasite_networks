
rm(list=ls())
setwd("~/University of Oregon Dropbox/Lauren Ponisio/")
setwd("parasite_networks")
source('src/writeResultsTable.R')
source("src/init.R")
library(brms)
library(lme4)
library(glmmTMB)
library(car)
library(DHARMa)

load(file="data/sp_mets.RData")
ncores <- 3

net.cols <- c("zweighted.betweenness",
              "zweighted.closeness",
              "zd",
              ## "znestedrank",
              "zdegree")

net.cols.scale <- paste0("scale(", net.cols, ")")

par.cols <- c("SpApicystisSpp", "SpCrithidiaPresence",
              "SpNosemaBombi", "SpNosemaCeranae")


par.formulas <- vector(mode="list", length=length(par.cols))
freq.par.formulas <- vector(mode="list", length=length(par.cols))
names(freq.par.formulas) <- par.cols
names(par.formulas) <- par.cols

for(par in par.cols){
  par.formulas[[par]] <- bf(formula(paste0(par, "| trials(SpScreened)~", 
                         paste(c(net.cols.scale, "(1|Site)",
                                 "(1|Year)", "ProjectSubProject",
                                 "(1|gr(GenusSpecies, cov = phylo_matrix))"),
                               collapse="+"))), family="zero_inflated_binomial")

  freq.par.formulas[[par]] <-
    formula(paste0("cbind(", par, ", SpScreened)~", 
                   paste(c(net.cols.scale,  "(1|Site)", 
                                 "(1|Year)", "(1|Obs)", "ProjectSubProject"),
                         collapse="+")))
}


sp.network.metrics <- sp.network.metrics[!is.na(sp.network.metrics$ProjectSubProject),]
sp.network.metrics$Obs <- 1:nrow(sp.network.metrics)
sp.network.metrics$Year <- as.factor(sp.network.metrics$Year)
sp.network.metrics$Site <- as.factor(sp.network.metrics$Site)
sp.network.metrics$ProjectSubProject <- as.factor(sp.network.metrics$ProjectSubProject)

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

write.csv(bombus[, c(net.cols, "Site", "Year", "ProjectSubProject",
                     "GenusSpecies", "SpScreened",
                     par.cols)], file="data/bombus_net.csv")

## check models
for(form in freq.par.formulas){
  print(form)
  mod <- glmmTMB(form, data=bombus,
                             family="binomial",
                 ziformula=~1)
  print(summary(mod))
  ## print(vif(mod))
}

## *******************************************************

bombus.CrithidiaPresence <- brm(par.formulas$SpCrithidiaPresence,
                        bombus,
                        data2=list(phylo_matrix=phylo_matrix),
                        cores=ncores,
                        iter = 10^4,
                        chains =1,
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

## Some outliers, both otherwise looks okay
checked.bombus.CrithidiaPresence <-
  check_brms(bombus.CrithidiaPresence)

## Not sig, but not perfect
testDispersion(checked.bombus.CrithidiaPresence)

## Year SD skewed towards zero, otherwise okay
plot.res(bombus.CrithidiaPresence, "bombus_CrithidiaPresence")

## All rhats and ESS good
summary(bombus.CrithidiaPresence)

## .65, good
bayes_R2(bombus.CrithidiaPresence)

## Looks good
plot(pp_check(bombus.CrithidiaPresence,
              resp="SpCrithidiaPresence", ndraws=10^3))

## *******************************************************

## SF not converging, not very many bombus screened.

bombus.ApicystisSpp <- brm(par.formulas$SpApicystisSpp,
                        bombus[bombus$ProjectSubProject != "SF",],
                        cores=ncores,
                        iter = 10^4,
                        chains =1,
                        data2=list(phylo_matrix=phylo_matrix),
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

## looks good!
checked.bombus.ApicystisSpp <-
  check_brms(bombus.ApicystisSpp)

## looks good!
testDispersion(checked.bombus.ApicystisSpp)

## ProjectSubProjectSF and Year SD skewed, but otherwise okay
plot.res(bombus.ApicystisSpp, "bombus_ApicystisSpp")

## SF had issues, resolved with SF removed
summary(bombus.ApicystisSpp)

## .61, good
bayes_R2(bombus.ApicystisSpp)

## Looks good
plot(pp_check(bombus.ApicystisSpp, resp="SpApicystisSpp", ndraws=10^3))

## *******************************************************

## SI not converging, very few positives

bombus.SpNosemaBombi <- brm(par.formulas$SpNosemaBombi,
                        bombus[bombus$ProjectSubProject != "SI",],
                        cores=ncores,
                        data2=list(phylo_matrix=phylo_matrix),
                        iter = 10^4,
                        chains =1,
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

## good after dropping SI
checked.bombus.SpNosemaBombi <-
  check_brms(bombus.SpNosemaBombi)

## good after dropping SI
testDispersion(checked.bombus.SpNosemaBombi)

## Random effect SD skewed
plot.res(bombus.SpNosemaBombi, "bombus_SpNosemaBombi")

## Rhat and ESS good
summary(bombus.SpNosemaBombi)

## .41, good
bayes_R2(bombus.SpNosemaBombi)

## Good
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


write.csv(melissodes[, c(net.cols, "Site", "Year", "ProjectSubProject",
                     "GenusSpecies", "SpScreened",
                     par.cols)], file="data/melissodes_net.csv")


## check models, drop Nosemas
for(form in freq.par.formulas[1:2]){
  mod <- glmmTMB(form, data=melissodes,
                             family="binomial",
                 ziformula=~1)
  print(summary(mod))
}


## *******************************************************

melissodes.CrithidiaPresence <- brm(par.formulas$SpCrithidiaPresence,
                        melissodes,
                        cores=ncores,
                        data2=list(phylo_matrix=phylo_matrix),
                        iter = 10^4,
                        chains =1,
                        thin=1,
                        init=0,
                        open_progress = FALSE,
                        control = list(adapt_delta = 0.999,
                                       stepsize = 0.001,
                                       max_treedepth = 20))
                  
write.ms.table(melissodes.CrithidiaPresence, "melissodes_CrithidiaPresence")
save(melissodes, melissodes.CrithidiaPresence,
     file="saved/melissodes_CrithidiaPresence.Rdata")

load(file="saved/melissodes_CrithidiaPresence.Rdata")

## 
checked.melissodes.CrithidiaPresence <-
  check_brms(melissodes.CrithidiaPresence)

## 
testDispersion(melissodes.CrithidiaPresence)

plot.res(melissodes.CrithidiaPresence, "melissodes_CrithidiaPresence")

summary(melissodes.CrithidiaPresence)

bayes_R2(melissodes.CrithidiaPresence)

plot(pp_check(melissodes.CrithidiaPresence,
              resp="SpCrithidiaPresence", ndraws=10^3))
x
## *******************************************************

melissodes.ApicystisSpp <- brm(par.formulas$SpApicystisSpp,
                        melissodes,
                        cores=ncores,
                        data2=list(phylo_matrix=phylo_matrix),
                        iter = 10^4,
                        chains =1,
                        thin=1,
                        init=0,
                        open_progress = FALSE,
                        control = list(adapt_delta = 0.999,
                                       stepsize = 0.001,
                                       max_treedepth = 20))
                  
write.ms.table(melissodes.ApicystisSpp, "melissodes_ApicystisSpp")
save(melissodes, melissodes.ApicystisSpp,
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
                                 "(1|Year)", "ProjectSubProject"),
                               collapse="+"))), family="binomial")
}


## *******************************************************

apis.CrithidiaPresence <- brm(par.formulas$SpCrithidiaPresence,
                        apis[apis$ProjectSubProject != "PN-CA-FIRE",],
                        cores=ncores,
                        iter = 10^4,
                        chains =1,
                        thin=1,
                        init=0,
                        open_progress = FALSE,
                        control = list(adapt_delta = 0.999,
                                       stepsize = 0.001,
                                       max_treedepth = 20))
                  
write.ms.table(apis.CrithidiaPresence, "apis_CrithidiaPresence")
save(apis, apis.CrithidiaPresence,
     file="saved/apis_CrithidiaPresence.Rdata")

load(file="saved/apis_CrithidiaPresence.Rdata")

plot.res(apis.CrithidiaPresence, "apis_CrithidiaPresence")

summary(apis.CrithidiaPresence)

bayes_R2(apis.CrithidiaPresence)

plot(pp_check(apis.CrithidiaPresence,
              resp="SpCrithidiaPresence", ndraws=10^3))

## *******************************************************

apis.ApicystisSpp <- brm(par.formulas$SpApicystisSpp,
                        apis[apis$ProjectSubProject != "PN-CA-FIRE",],
                        cores=ncores,
                        iter = 10^4,
                        chains =1,
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
                        apis[apis$ProjectSubProject != "PN-CA-FIRE",],
                        cores=ncores,
                        iter = 10^4,
                        chains =1,
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
