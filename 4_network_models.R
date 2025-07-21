
rm(list=ls())
setwd("~/University of Oregon Dropbox/Lauren Ponisio/")
setwd("parasite_networks")
source('src/writeResultsTable.R')
source("src/init.R")
library(brms)
library(lme4)
library(glmmTMB)
library(performance)
library(DHARMa)

load(file="data/network_mets.RData")
ncores <- 3


network.metrics <- network.metrics[
  !is.na(network.metrics$ProjectSubProject),]
network.metrics$Year <- as.factor(network.metrics$Year)
network.metrics$Site <- as.factor(network.metrics$Site)
network.metrics$ProjectSubProject <- as.factor(
  network.metrics$ProjectSubProject)

net.cols <- c("connectance",
              "zweighted.NODF",
              "zH2",
              "zweighted.cluster.coefficient.HL")

net.cols.scale <- paste0("scale(", net.cols, ")")

par.cols <- c("GenusApicystisSpp", "GenusCrithidiaPresence",
              "GenusNosemaBombi", "GenusNosemaCeranae")


par.formulas <- vector(mode="list", length=length(par.cols))
freq.par.formulas <- vector(mode="list", length=length(par.cols))
names(freq.par.formulas) <- par.cols
names(par.formulas) <- par.cols

## *******************************************************
## Genus Specific models
## *******************************************************

par.cols <- c("GenusApicystisSpp", "GenusCrithidiaPresence",
              "GenusNosemaBombi", "GenusNosemaCeranae")


par.formulas <- vector(mode="list", length=length(par.cols))
freq.par.formulas <- vector(mode="list", length=length(par.cols))
names(freq.par.formulas) <- par.cols
names(par.formulas) <- par.cols

for(par in par.cols){
  par.formulas[[par]] <- bf(
    formula(paste0(par, "| trials(GenusScreened)~", 
                   paste(c(net.cols.scale,
                           "(1|Site)",
                           "(1|Year)",
                           "ProjectSubProject"),
                         collapse="+"))),
    family="zero_inflated_binomial")

  freq.par.formulas[[par]] <-
    formula(paste0("cbind(", par, ", GenusScreened)~", 
                   paste(c(net.cols.scale,
                           "(1|Site)", 
                           "(1|Year)",
                           "ProjectSubProject"),
                         collapse="+")))
}

## *******************************************************
## Bombus
## *******************************************************

bombus <- network.metrics[network.metrics$Genus == "Bombus",]

write.csv(bombus[, c(net.cols, "Site", "Year", "ProjectSubProject",
                     "GenusScreened",
                     par.cols)], file="data/bombus_network.csv")

## There aren't enough screenings in SI for Apicystis or Crithidia
## model to run. 

sub.bombus <- list(
  GenusApicystisSpp=bombus[bombus$ProjectSubProject != "SF",], 
  GenusCrithidiaPresence=bombus[bombus$ProjectSubProject != "SF",],
  GenusNosemaBombi= bombus[bombus$ProjectSubProject != "SI" &
                           bombus$ProjectSubProject != "PN-CA-FIRE",],
  GenusNosemaCeranae=bombus[bombus$ProjectSubProject != "SI" &
                            bombus$ProjectSubProject != "PN-CA-FIRE",]
)

## check models
for(i in names(freq.par.formulas)[1:3]){
  print(i)
  mod <- glmmTMB(freq.par.formulas[[i]],
                 data=sub.bombus[[i]],
                 family="binomial",
                 ziformula=~1)
  print(summary(mod))
  print(check_collinearity(mod))
}

## *******************************************************

bombus.CrithidiaPresence <- brm(par.formulas$GenusCrithidiaPresence,
                        sub.bombus$GenusCrithidiaPresence
                        cores=ncores,
                        iter = 10^4,
                        chains =1,
                        thin=1,
                        init=0,
                        open_progress = FALSE,
                        control = list(adapt_delta = 0.999,
                                       stepsize = 0.001,
                                       max_treedepth = 20))
                  
write.ms.table(bombus.CrithidiaPresence,
               "network_bombus_CrithidiaPresence")
save(bombus,bombus.CrithidiaPresence,
     file="saved/network_bombus_CrithidiaPresence.Rdata")

load(file="saved/network_bombus_CrithidiaPresence.Rdata")

## Some outliers, both otherwise looks okay
checked.bombus.CrithidiaPresence <-
  check_brms(bombus.CrithidiaPresence)

## Not sig, but not perfect
testDispersion(checked.bombus.CrithidiaPresence)

## Year SD skewed towards zero, otherwise okay
plot.res(bombus.CrithidiaPresence, "network_bombus_CrithidiaPresence")

## All rhats and ESS good
summary(bombus.CrithidiaPresence)

## .67, good
bayes_R2(bombus.CrithidiaPresence)

## Looks good
plot(pp_check(bombus.CrithidiaPresence,
              resp="GenusCrithidiaPresence", ndraws=10^3))

## *******************************************************

## SF not converging, not very many bombus screened.

bombus.ApicystisSpp <- brm(par.formulas$GenusApicystisSpp,
                        sub.bombus$GenusApicystisSpp,
                        cores=ncores,
                        iter = 10^4,
                        chains =1,
                        thin=1,
                        init=0,
                        open_progress = FALSE,
                        control = list(adapt_delta = 0.999,
                                       stepsize = 0.001,
                                       max_treedepth = 20))
                  
write.ms.table(bombus.ApicystisSpp, "network_bombus_ApicystisSpp")
save(bombus, bombus.ApicystisSpp,
     file="saved/network_bombus_ApicystisSpp.Rdata")

load(file="saved/network_bombus_ApicystisSpp.Rdata")

## looks good!
checked.bombus.ApicystisSpp <-
  check_brms(bombus.ApicystisSpp)

## looks good!
testDispersion(checked.bombus.ApicystisSpp)

## ProjectSubProjectSF and Year SD skewed, but otherwise okay
plot.res(bombus.ApicystisSpp, "network_bombus_ApicystisSpp")

## SF had issues, resolved with SF removed
summary(bombus.ApicystisSpp)

## .76, good
bayes_R2(bombus.ApicystisSpp)

## Looks good
plot(pp_check(bombus.ApicystisSpp, resp="SpApicystisSpp", ndraws=10^3))

## *******************************************************

## SI not converging, very few positives

bombus.NosemaBombi <- brm(par.formulas$GenusNosemaBombi,
                        sub.bombus$GenusNosemaBombi,
                        cores=ncores,
                        iter = 10^4,
                        chains =1,
                        thin=1,
                        init=0,
                        open_progress = FALSE,
                        control = list(adapt_delta = 0.999,
                                       stepsize = 0.001,
                                       max_treedepth = 20))
                  
write.ms.table(bombus.NosemaBombi,
               "network_bombus_NosemaBombi")
save(bombus,bombus.NosemaBombi,
     file="saved/network_bombus_NosemaBombi.Rdata")

load(file="saved/network_bombus_NosemaBombi.Rdata")

## good after dropping SI
checked.bombus.NosemaBombi <-
  check_brms(bombus.NosemaBombi)

## good after dropping SI
testDispersion(checked.bombus.NosemaBombi)

## Random effect SD skewed
plot.res(bombus.NosemaBombi, "network_bombus_NosemaBombi")

## Rhat and ESS good
summary(bombus.NosemaBombi)

## .41, good
bayes_R2(bombus.NosemaBombi)

## Good
plot(pp_check(bombus.NosemaBombi,
              resp="GenusNosemaBombi", ndraws=10^3))

## *******************************************************
## Melissodes
## *******************************************************

melissodes <- network.metrics[network.metrics$Genus ==
                                 "Melissodes",]


write.csv(melissodes[, c(net.cols, "Site", "Year", "ProjectSubProject",
                     "Genus", "GenusScreened",
                     par.cols)], file="data/melissodes_net.csv")



## check models
for(i in names(freq.par.formulas)[1:2]){
  print(i)
  mod <- glmmTMB(freq.par.formulas[[i]],
                 data=melissodes,
                 family="binomial",
                 ziformula=~1)
  print(summary(mod))
  print(check_collinearity(mod))
}


## *******************************************************

melissodes.CrithidiaPresence <- brm(par.formulas$GenusCrithidiaPresence,
                        melissodes,
                        cores=ncores,
                        iter = 10^4,
                        chains =1,
                        thin=1,
                        init=0,
                        open_progress = FALSE,
                        control = list(adapt_delta = 0.999,
                                       stepsize = 0.001,
                                       max_treedepth = 20))
                  
write.ms.table(melissodes.CrithidiaPresence,
               "network_melissodes_CrithidiaPresence")
save(melissodes, melissodes.CrithidiaPresence,
     file="saved/network_melissodes_CrithidiaPresence.Rdata")

load(file="saved/network_melissodes_CrithidiaPresence.Rdata")

## 
checked.melissodes.CrithidiaPresence <-
  check_brms(melissodes.CrithidiaPresence)

## 
testDispersion(melissodes.CrithidiaPresence)

plot.res(melissodes.CrithidiaPresence,
         "network_melissodes_CrithidiaPresence")

summary(melissodes.CrithidiaPresence)

bayes_R2(melissodes.CrithidiaPresence)

plot(pp_check(melissodes.CrithidiaPresence,
              resp="GenusCrithidiaPresence", ndraws=10^3))

## *******************************************************

melissodes.ApicystisSpp <- brm(par.formulas$GenusApicystisSpp,
                        melissodes,
                        cores=ncores,
                        iter = 10^4,
                        chains =1,
                        thin=1,
                        init=0,
                        open_progress = FALSE,
                        control = list(adapt_delta = 0.999,
                                       stepsize = 0.001,
                                       max_treedepth = 20))
                  
write.ms.table(melissodes.ApicystisSpp, "network_melissodes_ApicystisSpp")
save(melissodes, melissodes.ApicystisSpp,
     file="saved/network_melissodes_ApicystisSpp.Rdata")

load(file="saved/network_melissodes_ApicystisSpp.Rdata")

plot.res(melissodes.ApicystisSpp, "network_melissodes_ApicystisSpp")

summary(melissodes.ApicystisSpp)

bayes_R2(melissodes.ApicystisSpp)

plot(pp_check(melissodes.ApicystisSpp, resp="GenusApicystisSpp", ndraws=10^3))

## *******************************************************
## Apis
## *******************************************************

apis <- network.metrics[network.metrics$Genus == "Apis",]

apis.sub <- apis[apis$ProjectSubProject != "PN-CA-FIRE" &
                apis$ProjectSubProject != "PN-COAST",]
                 
## check models
for(i in names(freq.par.formulas)[c(1,2,4)]){
  print(i)
  mod <- glmmTMB(freq.par.formulas[[i]],
                 data=apis.sub,
                 family="binomial",
                 ziformula=~1)
  print(summary(mod))
  print(check_collinearity(mod))
}

## *******************************************************

apis.CrithidiaPresence <- brm(par.formulas$GenusCrithidiaPresence,
                        apis.sub,
                        cores=ncores,
                        iter = 10^4,
                        chains =1,
                        thin=1,
                        init=0,
                        open_progress = FALSE,
                        control = list(adapt_delta = 0.999,
                                       stepsize = 0.001,
                                       max_treedepth = 20))
                  
write.ms.table(apis.CrithidiaPresence, "network_apis_CrithidiaPresence")
save(apis, apis.CrithidiaPresence,
     file="saved/network_apis_CrithidiaPresence.Rdata")

load(file="saved/network_apis_CrithidiaPresence.Rdata")

plot.res(apis.CrithidiaPresence, "network_apis_CrithidiaPresence")

summary(apis.CrithidiaPresence)

## .51
bayes_R2(apis.CrithidiaPresence)

plot(pp_check(apis.CrithidiaPresence,
              resp="GenusCrithidiaPresence", ndraws=10^3))

## *******************************************************

apis.ApicystisSpp <- brm(par.formulas$GenusApicystisSpp,
                        apis.sub,
                        cores=ncores,
                        iter = 10^4,
                        chains =1,
                        thin=1,
                        init=0,
                        open_progress = FALSE,
                        control = list(adapt_delta = 0.999,
                                       stepsize = 0.001,
                                       max_treedepth = 20))
                  
write.ms.table(apis.ApicystisSpp, "network_apis_ApicystisSpp")
save(apis,apis.ApicystisSpp, 
     file="saved/network_apis_ApicystisSpp.Rdata")

load(file="saved/network_apis_ApicystisSpp.Rdata")

plot.res(apis.ApicystisSpp, "network_apis_ApicystisSpp")

summary(apis.ApicystisSpp)

## .52
bayes_R2(apis.ApicystisSpp)

plot(pp_check(apis.ApicystisSpp, resp="GenusApicystisSpp", ndraws=10^3))

## *******************************************************

apis.GenusNosemaCeranae <- brm(par.formulas$GenusNosemaCeranae,
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
                  
write.ms.table(apis.GenusNosemaCeranae, "network_apis_GenusNosemaCeranae")
save(apis,apis.GenusNosemaCeranae,
     file="saved/network_apis_GenusNosemaCeranae.Rdata")

load(file="saved/network_apis_GenusNosemaCeranae.Rdata")

plot.res(apis.GenusNosemaCeranae, "network_apis_GenusNosemaCeranae")

summary(apis.GenusNosemaCeranae)

## .58
bayes_R2(apis.GenusNosemaCeranae)

plot(pp_check(apis.GenusNosemaCeranae,
              resp="GenusNosemaCeranae", ndraws=10^3))
