
rm(list=ls())
setwd("~/University of Oregon Dropbox/Lauren Ponisio/")
setwd("parasite_networks")
source('src/writeResultsTable.R')
source("src/init.R")
source("src/makeBayesR2Table.R")
source("src/makeModelOutputTables.R")

bayes_R2_table <- init_bayes_R2_table()
project_inclusion_table <- init_model_project_table()
library(brms)
library(lme4)
library(glmmTMB)
library(performance)
library(DHARMa)

load(file="data/network_mets.RData")
ncores <- 3
model_seed <- 123
set.seed(model_seed)


network.metrics <- network.metrics[
  !is.na(network.metrics$ProjectSubProject),]
network.metrics$Year <- as.factor(network.metrics$Year)
network.metrics$Site <- as.factor(network.metrics$Site)
network.metrics$ProjectSubProject <- as.factor(
  network.metrics$ProjectSubProject)

net.cols <- c("zweighted.NODF",
              "zH2",
              "zweighted.cluster.coefficient.HL",
              "number.of.species.HL",
              "number.of.species.LL")

net.cols.scale <- paste0("scale(", net.cols, ")")

par.cols <- c("GenusApicystisSpp", "GenusCrithidiaPresence")


par.formulas <- vector(mode="list", length=length(par.cols))
freq.par.formulas <- vector(mode="list", length=length(par.cols))
names(freq.par.formulas) <- par.cols
names(par.formulas) <- par.cols

## *******************************************************
## Genus Specific models
## *******************************************************

par.formulas <- vector(mode="list", length=length(par.cols))
freq.par.formulas <- vector(mode="list", length=length(par.cols))
names(freq.par.formulas) <- par.cols
names(par.formulas) <- par.cols

for(par in par.cols){
  par.formulas[[par]] <- bf(
    formula(paste0(par, "| trials(GenusScreened)~", 
                   paste(c(net.cols.scale, "ProjectSubProject",
                           "(1|Site)",
                           "(1|Year)"
                           ),
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

network.metrics$UniNetwork <- paste(network.metrics$Site,
                                    network.metrics$Year,
                                    network.metrics$SampleRound)

ms.sample.table <- network.metrics[network.metrics$Genus %in%
                                  c("Bombus", "Apis", "Melissodes"),]

table(ms.sample.table$Genus,
      ms.sample.table$ProjectSubProject)

tapply(ms.sample.table$UniNetwork,
       ms.sample.table$ProjectSubProject,
       function(x) length(unique(x)))

## *******************************************************
## Bombus
## *******************************************************
bombus <- network.metrics[network.metrics$Genus == "Bombus",]
bombus.projects <- sort(unique(as.character(bombus$ProjectSubProject)))

## reasonable sample sized
tapply(bombus$GenusScreened, bombus$ProjectSubProject, sum)
tapply(bombus$GenusCrithidiaPresence, bombus$ProjectSubProject, sum)
tapply(bombus$GenusApicystisSpp, bombus$ProjectSubProject, sum)

## Not enough positives to do anything with
tapply(bombus$GenusNosemaBombi, bombus$ProjectSubProject, sum)
tapply(bombus$GenusNosemaCeranae, bombus$ProjectSubProject, sum)

sub.bombus <- list(
  GenusApicystisSpp=bombus[!bombus$ProjectSubProject %in%
                           c("SF"),], 
  GenusCrithidiaPresence=bombus[!bombus$ProjectSubProject %in%
                                c("SF"),]
)

bombus.model.data <- do.call(
  rbind,
  lapply(names(sub.bombus), function(parasite_model){
    out <- sub.bombus[[parasite_model]]
    out$ParasiteModel <- parasite_model
    out
  })
)

write.csv(bombus.model.data[, c("ParasiteModel", net.cols, "Site", "Year",
                                "ProjectSubProject",
                                "GenusScreened",
                                par.cols)],
          file="data/bombus_network_model_data.csv",
          row.names = FALSE)


## check models
for(i in names(freq.par.formulas)){
  print(i)
  mod <- glmmTMB(freq.par.formulas[[i]],
                 data=sub.bombus[[i]],
                 family="binomial",
                 ziformula=~1)
  print(summary(mod))
  print(check_collinearity(mod))
}

## *******************************************************

print("Bombus Crithidia model")
bombus.CrithidiaPresence <- brm(par.formulas$GenusCrithidiaPresence,
                        sub.bombus$GenusCrithidiaPresence,
                        cores=ncores,
                        iter = 10^4,
                        chains =3,
                        thin=1,
                        init=0,
                        open_progress = FALSE,
                        seed = model_seed,
                        control = list(adapt_delta = 0.999,
                                       stepsize = 0.001,
                                       max_treedepth = 20))
                  
write.ms.table(bombus.CrithidiaPresence,
               "network_bombus_CrithidiaPresence")
save(bombus,bombus.CrithidiaPresence,
     file="saved/network_bombus_CrithidiaPresence.Rdata")

load(file="saved/network_bombus_CrithidiaPresence.Rdata")

## Some outliers, both otherwise looks okay

plot.res(bombus.CrithidiaPresence, "network_bombus_CrithidiaPresence")

summary(bombus.CrithidiaPresence)

bayes_R2_table <- record_bayes_R2(
  bayes_R2_table,
  bombus.CrithidiaPresence,
  model_id = "network_bombus_CrithidiaPresence"
)

project_inclusion_table <- record_model_projects(
  project_inclusion_table,
  model_id = "network_bombus_CrithidiaPresence",
  data = sub.bombus$GenusCrithidiaPresence,
  response_col = "GenusCrithidiaPresence",
  trials_col = "GenusScreened",
  all_projects = bombus.projects
)

save_model_diagnostics(
  bombus.CrithidiaPresence,
  model_id = "network_bombus_CrithidiaPresence",
  resp = "GenusCrithidiaPresence",
  ndraws = 10^3
)

## *******************************************************

print("Bombus ApicystisSpp model")
bombus.ApicystisSpp <- brm(par.formulas$GenusApicystisSpp,
                        sub.bombus$GenusApicystisSpp,
                        cores=ncores,
                        iter = 10^4,
                        chains =3,
                        thin=1,
                        init=0,
                        open_progress = FALSE,
                        seed = model_seed,
                        control = list(adapt_delta = 0.999,
                                       stepsize = 0.001,
                                       max_treedepth = 20))
                  
write.ms.table(bombus.ApicystisSpp, "network_bombus_ApicystisSpp")
save(bombus, bombus.ApicystisSpp,
     file="saved/network_bombus_ApicystisSpp.Rdata")

load(file="saved/network_bombus_ApicystisSpp.Rdata")

plot.res(bombus.ApicystisSpp, "network_bombus_ApicystisSpp")

summary(bombus.ApicystisSpp)

bayes_R2_table <- record_bayes_R2(
  bayes_R2_table,
  bombus.ApicystisSpp,
  model_id = "network_bombus_ApicystisSpp"
)

project_inclusion_table <- record_model_projects(
  project_inclusion_table,
  model_id = "network_bombus_ApicystisSpp",
  data = sub.bombus$GenusApicystisSpp,
  response_col = "GenusApicystisSpp",
  trials_col = "GenusScreened",
  all_projects = bombus.projects
)

save_model_diagnostics(
  bombus.ApicystisSpp,
  model_id = "network_bombus_ApicystisSpp",
  resp = "GenusApicystisSpp",
  ndraws = 10^3
)

## *******************************************************
## Melissodes
## *******************************************************

melissodes <- network.metrics[network.metrics$Genus ==
                                 "Melissodes",]
melissodes.projects <- sort(unique(as.character(melissodes$ProjectSubProject)))


write.csv(melissodes[, c(net.cols, "Site", "Year", "ProjectSubProject",
                     "Genus", "GenusScreened",
                     par.cols)], file="data/melissodes_net.csv")

tapply(melissodes$GenusScreened, melissodes$ProjectSubProject, sum)
tapply(melissodes$GenusCrithidiaPresence,
       melissodes$ProjectSubProject, sum)
tapply(melissodes$GenusApicystisSpp, melissodes$ProjectSubProject,
       sum)

## Only enough screened individuals in these projects
melissodes <- melissodes[melissodes$ProjectSubProject %in%
                         c("SF", "SI"),]

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
print("Melissodes Crithidia model")

melissodes.CrithidiaPresence <- brm(par.formulas$GenusCrithidiaPresence,
                        melissodes,
                        cores=ncores,
                        iter = 10^4,
                        chains =3,
                        thin=1,
                        init=0,
                        open_progress = FALSE,
                        seed = model_seed,
                        control = list(adapt_delta = 0.9999,
                                       stepsize = 0.0001,
                                       max_treedepth = 25))
                  
write.ms.table(melissodes.CrithidiaPresence,
               "network_melissodes_CrithidiaPresence")
save(melissodes, melissodes.CrithidiaPresence,
     file="saved/network_melissodes_CrithidiaPresence.Rdata")

load(file="saved/network_melissodes_CrithidiaPresence.Rdata")

plot.res(melissodes.CrithidiaPresence,
         "network_melissodes_CrithidiaPresence")

summary(melissodes.CrithidiaPresence)

bayes_R2_table <- record_bayes_R2(
  bayes_R2_table,
  melissodes.CrithidiaPresence,
  model_id = "network_melissodes_CrithidiaPresence"
)

project_inclusion_table <- record_model_projects(
  project_inclusion_table,
  model_id = "network_melissodes_CrithidiaPresence",
  data = melissodes,
  response_col = "GenusCrithidiaPresence",
  trials_col = "GenusScreened",
  all_projects = melissodes.projects
)

save_model_diagnostics(
  melissodes.CrithidiaPresence,
  model_id = "network_melissodes_CrithidiaPresence",
  resp = "GenusCrithidiaPresence",
  ndraws = 10^3
)

## *******************************************************

print("Melissodes ApicystisSpp model")
melissodes.ApicystisSpp <- brm(par.formulas$GenusApicystisSpp,
                        melissodes,
                        cores=ncores,
                        iter = 10^4,
                        chains =3,
                        thin=1,
                        init=0,
                        open_progress = FALSE,
                        seed = model_seed,
                        control = list(adapt_delta = 0.999,
                                       stepsize = 0.001,
                                       max_treedepth = 20))
                  
write.ms.table(melissodes.ApicystisSpp, "network_melissodes_ApicystisSpp")
save(melissodes, melissodes.ApicystisSpp,
     file="saved/network_melissodes_ApicystisSpp.Rdata")

load(file="saved/network_melissodes_ApicystisSpp.Rdata")

plot.res(melissodes.ApicystisSpp, "network_melissodes_ApicystisSpp")

summary(melissodes.ApicystisSpp)

bayes_R2_table <- record_bayes_R2(
  bayes_R2_table,
  melissodes.ApicystisSpp,
  model_id = "network_melissodes_ApicystisSpp"
)

project_inclusion_table <- record_model_projects(
  project_inclusion_table,
  model_id = "network_melissodes_ApicystisSpp",
  data = melissodes,
  response_col = "GenusApicystisSpp",
  trials_col = "GenusScreened",
  all_projects = melissodes.projects
)

save_model_diagnostics(
  melissodes.ApicystisSpp,
  model_id = "network_melissodes_ApicystisSpp",
  resp = "GenusApicystisSpp",
  ndraws = 10^3
)

## *******************************************************
## Apis
## *******************************************************

apis <- network.metrics[network.metrics$Genus == "Apis",]
apis.projects <- sort(unique(as.character(apis$ProjectSubProject)))

tapply(apis$GenusScreened, apis$ProjectSubProject, sum)
tapply(apis$GenusCrithidiaPresence, apis$ProjectSubProject, sum)
tapply(apis$GenusApicystisSpp, apis$ProjectSubProject, sum)

apis.sub <- apis[!apis$ProjectSubProject %in% c( "PN-COAST", "SF"),]
                 
## check models
for(i in names(freq.par.formulas)){
  print(i)
  mod <- glmmTMB(freq.par.formulas[[i]],
                 data=apis.sub,
                 family="binomial",
                 ziformula=~1)
  print(summary(mod))
  print(check_collinearity(mod))
}

## *******************************************************
print("Apis Crithidia model")
apis.CrithidiaPresence <- brm(par.formulas$GenusCrithidiaPresence,
                        apis.sub,
                        cores=ncores,
                        iter = 10^4,
                        chains =3,
                        thin=1,
                        init=0,
                        open_progress = FALSE,
                        seed = model_seed,
                        control = list(adapt_delta = 0.999,
                                       stepsize = 0.001,
                                       max_treedepth = 20))
                  
write.ms.table(apis.CrithidiaPresence, "network_apis_CrithidiaPresence")
save(apis, apis.CrithidiaPresence,
     file="saved/network_apis_CrithidiaPresence.Rdata")

load(file="saved/network_apis_CrithidiaPresence.Rdata")

plot.res(apis.CrithidiaPresence, "network_apis_CrithidiaPresence")

summary(apis.CrithidiaPresence)

bayes_R2_table <- record_bayes_R2(
  bayes_R2_table,
  apis.CrithidiaPresence,
  model_id = "network_apis_CrithidiaPresence"
)

project_inclusion_table <- record_model_projects(
  project_inclusion_table,
  model_id = "network_apis_CrithidiaPresence",
  data = apis.sub,
  response_col = "GenusCrithidiaPresence",
  trials_col = "GenusScreened",
  all_projects = apis.projects
)

save_model_diagnostics(
  apis.CrithidiaPresence,
  model_id = "network_apis_CrithidiaPresence",
  resp = "GenusCrithidiaPresence",
  ndraws = 10^3
)

## *******************************************************
print("Apis ApicystisSpp model")
apis.ApicystisSpp <- brm(par.formulas$GenusApicystisSpp,
                        apis.sub,
                        cores=ncores,
                        iter = 10^4,
                        chains =3,
                        thin=1,
                        init=0,
                        open_progress = FALSE,
                        seed = model_seed,
                        control = list(adapt_delta = 0.999,
                                       stepsize = 0.001,
                                       max_treedepth = 20))
                  
write.ms.table(apis.ApicystisSpp, "network_apis_ApicystisSpp")
save(apis,apis.ApicystisSpp, 
     file="saved/network_apis_ApicystisSpp.Rdata")

load(file="saved/network_apis_ApicystisSpp.Rdata")

plot.res(apis.ApicystisSpp, "network_apis_ApicystisSpp")

summary(apis.ApicystisSpp)

bayes_R2_table <- record_bayes_R2(
  bayes_R2_table,
  apis.ApicystisSpp,
  model_id = "network_apis_ApicystisSpp"
)

project_inclusion_table <- record_model_projects(
  project_inclusion_table,
  model_id = "network_apis_ApicystisSpp",
  data = apis.sub,
  response_col = "GenusApicystisSpp",
  trials_col = "GenusScreened",
  all_projects = apis.projects
)

save_model_diagnostics(
  apis.ApicystisSpp,
  model_id = "network_apis_ApicystisSpp",
  resp = "GenusApicystisSpp",
  ndraws = 10^3
)

## *******************************************************
