
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
library(DHARMa)
library(performance)
load(file="data/sp_mets.RData")
ncores <- 3
model_seed <- 123
set.seed(model_seed)

net.cols <- c("zweighted.betweenness",
              "zweighted.closeness",
              "zd",
              "zdegree",
              "zHBOverlap")

net.cols.scale <- paste0("scale(", net.cols, ")")

par.cols <- c("SpApicystisSpp", "SpCrithidiaPresence")

par.formulas <- vector(mode="list", length=length(par.cols))
freq.par.formulas <- vector(mode="list", length=length(par.cols))
names(freq.par.formulas) <- par.cols
names(par.formulas) <- par.cols

for(par in par.cols){
  par.formulas[[par]] <- bf(formula(
    paste0(par, "| trials(SpScreened)~", 
           paste(c(net.cols.scale, "(1|Site)",
                   "(1|Year)", "ProjectSubProject",
                   "(1|gr(GenusSpecies, cov = phylo_matrix))"),
                 collapse="+"))), family="zero_inflated_binomial")

  freq.par.formulas[[par]] <-
    formula(paste0("cbind(", par, ", SpScreened)~", 
                   paste(c(net.cols.scale,  "(1|Site)", 
                           "(1|Year)", "ProjectSubProject"),
                         collapse="+")))
}


sp.network.metrics <- sp.network.metrics[!is.na(sp.network.metrics$ProjectSubProject),]
sp.network.metrics$Obs <- 1:nrow(sp.network.metrics)
sp.network.metrics$Year <- as.factor(sp.network.metrics$Year)
sp.network.metrics$Site <- as.factor(sp.network.metrics$Site)
sp.network.metrics$ProjectSubProject <- as.factor(sp.network.metrics$ProjectSubProject)


ms.sample.table <- sp.network.metrics[sp.network.metrics$Genus %in%
                                  c("Bombus", "Apis", "Melissodes"),]

table(ms.sample.table$GenusSpecies,  ms.sample.table$ProjectSubProject)

## *******************************************************
## 
bombus <- sp.network.metrics[sp.network.metrics$Genus == "Bombus",]
bombus.projects <- sort(unique(as.character(bombus$ProjectSubProject)))

## vancouverensis was previously known as bifarius, change to match
## pylogeny before name change
bombus$GenusSpecies[bombus$GenusSpecies == "Bombus vancouverensis"] <-
  "Bombus bifarius"

load("data/Bombus_phylogeny.Rdata")
not.in.phylo <- unique(bombus$GenusSpecies[
  !bombus$GenusSpecies
  %in%
  phylo$tip.label])
not.in.phylo

write.csv(bombus[, c(net.cols, "Site", "Year", "ProjectSubProject",
                     "GenusSpecies", "SpScreened",
                     par.cols)], file="data/bombus_Spnet.csv")


## ## There aren't enough screenings/positives in SF for the
## ## moels to run
sub.bombus <- list(
  SpApicystisSpp=bombus[!bombus$ProjectSubProject %in% c("SF"),], 
  SpCrithidiaPresence=bombus[!bombus$ProjectSubProject %in% c("SF"),]
)

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


bombus.CrithidiaPresence <- brm(par.formulas$SpCrithidiaPresence,
                        sub.bombus$SpCrithidiaPresence,
                        data2=list(phylo_matrix=phylo_matrix),
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
                  
write.ms.table(bombus.CrithidiaPresence, "bombus_CrithidiaPresence")
save(bombus,bombus.CrithidiaPresence,
     file="saved/bombus_CrithidiaPresence.Rdata")

load(file="saved/bombus_CrithidiaPresence.Rdata")

## Year SD skewed towards zero, otherwise okay
plot.res(bombus.CrithidiaPresence, "bombus_CrithidiaPresence")

## All rhats and ESS good
summary(bombus.CrithidiaPresence)

bayes_R2_table <- record_bayes_R2(
  bayes_R2_table,
  bombus.CrithidiaPresence,
  model_id = "bombus_CrithidiaPresence"
)

## Looks good
project_inclusion_table <- record_model_projects(
  project_inclusion_table,
  model_id = "bombus_CrithidiaPresence",
  data = sub.bombus$SpCrithidiaPresence,
  response_col = "SpCrithidiaPresence",
  trials_col = "SpScreened",
  all_projects = bombus.projects
)

save_model_diagnostics(
  bombus.CrithidiaPresence,
  model_id = "bombus_CrithidiaPresence",
  resp = "SpCrithidiaPresence",
  ndraws = 10^3
)

## *******************************************************

bombus.ApicystisSpp <- brm(par.formulas$SpApicystisSpp,
                        sub.bombus$SpApicystisSpp,
                        cores=ncores,
                        iter = 10^4,
                        chains =3,
                        data2=list(phylo_matrix=phylo_matrix),
                        thin=1,
                        init=0,
                        open_progress = FALSE,
                        seed = model_seed,
                        control = list(adapt_delta = 0.999,
                                       stepsize = 0.001,
                                       max_treedepth = 20))
                  
write.ms.table(bombus.ApicystisSpp, "bombus_ApicystisSpp")
save(bombus, bombus.ApicystisSpp,
     file="saved/bombus_ApicystisSpp.Rdata")

load(file="saved/bombus_ApicystisSpp.Rdata")

plot.res(bombus.ApicystisSpp, "bombus_ApicystisSpp")

summary(bombus.ApicystisSpp)

bayes_R2_table <- record_bayes_R2(
  bayes_R2_table,
  bombus.ApicystisSpp,
  model_id = "bombus_ApicystisSpp"
)

## Looks good
project_inclusion_table <- record_model_projects(
  project_inclusion_table,
  model_id = "bombus_ApicystisSpp",
  data = sub.bombus$SpApicystisSpp,
  response_col = "SpApicystisSpp",
  trials_col = "SpScreened",
  all_projects = bombus.projects
)

save_model_diagnostics(
  bombus.ApicystisSpp,
  model_id = "bombus_ApicystisSpp",
  resp = "SpApicystisSpp",
  ndraws = 10^3
)

## *******************************************************
## Melissodes
## *******************************************************

melissodes <- sp.network.metrics[sp.network.metrics$Genus ==
                                 "Melissodes",]
melissodes.projects <- sort(unique(as.character(melissodes$ProjectSubProject)))
load("data/Melissodes_phylogeny.Rdata")
not.in.phylo <- unique(melissodes$GenusSpecies[
  !melissodes$GenusSpecies
  %in%
  phylo$tip.label])
not.in.phylo

melissodes <- melissodes[!melissodes$GenusSpecies %in% not.in.phylo,]


## Only enough screened individuals in these projects
melissodes <- melissodes[melissodes$ProjectSubProject %in% c("SF",
                                                             "SI"),]
melissodes$ProjectSubProject <- as.character(melissodes$ProjectSubProject)
table(melissodes$GenusSpecies,  melissodes$ProjectSubProject)
melissodes$Year <- as.character(melissodes$Year)
melissodes$Site <- as.character(melissodes$Site)
melissodes$GenusSpecies <- as.character(melissodes$GenusSpecies)

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
                        chains =3,
                        thin=1,
                        init=0,
                        open_progress = FALSE,
                        seed = model_seed,
                        control = list(adapt_delta = 0.999,
                                       stepsize = 0.001,
                                       max_treedepth = 20))
                  
write.ms.table(melissodes.CrithidiaPresence, "melissodes_CrithidiaPresence")
save(melissodes, melissodes.CrithidiaPresence,
     file="saved/melissodes_CrithidiaPresence.Rdata")

load(file="saved/melissodes_CrithidiaPresence.Rdata")

plot.res(melissodes.CrithidiaPresence, "melissodes_CrithidiaPresence")

summary(melissodes.CrithidiaPresence)

bayes_R2_table <- record_bayes_R2(
  bayes_R2_table,
  melissodes.CrithidiaPresence,
  model_id = "melissodes_CrithidiaPresence"
)

project_inclusion_table <- record_model_projects(
  project_inclusion_table,
  model_id = "melissodes_CrithidiaPresence",
  data = melissodes,
  response_col = "SpCrithidiaPresence",
  trials_col = "SpScreened",
  all_projects = melissodes.projects
)

save_model_diagnostics(
  melissodes.CrithidiaPresence,
  model_id = "melissodes_CrithidiaPresence",
  resp = "SpCrithidiaPresence",
  ndraws = 10^3
)

## *******************************************************

melissodes.ApicystisSpp <- brm(par.formulas$SpApicystisSpp,
                        melissodes,
                        cores=ncores,
                        data2=list(phylo_matrix=phylo_matrix),
                        iter = 10^4,
                        chains =3,
                        thin=1,
                        init=0,
                        open_progress = FALSE,
                        seed = model_seed,
                        control = list(adapt_delta = 0.999,
                                       stepsize = 0.001,
                                       max_treedepth = 20))
                  
write.ms.table(melissodes.ApicystisSpp, "melissodes_ApicystisSpp")
save(melissodes, melissodes.ApicystisSpp,
     file="saved/melissodes_ApicystisSpp.Rdata")

load(file="saved/melissodes_ApicystisSpp.Rdata")

plot.res(melissodes.ApicystisSpp, "melissodes_ApicystisSpp")

summary(melissodes.ApicystisSpp)

bayes_R2_table <- record_bayes_R2(
  bayes_R2_table,
  melissodes.ApicystisSpp,
  model_id = "melissodes_ApicystisSpp"
)

project_inclusion_table <- record_model_projects(
  project_inclusion_table,
  model_id = "melissodes_ApicystisSpp",
  data = melissodes,
  response_col = "SpApicystisSpp",
  trials_col = "SpScreened",
  all_projects = melissodes.projects
)

save_model_diagnostics(
  melissodes.ApicystisSpp,
  model_id = "melissodes_ApicystisSpp",
  resp = "SpApicystisSpp",
  ndraws = 10^3
)

## *******************************************************
## Apis
## *******************************************************

apis <- sp.network.metrics[sp.network.metrics$Genus == "Apis",]
apis.projects <- sort(unique(as.character(apis$ProjectSubProject)))
apis <- apis[!apis$ProjectSubProject %in% c( "PN-COAST", "SF"),]

## Drop phylogeny from models
par.formulas <- vector(mode="list", length=length(par.cols))
names(par.formulas) <- par.cols
for(par in par.cols){
  par.formulas[[par]] <- bf(formula(paste0(par, "| trials(SpScreened)~", 
                         paste(c(net.cols.scale[net.cols.scale!="scale(zHBOverlap)"],
                                 "(1|Site)",
                                 "(1|Year)", "ProjectSubProject"),
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
                        seed = model_seed,
                        control = list(adapt_delta = 0.999,
                                       stepsize = 0.001,
                                       max_treedepth = 20))
                  
write.ms.table(apis.CrithidiaPresence, "apis_CrithidiaPresence")
save(apis, apis.CrithidiaPresence,
     file="saved/apis_CrithidiaPresence.Rdata")

load(file="saved/apis_CrithidiaPresence.Rdata")

plot.res(apis.CrithidiaPresence, "apis_CrithidiaPresence")

summary(apis.CrithidiaPresence)

bayes_R2_table <- record_bayes_R2(
  bayes_R2_table,
  apis.CrithidiaPresence,
  model_id = "apis_CrithidiaPresence"
)

project_inclusion_table <- record_model_projects(
  project_inclusion_table,
  model_id = "apis_CrithidiaPresence",
  data = apis,
  response_col = "SpCrithidiaPresence",
  trials_col = "SpScreened",
  all_projects = apis.projects
)

save_model_diagnostics(
  apis.CrithidiaPresence,
  model_id = "apis_CrithidiaPresence",
  resp = "SpCrithidiaPresence",
  ndraws = 10^3
)

## *******************************************************

apis.ApicystisSpp <- brm(par.formulas$SpApicystisSpp,
                        apis,
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
                  
write.ms.table(apis.ApicystisSpp, "apis_ApicystisSpp")
save(apis,apis.ApicystisSpp, 
     file="saved/apis_ApicystisSpp.Rdata")

load(file="saved/apis_ApicystisSpp.Rdata")

plot.res(apis.ApicystisSpp, "apis_ApicystisSpp")

summary(apis.ApicystisSpp)

bayes_R2_table <- record_bayes_R2(
  bayes_R2_table,
  apis.ApicystisSpp,
  model_id = "apis_ApicystisSpp"
)

project_inclusion_table <- record_model_projects(
  project_inclusion_table,
  model_id = "apis_ApicystisSpp",
  data = apis,
  response_col = "SpApicystisSpp",
  trials_col = "SpScreened",
  all_projects = apis.projects
)

save_model_diagnostics(
  apis.ApicystisSpp,
  model_id = "apis_ApicystisSpp",
  resp = "SpApicystisSpp",
  ndraws = 10^3
)
