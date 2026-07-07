source("lab_paths.R")
setwd(local.path)
setwd("parasite_networks")

rm(list = ls())

source("src/writeResultsTable.R")
source("src/init.R")
source("src/makeBayesR2Table.R")
source("src/makeModelOutputTables.R")
source("src/runHostParasiteModels.R")

library(brms)
library(glmmTMB)
library(performance)
library(posterior)
library(bayesplot)

load(file = "data/network_mets.RData")

## Toggle these if you want to rerun only one part of the workflow.
run_models <- TRUE
run_post_checks <- TRUE
reset_output_tables <- TRUE

ncores <- 3
model_seed <- 123
pp_ndraws <- 10^3
set.seed(model_seed)

network.metrics <- network.metrics[!is.na(network.metrics$ProjectSubProject), ]
network.metrics$Year <- as.factor(network.metrics$Year)
network.metrics$Site <- as.factor(network.metrics$Site)
network.metrics$ProjectSubProject <- as.factor(network.metrics$ProjectSubProject)

network.metrics$UniNetwork <- paste(
  network.metrics$Site,
  network.metrics$Year,
  network.metrics$SampleRound
)

net.cols <- c(
  "zweighted.NODF",
  "zH2",
  "zweighted.cluster.coefficient.HL",
  "number.of.species.HL",
  "number.of.species.LL"
)

net.cols.scale <- paste0("scale(", net.cols, ")")

par.cols <- c("GenusApicystisSpp", "GenusCrithidiaPresence")

network_terms <- c(
  net.cols.scale,
  "ProjectSubProject",
  "(1|Site)",
  "(1|Year)"
)

make_genus_brms_formula <- function(response_col){
  bf(
    formula(
      paste0(
        response_col,
        " | trials(GenusScreened) ~ ",
        paste(network_terms, collapse = " + ")
      )
    ),
    family = "zero_inflated_binomial"
  )
}

make_genus_freq_formula <- function(response_col){
  make_binomial_check_formula(
    response_col = response_col,
    trials_col = "GenusScreened",
    terms = network_terms
  )
}

par.formulas <- setNames(
  lapply(par.cols, make_genus_brms_formula),
  par.cols
)

freq.par.formulas <- setNames(
  lapply(par.cols, make_genus_freq_formula),
  par.cols
)

## Basic sample summaries.
ms.sample.table <- network.metrics[
  network.metrics$Genus %in% c("Bombus", "Apis", "Melissodes"),
]

table(ms.sample.table$Genus, ms.sample.table$ProjectSubProject)

tapply(
  ms.sample.table$UniNetwork,
  ms.sample.table$ProjectSubProject,
  function(x) length(unique(x))
)

## *******************************************************
## Bombus
## *******************************************************

bombus <- network.metrics[network.metrics$Genus == "Bombus", ]
bombus.projects <- sort(unique(as.character(bombus$ProjectSubProject)))

tapply(bombus$GenusScreened, bombus$ProjectSubProject, sum)
tapply(bombus$GenusCrithidiaPresence, bombus$ProjectSubProject, sum)
tapply(bombus$GenusApicystisSpp, bombus$ProjectSubProject, sum)

## Not enough positives to model these parasites.
tapply(bombus$GenusNosemaBombi, bombus$ProjectSubProject, sum)
tapply(bombus$GenusNosemaCeranae, bombus$ProjectSubProject, sum)

sub.bombus <- list(
  GenusApicystisSpp = bombus[!bombus$ProjectSubProject %in% c("SF"), ],
  GenusCrithidiaPresence = bombus[!bombus$ProjectSubProject %in% c("SF"), ]
)

## *******************************************************
## Melissodes
## *******************************************************

melissodes.all <- network.metrics[network.metrics$Genus == "Melissodes", ]
melissodes.projects <- sort(unique(as.character(melissodes.all$ProjectSubProject)))

write.csv(
  melissodes.all[, c(net.cols, "Site", "Year", "ProjectSubProject",
                     "Genus", "GenusScreened", par.cols)],
  file = "data/melissodes_net.csv",
  row.names = FALSE
)

tapply(melissodes.all$GenusScreened, melissodes.all$ProjectSubProject, sum)
tapply(melissodes.all$GenusCrithidiaPresence, melissodes.all$ProjectSubProject, sum)
tapply(melissodes.all$GenusApicystisSpp, melissodes.all$ProjectSubProject, sum)

## Only enough screened individuals in these projects.
melissodes <- melissodes.all[
  melissodes.all$ProjectSubProject %in% c("SF", "SI"),
]

## *******************************************************
## Apis
## *******************************************************

apis.all <- network.metrics[network.metrics$Genus == "Apis", ]
apis.projects <- sort(unique(as.character(apis.all$ProjectSubProject)))

tapply(apis.all$GenusScreened, apis.all$ProjectSubProject, sum)
tapply(apis.all$GenusCrithidiaPresence, apis.all$ProjectSubProject, sum)
tapply(apis.all$GenusApicystisSpp, apis.all$ProjectSubProject, sum)

apis.sub <- apis.all[
  !apis.all$ProjectSubProject %in% c("PN-COAST", "SF"),
]

## *******************************************************
## Model specifications
## *******************************************************

model_specs <- list(
  list(
    model_id = "network_bombus_CrithidiaPresence",
    model_object_name = "bombus.CrithidiaPresence",
    data = sub.bombus$GenusCrithidiaPresence,
    brms_formula = par.formulas$GenusCrithidiaPresence,
    frequentist_formula = freq.par.formulas$GenusCrithidiaPresence,
    response_col = "GenusCrithidiaPresence",
    trials_col = "GenusScreened",
    all_projects = bombus.projects,
    control = list(adapt_delta = 0.999, stepsize = 0.001,
                   max_treedepth = 20)
  ),
  list(
    model_id = "network_bombus_ApicystisSpp",
    model_object_name = "bombus.ApicystisSpp",
    data = sub.bombus$GenusApicystisSpp,
    brms_formula = par.formulas$GenusApicystisSpp,
    frequentist_formula = freq.par.formulas$GenusApicystisSpp,
    response_col = "GenusApicystisSpp",
    trials_col = "GenusScreened",
    all_projects = bombus.projects,
    control = list(adapt_delta = 0.999, stepsize = 0.001,
                   max_treedepth = 20)
  ),
  list(
    model_id = "network_melissodes_CrithidiaPresence",
    model_object_name = "melissodes.CrithidiaPresence",
    data = melissodes,
    brms_formula = par.formulas$GenusCrithidiaPresence,
    frequentist_formula = freq.par.formulas$GenusCrithidiaPresence,
    response_col = "GenusCrithidiaPresence",
    trials_col = "GenusScreened",
    all_projects = melissodes.projects,
    control = list(adapt_delta = 0.9999, stepsize = 0.0001,
                   max_treedepth = 25)
  ),
  list(
    model_id = "network_melissodes_ApicystisSpp",
    model_object_name = "melissodes.ApicystisSpp",
    data = melissodes,
    brms_formula = par.formulas$GenusApicystisSpp,
    frequentist_formula = freq.par.formulas$GenusApicystisSpp,
    response_col = "GenusApicystisSpp",
    trials_col = "GenusScreened",
    all_projects = melissodes.projects,
    control = list(adapt_delta = 0.999, stepsize = 0.001,
                   max_treedepth = 20)
  ),
  list(
    model_id = "network_apis_CrithidiaPresence",
    model_object_name = "apis.CrithidiaPresence",
    data = apis.sub,
    brms_formula = par.formulas$GenusCrithidiaPresence,
    frequentist_formula = freq.par.formulas$GenusCrithidiaPresence,
    response_col = "GenusCrithidiaPresence",
    trials_col = "GenusScreened",
    all_projects = apis.projects,
    control = list(adapt_delta = 0.999, stepsize = 0.001,
                   max_treedepth = 20)
  ),
  list(
    model_id = "network_apis_ApicystisSpp",
    model_object_name = "apis.ApicystisSpp",
    data = apis.sub,
    brms_formula = par.formulas$GenusApicystisSpp,
    frequentist_formula = freq.par.formulas$GenusApicystisSpp,
    response_col = "GenusApicystisSpp",
    trials_col = "GenusScreened",
    all_projects = apis.projects,
    control = list(adapt_delta = 0.999, stepsize = 0.001,
                   max_treedepth = 20)
  )
)

## *******************************************************
## Fit models
## *******************************************************

if (isTRUE(run_models)) {
  for (spec in model_specs) {
    run_host_parasite_model(
      model_id = spec$model_id,
      model_object_name = spec$model_object_name,
      data = spec$data,
      brms_formula = spec$brms_formula,
      frequentist_formula = spec$frequentist_formula,
      response_col = spec$response_col,
      trials_col = spec$trials_col,
      all_projects = spec$all_projects,
      ncores = ncores,
      iter = 10^4,
      chains = 3,
      thin = 1,
      init = 0,
      seed = model_seed,
      control = spec$control,
      open_progress = FALSE
    )
  }
}

## *******************************************************
## Post-model checks and summary tables
## *******************************************************

if (isTRUE(run_post_checks)) {
  bayes_R2_table <- init_bayes_R2_table(reset = reset_output_tables)
  project_inclusion_table <- init_model_project_table(reset = reset_output_tables)

  for (spec in model_specs) {
    post <- check_saved_host_parasite_model(
      model_id = spec$model_id,
      model_file = file.path("saved", paste0(spec$model_id, ".Rdata")),
      resp = spec$response_col,
      ndraws = pp_ndraws,
      bayes_R2_table = bayes_R2_table,
      project_inclusion_table = project_inclusion_table,
      run_check_model = FALSE
    )

    bayes_R2_table <- post$bayes_R2_table
    project_inclusion_table <- post$project_inclusion_table
  }
}
