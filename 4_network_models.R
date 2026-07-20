rm(list = ls())

source("lab_paths.R")
setwd(local.path)
setwd("parasite_networks")

source("src/writeResultsTable.R")
source("src/init.R")
source("src/makeBayesR2Table.R")
source("src/makeModelOutputTables.R")
source("src/makeModelResultsTable.R")
source("src/runHostParasiteModels.R")
source("src/modelWorkflowHelpers.R")

library(brms)
library(glmmTMB)
library(performance)
library(posterior)
library(bayesplot)

load(file = "data/network_mets.RData")

## Toggle these if you want to rerun only one part of the workflow.
run_models <- FALSE
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
    prior = priors_standard,
    control = control_standard
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
    prior = priors_standard,
    control = control_standard
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
    prior = priors_sparse,
    control = control_sparse
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
    prior = priors_sparse,
    control = control_sparse
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
    prior = priors_standard,
    control = control_minor_issue
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
    prior = priors_standard,
    control = control_minor_issue
  )
)

## *******************************************************
## Fit models
## *******************************************************

if (isTRUE(run_models)) {
  for (spec in model_specs) {
    run_model_from_spec(spec, ncores = ncores, model_seed = model_seed)
  }
}

## *******************************************************
## Load saved models, then run post-model checks and summary tables
## *******************************************************

if (isTRUE(run_post_checks)) {

  ## When run_models = FALSE, this loads all saved brmsfit objects into the
  ## workspace without refitting. This is optional for the automated post-check
  ## loop below, but useful if you want to run make.ms.table() manually or
  ## inspect fitted models interactively.
  saved_model_objects <- load_saved_host_parasite_models(
    model_specs = model_specs,
    saved_dir = "saved",
    assign_to_environment = TRUE,
    envir = .GlobalEnv,
    stop_if_missing = TRUE
  )

  post_tables <- post_check_saved_host_parasite_models(
    model_specs = model_specs,
    ndraws = pp_ndraws,
    reset_output_tables = reset_output_tables,
    saved_dir = "saved",
    run_check_model = FALSE,
    output_prefix = "network"
  )

  bayes_R2_table <- post_tables$bayes_R2_table
  project_inclusion_table <- post_tables$project_inclusion_table
  model_results_table <- post_tables$model_results_table

  sampler_diagnostic_table <- make_sampler_diagnostic_table(
    model_specs = model_specs,
    saved_model_objects = saved_model_objects,
    out_file = file.path("output", "model_diagnostics", "network_sampler_diagnostics.csv")
  )
}

