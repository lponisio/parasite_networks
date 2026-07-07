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

load(file = "data/sp_mets.RData")

## Toggle these if you want to rerun only one part of the workflow.
run_models <- TRUE
run_post_checks <- TRUE
reset_output_tables <- FALSE

ncores <- 3
model_seed <- 123
pp_ndraws <- 10^3
set.seed(model_seed)

net.cols <- c(
  "zweighted.betweenness",
  "zweighted.closeness",
  "zd",
  "zdegree",
  "zHBOverlap"
)

net.cols.scale <- paste0("scale(", net.cols, ")")

par.cols <- c("SpApicystisSpp", "SpCrithidiaPresence")

sp.network.metrics <- sp.network.metrics[
  !is.na(sp.network.metrics$ProjectSubProject),
]

sp.network.metrics$Obs <- seq_len(nrow(sp.network.metrics))
sp.network.metrics$Year <- as.factor(sp.network.metrics$Year)
sp.network.metrics$Site <- as.factor(sp.network.metrics$Site)
sp.network.metrics$ProjectSubProject <- as.factor(
  sp.network.metrics$ProjectSubProject
)

ms.sample.table <- sp.network.metrics[
  sp.network.metrics$Genus %in% c("Bombus", "Apis", "Melissodes"),
]

table(ms.sample.table$GenusSpecies, ms.sample.table$ProjectSubProject)

phylo_terms <- c(
  net.cols.scale,
  "(1|Site)",
  "(1|Year)",
  "ProjectSubProject",
  "(1|gr(GenusSpecies, cov = phylo_matrix))"
)

freq_terms <- c(
  net.cols.scale,
  "(1|Site)",
  "(1|Year)",
  "ProjectSubProject"
)

make_species_phylo_brms_formula <- function(response_col){
  bf(
    formula(
      paste0(
        response_col,
        " | trials(SpScreened) ~ ",
        paste(phylo_terms, collapse = " + ")
      )
    ),
    family = "zero_inflated_binomial"
  )
}

make_species_freq_formula <- function(response_col, terms = freq_terms){
  make_binomial_check_formula(
    response_col = response_col,
    trials_col = "SpScreened",
    terms = terms
  )
}

phylo.par.formulas <- setNames(
  lapply(par.cols, make_species_phylo_brms_formula),
  par.cols
)

freq.par.formulas <- setNames(
  lapply(par.cols, make_species_freq_formula),
  par.cols
)

## *******************************************************
## Bombus
## *******************************************************

bombus <- sp.network.metrics[sp.network.metrics$Genus == "Bombus", ]
bombus.projects <- sort(unique(as.character(bombus$ProjectSubProject)))

## vancouverensis was previously known as bifarius; change to match
## phylogeny before the name change.
bombus$GenusSpecies[bombus$GenusSpecies == "Bombus vancouverensis"] <-
  "Bombus bifarius"

load("data/Bombus_phylogeny.Rdata")
bombus_phylo_matrix <- phylo_matrix
not.in.phylo <- unique(bombus$GenusSpecies[
  !bombus$GenusSpecies %in% phylo$tip.label
])
not.in.phylo

write.csv(
  bombus[, c(net.cols, "Site", "Year", "ProjectSubProject",
             "GenusSpecies", "SpScreened", par.cols)],
  file = "data/bombus_Spnet.csv",
  row.names = FALSE
)

## There are not enough screenings/positives in SF for the models to run.
sub.bombus <- list(
  SpApicystisSpp = bombus[!bombus$ProjectSubProject %in% c("SF"), ],
  SpCrithidiaPresence = bombus[!bombus$ProjectSubProject %in% c("SF"), ]
)

## *******************************************************
## Melissodes
## *******************************************************

melissodes.all <- sp.network.metrics[
  sp.network.metrics$Genus == "Melissodes",
]

melissodes.projects <- sort(unique(as.character(
  melissodes.all$ProjectSubProject
)))

load("data/Melissodes_phylogeny.Rdata")
melissodes_phylo_matrix <- phylo_matrix

not.in.phylo <- unique(melissodes.all$GenusSpecies[
  !melissodes.all$GenusSpecies %in% phylo$tip.label
])
not.in.phylo

melissodes <- melissodes.all[
  !melissodes.all$GenusSpecies %in% not.in.phylo,
]

## Only enough screened individuals in these projects.
melissodes <- melissodes[
  melissodes$ProjectSubProject %in% c("SF", "SI"),
]

melissodes$ProjectSubProject <- as.character(melissodes$ProjectSubProject)
melissodes$Year <- as.character(melissodes$Year)
melissodes$Site <- as.character(melissodes$Site)
melissodes$GenusSpecies <- as.character(melissodes$GenusSpecies)

table(melissodes$GenusSpecies, melissodes$ProjectSubProject)

write.csv(
  melissodes[, c(net.cols, "Site", "Year", "ProjectSubProject",
                 "GenusSpecies", "SpScreened", par.cols)],
  file = "data/melissodes_net.csv",
  row.names = FALSE
)

## *******************************************************
## Apis
## *******************************************************

apis.all <- sp.network.metrics[sp.network.metrics$Genus == "Apis", ]
apis.projects <- sort(unique(as.character(apis.all$ProjectSubProject)))

apis <- apis.all[
  !apis.all$ProjectSubProject %in% c("PN-COAST", "SF"),
]

## Drop phylogeny from Apis models.
apis_terms <- c(
  net.cols.scale[net.cols.scale != "scale(zHBOverlap)"],
  "(1|Site)",
  "(1|Year)",
  "ProjectSubProject"
)

make_apis_brms_formula <- function(response_col){
  bf(
    formula(
      paste0(
        response_col,
        " | trials(SpScreened) ~ ",
        paste(apis_terms, collapse = " + ")
      )
    ),
    family = "binomial"
  )
}

apis.par.formulas <- setNames(
  lapply(par.cols, make_apis_brms_formula),
  par.cols
)

apis.freq.formulas <- setNames(
  lapply(par.cols, function(par){
    make_species_freq_formula(par, terms = apis_terms)
  }),
  par.cols
)

## *******************************************************
## Model specifications
## *******************************************************

model_specs <- list(
  list(
    model_id = "bombus_CrithidiaPresence",
    model_object_name = "bombus.CrithidiaPresence",
    data = sub.bombus$SpCrithidiaPresence,
    brms_formula = phylo.par.formulas$SpCrithidiaPresence,
    frequentist_formula = freq.par.formulas$SpCrithidiaPresence,
    data2 = list(phylo_matrix = bombus_phylo_matrix),
    response_col = "SpCrithidiaPresence",
    trials_col = "SpScreened",
    all_projects = bombus.projects,
    control = list(adapt_delta = 0.999, stepsize = 0.001,
                   max_treedepth = 20)
  ),
  list(
    model_id = "bombus_ApicystisSpp",
    model_object_name = "bombus.ApicystisSpp",
    data = sub.bombus$SpApicystisSpp,
    brms_formula = phylo.par.formulas$SpApicystisSpp,
    frequentist_formula = freq.par.formulas$SpApicystisSpp,
    data2 = list(phylo_matrix = bombus_phylo_matrix),
    response_col = "SpApicystisSpp",
    trials_col = "SpScreened",
    all_projects = bombus.projects,
    control = list(adapt_delta = 0.999, stepsize = 0.001,
                   max_treedepth = 20)
  ),
  list(
    model_id = "melissodes_CrithidiaPresence",
    model_object_name = "melissodes.CrithidiaPresence",
    data = melissodes,
    brms_formula = phylo.par.formulas$SpCrithidiaPresence,
    frequentist_formula = freq.par.formulas$SpCrithidiaPresence,
    data2 = list(phylo_matrix = melissodes_phylo_matrix),
    response_col = "SpCrithidiaPresence",
    trials_col = "SpScreened",
    all_projects = melissodes.projects,
    control = list(adapt_delta = 0.999, stepsize = 0.001,
                   max_treedepth = 20)
  ),
  list(
    model_id = "melissodes_ApicystisSpp",
    model_object_name = "melissodes.ApicystisSpp",
    data = melissodes,
    brms_formula = phylo.par.formulas$SpApicystisSpp,
    frequentist_formula = freq.par.formulas$SpApicystisSpp,
    data2 = list(phylo_matrix = melissodes_phylo_matrix),
    response_col = "SpApicystisSpp",
    trials_col = "SpScreened",
    all_projects = melissodes.projects,
    control = list(adapt_delta = 0.999, stepsize = 0.001,
                   max_treedepth = 20)
  ),
  list(
    model_id = "apis_CrithidiaPresence",
    model_object_name = "apis.CrithidiaPresence",
    data = apis,
    brms_formula = apis.par.formulas$SpCrithidiaPresence,
    frequentist_formula = apis.freq.formulas$SpCrithidiaPresence,
    response_col = "SpCrithidiaPresence",
    trials_col = "SpScreened",
    all_projects = apis.projects,
    control = list(adapt_delta = 0.999, stepsize = 0.001,
                   max_treedepth = 20)
  ),
  list(
    model_id = "apis_ApicystisSpp",
    model_object_name = "apis.ApicystisSpp",
    data = apis,
    brms_formula = apis.par.formulas$SpApicystisSpp,
    frequentist_formula = apis.freq.formulas$SpApicystisSpp,
    response_col = "SpApicystisSpp",
    trials_col = "SpScreened",
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
      data2 = spec$data2,
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
