
rm(list=ls())
setwd('~/Dropbox (University of Oregon)/parasite_networks')

library(ggpubr)
library(tidyverse)
library(brms)
library(tidybayes)
library(viridis)

source("src/misc.R")
source("src/ggplotThemes.R")
source("src/genusNetworkPlotting.R")

## *************************************************************
## Bombus - Crithidia
## *************************************************************
## Load model & build CE object on demand inside helper
load("saved/bombus_CrithidiaPresence.Rdata")  
bombus<-model_data

## Restrict raw data as in your script
sub.bombus <- bombus %>%
  filter(!ProjectSubProject %in% c("SF", "PN-CA-FIRE"))


## Each row defines one panel.
## coef_term should match the corresponding row in fixef(fit).
metrics_spec <- tibble::tribble(
  ~effect_key,              ~coef_term,              ~xcol,                   ~xlab,
  "zweighted.betweenness",  "zweighted.betweenness",  "zweighted.betweenness",  "Betweenness centrality",
  "zweighted.closeness",    "zweighted.closeness",    "zweighted.closeness",    "Closeness centrality",
  "zdegree",               "zdegree",               "zdegree",               "Degree (number of plant partners)",
  "zd",                    "zd",                    "zd",                    "Partner selectivity d'",
  "zHBOverlap",            "zHBOverlap",            "zHBOverlap",            "Apis floral overlap"
)

bombus.CrithidiaPresence_fig <- make_network_figure(
  fit           = bombus.CrithidiaPresence,
  raw_df        = sub.bombus,
  outcome_label = "Crithidia infection probability",
  prop_col      = NULL,
  succ_col      = "SpCrithidiaPresence",
  trials_col    = "SpScreened",
  metrics_spec  = metrics_spec,
  ncol          = 3,
  show_points   = TRUE)

dir.create("figures", recursive = TRUE, showWarnings = FALSE)

ggsave(
  filename = "figures/bombus.CrithidiaPresence.pdf",
  plot     = bombus.CrithidiaPresence_fig,
  height   = 5,
  width    = 22,
  units    = "in",
  bg       = "black"
)



## *************************************************************
## Apis — Crithidia
## *************************************************************

load("saved/apis_CrithidiaPresence.Rdata")
apis <- model_data

sub.apis <- apis %>%
  filter(!ProjectSubProject %in% c("PN-CA-FIRE", "PN-COAST", "SF"))

apis_metrics_spec <- metrics_spec %>%
  filter(effect_key != "zHBOverlap")

apis.CrithidiaPresence_fig <- make_network_figure(
  fit           = apis.CrithidiaPresence,
  raw_df        = sub.apis,
  outcome_label = "Crithidia infection probability",
  prop_col      = NULL,
  succ_col      = "SpCrithidiaPresence",
  trials_col    = "SpScreened",
  metrics_spec  = apis_metrics_spec,
  ncol          = 4,
  show_points   = TRUE,
  theme         = "ms"
)

ggsave(
  filename = "figures/apis.CrithidiaPresence.pdf",
  plot     = apis.CrithidiaPresence_fig,
  height   = 5,
  width    = 20,
  units    = "in",
  bg       = "white"
)


## *************************************************************
## Bombus — Apicystis
## *************************************************************

load("saved/bombus_ApicystisSpp.Rdata")
bombus.apicystis_data <- model_data

sub.bombus.apicystis <- bombus.apicystis_data %>%
  filter(!ProjectSubProject %in% c("SF", "PN-CA-FIRE"))

bombus.ApicystisSpp_fig <- make_network_figure(
  fit           = bombus.ApicystisSpp,
  raw_df        = sub.bombus.apicystis,
  outcome_label = "Apicystis infection probability",
  prop_col      = NULL,
  succ_col      = "SpApicystisSpp",
  trials_col    = "SpScreened",
  metrics_spec  = metrics_spec,
  ncol          = 3,
  show_points   = TRUE,
  theme         = "ms"
)

ggsave(
  filename = "figures/bombus.ApicystisSpp.pdf",
  plot     = bombus.ApicystisSpp_fig,
  height   = 5,
  width    = 22,
  units    = "in",
  bg       = "white"
)

## *************************************************************
## Apis — Apicystis
## *************************************************************

load("saved/apis_ApicystisSpp.Rdata")
apis.apicystis_data <- model_data

sub.apis.apicystis <- apis.apicystis_data %>%
  filter(!ProjectSubProject %in% c("PN-CA-FIRE", "PN-COAST", "SF"))

apis.ApicystisSpp_fig <- make_network_figure(
  fit           = apis.ApicystisSpp,
  raw_df        = sub.apis.apicystis,
  outcome_label = "Apicystis infection probability",
  prop_col      = NULL,
  succ_col      = "SpApicystisSpp",
  trials_col    = "SpScreened",
  metrics_spec  = apis_metrics_spec,
  ncol          = 4,
  show_points   = TRUE,
  theme         = "ms"
)

ggsave(
  filename = "figures/apis.ApicystisSpp.pdf",
  plot     = apis.ApicystisSpp_fig,
  height   = 5,
  width    = 20,
  units    = "in",
  bg       = "white"
)





## *************************************************************
## Bombus — combined parasites: 2 rows × 5 columns
## *************************************************************

load("saved/bombus_CrithidiaPresence.Rdata")
bombus.crithidia_data <- model_data %>%
  filter(!ProjectSubProject %in% c("SF", "PN-CA-FIRE"))

load("saved/bombus_ApicystisSpp.Rdata")
bombus.apicystis_data <- model_data %>%
  filter(!ProjectSubProject %in% c("SF", "PN-CA-FIRE"))

bombus.parasites_fig <- make_stacked_parasite_network_figure(
  fit_top              = bombus.CrithidiaPresence,
  fit_bottom           = bombus.ApicystisSpp,
  raw_df_top           = bombus.crithidia_data,
  raw_df_bottom        = bombus.apicystis_data,
  top_outcome_label    = "Crithidia infection probability",
  bottom_outcome_label = "Apicystis infection probability",
  top_prop_col         = NULL,
  top_succ_col         = "SpCrithidiaPresence",
  top_trials_col       = "SpScreened",
  bottom_prop_col      = NULL,
  bottom_succ_col      = "SpApicystisSpp",
  bottom_trials_col    = "SpScreened",
  metrics_spec         = metrics_spec,
  ncol                 = 5,
  show_points          = TRUE,
  theme                = "ms"
)

dir.create("figures", recursive = TRUE, showWarnings = FALSE)

ggsave(
  filename = "figures/bombus_combined_parasites.pdf",
  plot     = bombus.parasites_fig,
  width    = 22,
  height   = 10,
  units    = "in",
  bg       = "white"
)

## *************************************************************
## Apis — combined parasites: 2 rows × 4 columns
## *************************************************************

apis_metrics_spec <- metrics_spec %>%
  filter(effect_key != "zHBOverlap")

load("saved/apis_CrithidiaPresence.Rdata")
apis.crithidia_data <- model_data %>%
  filter(!ProjectSubProject %in% c("PN-CA-FIRE", "PN-COAST", "SF"))

load("saved/apis_ApicystisSpp.Rdata")
apis.apicystis_data <- model_data %>%
  filter(!ProjectSubProject %in% c("PN-CA-FIRE", "PN-COAST", "SF"))

apis.parasites_fig <- make_stacked_parasite_network_figure(
  fit_top              = apis.CrithidiaPresence,
  fit_bottom           = apis.ApicystisSpp,
  raw_df_top           = apis.crithidia_data,
  raw_df_bottom        = apis.apicystis_data,
  top_outcome_label    = "Crithidia infection probability",
  bottom_outcome_label = "Apicystis infection probability",
  top_prop_col         = NULL,
  top_succ_col         = "SpCrithidiaPresence",
  top_trials_col       = "SpScreened",
  bottom_prop_col      = NULL,
  bottom_succ_col      = "SpApicystisSpp",
  bottom_trials_col    = "SpScreened",
  metrics_spec         = apis_metrics_spec,
  ncol                 = 4,
  show_points          = TRUE,
  theme                = "ms"
)

ggsave(
  filename = "figures/apis_combined_parasites.pdf",
  plot     = apis.parasites_fig,
  width    = 18,
  height   = 10,
  units    = "in",
  bg       = "white"
)






## *************************************************************
##  Bombus — ApicystisSpp
## *************************************************************

load("saved/bombus_ApicystisSpp.Rdata")
## provides bombus.ApicystisSpp
sub.bombus <- bombus |>
  filter(!ProjectSubProject %in% c("SF", "PN-CA-FIRE"))

                                        # panels (request the same set of network metrics; any missing are skipped)
pB_betweenness <- make_panel(bombus.ApicystisSpp,
                             "zweighted.betweenness",
                             sub.bombus, "zweighted.betweenness",
                             "SpApicystisSpp", "SpScreened",
                             "Betweenness centrality",
                             theme = "talk_black",
                             show_points = FALSE)

pB_closeness   <- make_panel(bombus.ApicystisSpp,
                             "zweighted.closeness",
                             sub.bombus, "zweighted.closeness",
                             "SpApicystisSpp", "SpScreened",
                             "Closeness centrality",
                             theme = "talk_black",
                             show_points = FALSE)

pB_degree      <- make_panel(bombus.ApicystisSpp, "zdegree",
                             sub.bombus, "zdegree",
                             "SpApicystisSpp", "SpScreened",
                             "Degree centrality",
                             theme = "talk_black",
                             show_points = FALSE)

pB_d           <- make_panel(bombus.ApicystisSpp, "zd",
                             sub.bombus, "zd",
                             "SpApicystisSpp", "SpScreened",
                             "d",
                             theme = "talk_black",
                             show_points = FALSE)

pB_overlap     <- make_panel(bombus.ApicystisSpp, "zHBOverlap",
                             sub.bombus, "zHBOverlap",
                             "SpApicystisSpp", "SpScreened",
                             "HB overlap",
                             theme = "talk_black",
                             show_points = FALSE)

bombus_panels <- list(pB_betweenness, pB_closeness,
                      pB_degree, pB_d, pB_overlap)
bombus_panels <- Filter(Negate(is.null), bombus_panels)

bombus.ApicystisSpp_fig <- ggarrange(
  plotlist = bombus_panels,
  labels = LETTERS[seq_along(bombus_panels)],
  ncol = 5, nrow = 1,
  common.legend = TRUE, legend = "bottom"
)

ggsave("figures/bombus.ApicystisSpp.pdf",
       plot = bombus.ApicystisSpp_fig, width = 22, height = 5)

## *************************************************************
## 2) Apis — ApicystisSpp
## *************************************************************

load("saved/apis_ApicystisSpp.Rdata")     # provides apis.ApicystisSpp

sub.apis <- apis[!apis$ProjectSubProject %in% c("PN-CA-FIRE",
                                               "PN-COAST", "SF"),]

pA_betweenness <- make_panel(apis.ApicystisSpp, "zweighted.betweenness",
                             sub.apis, "zweighted.betweenness",
                             "SpApicystisSpp", "SpScreened",
                             "Betweenness centrality")

pA_closeness   <- make_panel(apis.ApicystisSpp, "zweighted.closeness",
                             sub.apis, "zweighted.closeness",
                             "SpApicystisSpp", "SpScreened",
                             "Closeness centrality")

pA_degree      <- make_panel(apis.ApicystisSpp, "zdegree",
                             sub.apis, "zdegree",
                             "SpApicystisSpp", "SpScreened",
                             "Degree centrality")

pA_d           <- make_panel(apis.ApicystisSpp, "zd",
                             sub.apis, "zd",
                             "SpApicystisSpp", "SpScreened",
                             "d")

apis_panels <- list(pA_betweenness, pA_closeness, pA_degree, pA_d)
apis_panels <- Filter(Negate(is.null), apis_panels)

apis.ApicystisSpp_fig <- ggarrange(
  plotlist = apis_panels,
  labels = LETTERS[seq_along(apis_panels)],
  ncol = 4, nrow = 1,
  common.legend = TRUE, legend = "bottom"
)

ggsave("figures/apis.ApicystisSpp.pdf",
       plot = apis.ApicystisSpp_fig, width = 20, height = 5)


## *************************************************************
## Manuscript figures: stack parasites (CrithidiaPresence / ApicystisSpp) for each host
##   - columns: network metrics
##   - rows: parasites
##   - panel labels: A), B), C)... going left-to-right across rows
##   - one figure per host (Bombus / Apis / Melissodes)
## *************************************************************

metrics_species_all <- tibble::tibble(
  effect = c("zweighted.betweenness", "zweighted.closeness", "zdegree", "zd", "zHBOverlap"),
  x      = c("zweighted.betweenness", "zweighted.closeness", "zdegree", "zd", "zHBOverlap"),
  xlab   = c("Betweenness", "Closeness", "Degree", "d", "Apis overlap")
)

metrics_species_apis <- metrics_species_all |>
  dplyr::filter(effect != "zHBOverlap")  # Apis models typically omit overlap

parasite_order  <- c("CrithidiaPresence", "ApicystisSpp")
parasite_labels <- c("Crithidia", "Apicystis")
succ_cols       <- c(CrithidiaPresence = "SpCrithidiaPresence",
                     ApicystisSpp     = "SpApicystisSpp")
outcome_labels  <- c(CrithidiaPresence = "Crithidia prevalence",
                     ApicystisSpp     = "Apicystis prevalence")

# ---- 1) Bombus (stack parasites) ----
if (file.exists("saved/bombus_CrithidiaPresence.Rdata") &&
    file.exists("saved/bombus_ApicystisSpp.Rdata")) {

  load("saved/bombus_CrithidiaPresence.Rdata")
  load("saved/bombus_ApicystisSpp.Rdata")

  sub.bombus <- bombus %>%
    filter(!ProjectSubProject %in% c("SF", "PN-CA-FIRE"))

  bombus_species_parasites_stacked <- make_stacked_parasite_species_network_figure(
    fits_by_parasite = list(
      CrithidiaPresence = bombus.CrithidiaPresence,
      ApicystisSpp     = bombus.ApicystisSpp
    ),
    raw_df         = sub.bombus,
    outcome_labels = outcome_labels,
    succ_cols      = succ_cols,
    trials_cols    = "SpScreened",
    metrics_spec   = metrics_species_all,
    parasite_order = parasite_order,
    parasite_labels = parasite_labels,
    theme          = "ms",
    show_points    = FALSE,
    # Put metric names on the x-axis (not as top column titles)
    add_col_headers = FALSE,
    add_row_labels  = FALSE,
    label_style     = "paren",
    common_legend   = TRUE,
    legend          = "bottom"
  )

  ggsave("figures/species_network_bombus_parasites_stacked.pdf",
         plot = bombus_species_parasites_stacked,
         width = 22, height = 10)
}

# ---- 2) Apis (stack parasites) ----
if (file.exists("saved/apis_CrithidiaPresence.Rdata") &&
    file.exists("saved/apis_ApicystisSpp.Rdata")) {

  load("saved/apis_CrithidiaPresence.Rdata")
  load("saved/apis_ApicystisSpp.Rdata")

  sub.apis <- apis[!apis$ProjectSubProject %in% c("PN-CA-FIRE", "PN-COAST", "SF"),]

  apis_species_parasites_stacked <- make_stacked_parasite_species_network_figure(
    fits_by_parasite = list(
      CrithidiaPresence = apis.CrithidiaPresence,
      ApicystisSpp     = apis.ApicystisSpp
    ),
    raw_df         = sub.apis,
    outcome_labels = outcome_labels,
    succ_cols      = succ_cols,
    trials_cols    = "SpScreened",
    metrics_spec   = metrics_species_apis,
    parasite_order = parasite_order,
    parasite_labels = parasite_labels,
    theme          = "ms",
    show_points    = FALSE,
    # Put metric names on the x-axis (not as top column titles)
    add_col_headers = FALSE,
    add_row_labels  = FALSE,
    label_style     = "paren",
    common_legend   = TRUE,
    legend          = "bottom"
  )

  ggsave("figures/species_network_apis_parasites_stacked.pdf",
         plot = apis_species_parasites_stacked,
         width = 20, height = 10)
}

# ---- 3) Melissodes (stack parasites; optional) ----
if (file.exists("saved/melissodes_CrithidiaPresence.Rdata") &&
    file.exists("saved/melissodes_ApicystisSpp.Rdata")) {

  load("saved/melissodes_CrithidiaPresence.Rdata")
  load("saved/melissodes_ApicystisSpp.Rdata")

  sub.mel <- melissodes

  mel_species_parasites_stacked <- make_stacked_parasite_species_network_figure(
    fits_by_parasite = list(
      CrithidiaPresence = melissodes.CrithidiaPresence,
      ApicystisSpp     = melissodes.ApicystisSpp
    ),
    raw_df         = sub.mel,
    outcome_labels = outcome_labels,
    succ_cols      = succ_cols,
    trials_cols    = "SpScreened",
    metrics_spec   = metrics_species_all,
    parasite_order = parasite_order,
    parasite_labels = parasite_labels,
    theme          = "ms",
    show_points    = FALSE,
    # Put metric names on the x-axis (not as top column titles)
    add_col_headers = FALSE,
    add_row_labels  = FALSE,
    label_style     = "paren",
    common_legend   = TRUE,
    legend          = "bottom"
  )

  ggsave("figures/species_network_melissodes_parasites_stacked.pdf",
         plot = mel_species_parasites_stacked,
         width = 22, height = 10)
}

