
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

## Metric map Model formula uses scale(connectance),
## scale(zweighted.NODF), scale(zH2),
## scale(zweighted.cluster.coefficient.HL), so fixef names are
## "scaleconnectance", etc.

metrics <- tibble::tribble(
  ~effect_key,
  ~coef_term,
  ~xcol,
  ~xlab,
  "zweighted.NODF",
  "scalezweighted.NODF",
  "zweighted.NODF",
  "Nestedness (NODF)",
  "zweighted.cluster.coefficient.HL",
  "scalezweighted.cluster.coefficient.HL",
  "zweighted.cluster.coefficient.HL",
  "Clustering",
  "zH2",
  "scalezH2",
  "zH2",
  "Selectivity (H2)"
)

## ============================================================
## CRITHIDIA — Bombus / Apis / Melissodes
## ============================================================

##  Bombus 
load("saved/network_bombus_CrithidiaPresence.Rdata")
## object: bombus.CrithidiaPresence
bombus<-model_data

bombus_net <- bombus |> filter(!ProjectSubProject %in% c("SF", "PN-CA-FIRE"))

bombus_crith_fig <- make_network_figure(
  fit = bombus.CrithidiaPresence,
  raw_df = bombus_net,
  outcome_label = "Crithidia prevalence",
  prop_col = "PropGenusCrithidiaPresence",  ## if not present, it will compute successes/trials
  succ_col = "GenusCrithidiaPresence",
  trials_col = "GenusScreened",
  metrics_spec = metrics,
  ncol = 3,
  show_points = FALSE
)
ggsave("figures/network_bombus.CrithidiaPresence.pdf",
       plot = bombus_crith_fig, width = 20, height = 5)


## black talk theme + points off
bombus_crith_fig_dark <- make_network_figure(
  fit = bombus.CrithidiaPresence,
  raw_df = bombus_net,
  outcome_label = "Crithidia prevalence",
  prop_col = "PropGenusCrithidiaPresence",
  succ_col = "GenusCrithidiaPresence",
  trials_col = "GenusScreened",
  metrics_spec = metrics,
  ncol = 3,
  show_points = FALSE,
  theme = "talk_black"
)

ggsave("figures/network_bombus.CrithidiaPresence_dark.pdf",
       plot = bombus_crith_fig_dark, width = 20, height = 5)



##  Apis 
load("saved/network_apis_CrithidiaPresence.Rdata")
apis_net <- model_data[!model_data$ProjectSubProject %in% c("PN-CA-FIRE",
                                               "PN-COAST", "SF"),]
apis_crith_fig <- make_network_figure(
  fit = apis.CrithidiaPresence,
  raw_df = apis_net,
  outcome_label = "Crithidia prevalence",
  prop_col = "PropGenusCrithidiaPresence",
  succ_col = "GenusCrithidiaPresence",
  trials_col = "GenusScreened",
  metrics_spec = metrics,
  ncol=3,
  show_points = FALSE
)
ggsave("figures/network_apis.CrithidiaPresence.pdf",
       plot = apis_crith_fig,  width = 20, height = 5)

##  Melissodes 
load("saved/network_melissodes_CrithidiaPresence.Rdata")
## object: melissodes.CrithidiaPresence
mel_net <- model_data[model_data$ProjectSubProject %in% c("SF", "SI"),]
mel_crith_fig <- make_network_figure(
  fit = melissodes.CrithidiaPresence,
  raw_df = mel_net,
  outcome_label = "Crithidia prevalence",
  prop_col = "PropGenusCrithidiaPresence",
  succ_col = "GenusCrithidiaPresence",
  trials_col = "GenusScreened",
  metrics_spec = metrics,
  ncol=3
)
ggsave("figures/network_melissodes.CrithidiaPresence.pdf",
       plot = mel_crith_fig, width = 20, height = 5)

## ============================================================
## APICYSTIS — Bombus / Apis / Melissodes
## ============================================================

##  Bombus 
load("saved/network_bombus_ApicystisSpp.Rdata")


## object: bombus.ApicystisSpp
bombus_apic_fig <- make_network_figure(
  fit = bombus.ApicystisSpp,
  raw_df = bombus_net,  ## same subset rule as Crithidia
  outcome_label = "Apicystis prevalence",
  prop_col = "PropGenusApicystisSpp",
  succ_col = "GenusApicystisSpp",
  trials_col = "GenusScreened",
  metrics_spec = metrics,
  ncol=3,
  show_points = FALSE
)
ggsave("figures/network_bombus.ApicystisSpp.pdf",
       plot = bombus_apic_fig, width = 20, height = 5)



## black talk theme + points off
load("saved/network_bombus_ApicystisSpp.Rdata")
## object: bombus.ApicystisSpp
bombus_apic_fig <- make_network_figure(
  fit = bombus.ApicystisSpp,
  raw_df = bombus_net,  ## same subset rule as Crithidia
  outcome_label = "Apicystis prevalence",
  prop_col = "PropGenusApicystisSpp",
  succ_col = "GenusApicystisSpp",
  trials_col = "GenusScreened",
  metrics_spec = metrics,
  ncol=3,
  show_points = FALSE,
  theme = "talk_black"
)
ggsave("figures/network_bombus.ApicystisSpp_dark.pdf",
       plot = bombus_apic_fig, width = 20, height = 5)



##  Apis 
load("saved/network_apis_ApicystisSpp.Rdata")
## object: apis.ApicystisSpp
apis_apic_fig <- make_network_figure(
  fit = apis.ApicystisSpp,
  raw_df = apis_net,
  outcome_label = "Apicystis prevalence",
  prop_col = "PropGenusApicystisSpp",
  succ_col = "GenusApicystisSpp",
  trials_col = "GenusScreened",
  metrics_spec = metrics,
  ncol=3,
  show_points = FALSE
)
ggsave("figures/network_apis.ApicystisSpp.pdf",
       plot = apis_apic_fig, width = 20, height = 5)

##  Melissodes 
load("saved/network_melissodes_ApicystisSpp.Rdata")
## object: melissodes.ApicystisSpp
mel_apic_fig <- make_network_figure(
  fit = melissodes.ApicystisSpp,
  raw_df = mel_net,
  outcome_label = "Apicystis prevalence",
  prop_col = "PropGenusApicystisSpp",
  succ_col = "GenusApicystisSpp",
  trials_col = "GenusScreened",
  metrics_spec = metrics,
   ncol=3
)
ggsave("figures/network_melissodes.ApicystisSpp.pdf",
       plot = mel_apic_fig, width = 20, height = 5)


## ============================================================
## MANUSCRIPT STACKED FIGURES (per host)
## Rows = parasites (Crithidia, Apicystis); columns = network metrics
## Panel labels run A), B), C) across the first row, then D), E), F) on the second row.
## ============================================================

# Bombus
bombus_stack_fig <- make_stacked_parasite_network_figure(
  fit_top = bombus.CrithidiaPresence,
  fit_bottom = bombus.ApicystisSpp,
  raw_df_top = bombus_net,
  raw_df_bottom = bombus_net,
  top_outcome_label = "Crithidia prevalence",
  bottom_outcome_label = "Apicystis prevalence",
  top_prop_col = "PropGenusCrithidiaPresence",
  top_succ_col = "GenusCrithidiaPresence",
  top_trials_col = "GenusScreened",
  bottom_prop_col = "PropGenusApicystisSpp",
  bottom_succ_col = "GenusApicystisSpp",
  bottom_trials_col = "GenusScreened",
  metrics_spec = metrics,
  ncol = nrow(metrics),
  show_points = TRUE,
  theme = "ms"
)


ggsave("figures/network_bombus_parasites_stacked.pdf",
       plot = bombus_stack_fig, width = 20, height = 10)


# Apis
apis_stack_fig <- make_stacked_parasite_network_figure(
  fit_top = apis.CrithidiaPresence,
  fit_bottom = apis.ApicystisSpp,
  raw_df_top = apis_net,
  raw_df_bottom = apis_net,
  top_outcome_label = "Crithidia prevalence",
  bottom_outcome_label = "Apicystis prevalence",
  top_prop_col = "PropGenusCrithidiaPresence",
  top_succ_col = "GenusCrithidiaPresence",
  top_trials_col = "GenusScreened",
  bottom_prop_col = "PropGenusApicystisSpp",
  bottom_succ_col = "GenusApicystisSpp",
  bottom_trials_col = "GenusScreened",
  metrics_spec = metrics,
  ncol = nrow(metrics),
  show_points = TRUE,
  theme = "ms"
)

ggsave("figures/network_apis_parasites_stacked.pdf",
       plot = apis_stack_fig, width = 20, height = 10)


# Melissodes
mel_stack_fig <- make_stacked_parasite_network_figure(
  fit_top = melissodes.CrithidiaPresence,
  fit_bottom = melissodes.ApicystisSpp,
  raw_df_top = mel_net,
  raw_df_bottom = mel_net,
  top_outcome_label = "Crithidia prevalence",
  bottom_outcome_label = "Apicystis prevalence",
  top_prop_col = "PropGenusCrithidiaPresence",
  top_succ_col = "GenusCrithidiaPresence",
  top_trials_col = "GenusScreened",
  bottom_prop_col = "PropGenusApicystisSpp",
  bottom_succ_col = "GenusApicystisSpp",
  bottom_trials_col = "GenusScreened",
  metrics_spec = metrics,
  ncol = nrow(metrics),
  show_points = TRUE,
  theme = "ms"
)

ggsave("figures/network_melissodes_parasites_stacked.pdf",
       plot = mel_stack_fig, width = 20, height = 10)

