
rm(list=ls())
setwd('~/Dropbox (University of Oregon)/parasite_networks')

library(ggpubr)
library(tidyverse)
library(brms)
library(tidybayes)
library(viridis)

source("src/misc.R")
source("src/ggplotThemes.R")
source("src/spNetworkPlotting_fixed.R")

## *************************************************************
## Bombus - Crithidia
## *************************************************************
## Load model & build CE object on demand inside helper
load("saved/bombus_CrithidiaPresence.Rdata")  

## Restrict raw data as in your script
sub.bombus <- bombus %>%
  filter(!ProjectSubProject %in% c("SF", "PN-CA-FIRE"))

## Make the five plots in a compact, reproducible way
p_betweenness <- plot_ce_with_points(
  fit    = bombus.CrithidiaPresence,
  effect = "zweighted.betweenness",
  raw_df = sub.bombus,
  x      = "zweighted.betweenness",
  succ   = "SpCrithidiaPresence",
  trials = "SpScreened",
  xlab   = "Betweenness centrality",
  theme = "talk_black",
  show_points = FALSE
)

p_closeness <- plot_ce_with_points(
  bombus.CrithidiaPresence, "zweighted.closeness",
  sub.bombus, "zweighted.closeness", "SpCrithidiaPresence",
  "SpScreened",
  xlab = "Closeness centrality",
  theme = "talk_black",
  show_points = FALSE
)

p_degree <- plot_ce_with_points(
  bombus.CrithidiaPresence, "zdegree",
  sub.bombus, "zdegree", "SpCrithidiaPresence", "SpScreened",
  xlab = "Degree centrality",
  theme = "talk_black",
  show_points = FALSE
)

p_d <- plot_ce_with_points(
  bombus.CrithidiaPresence, "zd",
  sub.bombus, "zd", "SpCrithidiaPresence", "SpScreened",
  xlab = "d",
  theme = "talk_black",
  show_points = FALSE
)

p_overlap <- plot_ce_with_points(
  bombus.CrithidiaPresence, "zHBOverlap",
  sub.bombus, "zHBOverlap", "SpCrithidiaPresence", "SpScreened",
  xlab = "Apis overlap",
  theme = "talk_black",
  show_points = FALSE
)

## Arrange + save
bombus.CrithidiaPresence_fig <- ggarrange(
  p_betweenness, p_closeness, p_degree, p_d, p_overlap,
  labels = c("A","B","C","D","E"),
  ncol = 5, nrow = 1, common.legend = TRUE, legend = "bottom"
)

ggsave("figures/bombus.CrithidiaPresence.pdf",
       plot = bombus.CrithidiaPresence_fig, height = 5,
       width = 22)

## *************************************************************
## Apis - Crithidia
## *************************************************************

load("saved/apis_CrithidiaPresence.Rdata")  ## apis.CrithidiaPresence

sub.apis <- apis[!apis$ProjectSubProject %in% c("PN-CA-FIRE",
                                               "PN-COAST", "SF"),]


pA <- plot_ce_with_points(apis.CrithidiaPresence,
                          "zweighted.betweenness", sub.apis,
                          "zweighted.betweenness",
                          "SpCrithidiaPresence", "SpScreened",
                          xlab = "Betweenness centrality")

pB <- plot_ce_with_points(apis.CrithidiaPresence,
                          "zweighted.closeness", sub.apis,
                          "zweighted.closeness",
                          "SpCrithidiaPresence", "SpScreened",
                          xlab = "Closeness centrality")

pC <- plot_ce_with_points(apis.CrithidiaPresence,
                          "zdegree", sub.apis,
                          "zdegree", "SpCrithidiaPresence",
                          "SpScreened",
                          xlab = "Degree centrality")

pD <- plot_ce_with_points(apis.CrithidiaPresence,
                          "zd", sub.apis,
                          "zd", "SpCrithidiaPresence",
                          "SpScreened", xlab = "d")

apis.CrithidiaPresence_fig <- ggarrange(pA, pB, pC, pD,
                                        labels = c("A","B","C","D"),
                                        ncol = 4, nrow = 1,
                                        common.legend = TRUE,
                                        legend = "bottom"
                                        )
ggsave("figures/apis.CrithidiaPresence.pdf",
       plot = apis.CrithidiaPresence_fig, height = 5, width = 20)

## ## *************************************************************
## ## Melissodes - Crithidia
## ## *************************************************************
## load("saved/melissodes_CrithidiaPresence.Rdata")
## ## melissodes.CrithidiaPresence
## sub.mel <- melissodes

## mA <- plot_ce_with_points(melissodes.CrithidiaPresence,
##                           "zweighted.betweenness", sub.mel,
##                           "zweighted.betweenness",
##                           "SpCrithidiaPresence", "SpScreened",
##                           xlab = "Betweenness centrality")
## mB <- plot_ce_with_points(melissodes.CrithidiaPresence,
##                           "zweighted.closeness", sub.mel,
##                           "zweighted.closeness",
##                           "SpCrithidiaPresence", "SpScreened",
##                           xlab = "Closeness centrality")
## mC <- plot_ce_with_points(melissodes.CrithidiaPresence,
##                           "zdegree", sub.mel,
##                           "zdegree", "SpCrithidiaPresence",
##                           "SpScreened",
##                           xlab = "Degree centrality")
## mD <- plot_ce_with_points(melissodes.CrithidiaPresence,
##                           "zd", sub.mel,
##                           "zd", "SpCrithidiaPresence",
##                           "SpScreened", xlab = "d")

## melissodes.CrithidiaPresence_fig <- ggarrange(mA, mB, mC, mD,
##   labels = c("A","B","C","D"), ncol = 5, nrow = 1,
##   common.legend = TRUE, legend = "bottom"
## )
## ggsave("figures/melissodes.CrithidiaPresence.pdf",
##        plot = melissodes.CrithidiaPresence_fig, height = 5, width = 20)

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

## ## *************************************************************
## ## 3) Melissodes — SpApicystisSpp
## ## *************************************************************
## load("saved/melissodes_ApicystisSpp.Rdata")  # provides melissodes.ApicystisSpp
## sub.mel <- melissodes

## m_betweenness <- make_panel(melissodes.ApicystisSpp,
##                             "zweighted.betweenness",
##                             sub.mel, "zweighted.betweenness",
##                             "SpApicystisSpp", "SpScreened",
##                             "Betweenness centrality")

## m_closeness   <- make_panel(melissodes.ApicystisSpp,
##                             "zweighted.closeness",
##                             sub.mel, "zweighted.closeness",
##                             "SpApicystisSpp", "SpScreened",
##                             "Closeness centrality")

## m_degree      <- make_panel(melissodes.ApicystisSpp, "zdegree",
##                             sub.mel, "zdegree",
##                             "SpApicystisSpp", "SpScreened",
##                             "Degree centrality")

## m_d           <- make_panel(melissodes.ApicystisSpp, "zd",
##                             sub.mel, "zd",
##                             "SpApicystisSpp", "SpScreened",
##                             "d")

## m_overlap     <- make_panel(melissodes.ApicystisSpp, "zHBOverlap",
##                             sub.mel, "zHBOverlap",
##                             "SpApicystisSpp", "SpScreened",
##                             "HB overlap")

## mel_panels <- list(m_betweenness, m_closeness, m_degree, m_d, m_overlap)
## mel_panels <- Filter(Negate(is.null), mel_panels)

## melissodes.ApicystisSpp_fig <- ggarrange(
##   plotlist = mel_panels,
##   labels = LETTERS[seq_along(mel_panels)],
##   ncol = 5, nrow = 1,
##   common.legend = TRUE, legend = "bottom"
## )

## ggsave("figures/melissodes.ApicystisSpp.pdf",
##        plot = melissodes.ApicystisSpp_fig, width = 22, height = 5)

## *************************************************************
## Manuscript figures: stack hosts (Bombus / Apis / Melissodes) for each parasite
##   - columns: network metrics
##   - rows: hosts
##   - panel labels: A), B), C)... going left-to-right across rows
## *************************************************************

metrics_species <- tibble::tibble(
  effect = c("zweighted.betweenness", "zweighted.closeness", "zdegree", "zd", "zHBOverlap"),
  x      = c("zweighted.betweenness", "zweighted.closeness", "zdegree", "zd", "zHBOverlap"),
  xlab   = c("Betweenness", "Closeness", "Degree", "d", "Apis overlap")
)

# ---- 1) CrithidiaPresence (stack hosts) ----
fits_crith <- list()
raw_crith  <- list()
host_labs_crith <- c()

# Bombus
load("saved/bombus_CrithidiaPresence.Rdata")
sub.bombus <- bombus %>% filter(!ProjectSubProject %in% c("SF", "PN-CA-FIRE"))
fits_crith$bombus <- bombus.CrithidiaPresence
raw_crith$bombus  <- sub.bombus
host_labs_crith   <- c(host_labs_crith, "Bombus")

# Apis
load("saved/apis_CrithidiaPresence.Rdata")
sub.apis <- apis[!apis$ProjectSubProject %in% c("PN-CA-FIRE", "PN-COAST", "SF"),]
fits_crith$apis <- apis.CrithidiaPresence
raw_crith$apis  <- sub.apis
host_labs_crith <- c(host_labs_crith, "Apis")

# Melissodes (optional; included if the model file exists)
if (file.exists("saved/melissodes_CrithidiaPresence.Rdata")) {
  load("saved/melissodes_CrithidiaPresence.Rdata")
  sub.mel <- melissodes
  fits_crith$melissodes <- melissodes.CrithidiaPresence
  raw_crith$melissodes  <- sub.mel
  host_labs_crith <- c(host_labs_crith, "Melissodes")
}

species_Crithidia_hosts_stacked <- make_stacked_host_species_network_figure(
  fits_by_host   = fits_crith,
  raw_by_host    = raw_crith,
  outcome_label  = "Crithidia prevalence",
  succ_col       = "SpCrithidiaPresence",
  trials_col     = "SpScreened",
  metrics_spec   = metrics_species,
  host_order     = names(fits_crith),
  host_labels    = host_labs_crith,
  theme          = "ms",
  show_points    = FALSE,
  add_col_headers = TRUE,
  add_row_labels  = TRUE,
  label_style     = "paren",
  common_legend   = TRUE,
  legend          = "bottom"
)

ggsave(
  "figures/species_CrithidiaPresence_hosts_stacked.pdf",
  plot   = species_Crithidia_hosts_stacked,
  width  = 22,
  height = 3.6 * length(fits_crith) + 2
)

# ---- 2) ApicystisSpp (stack hosts) ----
fits_api <- list()
raw_api  <- list()
host_labs_api <- c()

# Bombus
load("saved/bombus_ApicystisSpp.Rdata")
sub.bombus <- bombus %>% filter(!ProjectSubProject %in% c("SF", "PN-CA-FIRE"))
fits_api$bombus <- bombus.ApicystisSpp
raw_api$bombus  <- sub.bombus
host_labs_api   <- c(host_labs_api, "Bombus")

# Apis
load("saved/apis_ApicystisSpp.Rdata")
sub.apis <- apis[!apis$ProjectSubProject %in% c("PN-CA-FIRE", "PN-COAST", "SF"),]
fits_api$apis <- apis.ApicystisSpp
raw_api$apis  <- sub.apis
host_labs_api <- c(host_labs_api, "Apis")

# Melissodes (optional)
if (file.exists("saved/melissodes_ApicystisSpp.Rdata")) {
  load("saved/melissodes_ApicystisSpp.Rdata")
  sub.mel <- melissodes
  fits_api$melissodes <- melissodes.ApicystisSpp
  raw_api$melissodes  <- sub.mel
  host_labs_api <- c(host_labs_api, "Melissodes")
}

species_Apicystis_hosts_stacked <- make_stacked_host_species_network_figure(
  fits_by_host   = fits_api,
  raw_by_host    = raw_api,
  outcome_label  = "Apicystis prevalence",
  succ_col       = "SpApicystisSpp",
  trials_col     = "SpScreened",
  metrics_spec   = metrics_species,
  host_order     = names(fits_api),
  host_labels    = host_labs_api,
  theme          = "ms",
  show_points    = FALSE,
  add_col_headers = TRUE,
  add_row_labels  = TRUE,
  label_style     = "paren",
  common_legend   = TRUE,
  legend          = "bottom"
)

ggsave(
  "figures/species_ApicystisSpp_hosts_stacked.pdf",
  plot   = species_Apicystis_hosts_stacked,
  width  = 22,
  height = 3.6 * length(fits_api) + 2
)
