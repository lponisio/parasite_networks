
rm(list=ls())
setwd('~/Dropbox (University of Oregon)/parasite_networks')

source("src/misc.R")
source("src/ggplotThemes.R")
source("src/genusNetworkPlotting.R")
source("src/ggplotThemes")

load(file="data/network_mets.RData")

library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)
library(ggnewscale)

# 1. Compute site-level prevalence and reshape
proj_long <- network.metrics %>%
  # drop rows with zero screened to avoid division by zero
  filter(GenusScreened > 0,
         Genus %in% c("Bombus", "Apis", "Melissodes"),
         ProjectSubProject != "PN-CA-FIRE") %>%
  mutate(
    prev_apicystis = GenusApicystisSpp    / GenusScreened,
    prev_crithidia = GenusCrithidiaPresence / GenusScreened
  ) %>%
  select(ProjectSubProject, Genus, GenusScreened,
         prev_apicystis, prev_crithidia) %>%
  pivot_longer(
    cols = starts_with("prev_"),
    names_to = "Parasite",
    values_to = "Prevalence"
  ) %>%
  mutate(
    Parasite = case_when(
      Parasite == "prev_apicystis" ~ "Apicystis",
      Parasite == "prev_crithidia" ~ "Crithidia",
      TRUE                         ~ Parasite
    )
  )

par.box <- ggplot(proj_long, aes(x = ProjectSubProject, y = Prevalence)) +
  # 1) Boxplots: fill by Genus
  geom_boxplot(
    aes(
      fill  = Genus,
      group = interaction(ProjectSubProject, Genus)
    ),
    outlier.shape = NA,
    alpha = 0.6,
    position = position_dodge(width = 0.8)
  ) +
  scale_fill_brewer(palette = "Dark2", name = "Genus") +

  # --- start a NEW fill scale for the points ---
  ggnewscale::new_scale_fill() +

  # 2) Points: fill by N screened (light → dark grey)
  geom_point(
    aes(
      fill  = log(GenusScreened),
      group = Genus
    ),
    position = position_jitterdodge(
      jitter.width = 0.15,
      dodge.width  = 0.8
    ),
    shape  = 21,        # filled circle: fill = N screened, black outline
    size   = 2.5,
    colour = "black",
    alpha  = 0.9
  ) +
  scale_fill_gradient(
    name = "N screened (log)",
    low  = "grey85",
    high = "grey10"
  ) +

  facet_wrap(~ Parasite, ncol = 1) +
  scale_y_continuous(
    labels = percent_format(accuracy = 1),
    limits = c(0, 1)
  ) +
  labs(
    y = "Infection prevalence (site-survey level)",
    x = ""
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.minor = element_blank()
  )

ggsave("figures/all_projects_boxplot.pdf",
       plot = par.box, width = 10, height = 15)


## Talk version of the plot -----------------------------------
par.box.talk <- ggplot(proj_long, aes(x = ProjectSubProject, y = Prevalence)) +
  # 1) Boxplots: fill by Genus
  geom_boxplot(
    aes(
      fill  = Genus,
      group = interaction(ProjectSubProject, Genus)
    ),
    outlier.shape = NA,
    alpha = 0.7,
    position = position_dodge(width = 0.8),
    linewidth = 0.7
  ) +
  scale_fill_brewer(palette = "Dark2", name = "Genus") +

  # --- new fill scale for points (N screened) ---
  ggnewscale::new_scale_fill() +

  # 2) Points: fill by log(N screened) (light → bright on black)
  geom_point(
    aes(
      fill  = log(GenusScreened),
      group = Genus
    ),
    position = position_jitterdodge(
      jitter.width = 0.2,
      dodge.width  = 0.8
    ),
    shape  = 21,        # filled circle
    size   = 4,         # bigger for talks
    colour = "black",
    alpha  = 0.95
  ) +
  scale_fill_gradient(
    name = "N screened (log)",
    low  = "grey30",
    high = "white"
  ) +

  facet_wrap(~ Parasite, ncol = 2) +
  scale_y_continuous(
    labels = percent_format(accuracy = 1),
    limits = c(0, 1)
  ) +
  labs(
    y = "Infection prevalence (site-survey level)",
    x = "",
  ) +
  theme_talk_black()

par.box.talk

## Save a talk-friendly version (PNG works well for slides)
ggsave(
  "figures/all_projects_boxplot_talk.pdf",
  plot   = par.box.talk,
  width  = 20,
  height = 10,
  dpi    = 300,
  bg     = "black"
)
