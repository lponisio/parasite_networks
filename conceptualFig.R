## ------------------------------------------------------------
## Conceptual figure: plant–pollinator networks & parasite spread
## Clean legends, black background, wider panels
## ------------------------------------------------------------

library(tidyverse)
library(patchwork)

set.seed(123)

## ------------------------------------------------------------
## Global black theme (applied at the end)
## ------------------------------------------------------------

theme_black <- theme_void(base_size = 11) +
  theme(
    plot.background   = element_rect(fill = "black", colour = NA),
    panel.background  = element_rect(fill = "black", colour = NA),
    text              = element_text(colour = "white"),
    plot.title        = element_text(colour = "white"),
    plot.subtitle     = element_text(colour = "white"),
    legend.text       = element_text(colour = "white"),
    legend.title      = element_text(colour = "white"),
    legend.background = element_rect(fill = "black", colour = NA),
    legend.key        = element_rect(fill = "black", colour = NA)
  )

## ------------------------------------------------------------
## Helper functions
## ------------------------------------------------------------

make_nodes <- function(plants, pollinators) {
  tibble(
    name  = c(plants, pollinators),
    guild = c(rep("Plant", length(plants)),
              rep("Pollinator", length(pollinators)))
  ) |>
    group_by(guild) |>
    mutate(
      x = seq(0, 1, length.out = n()),
      y = if_else(guild == "Plant", 0, 1)
    ) |>
    ungroup()
}

prep_edges <- function(edge_df, nodes_df,
                       plant_col = "plant",
                       poll_col  = "pollinator") {
  edge_df |>
    left_join(nodes_df |> select(name, x, y),
              by = setNames("name", plant_col)) |>
    rename(x_pl = x, y_pl = y) |>
    left_join(nodes_df |> select(name, x, y),
              by = setNames("name", poll_col)) |>
    rename(x_bee = x, y_bee = y)
}

## ============================================================
## PANEL A: Network structure → community parasite prevalence
## Strongly modular vs highly connected / nested
## ============================================================

plants_A <- paste0("P", 1:6)
polls_A  <- paste0("B", 1:7)

# Meadow: three tight modules, very few cross-links
edges_A_meadow <- tribble(
  ~pollinator, ~plant,
  # Module 1
  "B1", "P1",
  "B1", "P2",
  "B2", "P1",
  "B2", "P2",
  # Module 2
  "B3", "P3",
  "B3", "P4",
  "B4", "P3",
  "B5", "P4",
  # Module 3
  "B6", "P5",
  "B6", "P6",
  "B7", "P6"
)

# Managed: higher connectance + nested hub (B1) linking modules
edges_A_managed <- tribble(
  ~pollinator, ~plant,
  # Generalist hub B1 connects all plants
  "B1", "P1",
  "B1", "P2",
  "B1", "P3",
  "B1", "P4",
  "B1", "P5",
  "B1", "P6",
  # Others more generalized & cross-module
  "B2", "P1",
  "B2", "P2",
  "B2", "P3",
  "B3", "P2",
  "B3", "P3",
  "B3", "P4",
  "B4", "P3",
  "B4", "P4",
  "B4", "P5",
  "B5", "P3",
  "B5", "P4",
  "B5", "P5",
  "B6", "P4",
  "B6", "P5",
  "B6", "P6",
  "B7", "P5",
  "B7", "P6"
)

# Custom layout for meadow to show three modules clearly
nodes_A_meadow <- tribble(
  ~name, ~guild,       ~x,   ~y,
  "P1",  "Plant",      0.10, 0,
  "P2",  "Plant",      0.22, 0,
  "P3",  "Plant",      0.50, 0,
  "P4",  "Plant",      0.62, 0,
  "P5",  "Plant",      0.82, 0,
  "P6",  "Plant",      0.94, 0,
  "B1",  "Pollinator", 0.16, 1,
  "B2",  "Pollinator", 0.26, 1,
  "B3",  "Pollinator", 0.50, 1,
  "B4",  "Pollinator", 0.60, 1,
  "B5",  "Pollinator", 0.44, 1,
  "B6",  "Pollinator", 0.84, 1,
  "B7",  "Pollinator", 0.94, 1
)

# Same layout reused for managed network
nodes_A_managed <- nodes_A_meadow

# Conceptual infection patterns
infected_meadow  <- c("B2", "B3", "B4", "P2", "P3")      # localized
infected_managed <- c("B1", "B2", "B3", "B4", "B5",
                      "P1", "P2", "P3", "P4")            # widespread

nodes_A_meadow <- nodes_A_meadow |>
  mutate(
    infection = case_when(
      name %in% infected_meadow ~ "Moderate",
      guild == "Pollinator"     ~ "Low",
      TRUE                      ~ "None"
    )
  )

nodes_A_managed <- nodes_A_managed |>
  mutate(
    infection = case_when(
      name %in% infected_managed ~ "High",
      guild == "Pollinator"      ~ "Moderate",
      TRUE                       ~ "Low"
    )
  )

edges_A_meadow_plot  <- prep_edges(edges_A_meadow,  nodes_A_meadow)
edges_A_managed_plot <- prep_edges(edges_A_managed, nodes_A_managed)

# Background module bands
modules_bg <- tribble(
  ~xmin, ~xmax,
  0.03,  0.32,   # module 1 band
  0.40,  0.70,   # module 2 band
  0.78,  1.02    # module 3 band
)

pA_meadow <- ggplot() +
  geom_rect(
    data = modules_bg,
    aes(xmin = xmin, xmax = xmax,
        ymin = -0.25, ymax = 1.25),
    inherit.aes = FALSE,
    fill = c("#202020", "#303030", "#202020"),
    colour = NA
  ) +
  geom_segment(
    data = edges_A_meadow_plot,
    aes(x = x_pl, y = y_pl, xend = x_bee, yend = y_bee),
    colour = "grey60",
    alpha = 0.7
  ) +
  geom_point(
    data = nodes_A_meadow,
    aes(x = x, y = y, shape = guild, colour = infection),
    size = 4.2,
    fill = "black",
    stroke = 0.6
  ) +
  scale_shape_manual(values = c(Plant = 21, Pollinator = 24), guide = "none") +
  scale_color_manual(
    values = c(
      "None"     = "grey40",
      "Low"      = "grey70",
      "Moderate" = "orange",
      "High"     = "red3"
    ),
    guide = "none"
  ) +
  coord_fixed(
    ratio = 0.35,
    xlim = c(0, 1.05),
    ylim = c(-0.25, 1.25)
  ) +
  theme_void(base_size = 11) +
  theme(
    plot.subtitle = element_text(hjust = 0),
    legend.position = "none",
    plot.margin   = margin(5.5, 25, 5.5, 5.5)
  ) +
  labs(
    subtitle = "Natural meadow:\nstrong modularity and higher specialization (higher H2′)"
  )

pA_managed <- ggplot() +
  geom_segment(
    data = edges_A_managed_plot,
    aes(x = x_pl, y = y_pl, xend = x_bee, yend = y_bee),
    colour = "grey60",
    alpha = 0.7
  ) +
  geom_point(
    data = nodes_A_managed,
    aes(x = x, y = y, shape = guild, colour = infection),
    size = 4.2,
    fill = "black",
    stroke = 0.6
  ) +
  # emphasize B1 as generalist hub
  geom_point(
    data = dplyr::filter(nodes_A_managed, name == "B1"),
    aes(x = x, y = y),
    shape = 24,
    size  = 6,
    fill  = "black",
    colour = "red3",
    stroke = 0.9,
    inherit.aes = FALSE
  ) +
  geom_text(
    data = dplyr::filter(nodes_A_managed, name == "B1"),
    aes(x = x, y = y + 0.18),
    label = "Generalist hub",
    size  = 3
  ) +
  scale_shape_manual(values = c(Plant = 21, Pollinator = 24), guide = "none") +
  scale_color_manual(
    values = c(
      "None"     = "grey40",
      "Low"      = "grey70",
      "Moderate" = "orange",
      "High"     = "red3"
    ),
    guide = "none"
  ) +
  coord_fixed(
    ratio = 0.35,
    xlim = c(0, 1.05),
    ylim = c(-0.25, 1.25)
  ) +
  theme_void(base_size = 11) +
  theme(
    plot.subtitle  = element_text(hjust = 0),
    legend.position = "none",
    plot.margin    = margin(5.5, 5.5, 5.5, 25)
  ) +
  labs(
    subtitle = "Intensively managed habitat:\nlow modularity, high connectance & nestedness (hub B1)"
  )

panelA <- (pA_meadow | pA_managed) +
  plot_annotation(
    title = "A. Network structure \u2192 community-level parasite prevalence"
  )

## ============================================================
## PANEL B: Species roles → species parasite prevalence
## Fill = centrality; outline = infection; shapes = roles
## ============================================================

nodes_A <- make_nodes(plants_A, polls_A)

deg_df <- edges_A_meadow |>
  count(pollinator, name = "degree")

nodes_B <- nodes_A |>
  left_join(deg_df, by = c("name" = "pollinator")) |>
  mutate(
    degree = if_else(is.na(degree), 0L, degree),
    # toy d' to distinguish generalists vs specialists
    d_prime = case_when(
      guild == "Plant"      ~ NA_real_,
      degree <= 1           ~ 0.9,
      degree == 2           ~ 0.6,
      degree >= 3           ~ 0.2,
      TRUE                  ~ 0.5
    ),
    infection = case_when(
      guild == "Plant"                 ~ "Low",
      degree >= 3 & d_prime <= 0.3    ~ "High",
      degree == 2                     ~ "Moderate",
      degree <= 1                     ~ "Low",
      TRUE                            ~ "Low"
    ),
    role = case_when(
      guild == "Plant"       ~ "Plant",
      d_prime <= 0.3         ~ "Generalist (low d')",
      d_prime >= 0.7         ~ "Specialist (high d')",
      TRUE                   ~ "Other pollinator"   # diamond
    )
  )

# centrality scaled 0–1 for pollinators (plants set to 0)
nodes_B <- nodes_B |>
  mutate(
    centrality = if_else(
      guild == "Pollinator" & max(degree, na.rm = TRUE) > 0,
      degree / max(degree, na.rm = TRUE),
      0
    ),
    size_value = if_else(guild == "Plant", 3.0, 3.0 + degree)
  )

edges_B_plot <- prep_edges(edges_A_meadow, nodes_B)

pB <- ggplot() +
  geom_segment(
    data = edges_B_plot,
    aes(x = x_pl, y = y_pl, xend = x_bee, yend = y_bee),
    colour = "grey60",
    alpha = 0.5
  ) +
  geom_point(
    data = nodes_B,
    aes(
      x = x, y = y,
      shape = role,
      fill  = centrality,
      colour = infection,
      size  = size_value
    ),
    stroke = 0.8
  ) +
  scale_shape_manual(values = c(
    "Plant"               = 21,
    "Generalist (low d')" = 24,  # triangle
    "Specialist (high d')"= 22,  # square
    "Other pollinator"    = 23   # diamond
  )) +
  scale_fill_gradient(
    low = "grey20",
    high = "deepskyblue2",
    na.value = "grey20",
    name = "Centrality\n(degree, scaled)"
  ) +
  scale_color_manual(
    values = c(
      "None"     = "grey40",
      "Low"      = "grey70",
      "Moderate" = "orange",
      "High"     = "red3"
    ),
    guide = "none"   # no separate parasite legend
  ) +
  scale_size_continuous(range = c(2.5, 6), guide = "none") +
  coord_fixed(
    ratio = 0.35,
    xlim = c(-0.05, 1.05),
    ylim = c(-0.25, 1.25)
  ) +
  theme_void(base_size = 11) +
  theme(
    plot.title    = element_text(hjust = 0, face = "bold"),
    plot.subtitle = element_text(hjust = 0),
    legend.position = "bottom",
    legend.box     = "vertical"
  ) +
  labs(
    title    = "B. Species network roles \u2192 species-level parasite prevalence",
    subtitle = "Nodes filled by centrality (degree): highly central (bright) generalists\nhave higher parasite prevalence (red outlines)",
    shape    = "Node role"
  ) +
  # Make shape legend visible on black
  guides(
    shape = guide_legend(
      override.aes = list(
        fill   = "white",
        colour = "white",
        size   = 4
      )
    )
  )

## ============================================================
## PANEL C: Superspreader overlap with Apis mellifera
## ============================================================

plants_C  <- paste0("P", 1:4)
polls_C   <- c("Apis_mellifera", "W1", "W2", "W3")
nodes_C   <- make_nodes(plants_C, polls_C) |>
  mutate(
    guild2 = case_when(
      name == "Apis_mellifera" ~ "Apis mellifera",
      guild == "Plant"         ~ "Plant",
      TRUE                     ~ "Wild bee"
    )
  )

edges_C <- tribble(
  ~pollinator,      ~plant,
  "Apis_mellifera", "P1",
  "Apis_mellifera", "P2",
  "Apis_mellifera", "P3",
  "Apis_mellifera", "P4",
  "W1",             "P1",
  "W1",             "P2",
  "W2",             "P2",
  "W2",             "P3",
  "W3",             "P4"
)

edges_C_plot <- prep_edges(edges_C, nodes_C)

nodes_C <- nodes_C |>
  mutate(
    overlap = case_when(
      name == "W1" ~ "High overlap with Apis",
      name == "W2" ~ "Moderate overlap with Apis",
      name == "W3" ~ "Low overlap with Apis",
      TRUE         ~ NA_character_
    ),
    infection = case_when(
      name == "Apis_mellifera" ~ "High",
      name == "W1"             ~ "High",
      name == "W2"             ~ "Moderate",
      name == "W3"             ~ "Low",
      TRUE                     ~ "None"
    )
  )

pC <- ggplot() +
  # thicker edges for Apis
  geom_segment(
    data = dplyr::filter(edges_C_plot, pollinator == "Apis_mellifera"),
    aes(x = x_pl, y = y_pl, xend = x_bee, yend = y_bee),
    colour = "grey70",
    alpha  = 0.7,
    linewidth = 1.1
  ) +
  # thinner edges for wild bees
  geom_segment(
    data = dplyr::filter(edges_C_plot, pollinator != "Apis_mellifera"),
    aes(x = x_pl, y = y_pl, xend = x_bee, yend = y_bee),
    colour = "grey50",
    alpha  = 0.7,
    linewidth = 0.6
  ) +
  geom_point(
    data = nodes_C,
    aes(
      x = x, y = y,
      shape = guild2,
      colour = infection,
      size  = guild2
    ),
    fill = "black",
    stroke = 0.7
  ) +
  scale_shape_manual(values = c(
    "Plant"          = 21,
    "Wild bee"       = 24,
    "Apis mellifera" = 22
  ), guide = "none") +
  scale_size_manual(values = c(
    "Plant"          = 3.3,
    "Wild bee"       = 4,
    "Apis mellifera" = 5
  ), guide = "none") +
  scale_color_manual(
    values = c(
      "None"     = "grey40",
      "Low"      = "grey70",
      "Moderate" = "orange",
      "High"     = "red3"
    ),
    guide = "none"
  ) +
  coord_fixed(
    ratio = 0.35,
    xlim = c(-0.05, 1.05),
    ylim = c(-0.25, 1.25)
  ) +
  theme_void(base_size = 11) +
  theme(
    plot.title    = element_text(hjust = 0, face = "bold"),
    plot.subtitle = element_text(hjust = 0),
    legend.position = "none"
  ) +
  labs(
    title    = "C. Superspreader overlap with Apis mellifera",
    subtitle = "Higher floral overlap with abundant Apis mellifera increases infection risk\nfor co-foraging wild bees"
  )

## ============================================================
## Combine all panels and apply black theme globally
## ============================================================

final_fig <- panelA / pB / pC +
  plot_layout(heights = c(1.1, 1, 1))

final_fig <- final_fig & theme_black

final_fig

## Example export
# ggsave(
#   "network_epidemiology_conceptual.png",
#   final_fig,
#   width  = 10,
#   height = 7,
#   dpi    = 300
# )

## Example export
ggsave(
  "network_epidemiology_conceptual.png",
  final_fig,
  width  = 10,
  height = 7,
  dpi    = 300
)


## ============================================================
## legend
## ============================================================

library(ggplot2)
library(cowplot)   # install.packages("cowplot") if needed

## Colors match the main figure
infection_cols <- c(
  "None"     = "grey40",
  "Low"      = "grey70",
  "Moderate" = "orange",
  "High"     = "red3"
)

# Dummy data just to generate the legend
dummy <- data.frame(
  infection = factor(names(infection_cols),
                     levels = c("High", "Moderate", "Low", "None")),
  x = 1,
  y = 1
)

# A throwaway plot whose only purpose is to create the legend
legend_base <- ggplot(dummy, aes(x = x, y = y, colour = infection)) +
  geom_point(size = 4, fill = "black", shape = 21, stroke = 1) +
  scale_color_manual(
    values = infection_cols,
    name   = "Parasite prevalence"
  ) +
  theme_void() +
  theme(
    legend.position   = "bottom",
    legend.background = element_rect(fill = "black", colour = NA),
    legend.key        = element_rect(fill = "black", colour = NA),
    legend.title      = element_text(colour = "white"),
    legend.text       = element_text(colour = "white")
  )

# Extract just the legend and draw it as its own plot
legend_only  <- cowplot::get_legend(legend_base)
legend_plot  <- cowplot::ggdraw(legend_only)

legend_plot

ggsave(
  "parasite_prevalence_legend.png",
  legend_plot,
  width = 7,
  height = 1.5,
  dpi = 300,
  bg = "black"
)



library(ggplot2)
library(dplyr)

## Small bee–flower network icon -----------------------------

nodes_icon <- tibble::tibble(
  id    = c("bee1", "bee2", "bee3", "fl1",  "fl2",  "fl3"),
  type  = c(rep("bee", 3),  rep("flower", 3)),
  x     = c(0.2, 0.5, 0.8,  0.2,  0.5,  0.8),
  y     = c(0.80, 0.95, 0.80,  0.25,  0.10,  0.25),
  label = ifelse(type == "bee", "🐝", "🌸")
)

edges_icon <- tibble::tibble(
  from  = c("bee1", "bee1", "bee2", "bee2", "bee3", "bee3"),
  to    = c("fl1",  "fl2",  "fl2",  "fl3",  "fl1",  "fl3")
) |>
  left_join(nodes_icon |> select(id, x, y),
            by = c("from" = "id")) |>
  rename(x_bee = x, y_bee = y) |>
  left_join(nodes_icon |> select(id, x, y),
            by = c("to" = "id")) |>
  rename(x_fl = x, y_fl = y)

icon_network <- ggplot() +
  # edges
  geom_segment(
    data = edges_icon,
    aes(x = x_bee, y = y_bee, xend = x_fl, yend = y_fl),
    colour = "grey70",
    linewidth = 0.6
  ) +
  # nodes as emoji
  geom_text(
    data = nodes_icon,
    aes(x = x, y = y, label = label),
    size = 6
  ) +
  coord_fixed(xlim = c(0, 1), ylim = c(0, 1)) +
  theme_void(base_size = 10) +
  theme(
    plot.background  = element_rect(fill = "black", colour = NA),
    panel.background = element_rect(fill = "black", colour = NA)
  )

icon_network

ggsave(
  "bee_flower_network_icon.png",
  icon_network,
  width = 2.0,
  height = 2.0,
  dpi = 300,
  bg = "black"
)
