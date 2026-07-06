## Plot all plant-pollinator networks as Sankey/alluvial plots
## Input: all.nets = named list of adjacency matrices
##        rows = plants, columns = pollinators
## Output: PDF with 4 networks per page

library(ggplot2)
library(dplyr)
library(tidyr)
library(purrr)
library(ggalluvial)
library(gridExtra)
library(grid)

## ---- Settings ----

plots_per_page <- 4
nrow_page <- 2
ncol_page <- 2

## Optional: set to TRUE if you want to remove species with zero interactions
drop_zero_rows_cols <- TRUE


## ---- Helper function: convert matrix to edge list ----

network_to_edges <- function(mat) {
  
  mat <- as.matrix(mat)
  
  if (drop_zero_rows_cols) {
    mat <- mat[rowSums(mat, na.rm = TRUE) > 0, , drop = FALSE]
    mat <- mat[, colSums(mat, na.rm = TRUE) > 0, drop = FALSE]
  }
  
  if (nrow(mat) == 0 || ncol(mat) == 0) {
    return(NULL)
  }
  
  edges <- as.data.frame(as.table(mat), stringsAsFactors = FALSE)
  
  ## Robustly rename the first two columns, regardless of whether they are
  ## called Var1/Var2 or inherit dimnames like site/spp
  names(edges)[1:3] <- c("plant", "pollinator", "weight")
  
  edges <- edges %>%
    mutate(
      weight = as.numeric(weight),
      plant = as.character(plant),
      pollinator = as.character(pollinator)
    ) %>%
    filter(!is.na(weight), weight > 0)
  
  if (nrow(edges) == 0) return(NULL)
  
  edges
}

## ---- Helper function: make one Sankey plot ----

plot_network_sankey <- function(mat, net_name = "") {
  
  edges <- network_to_edges(mat)
  
  if (is.null(edges)) {
    return(
      ggplot() +
        theme_void() +
        labs(title = paste0(net_name, "\nNo interactions"))
    )
  }
  
  total_interactions <- sum(edges$weight)
  n_plants <- n_distinct(edges$plant)
  n_pollinators <- n_distinct(edges$pollinator)
  
  ggplot(
    edges,
    aes(
      axis1 = plant,
      axis2 = pollinator,
      y = weight
    )
  ) +
    geom_alluvium(
      aes(fill = plant),
      width = 1 / 12,
      alpha = 0.75,
      color = "gray70",
      linewidth = 0.15
    ) +
    geom_stratum(
      width = 1 / 10,
      fill = "gray95",
      color = "gray40",
      linewidth = 0.25
    ) +
    geom_text(
      stat = "stratum",
      aes(label = after_stat(stratum)),
      size = 2
    ) +
    scale_x_discrete(
      limits = c("Plants", "Pollinators"),
      expand = c(0.08, 0.08)
    ) +
    labs(
      title = net_name,
      subtitle = paste0(
        n_plants, " plants; ",
        n_pollinators, " pollinators; ",
        total_interactions, " interactions"
      ),
      x = NULL,
      y = "Interaction frequency"
    ) +
    theme_minimal(base_size = 9) +
    theme(
      legend.position = "none",
      panel.grid = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks = element_blank(),
      plot.title = element_text(face = "bold", size = 10),
      plot.subtitle = element_text(size = 8),
      axis.text.x = element_text(face = "bold")
    )
}


## ---- Build plots ----

if (is.null(names(all.nets))) {
  names(all.nets) <- paste0("network_", seq_along(all.nets))
}

network_plots <- map2(
  all.nets,
  names(all.nets),
  ~ plot_network_sankey(.x, .y)
)


## ---- Write PDF, 4 plots per page ----

pdf(out_file, width = 14, height = 10)

for (i in seq(1, length(network_plots), by = plots_per_page)) {
  
  plot_chunk <- network_plots[i:min(i + plots_per_page - 1, length(network_plots))]
  
  grid.arrange(
    grobs = plot_chunk,
    nrow = nrow_page,
    ncol = ncol_page
  )
}

dev.off()

message("Saved Sankey plots to: ", normalizePath(out_file))
