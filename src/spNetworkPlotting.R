# spNetworkPlotting_pointsToggle_themesAppended_darkFix.R
# Species-level network conditional-effect plotting helpers
#
# Adds the same "talk_black visibility" fixes you requested for genus-level plots:
#   - prediction / fitted line is white when theme = "talk_black"
#   - credible-interval gradient bands use lighter greys on black background
#
# Also supports:
#   - show_points argument to toggle raw point overlay on/off (default TRUE)
#   - theme argument to choose between "ms" and "talk_black"
#
# The theme functions (theme_ms, theme_talk_black) are appended at the end
# so this file is self-contained.
#
# Intended use:
#   source("src/spNetworkPlotting_pointsToggle_themesAppended_darkFix.R")
#
# Example:
#   p <- plot_ce_with_points(bombus.CrithidiaPresence,
#                            effect = "zdegree",
#                            raw_df = sub.bombus,
#                            x = "zdegree",
#                            succ = "SpCrithidiaPresence",
#                            trials = "SpScreened",
#                            xlab = "Degree centrality",
#                            ylab = "Crithidia prevalence",
#                            theme = "talk_black",
#                            show_points = FALSE)

# ---- internal theme resolver ----
.get_theme_fun <- function(theme = "ms") {
  if (is.null(theme)) return(NULL)
  if (is.function(theme)) return(theme)
  if (is.character(theme)) {
    theme <- tolower(theme[1])
    if (theme %in% c("ms", "manuscript", "white")) return(theme_ms)
    if (theme %in% c("talk_black", "black", "talk")) return(theme_talk_black)
    if (theme %in% c("none", "default")) return(NULL)
  }
  NULL
}

.is_talk_black <- function(theme) {
  is.character(theme) && tolower(theme[1]) %in% c("talk_black", "black", "talk")
}

# --- 1) Extract conditional_effects table for a given effect ---
get_effect_df <- function(fit, effect, re_formula = NA) {
  ce <- brms::conditional_effects(fit, re_formula = re_formula)
  keys <- names(ce)
  idx  <- which(grepl(effect, keys, fixed = TRUE))
  if (!length(idx)) {
    stop("Effect '", effect, "' not found in conditional_effects().")
  }
  hit <- keys[idx[1]]
  df  <- ce[[hit]]

  if (!effect %in% names(df)) {
    cand <- names(df)[grepl(effect, names(df), fixed = TRUE)]
    if (length(cand)) effect <- cand[1] else {
      stop("No matching x column found for effect: ", effect)
    }
  }
  list(df = df, xcol = effect)
}

# --- 2) Linetype from posterior sign probability ---
linetype_from_coef <- function(fit, term, prob = 0.95) {
  fe <- brms::fixef(fit, probs = c((1-prob)/2, 1-(1-prob)/2))
  row <- if (term %in% rownames(fe)) term else {
    rn <- rownames(fe); rn[which.max(grepl(term, rn, fixed = TRUE))]
  }
  ci_includes0 <- fe[row, "Q2.5"] <= 0 && fe[row, "Q97.5"] >= 0

  draws <- try(posterior::as_draws_df(fit), silent = TRUE)
  if (inherits(draws, "try-error")) return(if (ci_includes0) "dashed" else "solid")

  bname <- paste0("b_", row)
  if (!bname %in% names(draws)) return(if (ci_includes0) "dashed" else "solid")

  ppos  <- mean(draws[[bname]] > 0)
  psign <- max(ppos, 1-ppos)
  if (psign >= prob) "solid" else "dashed"
}

# --- 3) Raw points helper ---
prep_points <- function(dat, x, succ, trials, group = "ProjectSubProject") {
  dat |>
    dplyr::transmute(
      x    = .data[[x]],
      y    = .data[[succ]] / .data[[trials]],
      grp  = .data[[group]]
    )
}

# --- 4) Layered CI ribbons for a gradient look ---
add_ci_gradient_ce <- function(p, fit, effect, re_formula = NA,
                               probs  = c(0.50, 0.70, 0.85, 0.95),
                               fills  = rep("grey20", 4),
                               alphas = c(0.35, 0.25, 0.18, 0.12)) {
  stopifnot(length(probs) == length(fills), length(fills) == length(alphas))
  for (i in seq_along(probs)) {
    ce_i <- brms::conditional_effects(fit, re_formula = re_formula, prob = probs[i])
    keys <- names(ce_i)
    idx  <- which(grepl(effect, keys, fixed = TRUE))
    if (!length(idx)) next
    hit  <- keys[idx[1]]
    df_i <- ce_i[[hit]]

    xcol_i <- if (effect %in% names(df_i)) effect else {
      cand <- names(df_i)[grepl(effect, names(df_i), fixed = TRUE)]
      if (length(cand)) cand[1] else stop("No x column for effect in CE at prob=", probs[i])
    }

    p <- p + ggplot2::geom_ribbon(
      data = df_i,
      ggplot2::aes(x = .data[[xcol_i]], ymin = lower__, ymax = upper__),
      inherit.aes = FALSE,
      fill = fills[i], alpha = alphas[i]
    )
  }
  p
}

# --- 5) UPDATED: CE plot with optional points + theme + dark fixes ---
plot_ce_with_points <- function(fit, effect, raw_df, x, succ, trials,
                                xlab,
                                ylab = "Crithidia prevalence",
                                re_formula = NA,
                                prob = 0.95,
                                show_points = TRUE,
                                theme = "ms") {

  ef <- get_effect_df(fit, effect, re_formula = re_formula)
  df <- ef$df; x_in_df <- ef$xcol

  # line linetype from posterior
  lt  <- linetype_from_coef(fit, effect, prob = prob)

  is_black <- .is_talk_black(theme)
  line_col <- if (is_black) "white" else "black"

  # base: line (points added conditionally; no single flat ribbon here)
  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data[[x_in_df]], y = estimate__)) +
    ggplot2::geom_line(linewidth = 1.6, linetype = lt, colour = line_col)

  if (isTRUE(show_points)) {
    pts <- prep_points(raw_df, x = x, succ = succ, trials = trials)
    p <- p + ggplot2::geom_point(
      data = pts,
      ggplot2::aes(x = x, y = y, color = grp),
      inherit.aes = FALSE,
      alpha = 0.85
    ) +
      viridis::scale_color_viridis(discrete = TRUE)
  }

  p <- p + ggplot2::labs(x = xlab, y = ylab, color = "", fill = "Credible interval")

  # darker→lighter CI bands, adjusted for black slides
  fills_use  <- if (is_black) rep("grey70", 4) else rep("grey20", 4)
  alphas_use <- if (is_black) c(0.28, 0.22, 0.16, 0.10) else c(0.35, 0.25, 0.18, 0.12)

  p <- add_ci_gradient_ce(
    p, fit, effect, re_formula = re_formula,
    probs  = c(0.50, 0.70, 0.85, 0.95),
    fills  = fills_use,
    alphas = alphas_use
  )

  theme_fun <- .get_theme_fun(theme)
  if (!is.null(theme_fun)) {
    p <- p + theme_fun()
  }

  p
}

# --- 6) Convenience wrapper used in your viz script ---
make_panel <- function(fit, effect, raw_df, x, succ, trials, xlab,
                       ylab = "Apicystis prevalence",
                       show_points = TRUE,
                       theme = "ms") {
  out <- try({
    plot_ce_with_points(
      fit    = fit,
      effect = effect,
      raw_df = raw_df,
      x      = x,
      succ   = succ,
      trials = trials,
      xlab   = xlab,
      ylab   = ylab,
      re_formula = NA,
      show_points = show_points,
      theme  = theme
    )
  }, silent = TRUE)
  if (inherits(out, "try-error")) NULL else out
}

# --- 7) Manuscript helper: stack hosts into rows (one figure per parasite) ---
.blank_panel <- function() {
  ggplot2::ggplot() + ggplot2::theme_void()
}

make_stacked_host_species_network_figure <- function(
    fits_by_host,
    raw_by_host,
    outcome_label,
    succ_col,
    trials_col = "SpScreened",
    metrics_spec = NULL,
    host_order = names(fits_by_host),
    host_labels = host_order,
    theme = "ms",
    show_points = FALSE,
    add_col_headers = TRUE,
    add_row_labels = TRUE,
    label_style = c("paren", "plain"),
    verbose = FALSE,
    common_legend = TRUE,
    legend = "bottom"
) {
  label_style <- match.arg(label_style)

  if (is.null(metrics_spec)) {
    metrics_spec <- tibble::tibble(
      effect = c("zweighted.betweenness", "zweighted.closeness", "zdegree", "zd", "zHBOverlap"),
      x      = c("zweighted.betweenness", "zweighted.closeness", "zdegree", "zd", "zHBOverlap"),
      xlab   = c("Betweenness", "Closeness", "Degree", "d", "Apis overlap")
    )
  }

  n_hosts   <- length(host_order)
  n_metrics <- nrow(metrics_spec)

  # Build panels in row-major order: host1 across metrics, host2 across metrics, ...
  panels <- list()
  for (h in host_order) {
    fit_h <- fits_by_host[[h]]
    raw_h <- raw_by_host[[h]]

    for (i in seq_len(n_metrics)) {
      eff  <- metrics_spec$effect[i]
      xcol <- metrics_spec$x[i]
      # If we add column headers, keep x-axis titles blank; otherwise put metric name on the x-axis.
      xlab <- if (isTRUE(add_col_headers)) "" else metrics_spec$xlab[i]

      p <- try(
        plot_ce_with_points(
          fit    = fit_h,
          effect = eff,
          raw_df = raw_h,
          x      = xcol,
          succ   = succ_col,
          trials = trials_col,
          xlab   = xlab,
          ylab   = outcome_label,
          re_formula = NA,
          show_points = show_points,
          theme  = theme
        ) + {
          if (isTRUE(add_col_headers)) ggplot2::theme(axis.title.x = ggplot2::element_blank()) else ggplot2::theme()
        },
        silent = TRUE
      )

      if (inherits(p, "try-error")) {
        if (isTRUE(verbose)) {
          msg <- conditionMessage(attr(p, "condition"))
          message("\n❌ Host=", h, " effect=", eff, " failed:\n", msg, "\n")
        }
        p <- .blank_panel()
      }
      panels[[length(panels) + 1]] <- p
    }
  }

  # Panel labels A), B), ... along rows
  n_panels <- length(panels)
  labs <- LETTERS[seq_len(min(n_panels, length(LETTERS)))]
  if (n_panels > length(LETTERS)) {
    # fall back for >26 (unlikely here)
    labs <- c(labs, paste0("A", seq_len(n_panels - length(LETTERS))))
  }
  if (label_style == "paren") labs <- paste0(labs, ")")

  main_grid <- ggpubr::ggarrange(
    plotlist = panels,
    ncol = n_metrics, nrow = n_hosts,
    labels = labs[seq_len(n_panels)],
    font.label = list(face = "bold", size = 12),
    common.legend = common_legend,
    legend = legend
  )

  # Optional row labels (hosts)
  if (isTRUE(add_row_labels)) {
    row_grobs <- lapply(host_labels, function(h) {
      ggpubr::text_grob(h, rot = 90, face = "bold", size = 12)
    })
    left_col <- ggpubr::ggarrange(plotlist = row_grobs, ncol = 1, nrow = n_hosts)
    main_grid <- ggpubr::ggarrange(left_col, main_grid, ncol = 2, widths = c(0.05, 1))
  }

  # Optional column headers (metrics)
  if (isTRUE(add_col_headers)) {
    col_grobs <- lapply(metrics_spec$xlab, function(x) ggpubr::text_grob(x, face = "bold", size = 12))
    top_row <- ggpubr::ggarrange(plotlist = col_grobs, ncol = n_metrics, nrow = 1)
    main_grid <- ggpubr::ggarrange(top_row, main_grid, ncol = 1, heights = c(0.08, 1))
  }

  main_grid
}



# --- 8) Manuscript helper: stack parasites into rows (one figure per host) ---
make_stacked_parasite_species_network_figure <- function(
    fits_by_parasite,
    raw_df,
    outcome_labels,
    succ_cols,
    trials_cols = "SpScreened",
    metrics_spec = NULL,
    parasite_order = names(fits_by_parasite),
    parasite_labels = parasite_order,
    theme = "ms",
    show_points = FALSE,
    add_col_headers = TRUE,
    add_row_labels = TRUE,
    label_style = c("paren", "plain"),
    verbose = FALSE,
    common_legend = TRUE,
    legend = "bottom"
) {
  label_style <- match.arg(label_style)

  if (is.null(metrics_spec)) {
    metrics_spec <- tibble::tibble(
      effect = c("zweighted.betweenness", "zweighted.closeness", "zdegree", "zd", "zHBOverlap"),
      x      = c("zweighted.betweenness", "zweighted.closeness", "zdegree", "zd", "zHBOverlap"),
      xlab   = c("Betweenness", "Closeness", "Degree", "d", "Apis overlap")
    )
  }

  n_parasites <- length(parasite_order)
  n_metrics   <- nrow(metrics_spec)

  # normalize mapping inputs
  .get_map <- function(x, key, default = NULL) {
    if (is.list(x) || is.vector(x)) {
      if (!is.null(names(x)) && key %in% names(x)) return(x[[key]])
    }
    if (length(x) == 1) return(x[[1]])
    default
  }

  panels <- list()
  for (pname in parasite_order) {
    fit_p <- fits_by_parasite[[pname]]
    succ_p <- .get_map(succ_cols, pname, default = NULL)
    if (is.null(succ_p)) stop("No succ_col provided for parasite '", pname, "'.")
    trials_p <- .get_map(trials_cols, pname, default = trials_cols)
    ylab_p <- .get_map(outcome_labels, pname, default = as.character(outcome_labels[1]))

    for (i in seq_len(n_metrics)) {
      eff  <- metrics_spec$effect[i]
      xcol <- metrics_spec$x[i]
      xlab <- if (isTRUE(add_col_headers)) "" else metrics_spec$xlab[i]

      g <- try(
        plot_ce_with_points(
          fit    = fit_p,
          effect = eff,
          raw_df = raw_df,
          x      = xcol,
          succ   = succ_p,
          trials = trials_p,
          xlab   = xlab,
          ylab   = ylab_p,
          re_formula = NA,
          show_points = show_points,
          theme  = theme
        ) + {
          if (isTRUE(add_col_headers)) ggplot2::theme(axis.title.x = ggplot2::element_blank()) else ggplot2::theme()
        },
        silent = TRUE
      )

      if (inherits(g, "try-error")) {
        if (isTRUE(verbose)) {
          msg <- conditionMessage(attr(g, "condition"))
          message("\n❌ Parasite=", pname, " effect=", eff, " failed:\n", msg, "\n")
        }
        g <- .blank_panel()
      }

      panels[[length(panels) + 1]] <- g
    }
  }

  # Panel labels A), B), C)... along rows
  n_panels <- length(panels)
  labs <- LETTERS[seq_len(min(n_panels, length(LETTERS)))]
  if (n_panels > length(LETTERS)) {
    labs <- c(labs, paste0("A", seq_len(n_panels - length(LETTERS))))
  }
  if (label_style == "paren") labs <- paste0(labs, ")")

  main_grid <- ggpubr::ggarrange(
    plotlist = panels,
    ncol = n_metrics, nrow = n_parasites,
    labels = labs[seq_len(n_panels)],
    font.label = list(face = "bold", size = 12),
    common.legend = common_legend,
    legend = legend
  )

  # Optional row labels (parasites)
  if (isTRUE(add_row_labels)) {
    row_grobs <- lapply(parasite_labels, function(lbl) {
      ggpubr::text_grob(lbl, rot = 90, face = "bold", size = 12)
    })
    left_col <- ggpubr::ggarrange(plotlist = row_grobs, ncol = 1, nrow = n_parasites)
    main_grid <- ggpubr::ggarrange(left_col, main_grid, ncol = 2, widths = c(0.06, 1))
  }

  # Optional column headers (metrics)
  if (isTRUE(add_col_headers)) {
    col_grobs <- lapply(metrics_spec$xlab, function(x) ggpubr::text_grob(x, face = "bold", size = 12))
    top_row <- ggpubr::ggarrange(plotlist = col_grobs, ncol = n_metrics, nrow = 1)
    main_grid <- ggpubr::ggarrange(top_row, main_grid, ncol = 1, heights = c(0.08, 1))
  }

  main_grid
}

# =============================================================================
# Appended themes from ggplotThemes.txt (self-contained)
# =============================================================================

theme_talk_black <- function(base_size=14, base_family="sans") {
  # Robust: uses ggthemes if available; otherwise falls back to ggplot2 theming.
  if (!requireNamespace("ggthemes", quietly = TRUE)) {
    return(
      ggplot2::theme_minimal(base_size = base_size, base_family = base_family) +
        ggplot2::theme(
          plot.title = ggplot2::element_text(face = "bold", colour = "#ffffb3",
                                             size = ggplot2::rel(1.2), hjust = 0.5,
                                             margin = ggplot2::margin(0,0,20,0)),
          panel.background = ggplot2::element_rect(colour = NA, fill = "black"),
          plot.background  = ggplot2::element_rect(colour = NA, fill = "black"),
          axis.title = ggplot2::element_text(face = "bold", size = ggplot2::rel(1), colour = "white"),
          axis.text  = ggplot2::element_text(colour = "white"),
          axis.line  = ggplot2::element_line(colour = "white"),
          axis.ticks = ggplot2::element_line(colour = "white"),
          panel.grid.major = ggplot2::element_blank(),
          panel.grid.minor = ggplot2::element_blank(),
          legend.background = ggplot2::element_rect(fill = "black"),
          legend.text = ggplot2::element_text(color = "white"),
          legend.key  = ggplot2::element_rect(colour = NA, fill = "black"),
          legend.position = "bottom",
          legend.direction = "horizontal",
          legend.key.size = grid::unit(0.5, "cm"),
          legend.margin = grid::unit(0, "cm"),
          legend.title = ggplot2::element_text(face = "italic", colour = "white"),
          plot.margin = grid::unit(c(10,5,5,5), "mm"),
          strip.background = ggplot2::element_rect(colour = "#2D3A4C", fill = "black"),
          strip.text = ggplot2::element_text(face = "bold", colour = "white")
        )
    )
  }

  ggthemes::theme_foundation(base_size = base_size, base_family = base_family) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", colour = "#ffffb3",
                                         size = ggplot2::rel(1.2), hjust = 0.5,
                                         margin = ggplot2::margin(0,0,20,0)),
      panel.background = ggplot2::element_rect(colour = NA, fill = "black"),
      plot.background  = ggplot2::element_rect(colour = NA, fill = "black"),
      panel.border = ggplot2::element_rect(colour = NA),
      axis.title = ggplot2::element_text(face = "bold", size = ggplot2::rel(1), colour = "white"),
      axis.title.y = ggplot2::element_text(angle = 90, vjust = 2),
      axis.title.x = ggplot2::element_text(vjust = -0.2),
      axis.text = ggplot2::element_text(colour = "white"),
      axis.line.x = ggplot2::element_line(colour = "white"),
      axis.line.y = ggplot2::element_line(colour = "white"),
      axis.ticks  = ggplot2::element_line(colour = "white"),
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      legend.background = ggplot2::element_rect(fill = "black"),
      legend.text = ggplot2::element_text(color = "white"),
      legend.key = ggplot2::element_rect(colour = NA, fill = "black"),
      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.box = "vetical",
      legend.key.size = grid::unit(0.5, "cm"),
      legend.margin = grid::unit(0, "cm"),
      legend.title = ggplot2::element_text(face = "italic", colour = "white"),
      plot.margin = grid::unit(c(10,5,5,5), "mm"),
      strip.background = ggplot2::element_rect(colour = "#2D3A4C", fill = "black"),
      strip.text = ggplot2::element_text(face = "bold", colour = "white")
    )
}

theme_ms <- function(base_size=14, base_family="sans") {
  # Robust: uses ggthemes if available; otherwise falls back to ggplot2 theming.
  if (!requireNamespace("ggthemes", quietly = TRUE)) {
    return(
      ggplot2::theme_minimal(base_size = base_size, base_family = base_family) +
        ggplot2::theme(
          plot.title = ggplot2::element_text(face = "bold", colour = "#000000",
                                             size = ggplot2::rel(1.2), hjust = 0.5,
                                             margin = ggplot2::margin(0,0,20,0)),
          panel.background = ggplot2::element_rect(colour = NA, fill = "white"),
          plot.background  = ggplot2::element_rect(colour = NA, fill = "white"),
          axis.title = ggplot2::element_text(face = "bold", size = ggplot2::rel(1), colour = "black"),
          axis.text  = ggplot2::element_text(colour = "black"),
          axis.line  = ggplot2::element_line(colour = "black"),
          axis.ticks = ggplot2::element_line(colour = "black"),
          panel.grid.major = ggplot2::element_blank(),
          panel.grid.minor = ggplot2::element_blank(),
          legend.position = "none",
          legend.background = ggplot2::element_rect(fill = "white"),
          legend.text = ggplot2::element_text(color = "black"),
          legend.key  = ggplot2::element_rect(colour = NA, fill = "white"),
          legend.direction = "horizontal",
          legend.box = "vetical",
          legend.key.size = grid::unit(0.5, "cm"),
          legend.title = ggplot2::element_text(face = "italic", colour = "black"),
          plot.margin = grid::unit(c(5,5,5,5), "mm"),
          strip.background = ggplot2::element_rect(colour = "#2D3A4C", fill = "white"),
          strip.text = ggplot2::element_text(face = "bold", colour = "black")
        )
    )
  }

  ggthemes::theme_foundation(base_size = base_size, base_family = base_family) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", colour = "#000000",
                                         size = ggplot2::rel(1.2), hjust = 0.5,
                                         margin = ggplot2::margin(0,0,20,0)),
      panel.background = ggplot2::element_rect(colour = NA, fill = "white"),
      plot.background  = ggplot2::element_rect(colour = NA, fill = "white"),
      panel.border = ggplot2::element_rect(colour = NA),
      axis.title = ggplot2::element_text(face = "bold", size = ggplot2::rel(1), colour = "black"),
      axis.title.y = ggplot2::element_text(angle = 90, vjust = 2),
      axis.title.x = ggplot2::element_text(vjust = 0),
      axis.text = ggplot2::element_text(colour = "black"),
      axis.line.x = ggplot2::element_line(colour = "black"),
      axis.line.y = ggplot2::element_line(colour = "black"),
      axis.ticks  = ggplot2::element_line(colour = "black"),
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      legend.position = "none",
      legend.background = ggplot2::element_rect(fill = "white"),
      legend.text = ggplot2::element_text(color = "black"),
      legend.key = ggplot2::element_rect(colour = NA, fill = "white"),
      legend.direction = "horizontal",
      legend.box = "vetical",
      legend.key.size = grid::unit(0.5, "cm"),
      legend.title = ggplot2::element_text(face = "italic", colour = "black"),
      plot.margin = grid::unit(c(5,5,5,5), "mm"),
      strip.background = ggplot2::element_rect(colour = "#2D3A4C", fill = "white"),
      strip.text = ggplot2::element_text(face = "bold", colour = "black")
    )
}
