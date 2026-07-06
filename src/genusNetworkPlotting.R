# genusNetworkPlotting_pointsToggle_themesAppended_darkFix.R
# Adds:
#   1) show_points argument to toggle raw points on/off (default TRUE)
#   2) theme argument to select ggthemes-based themes from ggplotThemes.txt
#   3) Automatic dark-theme styling:
#        - prediction line becomes white when theme is talk_black
#        - CI "curves"/bands use lighter greys on black background
#
# The theme functions (theme_ms, theme_talk_black) are appended at the end
# so this file is self-contained.
#
# Intended use:
#   source("src/genusNetworkPlotting_pointsToggle_themesAppended_darkFix.R")
#
# Example:
#   fig <- make_network_figure(..., show_points = FALSE, theme = "talk_black")

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

# Find CE table for a given effect key (e.g., "connectance")
get_effect_df <- function(fit, effect, re_formula = NA) {
  ce <- brms::conditional_effects(fit, re_formula = re_formula)
  keys <- names(ce)

  matches <- grepl(effect, keys, fixed = TRUE)
  if (!any(matches)) {
    stop("Effect '", effect, "' not found among: ", paste(keys, collapse = ", "))
  }

  hit <- keys[which(matches)[1]]
  df  <- ce[[hit]]

  xcol <- if (effect %in% names(df)) {
    effect
  } else {
    nm_matches <- grepl(effect, names(df), fixed = TRUE)
    if (!any(nm_matches)) {
      stop("Could not find x column for effect '", effect, "' in conditional_effects table.")
    }
    names(df)[which(nm_matches)[1]]
  }

  list(df = df, xcol = xcol, key = hit)
}


# Linetype from coefficient: solid if posterior sign prob >= prob, dashed otherwise
linetype_from_coef <- function(fit, coef_term, prob = 0.95) {
  fe <- brms::fixef(fit, probs = c((1-prob)/2, 1-(1-prob)/2))
  row <- if (coef_term %in% rownames(fe)) {
    coef_term
  } else {
    rn <- rownames(fe)
    matches <- grepl(coef_term, rn, fixed = TRUE)
    if (!any(matches)) {
      stop("Coefficient term '", coef_term, "' not found among: ", paste(rn, collapse = ", "))
    }
    rn[which(matches)[1]]
  }
  draws <- try(posterior::as_draws_df(fit), silent = TRUE)
  if (inherits(draws, "try-error"))
    return(if (fe[row, "Q2.5"] <= 0 && fe[row, "Q97.5"] >= 0) "dashed" else "solid")
  bname <- paste0("b_", row)
  if (!bname %in% names(draws))
    return(if (fe[row, "Q2.5"] <= 0 && fe[row, "Q97.5"] >= 0) "dashed" else "solid")
  ppos  <- mean(draws[[bname]] > 0)
  psign <- max(ppos, 1 - ppos)
  if (psign >= prob) "solid" else "dashed"
}

# Prepare raw points: prefers a precomputed proportion if provided; else successes/trials
prep_points_network <- function(dat, x, prop_col = NULL, succ_col, trials_col, group = "ProjectSubProject") {
  if (!is.null(prop_col) && prop_col %in% names(dat)) {
    dat |>
      dplyr::transmute(x = .data[[x]], y = .data[[prop_col]], grp = .data[[group]])
  } else {
    dat |>
      dplyr::transmute(x = .data[[x]], y = .data[[succ_col]] / .data[[trials_col]], grp = .data[[group]])
  }
}

# --- helper to draw layered CIs (gives a dark→light gradient look) ---
add_ci_gradient <- function(p, fit, effect_key, re_formula = NA,
                            probs = c(0.50, 0.70, 0.85, 0.95),
                            fills = rep("grey20", 4),
                            alphas = c(0.35, 0.25, 0.18, 0.12)) {
  stopifnot(length(probs) == length(fills), length(fills) == length(alphas))
  for (i in seq_along(probs)) {
    ce_i <- brms::conditional_effects(fit, re_formula = re_formula, prob = probs[i])
    keys <- names(ce_i)
    matches <- grepl(effect_key, keys, fixed = TRUE)
    if (!any(matches)) {
      stop("Effect '", effect_key, "' not found among: ", paste(keys, collapse = ", "))
    }
    hit  <- keys[which(matches)[1]]
    df_i <- ce_i[[hit]]

    xcol <- if (effect_key %in% names(df_i)) {
      effect_key
    } else {
      nm_matches <- grepl(effect_key, names(df_i), fixed = TRUE)
      if (!any(nm_matches)) {
        stop("Could not find x column for effect '", effect_key, "' in conditional_effects table.")
      }
      names(df_i)[which(nm_matches)[1]]
    }
    p <- p + ggplot2::geom_ribbon(
      data = df_i,
      ggplot2::aes(x = .data[[xcol]], ymin = lower__, ymax = upper__),
      inherit.aes = FALSE,
      fill = fills[i], alpha = alphas[i]
    )
  }
  p
}

# Single panel for a network metric
plot_network_panel <- function(fit, effect_key, coef_term,
                               raw_df, xcol, prop_col, succ_col, trials_col,
                               xlab, ylab,
                               re_formula = NA,
                               show_points = TRUE,
                               theme = "ms") {

  ef <- get_effect_df(fit, effect = effect_key, re_formula = re_formula)
  df <- ef$df; x_in_df <- ef$xcol
  lt <- linetype_from_coef(fit, coef_term, prob = 0.95)


  is_black <- .is_talk_black(theme)

  # Line color fix for black slides
  line_col <- if (is_black) "white" else "black"

  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data[[x_in_df]], y = estimate__)) +
    ggplot2::geom_line(linewidth = 1.6, linetype = lt, colour = line_col)

  if (isTRUE(show_points)) {
    pts <- prep_points_network(
      raw_df, x = xcol, prop_col = prop_col,
      succ_col = succ_col, trials_col = trials_col
    )

    p <- p + ggplot2::geom_point(
      data = pts,
      ggplot2::aes(x = x, y = y, color = grp),
      inherit.aes = FALSE,
      alpha = 0.85
    )
  }

  p <- p +
    viridis::scale_color_viridis(discrete = TRUE) +
    ggplot2::labs(x = xlab, y = ylab, color = "", fill = "Credible interval")

  # Lighter CI bands on black background
  fills_use  <- if (is_black) rep("grey70", 4) else rep("grey20", 4)
  alphas_use <- if (is_black) c(0.28, 0.22, 0.16, 0.10) else c(0.35, 0.25, 0.18, 0.12)

  p <- add_ci_gradient(
    p, fit, effect_key, re_formula = re_formula,
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

# Build figure for a genus & outcome across common network metrics
make_network_figure <- function(fit, raw_df, outcome_label,
                                prop_col, succ_col, trials_col,
                                metrics_spec, # data.frame with effect_key, coef_term, xcol, xlab
                                ncol = 2,
                                show_points = TRUE,
                                theme = "ms") {

  panels <- purrr::pmap(metrics_spec, function(effect_key, coef_term, xcol, xlab) {

    out <- try(
      plot_network_panel(
        fit, effect_key, coef_term,
        raw_df, xcol, prop_col, succ_col, trials_col,
        xlab, ylab = outcome_label,
        re_formula = NA,
        show_points = show_points,
        theme = theme
      ),
      silent = TRUE
    )

    if (inherits(out, "try-error")) {
      msg <- conditionMessage(attr(out, "condition"))
      message("
❌ Panel failed: ", effect_key, " (", coef_term, ")
", msg, "
")
    }

    out
  })

  panels <- Filter(function(x) !inherits(x, "try-error"), panels)
  if (!length(panels)) stop("No panels could be constructed.")

  ggpubr::ggarrange(
    plotlist = panels,
    labels = LETTERS[seq_along(panels)],
    ncol = ncol, nrow = ceiling(length(panels)/ncol),
    common.legend = TRUE, legend = "bottom"
  )
}



# Build a list of ggplot panels for one outcome across metrics (internal utility)
make_network_panels <- function(fit, raw_df, outcome_label,
                                prop_col, succ_col, trials_col,
                                metrics_spec,
                                show_points = TRUE,
                                theme = "ms",
                                ylab_left_only = TRUE) {

  panels <- purrr::map(seq_len(nrow(metrics_spec)), function(i) {
    effect_key <- metrics_spec$effect_key[[i]]
    coef_term  <- metrics_spec$coef_term[[i]]
    xcol       <- metrics_spec$xcol[[i]]
    xlab       <- metrics_spec$xlab[[i]]

    ylab_i <- if (isTRUE(ylab_left_only) && i != 1) "" else outcome_label

    out <- try(
      plot_network_panel(
        fit, effect_key, coef_term,
        raw_df, xcol, prop_col, succ_col, trials_col,
        xlab, ylab = ylab_i,
        re_formula = NA,
        show_points = show_points,
        theme = theme
      ),
      silent = TRUE
    )

    if (inherits(out, "try-error")) {
      msg <- conditionMessage(attr(out, "condition"))
      message("\n❌ Panel failed: ", effect_key, " (", coef_term, ")\n", msg, "\n")
    }

    out
  })

  panels
}


# Manuscript-style figure: stack two parasite outcomes (rows) with network metrics as columns.
# - Panels are labeled A), B), C), ... going left-to-right across the first row, then the second row.
# - By default, only the leftmost panel in each row shows the y-axis label.
make_stacked_parasite_network_figure <- function(
  fit_top, fit_bottom,
  raw_df_top, raw_df_bottom,
  top_outcome_label, bottom_outcome_label,
  top_prop_col, top_succ_col, top_trials_col,
  bottom_prop_col, bottom_succ_col, bottom_trials_col,
  metrics_spec,
  ncol = nrow(metrics_spec),
  show_points = FALSE,
  theme = "ms",
  label_suffix = ")",
  allow_partial = FALSE
) {

  panels_top <- make_network_panels(
    fit = fit_top,
    raw_df = raw_df_top,
    outcome_label = top_outcome_label,
    prop_col = top_prop_col,
    succ_col = top_succ_col,
    trials_col = top_trials_col,
    metrics_spec = metrics_spec,
    show_points = show_points,
    theme = theme,
    ylab_left_only = TRUE
  )

  panels_bottom <- make_network_panels(
    fit = fit_bottom,
    raw_df = raw_df_bottom,
    outcome_label = bottom_outcome_label,
    prop_col = bottom_prop_col,
    succ_col = bottom_succ_col,
    trials_col = bottom_trials_col,
    metrics_spec = metrics_spec,
    show_points = show_points,
    theme = theme,
    ylab_left_only = TRUE
  )

  panels <- c(panels_top, panels_bottom)

  is_err <- vapply(panels, inherits, logical(1), what = "try-error")
  if (any(is_err) && !isTRUE(allow_partial)) {
    stop("One or more panels failed to build (", sum(is_err), " of ", length(panels), "). See printed errors above.")
  }

  panels_ok <- panels[!is_err]
  if (!length(panels_ok)) stop("No panels could be constructed.")

  labels <- paste0(LETTERS[seq_along(panels_ok)], label_suffix)

  ggpubr::ggarrange(
    plotlist = panels_ok,
    labels = labels,
    ncol = ncol,
    nrow = 2,
    common.legend = TRUE,
    legend = "bottom"
  )
}


# =============================================================================
# Appended themes from ggplotThemes.txt
# =============================================================================

theme_talk_black <- function(base_size = 14, base_family = "sans") {
  base <- if (requireNamespace("ggthemes", quietly = TRUE)) {
    ggthemes::theme_foundation(base_size = base_size, base_family = base_family)
  } else {
    ggplot2::theme_minimal(base_size = base_size, base_family = base_family)
  }

  base + ggplot2::theme(
    plot.title = ggplot2::element_text(
      face = "bold", colour = "#ffffb3",
      size = ggplot2::rel(1.2), hjust = 0.5,
      margin = ggplot2::margin(0, 0, 20, 0)
    ),
    text = ggplot2::element_text(),
    panel.background = ggplot2::element_rect(colour = NA, fill = "black"),
    plot.background = ggplot2::element_rect(colour = NA, fill = "black"),
    panel.border = ggplot2::element_rect(colour = NA),
    axis.title = ggplot2::element_text(face = "bold", size = ggplot2::rel(1), colour = "white"),
    axis.title.y = ggplot2::element_text(angle = 90, vjust = 2),
    axis.title.x = ggplot2::element_text(vjust = -0.2),
    axis.text = ggplot2::element_text(colour = "white"),
    axis.line.x = ggplot2::element_line(colour = "white"),
    axis.line.y = ggplot2::element_line(colour = "white"),
    axis.ticks = ggplot2::element_line(colour = "white"),
    panel.grid.major = ggplot2::element_blank(),
    panel.grid.minor = ggplot2::element_blank(),
    legend.background = ggplot2::element_rect(fill = "black"),
    legend.text = ggplot2::element_text(color = "white"),
    legend.key = ggplot2::element_rect(colour = NA, fill = "black"),
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.box = "vertical",
    legend.key.size = grid::unit(0.5, "cm"),
    legend.margin = grid::unit(0, "cm"),
    legend.title = ggplot2::element_text(face = "italic", colour = "white"),
    plot.margin = grid::unit(c(10, 5, 5, 5), "mm"),
    strip.background = ggplot2::element_rect(colour = "#2D3A4C", fill = "black"),
    strip.text = ggplot2::element_text(face = "bold", colour = "white")
  )
}


theme_ms <- function(base_size = 14, base_family = "sans") {
  base <- if (requireNamespace("ggthemes", quietly = TRUE)) {
    ggthemes::theme_foundation(base_size = base_size, base_family = base_family)
  } else {
    ggplot2::theme_minimal(base_size = base_size, base_family = base_family)
  }

  base + ggplot2::theme(
    plot.title = ggplot2::element_text(
      face = "bold", colour = "#ffffb3",
      size = ggplot2::rel(1.2), hjust = 0.5,
      margin = ggplot2::margin(0, 0, 20, 0)
    ),
    text = ggplot2::element_text(),
    panel.background = ggplot2::element_rect(colour = NA, fill = "white"),
    plot.background = ggplot2::element_rect(colour = NA, fill = "white"),
    panel.border = ggplot2::element_rect(colour = NA),
    axis.title = ggplot2::element_text(face = "bold", size = ggplot2::rel(1), colour = "black"),
    axis.title.y = ggplot2::element_text(angle = 90, vjust = 2),
    axis.title.x = ggplot2::element_text(vjust = 0),
    axis.text = ggplot2::element_text(colour = "black"),
    axis.line.x = ggplot2::element_line(colour = "black"),
    axis.line.y = ggplot2::element_line(colour = "black"),
    axis.ticks = ggplot2::element_line(colour = "black"),
    panel.grid.major = ggplot2::element_blank(),
    panel.grid.minor = ggplot2::element_blank(),
    legend.position = "none",
    legend.background = ggplot2::element_rect(fill = "white"),
    legend.text = ggplot2::element_text(color = "black"),
    legend.key = ggplot2::element_rect(colour = NA, fill = "white"),
    legend.direction = "horizontal",
    legend.box = "vertical",
    legend.key.size = grid::unit(0.5, "cm"),
    legend.title = ggplot2::element_text(face = "italic", colour = "black"),
    plot.margin = grid::unit(c(5, 5, 5, 5), "mm"),
    strip.background = ggplot2::element_rect(colour = "#2D3A4C", fill = "white"),
    strip.text = ggplot2::element_text(face = "bold", colour = "black")
  )
}
