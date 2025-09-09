
# Find CE table for a given effect key (e.g., "connectance")
get_effect_df <- function(fit, effect, re_formula = NA) {
  ce <- conditional_effects(fit, re_formula = re_formula)
  keys <- names(ce)
  hit  <- keys[which.max(grepl(effect, keys, fixed = TRUE))]
  if (!length(hit) || is.na(hit)) stop("Effect '", effect, "' not found.")
  df <- ce[[hit]]
  list(df = df, xcol = if (effect %in% names(df)) effect else names(df)[grepl(effect, names(df), fixed = TRUE)][1])
}

# Linetype from coefficient: solid if posterior sign prob >= prob, dashed otherwise
linetype_from_coef <- function(fit, coef_term, prob = 0.95) {
  fe <- fixef(fit, probs = c((1-prob)/2, 1-(1-prob)/2))
  row <- if (coef_term %in% rownames(fe)) coef_term else {
    rn <- rownames(fe); rn[which.max(grepl(coef_term, rn, fixed = TRUE))]
  }
  draws <- try(posterior::as_draws_df(fit), silent = TRUE)
  if (inherits(draws, "try-error")) return(if (fe[row, "Q2.5"] <= 0 && fe[row, "Q97.5"] >= 0) "dashed" else "solid")
  bname <- paste0("b_", row)
  if (!bname %in% names(draws)) return(if (fe[row, "Q2.5"] <= 0 && fe[row, "Q97.5"] >= 0) "dashed" else "solid")
  ppos  <- mean(draws[[bname]] > 0)
  psign <- max(ppos, 1 - ppos)
  if (psign >= prob) "solid" else "dashed"
}

# Prepare raw points: prefers a precomputed proportion if provided; else successes/trials
prep_points_network <- function(dat, x, prop_col = NULL, succ_col, trials_col, group = "ProjectSubProject") {
  if (!is.null(prop_col) && prop_col %in% names(dat)) {
    dat |>
      transmute(x = .data[[x]], y = .data[[prop_col]], grp = .data[[group]])
  } else {
    dat |>
      transmute(x = .data[[x]], y = .data[[succ_col]] / .data[[trials_col]], grp = .data[[group]])
  }
}

# Single panel for a network metric
plot_network_panel <- function(fit, effect_key, coef_term,
                               raw_df, xcol, prop_col, succ_col, trials_col,
                               xlab, ylab, re_formula = NA) {
  ef <- get_effect_df(fit, effect = effect_key, re_formula = re_formula)
  df <- ef$df; x_in_df <- ef$xcol
  lt <- linetype_from_coef(fit, coef_term, prob = 0.95)
  pts <- prep_points_network(raw_df, x = xcol, prop_col = prop_col,
                             succ_col = succ_col, trials_col = trials_col)

  # base plot: line & points (no single ribbon here)
  p <- ggplot(df, aes(x = .data[[x_in_df]], y = estimate__)) +
    geom_line(linewidth = 1.5, linetype = lt) +
    geom_point(data = pts, aes(x = x, y = y, color = grp), inherit.aes = FALSE) +
    scale_color_viridis(discrete = TRUE) +
    theme_ms() +
    labs(x = xlab, y = ylab, color = "", fill = "Credible interval")

  # add layered CI bands for gradient effect
  p <- add_ci_gradient(p, fit, effect_key, re_formula = re_formula,
                       probs  = c(0.50, 0.70, 0.85, 0.95),
                       fills  = rep("grey20", 4),
                       alphas = c(0.35, 0.25, 0.18, 0.12))

  p
}

# Build figure for a genus & outcome across common network metrics
make_network_figure <- function(fit, raw_df, outcome_label,
                                prop_col, succ_col, trials_col,
                                metrics_spec, # data.frame with effect_key, coef_term, xcol, xlab
                                ncol = 2) {
  panels <- purrr::pmap(metrics_spec, function(effect_key, coef_term, xcol, xlab) {
    try(
      plot_network_panel(fit, effect_key, coef_term,
                         raw_df, xcol, prop_col, succ_col, trials_col,
                         xlab, ylab = outcome_label, re_formula = NA),
      silent = TRUE
    )
  })
  panels <- Filter(function(x) !inherits(x, "try-error"), panels)

  if (!length(panels)) stop("No panels could be constructed.")
  ggarrange(plotlist = panels,
            labels = LETTERS[seq_along(panels)],
            ncol = ncol, nrow = ceiling(length(panels)/ncol),
            common.legend = TRUE, legend = "bottom")
}

# --- helper to draw layered CIs (gives a dark→light gradient look) ---
add_ci_gradient <- function(p, fit, effect_key, re_formula = NA,
                            probs = c(0.50, 0.70, 0.85, 0.95),
                            fills = rep("grey20", 4),
                            alphas = c(0.35, 0.25, 0.18, 0.12)) {
  stopifnot(length(probs) == length(fills), length(fills) == length(alphas))
  # We'll reuse the CE grid the same way you already do, but with different probs
  for (i in seq_along(probs)) {
    ce_i <- conditional_effects(fit, re_formula = re_formula, prob = probs[i])
    keys <- names(ce_i)
    hit  <- keys[which.max(grepl(effect_key, keys, fixed = TRUE))]
    df_i <- ce_i[[hit]]
    xcol <- if (effect_key %in% names(df_i)) effect_key else names(df_i)[grepl(effect_key, names(df_i), fixed = TRUE)][1]
    p <- p + geom_ribbon(data = df_i,
                         aes(x = .data[[xcol]], ymin = lower__, ymax = upper__),
                         inherit.aes = FALSE,
                         fill = fills[i], alpha = alphas[i])
  }
  p
}
