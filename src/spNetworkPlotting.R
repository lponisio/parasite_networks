# --- 1) Your extractor (kept as-is) ---
get_effect_df <- function(fit, effect, re_formula = NA) {
  ce <- conditional_effects(fit, re_formula = re_formula)
  keys <- names(ce)
  hit  <- keys[which.max(grepl(effect, keys, fixed = TRUE))]
  if (length(hit) == 0 || is.na(hit)) stop("Effect '", effect, "' not found in conditional_effects().")
  df <- ce[[hit]]

  if (!effect %in% names(df)) {
    cand <- names(df)[grepl(effect, names(df), fixed = TRUE)]
    if (length(cand)) effect <- cand[1] else stop("No matching x column found for effect: ", effect)
  }
  list(df = df, xcol = effect)
}

# --- 2) Linetype from posterior sign prob (kept, with a small guard) ---
linetype_from_coef <- function(fit, term, prob = 0.95) {
  fe <- fixef(fit, probs = c((1-prob)/2, 1-(1-prob)/2))
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

# --- 3) Raw points helper (kept) ---
prep_points <- function(dat, x, succ, trials, group = "ProjectSubProject") {
  dat |>
    transmute(
      x    = .data[[x]],
      y    = .data[[succ]] / .data[[trials]],
      grp  = .data[[group]]
    )
}

# --- 4) NEW: layered CI ribbons for a gradient look ---
add_ci_gradient_ce <- function(p, fit, effect, re_formula = NA,
                               probs  = c(0.50, 0.70, 0.85, 0.95),
                               fills  = rep("grey20", 4),
                               alphas = c(0.35, 0.25, 0.18, 0.12)) {
  stopifnot(length(probs) == length(fills), length(fills) == length(alphas))
  for (i in seq_along(probs)) {
    ce_i <- conditional_effects(fit, re_formula = re_formula, prob = probs[i])
    keys <- names(ce_i)
    hit  <- keys[which.max(grepl(effect, keys, fixed = TRUE))]
    if (!length(hit) || is.na(hit)) next
    df_i <- ce_i[[hit]]

    xcol_i <- if (effect %in% names(df_i)) effect else {
      cand <- names(df_i)[grepl(effect, names(df_i), fixed = TRUE)]
      if (length(cand)) cand[1] else stop("No x column for effect in CE at prob=", probs[i])
    }

    p <- p + geom_ribbon(
      data = df_i,
      aes(x = .data[[xcol_i]], ymin = lower__, ymax = upper__),
      inherit.aes = FALSE,
      fill = fills[i], alpha = alphas[i]
    )
  }
  p
}

# --- 5) UPDATED: species-level CE plot with gradient shading ---
plot_ce_with_points <- function(fit, effect, raw_df, x, succ, trials,
                                xlab, ylab = "Crithidia prevalence",
                                re_formula = NA, prob = 0.95) {
  ef <- get_effect_df(fit, effect, re_formula = re_formula)
  df <- ef$df; x_in_df <- ef$xcol

  # line linetype from posterior
  lt  <- linetype_from_coef(fit, effect, prob = prob)

  # raw points
  pts <- prep_points(raw_df, x = x, succ = succ, trials = trials)

  # base: line + points (no single flat ribbon here)
  p <- ggplot(df, aes(x = .data[[x_in_df]], y = estimate__)) +
    geom_line(linewidth = 1.5, linetype = lt) +
    geom_point(data = pts, aes(x = x, y = y, color = grp), inherit.aes = FALSE) +
    scale_color_viridis(discrete = TRUE) +
    theme_ms() +
    labs(x = xlab, y = ylab, color = "", fill = "Credible interval")

  # layered CI bands for a dark→light fade
  p <- add_ci_gradient_ce(
    p, fit, effect, re_formula = re_formula,
    probs  = c(0.50, 0.70, 0.85, 0.95),
    fills  = rep("grey20", 4),
    alphas = c(0.35, 0.25, 0.18, 0.12)
  )

  p
}

# --- 6) Your make_panel() stays the same (it calls plot_ce_with_points) ---
make_panel <- function(fit, effect, raw_df, x, succ, trials, xlab) {
  out <- try({
    plot_ce_with_points(
      fit    = fit,
      effect = effect,
      raw_df = raw_df,
      x      = x,
      succ   = succ,
      trials = trials,
      xlab   = xlab,
      ylab   = "Apicystis prevalence",   # or set per-outcome
      re_formula = NA
    )
  }, silent = TRUE)
  if (inherits(out, "try-error")) NULL else out
}
