## Helper functions for posterior coefficient tables with directional
## posterior probabilities.
##
## Source after src/makeBayesR2Table.R because this file reuses
## parse_bayes_R2_model_id() and latex_escape().

## *******************************************************
## Posterior-results tables for manuscript/supplement
## *******************************************************

model_results_dir <- file.path("tables", "model_results")
model_results_csv_file <- file.path(model_results_dir, "model_results.csv")
model_results_tex_file <- file.path(model_results_dir, "model_results.tex")
model_results_txt_file <- file.path(model_results_dir, "model_results_latex_rows.txt")

empty_model_results_table <- function(){
  data.frame(
    ModelID = character(),
    ModelSet = character(),
    Host = character(),
    Parasite = character(),
    Term = character(),
    Estimate = numeric(),
    Est.Error = numeric(),
    Q2.5 = numeric(),
    Q97.5 = numeric(),
    P.gt.0 = numeric(),
    P.lt.0 = numeric(),
    P.direction = numeric(),
    Direction = character(),
    Evidence = character(),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

init_model_results_table <- function(table_file = model_results_csv_file,
                                     reset = FALSE){
  if (!reset && file.exists(table_file)) {
    out <- utils::read.csv(table_file, stringsAsFactors = FALSE,
                           check.names = FALSE)
    expected_cols <- names(empty_model_results_table())
    missing_cols <- setdiff(expected_cols, names(out))
    if (length(missing_cols) > 0) {
      for (cc in missing_cols) out[[cc]] <- NA
    }
    out <- out[, expected_cols]
    return(out)
  }

  empty_model_results_table()
}

posterior_evidence_label <- function(p_direction,
                                     strong = 0.975,
                                     weak = 0.95){
  ifelse(
    is.na(p_direction),
    NA_character_,
    ifelse(p_direction >= strong,
           "strong",
           ifelse(p_direction >= weak, "weak", "little"))
  )
}

clean_model_term <- function(term){
  out <- term

  ## brms often removes punctuation from scale(...) terms in coefficient names.
  out <- gsub("scale\\(([^)]+)\\)", "scaled \\1", out)
  out <- gsub("^scale", "scaled ", out)
  out <- gsub("[()]", "", out)

  replacements <- c(
    "zweighted.NODF" = "weighted NODF",
    "zH2" = "H2",
    "zweighted.cluster.coefficient.HL" = "weighted cluster coefficient",
    "number.of.species.HL" = "higher-level species richness",
    "number.of.species.LL" = "lower-level species richness",
    "zweighted.betweenness" = "weighted betweenness",
    "zweighted.closeness" = "weighted closeness",
    "zd" = "d'",
    "zdegree" = "degree",
    "zHBOverlap" = "honey bee overlap",
    "zproportional.generality" = "proportional generality",
    "znormalised.degree" = "normalized degree"
  )

  for (pat in names(replacements)) {
    out <- gsub(pat, replacements[[pat]], out, fixed = TRUE)
  }

  out <- gsub("^ProjectSubProject", "Project: ", out)
  out <- gsub("_", " ", out, fixed = TRUE)
  out <- gsub("\\.", " ", out)
  out <- gsub("  +", " ", out)
  trimws(out)
}

extract_model_results <- function(model,
                                  model_id,
                                  ci_probs = c(0.025, 0.975),
                                  exclude_intercept = TRUE,
                                  term_labels = NULL){
  fx <- brms::fixef(model, summary = FALSE)
  fx <- as.matrix(fx)

  if (ncol(fx) == 0) {
    return(empty_model_results_table())
  }

  terms <- colnames(fx)
  keep <- rep(TRUE, length(terms))

  if (isTRUE(exclude_intercept)) {
    keep <- keep & !grepl("Intercept", terms)
  }

  terms <- terms[keep]
  fx <- fx[, keep, drop = FALSE]

  if (length(terms) == 0) {
    return(empty_model_results_table())
  }

  meta <- if (exists("parse_bayes_R2_model_id", mode = "function")) {
    parse_bayes_R2_model_id(model_id)
  } else {
    list(ModelSet = NA_character_, Host = NA_character_, Parasite = NA_character_)
  }

  out <- lapply(terms, function(tt){
    vals <- fx[, tt]
    vals <- vals[is.finite(vals)]

    estimate <- mean(vals, na.rm = TRUE)
    est_error <- stats::sd(vals, na.rm = TRUE)
    qs <- stats::quantile(vals, probs = ci_probs, na.rm = TRUE,
                          names = FALSE)
    p_gt_0 <- mean(vals > 0, na.rm = TRUE)
    p_lt_0 <- mean(vals < 0, na.rm = TRUE)
    direction <- ifelse(estimate >= 0, "positive", "negative")
    p_direction <- max(p_gt_0, p_lt_0, na.rm = TRUE)

    label <- if (!is.null(term_labels) && tt %in% names(term_labels)) {
      unname(term_labels[[tt]])
    } else {
      clean_model_term(tt)
    }

    data.frame(
      ModelID = model_id,
      ModelSet = meta$ModelSet,
      Host = meta$Host,
      Parasite = meta$Parasite,
      Term = label,
      Estimate = estimate,
      Est.Error = est_error,
      Q2.5 = qs[1],
      Q97.5 = qs[2],
      P.gt.0 = p_gt_0,
      P.lt.0 = p_lt_0,
      P.direction = p_direction,
      Direction = direction,
      Evidence = posterior_evidence_label(p_direction),
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  })

  do.call(rbind, out)
}

format_model_results_table <- function(results_table, digits = 3){
  if (nrow(results_table) == 0) {
    return(results_table)
  }

  results_table <- results_table[order(results_table$Host,
                                       results_table$Parasite,
                                       results_table$ModelSet,
                                       results_table$Term), ]

  num_cols <- c("Estimate", "Est.Error", "Q2.5", "Q97.5",
                "P.gt.0", "P.lt.0", "P.direction")
  for (cc in intersect(num_cols, names(results_table))) {
    results_table[[cc]] <- round(results_table[[cc]], digits)
  }

  results_table
}

make_result_latex_rows <- function(results_table,
                                   digits = 3,
                                   include_model_cols = TRUE){
  if (nrow(results_table) == 0) {
    return(character(0))
  }

  fmt <- format_model_results_table(results_table, digits = digits)

  rows <- apply(fmt, 1, function(rr){
    base_cols <- if (isTRUE(include_model_cols)) {
      c(latex_escape(rr[["ModelSet"]]),
        latex_escape(rr[["Host"]]),
        latex_escape(rr[["Parasite"]]))
    } else {
      character(0)
    }

    cols <- c(
      base_cols,
      latex_escape(rr[["Term"]]),
      sprintf(paste0("%.", digits, "f"), as.numeric(rr[["Estimate"]])),
      sprintf(paste0("%.", digits, "f"), as.numeric(rr[["Q2.5"]])),
      sprintf(paste0("%.", digits, "f"), as.numeric(rr[["Q97.5"]])),
      sprintf(paste0("%.", digits, "f"), as.numeric(rr[["P.gt.0"]])),
      sprintf(paste0("%.", digits, "f"), as.numeric(rr[["P.lt.0"]])),
      sprintf(paste0("%.", digits, "f"), as.numeric(rr[["P.direction"]])),
      latex_escape(rr[["Evidence"]])
    )

    paste(cols, collapse = " & ")
  })

  paste0(rows, " \\\\")
}

write_model_results_latex <- function(results_table,
                                      tex_file,
                                      digits = 3,
                                      include_model_cols = TRUE,
                                      caption = "Posterior summaries for host-parasite network models.",
                                      label = "tab:model-results"){
  make_dir(dirname(tex_file))

  if (nrow(results_table) == 0) {
    writeLines("% No model-result rows have been recorded yet.", tex_file)
    return(invisible(tex_file))
  }

  rows <- make_result_latex_rows(results_table,
                                 digits = digits,
                                 include_model_cols = include_model_cols)

  if (isTRUE(include_model_cols)) {
    tabular <- "lll l rrrrrrl"
    header <- paste(
      "Model set & Host & Parasite & Term & Estimate & Q2.5 & Q97.5 &",
      "$P(>0)$ & $P(<0)$ & $P(direction)$ & Evidence \\\\",
      sep = " "
    )
  } else {
    tabular <- "l rrrrrrl"
    header <- paste(
      "Term & Estimate & Q2.5 & Q97.5 & $P(>0)$ & $P(<0)$ &",
      "$P(direction)$ & Evidence \\\\",
      sep = " "
    )
  }

  latex_lines <- c(
    "% Generated by src/runHostParasiteModels.R",
    "\\begin{table}[!ht]",
    "\\centering",
    "\\small",
    paste0("\\caption{", latex_escape(caption), "}"),
    paste0("\\label{", latex_escape(label), "}"),
    paste0("\\begin{tabular}{", tabular, "}"),
    "\\hline",
    header,
    "\\hline",
    rows,
    "\\hline",
    "\\end{tabular}",
    "\\end{table}"
  )

  writeLines(latex_lines, tex_file)
  invisible(tex_file)
}

write_model_results_outputs <- function(results_table,
                                        csv_file,
                                        tex_file,
                                        txt_file,
                                        digits = 3,
                                        include_model_cols = TRUE,
                                        caption = "Posterior summaries for host-parasite network models.",
                                        label = "tab:model-results"){
  make_dir(dirname(csv_file))

  utils::write.csv(format_model_results_table(results_table, digits = digits),
                   file = csv_file,
                   row.names = FALSE)

  write_model_results_latex(results_table,
                            tex_file = tex_file,
                            digits = digits,
                            include_model_cols = include_model_cols,
                            caption = caption,
                            label = label)

  rows <- make_result_latex_rows(results_table,
                                 digits = digits,
                                 include_model_cols = include_model_cols)
  writeLines(rows, txt_file)

  invisible(results_table)
}

make_ms_table <- function(model,
                          model_id,
                          out_dir = model_results_dir,
                          digits = 3,
                          ci_probs = c(0.025, 0.975),
                          exclude_intercept = TRUE,
                          term_labels = NULL){
  make_dir(out_dir)

  one_model <- extract_model_results(model = model,
                                     model_id = model_id,
                                     ci_probs = ci_probs,
                                     exclude_intercept = exclude_intercept,
                                     term_labels = term_labels)

  prefix <- file.path(out_dir, model_id)

  write_model_results_outputs(
    one_model,
    csv_file = paste0(prefix, "_results.csv"),
    tex_file = paste0(prefix, "_results.tex"),
    txt_file = paste0(prefix, "_latex_rows.txt"),
    digits = digits,
    include_model_cols = FALSE,
    caption = paste("Posterior summaries for", model_id, "."),
    label = paste0("tab:", gsub("_", "-", model_id), "-results")
  )

  invisible(one_model)
}

## Dot-name alias, because older scripts/functions used dot-style names.
make.ms.table <- make_ms_table

record_model_results <- function(results_table,
                                 model,
                                 model_id,
                                 csv_file = model_results_csv_file,
                                 tex_file = model_results_tex_file,
                                 txt_file = model_results_txt_file,
                                 digits = 3,
                                 ci_probs = c(0.025, 0.975),
                                 exclude_intercept = TRUE,
                                 term_labels = NULL,
                                 caption = "Posterior summaries and directional posterior probabilities for host-parasite network models.",
                                 label = "tab:model-results"){
  if (is.null(results_table)) {
    results_table <- init_model_results_table(table_file = csv_file)
  }

  one_model <- make_ms_table(model = model,
                             model_id = model_id,
                             digits = digits,
                             ci_probs = ci_probs,
                             exclude_intercept = exclude_intercept,
                             term_labels = term_labels)

  results_table <- results_table[results_table$ModelID != model_id, , drop = FALSE]
  results_table <- rbind(results_table, one_model)

  write_model_results_outputs(
    results_table,
    csv_file = csv_file,
    tex_file = tex_file,
    txt_file = txt_file,
    digits = digits,
    include_model_cols = TRUE,
    caption = caption,
    label = label
  )

  invisible(results_table)
}

