library(ggplot2)
library(ggpubr)
rm(list=ls())
setwd("~/")
source("lab_paths.R")
setwd(local.path)
setwd("parasite_networks")


library(tidyverse)

# ============================================================
# 1. Load model summary CSVs only
# ============================================================

table_path <- "saved/tables"

# Only true model-summary CSVs:
# roles:   apis_CrithidiaPresence.csv
#          bombus_ApicystisSpp.csv
# network: network_apis_CrithidiaPresence.csv
#          network_bombus_ApicystisSpp.csv
model_csv_pattern <- "^(network_)?(apis|bombus|melissodes)_(crithidiapresence|apicystisspp)\\.csv$"

csv_files <- list.files(
  table_path,
  pattern = model_csv_pattern,
  full.names = TRUE,
  ignore.case = TRUE
)

if (length(csv_files) == 0) {
  stop("No model-summary CSVs found. Check saved/tables and file names.")
}

model_key <- tibble(
  path = csv_files,
  file = basename(csv_files),
  file_lower = str_to_lower(file)
) %>%
  mutate(
    scope = if_else(str_detect(file_lower, "^network_"), "network", "roles"),
    
    host_plain = case_when(
      str_detect(file_lower, "(^|_)apis_") ~ "Apis",
      str_detect(file_lower, "(^|_)bombus_") ~ "Bombus",
      str_detect(file_lower, "(^|_)melissodes_") ~ "Melissodes",
      TRUE ~ NA_character_
    ),
    
    parasite_plain = case_when(
      str_detect(file_lower, "crithidiapresence") ~ "Crithidia",
      str_detect(file_lower, "apicystisspp") ~ "Apicystis",
      TRUE ~ NA_character_
    ),
    
    host = paste0("\\textit{", host_plain, "}"),
    parasite = paste0("\\textit{", parasite_plain, "}")
  )

cat("\nUsing model-summary CSVs:\n")
print(model_key$file)

ignored_csvs <- setdiff(
  basename(list.files(table_path, pattern = "\\.csv$", full.names = TRUE)),
  model_key$file
)

cat("\nIgnoring non-model/generated CSVs:\n")
print(ignored_csvs)

tables_raw <- map2(
  model_key$path,
  model_key$file,
  ~ read.csv(.x, stringsAsFactors = FALSE) %>%
    mutate(file = .y)
)

all_results <- bind_rows(tables_raw) %>%
  left_join(
    model_key %>%
      select(file, scope, host, parasite, host_plain, parasite_plain),
    by = "file"
  )

# ============================================================
# 2. Support symbols
# ============================================================

make_support <- function(pgt0, plt0) {
  case_when(
    pgt0 >= 0.975 ~ "\\texttt{++}",
    plt0 >= 0.975 ~ "\\texttt{--}",
    pgt0 >= 0.950 ~ "\\texttt{+}",
    plt0 >= 0.950 ~ "\\texttt{-}",
    TRUE ~ ""
  )
}

all_results <- all_results %>%
  mutate(support = make_support(Pgt0, Plt0))
# ============================================================
# 3. General model summary tables
# Combined table for species-role models
# Combined table for network-topology models
# ============================================================

clean_term <- function(x) {
  recode(
    x,
    "Intercept" = "Intercept",
    
    # Species role predictors
    "scalezdegree" = "Degree",
    "scalezweighted.closeness" = "Weighted closeness",
    "scalezweighted.betweenness" = "Weighted betweenness",
    "scalezd" = "Specialization ($d'$)",
    "scalezHBOverlap" = "\\textit{Apis} overlap",
    
    # Network predictors
    "scalezweighted.NODF" = "Weighted NODF",
    "scalezH2" = "H2'",
    "scalezweighted.cluster.coefficient.HL" = "Cluster coefficient",
    "scalenumber.of.species.HL" = "Bee richness",
    "scalenumber.of.species.LL" = "Plant richness",
    
    # Project contrasts
    "ProjectPN" = "PNW Coast Range",
    "ProjectSF" = "PNW Cascades post-fire",
    "ProjectSI" = "Sky Islands",
    "ProjectSubProjectPNMCOAST" = "PNW Coast Range",
    "ProjectSubProjectPNMORMFIRE" = "PNW Cascades post-fire",
    "ProjectSubProjectSI" = "Sky Islands",
    
    .default = x
  )
}

format_num <- function(x, digits = 2) {
  x <- suppressWarnings(as.numeric(x))
  ifelse(is.na(x), "", formatC(x, format = "f", digits = digits))
}

get_first_col <- function(dat, possible_names) {
  nm <- possible_names[possible_names %in% names(dat)]
  
  if (length(nm) == 0) {
    rep(NA_real_, nrow(dat))
  } else {
    dat[[nm[1]]]
  }
}

traditional_results_raw <- all_results

traditional_results_raw$CI_low <- get_first_col(
  traditional_results_raw,
  c("l.95..CI", "l.95.CI", "Q2.5", "X2.5.")
)

traditional_results_raw$CI_high <- get_first_col(
  traditional_results_raw,
  c("u.95..CI", "u.95.CI", "Q97.5", "X97.5.")
)

traditional_results_raw$Rhat_safe <- get_first_col(
  traditional_results_raw,
  c("Rhat")
)

traditional_results <- traditional_results_raw %>%
  mutate(
    Model = if_else(scope == "network", "Network topology", "Species roles"),
    Host = host_plain,
    Parasite = parasite_plain,
    Term = clean_term(X),
    
    Estimate_chr = format_num(Estimate, 2),
    SE_chr = format_num(Est.Error, 2),
    CI_chr = paste0("[", format_num(CI_low, 2), ", ", format_num(CI_high, 2), "]"),
    Pgt0_chr = format_num(Pgt0, 3),
    Plt0_chr = format_num(Plt0, 3),
    Rhat_chr = format_num(Rhat_safe, 2),
    Support = support
  ) %>%
  filter(Term != "Intercept") %>%
  filter(!is.na(Host), !is.na(Parasite)) %>%
  select(
    Model, Host, Parasite, Term,
    Estimate_chr, SE_chr, CI_chr,
    Pgt0_chr, Plt0_chr, Support, Rhat_chr
  )

make_general_model_table <- function(dat, model_type, label, caption) {
  
  dat_use <- dat %>%
    filter(Model == model_type) %>%
    mutate(
      Host = factor(Host, levels = c("Bombus", "Apis", "Melissodes")),
      Parasite = factor(Parasite, levels = c("Crithidia", "Apicystis"))
    ) %>%
    arrange(Host, Parasite, Term)
  
  if (nrow(dat_use) == 0) {
    stop(paste("No rows found for", model_type))
  }
  
  rows <- dat_use %>%
    mutate(
      row = paste0(
        "\\textit{", Host, "} & ",
        "\\textit{", Parasite, "} & ",
        Term, " & ",
        Estimate_chr, " & ",
        SE_chr, " & ",
        CI_chr, " & ",
        Pgt0_chr, " & ",
        Plt0_chr, " & ",
        Support, " & ",
        Rhat_chr, " \\\\"
      )
    ) %>%
    pull(row) %>%
    paste(collapse = "\n")
  
  paste(
    "\\begin{landscape}",
    "\\begin{table}[ht]",
    "\\centering",
    "\\scriptsize",
    "\\setlength{\\tabcolsep}{3.5pt}",
    "\\renewcommand{\\arraystretch}{1.08}",
    paste0("\\caption{", caption, "}"),
    paste0("\\label{tab:", label, "}"),
    "\\begin{tabular}{llp{3.2cm}ccccccc}",
    "\\toprule",
    "Host & Parasite & Predictor & Estimate & SE & 95\\% CI & $P(>0)$ & $P(<0)$ & Support & Rhat \\\\",
    "\\midrule",
    rows,
    "\\bottomrule",
    "\\end{tabular}",
    "\\end{table}",
    "\\end{landscape}",
    sep = "\n"
  )
}

roles_traditional_latex <- make_general_model_table(
  traditional_results,
  model_type = "Species roles",
  label = "species_roles_model_summary",
  caption = "General model summaries for species-role models. Estimates are posterior means with posterior standard errors and 95\\% credible intervals. Support indicates posterior directional support: \\texttt{++} or \\texttt{--} for posterior probability $\\geq$ 0.975, and \\texttt{+} or \\texttt{-} for posterior probability $\\geq$ 0.950."
)

network_traditional_latex <- make_general_model_table(
  traditional_results,
  model_type = "Network topology",
  label = "network_topology_model_summary",
  caption = "General model summaries for network-topology models. Estimates are posterior means with posterior standard errors and 95\\% credible intervals. Support indicates posterior directional support: \\texttt{++} or \\texttt{--} for posterior probability $\\geq$ 0.975, and \\texttt{+} or \\texttt{-} for posterior probability $\\geq$ 0.950."
)

cat(roles_traditional_latex)
cat("\n\n\\clearpage\n\n")
cat(network_traditional_latex)


# ------------------------------------------------------------
# 3. Rename predictors for table columns
# ------------------------------------------------------------

network_terms <- c(
  "scalezweighted.NODF" = "Nestedness",
  "scalezH2" = "Selectivity",
  "scalezweighted.cluster.coefficient.HL" = "Modularity",
  "scalenumber.of.species.HL" = "BeeRichness",
  "scalenumber.of.species.LL" = "PlantRichness",
  "ProjectPN" = "PNWCoast",
  "ProjectSF" = "PNWCascades",
  "ProjectSI" = "SkyIslands"
)

role_terms <- c(
  "scalezdegree" = "Degree",
  "scalezweighted.closeness" = "Closeness",
  "scalezweighted.betweenness" = "Betweenness",
  "scalezd" = "Selectivity_d",
  "scalezHBOverlap" = "ApisOverlap",
  "ProjectPN" = "PNWCoast",
  "ProjectSF" = "PNWCascades",
  "ProjectSI" = "SkyIslands"
)

# ------------------------------------------------------------
# 4. Make wide tables
# ------------------------------------------------------------

make_wide_table <- function(dat, terms) {
  dat %>%
    filter(X %in% names(terms)) %>%
    mutate(term = recode(X, !!!terms)) %>%
    select(host, parasite, term, support) %>%
    distinct(host, parasite, term, .keep_all = TRUE) %>%
    pivot_wider(
      names_from = term,
      values_from = support,
      values_fill = list(support = ""),
      values_fn = list(support = dplyr::first)
    )
}

network_tab <- all_results %>%
  filter(scope == "network") %>%
  make_wide_table(network_terms)

roles_tab <- all_results %>%
  filter(scope == "roles") %>%
  make_wide_table(role_terms)

# ------------------------------------------------------------
# 5. Function to write LaTeX rows
# ------------------------------------------------------------

make_latex_rows <- function(tab, cols) {
  host_order <- c("\\textit{Bombus}", "\\textit{Apis}", "\\textit{Melissodes}")
  
  rows <- c()
  
  for (h in host_order) {
    sub <- tab %>% filter(host == h)
    
    crith <- sub %>% filter(parasite == "\\textit{Crithidia}")
    api   <- sub %>% filter(parasite == "\\textit{Apicystis}")
    
    crith_vals <- sapply(cols, \(cc) ifelse(cc %in% names(crith), crith[[cc]], ""))
    api_vals   <- sapply(cols, \(cc) ifelse(cc %in% names(api), api[[cc]], ""))
    
    rows <- c(
      rows,
      paste0("\\multirow{2}{*}{", h, "} & \\textit{Crithidia}  & ",
             paste(crith_vals, collapse = " & "), " \\\\"),
      "\\cline{2-10}",
      paste0("& \\textit{Apicystis} & ",
             paste(api_vals, collapse = " & "), " \\\\"),
      "\\hline"
    )
  }
  
  paste(rows, collapse = "\n")
}

network_cols <- c(
  "Nestedness", "Selectivity", "Modularity", "BeeRichness", "PlantRichness",
  "PNWCoast", "PNWCascades", "SkyIslands"
)

roles_cols <- c(
  "Degree", "Closeness", "Betweenness", "Selectivity_d", "ApisOverlap",
  "PNWCoast", "PNWCascades", "SkyIslands"
)

network_rows <- make_latex_rows(network_tab, network_cols)
roles_rows   <- make_latex_rows(roles_tab, roles_cols)

cat(network_rows)
cat("\n\n")
cat(roles_rows)




# ============================================================
# Load individual model dataframes only
# ============================================================

saved_path <- "saved"

rdata_files <- list.files(
  saved_path,
  pattern = "^(apis|bombus|melissodes)_.*\\.Rdata$",
  full.names = TRUE
)

load_one_rdata <- function(path) {
  env <- new.env(parent = emptyenv())
  obj_names <- load(path, envir = env)
  
  objs <- mget(obj_names, envir = env)
  df_objs <- objs[vapply(objs, is.data.frame, logical(1))]
  
  df <- df_objs[[which.max(vapply(df_objs, nrow, integer(1)))]]
  
  df %>%
    mutate(source_file = basename(path))
}

model_dat <- bind_rows(lapply(rdata_files, load_one_rdata)) %>%
  mutate(
    host_model = case_when(
      str_detect(source_file, "^apis_") ~ "Apis",
      str_detect(source_file, "^bombus_") ~ "Bombus",
      str_detect(source_file, "^melissodes_") ~ "Melissodes"
    ),
    parasite_model = case_when(
      str_detect(source_file, "CrithidiaPresence") ~ "Crithidia",
      str_detect(source_file, "ApicystisSpp") ~ "Apicystis"
    )
  )



# ============================================================
# Supplement table: stacked study-system effort table
# Project names updated for actual project codes:
# SF = Sunflower Fields
# SI = Sky Islands Meadows
# PN = Harvested Forests of PNW
# HJ = Cascades Meadows
# ============================================================

# ------------------------------------------------------------
# 1. Manual metadata
# ------------------------------------------------------------
# Only fill fields that cannot be calculated from network.metrics.
# Survey time/round and total hours are calculated from SurveyMin.

project_metadata <- tibble::tribble(
  ~ProjectSubProject, ~StudySystem, ~SurveyLayout, ~SurveyMonths,
  
  "HJA",
  "Cascades Meadows",
  "Ten 3 m $\\times$ 3 m plots per site, for 90 $m^2$ total surveyed area per site.",
  "June--September",
  
  "SI",
  "Sky Islands Meadows",
  "Three 50 m $\\times$ 50 m plots per site, for 7500 $m^2$ total surveyed area per site.",
  "July--September",
  
  "SF",
  "California Sunflower Fields",
  "Two 50 m transects per sampled habitat type.Sunflower transects included an additional 30 min rare-species collection not included in richness or abundance calculations",
  "June--August",
  
  "PN-OR-FIRE",
  "Salvage logged, Post-fire Timber: Cascade Range of Oregon",
  "even 32 m transects per site, for a distance of 224 m with a breadth of 2 m on each side per site of surveying",
  "May--September",
  
  "PN-COAST",
  "Harvested Forests: Coast Range of Oregon",
  "even 32 m transects per site, for a distance of 224 m with a breadth of 2 m on each side per site of surveying",
  "May--September",
  
  "PN-CA-FIRE",
  "Salvage logged, Post-fire Timber: Northern California",
  "even 32 m transects per site, for a distance of 224 m with a breadth of 2 m on each side per site of surveying",
  "May--September"
  
)

# ------------------------------------------------------------
# 2. Collapse to one row per site-year-round survey
# ------------------------------------------------------------
# This avoids inflation from genus-level replicated rows.

screened_individuals <- model_dat %>%
  distinct(Project,ProjectSubProject,
           Year,
           SampleRound,
           Site,
           SurveyMin, .keep_all = TRUE)
# ------------------------------------------------------------
# 3. Formatting helper functions
# ------------------------------------------------------------

format_minutes <- function(x) {
  x <- as.numeric(x)
  
  ifelse(
    is.na(x),
    "",
    ifelse(
      x >= 60,
      paste0(round(x / 60, 2), " hrs"),
      paste0(round(x, 1), " min")
    )
  )
}

format_survey_time <- function(x) {
  x <- sort(unique(na.omit(as.numeric(x))))
  
  if (length(x) == 0) {
    ""
  } else if (length(x) == 1) {
    format_minutes(x)
  } else if (length(x) <= 3) {
    paste(format_minutes(x), collapse = ", ")
  } else {
    paste0(format_minutes(min(x)), "--", format_minutes(max(x)))
  }
}

format_total_hours <- function(x) {
  total_hours <- sum(as.numeric(x), na.rm = TRUE) / 60
  
  ifelse(
    is.na(total_hours) | total_hours == 0,
    "",
    as.character(round(total_hours, 1))
  )
}

# ------------------------------------------------------------
# 4. Calculate rounds per year
# ------------------------------------------------------------

rounds_by_year <- model_dat %>%
  group_by(Project, Year,ProjectSubProject) %>%
  summarise(
    n_rounds = n_distinct(SampleRound),
    .groups = "drop"
  ) %>%
  group_by(ProjectSubProject) %>%
  summarise(
    rounds_per_year = ifelse(
      min(n_rounds, na.rm = TRUE) == max(n_rounds, na.rm = TRUE),
      as.character(max(n_rounds, na.rm = TRUE)),
      paste0(min(n_rounds, na.rm = TRUE), "--", max(n_rounds, na.rm = TRUE))
    ),
    .groups = "drop"
  )

# ------------------------------------------------------------
# 5. Summarize effort by project
# ------------------------------------------------------------

project_effort <- model_dat %>%
  group_by(Project,ProjectSubProject) %>%
  summarise(
    n_sites_num = n_distinct(Site),
    
    survey_time_round = format_survey_time(SurveyMin),
    
    years_chr = ifelse(
      n_distinct(Year) == 1,
      as.character(unique(Year)),
      paste0(min(Year, na.rm = TRUE), "--", max(Year, na.rm = TRUE))
    ),
    
    total_surveys_num = n_distinct(paste(Site, Year, SampleRound, sep = "_")),
    
    total_hours_chr = format_total_hours(SurveyMin),
    
    .groups = "drop"
  ) %>%
  left_join(rounds_by_year, by = "ProjectSubProject")

# ------------------------------------------------------------
# 6. Join metadata and format blanks
# ------------------------------------------------------------

supp_effort <- project_metadata %>%
  left_join(project_effort, by = "ProjectSubProject") %>%
  mutate(
    n_sites = ifelse(is.na(n_sites_num), "", as.character(n_sites_num)),
    survey_time_round = ifelse(is.na(survey_time_round), "", survey_time_round),
    rounds_per_year = ifelse(is.na(rounds_per_year), "", rounds_per_year),
    years = ifelse(is.na(years_chr), "", years_chr),
    total_surveys = ifelse(is.na(total_surveys_num), "", as.character(total_surveys_num)),
    total_hours = ifelse(is.na(total_hours_chr), "", total_hours_chr)
  )

# ------------------------------------------------------------
# 7. Function to make each study-system block
# Cleaner spacing + smaller within-project column headers
# ------------------------------------------------------------

make_project_block <- function(x) {
  
  layout_line <- if (x$SurveyLayout == "") {
    ""
  } else {
    paste0(
      "\\multicolumn{7}{@{}p{0.98\\textwidth}@{}}{",
      "\\footnotesize \\textit{Survey layout:} ",
      x$SurveyLayout,
      "} \\\\"
    )
  }
  
  paste(
    "\\midrule",
    paste0("\\multicolumn{7}{@{}l}{\\large\\textbf{", x$StudySystem, "}} \\\\"),
    layout_line,
    "\\addlinespace[0.10cm]",
    
    paste0(
      "\\footnotesize ",
      "\\textbf{Sites} & ",
      "\\textbf{Months} & ",
      "\\textbf{Time/round} & ",
      "\\textbf{Rounds/yr} & ",
      "\\textbf{Years} & ",
      "\\textbf{Surveys} & ",
      "\\textbf{Hours} \\\\"
    ),
    
    paste0(
      x$n_sites, " & ",
      x$SurveyMonths, " & ",
      x$survey_time_round, " & ",
      x$rounds_per_year, " & ",
      x$years, " & ",
      x$total_surveys, " & ",
      x$total_hours, " \\\\"
    ),
    
    "\\addlinespace[0.08cm]",
    sep = "\n"
  )
}

# ------------------------------------------------------------
# 8. Build table body
# ------------------------------------------------------------

table_body <- supp_effort %>%
  split(seq_len(nrow(.))) %>%
  map_chr(make_project_block) %>%
  paste(collapse = "\n")

# ------------------------------------------------------------
# 9. Full LaTeX table
# ------------------------------------------------------------

supp_effort_latex <- paste(
  "\\begin{table}[ht]",
  "\\centering",
  "\\small",
  "\\setlength{\\tabcolsep}{7pt}",
  "\\renewcommand{\\arraystretch}{1.15}",
  "\\caption{Survey design and sampling effort across study systems. Survey layout is summarized below each study-system heading when applicable. Total hours are calculated from summed survey minutes across all site-by-year-by-round surveys.}",
  "\\label{tab:supp_survey_effort}",
  "\\begin{tabular}{@{}",
  ">{\\raggedright\\arraybackslash}p{1.1cm}",
  ">{\\raggedright\\arraybackslash}p{2.1cm}",
  ">{\\raggedright\\arraybackslash}p{1.9cm}",
  ">{\\raggedright\\arraybackslash}p{1.5cm}",
  ">{\\raggedright\\arraybackslash}p{1.7cm}",
  ">{\\raggedright\\arraybackslash}p{1.4cm}",
  ">{\\raggedright\\arraybackslash}p{1.3cm}",
  "@{}}",
  "\\toprule",
  table_body,
  "\\bottomrule",
  "\\end{tabular}",
  "\\end{table}",
  sep = "\n"
)

cat(supp_effort_latex)

# ============================================================
# Supplement table: parasite groups not modeled
# Style: project headers span columns
# Columns: Genus, Crithidia, Apicystis, Nosema bombi, Nosema ceranae
# ============================================================

# ------------------------------------------------------------
# 1. Project labels
# ------------------------------------------------------------

project_labels <- tibble::tribble(
  ~ProjectSubProject, ~StudySystem,
  "SF", "California Sunflower Fields",
  "SI", "Sky Islands Meadows",
  "PN-OR-FIRE", "Salvage logged, Post-fire Timber: Cascade Range of Oregon",
  "PN-COAST", "Harvested Forests: Coast Range of Oregon",  
  "PN-CA-FIRE", "Salvage logged, Post-fire Timber: Northern California",  
  "HJA", "Cascades Meadows"
)

project_order <- c(
  "Cascades Meadows",
  "Sky Islands Meadows",
  "California Sunflower Fields",
  "Salvage logged, Post-fire Timber: Cascade Range of Oregon",
  "Harvested Forests: Coast Range of Oregon",
  "Salvage logged, Post-fire Timber: Northern California" 
)

# ------------------------------------------------------------
# 2. Helper functions
# ------------------------------------------------------------

format_rate <- function(pos, screened, digits = 1) {
  if (is.na(screened) || screened == 0) return("")
  paste0(pos, "/", screened, " (", round(100 * pos / screened, digits), "\\%)")
}

escape_latex <- function(x) {
  x <- as.character(x)
  x <- gsub("&", "\\\\&", x)
  x <- gsub("_", "\\\\_", x)
  x
}

# ============================================================
# A. Bee parasite rates by project × genus
# ============================================================
# Genus-level parasite columns are repeated across species rows,
# so collapse first to one row per project-site-year-round-genus.

genus_events <- model_dat %>%
  distinct(
    ProjectSubProject,
    Site,
    Year,
    SampleRound,
    Genus,
    GenusScreened,
    GenusCrithidiaPresence,
    GenusApicystisSpp,
    GenusNosemaBombi,
    GenusNosemaCeranae
  )

bee_genus_summary <- genus_events %>%
  group_by(ProjectSubProject, Genus) %>%
  summarise(
    screened = sum(GenusScreened, na.rm = TRUE),
    crithidia = sum(GenusCrithidiaPresence, na.rm = TRUE),
    apicystis = sum(GenusApicystisSpp, na.rm = TRUE),
    nosema_bombi = sum(GenusNosemaBombi, na.rm = TRUE),
    nosema_ceranae = sum(GenusNosemaCeranae, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  filter(screened > 0) %>%
  left_join(project_labels, by = "ProjectSubProject") %>%
  rowwise() %>%
  mutate(
    Crithidia = format_rate(crithidia, screened),
    Apicystis = format_rate(apicystis, screened),
    `Nosema bombi` = format_rate(nosema_bombi, screened),
    `Nosema ceranae` = format_rate(nosema_ceranae, screened)
  ) %>%
  ungroup() %>%
  select(
    StudySystem,
    Genus,
    Crithidia,
    Apicystis,
    `Nosema bombi`,
    `Nosema ceranae`
  )

# ============================================================
# B. Syrphid parasite rates from raw sunflower specimen data
# ============================================================
load("../sunflower/data/spec_RBCL_16s.Rdata")

spec_individuals <- spec %>%
  mutate(
    specimen_id = case_when(
      "UniqueID" %in% names(.) & !is.na(UniqueID) ~ as.character(UniqueID),
      "TempID" %in% names(.) & !is.na(TempID) ~ paste0("TempID_", TempID),
      TRUE ~ paste0("row_", row_number())
    )
  ) %>%
  distinct(specimen_id, .keep_all = TRUE)

syrphid_spec <- spec_individuals %>%
  filter(Family == "Syrphidae")

crithidia_cols <- intersect(
  c("CrithidiaSpp", "CrithidiaBombi", "CrithidiaExpoeki"),
  names(syrphid_spec)
)

nosema_bombi_cols <- intersect(
  c("NosemaBombi", "VarimorphaBombi"),
  names(syrphid_spec)
)

nosema_ceranae_cols <- intersect(
  c("NosemaCeranae", "VarimorphaCeranae"),
  names(syrphid_spec)
)

any_screened <- function(dat, cols) {
  if (length(cols) == 0) return(rep(FALSE, nrow(dat)))
  rowSums(!is.na(dat[, cols, drop = FALSE])) > 0
}

any_positive <- function(dat, cols) {
  if (length(cols) == 0) return(rep(FALSE, nrow(dat)))
  rowSums(dat[, cols, drop = FALSE] == 1, na.rm = TRUE) > 0
}

syrphid_parasites <- syrphid_spec %>%
  mutate(
    crithidia_screened = any_screened(., crithidia_cols),
    crithidia_positive = any_positive(., crithidia_cols),
    
    apicystis_screened = !is.na(Apicystis),
    apicystis_positive = Apicystis == 1,
    
    nosema_bombi_screened = any_screened(., nosema_bombi_cols),
    nosema_bombi_positive = any_positive(., nosema_bombi_cols),
    
    nosema_ceranae_screened = any_screened(., nosema_ceranae_cols),
    nosema_ceranae_positive = any_positive(., nosema_ceranae_cols)
  )

syrphid_summary <- tibble(
  StudySystem = "California Sunflower Fields",
  Genus = "Syrphidae",
  Crithidia = format_rate(
    sum(syrphid_parasites$crithidia_positive, na.rm = TRUE),
    sum(syrphid_parasites$crithidia_screened, na.rm = TRUE)
  ),
  Apicystis = format_rate(
    sum(syrphid_parasites$apicystis_positive, na.rm = TRUE),
    sum(syrphid_parasites$apicystis_screened, na.rm = TRUE)
  ),
  `Nosema bombi` = format_rate(
    sum(syrphid_parasites$nosema_bombi_positive, na.rm = TRUE),
    sum(syrphid_parasites$nosema_bombi_screened, na.rm = TRUE)
  ),
  `Nosema ceranae` = format_rate(
    sum(syrphid_parasites$nosema_ceranae_positive, na.rm = TRUE),
    sum(syrphid_parasites$nosema_ceranae_screened, na.rm = TRUE)
  )
)

# ============================================================
# C. Final table data
# ============================================================

parasite_not_modeled_table <- bind_rows(
  bee_genus_summary,
  syrphid_summary
) %>%
  mutate(
    StudySystem = factor(StudySystem, levels = project_order)
  ) %>%
  arrange(StudySystem, Genus) %>%
  mutate(
    StudySystem = as.character(StudySystem)
  )

parasite_not_modeled_table
# ============================================================
# D. Manual LaTeX table with project headers
# ============================================================

make_project_parasite_block <- function(dat, study_system) {
  
  dat_project <- dat %>%
    filter(StudySystem == study_system) %>%
    select(Genus, Crithidia, Apicystis, `Nosema bombi`, `Nosema ceranae`)
  
  if (nrow(dat_project) == 0) return("")
  
  rows <- dat_project %>%
    mutate(
      Genus = escape_latex(Genus)
    ) %>%
    pmap_chr(function(Genus, Crithidia, Apicystis, `Nosema bombi`, `Nosema ceranae`) {
      paste0(
        Genus, " & ",
        Crithidia, " & ",
        Apicystis, " & ",
        `Nosema bombi`, " & ",
        `Nosema ceranae`, " \\\\"
      )
    }) %>%
    paste(collapse = "\n")
  
  paste(
    "\\midrule",
    paste0("\\multicolumn{5}{@{}l}{\\normalsize\\textbf{", study_system, "}} \\\\"),
    "\\addlinespace[0.03cm]",
    "\\scriptsize \\textbf{Genus} & \\textbf{Crithidia} & \\textbf{Apicystis} & \\textbf{Nosema bombi} & \\textbf{Nosema ceranae} \\\\",
    rows,
    "\\addlinespace[0.03cm]",
    sep = "\n"
  )
}
table_body <- project_order %>%
  map_chr(~ make_project_parasite_block(parasite_not_modeled_table, .x)) %>%
  paste(collapse = "\n")

parasite_not_modeled_latex <- paste(
  "\\clearpage",
  "\\begin{table}[ht]",
  "\\centering",
  "\\scriptsize",
  "\\setlength{\\tabcolsep}{5pt}",
  "\\renewcommand{\\arraystretch}{1.02}",
  "\\caption{Parasite groups summarized descriptively but not included in formal models. Parasite summaries are reported by study system and genus using genus-level screening totals. Syrphid fly parasite rates are reported from the raw sunflower specimen data.}",
  "\\label{tab:parasites_not_modeled}",
  "\\begin{tabular}{@{}",
  ">{\\raggedright\\arraybackslash}p{2.4cm}",
  ">{\\raggedright\\arraybackslash}p{2.3cm}",
  ">{\\raggedright\\arraybackslash}p{2.3cm}",
  ">{\\raggedright\\arraybackslash}p{2.5cm}",
  ">{\\raggedright\\arraybackslash}p{2.5cm}",
  "@{}}",
  "\\toprule",
  table_body,
  "\\bottomrule",
  "\\end{tabular}",
  "\\end{table}",
  sep = "\n"
)

cat(parasite_not_modeled_latex)

# ============================================================
# Supplement table: individuals screened by species x project
# Uses SpScreened, not Obs
# ============================================================

species_events <- model_dat %>%
  distinct(
    ProjectSubProject,
    Site,
    Year,
    SampleRound,
    GenusSpecies,
    SpScreened,
    .keep_all = TRUE
  )

screened_by_species <- species_events %>%
  group_by(GenusSpecies, ProjectSubProject) %>%
  summarise(
    n = sum(SpScreened, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  pivot_wider(
    names_from = ProjectSubProject,
    values_from = n,
    values_fill = 0
  )

for (cc in c("HJA", "SI", "PN-COAST", "PN-OR-FIRE", "SF")) {
  if (!cc %in% names(screened_by_species)) screened_by_species[[cc]] <- 0
}

screened_by_species <- screened_by_species %>%
  mutate(
    Total = HJA + SI + `PN-COAST` + `PN-OR-FIRE` + SF,
    Genus = stringr::word(GenusSpecies, 1)
  ) %>%
  filter(Total > 0) %>%
  arrange(factor(Genus, levels = c("Apis", "Bombus", "Melissodes")), GenusSpecies)


species_rows <- screened_by_species %>%
  mutate(
    SpeciesLatex = paste0("\\textit{", GenusSpecies, "}"),
    row = paste0(
      SpeciesLatex, " & ",
      HJA, " & ",
      SI, " & ",
      `PN-COAST`, " & ",
      `PN-OR-FIRE`, " & ",
      SF, " & ",
      Total, " \\\\"
    )
  ) %>%
  group_by(Genus) %>%
  summarise(
    rows = paste(row, collapse = "\n"),
    .groups = "drop"
  ) %>%
  mutate(
    rows = ifelse(
      Genus %in% c("Bombus", "Melissodes"),
      paste0("\\hline\n", rows),
      rows
    )
  ) %>%
  pull(rows) %>%
  paste(collapse = "\n")

totals <- screened_by_species %>%
  summarise(
    HJA = sum(HJA),
    SI = sum(SI),
    PN_COAST = sum(`PN-COAST`),
    PN_OR_FIRE = sum(`PN-OR-FIRE`),
    SF = sum(SF),
    Total = sum(Total)
  )

total_row <- paste0(
  "\\hline\n",
  "\\textbf{Total} & ",
  "\\textbf{", totals$HJA, "} & ",
  "\\textbf{", totals$SI, "} & ",
  "\\textbf{", totals$PN_COAST, "} & ",
  "\\textbf{", totals$PN_OR_FIRE, "} & ",
  "\\textbf{", totals$SF, "} & ",
  "\\textbf{", totals$Total, "} \\\\"
)

screened_species_latex <- paste(
  "\\begin{table}[t]",
  "\\centering",
  "\\scriptsize",
  "\\setlength{\\tabcolsep}{4pt}",
  "\\renewcommand{\\arraystretch}{1.08}",
  "\\caption{Number of individuals screened for each bee species across study systems included in the models. Columns are grouped into natural meadow systems (HJ Andrews = HJA, Sky Islands = SI), Pacific Northwest forest systems (Coast Range, Cascades), and sunflower fields (SF).}",
  "\\label{tab:screened_by_species_projects}",
  "\\begin{tabular}{|>{\\raggedright\\arraybackslash}p{4.8cm}|c|c|c|c|c|c|}",
  "\\hline",
  "\\textbf{Species} &",
  "\\multicolumn{2}{c|}{\\textbf{Natural meadows}} &",
  "\\multicolumn{2}{c|}{\\textbf{PNW harvested forests}} &",
  "\\textbf{SF} &",
  "\\textbf{Total} \\\\",
  "\\hline",
  "& \\textbf{HJA} & \\textbf{SI} &",
  "\\makecell{\\textbf{Coast Range}} &",
  "\\makecell{\\textbf{Cascades}} &",
  "& \\\\",
  "\\hline",
  species_rows,
  total_row,
  "\\hline",
  "\\end{tabular}",
  "\\end{table}",
  sep = "\n"
)

cat(screened_species_latex)
