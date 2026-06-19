library(ggplot2)
library(ggpubr)
rm(list=ls())
setwd("~/")
source("lab_paths.R")
setwd(local.path)
setwd("parasite_networks")


load(file="data/network_mets.RData")
load(file="data/sp_mets.RData")
load(file="../sunflower/data/spec_RBCL_16s.Rdata")

# ============================================================
# Supplement table: stacked study-system effort table
# Project names updated for actual project codes:
# SF = Sunflower Fields
# SI = Sky Islands Meadows
# PN = Harvested Forests of PNW
# HJ = Cascades Meadows
# ============================================================

library(dplyr)
library(purrr)

# ------------------------------------------------------------
# 1. Manual metadata
# ------------------------------------------------------------
# Only fill fields that cannot be calculated from network.metrics.
# Survey time/round and total hours are calculated from SurveyMin.

project_metadata <- tibble::tribble(
  ~Project, ~StudySystem, ~SurveyLayout, ~SurveyMonths,
  
  "HJ",
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
  
  "PN",
  "Harvested Forests of the Pacific Northwest",
  "even 32 m transects per site, for a distance of 224 m with a breadth of 2 m on each side per site of surveying",
  "May--September"
)

# ------------------------------------------------------------
# 2. Collapse to one row per site-year-round survey
# ------------------------------------------------------------
# This avoids inflation from genus-level replicated rows.

survey_events <- network.metrics %>%
  distinct(
    Project,
    Year,
    SampleRound,
    Site,
    SurveyMin
  )

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

rounds_by_year <- survey_events %>%
  group_by(Project, Year) %>%
  summarise(
    n_rounds = n_distinct(SampleRound),
    .groups = "drop"
  ) %>%
  group_by(Project) %>%
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

project_effort <- survey_events %>%
  group_by(Project) %>%
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
  left_join(rounds_by_year, by = "Project")

# ------------------------------------------------------------
# 6. Join metadata and format blanks
# ------------------------------------------------------------

supp_effort <- project_metadata %>%
  left_join(project_effort, by = "Project") %>%
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


writeLines(
  supp_effort_latex,
  "figures/pub_tables/supp_survey_effort_stacked.tex"
)

# ============================================================
# Supplement table: parasite groups not modeled
# Style: project headers span columns
# Columns: Genus, Crithidia, Apicystis, Nosema bombi, Nosema ceranae
# ============================================================

library(dplyr)
library(tidyr)
library(purrr)

# ------------------------------------------------------------
# 1. Project labels
# ------------------------------------------------------------

project_labels <- tibble::tribble(
  ~Project, ~StudySystem,
  "SF", "California Sunflower Fields",
  "SI", "Sky Islands Meadows",
  "PN", "Harvested Forests of the Pacific Northwest",
  "HJ", "Cascades Meadows"
)

project_order <- c(
  "Cascades Meadows",
  "Sky Islands Meadows",
  "California Sunflower Fields",
  "Harvested Forests of the Pacific Northwest"
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

genus_events <- sp.network.metrics %>%
  distinct(
    Project,
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
  group_by(Project, Genus) %>%
  summarise(
    screened = sum(GenusScreened, na.rm = TRUE),
    crithidia = sum(GenusCrithidiaPresence, na.rm = TRUE),
    apicystis = sum(GenusApicystisSpp, na.rm = TRUE),
    nosema_bombi = sum(GenusNosemaBombi, na.rm = TRUE),
    nosema_ceranae = sum(GenusNosemaCeranae, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  filter(screened > 0) %>%
  left_join(project_labels, by = "Project") %>%
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

syrphid_parasites <- syrphid_spec %>%
  mutate(
    crithidia_screened = if (length(crithidia_cols) > 0) {
      rowSums(!is.na(across(all_of(crithidia_cols)))) > 0
    } else {
      FALSE
    },
    
    crithidia_positive = if (length(crithidia_cols) > 0) {
      rowSums(across(all_of(crithidia_cols), ~ .x == 1), na.rm = TRUE) > 0
    } else {
      FALSE
    },
    
    apicystis_screened = !is.na(Apicystis),
    apicystis_positive = Apicystis == 1
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
  `Nosema bombi` = "",
  `Nosema ceranae` = ""
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

dir.create("figures/pub_tables", recursive = TRUE, showWarnings = FALSE)

writeLines(
  parasite_not_modeled_latex,
  "figures/pub_tables/supp_parasites_not_modeled.tex"
)