# ============================================================
# Descriptive supplementary tables
#   1. Survey effort and total hours
#   2. Individuals screened by species
#   3. Parasite screening rates
#
# Post-fire datasets from Oregon and California are combined.
# ============================================================

library(tidyverse)

saved_path <- "saved"

# ============================================================
# 1. Load the current model data frames
# ============================================================

rdata_files <- list.files(
  saved_path,
  pattern = "^(apis|bombus|melissodes)_.*\\.Rdata$",
  full.names = TRUE,
  ignore.case = TRUE
)

if (length(rdata_files) == 0) {
  stop("No host model-data Rdata files were found in saved/.")
}

load_one_rdata <- function(path) {
  
  env <- new.env(parent = emptyenv())
  object_names <- load(path, envir = env)
  objects <- mget(object_names, envir = env)
  
  data_objects <- objects[
    vapply(objects, is.data.frame, logical(1))
  ]
  
  if (length(data_objects) == 0) {
    stop("No data frame found in: ", basename(path))
  }
  
  # Retain the largest data frame in each saved file.
  dat <- data_objects[[which.max(vapply(data_objects, nrow, integer(1)))]]
  
  dat %>%
    mutate(source_file = basename(path))
}

model_dat <- map_dfr(rdata_files, load_one_rdata)


# ============================================================
# 2. Standardize study-system labels
# ============================================================

model_dat <- model_dat %>%
  mutate(
    StudySystemCode = case_when(
      ProjectSubProject == "HJA" ~ "HJA",
      ProjectSubProject == "SI" ~ "SI",
      ProjectSubProject == "SF" ~ "SF",
      ProjectSubProject == "PN-COAST" ~ "PN-COAST",
      
      # Combine Oregon and California post-fire datasets.
      ProjectSubProject %in% c(
        "PN-OR-FIRE",
        "PN-CA-FIRE"
      ) ~ "PN-FIRE",
      
      TRUE ~ as.character(ProjectSubProject)
    ),
    
    StudySystem = recode(
      StudySystemCode,
      "HJA" = "Cascades Meadows",
      "SI" = "Sky Islands Meadows",
      "SF" = "California Sunflower Fields",
      "PN-COAST" = "Coast Range Harvested Forests",
      "PN-FIRE" = "Post-fire Harvested Forests"
    )
  )

study_system_order <- c(
  "Cascades Meadows",
  "Sky Islands Meadows",
  "California Sunflower Fields",
  "Coast Range Harvested Forests",
  "Post-fire Harvested Forests"
)

model_dat <- model_dat %>%
  mutate(
    StudySystem = factor(
      StudySystem,
      levels = study_system_order
    )
  )

# ============================================================
# 3. Survey effort and total hours
# ============================================================

survey_events <- model_dat %>%
  mutate(
    # Preserve project identity so sites with the same label
    # in different datasets are not merged.
    SiteID = interaction(
      ProjectSubProject,
      Site,
      drop = TRUE
    ),
    
    SurveyID = interaction(
      ProjectSubProject,
      Site,
      Year,
      SampleRound,
      drop = TRUE
    )
  ) %>%
  distinct(
    SurveyID,
    StudySystem,
    StudySystemCode,
    SiteID,
    Year,
    SampleRound,
    SurveyMin,
    .keep_all = TRUE
  )

format_minutes <- function(x) {
  
  x <- sort(unique(na.omit(as.numeric(x))))
  
  if (length(x) == 0) {
    ""
  } else if (length(x) == 1) {
    if (x >= 60) {
      paste0(round(x / 60, 2), " hr")
    } else {
      paste0(round(x, 1), " min")
    }
  } else {
    paste0(
      round(min(x), 1),
      "--",
      round(max(x), 1),
      " min"
    )
  }
}

rounds_by_year <- survey_events %>%
  group_by(StudySystem, Year) %>%
  summarise(
    rounds = n_distinct(SampleRound),
    .groups = "drop"
  ) %>%
  group_by(StudySystem) %>%
  summarise(
    rounds_per_year = case_when(
      n() == 0 ~ "",
      min(rounds) == max(rounds) ~ as.character(max(rounds)),
      TRUE ~ paste0(min(rounds), "--", max(rounds))
    ),
    .groups = "drop"
  )

survey_effort <- survey_events %>%
  group_by(StudySystem) %>%
  summarise(
    sites = n_distinct(SiteID),
    
    years = if_else(
      n_distinct(Year) == 1,
      as.character(first(Year)),
      paste0(
        min(Year, na.rm = TRUE),
        "--",
        max(Year, na.rm = TRUE)
      )
    ),
    
    time_per_round = format_minutes(SurveyMin),
    
    surveys = n_distinct(SurveyID),
    
    total_hours = round(
      sum(as.numeric(SurveyMin), na.rm = TRUE) / 60,
      1
    ),
    
    .groups = "drop"
  ) %>%
  left_join(
    rounds_by_year,
    by = "StudySystem"
  ) %>%
  arrange(StudySystem)

survey_effort

# ============================================================
# Combined table:
# individuals screened by species and parasite rates
# ============================================================

first_existing_column <- function(dat, candidates) {
  hit <- candidates[candidates %in% names(dat)]
  
  if (length(hit) == 0) {
    return(NULL)
  }
  
  hit[[1]]
}

format_rate <- function(positive, screened, digits = 1) {
  if (is.na(screened) || screened == 0) {
    return("")
  }
  
  paste0(
    positive,
    "/",
    screened,
    " (",
    round(100 * positive / screened, digits),
    "\\%)"
  )
}

crithidia_positive_col <- first_existing_column(
  model_dat,
  c(
    "SpCrithidiaPresence",
    "SpeciesCrithidiaPresence",
    "CrithidiaPresence"
  )
)

apicystis_positive_col <- first_existing_column(
  model_dat,
  c(
    "SpApicystisSpp",
    "SpeciesApicystisSpp",
    "ApicystisSpp"
  )
)

varimorpha_bombi_positive_col <- first_existing_column(
  model_dat,
  c(
    "SpVarimorphaBombi",
    "SpNosemaBombi",
    "SpeciesVarimorphaBombi",
    "SpeciesNosemaBombi"
  )
)

varimorpha_ceranae_positive_col <- first_existing_column(
  model_dat,
  c(
    "SpVarimorphaCeranae",
    "SpNosemaCeranae",
    "SpeciesVarimorphaCeranae",
    "SpeciesNosemaCeranae"
  )
)

list(
  Crithidia = crithidia_positive_col,
  Apicystis = apicystis_positive_col,
  Varimorpha_bombi = varimorpha_bombi_positive_col,
  Varimorpha_ceranae = varimorpha_ceranae_positive_col
)

species_events <- model_dat %>%
  mutate(
    SpeciesEventID = interaction(
      ProjectSubProject,
      Site,
      Year,
      SampleRound,
      GenusSpecies,
      drop = TRUE
    )
  ) %>%
  transmute(
    SpeciesEventID,
    StudySystem,
    GenusSpecies,
    SpScreened,
    
    CrithidiaPositive =
      .data[[crithidia_positive_col]],
    
    ApicystisPositive =
      .data[[apicystis_positive_col]],
    
    VarimorphaBombiPositive =
      .data[[varimorpha_bombi_positive_col]],
    
    VarimorphaCeranaePositive =
      .data[[varimorpha_ceranae_positive_col]]
  ) %>%
  distinct(
    SpeciesEventID,
    .keep_all = TRUE
  )

screened_by_system <- species_events %>%
  group_by(
    GenusSpecies,
    StudySystem
  ) %>%
  summarise(
    screened = sum(SpScreened, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    StudySystem = as.character(StudySystem)
  ) %>%
  pivot_wider(
    names_from = StudySystem,
    values_from = screened,
    values_fill = 0
  )

expected_systems <- c(
  "Cascades Meadows",
  "Sky Islands Meadows",
  "California Sunflower Fields",
  "Coast Range Harvested Forests",
  "Post-fire Harvested Forests"
)

for (system_name in expected_systems) {
  if (!system_name %in% names(screened_by_system)) {
    screened_by_system[[system_name]] <- 0
  }
}


species_parasite_rates <- species_events %>%
  group_by(GenusSpecies) %>%
  summarise(
    Screened = sum(SpScreened, na.rm = TRUE),
    
    CrithidiaPositive =
      sum(CrithidiaPositive, na.rm = TRUE),
    
    ApicystisPositive =
      sum(ApicystisPositive, na.rm = TRUE),
    
    VarimorphaBombiPositive =
      sum(VarimorphaBombiPositive, na.rm = TRUE),
    
    VarimorphaCeranaePositive =
      sum(VarimorphaCeranaePositive, na.rm = TRUE),
    
    .groups = "drop"
  ) %>%
  rowwise() %>%
  mutate(
    Crithidia = format_rate(
      CrithidiaPositive,
      Screened
    ),
    
    Apicystis = format_rate(
      ApicystisPositive,
      Screened
    ),
    
    `Varimorpha bombi` = format_rate(
      VarimorphaBombiPositive,
      Screened
    ),
    
    `Varimorpha ceranae` = format_rate(
      VarimorphaCeranaePositive,
      Screened
    )
  ) %>%
  ungroup() %>%
  select(
    GenusSpecies,
    Crithidia,
    Apicystis,
    `Varimorpha bombi`,
    `Varimorpha ceranae`
  )


combined_species_summary <- screened_by_system %>%
  mutate(
    Total = rowSums(
      across(all_of(expected_systems)),
      na.rm = TRUE
    ),
    
    Genus = stringr::word(GenusSpecies, 1)
  ) %>%
  left_join(
    species_parasite_rates,
    by = "GenusSpecies"
  ) %>%
  filter(Total > 0) %>%
  arrange(
    factor(
      Genus,
      levels = c("Apis", "Bombus", "Melissodes")
    ),
    GenusSpecies
  ) %>%
  select(
    GenusSpecies,
    all_of(expected_systems),
    Total,
    Crithidia,
    Apicystis,
    `Varimorpha bombi`,
    `Varimorpha ceranae`
  )

combined_species_summary

# ============================================================
# LaTeX table:
# species × study-system totals × parasite positives
# ============================================================

latex_escape_species <- function(x) {
  paste0("\\textit{", gsub("_", "\\\\_", x), "}")
}

table_dat <- combined_species_summary %>%
  transmute(
    Species = latex_escape_species(GenusSpecies),
    
    HJA = `Cascades Meadows`,
    SI = `Sky Islands Meadows`,
    SF = `California Sunflower Fields`,
    Coast = `Coast Range Harvested Forests`,
    PostFire = `Post-fire Harvested Forests`,
    Total = Total,
    
    Crithidia = Crithidia,
    Apicystis = Apicystis,
    V_bombi = `Varimorpha bombi`,
    V_ceranae = `Varimorpha ceranae`,
    
    Genus = stringr::word(GenusSpecies, 1)
  )

species_rows <- table_dat %>%
  group_by(Genus) %>%
  summarise(
    rows = paste0(
      Species, " & ",
      HJA, " & ",
      SI, " & ",
      SF, " & ",
      Coast, " & ",
      PostFire, " & ",
      Total, " & ",
      Crithidia, " & ",
      Apicystis, " & ",
      V_bombi, " & ",
      V_ceranae, " \\\\",
      collapse = "\n"
    ),
    .groups = "drop"
  ) %>%
  arrange(
    factor(
      Genus,
      levels = c("Apis", "Bombus", "Melissodes")
    )
  ) %>%
  mutate(
    rows = if_else(
      row_number() == 1,
      rows,
      paste0("\\midrule\n", rows)
    )
  ) %>%
  pull(rows) %>%
  paste(collapse = "\n")

table_totals <- table_dat %>%
  summarise(
    HJA = sum(HJA, na.rm = TRUE),
    SI = sum(SI, na.rm = TRUE),
    SF = sum(SF, na.rm = TRUE),
    Coast = sum(Coast, na.rm = TRUE),
    PostFire = sum(PostFire, na.rm = TRUE),
    Total = sum(Total, na.rm = TRUE)
  )

total_row <- paste0(
  "\\midrule\n",
  "\\textbf{Total} & ",
  "\\textbf{", table_totals$HJA, "} & ",
  "\\textbf{", table_totals$SI, "} & ",
  "\\textbf{", table_totals$SF, "} & ",
  "\\textbf{", table_totals$Coast, "} & ",
  "\\textbf{", table_totals$PostFire, "} & ",
  "\\textbf{", table_totals$Total, "} & ",
  " &  &  &  \\\\"
)

combined_species_latex <- paste(
  "\\begin{landscape}",
  "\\begin{table}[ht]",
  "\\centering",
  "\\scriptsize",
  "\\setlength{\\tabcolsep}{3.5pt}",
  "\\renewcommand{\\arraystretch}{1.08}",
  
  "\\caption{Numbers of focal bee specimens screened across study systems and parasite detection rates by species. Columns under ``Total collected per project'' report the number of individuals screened in each study system. Parasite columns report the total number of positive individuals divided by the total number screened, with prevalence in parentheses. Oregon and California post-fire sites were combined into one post-fire harvested-forest category.}",
  
  "\\label{tab:species_screening_parasite_rates}",
  
  "\\begin{tabular}{@{}",
  ">{\\raggedright\\arraybackslash}p{3.4cm}",
  "rrrrrr",
  ">{\\raggedright\\arraybackslash}p{2.1cm}",
  ">{\\raggedright\\arraybackslash}p{2.1cm}",
  ">{\\raggedright\\arraybackslash}p{2.1cm}",
  ">{\\raggedright\\arraybackslash}p{2.1cm}",
  "@{}}",
  
  "\\toprule",
  
  "\\multicolumn{1}{c}{\\textbf{Species}} &",
  "\\multicolumn{6}{c}{\\textbf{Total collected per project}} &",
  "\\multicolumn{4}{c}{\\textbf{Total positive (\\%)}} \\\\",
  
  "\\cmidrule(lr){2-7}",
  "\\cmidrule(lr){8-11}",
  
  "&",
  "\\makecell{\\textbf{Cascades}\\\\\\textbf{Meadows}} &",
  "\\makecell{\\textbf{Sky}\\\\\\textbf{Islands}} &",
  "\\makecell{\\textbf{Sunflower}\\\\\\textbf{Fields}} &",
  "\\makecell{\\textbf{Coast Range}\\\\\\textbf{Forests}} &",
  "\\makecell{\\textbf{Post-fire}\\\\\\textbf{Forests}} &",
  "\\textbf{Total} &",
  "\\textit{Crithidia} &",
  "\\textit{Apicystis} &",
  "\\makecell{\\textit{V. bombi}} &",
  "\\makecell{\\textit{V. ceranae}} \\\\",
  
  "\\midrule",
  species_rows,
  total_row,
  "\\bottomrule",
  
  "\\end{tabular}",
  "\\end{table}",
  "\\end{landscape}",
  sep = "\n"
)

cat(combined_species_latex)