# ============================================================
# Survey effort for networks represented in model data
# Each project-site-year-round is counted once across hosts
# and parasite files.
# ============================================================

rm(list = ls())
setwd("~/Dropbox (University of Oregon)/parasite_networks")
library(tidyverse)

# 1. Load network model data ----------------------------------

files <- list.files(
  "saved",
  pattern = "^network_(apis|bombus|melissodes)_ApicystisSpp\\.Rdata$",
  full.names = TRUE,
  ignore.case = TRUE
)

if (!length(files)) stop("No network model-data files found.")

model_dat <- map_dfr(files, function(path) {
  env <- new.env()
  load(path, envir = env)
  
  if (!is.data.frame(env$model_data)) {
    stop("Missing model_data in ", basename(path))
  }
  
  as_tibble(env$model_data) %>%
    mutate(source_file = basename(path))
}) %>%
  mutate(
    across(
      c(ProjectSubProject, Site, SampleRound, UniNetwork),
      ~ trimws(as.character(.x))
    ),
    Year = as.integer(as.character(Year)),
    SurveyMin = as.numeric(as.character(SurveyMin))
  )

# 2. Set reporting labels; retain original project identities --

labels <- c(
  "HJA" = "Cascades Meadows",
  "SI" = "Sky Islands Meadows",
  "SF" = "California Sunflower Fields",
  "PN-COAST" = "Coast Range Harvested Forests",
  "PN-OR-FIRE" = "Oregon Post-fire Harvested Forests",
  "PN-CA-FIRE" = "California Post-fire Harvested Forests"
)

model_dat <- model_dat %>%
  mutate(
    StudySystemCode = ProjectSubProject,
    StudySystem = factor(
      unname(labels[StudySystemCode]),
      levels = unname(labels)
    )
  )

# 3. Validate identifiers and effort before deduplicating ------

key <- c("ProjectSubProject", "Site", "Year", "SampleRound")

if (any(model_dat %>% transmute(
  bad = if_any(
    all_of(c(key, "UniNetwork")),
    ~ is.na(.x) | as.character(.x) == ""
  ) | is.na(StudySystem)
) %>% pull(bad))) {
  stop("Missing identifiers or unknown study-system labels.")
}

identity_map <- model_dat %>%
  distinct(across(all_of(c(key, "UniNetwork"))))

if (
  any(count(identity_map, across(all_of(key)))$n > 1) ||
  any(count(identity_map, ProjectSubProject, UniNetwork)$n > 1)
) {
  stop("Survey and network identifiers are not one-to-one.")
}

conflicts <- model_dat %>%
  group_by(across(all_of(key))) %>%
  summarise(
    effort_values = n_distinct(SurveyMin, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  filter(effort_values > 1)

if (nrow(conflicts)) {
  print(conflicts, n = Inf)
  stop("Conflicting effort values for the same survey.")
}

# Collapse copies, using the available effort value.
survey_events <- model_dat %>%
  group_by(across(all_of(key))) %>%
  summarise(
    UniNetwork = first(UniNetwork),
    StudySystem = first(StudySystem),
    SurveyMin = first(SurveyMin[!is.na(SurveyMin)],
                      default = NA_real_),
    .groups = "drop"
  )

if (any(!is.finite(survey_events$SurveyMin) |
        survey_events$SurveyMin < 0)) {
  stop("Missing or invalid survey effort.")
}

# 4. Summarize unique surveys ---------------------------------

range_text <- function(x) {
  if (min(x) == max(x)) as.character(min(x)) else
    paste0(min(x), "--", max(x))
}

rounds <- survey_events %>%
  count(StudySystem, ProjectSubProject, Site, Year) %>%
  group_by(StudySystem) %>%
  summarise(rounds_per_year = range_text(n), .groups = "drop")

survey_effort <- survey_events %>%
  group_by(StudySystem) %>%
  summarise(
    sites = nrow(distinct(pick(ProjectSubProject, Site))),
    years = paste(sort(unique(Year)), collapse = ", "),
    time_per_round = paste0(range_text(SurveyMin), " min"),
    surveys = n(),
    total_hours = round(sum(SurveyMin) / 60, 1),
    .groups = "drop"
  ) %>%
  left_join(rounds, by = "StudySystem") %>%
  arrange(StudySystem)

print(survey_effort, n = Inf, width = Inf)

# ============================================================
# 5. Network richness and distinct interaction links
# Assumes network rows = plants, columns = pollinators.
# Links count unique plant-pollinator pairs, not total visits.
# ============================================================

metrics <- c(
  "number.of.species.HL",
  "number.of.species.LL",
  "connectance"
)

# Check that repeated copies of each network agree.
metric_conflicts <- model_dat %>%
  group_by(ProjectSubProject, UniNetwork) %>%
  summarise(
    across(all_of(metrics), ~ n_distinct(.x, na.rm = TRUE)),
    .groups = "drop"
  ) %>%
  filter(if_any(all_of(metrics), ~ .x > 1))

if (nrow(metric_conflicts)) {
  print(metric_conflicts, n = Inf)
  stop("Conflicting metrics for the same network.")
}

# Collapse repeated records, retaining available metric values.
network_dat <- model_dat %>%
  group_by(ProjectSubProject, UniNetwork) %>%
  summarise(
    StudySystem = first(StudySystem),
    across(
      all_of(metrics),
      ~ first(.x[!is.na(.x)], default = NA_real_)
    ),
    .groups = "drop"
  ) %>%
  rename(
    pollinator_richness = number.of.species.HL,
    plant_richness = number.of.species.LL
  ) %>%
  mutate(
    interaction_links =
      connectance * pollinator_richness * plant_richness
  )

# Inspect networks with unavailable metrics.
network_dat %>%
  filter(if_any(
    c(plant_richness, pollinator_richness, interaction_links),
    ~ !is.finite(.x)
  )) %>%
  print(n = Inf, width = Inf)

# Treat unavailable values as NA, not zero.
network_dat <- network_dat %>%
  mutate(across(
    c(plant_richness, pollinator_richness, interaction_links),
    ~ replace(.x, !is.finite(.x), NA_real_)
  ))

# Check reconstructed links where values are available.
if (any(
  abs(network_dat$interaction_links -
      round(network_dat$interaction_links)) > 1e-6,
  na.rm = TRUE
)) {
  stop("Non-integer link counts: check raw connectance and richness.")
}

network_dat <- network_dat %>%
  mutate(interaction_links = round(interaction_links))

safe_mean <- function(x) {
  if (all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE)
}

safe_median <- function(x) {
  if (all(is.na(x))) NA_real_ else median(x, na.rm = TRUE)
}

network_summary <- network_dat %>%
  group_by(StudySystem) %>%
  summarise(
    networks = n(),
    across(
      c(plant_richness, pollinator_richness, interaction_links),
      list(
        n = ~ sum(!is.na(.x)),
        mean = safe_mean,
        median = safe_median
      ),
      .names = "{.col}_{.fn}"
    ),
    .groups = "drop"
  ) %>%
  mutate(across(
    ends_with("_mean") | ends_with("_median"),
    ~ round(.x, 1)
  ))

survey_effort_extended <- survey_effort %>%
  left_join(network_summary, by = "StudySystem")

print(survey_effort_extended, n = Inf, width = Inf)







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
  filter(Family == "Syrphidae",
         Pathogen==1)


crithidia_cols <- intersect(
  c("CrithidiaSpp", "CrithidiaBombi", "CrithidiaExpoeki"),
  names(syrphid_spec)
)


nosema_bombi_cols <- intersect(
  c("NosemaBombi"),
  names(syrphid_spec)
)


nosema_ceranae_cols <- intersect(
  c("NosemaCeranae"),
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

format_rate <- function(pos, screened, digits = 1) {
  if (is.na(screened) || screened == 0) return("")
  paste0(pos, "/", screened, " (", round(100 * pos / screened, digits), "\\%)")
}


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
