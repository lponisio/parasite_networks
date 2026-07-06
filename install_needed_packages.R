## install_needed_packages.R
## Install packages needed for the parasite_networks analysis scripts.
## Run from any R session with:
## source("install_needed_packages.R")

## Use a stable CRAN mirror
options(repos = c(CRAN = "https://cloud.r-project.org"))

## Packages used directly by the current analysis/helper scripts
cran_packages <- c(
  ## General data manipulation / plotting
  "tidyverse",
  "ggplot2",

  ## PCA/correlation scripts
  "corrr",
  "ggcorrplot",
  "FactoMineR",
  "factoextra",

  ## Phylogeny scripts
  "ape",

  ## Model fitting
  "brms",
  "rstan",
  "lme4",
  "glmmTMB",

  ## Model diagnostics / posterior checks
  "performance",
  "DHARMa",
  "bayesplot",
  "posterior",
  "loo",
  "see",

  ## Table/output helpers
  "knitr",
  "kableExtra"
)

install_if_missing <- function(pkgs) {
  installed <- rownames(installed.packages())
  missing <- setdiff(pkgs, installed)

  if (length(missing) == 0) {
    message("All CRAN packages are already installed.")
    return(invisible(TRUE))
  }

  message("Installing missing CRAN packages: ", paste(missing, collapse = ", "))
  install.packages(missing, dependencies = TRUE)

  still_missing <- setdiff(pkgs, rownames(installed.packages()))
  if (length(still_missing) > 0) {
    stop(
      "These packages were not installed successfully: ",
      paste(still_missing, collapse = ", "),
      call. = FALSE
    )
  }

  invisible(TRUE)
}

install_if_missing(cran_packages)

## Optional but recommended for brms on many systems.
## The current scripts will run with the default brms backend, but cmdstanr is
## often easier to use/reproduce than rstan once CmdStan is installed.
install_cmdstanr <- TRUE
install_cmdstan <- FALSE  # Set TRUE if you also want to install CmdStan itself.

if (isTRUE(install_cmdstanr)) {
  if (!requireNamespace("cmdstanr", quietly = TRUE)) {
    message("Installing cmdstanr from the Stan R package repository.")
    install.packages(
      "cmdstanr",
      repos = c("https://mc-stan.org/r-packages/", getOption("repos"))
    )
  }

  if (!requireNamespace("cmdstanr", quietly = TRUE)) {
    warning("cmdstanr did not install successfully. The scripts can still use the default brms backend.")
  }
}

if (isTRUE(install_cmdstan) && requireNamespace("cmdstanr", quietly = TRUE)) {
  message("Installing CmdStan. This can take a while.")
  cmdstanr::install_cmdstan()
}

## Final check: make sure all directly required CRAN packages load
message("Checking required packages...")
load_check <- vapply(
  cran_packages,
  requireNamespace,
  logical(1),
  quietly = TRUE
)

if (!all(load_check)) {
  stop(
    "The following required packages could not be loaded: ",
    paste(names(load_check)[!load_check], collapse = ", "),
    call. = FALSE
  )
}

message("Package installation/check complete.")

## Print session info for reproducibility
sessionInfo()
