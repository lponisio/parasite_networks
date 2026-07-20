## Shared helper functions and model settings for host-parasite brms workflows.
##
## Source this file after:
##   src/runHostParasiteModels.R
##
## It keeps the main model scripts focused on:
##   1. loading/preparing data,
##   2. defining model_specs,
##   3. fitting and/or post-checking models.

## *******************************************************
## Model-fitting tweaks for sparse/problematic binomial models
## *******************************************************

## These priors are intentionally modest regularizers.
priors_standard <- c(
  brms::prior(normal(0, 1), class = "b"),
  brms::prior(student_t(3, 0, 2.5), class = "Intercept"),
  brms::prior(exponential(1), class = "sd")
)

## Slightly stronger regularization for sparse Melissodes models.
priors_sparse <- c(
  brms::prior(normal(0, 0.75), class = "b"),
  brms::prior(student_t(3, 0, 2.5), class = "Intercept"),
  brms::prior(exponential(2), class = "sd")
)

control_standard <- list(
  adapt_delta = 0.999,
  max_treedepth = 20
)

control_minor_issue <- list(
  adapt_delta = 0.9995,
  max_treedepth = 22
)

control_sparse <- list(
  adapt_delta = 0.9999,
  max_treedepth = 25
)

## *******************************************************
## Formula helpers
## *******************************************************

make_genus_brms_formula <- function(response_col,
                                    trials_col = "GenusScreened",
                                    terms = network_terms,
                                    family = "zero_inflated_binomial"){
  brms::bf(
    formula(
      paste0(
        response_col,
        " | trials(", trials_col, ") ~ ",
        paste(terms, collapse = " + ")
      )
    ),
    family = family
  )
}

make_genus_freq_formula <- function(response_col,
                                    trials_col = "GenusScreened",
                                    terms = network_terms){
  make_binomial_check_formula(
    response_col = response_col,
    trials_col = trials_col,
    terms = terms
  )
}

make_species_phylo_brms_formula <- function(response_col,
                                            trials_col = "SpScreened",
                                            terms = phylo_terms,
                                            family = "zero_inflated_binomial"){
  brms::bf(
    formula(
      paste0(
        response_col,
        " | trials(", trials_col, ") ~ ",
        paste(terms, collapse = " + ")
      )
    ),
    family = family
  )
}

make_species_freq_formula <- function(response_col,
                                      trials_col = "SpScreened",
                                      terms = freq_terms){
  make_binomial_check_formula(
    response_col = response_col,
    trials_col = trials_col,
    terms = terms
  )
}

make_apis_brms_formula <- function(response_col,
                                   trials_col = "SpScreened",
                                   terms = apis_terms,
                                   family = "binomial"){
  brms::bf(
    formula(
      paste0(
        response_col,
        " | trials(", trials_col, ") ~ ",
        paste(terms, collapse = " + ")
      )
    ),
    family = family
  )
}

## *******************************************************
## Model-running wrapper
## *******************************************************

run_model_from_spec <- function(spec,
                                ncores = get("ncores", envir = .GlobalEnv),
                                model_seed = get("model_seed", envir = .GlobalEnv),
                                iter = 10^4,
                                chains = 3,
                                thin = 1,
                                init = 0,
                                open_progress = FALSE){

  run_args <- list(
    model_id = spec$model_id,
    model_object_name = spec$model_object_name,
    data = spec$data,
    brms_formula = spec$brms_formula,
    frequentist_formula = spec$frequentist_formula,
    response_col = spec$response_col,
    trials_col = spec$trials_col,
    all_projects = spec$all_projects,
    ncores = ncores,
    iter = iter,
    chains = chains,
    thin = thin,
    init = init,
    seed = model_seed,
    control = spec$control,
    open_progress = open_progress
  )

  if (!is.null(spec$data2)) {
    run_args$data2 <- spec$data2
  }

  helper_args <- names(formals(run_host_parasite_model))
  helper_accepts_prior <- any(helper_args %in% c("prior", "..."))

  if (helper_accepts_prior && !is.null(spec$prior)) {
    run_args$prior <- spec$prior
  }

  do.call(run_host_parasite_model, run_args)
}

## *******************************************************
## Sampler diagnostic summary table
## *******************************************************

make_sampler_diagnostic_table <- function(model_specs,
                                          saved_model_objects = NULL,
                                          out_file = NULL,
                                          envir = .GlobalEnv){

  get_fit <- function(spec){
    if (!is.null(saved_model_objects)) {
      if (!is.null(saved_model_objects[[spec$model_object_name]])) {
        return(saved_model_objects[[spec$model_object_name]])
      }
      if (!is.null(saved_model_objects[[spec$model_id]])) {
        return(saved_model_objects[[spec$model_id]])
      }
    }

    if (exists(spec$model_object_name, envir = envir)) {
      return(get(spec$model_object_name, envir = envir))
    }

    if (exists(spec$model_id, envir = envir)) {
      return(get(spec$model_id, envir = envir))
    }

    NULL
  }

  diag_list <- lapply(model_specs, function(spec){
    fit <- get_fit(spec)

    if (is.null(fit)) {
      return(data.frame(
        model_id = spec$model_id,
        divergences = NA_integer_,
        treedepth_hits = NA_integer_,
        max_rhat = NA_real_,
        min_bulk_ess = NA_real_,
        min_tail_ess = NA_real_,
        diagnostic_flag = "model object not found"
      ))
    }

    np <- brms::nuts_params(fit)
    n_div <- sum(np$Value[np$Parameter == "divergent__"], na.rm = TRUE)

    max_td <- spec$control$max_treedepth
    n_td <- if (is.null(max_td)) {
      NA_integer_
    } else {
      sum(np$Value[np$Parameter == "treedepth__"] >= max_td, na.rm = TRUE)
    }

    draws_sum <- posterior::summarise_draws(posterior::as_draws_df(fit))
    max_rhat <- suppressWarnings(max(draws_sum$rhat, na.rm = TRUE))
    min_bulk <- suppressWarnings(min(draws_sum$ess_bulk, na.rm = TRUE))
    min_tail <- suppressWarnings(min(draws_sum$ess_tail, na.rm = TRUE))

    flag <- "looks okay"
    if (!is.na(n_div) && n_div > 0) {
      flag <- "divergences: refit or interpret cautiously"
    } else if (!is.na(n_td) && n_td > 0) {
      flag <- "max treedepth hits: consider higher max_treedepth or simplification"
    } else if (is.finite(max_rhat) && max_rhat > 1.01) {
      flag <- "Rhat > 1.01"
    } else if (is.finite(min_bulk) && min_bulk < 400) {
      flag <- "low bulk ESS"
    } else if (is.finite(min_tail) && min_tail < 400) {
      flag <- "low tail ESS"
    }

    data.frame(
      model_id = spec$model_id,
      divergences = n_div,
      treedepth_hits = n_td,
      max_rhat = max_rhat,
      min_bulk_ess = min_bulk,
      min_tail_ess = min_tail,
      diagnostic_flag = flag
    )
  })

  diag_table <- do.call(rbind, diag_list)

  if (!is.null(out_file)) {
    dir.create(dirname(out_file), showWarnings = FALSE, recursive = TRUE)
    write.csv(diag_table, file = out_file, row.names = FALSE)
  }

  diag_table
}
