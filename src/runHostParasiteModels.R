## Helper functions for fitting and checking host-parasite brms models.
##
## Source after:
##   src/writeResultsTable.R
##   src/makeBayesR2Table.R
##   src/makeModelOutputTables.R
##
## Main exported functions:
##   run_host_parasite_model()
##   check_saved_host_parasite_model()

model_input_dir <- file.path("data", "model_inputs")
pre_model_check_dir <- file.path("tables", "pre_model_checks")
post_model_check_dir <- file.path("tables", "post_model_checks")
native_model_check_dir <- file.path("figures", "model_checks")
pp_check_dir <- file.path("figures", "pp_checks")

make_dir <- function(path){
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
  invisible(path)
}

`%||%` <- function(x, y){
  if (is.null(x)) y else x
}

drop_nulls <- function(x){
  x[!vapply(x, is.null, logical(1))]
}

safe_one_line <- function(x){
  if (length(x) == 0 || all(is.na(x))) return("None")
  x <- sort(unique(as.character(x[!is.na(x)])))
  if (length(x) == 0) return("None")
  paste(x, collapse = "; ")
}

safe_write_lines <- function(x, file, append = FALSE){
  make_dir(dirname(file))
  cat(paste0(x, collapse = "\n"), "\n", file = file, append = append)
  invisible(file)
}

capture_output_warnings <- function(expr){
  warnings <- character()
  messages <- character()

  value <- tryCatch(
    withCallingHandlers(
      expr,
      warning = function(w){
        warnings <<- c(warnings, conditionMessage(w))
        invokeRestart("muffleWarning")
      },
      message = function(m){
        messages <<- c(messages, conditionMessage(m))
        invokeRestart("muffleMessage")
      }
    ),
    error = function(e){
      structure(list(error = conditionMessage(e)), class = "captured_error")
    }
  )

  list(value = value, warnings = warnings, messages = messages)
}

make_binomial_check_formula <- function(response_col, trials_col, terms){
  stats::as.formula(
    paste0(
      "cbind(", response_col, ", ", trials_col, " - ", response_col, ") ~ ",
      paste(terms, collapse = " + ")
    )
  )
}

make_project_summary <- function(data,
                                 response_col,
                                 trials_col,
                                 project_col = "ProjectSubProject",
                                 site_col = "Site",
                                 year_col = "Year"){
  if (!project_col %in% names(data)) {
    return(data.frame())
  }

  projects <- sort(unique(as.character(data[[project_col]][!is.na(data[[project_col]])])))

  out <- lapply(projects, function(pr){
    sub <- data[as.character(data[[project_col]]) == pr, , drop = FALSE]

    n_screened <- if (trials_col %in% names(sub)) sum(sub[[trials_col]], na.rm = TRUE) else NA_real_
    n_positive <- if (response_col %in% names(sub)) sum(sub[[response_col]], na.rm = TRUE) else NA_real_

    data.frame(
      ProjectSubProject = pr,
      NRows = nrow(sub),
      NSites = if (site_col %in% names(sub)) length(unique(sub[[site_col]][!is.na(sub[[site_col]])])) else NA_integer_,
      NYears = if (year_col %in% names(sub)) length(unique(sub[[year_col]][!is.na(sub[[year_col]])])) else NA_integer_,
      NScreened = n_screened,
      NPositive = n_positive,
      Prevalence = if (isTRUE(n_screened > 0)) n_positive / n_screened else NA_real_,
      stringsAsFactors = FALSE
    )
  })

  do.call(rbind, out)
}

validate_model_data <- function(data,
                                response_col,
                                trials_col,
                                required_cols = character()){
  missing_cols <- setdiff(unique(c(response_col, trials_col, required_cols)),
                          names(data))
  if (length(missing_cols) > 0) {
    stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
  }

  if (any(is.na(data[[response_col]]))) {
    warning("Response column contains NA values: ", response_col)
  }

  if (any(is.na(data[[trials_col]]))) {
    warning("Trials column contains NA values: ", trials_col)
  }

  bad_trials <- which(data[[trials_col]] < data[[response_col]])
  if (length(bad_trials) > 0) {
    stop(
      "Some rows have trials < positives for ",
      response_col,
      ". First bad row index: ",
      bad_trials[1]
    )
  }

  invisible(TRUE)
}

pre_model_check <- function(model_id,
                            data,
                            response_col,
                            trials_col,
                            frequentist_formula = NULL,
                            project_col = "ProjectSubProject",
                            site_col = "Site",
                            year_col = "Year",
                            all_projects = NULL,
                            run_glmmTMB = TRUE,
                            ziformula = ~1,
                            out_dir = pre_model_check_dir){
  make_dir(out_dir)

  log_file <- file.path(out_dir, paste0(model_id, "_precheck_log.txt"))
  project_file <- file.path(out_dir, paste0(model_id, "_project_summary.csv"))
  collinearity_file <- file.path(out_dir, paste0(model_id, "_collinearity.txt"))

  data <- droplevels(as.data.frame(data))

  validate_model_data(
    data = data,
    response_col = response_col,
    trials_col = trials_col,
    required_cols = c(project_col, site_col, year_col)
  )

  project_summary <- make_project_summary(
    data = data,
    response_col = response_col,
    trials_col = trials_col,
    project_col = project_col,
    site_col = site_col,
    year_col = year_col
  )

  utils::write.csv(project_summary, project_file, row.names = FALSE)

  included_projects <- if (project_col %in% names(data)) {
    sort(unique(as.character(data[[project_col]])))
  } else {
    character(0)
  }

  excluded_projects <- if (!is.null(all_projects)) {
    setdiff(sort(unique(as.character(all_projects))), included_projects)
  } else {
    character(0)
  }

  n_screened <- sum(data[[trials_col]], na.rm = TRUE)
  n_positive <- sum(data[[response_col]], na.rm = TRUE)

  log_lines <- c(
    paste("Model ID:", model_id),
    paste("Rows:", nrow(data)),
    paste("Response column:", response_col),
    paste("Trials column:", trials_col),
    paste("Total screened:", n_screened),
    paste("Total positives:", n_positive),
    paste("Overall prevalence:", if (isTRUE(n_screened > 0)) round(n_positive / n_screened, 4) else NA),
    paste("Included projects:", safe_one_line(included_projects)),
    paste("Excluded projects:", safe_one_line(excluded_projects)),
    paste("Project summary CSV:", project_file),
    "",
    "Missing values by key column:",
    paste(response_col, sum(is.na(data[[response_col]]))),
    paste(trials_col, sum(is.na(data[[trials_col]])))
  )

  safe_write_lines(log_lines, log_file, append = FALSE)

  if (isTRUE(run_glmmTMB) && !is.null(frequentist_formula)) {
    glmm_data <- data

    fit_result <- capture_output_warnings({
      glmmTMB::glmmTMB(
        formula = frequentist_formula,
        data = glmm_data,
        family = stats::binomial(),
        ziformula = ziformula
      )
    })

    safe_write_lines(c("", "glmmTMB pre-check:"), log_file, append = TRUE)

    if (inherits(fit_result$value, "captured_error")) {
      safe_write_lines(
        paste("glmmTMB failed:", fit_result$value$error),
        log_file,
        append = TRUE
      )
    } else {
      mod <- fit_result$value
      safe_write_lines(
        capture.output(print(summary(mod))),
        log_file,
        append = TRUE
      )

      coll_result <- capture_output_warnings({
        performance::check_collinearity(mod)
      })

      if (inherits(coll_result$value, "captured_error")) {
        safe_write_lines(
          paste("check_collinearity failed:", coll_result$value$error),
          log_file,
          append = TRUE
        )
      } else {
        safe_write_lines(capture.output(print(coll_result$value)),
                         collinearity_file,
                         append = FALSE)
        safe_write_lines(
          paste("Collinearity output:", collinearity_file),
          log_file,
          append = TRUE
        )
      }

      if (length(coll_result$warnings) > 0) {
        safe_write_lines(
          c("check_collinearity warnings:", unique(coll_result$warnings)),
          log_file,
          append = TRUE
        )
      }
    }

    if (length(fit_result$warnings) > 0) {
      safe_write_lines(
        c("glmmTMB warnings:", unique(fit_result$warnings)),
        log_file,
        append = TRUE
      )
    }

    if (length(fit_result$messages) > 0) {
      safe_write_lines(
        c("glmmTMB messages:", unique(fit_result$messages)),
        log_file,
        append = TRUE
      )
    }
  }

  invisible(list(
    log_file = log_file,
    project_file = project_file,
    project_summary = project_summary
  ))
}

save_model_input_csv <- function(data,
                                 model_id,
                                 out_dir = model_input_dir){
  make_dir(out_dir)
  out_file <- file.path(out_dir, paste0(model_id, "_data.csv"))
  utils::write.csv(data, out_file, row.names = FALSE)
  invisible(out_file)
}

run_host_parasite_model <- function(model_id,
                                    data,
                                    brms_formula,
                                    response_col,
                                    trials_col,
                                    frequentist_formula = NULL,
                                    host = NULL,
                                    parasite = NULL,
                                    model_set = NULL,
                                    all_projects = NULL,
                                    data2 = NULL,
                                    model_file = file.path("saved", paste0(model_id, ".Rdata")),
                                    model_object_name = NULL,
                                    ncores = 3,
                                    iter = 10000,
                                    chains = 3,
                                    thin = 1,
                                    init = 0,
                                    seed = 123,
                                    open_progress = FALSE,
                                    control = list(adapt_delta = 0.999,
                                                   stepsize = 0.001,
                                                   max_treedepth = 20),
                                    backend = NULL,
                                    run_precheck = TRUE,
                                    run_glmmTMB_precheck = TRUE,
                                    write_ms_table = TRUE,
                                    ...){
  make_dir(dirname(model_file))

  model_data <- droplevels(as.data.frame(data))

  if (isTRUE(run_precheck)) {
    pre_model_check(
      model_id = model_id,
      data = model_data,
      response_col = response_col,
      trials_col = trials_col,
      frequentist_formula = frequentist_formula,
      all_projects = all_projects,
      run_glmmTMB = run_glmmTMB_precheck
    )
  }

  input_file <- save_model_input_csv(model_data, model_id)

  meta <- if (exists("parse_bayes_R2_model_id", mode = "function")) {
    parse_bayes_R2_model_id(model_id)
  } else {
    list(ModelSet = model_set %||% NA_character_,
         Host = host %||% NA_character_,
         Parasite = parasite %||% NA_character_)
  }

  model_metadata <- list(
    model_id = model_id,
    model_set = model_set %||% meta$ModelSet,
    host = host %||% meta$Host,
    parasite = parasite %||% meta$Parasite,
    response_col = response_col,
    trials_col = trials_col,
    all_projects = all_projects,
    model_input_csv = input_file,
    model_file = model_file,
    seed = seed,
    iter = iter,
    chains = chains,
    thin = thin,
    control = control
  )

  brm_args <- list(
    formula = brms_formula,
    data = model_data,
    data2 = data2,
    cores = ncores,
    iter = iter,
    chains = chains,
    thin = thin,
    init = init,
    open_progress = open_progress,
    seed = seed,
    control = control,
    backend = backend,
    ...
  )
  brm_args <- drop_nulls(brm_args)

  message("Fitting ", model_id)
  model <- do.call(brms::brm, brm_args)

  if (isTRUE(write_ms_table) && exists("write.ms.table", mode = "function")) {
    tryCatch(
      write.ms.table(model, model_id),
      error = function(e) warning("write.ms.table failed for ", model_id, ": ",
                                  conditionMessage(e))
    )
  }

  save_env <- new.env(parent = emptyenv())
  save_env$model <- model
  save_env$model_data <- model_data
  save_env$model_metadata <- model_metadata

  if (!is.null(model_object_name)) {
    assign(model_object_name, model, envir = save_env)
  }

  save(list = ls(save_env), file = model_file, envir = save_env)

  invisible(list(
    model = model,
    model_file = model_file,
    model_metadata = model_metadata
  ))
}

load_saved_brms_model <- function(model_file){
  load_env <- new.env(parent = emptyenv())
  load(model_file, envir = load_env)

  if (exists("model", envir = load_env)) {
    model <- get("model", envir = load_env)
  } else {
    obj_names <- ls(load_env)
    is_brmsfit <- vapply(obj_names, function(nm){
      inherits(get(nm, envir = load_env), "brmsfit")
    }, logical(1))
    if (!any(is_brmsfit)) {
      stop("No brmsfit object found in ", model_file)
    }
    model <- get(obj_names[which(is_brmsfit)[1]], envir = load_env)
  }

  model_data <- if (exists("model_data", envir = load_env)) {
    get("model_data", envir = load_env)
  } else {
    NULL
  }

  model_metadata <- if (exists("model_metadata", envir = load_env)) {
    get("model_metadata", envir = load_env)
  } else {
    list()
  }

  list(model = model, model_data = model_data, model_metadata = model_metadata)
}

save_pp_check_pdf <- function(model,
                              model_id,
                              resp = NULL,
                              ndraws = 1000,
                              width = 7,
                              height = 5,
                              out_dir = pp_check_dir){
  make_dir(out_dir)
  out_file <- file.path(out_dir, paste0(model_id, ".pdf"))

  result <- capture_output_warnings({
    if (is.null(resp)) {
      brms::pp_check(model, ndraws = ndraws)
    } else {
      brms::pp_check(model, resp = resp, ndraws = ndraws)
    }
  })

  log_file <- file.path(out_dir, paste0(model_id, "_pp_check_log.txt"))

  if (inherits(result$value, "captured_error")) {
    safe_write_lines(
      paste("pp_check failed for", model_id, ":", result$value$error),
      log_file,
      append = FALSE
    )
    warning("pp_check failed for ", model_id, ": ", result$value$error)
    return(invisible(out_file))
  }

  grDevices::pdf(out_file, width = width, height = height)
  print(result$value)
  grDevices::dev.off()

  if (length(result$warnings) > 0 || length(result$messages) > 0) {
    safe_write_lines(
      c("pp_check messages/warnings:",
        unique(result$warnings),
        unique(result$messages)),
      log_file,
      append = FALSE
    )
  }

  invisible(out_file)
}

save_model_summary_txt <- function(model,
                                   model_id,
                                   out_dir = file.path("tables", "model_summaries")){
  make_dir(out_dir)
  out_file <- file.path(out_dir, paste0(model_id, "_summary.txt"))
  capture.output(summary(model), file = out_file)
  invisible(out_file)
}

save_draw_summary_csv <- function(model,
                                  model_id,
                                  out_dir = post_model_check_dir){
  make_dir(out_dir)
  out_file <- file.path(out_dir, paste0(model_id, "_draw_summary.csv"))

  result <- capture_output_warnings({
    draws <- brms::as_draws_df(model)
    posterior::summarise_draws(draws)
  })

  log_file <- file.path(out_dir, paste0(model_id, "_draw_summary_log.txt"))

  if (inherits(result$value, "captured_error")) {
    safe_write_lines(
      paste("posterior::summarise_draws failed:", result$value$error),
      log_file,
      append = FALSE
    )
    return(invisible(NULL))
  }

  utils::write.csv(as.data.frame(result$value), out_file, row.names = FALSE)

  if (length(result$warnings) > 0 || length(result$messages) > 0) {
    safe_write_lines(
      c("posterior::summarise_draws messages/warnings:",
        unique(result$warnings),
        unique(result$messages)),
      log_file,
      append = FALSE
    )
  }

  invisible(out_file)
}

save_nuts_params_csv <- function(model,
                                 model_id,
                                 out_dir = post_model_check_dir){
  make_dir(out_dir)
  out_file <- file.path(out_dir, paste0(model_id, "_nuts_params.csv"))
  summary_file <- file.path(out_dir, paste0(model_id, "_nuts_summary.txt"))

  result <- capture_output_warnings({
    brms::nuts_params(model)
  })

  if (inherits(result$value, "captured_error") || is.null(result$value)) {
    safe_write_lines(
      paste("nuts_params unavailable for", model_id),
      summary_file,
      append = FALSE
    )
    return(invisible(NULL))
  }

  np <- result$value
  utils::write.csv(np, out_file, row.names = FALSE)

  divergent <- if (all(c("Parameter", "Value") %in% names(np))) {
    sum(np$Parameter == "divergent__" & np$Value == 1, na.rm = TRUE)
  } else {
    NA_integer_
  }

  max_treedepth <- if (all(c("Parameter", "Value") %in% names(np))) {
    suppressWarnings(max(np$Value[np$Parameter == "treedepth__"], na.rm = TRUE))
  } else {
    NA_real_
  }

  safe_write_lines(
    c(
      paste("Model ID:", model_id),
      paste("Divergent transitions:", divergent),
      paste("Maximum treedepth reached:", max_treedepth),
      paste("NUTS parameters CSV:", out_file)
    ),
    summary_file,
    append = FALSE
  )

  invisible(out_file)
}

save_brms_mcmc_plot_pdf <- function(model,
                                    model_id,
                                    type,
                                    out_dir = native_model_check_dir,
                                    width = 11,
                                    height = 8.5){
  make_dir(out_dir)
  out_file <- file.path(out_dir, paste0(model_id, "_", type, ".pdf"))
  log_file <- file.path(out_dir, paste0(model_id, "_", type, "_log.txt"))

  result <- capture_output_warnings({
    brms::mcmc_plot(model, type = type)
  })

  if (inherits(result$value, "captured_error")) {
    safe_write_lines(
      paste("mcmc_plot type", type, "failed for", model_id, ":", result$value$error),
      log_file,
      append = FALSE
    )
    return(invisible(NULL))
  }

  grDevices::pdf(out_file, width = width, height = height)
  print(result$value)
  grDevices::dev.off()

  if (length(result$warnings) > 0 || length(result$messages) > 0) {
    safe_write_lines(
      c(paste("mcmc_plot", type, "messages/warnings:"),
        unique(result$warnings),
        unique(result$messages)),
      log_file,
      append = FALSE
    )
  }

  invisible(out_file)
}

save_brms_native_diagnostics <- function(model,
                                         model_id,
                                         mcmc_types = c("trace", "rhat", "neff"),
                                         out_dir = native_model_check_dir){
  make_dir(out_dir)

  files <- lapply(mcmc_types, function(tp){
    save_brms_mcmc_plot_pdf(model, model_id, type = tp, out_dir = out_dir)
  })

  save_draw_summary_csv(model, model_id)
  save_nuts_params_csv(model, model_id)

  invisible(files)
}

save_optional_check_model_pdf <- function(model,
                                          model_id,
                                          width = 11,
                                          height = 8.5,
                                          out_dir = native_model_check_dir){
  make_dir(out_dir)
  out_file <- file.path(out_dir, paste0(model_id, "_performance_check_model.pdf"))
  log_file <- file.path(out_dir, paste0(model_id, "_performance_check_model_log.txt"))

  result <- capture_output_warnings({
    performance::check_model(
      model,
      residual_type = "normal",
      panel = TRUE,
      verbose = FALSE
    )
  })

  if (inherits(result$value, "captured_error")) {
    safe_write_lines(
      paste("performance::check_model failed for", model_id, ":", result$value$error),
      log_file,
      append = FALSE
    )
    return(invisible(NULL))
  }

  grDevices::pdf(out_file, width = width, height = height)
  print(result$value)
  grDevices::dev.off()

  safe_write_lines(
    c("performance::check_model completed.",
      paste("Output:", out_file),
      "Messages/warnings:",
      unique(result$warnings),
      unique(result$messages)),
    log_file,
    append = FALSE
  )

  invisible(out_file)
}

check_saved_host_parasite_model <- function(model_id,
                                            model_file = file.path("saved", paste0(model_id, ".Rdata")),
                                            resp = NULL,
                                            ndraws = 1000,
                                            bayes_R2_table = NULL,
                                            project_inclusion_table = NULL,
                                            run_check_model = FALSE,
                                            mcmc_types = c("trace", "rhat", "neff")){
  saved <- load_saved_brms_model(model_file)
  model <- saved$model
  model_data <- saved$model_data
  meta <- saved$model_metadata

  resp <- resp %||% meta$response_col
  trials_col <- meta$trials_col %||% NA_character_
  all_projects <- meta$all_projects %||% NULL

  save_model_summary_txt(model, model_id)
  save_brms_native_diagnostics(model, model_id, mcmc_types = mcmc_types)
  save_pp_check_pdf(model, model_id, resp = resp, ndraws = ndraws)

  if (isTRUE(run_check_model)) {
    save_optional_check_model_pdf(model, model_id)
  }

  if (is.null(bayes_R2_table)) {
    bayes_R2_table <- init_bayes_R2_table()
  }

  bayes_R2_table <- record_bayes_R2(
    bayes_R2_table,
    model,
    model_id = model_id
  )

  if (is.null(project_inclusion_table)) {
    project_inclusion_table <- init_model_project_table()
  }

  if (!is.null(model_data) && !is.na(trials_col)) {
    project_inclusion_table <- record_model_projects(
      project_inclusion_table,
      model_id = model_id,
      data = model_data,
      response_col = resp,
      trials_col = trials_col,
      all_projects = all_projects
    )
  }

  invisible(list(
    bayes_R2_table = bayes_R2_table,
    project_inclusion_table = project_inclusion_table
  ))
}
