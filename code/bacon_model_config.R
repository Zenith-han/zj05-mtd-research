bacon_model_config <- list(
  core = "ZJ05_formal_chronology",
  thick = 20,
  d.min = 247,
  d.max = 7880,
  d.by = 1,
  acc.mean = 28,
  acc.shape = 1,
  mem.mean = 0.95,
  mem.strength = 20,
  ssize = 18000L,
  burnin = 500L,
  seed = 20264621L,
  normal = FALSE,
  archived_relative_path = file.path("data", "processed", "bacon_model.rds"),
  rebuilt_relative_path = file.path(
    "bacon_model_build", "ZJ05_formal_chronology", "bacon_model.rds"
  ),
  output_columns = 385L,
  range_rows = 7634L
)

bacon_result_relative_paths <- c(
  file.path("bacon", "age_depth_model_1cm.csv"),
  file.path("bacon", "event_age_posterior_draws_wide.csv"),
  file.path("bacon", "mtd_linked_sedimentation_rate_posterior_summary.csv"),
  file.path("bacon", "stage_sedimentation_rate_posterior_summary.csv"),
  file.path("bacon", "ZJ05_LSR_step_10ka_pseudoPosterior.csv"),
  "bivariate_enrichment_summary.csv",
  "categorized_timeseries.csv",
  "circular_shift_by_offset_midRSL_bivariate.csv",
  "circular_shift_by_offset_univariate.csv",
  "circular_shift_robustness_midRSL_bivariate.csv",
  "circular_shift_robustness_univariate.csv",
  "event_age_summary.csv",
  "midRSL_bivariate_draws_ALL.csv",
  "posterior_mean_phase_shift_null_long.csv",
  "posterior_mean_phase_shift_null_summary.csv",
  "univariate_draws_ALL.csv",
  "univariate_enrichment_summary.csv"
)

bacon_model_file <- function(
  public_root,
  model_source,
  must_exist = TRUE,
  results_root = file.path(public_root, "results")
) {
  public_root <- normalizePath(public_root, mustWork = TRUE)
  if (!model_source %in% c("archived", "rebuilt")) {
    stop("--model-source must be archived or rebuilt.")
  }
  model_file <- if (model_source == "archived") {
    file.path(public_root, bacon_model_config$archived_relative_path)
  } else {
    file.path(
      normalizePath(results_root, mustWork = FALSE),
      bacon_model_config$rebuilt_relative_path
    )
  }
  if (must_exist && !file.exists(model_file)) {
    label <- if (model_source == "archived") "Archived" else "Rebuilt"
    stop(label, " Bacon model object not found: ", model_file)
  }
  normalizePath(model_file, mustWork = must_exist)
}

validate_bacon_run_id <- function(run_id) {
  if (!is.character(run_id) || length(run_id) != 1L || is.na(run_id) ||
      !grepl("^[A-Za-z0-9][A-Za-z0-9._-]{0,63}$", run_id)) {
    stop(
      "Invalid --run-id. Use 1-64 ASCII letters, digits, dots, underscores, ",
      "or hyphens, starting with a letter or digit."
    )
  }
  run_id
}

timestamp_bacon_run_id <- function(now = Sys.time()) {
  time_zone <- attr(now, "tzone")
  if (is.null(time_zone) || !length(time_zone)) time_zone <- ""
  format(now, "%Y%m%d-%H%M%S", tz = time_zone)
}

assert_isolated_bacon_build_root <- function(public_root, selected_results_root) {
  public_root <- normalizePath(public_root, winslash = "/", mustWork = TRUE)
  selected_results_root <- normalizePath(
    selected_results_root,
    winslash = "/",
    mustWork = FALSE
  )
  formal_results_root <- normalizePath(
    file.path(public_root, "results"),
    winslash = "/",
    mustWork = FALSE
  )
  rebuilt_parent <- normalizePath(
    file.path(public_root, "results", "rebuilt_runs"),
    winslash = "/",
    mustWork = FALSE
  )
  parent_is_isolated <- identical(
    tolower(dirname(rebuilt_parent)),
    tolower(formal_results_root)
  )
  is_direct_child <- identical(
    tolower(dirname(selected_results_root)),
    tolower(rebuilt_parent)
  )
  run_id <- basename(selected_results_root)
  valid_run_id <- grepl(
    "^[A-Za-z0-9][A-Za-z0-9._-]{0,63}$",
    run_id
  )
  if (!parent_is_isolated || !is_direct_child || !valid_run_id) {
    stop(
      "Bacon model builds require an isolated rebuilt run result root under ",
      "results/rebuilt_runs/<run-id>; refusing to write to: ",
      selected_results_root
    )
  }
  invisible(TRUE)
}

reserve_bacon_run_directory <- function(run_root) {
  run_root <- normalizePath(run_root, winslash = "/", mustWork = FALSE)
  parent <- dirname(run_root)
  if (!dir.exists(parent) &&
      !dir.create(parent, recursive = TRUE, showWarnings = FALSE) &&
      !dir.exists(parent)) {
    stop("Could not create rebuilt-runs directory: ", parent)
  }
  if (file.exists(run_root)) {
    stop("Rebuilt run directory already exists: ", run_root)
  }
  if (!dir.create(run_root, recursive = FALSE, showWarnings = FALSE)) {
    if (file.exists(run_root)) {
      stop("Rebuilt run directory already exists: ", run_root)
    }
    stop("Could not reserve rebuilt run directory: ", run_root)
  }
  invisible(normalizePath(run_root, winslash = "/", mustWork = TRUE))
}

reserve_isolated_bacon_run_directory <- function(public_root, run_root) {
  assert_isolated_bacon_build_root(public_root, run_root)
  reserved_root <- reserve_bacon_run_directory(run_root)
  assert_isolated_bacon_build_root(public_root, reserved_root)
  invisible(reserved_root)
}

acquire_bacon_run_lock <- function(run_root) {
  run_root <- normalizePath(run_root, winslash = "/", mustWork = TRUE)
  lock_dir <- file.path(run_root, ".downstream-results.lock")
  if (file.exists(lock_dir) ||
      !dir.create(lock_dir, recursive = FALSE, showWarnings = FALSE)) {
    if (file.exists(lock_dir)) {
      stop("Rebuilt run is already being updated: ", run_root)
    }
    stop("Could not acquire rebuilt-run lock: ", lock_dir)
  }
  invisible(normalizePath(lock_dir, winslash = "/", mustWork = TRUE))
}

release_bacon_run_lock <- function(lock_dir) {
  lock_dir <- normalizePath(lock_dir, winslash = "/", mustWork = FALSE)
  if (!identical(basename(lock_dir), ".downstream-results.lock")) {
    stop("Refusing to release an invalid rebuilt-run lock path: ", lock_dir)
  }
  if (file.exists(lock_dir) &&
      unlink(lock_dir, recursive = TRUE, force = TRUE) != 0L) {
    stop("Could not release rebuilt-run lock: ", lock_dir)
  }
  invisible(TRUE)
}

parse_bacon_run_mode <- function(
  args = commandArgs(trailingOnly = TRUE),
  public_root,
  now = Sys.time()
) {
  public_root <- normalizePath(public_root, mustWork = TRUE)
  recognized <- args == "--build-bacon-model" |
    startsWith(args, "--model-source=") |
    startsWith(args, "--run-id=")
  unknown <- args[!recognized]
  if (length(unknown) > 0L) {
    stop(
      "Unknown command-line argument: ",
      paste(unique(unknown), collapse = ", ")
    )
  }
  build_flags <- args[args == "--build-bacon-model"]
  if (length(build_flags) > 1L) {
    stop("--build-bacon-model may be supplied only once.")
  }
  build_model <- length(build_flags) == 1L

  source_args <- args[startsWith(args, "--model-source=")]
  if (length(source_args) > 1L) {
    stop("--model-source may be supplied only once.")
  }
  source_is_explicit <- length(source_args) == 1L

  run_id_args <- args[startsWith(args, "--run-id=")]
  if (length(run_id_args) > 1L) {
    stop("--run-id may be supplied only once.")
  }
  supplied_run_id <- if (length(run_id_args) == 1L) {
    sub("--run-id=", "", run_id_args, fixed = TRUE)
  } else {
    NULL
  }

  model_source <- if (source_is_explicit) {
    sub("--model-source=", "", source_args, fixed = TRUE)
  } else if (build_model) {
    "rebuilt"
  } else {
    "archived"
  }
  if (!model_source %in% c("archived", "rebuilt")) {
    stop("--model-source must be archived or rebuilt.")
  }
  if (build_model && source_is_explicit && model_source == "archived") {
    stop("--build-bacon-model cannot be combined with --model-source=archived.")
  }

  if (model_source == "archived" && !is.null(supplied_run_id)) {
    stop("--run-id is only valid for rebuilt Bacon runs.")
  }

  run_id <- NULL
  new_run <- FALSE
  selected_results_root <- file.path(public_root, "results")
  if (model_source == "rebuilt") {
    if (build_model) {
      run_id <- if (is.null(supplied_run_id)) {
        timestamp_bacon_run_id(now)
      } else {
        validate_bacon_run_id(supplied_run_id)
      }
      new_run <- TRUE
    } else {
      if (is.null(supplied_run_id)) {
        stop("--run-id is required when reusing an existing rebuilt Bacon model.")
      }
      run_id <- validate_bacon_run_id(supplied_run_id)
    }
    selected_results_root <- file.path(
      public_root, "results", "rebuilt_runs", run_id
    )
    if (new_run && file.exists(selected_results_root)) {
      stop("Rebuilt run directory already exists: ", selected_results_root)
    }
  }

  run_downstream <- !(build_model && !source_is_explicit)
  require_existing <- !build_model
  list(
    build_model = build_model,
    run_downstream = run_downstream,
    model_source = model_source,
    run_id = run_id,
    results_root = normalizePath(selected_results_root, mustWork = FALSE),
    new_run = new_run,
    model_file = bacon_model_file(
      public_root,
      model_source,
      must_exist = require_existing,
      results_root = selected_results_root
    )
  )
}

normalize_bacon_model_output <- function(model) {
  if (is.null(model$output) || is.null(dim(model$output))) {
    stop("The selected Bacon model object has no output matrix.")
  }
  output <- as.matrix(model$output)
  expected_columns <- bacon_model_config$output_columns
  if (ncol(output) == expected_columns) {
    return(model)
  }
  if (ncol(output) != expected_columns + 1L) {
    stop(
      "Unexpected Bacon model output width: ", ncol(output),
      "; expected ", expected_columns, " model columns."
    )
  }
  diagnostic <- output[, ncol(output)]
  if (is.null(model$Us) || length(model$Us) != nrow(output) ||
      !isTRUE(all.equal(
        as.numeric(diagnostic),
        as.numeric(model$Us),
        tolerance = 0,
        check.attributes = FALSE
      ))) {
    stop("The trailing Bacon diagnostic column does not match model$Us.")
  }
  model$output <- output[, -ncol(output), drop = FALSE]
  model
}

assert_bacon_model_configuration <- function(model) {
  required <- c(
    "core", "thick", "d.min", "d.max", "d.by", "acc.mean", "acc.shape",
    "mem.mean", "mem.strength", "ssize", "seed", "normal", "output", "ranges"
  )
  missing <- required[!vapply(required, function(name) !is.null(model[[name]]), logical(1))]
  if (length(missing) > 0L) {
    stop(
      "The selected Bacon model object is missing fields: ",
      paste(missing, collapse = ", ")
    )
  }

  numeric_fields <- c(
    "thick", "d.min", "d.max", "d.by", "acc.mean", "acc.shape",
    "mem.mean", "mem.strength", "ssize", "seed"
  )
  numeric_mismatch <- numeric_fields[!vapply(
    numeric_fields,
    function(name) isTRUE(all.equal(
      as.numeric(model[[name]]),
      as.numeric(bacon_model_config[[name]]),
      tolerance = 0
    )),
    logical(1)
  )]
  if (!identical(as.character(model$core), bacon_model_config$core) ||
      length(model$normal) != 1L ||
      !identical(as.logical(model$normal), bacon_model_config$normal) ||
      length(numeric_mismatch) > 0L) {
    stop("The selected Bacon model object does not match the final parameter configuration.")
  }

  if (!identical(
    as.integer(dim(model$output)),
    c(bacon_model_config$ssize, bacon_model_config$output_columns)
  )) {
    stop("The selected Bacon model output matrix must be 18,000 x 385.")
  }
  if (!identical(
    as.integer(dim(model$ranges)),
    c(bacon_model_config$range_rows, 5L)
  )) {
    stop("The selected Bacon model 1-cm summary must be 7,634 x 5.")
  }
  invisible(TRUE)
}
