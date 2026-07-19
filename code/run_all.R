script_dir <- function() {
  args <- commandArgs(FALSE)
  hit <- grep("^--file=", args, value = TRUE)
  if (length(hit) > 0) return(dirname(normalizePath(sub("^--file=", "", hit[1]))))
  if (basename(getwd()) == "code") return(normalizePath(getwd()))
  normalizePath(file.path(getwd(), "code"), mustWork = TRUE)
}
code_dir <- script_dir()
source(file.path(code_dir, "00_paths.R"))
source(file.path(code_dir, "bacon_model_config.R"))

run_all_pipeline <- function(args = commandArgs(trailingOnly = TRUE)) {
  mode <- parse_bacon_run_mode(args, root_dir)
  if (mode$new_run) {
    reserve_isolated_bacon_run_directory(root_dir, mode$results_root)
  } else if (mode$model_source == "rebuilt") {
    assert_isolated_bacon_build_root(root_dir, mode$results_root)
  }

  run_lock <- NULL
  if (mode$model_source == "rebuilt" && mode$run_downstream) {
    assert_isolated_bacon_build_root(root_dir, mode$results_root)
    run_lock <- acquire_bacon_run_lock(mode$results_root)
    on.exit({
      if (!is.null(run_lock)) release_bacon_run_lock(run_lock)
    }, add = TRUE)
  }

  previous_model_file <- getOption("zj05.bacon_model_file", NULL)
  previous_results_root <- getOption("zj05.results_root", NULL)
  options(
    zj05.bacon_model_file = mode$model_file,
    zj05.results_root = mode$results_root
  )
  on.exit(options(
    zj05.bacon_model_file = previous_model_file,
    zj05.results_root = previous_results_root
  ), add = TRUE)

  message("Bacon model source: ", mode$model_source)
  if (!is.null(mode$run_id)) message("Rebuilt run ID: ", mode$run_id)
  message("Bacon model file: ", mode$model_file)
  message("Result root: ", mode$results_root)

  if (mode$build_model) {
    message("\n=== Running 01_build_bacon_chronology.R ===")
    source(file.path(code_dir, "01_build_bacon_chronology.R"), chdir = FALSE)
    mode$model_file <- bacon_model_file(
      root_dir,
      "rebuilt",
      must_exist = TRUE,
      results_root = mode$results_root
    )
    options(zj05.bacon_model_file = mode$model_file)
    if (!mode$run_downstream) {
      message("Bacon model build completed; downstream analysis was not requested.")
      return(invisible(mode))
    }
  }

  staging_root <- NULL
  if (mode$model_source == "rebuilt") {
    staging_root <- create_result_staging_root(mode$results_root)
    options(zj05.results_root = staging_root)
    on.exit({
      if (!is.null(staging_root) && dir.exists(staging_root)) {
        unlink(staging_root, recursive = TRUE, force = TRUE)
      }
    }, add = TRUE)
    message("Staging rebuilt results in: ", staging_root)
  }

  for (script in c(
    "02_extract_mtd_age_posteriors.R",
    "03_summarize_chronology.R",
    "04_build_environment_timeseries.R",
    "05_calculate_enrichment_lambda.R",
    "06_posterior_mean_phase_shift_null.R",
    "07_circular_shift_robustness.R"
  )) {
    message("\n=== Running ", script, " ===")
    source(file.path(code_dir, script), chdir = FALSE)
  }
  if (mode$model_source == "rebuilt") {
    commit_staged_result_files(
      staging_root,
      mode$results_root,
      bacon_result_relative_paths
    )
    staging_root <- NULL
    options(zj05.results_root = mode$results_root)
    message("Committed rebuilt results to: ", mode$results_root)
  }
  invisible(mode)
}

run_all_pipeline()
