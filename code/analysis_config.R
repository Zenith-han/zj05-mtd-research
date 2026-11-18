if (!exists("bacon_model_config", inherits = TRUE)) {
  if (!exists("script_dir", mode = "function")) {
    stop("script_dir() is required before sourcing analysis_config.R.")
  }
  source(file.path(script_dir(), "bacon_model_config.R"))
}

analysis_config <- list(
  analysis_start_ka = 0,
  analysis_end_ka = 237,
  step_ka = 0.2,
  variability_half_window_ka = 5,
  variability_half_window_points = 25L,
  variability_window_points = 51L,
  n_events = 11L,
  n_enrichment_iter = 5000L,
  enrichment_seed = 1L,
  shift_maxima_ka = c(3, 5, 10),
  bacon_model_id = bacon_model_config$core,
  posterior_iterations = bacon_model_config$ssize
)

stage_sedimentation_intervals <- data.frame(
  stage_id = c("early_middle_MIS_3", "MIS_5e", "MIS_2", "MIS_6"),
  stage_label = c("early-middle MIS 3", "MIS 5e", "MIS 2", "MIS 6"),
  age_left_ka = c(35, 116, 14, 130),
  age_right_ka = c(57, 130, 29, 191),
  stringsAsFactors = FALSE
)

sedimentation_rate_10ka_edges_ka <- seq(0, 230, by = 10)

mtd_sedimentation_rate_groups <- list(
  `MTD1-MTD3` = paste0("MTD", 1:3),
  `MTD4-MTD7` = paste0("MTD", 4:7),
  `MTD8-MTD11` = paste0("MTD", 8:11)
)

analysis_duration_ka <- with(
  analysis_config,
  analysis_end_ka - analysis_start_ka
)

analysis_grid_ka <- with(
  analysis_config,
  seq(analysis_start_ka, analysis_end_ka, by = step_ka)
)

if (length(analysis_grid_ka) != 1186L ||
    tail(analysis_grid_ka, 1) != analysis_config$analysis_end_ka) {
  stop("The configured 0-237 ka grid must contain exactly 1,186 points.")
}

if (analysis_config$variability_half_window_points * analysis_config$step_ka !=
      analysis_config$variability_half_window_ka ||
    analysis_config$variability_window_points !=
      2L * analysis_config$variability_half_window_points + 1L) {
  stop("The variability window must be a symmetric 51-point window spanning +/-5 ka.")
}
