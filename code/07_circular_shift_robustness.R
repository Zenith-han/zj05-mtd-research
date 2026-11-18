script_dir <- function() {
  args <- commandArgs(FALSE)
  hit <- grep("^--file=", args, value = TRUE)
  if (length(hit) > 0) return(dirname(normalizePath(sub("^--file=", "", hit[1]))))
  if (basename(getwd()) == "code") return(normalizePath(getwd()))
  normalizePath(file.path(getwd(), "code"), mustWork = TRUE)
}
source(file.path(script_dir(), "00_paths.R"))
source(file.path(script_dir(), "analysis_config.R"))
source(file.path(script_dir(), "environment_classification.R"))

vars <- c(
  "RSL_State", "RSL_Variability", "RSL_RateofChange",
  "EASM_Intensity", "EASM_Variability", "EASM_RateofChange"
)
env <- read_result_csv_required("categorized_timeseries.csv")
env <- env[order(env$age_ka), ]
cell_duration_ka <- grid_exposure_weights(env$age_ka, 0, 237)
total_time_ka <- sum(cell_duration_ka)
posterior <- read_result_csv_required(file.path(
  "bacon", "event_age_posterior_draws_wide.csv"
))
events <- read_csv_required(path("data", "mtd_events.csv"))
event_ids <- events$event_id
main_uni <- read_result_csv_required("univariate_enrichment_summary.csv")
main_bi <- read_result_csv_required("bivariate_enrichment_summary.csv")
main_mid <- main_bi[main_bi$varA == "RSL_State" & main_bi$A_level == "Mid", ]
if (nrow(main_mid) != 21L) {
  stop("Expected 21 Mid-RSL bivariate rows in bivariate_enrichment_summary.csv.")
}

set.seed(analysis_config$enrichment_seed)
sample_rows <- sample(
  seq_len(nrow(posterior)),
  analysis_config$n_enrichment_iter,
  replace = FALSE
)
selected_iterations <- posterior$iteration[sample_rows]
uni_draws <- read_result_csv_required("univariate_draws_ALL.csv")
check_draws <- uni_draws[
  uni_draws$variable == "RSL_State" & uni_draws$level == "High",
  c("iteration", "bacon_iteration")
]
check_draws <- check_draws[order(check_draws$iteration), ]
if (!identical(as.integer(check_draws$bacon_iteration), as.integer(selected_iterations))) {
  stop("The fixed posterior rows do not match the canonical enrichment draws.")
}

selected_ages <- as.matrix(posterior[sample_rows, event_ids, drop = FALSE])
storage.mode(selected_ages) <- "double"
nearest_index <- function(age_ka) {
  idx <- round((age_ka - env$age_ka[1]) / analysis_config$step_ka) + 1L
  pmin(pmax(idx, 1L), nrow(env))
}
event_idx <- nearest_index(selected_ages)
shift_series <- function(x, steps) {
  x[((seq_along(x) - 1L - steps) %% length(x)) + 1L]
}
safe_r <- function(count, tc) {
  (count / analysis_config$n_events) / (tc / total_time_ka)
}
summarise_draws <- function(values) {
  data.frame(
    R_mean = mean(values),
    R_p05 = stats::quantile(values, 0.05, names = FALSE),
    R_p50 = stats::quantile(values, 0.50, names = FALSE),
    R_p95 = stats::quantile(values, 0.95, names = FALSE)
  )
}
summarise_range <- function(values, obs, max_ka, n_offsets) {
  shift_means <- rowMeans(values)
  data.frame(
    shiftmax_ka = max_ka,
    n_offsets = n_offsets,
    mean_R = mean(values),
    R_p05 = stats::quantile(values, 0.05, names = FALSE),
    R_p50 = stats::quantile(values, 0.50, names = FALSE),
    R_p95 = stats::quantile(values, 0.95, names = FALSE),
    shift_mean_R_min = min(shift_means),
    shift_mean_R_max = max(shift_means),
    shift_retention_rate = if (obs >= 1) mean(shift_means >= 1) else mean(shift_means < 1)
  )
}

uni_specs <- data.frame(
  canonical_key = paste("univariate", main_uni$var, main_uni$level, sep = "|"),
  variable = main_uni$var,
  level = main_uni$level,
  Tc_ka = main_uni$Tc_ka,
  R_enrichment = main_uni$R_mean,
  stringsAsFactors = FALSE
)
bi_specs <- data.frame(
  canonical_key = paste(
    "bivariate", main_mid$varA, main_mid$A_level,
    main_mid$varB, main_mid$B_level, sep = "|"
  ),
  RSL_level = "Mid",
  y_var = main_mid$varB,
  y_level = main_mid$B_level,
  Tc_ka = main_mid$Tc_ka,
  R_enrichment = main_mid$R_mean,
  stringsAsFactors = FALSE
)

max_steps <- round(max(analysis_config$shift_maxima_ka) / analysis_config$step_ka)
steps <- seq.int(-max_steps, max_steps)
shift_ka <- steps * analysis_config$step_ka
uni_values <- lapply(seq_len(nrow(uni_specs)), function(i) {
  matrix(NA_real_, nrow = length(steps), ncol = analysis_config$n_enrichment_iter)
})
bi_values <- lapply(seq_len(nrow(bi_specs)), function(i) {
  matrix(NA_real_, nrow = length(steps), ncol = analysis_config$n_enrichment_iter)
})

for (sidx in seq_along(steps)) {
  shifted <- lapply(vars, function(v) shift_series(env[[v]], steps[sidx]))
  names(shifted) <- vars
  shifted_cell_duration_ka <- shift_series(cell_duration_ka, steps[sidx])

  for (j in seq_len(nrow(uni_specs))) {
    in_cat <- shifted[[uni_specs$variable[j]]] == uni_specs$level[j]
    tc <- sum(shifted_cell_duration_ka[in_cat])
    mapped <- matrix(in_cat[event_idx], nrow = analysis_config$n_enrichment_iter)
    uni_values[[j]][sidx, ] <- safe_r(rowSums(mapped), tc)
  }
  for (j in seq_len(nrow(bi_specs))) {
    in_cat <- shifted$RSL_State == "Mid" &
      shifted[[bi_specs$y_var[j]]] == bi_specs$y_level[j]
    tc <- sum(shifted_cell_duration_ka[in_cat])
    mapped <- matrix(in_cat[event_idx], nrow = analysis_config$n_enrichment_iter)
    bi_values[[j]][sidx, ] <- safe_r(rowSums(mapped), tc)
  }
}

zero_idx <- which(steps == 0L)
uni_zero <- vapply(uni_values, function(x) mean(x[zero_idx, ]), numeric(1))
bi_zero <- vapply(bi_values, function(x) mean(x[zero_idx, ]), numeric(1))
if (max(abs(uni_zero - uni_specs$R_enrichment)) > 1e-10 ||
    max(abs(bi_zero - bi_specs$R_enrichment)) > 1e-10) {
  stop("Zero-shift R does not match the canonical enrichment summaries.")
}

uni_by_offset <- do.call(rbind, lapply(seq_len(nrow(uni_specs)), function(j) {
  summaries <- do.call(rbind, lapply(seq_along(steps), function(i) {
    summarise_draws(uni_values[[j]][i, ])
  }))
  data.frame(
    canonical_key = uni_specs$canonical_key[j],
    variable = uni_specs$variable[j],
    level = uni_specs$level[j],
    Tc_ka = uni_specs$Tc_ka[j],
    R_enrichment = uni_specs$R_enrichment[j],
    shift_ka = shift_ka,
    is_zero_shift_reference = steps == 0L,
    summaries,
    stringsAsFactors = FALSE
  )
}))
bi_by_offset <- do.call(rbind, lapply(seq_len(nrow(bi_specs)), function(j) {
  summaries <- do.call(rbind, lapply(seq_along(steps), function(i) {
    summarise_draws(bi_values[[j]][i, ])
  }))
  data.frame(
    canonical_key = bi_specs$canonical_key[j],
    RSL_level = bi_specs$RSL_level[j],
    y_var = bi_specs$y_var[j],
    y_level = bi_specs$y_level[j],
    Tc_ka = bi_specs$Tc_ka[j],
    R_enrichment = bi_specs$R_enrichment[j],
    shift_ka = shift_ka,
    is_zero_shift_reference = steps == 0L,
    summaries,
    stringsAsFactors = FALSE
  )
}))

uni_summary <- do.call(rbind, lapply(seq_len(nrow(uni_specs)), function(j) {
  do.call(rbind, lapply(analysis_config$shift_maxima_ka, function(max_ka) {
    keep <- shift_ka != 0 & abs(shift_ka) <= max_ka
    data.frame(
      canonical_key = uni_specs$canonical_key[j],
      variable = uni_specs$variable[j],
      level = uni_specs$level[j],
      Tc_ka = uni_specs$Tc_ka[j],
      R_enrichment = uni_specs$R_enrichment[j],
      summarise_range(
        uni_values[[j]][keep, , drop = FALSE],
        uni_specs$R_enrichment[j],
        max_ka,
        sum(keep)
      ),
      stringsAsFactors = FALSE
    )
  }))
}))
bi_summary <- do.call(rbind, lapply(seq_len(nrow(bi_specs)), function(j) {
  do.call(rbind, lapply(analysis_config$shift_maxima_ka, function(max_ka) {
    keep <- shift_ka != 0 & abs(shift_ka) <= max_ka
    data.frame(
      canonical_key = bi_specs$canonical_key[j],
      RSL_level = bi_specs$RSL_level[j],
      y_var = bi_specs$y_var[j],
      y_level = bi_specs$y_level[j],
      Tc_ka = bi_specs$Tc_ka[j],
      R_enrichment = bi_specs$R_enrichment[j],
      summarise_range(
        bi_values[[j]][keep, , drop = FALSE],
        bi_specs$R_enrichment[j],
        max_ka,
        sum(keep)
      ),
      stringsAsFactors = FALSE
    )
  }))
}))

write_result_csv(
  uni_summary,
  "circular_shift_robustness_univariate.csv"
)
write_result_csv(
  bi_summary,
  "circular_shift_robustness_midRSL_bivariate.csv"
)
write_result_csv(
  uni_by_offset,
  "circular_shift_by_offset_univariate.csv"
)
write_result_csv(
  bi_by_offset,
  "circular_shift_by_offset_midRSL_bivariate.csv"
)

message("Exhaustive bounded circular-shift summaries complete.")
