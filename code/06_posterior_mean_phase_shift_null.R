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
if (!all(c("age_ka", vars) %in% names(env))) {
  stop("Environmental grid is missing required classified columns.")
}
if (nrow(env) != length(analysis_grid_ka) ||
    abs(min(env$age_ka) - analysis_config$analysis_start_ka) > 1e-12 ||
    abs(max(env$age_ka) - analysis_config$analysis_end_ka) > 1e-12) {
  stop("Environmental grid is not the required 0-237 ka, 0.2-ka grid.")
}
weights <- grid_exposure_weights(
  env$age_ka,
  analysis_config$analysis_start_ka,
  analysis_config$analysis_end_ka
)
total_time_ka <- sum(weights)

posterior <- read_result_csv_required(file.path(
  "bacon", "event_age_posterior_draws_wide.csv"
))
events <- read_csv_required(path("data", "mtd_events.csv"))
event_ids <- events$event_id
if (length(event_ids) != analysis_config$n_events ||
    nrow(posterior) != analysis_config$posterior_iterations ||
    !all(c("iteration", event_ids) %in% names(posterior))) {
  stop("Bacon posterior event-age table has invalid dimensions or event columns.")
}

main_uni <- read_result_csv_required("univariate_enrichment_summary.csv")
main_bi <- read_result_csv_required("bivariate_enrichment_summary.csv")
main_mid <- main_bi[main_bi$varA == "RSL_State" & main_bi$A_level == "Mid", ]
if (nrow(main_mid) != 21L) stop("Expected 21 Mid-RSL bivariate rows.")

uni_draws <- read_result_csv_required("univariate_draws_ALL.csv")
set.seed(analysis_config$enrichment_seed)
sample_rows <- sample(seq_len(nrow(posterior)), analysis_config$n_enrichment_iter, replace = FALSE)
selected_iterations <- posterior$iteration[sample_rows]
check_draws <- uni_draws[
  uni_draws$variable == "RSL_State" & uni_draws$level == "High",
  c("iteration", "bacon_iteration")
]
check_draws <- check_draws[order(check_draws$iteration), ]
if (!identical(as.integer(check_draws$bacon_iteration), as.integer(selected_iterations))) {
  stop("Reconstructed seed=1 posterior rows do not match submitted draw files.")
}

selected_ages <- as.matrix(posterior[sample_rows, event_ids, drop = FALSE])
storage.mode(selected_ages) <- "double"
age_vec <- as.vector(selected_ages)
n_draw <- nrow(selected_ages)
n_events <- length(event_ids)

wrap_age <- function(age_ka, shift_ka) {
  duration <- analysis_config$analysis_end_ka - analysis_config$analysis_start_ka
  ((age_ka - analysis_config$analysis_start_ka + shift_ka) %% duration) +
    analysis_config$analysis_start_ka
}

nearest_index <- function(age_ka) {
  idx <- round((age_ka - env$age_ka[1]) / analysis_config$step_ka) + 1L
  pmin(pmax(idx, 1L), nrow(env))
}

make_codes <- function(x) {
  lv <- unique(x[!is.na(x)])
  list(code = match(x, lv), levels = lv)
}

coded <- lapply(vars, function(v) make_codes(env[[v]]))
names(coded) <- vars

level_code <- function(var, level) {
  out <- match(level, coded[[var]]$levels)
  if (is.na(out)) stop("Level not found in environmental grid: ", var, " / ", level)
  out
}

weighted_duration <- function(mask) sum(weights[mask], na.rm = TRUE)

uni_specs <- data.frame(
  spec_id = paste("univariate", main_uni$var, main_uni$level, sep = "|"),
  analysis_type = "univariate",
  condition_variable = "",
  condition_level = "",
  variable = main_uni$var,
  level = main_uni$level,
  Tc_ka = main_uni$Tc_ka,
  R_mean = main_uni$R_mean,
  stringsAsFactors = FALSE
)

mid_specs <- data.frame(
  spec_id = paste("midRSL_bivariate", "RSL_State", "Mid", main_mid$varB, main_mid$B_level, sep = "|"),
  analysis_type = "midRSL_bivariate",
  condition_variable = "RSL_State",
  condition_level = "Mid",
  variable = main_mid$varB,
  level = main_mid$B_level,
  Tc_ka = main_mid$Tc_ka,
  R_mean = main_mid$R_mean,
  stringsAsFactors = FALSE
)

specs <- rbind(uni_specs, mid_specs)
row.names(specs) <- NULL

validate_exposure <- function() {
  for (i in seq_len(nrow(main_uni))) {
    mask <- env[[main_uni$var[i]]] == main_uni$level[i]
    if (abs(weighted_duration(mask) - main_uni$Tc_ka[i]) > 1e-10) {
      stop("Univariate exposure duration mismatch: ", main_uni$var[i], " / ", main_uni$level[i])
    }
  }
  for (i in seq_len(nrow(main_mid))) {
    mask <- env$RSL_State == "Mid" & env[[main_mid$varB[i]]] == main_mid$B_level[i]
    if (abs(weighted_duration(mask) - main_mid$Tc_ka[i]) > 1e-10) {
      stop("Mid-RSL bivariate exposure duration mismatch: ", main_mid$varB[i], " / ", main_mid$B_level[i])
    }
  }
}
validate_exposure()

compute_stats_for_index <- function(idx) {
  total_event_count <- length(idx)
  if (total_event_count != n_draw * n_events) stop("Unexpected event-index vector length.")

  out <- data.frame(
    spec_id = specs$spec_id,
    mean_count_numerator = NA_real_,
    mean_count_denominator = n_events,
    R = NA_real_,
    stringsAsFactors = FALSE
  )

  for (i in seq_len(nrow(main_uni))) {
    code <- level_code(main_uni$var[i], main_uni$level[i])
    count <- sum(coded[[main_uni$var[i]]]$code[idx] == code, na.rm = TRUE)
    out$mean_count_numerator[i] <- count / n_draw
    out$R[i] <- (count / total_event_count) / (main_uni$Tc_ka[i] / total_time_ka)
  }

  mid_code <- level_code("RSL_State", "Mid")
  in_mid <- coded$RSL_State$code[idx] == mid_code
  offset <- nrow(main_uni)
  for (i in seq_len(nrow(main_mid))) {
    code <- level_code(main_mid$varB[i], main_mid$B_level[i])
    joint <- sum(in_mid & coded[[main_mid$varB[i]]]$code[idx] == code, na.rm = TRUE)
    row_i <- offset + i
    out$mean_count_numerator[row_i] <- joint / n_draw
    out$R[row_i] <- (joint / total_event_count) / (main_mid$Tc_ka[i] / total_time_ka)
  }
  out
}

base_idx <- nearest_index(age_vec)
observed_stats <- compute_stats_for_index(base_idx)
obs_match <- match(specs$spec_id, observed_stats$spec_id)
max_obs_diff <- max(abs(observed_stats$R[obs_match] - specs$R_mean))
if (!is.finite(max_obs_diff) || max_obs_diff > 1e-10) {
  stop("Rebuilt R does not match source R_mean. Max diff = ", max_obs_diff)
}

phase_steps <- 0:as.integer(round(
  (analysis_config$analysis_end_ka - analysis_config$analysis_start_ka) / analysis_config$step_ka
) - 1L)
shift_values_ka <- phase_steps * analysis_config$step_ka
if (length(shift_values_ka) != 1185L || shift_values_ka[1] != 0 ||
    abs(tail(shift_values_ka, 1) - 236.8) > 1e-12) {
  stop("Exact phase grid is not 0, 0.2, ..., 236.8 ka.")
}

message("Running posterior-mean exact circular-shift null over ", length(shift_values_ka), " phases.")
null_r <- matrix(NA_real_, nrow = length(shift_values_ka), ncol = nrow(specs))
for (s in seq_along(shift_values_ka)) {
  idx <- nearest_index(wrap_age(age_vec, shift_values_ka[s]))
  null_r[s, ] <- compute_stats_for_index(idx)$R
  if (s %% 100L == 0L) {
    message("Completed ", s, " / ", length(shift_values_ka), " exact phases")
  }
}
if (max(abs(null_r[1, ] - specs$R_mean)) > 1e-10) {
  stop("Zero-phase R does not match the observed posterior-mean R.")
}

phase_long <- do.call(rbind, lapply(seq_len(nrow(specs)), function(i) {
  data.frame(
    spec_id = specs$spec_id[i],
    analysis_type = specs$analysis_type[i],
    condition_variable = specs$condition_variable[i],
    condition_level = specs$condition_level[i],
    variable = specs$variable[i],
    level = specs$level[i],
    Tc_ka = specs$Tc_ka[i],
    R_observed = specs$R_mean[i],
    shift_ka = shift_values_ka,
    R_phase_mean = null_r[, i],
    stringsAsFactors = FALSE
  )
}))
row.names(phase_long) <- NULL

observed <- merge(
  specs,
  observed_stats[, c("spec_id", "R")],
  by = "spec_id",
  all.x = TRUE,
  sort = FALSE
)
names(observed)[names(observed) == "R"] <- "rebuilt_R_mean"

summary_rows <- vector("list", nrow(specs))
for (i in seq_len(nrow(specs))) {
  obs <- observed[observed$spec_id == specs$spec_id[i], ]
  vals <- null_r[, i]
  if (nrow(obs) != 1L || length(vals) != length(shift_values_ka)) {
    stop("Invalid summary dimensions for ", specs$spec_id[i])
  }
  observed_threshold <- obs$rebuilt_R_mean
  p_enrichment <- sum(vals >= observed_threshold) / length(vals)
  p_depletion <- sum(vals <= observed_threshold) / length(vals)
  effect_direction <- if (obs$R_mean >= 1) "enrichment" else "depletion"
  p_one_sided <- if (effect_direction == "enrichment") p_enrichment else p_depletion

  summary_rows[[i]] <- data.frame(
    analysis_type = obs$analysis_type,
    condition_variable = obs$condition_variable,
    condition_level = obs$condition_level,
    variable = obs$variable,
    level = obs$level,
    R_mean = obs$R_mean,
    Tc_ka = obs$Tc_ka,
    null_mean = mean(vals),
    null_p05 = stats::quantile(vals, 0.05, names = FALSE),
    null_p50 = stats::quantile(vals, 0.50, names = FALSE),
    null_p95 = stats::quantile(vals, 0.95, names = FALSE),
    effect_direction = effect_direction,
    p_one_sided = p_one_sided,
    stringsAsFactors = FALSE
  )
}

summary <- do.call(rbind, summary_rows)
if (nrow(summary[summary$analysis_type == "univariate", ]) != 24L) {
  stop("Expected 24 univariate summary rows.")
}
if (nrow(summary[summary$analysis_type == "midRSL_bivariate", ]) != 21L) {
  stop("Expected 21 Mid-RSL bivariate summary rows.")
}
if (anyNA(summary$p_one_sided)) stop("one-sided p values contain NA.")
if (!all(summary$effect_direction == ifelse(summary$R_mean >= 1, "enrichment", "depletion"))) {
  stop("Effect direction does not match R_mean.")
}

invisible(write_result_csv(summary, "posterior_mean_phase_shift_null_summary.csv"))
invisible(write_result_csv(phase_long, "posterior_mean_phase_shift_null_long.csv"))

message("Posterior-mean exact phase-shift null summary complete.")
