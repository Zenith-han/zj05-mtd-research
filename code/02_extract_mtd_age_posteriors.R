script_dir <- function() {
  args <- commandArgs(FALSE)
  hit <- grep("^--file=", args, value = TRUE)
  if (length(hit) > 0) return(dirname(normalizePath(sub("^--file=", "", hit[1]))))
  if (basename(getwd()) == "code") return(normalizePath(getwd()))
  normalizePath(file.path(getwd(), "code"), mustWork = TRUE)
}
source(file.path(script_dir(), "bacon_helpers.R"))
source(file.path(script_dir(), "bacon_model_config.R"))
source(file.path(script_dir(), "analysis_config.R"))
require_bacon()

info_file <- getOption("zj05.bacon_model_file", NULL)
if (is.null(info_file)) {
  info_file <- bacon_model_file(root_dir, "archived", must_exist = TRUE)
} else {
  if (length(info_file) != 1L || !file.exists(info_file)) {
    stop("Selected Bacon model object not found: ", info_file)
  }
  info_file <- normalizePath(info_file, mustWork = TRUE)
}
set <- readRDS(info_file)
assert_bacon_model_configuration(set)

events <- event_table()
posterior <- event_draws_and_summary(set)
draws <- posterior$draws_years / 1000
if (nrow(draws) < bacon_model_config$ssize) {
  stop("The selected Bacon model has fewer retained posterior iterations than configured.")
}
if (nrow(draws) > bacon_model_config$ssize) {
  set.seed(20260619)
  keep <- sort(sample(seq_len(nrow(draws)), bacon_model_config$ssize, replace = FALSE))
  draws <- draws[keep, , drop = FALSE]
}
if (!all(apply(draws, 1, function(x) all(diff(x) > 0)))) {
  stop("At least one retained model iteration violates event-age monotonicity.")
}

wide <- data.frame(iteration = seq_len(nrow(draws)), draws, check.names = FALSE)
write_result_csv(wide, file.path("bacon", "event_age_posterior_draws_wide.csv"))

q <- t(apply(draws, 2, stats::quantile, probs = c(0.025, 0.5, 0.975), na.rm = TRUE))
summary <- data.frame(
  event_id = events$event_id,
  original_top_mbsf = events$original_top_mbsf,
  original_base_mbsf = events$original_base_mbsf,
  event_free_depth_cm = events$event_free_depth_cm,
  q2.5_ka = q[, 1],
  median_ka = q[, 2],
  mean_ka = colMeans(draws),
  q97.5_ka = q[, 3],
  stringsAsFactors = FALSE
)
write_result_csv(summary, "event_age_summary.csv")

official <- set$ranges
if (is.null(official) || ncol(official) != 5L) {
  stop("The selected Bacon model object does not contain a valid 1-cm age-depth table.")
}
names(official) <- c("depth_cm", "q2.5_year", "q97.5_year", "median_year", "mean_year")
write_result_csv(official, file.path("bacon", "age_depth_model_1cm.csv"))

idx <- match(events$event_free_depth_cm, official$depth_cm)
if (anyNA(idx)) stop("Event depths are missing from the official 1-cm Bacon summary.")
validation <- data.frame(
  event_id = events$event_id,
  depth_cm = events$event_free_depth_cm,
  draw_q2.5_year = q[, 1] * 1000,
  official_q2.5_year = official$q2.5_year[idx],
  draw_median_year = q[, 2] * 1000,
  official_median_year = official$median_year[idx],
  draw_mean_year = colMeans(draws) * 1000,
  official_mean_year = official$mean_year[idx],
  draw_q97.5_year = q[, 3] * 1000,
  official_q97.5_year = official$q97.5_year[idx],
  stringsAsFactors = FALSE
)
validation$max_abs_difference_year <- apply(abs(cbind(
  validation$draw_q2.5_year - validation$official_q2.5_year,
  validation$draw_median_year - validation$official_median_year,
  validation$draw_mean_year - validation$official_mean_year,
  validation$draw_q97.5_year - validation$official_q97.5_year
)), 1, max)
if (max(validation$max_abs_difference_year) > 2) {
  stop("Direct posterior summaries differ from Bacon's official rounded summary by >2 years.")
}
message("Wrote 18,000 complete event-age sets from the retained Bacon model iterations.")
