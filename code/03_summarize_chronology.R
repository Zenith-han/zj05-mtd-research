script_dir <- function() {
  args <- commandArgs(FALSE)
  hit <- grep("^--file=", args, value = TRUE)
  if (length(hit) > 0) return(dirname(normalizePath(sub("^--file=", "", hit[1]))))
  if (basename(getwd()) == "code") return(normalizePath(getwd()))
  normalizePath(file.path(getwd(), "code"), mustWork = TRUE)
}
source(file.path(script_dir(), "bacon_helpers.R"))
source(file.path(script_dir(), "analysis_config.R"))

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

summarize_rate_draws <- function(values, include_coverage = TRUE) {
  ok <- is.finite(values)
  q <- if (sum(ok) > 0L) {
    as.numeric(stats::quantile(
      values[ok],
      probs = c(0.025, 0.5, 0.975),
      names = FALSE
    ))
  } else {
    rep(NA_real_, 3L)
  }
  summary <- c(
    mean_cm_per_kyr = if (sum(ok) > 0L) mean(values[ok]) else NA_real_,
    q2.5_cm_per_kyr = q[1],
    median_cm_per_kyr = q[2],
    q97.5_cm_per_kyr = q[3],
    n_draws_used = sum(ok)
  )
  if (include_coverage) summary <- c(summary, coverage = mean(ok))
  summary
}

stage_rates <- bacon_posterior_interval_rates(
  set,
  stage_sedimentation_intervals$age_left_ka,
  stage_sedimentation_intervals$age_right_ka
)
stage_statistics <- t(apply(stage_rates, 2, summarize_rate_draws))
stage_summary <- data.frame(
  stage_sedimentation_intervals,
  duration_kyr = stage_sedimentation_intervals$age_right_ka -
    stage_sedimentation_intervals$age_left_ka,
  stage_statistics,
  check.names = FALSE,
  stringsAsFactors = FALSE
)
write_result_csv(
  stage_summary,
  file.path("bacon", "stage_sedimentation_rate_posterior_summary.csv")
)

age_edges_ka <- sedimentation_rate_10ka_edges_ka
lsr_draws <- bacon_posterior_interval_rates(
  set,
  head(age_edges_ka, -1),
  tail(age_edges_ka, -1)
)
lsr_rows <- lapply(seq_len(ncol(lsr_draws)), function(i) {
  rates <- lsr_draws[, i]
  ok <- is.finite(rates)
  q <- if (sum(ok) > 0L) {
    as.numeric(stats::quantile(
      rates[ok],
      probs = c(0.025, 0.25, 0.5, 0.75, 0.975),
      na.rm = TRUE
    ))
  } else {
    rep(NA_real_, 5L)
  }
  data.frame(
    age_left_ka = age_edges_ka[i],
    age_right_ka = age_edges_ka[i + 1L],
    LSR_p02_5 = q[1],
    LSR_p25 = q[2],
    LSR_p50 = q[3],
    LSR_p75 = q[4],
    LSR_p97_5 = q[5],
    n_draws_used = sum(ok),
    coverage = mean(ok),
    stringsAsFactors = FALSE
  )
})
lsr_step <- do.call(rbind, lsr_rows)
write_result_csv(
  lsr_step,
  file.path("bacon", "ZJ05_LSR_step_10ka_pseudoPosterior.csv")
)

wide <- read_result_csv_required(file.path(
  "bacon", "event_age_posterior_draws_wide.csv"
))
events <- event_table()
event_ids <- events$event_id
if (!identical(as.integer(wide$iteration), seq_len(nrow(wide))) ||
    nrow(wide) != nrow(lsr_draws) || !all(event_ids %in% names(wide))) {
  stop("Event-age posterior rows do not align with the retained Bacon curves.")
}
event_age_draws <- as.matrix(wide[, event_ids, drop = FALSE])
linked_event_rates <- link_same_row_event_rates(
  lsr_draws,
  event_age_draws,
  age_edges_ka
)

event_statistics <- t(apply(
  linked_event_rates,
  2,
  summarize_rate_draws,
  include_coverage = FALSE
))
event_summary <- data.frame(
  summary_level = "event",
  summary_id = event_ids,
  member_events = event_ids,
  event_count = 1L,
  window_rule = "10-kyr bins; left-closed/right-open",
  aggregation_rule = "same-row event age and Bacon-curve rate",
  event_statistics,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

group_rates <- vapply(
  mtd_sedimentation_rate_groups,
  function(members) rowMeans(linked_event_rates[, members, drop = FALSE]),
  numeric(nrow(linked_event_rates))
)
group_statistics <- t(apply(
  group_rates,
  2,
  summarize_rate_draws,
  include_coverage = FALSE
))
group_summary <- data.frame(
  summary_level = "group",
  summary_id = names(mtd_sedimentation_rate_groups),
  member_events = vapply(
    mtd_sedimentation_rate_groups,
    paste,
    collapse = ",",
    character(1)
  ),
  event_count = lengths(mtd_sedimentation_rate_groups),
  window_rule = "10-kyr bins; left-closed/right-open",
  aggregation_rule = "equal-weight arithmetic mean within posterior iteration",
  group_statistics,
  check.names = FALSE,
  stringsAsFactors = FALSE
)
linked_summary <- rbind(event_summary, group_summary)
write_result_csv(
  linked_summary,
  file.path("bacon", "mtd_linked_sedimentation_rate_posterior_summary.csv")
)
