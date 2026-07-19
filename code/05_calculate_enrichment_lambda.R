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
if (nrow(env) != 1181L || min(env$age_ka) != 0 || max(env$age_ka) != 236) {
  stop("Environmental grid is not the required 1,181-point 0-236 ka grid.")
}
weights <- grid_exposure_weights(
  env$age_ka,
  analysis_config$analysis_start_ka,
  analysis_config$analysis_end_ka
)
total_time_ka <- sum(weights)
dt_ka <- analysis_config$step_ka

posterior <- read_result_csv_required(file.path(
  "bacon", "event_age_posterior_draws_wide.csv"
))
events <- read_csv_required(path("data", "mtd_events.csv"))
event_ids <- events$event_id
if (length(event_ids) != analysis_config$n_events ||
    nrow(posterior) != analysis_config$posterior_iterations ||
    !all(c("iteration", event_ids) %in% names(posterior))) {
  stop("The retained Bacon iteration table has invalid dimensions or event columns.")
}

set.seed(analysis_config$enrichment_seed)
sample_rows <- sample(seq_len(nrow(posterior)), analysis_config$n_enrichment_iter, replace = FALSE)
ages_mat <- t(as.matrix(posterior[sample_rows, event_ids, drop = FALSE]))
nearest_index <- function(age_ka) {
  idx <- round((age_ka - env$age_ka[1]) / dt_ka) + 1L
  pmin(pmax(idx, 1L), nrow(env))
}
summarise_draws <- function(x) {
  c(
    mean = mean(x),
    p05 = stats::quantile(x, 0.05, names = FALSE),
    p50 = stats::quantile(x, 0.50, names = FALSE),
    p95 = stats::quantile(x, 0.95, names = FALSE),
    prob_gt_1 = mean(x > 1)
  )
}

calc_counts <- function(in_category) {
  vapply(seq_len(ncol(ages_mat)), function(j) {
    sum(in_category[nearest_index(ages_mat[, j])])
  }, numeric(1))
}

univariate_rows <- list()
univariate_draw_rows <- list()
k <- 1L
for (var in vars) {
  levels <- unique(env[[var]])
  levels <- levels[!is.na(levels)]
  tc_values <- vapply(levels, function(level) sum(weights[env[[var]] == level]), numeric(1))
  if (abs(sum(tc_values) - total_time_ka) > 1e-10) stop("Univariate Tc does not sum to 236 ka: ", var)
  for (level in levels) {
    in_category <- env[[var]] == level
    tc <- sum(weights[in_category])
    counts <- calc_counts(in_category)
    lambda <- counts / tc
    r <- (counts / analysis_config$n_events) / (tc / total_time_ka)
    rs <- summarise_draws(r)
    ls <- summarise_draws(lambda)
    univariate_rows[[k]] <- data.frame(
      var = var, level = level, Tc_ka = tc,
      R_mean = rs["mean"], R_p05 = rs["p05"], R_p50 = rs["p50"], R_p95 = rs["p95"],
      P_R_gt_1 = rs["prob_gt_1"],
      lambda_mean = ls["mean"], lambda_p05 = ls["p05"],
      lambda_p50 = ls["p50"], lambda_p95 = ls["p95"],
      row.names = NULL
    )
    univariate_draw_rows[[k]] <- data.frame(
      variable = var, level = level,
      iteration = seq_len(analysis_config$n_enrichment_iter),
      bacon_iteration = posterior$iteration[sample_rows],
      count = counts, Tc_ka = tc, R = r, lambda = lambda
    )
    k <- k + 1L
  }
}
univariate <- do.call(rbind, univariate_rows)
univariate_draws <- do.call(rbind, univariate_draw_rows)

bivariate_rows <- list()
midrsl_draw_rows <- list()
pairs <- utils::combn(vars, 2, simplify = FALSE)
k <- 1L
mid_k <- 1L
for (pair in pairs) {
  var_a <- pair[1]
  var_b <- pair[2]
  for (a in unique(env[[var_a]])) {
    for (b in unique(env[[var_b]])) {
      in_cell <- env[[var_a]] == a & env[[var_b]] == b
      tc <- sum(weights[in_cell])
      if (tc <= 0) next
      counts <- calc_counts(in_cell)
      lambda <- counts / tc
      r <- (counts / analysis_config$n_events) / (tc / total_time_ka)
      rs <- summarise_draws(r)
      ls <- summarise_draws(lambda)
      bivariate_rows[[k]] <- data.frame(
        varA = var_a, A_level = a, varB = var_b, B_level = b, Tc_ka = tc,
        R_mean = rs["mean"], R_p05 = rs["p05"], R_p50 = rs["p50"], R_p95 = rs["p95"],
        P_R_gt_1 = rs["prob_gt_1"],
        lambda_mean = ls["mean"], lambda_p05 = ls["p05"],
        lambda_p50 = ls["p50"], lambda_p95 = ls["p95"],
        row.names = NULL
      )
      if (var_a == "RSL_State" && a == "Mid") {
        midrsl_draw_rows[[mid_k]] <- data.frame(
          RSL_level = a,
          y_var = var_b,
          y_level = b,
          iteration = seq_len(analysis_config$n_enrichment_iter),
          bacon_iteration = posterior$iteration[sample_rows],
          count = counts,
          Tc_ka = tc,
          R = r,
          lambda = lambda
        )
        mid_k <- mid_k + 1L
      }
      k <- k + 1L
    }
  }
  pair_summary <- do.call(rbind, bivariate_rows[vapply(
    bivariate_rows,
    function(x) x$varA[1] == var_a && x$varB[1] == var_b,
    logical(1)
  )])
  if (abs(sum(pair_summary$Tc_ka) - total_time_ka) > 1e-10) {
    stop("Bivariate Tc does not sum to 236 ka: ", var_a, " x ", var_b)
  }
}
bivariate <- do.call(rbind, bivariate_rows)
midrsl <- bivariate[bivariate$varA == "RSL_State" & bivariate$A_level == "Mid", ]
midrsl_draws <- do.call(rbind, midrsl_draw_rows)
if (nrow(midrsl) != 21L ||
    nrow(midrsl_draws) != 21L * analysis_config$n_enrichment_iter) {
  stop("Unexpected Mid-RSL bivariate summary or draw dimensions.")
}

write_result_csv(univariate, "univariate_enrichment_summary.csv")
write_result_csv(bivariate, "bivariate_enrichment_summary.csv")
write_result_csv(univariate_draws, "univariate_draws_ALL.csv")
write_result_csv(midrsl_draws, "midRSL_bivariate_draws_ALL.csv")
