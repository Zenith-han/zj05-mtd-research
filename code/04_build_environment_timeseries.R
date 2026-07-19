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

raw_env <- read_csv_required(path("data", "environment_proxy_metrics.csv"))
required <- c(
  "age_ka",
  "RSL_State_Value",
  "EASM_Intensity_Value",
  "RSL_RateofChange_dka_h1.2",
  "EASM_RateofChange_dka_h1.2"
)
if (!all(required %in% names(raw_env))) {
  stop("Environment input is missing: ", paste(setdiff(required, names(raw_env)), collapse = ", "))
}
raw_env <- raw_env[order(raw_env$age_ka), ]
if (min(raw_env$age_ka, na.rm = TRUE) > analysis_config$analysis_start_ka ||
    max(raw_env$age_ka, na.rm = TRUE) < analysis_config$analysis_end_ka) {
  stop("Environment input does not span the configured 0-236 ka interval.")
}

grid <- data.frame(age_ka = analysis_grid_ka)
interp <- function(column) {
  stats::approx(raw_env$age_ka, raw_env[[column]], grid$age_ka, rule = 2)$y
}
grid$RSL_State_Value <- interp("RSL_State_Value")
grid$EASM_Intensity_Value <- interp("EASM_Intensity_Value")

grid$RSL_RateofChange_dka_h1.2 <- environment_rate_series(
  raw_env,
  "RSL_RateofChange_dka_h1.2",
  grid$age_ka
)
grid$EASM_RateofChange_dka_h1.2 <- environment_rate_series(
  raw_env,
  "EASM_RateofChange_dka_h1.2",
  grid$age_ka
)
grid$RSL_Variability_10ka_halfwin <- rolling_sd_grid(
  grid$RSL_State_Value,
  analysis_config$variability_half_window_points
)
grid$EASM_Variability_10ka_halfwin <- rolling_sd_grid(
  grid$EASM_Intensity_Value,
  analysis_config$variability_half_window_points
)

rsl <- classify_tertile(grid$RSL_State_Value, c("Low", "Mid", "High"))
easm <- classify_tertile(grid$EASM_Intensity_Value, c("Weak", "Moderate", "Strong"))
rsl_sd <- classify_tertile(grid$RSL_Variability_10ka_halfwin, c("Low", "Moderate", "High"))
easm_sd <- classify_tertile(grid$EASM_Variability_10ka_halfwin, c("Low", "Moderate", "High"))
rsl_rate <- classify_signed_rate(
  grid$RSL_RateofChange_dka_h1.2,
  c("Slow Rise", "Moderate Rise", "Rapid Rise"),
  c("Slow Fall", "Moderate Fall", "Rapid Fall")
)
easm_rate <- classify_signed_rate(
  grid$EASM_RateofChange_dka_h1.2,
  c("Slow Strengthening", "Moderate Strengthening", "Rapid Strengthening"),
  c("Slow Weakening", "Moderate Weakening", "Rapid Weakening")
)

grid$RSL_State <- rsl$values
grid$EASM_Intensity <- easm$values
grid$RSL_Variability <- rsl_sd$values
grid$EASM_Variability <- easm_sd$values
grid$RSL_RateofChange <- rsl_rate$values
grid$EASM_RateofChange <- easm_rate$values

vars <- c(
  "RSL_State", "RSL_Variability", "RSL_RateofChange",
  "EASM_Intensity", "EASM_Variability", "EASM_RateofChange"
)
if (anyNA(grid[, vars])) stop("The 0-236 ka environmental classification contains NA values.")

write_result_csv(grid, "categorized_timeseries.csv")
