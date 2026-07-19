script_dir <- function() {
  args <- commandArgs(FALSE)
  hit <- grep("^--file=", args, value = TRUE)
  if (length(hit) > 0) return(dirname(normalizePath(sub("^--file=", "", hit[1]))))
  if (basename(getwd()) == "code") return(normalizePath(getwd()))
  normalizePath(file.path(getwd(), "code"), mustWork = TRUE)
}
source(file.path(script_dir(), "00_paths.R"))
source(file.path(script_dir(), "bacon_model_config.R"))
source(file.path(script_dir(), "bacon_control_utils.R"))

selected_results_root <- results_root()
assert_isolated_bacon_build_root(root_dir, selected_results_root)
rebuilt_model_file <- bacon_model_file(
  root_dir,
  "rebuilt",
  must_exist = FALSE,
  results_root = selected_results_root
)
core_dir <- dirname(rebuilt_model_file)
build_root <- dirname(core_dir)

require_bacon <- function() {
  if (!requireNamespace("rbacon", quietly = TRUE)) {
    stop("The rbacon package is required.")
  }
  version <- as.character(utils::packageVersion("rbacon"))
  if (version != "3.5.2") stop("Expected rbacon 3.5.2; found ", version)
}
require_bacon()

bacon_params <- bacon_model_config

controls <- prepare_bacon_controls(
  read_csv_required(path("data", "age_control.csv")),
  expected_rows = 31L
)
age_reversals <- central_age_reversals(controls)
if (length(age_reversals) > 0) {
  reversal_pairs <- vapply(
    age_reversals,
    function(i) paste0(controls$LabID[i], " -> ", controls$LabID[i + 1L]),
    character(1)
  )
  warning(
    "Central ages decrease with depth at: ",
    paste(reversal_pairs, collapse = "; "),
    ". The input remains ordered by depth; review these controls before interpreting the model."
  )
}

core <- bacon_params$core
if (!dir.exists(core_dir)) dir.create(core_dir, recursive = TRUE, showWarnings = FALSE)
utils::write.csv(controls, file.path(core_dir, paste0(core, ".csv")), row.names = FALSE)

depths <- seq(
  bacon_params$d.min,
  bacon_params$d.max,
  by = bacon_params$d.by
)
set.seed(bacon_params$seed)
info <- rbacon::Bacon(
  core = core,
  coredir = build_root,
  thick = bacon_params$thick,
  d.min = min(depths),
  d.max = max(depths),
  d.by = bacon_params$d.by,
  depths = depths,
  acc.mean = bacon_params$acc.mean,
  acc.shape = bacon_params$acc.shape,
  mem.mean = bacon_params$mem.mean,
  mem.strength = bacon_params$mem.strength,
  ssize = bacon_params$ssize,
  burnin = bacon_params$burnin,
  seed = bacon_params$seed,
  normal = bacon_params$normal,
  suggest = FALSE,
  ask = FALSE,
  run = TRUE,
  plot.pdf = TRUE,
  save.info = TRUE,
  save.elbowages = TRUE,
  verbose = TRUE
)

info <- normalize_bacon_model_output(info)
assert_bacon_model_configuration(info)
saveRDS(info, rebuilt_model_file)
message("Wrote Bacon chronology model build outputs to ", core_dir)
