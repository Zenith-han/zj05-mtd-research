prepare_bacon_controls <- function(controls, expected_rows = 31L) {
  required_cols <- c("LabID", "age", "error", "depth", "cc", "delta.R", "delta.STD")

  if (!is.data.frame(controls)) {
    stop("Bacon controls must be supplied as a data frame.")
  }
  if (!identical(names(controls), required_cols)) {
    stop(
      "data/age_control.csv must contain columns: ",
      paste(required_cols, collapse = ", ")
    )
  }
  if (nrow(controls) != expected_rows || anyNA(controls)) {
    stop(
      "data/age_control.csv must contain ",
      expected_rows,
      " complete age-control rows."
    )
  }
  if (!is.character(controls$LabID) || any(!nzchar(trimws(controls$LabID)))) {
    stop("LabID values must be non-empty text.")
  }
  if (anyDuplicated(controls$LabID)) {
    stop("LabID values must be unique.")
  }

  numeric_cols <- c("age", "error", "depth", "cc", "delta.R", "delta.STD")
  non_numeric <- numeric_cols[!vapply(controls[numeric_cols], is.numeric, logical(1))]
  if (length(non_numeric) > 0) {
    stop("These Bacon control fields must be numeric: ", paste(non_numeric, collapse = ", "))
  }
  if (any(!is.finite(as.matrix(controls[numeric_cols])))) {
    stop("Bacon control values must be finite.")
  }
  if (any(controls$age <= 0)) stop("Control ages must be greater than zero.")
  if (any(controls$error <= 0)) stop("Control errors must be greater than zero.")
  if (any(controls$depth < 0)) stop("Control depths must be non-negative.")
  if (any(controls$cc != round(controls$cc)) || any(!controls$cc %in% 0:4)) {
    stop("Calibration-curve codes (cc) must be integers from 0 to 4.")
  }
  # rbacon requires dated levels to be listed in non-decreasing depth order;
  # multiple independent controls at the same depth are valid.
  # Age is deliberately not used as a sorting key because age reversals can be
  # legitimate observations that the age-depth model must evaluate.
  controls <- controls[order(controls$depth, controls$LabID, method = "radix"), required_cols]
  rownames(controls) <- NULL
  if (any(diff(controls$depth) < 0)) {
    stop("Control depths must be non-decreasing after sorting.")
  }

  controls
}

central_age_reversals <- function(controls) {
  which(diff(controls$depth) > 0 & diff(controls$age) < 0)
}
