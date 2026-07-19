classify_tertile <- function(x, labels) {
  qs <- stats::quantile(
    x,
    probs = c(1 / 3, 2 / 3),
    na.rm = TRUE,
    type = 7,
    names = FALSE
  )
  list(
    values = as.character(cut(
      x,
      breaks = c(-Inf, qs, Inf),
      labels = labels,
      include.lowest = TRUE
    )),
    thresholds = qs
  )
}

classify_signed_rate <- function(x, rise_labels, fall_labels) {
  out <- rep(NA_character_, length(x))
  pos <- which(is.finite(x) & x > 0)
  neg <- which(is.finite(x) & x < 0)
  zero <- which(is.finite(x) & x == 0)
  if (!length(pos) || !length(neg)) {
    stop("Signed-rate classification requires both positive and negative values.")
  }
  pos_q <- stats::quantile(abs(x[pos]), c(1 / 3, 2 / 3), na.rm = TRUE, names = FALSE)
  neg_q <- stats::quantile(abs(x[neg]), c(1 / 3, 2 / 3), na.rm = TRUE, names = FALSE)
  out[pos] <- as.character(cut(
    abs(x[pos]), c(-Inf, pos_q, Inf),
    labels = rise_labels, include.lowest = TRUE
  ))
  out[neg] <- as.character(cut(
    abs(x[neg]), c(-Inf, neg_q, Inf),
    labels = fall_labels, include.lowest = TRUE
  ))
  if (length(zero)) {
    stop("Exact zero rate values are not classifiable under the six-class scheme.")
  }
  list(
    values = out,
    positive_thresholds = pos_q,
    negative_abs_thresholds = neg_q
  )
}

grid_exposure_weights <- function(age_ka, start_ka, end_ka) {
  if (is.unsorted(age_ka, strictly = TRUE)) stop("Age grid must be strictly increasing.")
  if (age_ka[1] != start_ka || tail(age_ka, 1) != end_ka) {
    stop("Age grid does not match configured analysis bounds.")
  }
  boundaries <- c(
    start_ka,
    (age_ka[-length(age_ka)] + age_ka[-1]) / 2,
    end_ka
  )
  weights <- diff(boundaries)
  if (any(weights <= 0) || abs(sum(weights) - (end_ka - start_ka)) > 1e-12) {
    stop("Exposure weights do not exactly cover the analysis duration.")
  }
  weights
}

rolling_sd_grid <- function(value, half_window_points = 25L) {
  n <- length(value)
  if (n < 2L || half_window_points < 1L) {
    stop("A rolling standard deviation requires at least two values and a positive half-window.")
  }
  vapply(seq_len(n), function(i) {
    idx <- seq.int(max(1L, i - half_window_points), min(n, i + half_window_points))
    stats::sd(value[idx], na.rm = TRUE)
  }, numeric(1))
}

environment_rate_series <- function(
    raw_environment,
    rate_column,
    target_age_ka) {
  if (!is.data.frame(raw_environment) || !"age_ka" %in% names(raw_environment)) {
    stop("Environment input is missing the required age_ka column.")
  }
  if (!is.character(rate_column) || length(rate_column) != 1L ||
      is.na(rate_column) || !nzchar(rate_column)) {
    stop("rate_column must name one required precomputed rate column.")
  }
  if (!rate_column %in% names(raw_environment)) {
    stop("Environment input is missing required precomputed rate column: ", rate_column)
  }

  age_ka <- raw_environment$age_ka
  rate <- raw_environment[[rate_column]]
  if (!is.numeric(age_ka) || !is.numeric(rate) ||
      length(age_ka) != length(rate) || length(age_ka) < 2L) {
    stop("Environment ages and precomputed rates must be numeric vectors of equal length.")
  }
  if (any(!is.finite(age_ka)) || is.unsorted(age_ka, strictly = TRUE)) {
    stop("Environment ages must be finite and strictly increasing.")
  }
  if (any(is.infinite(rate)) || sum(is.finite(rate)) < 2L) {
    stop("The required precomputed rate column must contain at least two finite values.")
  }
  finite_span <- range(which(is.finite(rate)))
  if (any(!is.finite(rate[seq.int(finite_span[1], finite_span[2])]))) {
    stop("The required precomputed rate column contains an internal missing value.")
  }
  if (!is.numeric(target_age_ka) || !length(target_age_ka) ||
      any(!is.finite(target_age_ka))) {
    stop("Target ages must be a non-empty finite numeric vector.")
  }

  stats::approx(
    age_ka,
    rate,
    target_age_ka,
    rule = 2,
    ties = "ordered",
    na.rm = TRUE
  )$y
}
