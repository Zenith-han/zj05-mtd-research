script_dir <- function() {
  args <- commandArgs(FALSE)
  hit <- grep("^--file=", args, value = TRUE)
  if (length(hit) > 0) return(dirname(normalizePath(sub("^--file=", "", hit[1]))))
  if (basename(getwd()) == "code") return(normalizePath(getwd()))
  normalizePath(file.path(getwd(), "code"), mustWork = TRUE)
}
source(file.path(script_dir(), "00_paths.R"))

event_table <- function() {
  events <- read_csv_required(path("data", "mtd_events.csv"))
  required <- c("event_id", "original_top_mbsf", "original_base_mbsf", "event_free_depth_cm")
  if (!all(required %in% names(events))) {
    stop("mtd_events.csv must contain: ", paste(required, collapse = ", "))
  }
  if (nrow(events) != 11L) stop("Expected 11 MTD events; found ", nrow(events))
  events
}

posterior_at_depths <- function(set, depths) {
  age_curves_ka <- bacon_posterior_age_curves_ka(set)
  out <- vapply(
    depths,
    function(d) vapply(
      seq_len(nrow(age_curves_ka)),
      function(i) {
        stats::approx(
          x = set$elbows,
          y = age_curves_ka[i, ],
          xout = d,
          rule = 2,
          ties = "ordered"
        )$y
      },
      numeric(1)
    ) * 1000,
    numeric(nrow(age_curves_ka))
  )
  if (is.null(dim(out))) out <- matrix(out, ncol = 1)
  out
}

bacon_posterior_age_curves_ka <- function(set) {
  if (is.null(set$output) || is.null(set$elbows) || is.null(set$K) ||
      ncol(set$output) < set$K + 1L || length(set$elbows) != set$K) {
    stop("The selected Bacon model object lacks the retained elbow-age inputs.")
  }
  output <- as.matrix(set$output)
  storage.mode(output) <- "numeric"
  top_age <- output[, 1]
  accumulation <- output[, 2:(set$K + 1L), drop = FALSE]
  age_curves_ka <- cbind(
    top_age,
    top_age + t(apply(
      accumulation[, -ncol(accumulation), drop = FALSE],
      1,
      cumsum
    )) * set$thick
  ) / 1000
  if (any(apply(age_curves_ka, 1, function(x) any(diff(x) <= 0)))) {
    stop("At least one retained Bacon age-depth curve is not strictly monotonic.")
  }
  age_curves_ka
}

bacon_depth_at_age_ka <- function(age_curve_ka, elbows_cm, target_age_ka) {
  depth <- stats::approx(
    x = age_curve_ka,
    y = elbows_cm,
    xout = target_age_ka,
    rule = 1,
    ties = "ordered"
  )$y
  below_model_top <- target_age_ka < age_curve_ka[1]
  if (any(below_model_top)) {
    top_slope <- (elbows_cm[2] - elbows_cm[1]) /
      (age_curve_ka[2] - age_curve_ka[1])
    depth[below_model_top] <- elbows_cm[1] +
      (target_age_ka[below_model_top] - age_curve_ka[1]) * top_slope
  }
  depth
}

bacon_posterior_interval_rates <- function(set, age_left_ka, age_right_ka) {
  if (!is.numeric(age_left_ka) || !is.numeric(age_right_ka) ||
      length(age_left_ka) != length(age_right_ka) ||
      length(age_left_ka) == 0L || anyNA(age_left_ka) || anyNA(age_right_ka) ||
      any(!is.finite(age_left_ka)) || any(!is.finite(age_right_ka)) ||
      any(age_right_ka <= age_left_ka)) {
    stop("Posterior interval bounds must be finite paired values with right > left.")
  }
  age_curves_ka <- bacon_posterior_age_curves_ka(set)
  target_age_ka <- unique(c(age_left_ka, age_right_ka))
  depth_at_target <- t(apply(
    age_curves_ka,
    1,
    bacon_depth_at_age_ka,
    elbows_cm = set$elbows,
    target_age_ka = target_age_ka
  ))
  rates <- sweep(
    depth_at_target[, match(age_right_ka, target_age_ka), drop = FALSE] -
      depth_at_target[, match(age_left_ka, target_age_ka), drop = FALSE],
    2,
    age_right_ka - age_left_ka,
    "/"
  )
  colnames(rates) <- paste0(age_left_ka, "-", age_right_ka, "_ka")
  rates
}

left_closed_interval_index <- function(age_ka, age_edges_ka) {
  if (!is.numeric(age_edges_ka) || length(age_edges_ka) < 2L ||
      anyNA(age_edges_ka) || any(!is.finite(age_edges_ka)) ||
      is.unsorted(age_edges_ka, strictly = TRUE)) {
    stop("Age edges must be a strictly increasing finite numeric vector.")
  }
  index <- findInterval(age_ka, age_edges_ka)
  index[!is.finite(age_ka) | age_ka < age_edges_ka[1] |
    age_ka >= tail(age_edges_ka, 1)] <- NA_integer_
  as.integer(index)
}

link_same_row_event_rates <- function(rate_by_bin, event_age_draws_ka, age_edges_ka) {
  rate_by_bin <- as.matrix(rate_by_bin)
  event_age_draws_ka <- as.matrix(event_age_draws_ka)
  storage.mode(rate_by_bin) <- "numeric"
  storage.mode(event_age_draws_ka) <- "numeric"
  if (nrow(rate_by_bin) != nrow(event_age_draws_ka)) {
    stop("Rate bins and event ages must contain the same posterior rows.")
  }
  if (ncol(rate_by_bin) != length(age_edges_ka) - 1L) {
    stop("Rate-bin columns must correspond one-to-one with consecutive age edges.")
  }
  linked <- matrix(
    NA_real_,
    nrow = nrow(event_age_draws_ka),
    ncol = ncol(event_age_draws_ka),
    dimnames = dimnames(event_age_draws_ka)
  )
  for (j in seq_len(ncol(event_age_draws_ka))) {
    bin_index <- left_closed_interval_index(event_age_draws_ka[, j], age_edges_ka)
    ok <- !is.na(bin_index)
    linked[ok, j] <- rate_by_bin[cbind(which(ok), bin_index[ok])]
  }
  linked
}

event_draws_and_summary <- function(set) {
  events <- event_table()
  draws <- posterior_at_depths(set, events$event_free_depth_cm)
  colnames(draws) <- events$event_id
  q <- t(apply(
    draws, 2, stats::quantile,
    probs = c(0.025, 0.5, 0.975), na.rm = TRUE
  ))
  summary <- data.frame(
    event_id = events$event_id,
    original_top_mbsf = events$original_top_mbsf,
    original_base_mbsf = events$original_base_mbsf,
    event_free_depth_cm = events$event_free_depth_cm,
    q2.5_ka = q[, 1] / 1000,
    median_ka = q[, 2] / 1000,
    mean_ka = colMeans(draws, na.rm = TRUE) / 1000,
    q97.5_ka = q[, 3] / 1000,
    stringsAsFactors = FALSE
  )
  list(draws_years = draws, summary = summary)
}
