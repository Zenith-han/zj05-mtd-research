package_root <- function() {
  configured_root <- getOption("zj05.public_root", NULL)
  if (!is.null(configured_root)) {
    if (length(configured_root) != 1L || !nzchar(configured_root)) {
      stop("Option zj05.public_root must contain one non-empty path.")
    }
    return(normalizePath(configured_root, mustWork = TRUE))
  }
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- "--file="
  hit <- grep(paste0("^", file_arg), args, value = TRUE)
  if (length(hit) > 0) {
    script_path <- normalizePath(sub(file_arg, "", hit[1]), mustWork = TRUE)
    return(normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE))
  }
  wd <- normalizePath(getwd(), mustWork = TRUE)
  if (basename(wd) == "code") {
    return(normalizePath(file.path(wd, ".."), mustWork = TRUE))
  }
  wd
}

root_dir <- package_root()
path <- function(...) file.path(root_dir, ...)

results_root <- function() {
  configured_root <- getOption("zj05.results_root", NULL)
  if (is.null(configured_root)) {
    return(normalizePath(path("results"), mustWork = FALSE))
  }
  if (length(configured_root) != 1L || is.na(configured_root) ||
      !is.character(configured_root) || !nzchar(configured_root)) {
    stop("Option zj05.results_root must contain one non-empty path.")
  }
  normalizePath(configured_root, mustWork = FALSE)
}

result_path <- function(...) file.path(results_root(), ...)

repo_rel <- function(file) {
  normalized <- normalizePath(file, winslash = "/", mustWork = FALSE)
  root <- normalizePath(root_dir, winslash = "/", mustWork = TRUE)
  prefix <- paste0(root, "/")
  if (normalized == root) return(".")
  if (startsWith(normalized, prefix)) return(substr(normalized, nchar(prefix) + 1L, nchar(normalized)))
  normalized
}

local_r_lib <- file.path(root_dir, "R_libs")
if (dir.exists(local_r_lib)) {
  .libPaths(c(normalizePath(local_r_lib), .libPaths()))
}

cli_value <- function(name, default = NULL, args = commandArgs(trailingOnly = TRUE)) {
  prefix <- paste0("--", name, "=")
  hit <- args[startsWith(args, prefix)]
  if (length(hit) == 0) return(default)
  sub(prefix, "", hit[1], fixed = TRUE)
}

cli_flag <- function(name, args = commandArgs(trailingOnly = TRUE)) {
  paste0("--", name) %in% args
}

as_cli_numeric <- function(name, default = NULL) {
  value <- cli_value(name, default)
  if (is.null(value)) return(NULL)
  value <- as.numeric(value)
  if (!is.finite(value)) stop("Invalid numeric --", name)
  value
}

ensure_dir <- function(...) {
  d <- path(...)
  if (!dir.exists(d)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
  d
}

ensure_result_dir <- function(...) {
  d <- result_path(...)
  if (!dir.exists(d)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
  d
}

read_csv_required <- function(file, ...) {
  if (!file.exists(file)) stop("Required file not found: ", file)
  utils::read.csv(file, stringsAsFactors = FALSE, check.names = FALSE, ...)
}

write_csv <- function(x, rel, ...) {
  file <- path(rel)
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(x, file, row.names = FALSE, ...)
  message("Wrote ", file)
  file
}

read_result_csv_required <- function(rel, ...) {
  read_csv_required(result_path(rel), ...)
}

write_result_csv <- function(x, rel, ...) {
  file <- result_path(rel)
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(x, file, row.names = FALSE, ...)
  message("Wrote ", file)
  file
}

create_result_staging_root <- function(target_root) {
  target_root <- normalizePath(target_root, winslash = "/", mustWork = TRUE)
  staging_root <- tempfile(
    pattern = paste0(".", basename(target_root), "-staging-"),
    tmpdir = dirname(target_root)
  )
  if (!dir.create(staging_root, recursive = FALSE, showWarnings = FALSE)) {
    stop("Could not create rebuilt-result staging directory: ", staging_root)
  }
  normalizePath(staging_root, winslash = "/", mustWork = TRUE)
}

assert_result_target_containment <- function(target_root, relative_paths) {
  target_root <- normalizePath(target_root, winslash = "/", mustWork = TRUE)
  root_prefix <- paste0(target_root, "/")
  for (relative_path in relative_paths) {
    existing_parent <- dirname(file.path(target_root, relative_path))
    while (!file.exists(existing_parent)) {
      parent <- dirname(existing_parent)
      if (identical(parent, existing_parent)) {
        stop("Could not resolve result target parent: ", relative_path)
      }
      existing_parent <- parent
    }
    resolved_parent <- normalizePath(
      existing_parent,
      winslash = "/",
      mustWork = TRUE
    )
    if (!identical(resolved_parent, target_root) &&
        !startsWith(resolved_parent, root_prefix)) {
      stop(
        "Result path resolves outside the selected result root: ",
        relative_path
      )
    }
  }
  invisible(TRUE)
}

rollback_staged_result_files <- function(
  target_files,
  backup_files,
  installed,
  backed_up,
  relative_paths
) {
  rollback_failures <- character()
  if (length(installed) > 0L) {
    removal_failures <- vapply(installed, function(i) {
      if (!file.exists(target_files[i])) return(FALSE)
      suppressWarnings(file.remove(target_files[i]))
      file.exists(target_files[i])
    }, logical(1))
    if (any(removal_failures)) {
      rollback_failures <- c(
        rollback_failures,
        relative_paths[installed[removal_failures]]
      )
    }
  }
  if (length(backed_up) > 0L) {
    restore_failures <- vapply(backed_up, function(i) {
      dir.create(
        dirname(target_files[i]),
        recursive = TRUE,
        showWarnings = FALSE
      )
      !suppressWarnings(file.rename(backup_files[i], target_files[i]))
    }, logical(1))
    if (any(restore_failures)) {
      rollback_failures <- c(
        rollback_failures,
        relative_paths[backed_up[restore_failures]]
      )
    }
  }
  list(
    complete = length(rollback_failures) == 0L,
    failures = unique(rollback_failures)
  )
}

commit_staged_result_files <- function(staging_root, target_root, relative_paths) {
  staging_root <- normalizePath(staging_root, winslash = "/", mustWork = TRUE)
  target_root <- normalizePath(target_root, winslash = "/", mustWork = TRUE)
  relative_paths <- unique(as.character(relative_paths))
  portable_paths <- gsub("\\\\", "/", relative_paths)
  invalid_paths <- !nzchar(relative_paths) |
    grepl("^([A-Za-z]:)?/", portable_paths) |
    vapply(
      strsplit(portable_paths, "/", fixed = TRUE),
      function(parts) any(parts %in% c("", ".", "..")),
      logical(1)
    )
  if (length(relative_paths) == 0L || any(invalid_paths)) {
    stop("Result commit requires non-empty relative file paths.")
  }

  staged_files <- file.path(staging_root, relative_paths)
  missing <- relative_paths[!file.exists(staged_files)]
  if (length(missing) > 0L) {
    stop(
      "Rebuilt analysis did not produce the complete result contract: ",
      paste(missing, collapse = ", ")
    )
  }
  staged_info <- file.info(staged_files)
  if (any(staged_info$isdir)) {
    stop("Rebuilt result contract contains a directory instead of a file.")
  }
  assert_result_target_containment(target_root, relative_paths)

  backup_root <- tempfile(
    pattern = paste0(".", basename(target_root), "-backup-"),
    tmpdir = dirname(target_root)
  )
  if (!dir.create(backup_root, recursive = FALSE, showWarnings = FALSE)) {
    stop("Could not create rebuilt-result backup directory: ", backup_root)
  }

  target_files <- file.path(target_root, relative_paths)
  backup_files <- file.path(backup_root, relative_paths)
  backed_up <- integer()
  installed <- integer()
  committed <- FALSE
  on.exit({
    rollback <- if (committed) {
      list(complete = TRUE, failures = character())
    } else {
      rollback_staged_result_files(
        target_files,
        backup_files,
        installed,
        backed_up,
        relative_paths
      )
    }
    if (rollback$complete) {
      unlink(backup_root, recursive = TRUE, force = TRUE)
    } else {
      warning(
        "Could not fully restore rebuilt results after a failed commit: ",
        paste(rollback$failures, collapse = ", "),
        ". Backup files were preserved at: ", backup_root
      )
    }
  }, add = TRUE)

  for (i in which(file.exists(target_files))) {
    dir.create(dirname(backup_files[i]), recursive = TRUE, showWarnings = FALSE)
    if (!suppressWarnings(file.rename(target_files[i], backup_files[i]))) {
      stop("Could not back up rebuilt result before commit: ", relative_paths[i])
    }
    backed_up <- c(backed_up, i)
  }
  for (i in seq_along(relative_paths)) {
    dir.create(dirname(target_files[i]), recursive = TRUE, showWarnings = FALSE)
    if (!suppressWarnings(file.rename(staged_files[i], target_files[i]))) {
      stop("Could not commit rebuilt result: ", relative_paths[i])
    }
    installed <- c(installed, i)
  }

  committed <- TRUE
  unlink(staging_root, recursive = TRUE, force = TRUE)
  invisible(target_files)
}

append_csv <- function(x, rel) {
  file <- path(rel)
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  exists <- file.exists(file)
  utils::write.table(
    x, file, sep = ",", row.names = FALSE,
    col.names = !exists, append = exists, quote = TRUE
  )
  invisible(file)
}
