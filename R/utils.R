`%||%` <- function(x, y) if (is.null(x) || length(x) == 0L) y else x

script_path <- function() {
  hit <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (!length(hit)) return(normalizePath(getwd(), winslash = "/", mustWork = TRUE))
  normalizePath(sub("^--file=", "", hit[[1L]]), winslash = "/", mustWork = TRUE)
}

find_repo_root <- function(start = script_path()) {
  p <- if (dir.exists(start)) start else dirname(start)
  p <- normalizePath(p, winslash = "/", mustWork = TRUE)
  repeat {
    marker <- file.path(p, "config", "example_config.R")
    if (file.exists(marker) && dir.exists(file.path(p, "R"))) return(p)
    parent <- dirname(p)
    if (identical(parent, p)) stop("Could not locate repository root from: ", start)
    p <- parent
  }
}

parse_cli_args <- function(args = commandArgs(trailingOnly = TRUE)) {
  out <- list(config = NULL, force = FALSE, only = NULL)
  for (arg in args) {
    if (grepl("^--config=", arg)) out$config <- sub("^--config=", "", arg)
    if (identical(arg, "--force")) out$force <- TRUE
    if (grepl("^--only=", arg)) out$only <- sub("^--only=", "", arg)
  }
  out
}

resolve_repo_path <- function(path, root) {
  if (is.null(path) || !nzchar(path)) return(path)
  if (grepl("^[A-Za-z]:[/\\\\]", path) || grepl("^/", path)) {
    return(normalizePath(path, winslash = "/", mustWork = FALSE))
  }
  normalizePath(file.path(root, path), winslash = "/", mustWork = FALSE)
}

load_config <- function(root = find_repo_root(), config_path = NULL) {
  path <- config_path %||% file.path(root, "config", "local_config.R")
  if (!file.exists(path)) path <- file.path(root, "config", "example_config.R")
  env <- new.env(parent = baseenv())
  sys.source(path, envir = env)
  if (!exists("config", envir = env, inherits = FALSE) || !is.list(env$config)) {
    stop("Config file must define a list named `config`: ", path)
  }
  cfg <- env$config
  path_keys <- c("primary_long", "revision_primary_long", "model3_mi_long",
                 "external_long", "mediation_patient", "equal_lag_long", "patient_flags",
                 "revision_covariates", "output_dir", "figure_dir")
  for (key in intersect(path_keys, names(cfg))) cfg[[key]] <- resolve_repo_path(cfg[[key]], root)
  cfg$repo_root <- root
  cfg$config_file <- normalizePath(path, winslash = "/", mustWork = TRUE)
  cfg
}

require_packages <- function(packages) {
  missing <- packages[!vapply(packages, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))]
  if (length(missing)) {
    stop("Missing R package(s): ", paste(missing, collapse = ", "),
         ". Install them explicitly before running the analysis.")
  }
  invisible(TRUE)
}

ensure_dir <- function(path) {
  if (!dir.exists(path) && !dir.create(path, recursive = TRUE, showWarnings = FALSE)) {
    stop("Could not create output directory: ", path)
  }
  normalizePath(path, winslash = "/", mustWork = TRUE)
}

analysis_output_dir <- function(cfg, analysis_id) {
  ensure_dir(file.path(cfg$output_dir, analysis_id))
}

read_analysis_file <- function(path, sheet = 1L) {
  if (!file.exists(path)) stop("Input file does not exist: ", path)
  ext <- tolower(tools::file_ext(path))
  if (ext == "csv") return(read.csv(path, check.names = FALSE, stringsAsFactors = FALSE))
  if (ext == "tsv") return(read.delim(path, check.names = FALSE, stringsAsFactors = FALSE))
  if (ext %in% c("xlsx", "xls")) {
    require_packages("readxl")
    return(as.data.frame(readxl::read_excel(path, sheet = sheet), check.names = FALSE))
  }
  if (ext == "rds") return(readRDS(path))
  stop("Unsupported input format: .", ext)
}

rename_aliases <- function(data, aliases) {
  for (canonical in names(aliases)) {
    if (canonical %in% names(data)) next
    hit <- intersect(aliases[[canonical]], names(data))
    if (length(hit)) names(data)[match(hit[[1L]], names(data))] <- canonical
  }
  data
}

assert_columns <- function(data, required, label = "data") {
  missing <- setdiff(required, names(data))
  if (length(missing)) stop(label, " is missing required column(s): ", paste(missing, collapse = ", "))
  invisible(TRUE)
}

primary_aliases <- list(
  subject_id = c("subject_id", "analysis_id", "id", "ID(study group)"),
  day = c("day", "start_day"),
  fluid_intake_ml = c("fluid_intake_ml", "fluid_t", "fluid_input_ml"),
  urine_output_ml = c("urine_output_ml", "urine_t", "urine_ml"),
  estimated_fluid_balance_ml = c("estimated_fluid_balance_ml", "fluid_balance_t", "fluid_balance_ml"),
  iap_current = c("iap_current", "iap_t", "iap_start"),
  iap_next = c("iap_next", "iap_end"),
  apache_ii = c("apache_ii", "apache_t", "apache"),
  creatinine = c("creatinine", "cr_t", "cr"),
  map = c("map", "map_t", "map_repaired"),
  sbp = c("sbp", "sbp_t"), dbp = c("dbp", "dbp_t"),
  age = c("age"), sex = c("sex", "sex_code"), etiology = c("etiology", "etiology_code"),
  weight_kg = c("weight_kg", "weight"), heart_rate = c("heart_rate", "pulse_t", "pulse"),
  albumin = c("albumin", "alb_t", "alb"), ph = c("ph", "ph_t"), spo2 = c("spo2", "spo2_t"),
  smoking = c("smoking", "Smoking status"), drinking = c("drinking", "Drinking status"),
  hypertension = c("hypertension", "Hypertension"), diabetes = c("diabetes", "History of Diabetes"),
  hyperlipidemia_history = c("hyperlipidemia_history", "History of hyperlipidemia"), copd = c("copd", "COPD"),
  shock = c("shock", "shock_status", "Shock"), pcd = c("pcd", "pcd_status", "PCD", "PCD_ever"),
  crrt = c("crrt", "crrt_status", "CRRT", "CRRT_ever")
)

prepare_primary <- function(data) {
  d <- rename_aliases(data, primary_aliases)
  required <- c("subject_id", "day", "fluid_intake_ml", "estimated_fluid_balance_ml",
                "iap_current", "iap_next", "apache_ii", "creatinine", "map",
                "age", "sex", "etiology")
  assert_columns(d, required, "primary longitudinal data")
  d$subject_id <- as.character(d$subject_id)
  numeric_vars <- intersect(c("day", "fluid_intake_ml", "urine_output_ml",
                              "estimated_fluid_balance_ml", "iap_current", "iap_next",
                              "apache_ii", "creatinine", "map", "sbp", "dbp", "age",
                              "weight_kg", "heart_rate", "albumin", "ph", "spo2", "shock", "pcd", "crrt"), names(d))
  for (v in numeric_vars) d[[v]] <- suppressWarnings(as.numeric(d[[v]]))
  d$sex <- factor(d$sex)
  d$etiology <- factor(d$etiology)
  d$iap12_next <- as.integer(d$iap_next >= 12)
  d$iap15_next <- as.integer(d$iap_next >= 15)
  d$iap20_next <- as.integer(d$iap_next >= 20)
  d$fluid_intake_l <- d$fluid_intake_ml / 1000
  d$estimated_fluid_balance_l <- d$estimated_fluid_balance_ml / 1000
  d[order(d$subject_id, d$day), , drop = FALSE]
}

prepare_external <- function(data) {
  d <- rename_aliases(data, primary_aliases)
  if (!"subject_id" %in% names(d) && "id" %in% names(d)) d$subject_id <- d$id
  required <- c("subject_id", "day", "fluid_intake_ml", "iap_current", "iap_next",
                "map", "age", "sex", "etiology")
  assert_columns(d, required, "external longitudinal data")
  d$subject_id <- as.character(d$subject_id)
  numeric_vars <- c("day", "fluid_intake_ml", "iap_current", "iap_next", "map", "age")
  for (v in numeric_vars) d[[v]] <- suppressWarnings(as.numeric(d[[v]]))
  d$sex <- factor(d$sex)
  d$etiology <- factor(d$etiology)
  if ("IAH15_next" %in% names(d)) {
    d$iap15_next <- as.integer(d$IAH15_next)
  } else {
    d$iap15_next <- as.integer(d$iap_next >= 15)
  }
  if ("IAH20_next" %in% names(d)) {
    d$iap20_next <- as.integer(d$IAH20_next)
  } else {
    d$iap20_next <- as.integer(d$iap_next >= 20)
  }
  d$iap12_next <- as.integer(d$iap_next >= 12)
  d$fluid_intake_l <- d$fluid_intake_ml / 1000
  d[order(d$subject_id, d$day), , drop = FALSE]
}

merge_optional_fields <- function(data, path, fields, label = "supplemental data") {
  if (is.null(path) || !nzchar(path) || !file.exists(path)) return(data)
  supplemental <- rename_aliases(read_analysis_file(path), primary_aliases)
  assert_columns(supplemental, "subject_id", label)
  supplemental$subject_id <- as.character(supplemental$subject_id)
  keys <- "subject_id"
  if ("day" %in% names(data) && "day" %in% names(supplemental)) {
    supplemental$day <- suppressWarnings(as.numeric(supplemental$day))
    keys <- c(keys, "day")
  }
  available <- intersect(fields, names(supplemental))
  if (!length(available)) return(data)
  supplemental <- supplemental[c(keys, available)]
  if (anyDuplicated(supplemental[keys])) stop(label, " has duplicate rows for key: ", paste(keys, collapse = ", "))
  names(supplemental)[match(available, names(supplemental))] <- paste0(available, "__supplemental")
  merged <- merge(data, supplemental, by = keys, all.x = TRUE, sort = FALSE)
  for (field in available) {
    supplied <- paste0(field, "__supplemental")
    if (!field %in% names(merged)) merged[[field]] <- NA
    use <- !is.na(merged[[supplied]])
    merged[[field]][use] <- merged[[supplied]][use]
    merged[[supplied]] <- NULL
  }
  merged
}
validate_panel <- function(data, expected_patients = NULL, expected_windows = NULL,
                           expected_days = NULL, strict = TRUE) {
  problems <- character()
  n_patients <- length(unique(data$subject_id))
  if (!is.null(expected_patients) && n_patients != expected_patients) {
    problems <- c(problems, sprintf("patients=%s (expected %s)", n_patients, expected_patients))
  }
  if (!is.null(expected_windows) && nrow(data) != expected_windows) {
    problems <- c(problems, sprintf("windows=%s (expected %s)", nrow(data), expected_windows))
  }
  if (anyDuplicated(data[c("subject_id", "day")])) problems <- c(problems, "duplicate subject-day rows")
  if (!is.null(expected_days) && !identical(sort(unique(data$day)), sort(expected_days))) {
    problems <- c(problems, "unexpected day values")
  }
  if (length(problems) && strict) stop("Panel validation failed: ", paste(problems, collapse = "; "))
  data.frame(patients = n_patients, windows = nrow(data), problems = paste(problems, collapse = "; "))
}

write_csv_atomic <- function(data, path) {
  ensure_dir(dirname(path))
  tmp <- tempfile(pattern = "csv_", tmpdir = dirname(path), fileext = ".tmp")
  on.exit(unlink(tmp), add = TRUE)
  write.csv(data, tmp, row.names = FALSE, na = "")
  if (file.exists(path)) unlink(path)
  if (!file.rename(tmp, path)) stop("Could not finalize output: ", path)
  invisible(path)
}

write_run_metadata <- function(cfg, input_paths, output_dir, extra = list()) {
  paths <- input_paths[file.exists(input_paths)]
  info <- data.frame(
    item = c("timestamp_utc", "config_file", paste0("input_", seq_along(paths)),
             paste0("input_md5_", seq_along(paths)), names(extra)),
    value = c(format(Sys.time(), tz = "UTC", usetz = TRUE), cfg$config_file,
              normalizePath(paths, winslash = "/"), unname(tools::md5sum(paths)),
              unlist(extra, use.names = FALSE)),
    stringsAsFactors = FALSE
  )
  write_csv_atomic(info, file.path(output_dir, "run_metadata.csv"))
  capture.output(sessionInfo(), file = file.path(output_dir, "session_info.txt"))
  invisible(info)
}

skip_unless_enabled <- function(flag, label) {
  if (!isTRUE(flag)) {
    message("SKIP: ", label, " is disabled in the selected configuration.")
    return(TRUE)
  }
  FALSE
}
