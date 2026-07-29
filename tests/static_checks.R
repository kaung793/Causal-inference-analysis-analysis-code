file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
root <- normalizePath(file.path(dirname(sub("^--file=", "", file_arg[[1L]])), ".."),
                      winslash = "/", mustWork = TRUE)

fail <- function(...) stop(..., call. = FALSE)

required <- c(
  "README.md", "config/example_config.R", "data/README.md",
  "docs/analysis_crosswalk.md", "docs/reconciliation_notes.md",
  "docs/software_environment.md", "run_all.R"
)
missing <- required[!file.exists(file.path(root, required))]
if (length(missing)) fail("Missing required repository file(s): ", paste(missing, collapse = ", "))

analysis_files <- sort(list.files(file.path(root, "analysis"), pattern = "[.]R$", full.names = TRUE))
expected_scripts <- sprintf("%02d", 0:16)
observed_scripts <- substr(basename(analysis_files), 1, 2)
if (!identical(observed_scripts, expected_scripts)) {
  fail("Expected one ordered analysis script for every ID 00-16")
}

root_r <- list.files(root, pattern = "[.]R$", full.names = TRUE)
r_files <- c(
  list.files(file.path(root, "R"), pattern = "[.]R$", full.names = TRUE),
  analysis_files, root_r, file.path(root, "tests", "static_checks.R")
)
for (path in unique(r_files)) {
  tryCatch(parse(file = path), error = function(e) {
    fail("R parse failure in ", basename(path), ": ", conditionMessage(e))
  })
}

scan_files <- c(
  list.files(file.path(root, "analysis"), pattern = "[.]R$", full.names = TRUE),
  root_r, file.path(root, "config", "example_config.R")
)
forbidden <- c(
  "setwd[[:space:]]*[(]",
  "install[.]packages[[:space:]]*[(]",
  "PROJECT_ROOT_PLACEHOLDER",
  "[A-Za-z]:/",
  "/Users/",
  "/home/"
)
for (path in scan_files) {
  text <- paste(readLines(path, warn = FALSE, encoding = "UTF-8"), collapse = "\n")
  hits <- forbidden[vapply(forbidden, grepl, logical(1), x = text, perl = TRUE)]
  if (length(hits)) fail("Nonportable pattern in ", basename(path), ": ", paste(hits, collapse = ", "))
}

expected_dir <- file.path(root, "results", "expected")
expected_files <- list.files(expected_dir, pattern = "[.]csv$", full.names = TRUE)
if (length(expected_files) < 6L) fail("Expected aggregate reconciliation CSVs are missing")
for (path in expected_files) {
  x <- tryCatch(read.csv(path, check.names = FALSE), error = function(e) {
    fail("Could not parse aggregate checkpoint ", basename(path), ": ", conditionMessage(e))
  })
  if (!nrow(x)) fail("Aggregate checkpoint is empty: ", basename(path))
  if ("subject_id" %in% names(x)) fail("Participant-level identifier found in expected result: ", basename(path))
}

data_files <- list.files(file.path(root, "data"), recursive = TRUE, full.names = TRUE)
data_files <- data_files[!dir.exists(data_files)]
allowed_data <- normalizePath(file.path(root, "data", "README.md"), winslash = "/", mustWork = TRUE)
unexpected_data <- setdiff(normalizePath(data_files, winslash = "/", mustWork = TRUE), allowed_data)
if (length(unexpected_data)) fail("Unexpected file under data/: ", paste(basename(unexpected_data), collapse = ", "))

cat("PASS: R syntax, portable-path, repository-structure, and aggregate-data checks\n")
