file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
root <- normalizePath(file.path(dirname(sub("^--file=", "", file_arg[[1L]])), ".."),
                      winslash = "/", mustWork = TRUE)

fail <- function(...) stop(..., call. = FALSE)

required <- c(
  "README.md", "config/example_config.R", "data/README.md",
  "docs/analysis_crosswalk.md", "docs/methodological_notes.md",
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

deprecated_root_scripts <- c(
  "aipw_analysis.R", "aipw_analysis_extended.R", "causal_mediation_analysis.R",
  "main_association_analysis_per1000ml.R", "msm_step1_iptw_construction.R",
  "msm_step3_outcome_models.R", "rcs_revised_threshold_analysis.R"
)
remaining_deprecated <- deprecated_root_scripts[file.exists(file.path(root, deprecated_root_scripts))]
if (length(remaining_deprecated)) {
  fail("Deprecated compatibility entry point(s) remain: ", paste(remaining_deprecated, collapse = ", "))
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

analysis_text <- paste(vapply(analysis_files, function(path) {
  paste(readLines(path, warn = FALSE, encoding = "UTF-8"), collapse = "\n")
}, character(1)), collapse = "\n")
if (grepl("estimated_fluid_balance_(ml|l)", analysis_text, perl = TRUE)) {
  fail("Analysis scripts must use canonical fluid_balance_ml/fluid_balance_l names")
}
utils_text <- paste(readLines(file.path(root, "R", "utils.R"), warn = FALSE, encoding = "UTF-8"), collapse = "\n")
if (!grepl('fluid_balance_ml = c\\("fluid_balance_ml", "estimated_fluid_balance_ml", "fluid_balance_t"\\)', utils_text)) {
  fail("Historical fluid-balance aliases are not preserved in the input layer")
}

data_files <- list.files(file.path(root, "data"), recursive = TRUE, full.names = TRUE)
data_files <- data_files[!dir.exists(data_files)]
allowed_data <- normalizePath(file.path(root, "data", "README.md"), winslash = "/", mustWork = TRUE)
unexpected_data <- setdiff(normalizePath(data_files, winslash = "/", mustWork = TRUE), allowed_data)
if (length(unexpected_data)) fail("Unexpected file under data/: ", paste(basename(unexpected_data), collapse = ", "))

cat("PASS: R syntax, canonical naming, portable paths, repository structure, and data-boundary checks\n")
