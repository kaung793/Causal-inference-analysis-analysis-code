file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
root <- normalizePath(dirname(sub("^--file=", "", file_arg[[1L]])), winslash = "/")
source(file.path(root, "R", "compatibility_wrapper.R"))
run_compatibility(root, "analysis/04_rcs_dose_response.R", force = FALSE)
