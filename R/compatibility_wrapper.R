run_compatibility <- function(root, target, force = FALSE) {
  args <- commandArgs(trailingOnly = TRUE)
  if (force && !"--force" %in% args) args <- c("--force", args)
  target_path <- file.path(root, target)
  if (!file.exists(target_path)) stop("Compatibility target does not exist: ", target_path)
  message("Compatibility entry point; dispatching to ", target)
  quoted_args <- c(shQuote(target_path), vapply(args, shQuote, character(1)))
  status <- system2(file.path(R.home("bin"), "Rscript"), args = quoted_args)
  quit(save = "no", status = status)
}
