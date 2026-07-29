file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
root <- normalizePath(file.path(dirname(sub("^--file=", "", file_arg[[1L]])), ".."), winslash = "/")
source(file.path(root, "R", "utils.R"))
cli <- parse_cli_args(); cfg <- load_config(root, cli$config)
if (skip_unless_enabled(isTRUE(cfg$run_figures) || cli$force, "figure generation")) quit(save = "no", status = 0)
require_packages("ggplot2")
out <- ensure_dir(cfg$figure_dir)
produced <- character()
save_plot <- function(plot, stem, width, height) {
  png <- file.path(out, paste0(stem, ".png")); pdf <- file.path(out, paste0(stem, ".pdf"))
  ggplot2::ggsave(png, plot, width = width, height = height, dpi = 300, bg = "white")
  ggplot2::ggsave(pdf, plot, width = width, height = height, device = cairo_pdf, bg = "white")
  produced <<- c(produced, png, pdf)
}

primary_path <- file.path(cfg$output_dir, "02_primary_mixed_models", "primary_mixed_model_results.csv")
if (file.exists(primary_path)) {
  x <- read.csv(primary_path, check.names = FALSE); x <- x[x$model == "Model 3", ]
  x$label <- paste(x$exposure, x$outcome, sep = " | ")
  p <- ggplot2::ggplot(x, ggplot2::aes(estimate, reorder(label, estimate), colour = exposure)) +
    ggplot2::geom_errorbar(ggplot2::aes(xmin = ci_lower, xmax = ci_upper), orientation = "y", width = .18) +
    ggplot2::geom_point(size = 2.5) + ggplot2::facet_wrap(~effect_measure, scales = "free_x", ncol = 1) +
    ggplot2::theme_bw(base_size = 10) + ggplot2::labs(x = NULL, y = NULL, colour = NULL)
  save_plot(p, "primary_model3_forest", 7, 6)
}

rcs_path <- file.path(cfg$output_dir, "04_rcs_dose_response", "rcs_prediction_curves.csv")
if (file.exists(rcs_path)) {
  x <- read.csv(rcs_path, check.names = FALSE); x$exposure_label <- ifelse(x$exposure == "fluid_intake_ml", "Fluid intake", "Fluid balance")
  p <- ggplot2::ggplot(x, ggplot2::aes(exposure_ml, effect, colour = exposure_label, fill = exposure_label)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = ci_lower, ymax = ci_upper), alpha = .15, colour = NA) +
    ggplot2::geom_line(linewidth = .8) + ggplot2::facet_wrap(~outcome, scales = "free_y") +
    ggplot2::theme_bw(base_size = 10) + ggplot2::labs(x = "Exposure (mL)", y = "Adjusted effect relative to the median", colour = NULL, fill = NULL)
  save_plot(p, "rcs_dose_response", 9, 6)
}

time_path <- file.path(cfg$output_dir, "07_time_window_heterogeneity", "time_window_results.csv")
if (file.exists(time_path)) {
  x <- read.csv(time_path, check.names = FALSE); x$window_label <- factor(x$window_label, levels = unique(x$window_label))
  p <- ggplot2::ggplot(x, ggplot2::aes(estimate, window_label)) +
    ggplot2::geom_errorbar(ggplot2::aes(xmin = ci_lower, xmax = ci_upper), orientation = "y", width = .18) + ggplot2::geom_point(size = 2.5) +
    ggplot2::facet_wrap(~outcome, scales = "free_x") + ggplot2::theme_bw(base_size = 10) + ggplot2::labs(x = NULL, y = NULL)
  save_plot(p, "time_window_heterogeneity", 8, 4.5)
}

med_path <- file.path(cfg$output_dir, "09_mediation_mi_rubin", "mediation_mi_rubin_pooled_results.csv")
if (file.exists(med_path)) {
  x <- read.csv(med_path, check.names = FALSE); x <- x[x$effect != "Proportion mediated", ]
  p <- ggplot2::ggplot(x, ggplot2::aes(estimate, effect)) + ggplot2::geom_vline(xintercept = 0, linetype = 2, colour = "grey60") +
    ggplot2::geom_errorbar(ggplot2::aes(xmin = ci_lower, xmax = ci_upper), orientation = "y", width = .18) + ggplot2::geom_point(size = 2.5) +
    ggplot2::facet_wrap(~exposure, scales = "free_x", ncol = 1) + ggplot2::theme_bw(base_size = 10) +
    ggplot2::labs(x = "Absolute change in in-hospital mortality risk", y = NULL)
  save_plot(p, "mediation_mi_rubin", 7, 5.5)
}

writeLines(normalizePath(produced, winslash = "/", mustWork = FALSE), file.path(out, "figure_manifest.txt"))
message("Generated ", length(produced), " figure files")
