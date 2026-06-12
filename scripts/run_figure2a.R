#!/usr/bin/env Rscript

`%||%` <- function(x, y) if (is.null(x)) y else x

args <- commandArgs(trailingOnly = TRUE)
script_path <- tryCatch(normalizePath(sys.frame(1)$ofile), error = function(e) NULL)
root <- normalizePath(file.path(dirname(script_path %||% "scripts/run_figure2a.R"), ".."), mustWork = FALSE)
source(file.path(root, "R", "figure2a.R"))

input_xlsx <- if (length(args) >= 1) {
  args[[1]]
} else {
  file.path(root, "inst", "extdata", "pgen.1006453.s002-3.xlsx")
}

out_dir <- if (length(args) >= 2) args[[2]] else file.path(root, "outputs")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

prepared <- prepare_figure2a_data(input_xlsx)
export_figure2a_matrix(prepared, file.path(out_dir, "figure_2A_matrix.csv"))
plot_figure2a_heatmap(prepared, file.path(out_dir, "figure_2A_recreated.pdf"))

message("Outputs written to: ", out_dir)
