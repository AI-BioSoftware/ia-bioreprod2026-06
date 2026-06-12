#' Access the bundled supplementary Excel file
#'
#' @return A path to the bundled supplementary Excel file.
#' @export
default_excel_path <- function() {
  system.file("extdata", "pgen.1006453.s002-3.xlsx", package = "fig2aheatmap")
}

validate_required_columns <- function(data, required_columns) {
  missing_columns <- setdiff(required_columns, names(data))
  if (length(missing_columns) > 0) {
    stop(
      "Missing required columns: ",
      paste(missing_columns, collapse = ", "),
      call. = FALSE
    )
  }
  invisible(TRUE)
}

extract_time_columns <- function(data) {
  time_cols <- names(data)[stringr::str_detect(names(data), "^X?\\d+\\.?\\d*$")]
  if (length(time_cols) == 0) {
    stop("No expression time columns were detected.", call. = FALSE)
  }
  time_cols
}

coerce_numeric_matrix <- function(data, time_cols) {
  expr_tbl <- dplyr::select(data, dplyr::all_of(time_cols))
  expr_tbl[] <- lapply(expr_tbl, as.numeric)
  as.matrix(expr_tbl)
}

z_score_rows <- function(mat) {
  z <- t(scale(t(mat)))
  z[is.na(z)] <- 0
  z
}

#' Prepare the data needed to reproduce Figure 2A
#'
#' This function reads the supplementary Excel file, validates the expected
#' columns, filters genes with `LS_cutoff == "Yes"`, orders them by
#' `Figure2A_order_peaktime`, computes gene-wise z-scores, and returns the
#' processed objects needed for export or plotting.
#'
#' @param xlsx_file Path to the supplementary Excel file.
#' @param sheet_name Worksheet name. Defaults to `"Sheet1"`.
#' @param skip Number of leading rows to skip before the real header.
#' @param top_n Number of genes to retain after ordering.
#' @param reverse_rows Whether to reverse row order for plotting.
#'
#' @return A named list with filtered data, z-score matrix, time values and gene IDs.
#' @export
#'
#' @examples
#' prepared <- prepare_figure2a_data(default_excel_path())
prepare_figure2a_data <- function(
  xlsx_file,
  sheet_name = "Sheet1",
  skip = 2,
  top_n = 1246,
  reverse_rows = TRUE
) {
  if (!file.exists(xlsx_file)) {
    stop("Input file does not exist: ", xlsx_file, call. = FALSE)
  }

  raw <- readxl::read_excel(xlsx_file, sheet = sheet_name, skip = skip)
  names(raw) <- make.names(names(raw), unique = TRUE)

  required_columns <- c("gene_ID", "LS_cutoff", "Figure2A_order_peaktime")
  validate_required_columns(raw, required_columns)
  time_cols <- extract_time_columns(raw)

  df <- raw %>%
    dplyr::mutate(
      LS_cutoff = as.character(.data$LS_cutoff),
      Figure2A_order_peaktime = suppressWarnings(as.numeric(.data$Figure2A_order_peaktime))
    ) %>%
    dplyr::filter(.data$LS_cutoff == "Yes") %>%
    dplyr::filter(!is.na(.data$Figure2A_order_peaktime)) %>%
    dplyr::arrange(.data$Figure2A_order_peaktime) %>%
    dplyr::slice_head(n = top_n)

  if (nrow(df) == 0) {
    stop("Filtering produced an empty data set.", call. = FALSE)
  }

  expr <- coerce_numeric_matrix(df, time_cols)
  rownames(expr) <- df$gene_ID

  expr_z <- z_score_rows(expr)
  time_vals <- as.numeric(stringr::str_replace(time_cols, "^X", ""))
  col_order <- order(time_vals)
  expr_z <- expr_z[, col_order, drop = FALSE]
  time_vals <- time_vals[col_order]

  if (reverse_rows) {
    expr_z <- expr_z[nrow(expr_z):1, , drop = FALSE]
  }

  list(
    data = tibble::as_tibble(df),
    expr_z = expr_z,
    time_vals = time_vals,
    gene_ids = rownames(expr_z)
  )
}

#' Export the processed z-score matrix to CSV
#'
#' @param prepared Output of [prepare_figure2a_data()].
#' @param output_csv Path to the CSV file to create.
#'
#' @return The output path, invisibly.
#' @export
export_figure2a_matrix <- function(prepared, output_csv) {
  stopifnot(is.list(prepared), !is.null(prepared$expr_z))
  out <- data.frame(
    gene_ID = rownames(prepared$expr_z),
    prepared$expr_z,
    check.names = FALSE
  )
  utils::write.csv(out, output_csv, row.names = FALSE)
  invisible(output_csv)
}

#' Draw a Figure 2A-style heatmap and save it to PDF
#'
#' @param prepared Output of [prepare_figure2a_data()].
#' @param output_pdf Path to the PDF file to create.
#' @param tick_times Numeric vector of tick labels in minutes.
#' @param width Width of the output PDF in inches.
#' @param height Height of the output PDF in inches.
#'
#' @return The output path, invisibly.
#' @export
plot_figure2a_heatmap <- function(
  prepared,
  output_pdf,
  tick_times = c(0, 50, 100, 150, 200),
  width = 9,
  height = 12
) {
  stopifnot(is.list(prepared), !is.null(prepared$expr_z), !is.null(prepared$time_vals))

  keep <- tick_times %in% prepared$time_vals
  tick_labels <- tick_times[keep]
  tick_at <- match(tick_labels, prepared$time_vals)

  bottom_ha <- ComplexHeatmap::HeatmapAnnotation(
    time = ComplexHeatmap::anno_empty(border = FALSE, height = grid::unit(8, "mm")),
    show_annotation_name = FALSE
  )

  ht <- ComplexHeatmap::Heatmap(
    prepared$expr_z,
    name = "z-score",
    col = circlize::colorRamp2(c(-1.5, 0, 1.5), c("cyan", "black", "yellow")),
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    show_row_names = FALSE,
    show_column_names = FALSE,
    row_title = "Periodic genes",
    column_title = "Time (minutes)",
    use_raster = FALSE,
    bottom_annotation = bottom_ha
  )

  grDevices::pdf(output_pdf, width = width, height = height)
  on.exit(grDevices::dev.off(), add = TRUE)
  ComplexHeatmap::draw(ht)
  ComplexHeatmap::decorate_annotation("time", {
    x <- (tick_at - 0.5) / ncol(prepared$expr_z)
    grid::grid.segments(
      x0 = grid::unit(x, "npc"),
      y0 = grid::unit(1, "npc"),
      x1 = grid::unit(x, "npc"),
      y1 = grid::unit(0.72, "npc"),
      gp = grid::gpar(lwd = 1)
    )
    grid::grid.text(
      label = tick_labels,
      x = grid::unit(x, "npc"),
      y = grid::unit(0.25, "npc"),
      gp = grid::gpar(fontsize = 10)
    )
  })

  invisible(output_pdf)
}
