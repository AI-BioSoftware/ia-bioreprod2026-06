if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("ComplexHeatmap")
install.packages("readxl")

library(readxl)
library(dplyr)
library(tidyr)
library(ComplexHeatmap)
library(circlize)
library(stringr)
library(grid)

xlsx_file <- "pgen.1006453.s002.xlsx"
sheet_name <- "Sheet1"

raw <- read_excel(
  xlsx_file,
  sheet = sheet_name,
  skip = 2
)

names(raw) <- make.names(names(raw), unique = TRUE)

# Colonnes attendues après make.names()
# gene_ID, normalized_per_rank, LS_cutoff, combined_per_rank, ...
# Figure2A_order_peaktime, Sparklines, X0.0, X5.0, ...

gene_col <- "gene_ID"
ls_cutoff_col <- "LS_cutoff"
order_col <- "Figure2A_order_peaktime"

# Colonnes d'expression = colonnes numériques de temps
time_cols <- names(raw)[str_detect(names(raw), "^X?\\d+\\.?\\d*$")]

# Nettoyage et sélection
df <- raw %>%
  mutate(
    LS_cutoff = as.character(.data[[ls_cutoff_col]]),
    Figure2A_order_peaktime = suppressWarnings(as.numeric(.data[[order_col]])),
    combined_per_rank = suppressWarnings(as.numeric(combined_per_rank))
  ) %>%
  filter(LS_cutoff == "Yes") %>%
  filter(!is.na(Figure2A_order_peaktime)) %>%
  arrange(Figure2A_order_peaktime) %>%
  slice_head(n = 1246) %>%
  arrange(desc(Figure2A_order_peaktime))


# Si vous voulez reproduire strictement l'étape "top 1600"
#avant le filtre Lomb-Scargle, décommentez :
#df <- raw %>%
#  mutate(
#    LS_cutoff = as.character(.data[[ls_cutoff_col]]),
 #   Figure2A_order_peaktime = suppressWarnings(as.numeric(.data[[order_col]])),
#    combined_per_rank = suppressWarnings(as.numeric(combined_per_rank))
#  ) %>%
#  filter(!is.na(combined_per_rank)) %>%
#  arrange(combined_per_rank) %>%
#  slice_head(n = 1600) %>%
#  filter(LS_cutoff == "Yes") %>%
#  arrange(Figure2A_order_peaktime)

#expr <- df %>%
#  select(all_of(time_cols)) %>%
#  mutate(across(everything(), as.numeric)) %>%
#  as.matrix()

rownames(expr) <- df[[gene_col]]

# Z-score par gène
expr_z <- t(scale(t(expr)))
expr_z[is.na(expr_z)] <- 0

# Ordonner les colonnes par temps
time_vals <- as.numeric(str_replace(time_cols, "^X", ""))
col_order <- order(time_vals)
expr_z <- expr_z[, col_order, drop = FALSE]
time_vals <- time_vals[col_order]

# CSV de sortie
out_csv <- "figure_2A_matrix.csv"
write.csv(
  data.frame(gene_ID = rownames(expr_z), expr_z, check.names = FALSE),
  out_csv,
  row.names = FALSE
)


# Heatmap
# inversion matrice 
expr_z <- expr_z[nrow(expr_z):1, , drop = FALSE]


# positions des temps à afficher
tick_times <- c(0, 50, 100, 150, 200)

# indices des colonnes correspondantes dans la matrice
tick_at <- match(tick_times, time_vals)

# annotation d’axe en bas
bottom_ha <- HeatmapAnnotation(
  time = anno_empty(
    border = FALSE,
    height = unit(8, "mm")
  ),
  show_annotation_name = FALSE
)

ht <- Heatmap(
  expr_z,
  name = "z-score",
  col = colorRamp2(c(-1.5, 0, 1.5), c("cyan", "black", "yellow")),
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_names = FALSE,
  show_column_names = FALSE,
  title = "Saccharomyces cerevisiae",
  row_title = "Periodic genes",
  column_title = "Time (minutes)",
  use_raster = FALSE,
  bottom_annotation = bottom_ha
)

pdf("figure_2A_recreated-CBY-reverse-1246-ticks.pdf", width = 9, height = 12)
draw(ht)

decorate_annotation("time", {
  x <- (tick_at - 0.5) / ncol(expr_z)
  
  # tick marks
  grid.segments(
    x0 = unit(x, "npc"),
    y0 = unit(1, "npc"),
    x1 = unit(x, "npc"),
    y1 = unit(0.72, "npc"),
    gp = gpar(lwd = 1)
  )
  
  # labels
  grid.text(
    label = tick_times,
    x = unit(x, "npc"),
    y = unit(0.25, "npc"),
    gp = gpar(fontsize = 10)
  )
})

dev.off()







################### 
# anciennes versions

# Heatmap
# inversion matrice 
expr_z <- expr_z[nrow(expr_z):1, , drop = FALSE]

ht <- Heatmap(
  expr_z,
  name = "z-score",
  col = colorRamp2(c(-1.5, 0, 1.5), c("cyan", "black", "yellow")),
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_names = FALSE,
  show_column_names = TRUE,
  column_labels = as.character(time_vals),
  row_title = "Top Periodic genes (1246)",
  column_title = "Time (minutes)",
  use_raster = FALSE
#  raster_by_magick = FALSE,
#  raster_quality = 1
)
#png("figure_2A_recreated-vFiltered.png", width = 1800, height = 2400, res = 200)
#draw(ht)
#dev.off()

pdf("figure_2A_recreatedvFiltered-CBY-select-reverse.pdf", width = 9, height = 12)
draw(ht)
dev.off()

#________ essai 1, gén,ération d'une erreur due au raster
# Heatmap
#png("figure_2A_recreated.png", width = 1800, height = 2400, res = 200)

#Heatmap(
#  expr_z,
#  name = "z-score",
#  col = colorRamp2(c(-1.5, 0, 1.5), c("navy", "white", "firebrick")),
#  cluster_rows = FALSE,
#  cluster_columns = FALSE,
#  show_row_names = FALSE,
#  show_column_names = TRUE,
#  column_labels = as.character(time_vals),
#  row_title = "Periodic genes",
#  column_title = "Time (minutes)",
#  use_raster = TRUE,
#  raster_quality = 2
#)

#dev.off()
