library(Seurat)
library(pheatmap)

input_rds <- file.path("data", "pbmc_qc.rds")
output_dir <- file.path("outputs", "Figure_1")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Figure 1 ------------------------------------------
pbmc <- readRDS(input_rds)
hvgs <- VariableFeatures(pbmc)
exp_mat <- as.matrix(pbmc[["RNA"]]$data[hvgs, ])

stage_factors <- factor(pbmc$Stage, levels = c("MII_Oocyte", "1_cell", "2_cell", "4_cell", "8_cell", "Morula", "ICM", "TE", "EPI", "PE"))

pseudo_bulk <- sapply(levels(stage_factors), function(stage) {
  cells_in_stage <- colnames(pbmc)[stage_factors == stage]
  if (length(cells_in_stage) > 1) {
    return(rowMeans(exp_mat[, cells_in_stage]))
  } else if (length(cells_in_stage) == 1) {
    return(exp_mat[, cells_in_stage])
  } else {
    return(numeric(length(hvgs)))
  }
})

cor_matrix_stages <- cor(pseudo_bulk, method = "pearson")

annotation_data <- data.frame(Stage = rownames(cor_matrix_stages))
rownames(annotation_data) <- rownames(cor_matrix_stages)

heatmap_stage_colors <- list(
  Stage = c(
    "MII_Oocyte" = "#E41A1C", "1_cell" = "#FF7F00", "2_cell" = "#377EB8", "4_cell" = "#4DAF4A",
    "8_cell" = "#984EA3", "Morula" = "#A65628", "ICM" = "#F781BF",
    "TE" = "#999999", "EPI" = "#17BECF", "PE" = "#BCBD22"
  )
)

pheatmap::pheatmap(
  cor_matrix_stages,
  annotation_col = annotation_data,
  annotation_row = annotation_data,
  annotation_colors = heatmap_stage_colors,
  show_rownames = TRUE,
  show_colnames = TRUE,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  legend = TRUE,
  annotation_legend = FALSE,
  color = colorRampPalette(c("#3A7BD5", "#FFFFFF", "#D32F2F"))(100),
  border_color = "white",
  main = "Developmental Stage Correlation Matrix (Pseudo-bulk)",
  filename = file.path(output_dir, "Figure_1_Stage_Correlation.pdf")
)
