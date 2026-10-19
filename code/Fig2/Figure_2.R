library(Seurat)
library(harmony)
library(dplyr)
library(stringr)
library(ggplot2)
library(patchwork)
library(tidyr)

input_rds <- file.path("data", "pbmc_qc.rds")
heatmap_genes_file <- file.path("data", "Figure_2E_genes.txt")
output_dir <- file.path("outputs", "Figure_2")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Figure 2A ------------------------------------------
pbmc <- readRDS(input_rds)

stage_order <- c("MII_Oocyte", "1_cell", "2_cell", "4_cell", "8_cell", "Morula", "ICM", "TE", "EPI", "PE")
pbmc$Stage <- factor(pbmc$Stage, levels = stage_order)

pbmc$Batch <- case_when(
  str_detect(rownames(pbmc@meta.data), "^SRR") ~ str_extract(rownames(pbmc@meta.data), "^SRR\\d{5}"),
  TRUE                                         ~ "Custom_Batch_Cell1"
)

pbmc <- JoinLayers(pbmc)
pbmc <- NormalizeData(pbmc)

pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 5000)
pbmc <- ScaleData(pbmc, features = VariableFeatures(pbmc), verbose = FALSE)

pbmc <- RunPCA(pbmc, features = VariableFeatures(pbmc), npcs = 30, verbose = FALSE)

pbmc <- RunHarmony(pbmc,
                   group.by.vars = "Batch",
                   reduction.use = "pca",
                   dims.use = 1:30,
                   theta = 0.8,
                   lambda = 1,
                   max_iter = 25,
                   verbose = TRUE)

pbmc <- FindNeighbors(pbmc, reduction = "harmony", dims = 1:30, verbose = FALSE)

pbmc <- FindClusters(pbmc, resolution = 0.8, verbose = FALSE)

pbmc <- RunUMAP(pbmc,
                reduction = "harmony",
                dims = 1:30,
                n.neighbors = 50L,
                min.dist = 0.4,
                verbose = FALSE)

stage_colors <- c(
  "MII_Oocyte" = "#E41A1C", "1_cell" = "#FF7F00", "2_cell" = "#377EB8",
  "4_cell"     = "#4DAF4A", "8_cell" = "#984EA3", "Morula" = "#A65628",
  "ICM"        = "#F781BF", "TE"     = "#999999", "EPI"    = "#17BECF", "PE" = "#BCBD22"
)

stage_labels <- c(
  "MII_Oocyte" = "MII", "1_cell" = "1-cell", "2_cell" = "2-cell",
  "4_cell" = "4-cell", "8_cell" = "8-cell", "Morula" = "Morula",
  "ICM" = "ICM", "TE" = "TE", "EPI" = "EPI", "PE" = "PE"
)

p2_stage <- DimPlot(pbmc, reduction = "umap", group.by = "Stage", pt.size = 1.5) +
  scale_color_manual(
    values = stage_colors,
    breaks = stage_order,
    labels = unname(stage_labels[stage_order])
  ) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
    legend.title = element_blank()
  ) +
  labs(title = NULL, x = "UMAP-1", y = "UMAP-2")

print(p2_stage)
Idents(pbmc) <- "Stage"
saveRDS(pbmc, file = file.path(output_dir, "pbmc_integrated.rds"))
ggsave(file.path(output_dir, "Figure_2A_Stages.pdf"), plot = p2_stage, width = 7, height = 5)
ggsave(file.path(output_dir, "Figure_2A_Stages.png"), plot = p2_stage, width = 7, height = 5, dpi = 300)

# Figure 2B ------------------------------------------
features_figure2 <- c(
  "ZAR1", "DPPA3", "C3orf70", "PALS1", "PLEK", "MPIG6B",
  "LEUTX", "ZSCAN4", "DUXA", "NR5A2", "CD37", "KLF4",
  "AQP1", "C1QA", "POU5F1", "CDX2", "GATA3", "GATA2",
  "SOX2", "NANOG", "ERAS", "GATA6", "PDGFRA", "LAMA1"
)

existing_features <- features_figure2[features_figure2 %in% rownames(pbmc)]
missing_features <- setdiff(features_figure2, existing_features)
if (length(missing_features) > 0) {
  warning("Genes absent from the expression matrix: ", paste(missing_features, collapse = ", "), call. = FALSE)
}
if (length(existing_features) == 0) stop("No Figure 2B genes match the expression matrix.")

p_dotplot <- DotPlot(pbmc, features = existing_features, assay = "RNA") +
  coord_flip() +
  scale_y_discrete(labels = stage_labels) +
  scale_color_gradientn(colors = c("#F5F5F5", "#FFB90F", "#D62728")) +
  theme_bw() +
  theme(
    panel.grid.major = element_line(color = "grey96"),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10, color = "black"),
    axis.text.y = element_text(size = 10, color = "black"),
    legend.position = "right"
  ) +
  labs(x = "Lineage and stage marker genes", y = NULL, color = "Relative expression", size = "Pct. expressed")

print(p_dotplot)
ggsave(file.path(output_dir, "Figure_2B_DotPlot.pdf"), plot = p_dotplot, width = 9.5, height = 7.5)
ggsave(file.path(output_dir, "Figure_2B_DotPlot.png"), plot = p_dotplot, width = 9.5, height = 7.5, dpi = 300)

# Figure 2C ------------------------------------------
gene_stages <- c(
  "ZAR1" = "MII_Oocyte",
  "PALS1" = "1_cell",
  "PLEK" = "2_cell",
  "LEUTX" = "4_cell",
  "CD37" = "8_cell",
  "AQP1" = "Morula",
  "C1QA" = "ICM",
  "GATA2" = "TE",
  "ERAS" = "EPI",
  "LAMA1" = "PE"
)

plot_list <- list()
for (gene in names(gene_stages)) {
  if (gene %in% rownames(pbmc)) {
    stage_label <- stage_labels[gene_stages[gene]]
    p <- FeaturePlot(
      pbmc,
      features = gene,
      reduction = "umap",
      cols = c("#EFEFEF", "#D62728"),
      order = TRUE,
      pt.size = 0.4
    ) +
      theme_bw() +
      theme(
        panel.grid = element_blank(),
        plot.title = element_text(size = 11, hjust = 0.5),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.border = element_rect(color = "grey80", fill = NA, linewidth = 0.6),
        legend.position = "right"
      ) +
      labs(title = paste0(gene, " (", stage_label, ")"))
    plot_list[[gene]] <- p
  }
}
if (length(plot_list) == 0) stop("No Figure 2C genes match the expression matrix.")

combined_feature_plots <- wrap_plots(plot_list, ncol = 2)
plot_height <- ceiling(length(plot_list) / 2) * 3.2
print(combined_feature_plots)
ggsave(
  file.path(output_dir, "Figure_2C_FeaturePlot.pdf"),
  plot = combined_feature_plots,
  width = 8.5,
  height = plot_height,
  device = cairo_pdf,
  limitsize = FALSE
)
ggsave(
  file.path(output_dir, "Figure_2C_FeaturePlot.png"),
  plot = combined_feature_plots,
  width = 8.5,
  height = plot_height,
  dpi = 300,
  limitsize = FALSE
)

# Figure 2E ------------------------------------------
if (!file.exists(heatmap_genes_file)) {
  stop("Provide the original Figure 2E gene list in data/Figure_2E_genes.txt (one gene per line).")
}
final_features <- scan(heatmap_genes_file, what = character(), quiet = TRUE)
if (length(final_features) == 0L) stop("The Figure 2E gene list is empty.")

pbmc$Stage <- factor(pbmc$Stage, levels = stage_order)
Idents(pbmc) <- "Stage"

cells_ordered <- colnames(pbmc)[order(match(pbmc$Stage, stage_order))]
heatmap_genes <- final_features[final_features %in% rownames(pbmc)]
if (length(heatmap_genes) == 0L) stop("No Figure 2E genes match the expression matrix.")

pbmc <- ScaleData(pbmc, features = heatmap_genes)

mat <- GetAssayData(pbmc, assay = "RNA", layer = "scale.data")[heatmap_genes, cells_ordered]
mat[mat > 2.0] <- 2.0
mat[mat < -2.0] <- -2.0

cell_counts <- table(pbmc$Stage)[stage_order]
v_lines <- cumsum(cell_counts)[-length(cell_counts)]

ends <- cumsum(cell_counts)
starts <- c(1, ends[-length(ends)] + 1)
center_indices <- round((starts + ends) / 2)
center_cell_ids <- cells_ordered[center_indices]

label_mapping <- data.frame(
  Cell = cells_ordered,
  CleanLabel = ""
)
label_mapping$CleanLabel[center_indices] <- stage_labels[stage_order]

df <- as.data.frame(mat) %>%
  dplyr::mutate(Gene = rownames(.)) %>%
  pivot_longer(cols = -Gene, names_to = "Cell", values_to = "Expression") %>%
  dplyr::left_join(label_mapping, by = "Cell")

df$Gene <- factor(df$Gene, levels = rev(heatmap_genes))
df$Cell <- factor(df$Cell, levels = cells_ordered)

p2_heatmap <- ggplot(df) +
  geom_tile(aes(x = Cell, y = Gene, fill = Expression)) +
  geom_vline(xintercept = v_lines + 0.5, color = "white", linewidth = 0.5) +
  scale_fill_gradientn(
    colors = c("#313695", "#4575b4", "#74add1", "#abd9e9", "#ffffff", "#fee090", "#fdae61", "#f46d43", "#d73027"),
    name = "Expression"
  ) +
  scale_x_discrete(
    breaks = center_cell_ids,
    labels = label_mapping$CleanLabel[center_indices],
    expand = c(0, 0)
  ) +
  scale_y_discrete(expand = c(0, 0)) +
  theme_minimal() +
  theme(

    axis.text.x = element_text(face = "bold", size = 10, color = "black",
                               angle = 45, hjust = 1, vjust = 1),
    axis.ticks.x = element_blank(),
    axis.title.x = element_blank(),

    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title.y = element_blank(),
    panel.grid = element_blank(),

    plot.margin = margin(t = 15, r = 15, b = 60, l = 15, unit = "pt"),
    plot.title = element_text(face = "bold", size = 13, hjust = 0.5, vjust = 1)
  )

print(p2_heatmap)

out_width <- 12.5
out_height <- 9.0

ggsave(file.path(output_dir, "Figure_2E_SingleCell_Heatmap.pdf"), plot = p2_heatmap, width = out_width, height = out_height, device = cairo_pdf)
ggsave(file.path(output_dir, "Figure_2E_SingleCell_Heatmap.png"), plot = p2_heatmap, width = out_width, height = out_height, dpi = 300)
