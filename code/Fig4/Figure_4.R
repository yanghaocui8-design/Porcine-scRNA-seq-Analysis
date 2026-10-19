library(Seurat)
library(monocle3)
library(ggplot2)
library(ComplexHeatmap)
library(circlize)
library(SingleCellExperiment)
library(grid)
library(dplyr)

input_rds <- file.path("outputs", "Figure_2", "pbmc_integrated.rds")
output_dir <- file.path("outputs", "Figure_4")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Figure 4A ------------------------------------------
pbmc <- readRDS(input_rds)
seurat_sub <- subset(pbmc, subset = Stage %in% c("Morula", "TE", "ICM", "EPI", "PE"))

stage_order <- c("Morula", "ICM", "TE", "EPI", "PE")
seurat_sub$Stage <- factor(seurat_sub$Stage, levels = stage_order)
stage_colors <- c(
  "Morula" = "#4D69AF",
  "ICM" = "#60C261",
  "TE" = "#1CAFC9",
  "EPI" = "#837FBD",
  "PE" = "#B76C9F"
)

expression_matrix <- GetAssayData(seurat_sub, assay = "RNA", layer = "counts")
cell_metadata <- seurat_sub@meta.data
gene_annotation <- data.frame(gene_short_name = rownames(expression_matrix), row.names = rownames(expression_matrix))
cds <- new_cell_data_set(expression_matrix, cell_metadata = cell_metadata, gene_metadata = gene_annotation)

cds <- preprocess_cds(cds, num_dim = 15)
cds <- reduce_dimension(cds, reduction_method = "UMAP", umap.min_dist = 0.5, umap.n_neighbors = 30)
cds <- cluster_cells(cds)
cds <- learn_graph(cds, use_partition = FALSE, learn_graph_control = list(minimal_branch_len = 4, rann.k = 35))

get_earliest_principal_node <- function(cds, time_bin="Morula") {
  cell_ids <- which(colData(cds)$Stage == time_bin)
  closest_vertex <- cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex)
  root_pr_nodes <- igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names(which.max(table(closest_vertex[cell_ids,]))))]
  return(root_pr_nodes)
}
nodes <- get_earliest_principal_node(cds, time_bin="Morula")
cds <- order_cells(cds, root_pr_nodes = nodes)

saveRDS(cds, file = file.path(output_dir, "monocle3_ordered_cds.rds"))

p4A <- plot_cells(
  cds,
  color_cells_by = "pseudotime",
  cell_size = 1.5,
  label_cell_groups = FALSE,
  label_leaves = FALSE,
  label_branch_points = FALSE
) +
  scale_color_viridis_c(option = "plasma", name = "Pseudotime") +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  labs(title = NULL, x = "UMAP 1", y = "UMAP 2")

print(p4A)
ggsave(file.path(output_dir, "Figure_4A_Pseudotime_Tree.pdf"), plot = p4A, width = 6.5, height = 4.5)
ggsave(file.path(output_dir, "Figure_4A_Pseudotime_Tree.png"), plot = p4A, width = 6.5, height = 4.5, dpi = 300)

# Figure 4B ------------------------------------------
p4B <- plot_cells(
  cds,
  color_cells_by = "Stage",
  cell_size = 1.5,
  label_cell_groups = FALSE,
  label_leaves = FALSE,
  label_branch_points = FALSE
) +
  scale_color_manual(values = stage_colors, breaks = stage_order, name = "Stage") +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  labs(title = NULL, x = "UMAP 1", y = "UMAP 2")

print(p4B)
ggsave(file.path(output_dir, "Figure_4B_Celltype_Tree.pdf"), plot = p4B, width = 6.5, height = 4.5)
ggsave(file.path(output_dir, "Figure_4B_Celltype_Tree.png"), plot = p4B, width = 6.5, height = 4.5, dpi = 300)

# Figure 4C ------------------------------------------
umap_coords <- as.data.frame(SingleCellExperiment::reducedDim(cds, "UMAP"))
colnames(umap_coords) <- c("UMAP_1", "UMAP_2")
cell_metadata <- as.data.frame(colData(cds))
plot_data <- cbind(
  umap_coords,
  Stage = cell_metadata[rownames(umap_coords), "Stage"]
)
plot_data <- plot_data %>% dplyr::filter(!is.na(Stage))
plot_data$Stage <- factor(plot_data$Stage, levels = stage_order)

set.seed(42)
plot_data_shuffled <- plot_data[sample(nrow(plot_data)), ]

p4C <- ggplot(plot_data_shuffled, aes(x = UMAP_1, y = UMAP_2, color = Stage)) +
  geom_point(
    size = 2.6,
    position = position_jitter(width = 0.28, height = 0.28, seed = 42)
  ) +
  scale_color_manual(values = stage_colors, breaks = stage_order, name = "Stage") +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.text = element_text(size = 12, color = "black"),
    axis.title = element_text(size = 14),
    legend.text = element_text(size = 11)
  ) +
  labs(title = NULL, x = "UMAP 1", y = "UMAP 2") +
  guides(color = guide_legend(override.aes = list(size = 4, alpha = 1)))
p4C$layers[[1]]$aes_params$alpha <- 0.85

print(p4C)
ggsave(file.path(output_dir, "Figure_4C_Cell_Type_Annotation.pdf"), plot = p4C, width = 8.5, height = 7)
ggsave(file.path(output_dir, "Figure_4C_Cell_Type_Annotation.png"), plot = p4C, width = 8.5, height = 7, dpi = 300)

# Figure 4D ------------------------------------------
essential_markers <- c("POU5F1", "SOX2", "NANOG", "NANOG-1", "ERAS",
                       "CDX2", "GATA3", "KRT8", "TEAD4", "BMP4",
                       "GATA6", "SOX17", "PDGFRA", "GATA4")
existing_markers <- essential_markers[essential_markers %in% rownames(cds)]

track_genes <- graph_test(cds, neighbor_graph="principal_graph", cores=4)
valid_bg_genes <- track_genes %>%
  dplyr::filter(!grepl("^LOC|^Gm|^ENSSSCG|^2310|^Format|^Rik", gene_short_name)) %>%
  dplyr::filter(!gene_short_name %in% existing_markers) %>%
  dplyr::filter(q_value < 1e-12) %>%
  dplyr::arrange(desc(morans_I)) %>%
  head(60) %>%
  pull(gene_short_name)

genes_subset <- unique(c(existing_markers, valid_bg_genes))

cell_ptime <- pseudotime(cds)
valid_cells <- names(cell_ptime)[!is.infinite(cell_ptime)]
ordered_cells <- names(sort(cell_ptime[valid_cells]))

raw_counts <- normalized_counts(cds)[genes_subset, ordered_cells]

pt.matrix <- t(apply(raw_counts, 1, function(x) {
  smooth.spline(x, df = 3)$y
}))

pt.matrix <- t(apply(pt.matrix, 1, function(x) {
  if (sd(x) == 0) return(x)
  (x - mean(x)) / sd(x)
}))

rownames(pt.matrix) <- genes_subset
colnames(pt.matrix) <- ordered_cells

ha_top <- HeatmapAnnotation(
  Pseudotime = cell_ptime[ordered_cells],
  col = list(Pseudotime = colorRamp2(
    c(0, max(cell_ptime[ordered_cells]) / 2, max(cell_ptime[ordered_cells])),
    c("#2C53A0", "#A52C91", "#F22232")
  )),
  show_legend = TRUE,
  annotation_name_side = "left"
)

col_fun <- colorRamp2(c(-2, 0, 2), c("#2166ac", "#f7f7f7", "#b2182b"))

htkm <- Heatmap(
  pt.matrix,
  name = "z-score",
  col = col_fun,
  top_annotation = ha_top,
  show_row_names = TRUE,
  show_column_names = FALSE,
  row_names_gp = gpar(fontsize = 6, fontface = "italic"),
  row_title = NULL,
  row_title_rot = 0,
  show_row_dend = TRUE,
  cluster_columns = FALSE,
  cluster_rows = TRUE,
  row_km = 3,
  row_gap = unit(3, "mm"),
  border = FALSE
)

png(file.path(output_dir, "Figure_4D_Pseudotime_Heatmap.png"), width = 7.5, height = 10, units = "in", res = 300)
draw(htkm)
dev.off()

pdf(file.path(output_dir, "Figure_4D_Pseudotime_Heatmap.pdf"), width = 7.5, height = 10)
draw(htkm)
dev.off()
