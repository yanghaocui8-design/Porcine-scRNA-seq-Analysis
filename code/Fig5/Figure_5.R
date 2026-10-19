library(Seurat)
library(ComplexHeatmap)
library(circlize)
library(openxlsx)
library(Mfuzz)
library(Biobase)
library(clusterProfiler)
library(org.Ss.eg.db)
library(grid)
library(dplyr)

input_rds <- file.path("outputs", "Figure_2", "pbmc_integrated.rds")
output_dir <- file.path("outputs", "Figure_5")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Figure 5A ------------------------------------------
pbmc_lazy <- readRDS(input_rds)

stage_order_5s <- c("Morula", "ICM", "TE", "EPI", "PE")

pbmc_5stages <- subset(pbmc_lazy, subset = Stage %in% stage_order_5s)
pbmc_5stages$Stage <- factor(pbmc_5stages$Stage, levels = stage_order_5s)
Idents(pbmc_5stages) <- "Stage"
rm(pbmc_lazy)

zga_all_markers <- FindAllMarkers(pbmc_5stages, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.5)

strict_rescue <- FindMarkers(pbmc_5stages, ident.1 = "Morula", ident.2 = "PE", only.pos = TRUE, logfc.threshold = 0.25)
rescue_genes_all <- rownames(strict_rescue)[which(strict_rescue$p_val_adj < 0.05)]
extra_genes <- head(rescue_genes_all, min(30, length(rescue_genes_all)))

heatmap_genes <- unique(c(zga_all_markers$gene[zga_all_markers$p_val_adj < 0.05], extra_genes))
heatmap_genes <- heatmap_genes[!grepl("^LOC|^ENSSSCG", heatmap_genes)]
heatmap_genes <- unique(heatmap_genes[heatmap_genes %in% rownames(pbmc_5stages)])

exp_mtx <- as.matrix(GetAssayData(pbmc_5stages, assay = "RNA", layer = "data")[heatmap_genes, ])
avg_mtx <- sapply(levels(pbmc_5stages$Stage), function(st) rowMeans(exp_mtx[, which(Idents(pbmc_5stages) == st)]))
plot_mtx <- t(scale(t(avg_mtx)))

peak_st <- max.col(plot_mtx, ties.method = "first")

sort_idx <- order(peak_st)
plot_mtx_final <- plot_mtx[sort_idx, , drop = FALSE]
row_split_final <- peak_st[sort_idx]

palette_colors_5s <- c("#4D69AF", "#60C261", "#1CAFC9", "#837FBD", "#B76C9F")
names(palette_colors_5s) <- stage_order_5s

ht <- Heatmap(plot_mtx_final,
              name = "Z-Score",
              col = colorRamp2(c(-2, 0, 2), c("#3552A0", "#FFFFFF", "#D73027")),

              top_annotation = HeatmapAnnotation(
                Stage = levels(pbmc_5stages$Stage),
                col = list(Stage = palette_colors_5s),
                gp = gpar(col = NA),
                show_annotation_name = FALSE
              ),

              cluster_rows = TRUE,
              cluster_row_slices = FALSE,
              row_dend_width = unit(15, "mm"),

              cluster_columns = FALSE,

              row_split = row_split_final,
              rect_gp = gpar(col = NA),
              border = TRUE,
              row_gap = unit(2.5, "mm"),

              row_title = "Cluster %s",
              row_title_rot = 90,
              row_title_gp = gpar(fontsize = 11, fontface = "plain"),

              column_names_rot = 0,
              column_names_centered = TRUE,
              column_names_gp = gpar(fontsize = 11, fontface = "bold"),

              show_row_names = FALSE,
              use_raster = TRUE,
              column_title = NULL,
              column_title_gp = gpar(fontsize = 12, fontface = "bold")
)

pdf(file.path(output_dir, "Figure_5A_Lineage_Cascade.pdf"), width = 6.5, height = 8.8)
draw(ht, merge_legend = TRUE)
dev.off()

png(file.path(output_dir, "Figure_5A_Lineage_Cascade.png"), width = 6.5, height = 8.8, units = "in", res = 300)
draw(ht, merge_legend = TRUE)
dev.off()

draw(ht, merge_legend = TRUE)

avg_df <- as.data.frame(avg_mtx) %>%
  tibble::rownames_to_column("gene")

cascade_info <- data.frame(
  gene = rownames(plot_mtx_final),
  Peak_Stage_Num = row_split_final,
  Peak_Stage_Name = stage_order_5s[row_split_final],
  stringsAsFactors = FALSE
)

final_excel_data <- zga_all_markers %>%
  dplyr::filter(gene %in% heatmap_genes) %>%
  dplyr::left_join(cascade_info, by = "gene") %>%
  dplyr::left_join(avg_df, by = "gene") %>%
  dplyr::select(gene, cluster, p_val, avg_log2FC, pct.1, pct.2, p_val_adj,
         Peak_Stage_Num, Peak_Stage_Name, everything()) %>%
  dplyr::arrange(Peak_Stage_Num, desc(avg_log2FC))

write.xlsx(final_excel_data, file = file.path(output_dir, "Lineage_Table4_Markers.xlsx"), overwrite = TRUE)

# Figure 5B ------------------------------------------
time_series_5s <- avg_mtx

eset <- ExpressionSet(assayData = as.matrix(time_series_5s))

eset <- filter.std(eset, min.std = 0, visu = FALSE)
eset <- standardise(eset)

m1 <- mestimate(eset)
num_clusters <- 5

set.seed(123)
cl_5s <- mfuzz(eset, c = num_clusters, m = m1)

plot_mfuzz <- function() {
  par(mar = c(4, 4, 2, 1))
  mfuzz.plot2(
    eset, cl_5s,
    mfrow = c(3, 2),
    time.labels = colnames(time_series_5s),
    xlab = "Developmental Stages",
    ylab = "Expression Change (Z-score)",
    colo = "fancy",
    x11 = FALSE,
    lwd = 1.8
  )
}

pdf(file.path(output_dir, "Figure_5B_Mfuzz_5_Clusters.pdf"), width = 8, height = 9)
tryCatch(plot_mfuzz(), finally = dev.off())
png(file.path(output_dir, "Figure_5B_Mfuzz_5_Clusters.png"), width = 8, height = 9, units = "in", res = 300)
tryCatch(plot_mfuzz(), finally = dev.off())

cluster_df <- data.frame(
  Gene = rownames(Biobase::exprs(eset)),
  Mfuzz_Cluster = cl_5s$cluster
)
write.csv(cluster_df, file.path(output_dir, "Figure_5B_Mfuzz_5_Clusters_Genes.csv"), row.names = FALSE)

# Figure 5C ------------------------------------------
base_dir <- output_dir
out_dir  <- file.path(base_dir, "GO_Results")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

mfuzz_file <- file.path(base_dir, "Figure_5B_Mfuzz_5_Clusters_Genes.csv")
mfuzz_df <- read.csv(mfuzz_file, stringsAsFactors = FALSE)

cluster_map <- c("3" = "C1", "2" = "C2", "5" = "C3", "4" = "C4", "1" = "C5")
mfuzz_df$Ordered_Cluster <- cluster_map[as.character(mfuzz_df$Mfuzz_Cluster)]

go_results_list <- list()

for (cl in c("C1", "C2", "C3", "C4", "C5")) {
  genes <- mfuzz_df %>% dplyr::filter(Ordered_Cluster == cl) %>% pull(Gene)

  if (length(genes) > 0) {

    ego <- enrichGO(gene          = genes,
                    OrgDb         = org.Ss.eg.db,
                    keyType       = 'SYMBOL',
                    ont           = "BP",
                    pAdjustMethod = "BH",
                    pvalueCutoff  = 0.2,
                    qvalueCutoff  = 0.5)

    if (!is.null(ego) && nrow(as.data.frame(ego)) > 0) {
      res_df <- as.data.frame(ego)
      res_df$Cluster <- cl
      go_results_list[[cl]] <- res_df
    } else {

      ego_raw <- enrichGO(gene          = genes,
                          OrgDb         = org.Ss.eg.db,
                          keyType       = 'SYMBOL',
                          ont           = "BP",
                          pvalueCutoff  = 1,
                          qvalueCutoff  = 1)
      if (!is.null(ego_raw)) {
        res_df <- as.data.frame(ego_raw)
        res_df$Cluster <- cl
        go_results_list[[cl]] <- res_df
      }
    }
  }
}

all_go_res <- do.call(rbind, go_results_list)
write.xlsx(all_go_res, file.path(out_dir, "Figure_5C_GO_All_Results.xlsx"), overwrite = TRUE)

top5_go <- all_go_res %>%
  dplyr::group_by(Cluster) %>%
  dplyr::arrange(pvalue) %>%
  slice_head(n = 5) %>%
  dplyr::ungroup()

selected_terms <- unique(top5_go$Description)

go_matrix <- matrix(0, nrow = length(selected_terms), ncol = 5,
                    dimnames = list(selected_terms, c("C1", "C2", "C3", "C4", "C5")))

for (i in seq_len(nrow(all_go_res))) {
  term <- all_go_res$Description[i]
  cl <- all_go_res$Cluster[i]
  if (term %in% selected_terms) {

    p_val <- all_go_res$pvalue[i]
    go_matrix[term, cl] <- -log10(p_val)
  }
}

term_split <- factor(top5_go$Cluster[match(rownames(go_matrix), top5_go$Description)],
                     levels = c("C1", "C2", "C3", "C4", "C5"))

go_display <- matrix(
  "Not enriched",
  nrow = nrow(go_matrix), ncol = ncol(go_matrix),
  dimnames = dimnames(go_matrix)
)
for (cl in colnames(go_matrix)) {
  enriched <- which(go_matrix[, cl] >= 1.3)
  go_display[enriched, cl] <- cl
}
go_colors <- c(
  "Not enriched" = "#D9D9D9",
  setNames(unname(palette_colors_5s), c("C1", "C2", "C3", "C4", "C5"))
)

column_labels <- c("C1\nMorula", "C2\nICM", "C3\nTE", "C4\nEPI", "C5\nPE")

ht_go <- Heatmap(go_display,
                 name = "Enriched",
                 col = go_colors,
                 cluster_rows = FALSE,
                 cluster_columns = FALSE,
                 row_split = term_split,
                 row_title = NULL,
                 row_title_rot = 0,
                 row_title_gp = gpar(fontsize = 10, fontface = "bold"),
                 row_gap = unit(2, "mm"),
                 rect_gp = gpar(col = "black", lwd = 0.8),
                 column_labels = column_labels,
                 column_names_rot = 0,
                 column_names_centered = TRUE,
                 column_names_gp = gpar(fontsize = 10, fontface = "bold"),
                 row_names_gp = gpar(fontsize = 9),
                 width = unit(6, "cm"),
                 show_heatmap_legend = FALSE,
                 column_title = NULL,
                 column_title_gp = gpar(fontsize = 12, fontface = "bold")
)

pdf_path <- file.path(out_dir, "Figure_5C_GO_Enrichment_Heatmap.pdf")
pdf(pdf_path, width = 10, height = 11)
draw(ht_go)
dev.off()

png(file.path(out_dir, "Figure_5C_GO_Enrichment_Heatmap.png"), width = 10, height = 11, units = "in", res = 300)
draw(ht_go)
dev.off()
