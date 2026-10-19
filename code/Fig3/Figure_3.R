library(Seurat)
library(ggplot2)
library(patchwork)
library(ggrepel)
library(readxl)
library(clusterProfiler)
library(org.Ss.eg.db)
library(ComplexHeatmap)
library(circlize)
library(grid)
library(dplyr)

input_rds <- file.path("outputs", "Figure_2", "pbmc_integrated.rds")
target_excel_path <- file.path("data", "ZGA_Table4_Markers.xlsx")
output_dir <- file.path("outputs", "Figure_3")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
kegg_top_terms <- 3L
go_top_terms <- 5L
chord_top_genes <- 20L

# Figure 3A ------------------------------------------
pbmc <- readRDS(input_rds)

possible_cols <- grep("stage", colnames(pbmc@meta.data), ignore.case = TRUE, value = TRUE)
actual_stage_col <- if (length(possible_cols) > 0) possible_cols[1] else "active.ident"

if (actual_stage_col != "active.ident") {
  pbmc$Stage_Cleaned <- pbmc@meta.data[[actual_stage_col]]
  Idents(pbmc) <- "Stage_Cleaned"
} else {
  pbmc$Stage_Cleaned <- Idents(pbmc)
}

all_real_names <- as.character(unique(Idents(pbmc)))

c1_name <- grep("1.*cell|zygote", all_real_names, ignore.case = TRUE, value = TRUE)[1]
c2_name <- grep("2.*cell|2C", all_real_names, ignore.case = TRUE, value = TRUE)[1]
c4_name <- grep("4.*cell|4C", all_real_names, ignore.case = TRUE, value = TRUE)[1]
c8_name <- grep("8.*cell|8C", all_real_names, ignore.case = TRUE, value = TRUE)[1]

comparisons <- list(
  c(c1_name, c2_name, "2C vs 1C"),
  c(c2_name, c4_name, "4C vs 2C"),
  c(c4_name, c8_name, "8C vs 4C")
)

plot_list <- list()

FC_THRESHOLD   <- 1.0
PVAL_THRESHOLD <- 0.01

for (comp in comparisons) {
  g1 <- comp[1]
  g2 <- comp[2]
  title_name <- comp[3]

  if (!is.na(g1) && !is.na(g2) && g1 %in% Idents(pbmc) && g2 %in% Idents(pbmc)) {

    deg <- FindMarkers(pbmc, ident.1 = g2, ident.2 = g1, logfc.threshold = 0, min.pct = 0.05)
    deg$Gene <- rownames(deg)

    deg$p_val[deg$p_val == 0] <- 1e-300
    deg$neg_log10_pval <- -log10(deg$p_val)

    deg$Status <- "Unchanged"
    deg$Status[deg$avg_log2FC > FC_THRESHOLD & deg$p_val < PVAL_THRESHOLD] <- "Up (ZGA)"
    deg$Status[deg$avg_log2FC < -FC_THRESHOLD & deg$p_val < PVAL_THRESHOLD] <- "Down (Decay)"

    n_up   <- sum(deg$Status == "Up (ZGA)")
    n_down <- sum(deg$Status == "Down (Decay)")

    up_label   <- paste0("Up (ZGA): ", n_up)
    down_label <- paste0("Down (Decay): ", n_down)

    deg$Status_Legend <- "Unchanged"
    deg$Status_Legend[deg$Status == "Up (ZGA)"] <- up_label
    deg$Status_Legend[deg$Status == "Down (Decay)"] <- down_label

    mt_ribo_pattern <- "^loc|^mt-?|^rp|^cox|^nd[1-6]|^atp[68]|^cytb"

    clean_deg <- dplyr::filter(deg, !grepl(mt_ribo_pattern, Gene, ignore.case = TRUE))

    top_up   <- dplyr::filter(clean_deg, Status == "Up (ZGA)") %>%
      dplyr::arrange(p_val, desc(avg_log2FC)) %>%
      head(3)

    top_down <- dplyr::filter(clean_deg, Status == "Down (Decay)") %>%
      dplyr::arrange(p_val, avg_log2FC) %>%
      head(3)

    label_data <- rbind(top_up, top_down)

    color_map <- c("#D73027", "#4575B4", "#CCCCCC")
    names(color_map) <- c(up_label, down_label, "Unchanged")

    p <- ggplot(deg, aes(x = avg_log2FC, y = neg_log10_pval)) +
      geom_point(aes(color = Status_Legend), alpha = 0.35, size = 0.5) +
      scale_color_manual(values = color_map) +
      geom_vline(xintercept = c(-FC_THRESHOLD, FC_THRESHOLD), linetype = "dashed", color = "grey50", linewidth = 0.3) +
      geom_hline(yintercept = -log10(PVAL_THRESHOLD), linetype = "dashed", color = "grey50", linewidth = 0.3) +
      geom_text_repel(
        data = label_data, aes(label = Gene), size = 3, fontface = "bold.italic",
        color = "black", max.overlaps = 30, box.padding = 0.5, point.padding = 0.2,
        segment.size = 0.3, segment.color = "grey30"
      ) +
      labs(title = title_name, x = "Log2 Fold Change", y = "-Log10 P-value") +
      theme_classic() +
      theme(
        aspect.ratio = 1,
        plot.title = element_text(face = "bold", size = 11, hjust = 0.5),
        axis.title = element_text(face = "bold", size = 9.5),
        axis.text = element_text(size = 8.5, color = "black"),
        axis.line = element_line(linewidth = 0.4, color = "black"),
        legend.position = "bottom",
        legend.title = element_blank(),
        legend.text = element_text(size = 8)
      )

    plot_list[[title_name]] <- p
  }
}

if (length(plot_list) > 0) {
  combined_volcano <- wrap_plots(plot_list, ncol = 3)

  print(combined_volcano)
  ggsave(file.path(output_dir, "Figure_3A_Volcano.pdf"), plot = combined_volcano, width = 14, height = 5)
  ggsave(file.path(output_dir, "Figure_3A_Volcano.png"), plot = combined_volcano, width = 14, height = 5, dpi = 300)
}

# Figure 3D ------------------------------------------
pbmc_4stages <- subset(pbmc, subset = Stage %in% c("1_cell", "2_cell", "4_cell", "8_cell"))

pbmc_4stages <- RunUMAP(pbmc_4stages, dims = 1:30, reduction = "harmony", reduction.name = "umap")

pbmc_4stages$Stage <- factor(pbmc_4stages$Stage, levels = c("1_cell", "2_cell", "4_cell", "8_cell"))

zga_4_colors <- c(
  "1_cell" = "#FF7F00",
  "2_cell" = "#4DAF4A",
  "4_cell" = "#984EA3",
  "8_cell" = "#377EB8"
)

p3D <- DimPlot(
  pbmc_4stages,
  reduction = "umap",
  group.by = "Stage",
  cols = zga_4_colors,
  pt.size = 3.0,
  shuffle = TRUE,
  raster = FALSE
) +
  theme_classic() +
  labs(
    x = "UMAP_1",
    y = "UMAP_2",
    title = NULL
  ) +
  theme(
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    axis.text = element_text(size = 12, color = "black"),
    axis.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 13, face = "bold"),
    legend.spacing.y = unit(0.3, 'cm'),
    legend.title = element_blank()
  ) +
  guides(Stage = guide_legend(override.aes = list(size = 6, alpha = 1)))

p3D$layers[[1]]$aes_params$alpha <- 0.85

print(p3D)

ggsave(
  file.path(output_dir, "Figure_3D_ZGA_UMAP.pdf"),
  plot = p3D,
  width = 8.5,
  height = 7.0,
  dpi = 300
)

# Figure 3E ------------------------------------------
Idents(pbmc) <- "Stage_Cleaned"

deg_4C_vs_2C <- FindMarkers(
  object = pbmc,
  ident.1 = "4_cell",
  ident.2 = "2_cell",
  logfc.threshold = 0,
  min.pct = 0
)
deg_4C_vs_2C$Gene <- rownames(deg_4C_vs_2C)
deg_4C_vs_2C$Delta_Percent <- deg_4C_vs_2C$pct.1 - deg_4C_vs_2C$pct.2

FC_CUTOFF <- 1.0
DP_OFF <- 0.05

gene_sets <- list(
  ZGA_Activation_Up = deg_4C_vs_2C %>%
    dplyr::filter(avg_log2FC > FC_CUTOFF, Delta_Percent > DP_OFF),
  Maternal_Decay_Down = deg_4C_vs_2C %>%
    dplyr::filter(avg_log2FC < -FC_CUTOFF, Delta_Percent < -DP_OFF)
)

gene_inputs <- lapply(gene_sets, function(deg_df) {
  if (nrow(deg_df) == 0) {
    return(list(entrez = character(), deg_df = deg_df))
  }
  gene_map <- bitr(
    unique(deg_df$Gene),
    fromType = "SYMBOL",
    toType = "ENTREZID",
    OrgDb = org.Ss.eg.db
  )
  list(entrez = unique(na.omit(gene_map$ENTREZID)), deg_df = deg_df)
})

prepare_chord_data <- function(result_df, deg_df, n_terms, n_genes) {
  if (is.null(result_df) || nrow(result_df) == 0 || nrow(deg_df) == 0) return(NULL)
  result_df <- result_df[
    is.finite(result_df$p.adjust) & result_df$p.adjust >= 0 & result_df$p.adjust < 1,
    , drop = FALSE
  ]
  if (nrow(result_df) == 0) return(NULL)

  result_df <- result_df[order(result_df$p.adjust, result_df$pvalue), , drop = FALSE]
  top_terms <- head(result_df, n_terms)
  top_terms$Term <- paste0("P", seq_len(nrow(top_terms)))
  pathway_genes <- unique(unlist(strsplit(as.character(top_terms$geneID), "/", fixed = TRUE)))
  candidates <- deg_df[
    !is.na(deg_df$Gene) & deg_df$Gene %in% pathway_genes & is.finite(deg_df$avg_log2FC),
    , drop = FALSE
  ]
  candidates <- candidates[order(abs(candidates$avg_log2FC), decreasing = TRUE), , drop = FALSE]
  genes <- head(unique(as.character(candidates$Gene)), n_genes)
  if (length(genes) == 0) return(NULL)

  weights <- matrix(
    0, nrow = length(genes), ncol = nrow(top_terms),
    dimnames = list(genes, top_terms$Term)
  )
  for (i in seq_len(nrow(top_terms))) {
    term_genes <- strsplit(as.character(top_terms$geneID[i]), "/", fixed = TRUE)[[1]]
    matched_genes <- intersect(term_genes, genes)
    weights[matched_genes, i] <- -log10(max(top_terms$p.adjust[i], .Machine$double.xmin))
  }
  weights <- weights[rowSums(weights) > 0, colSums(weights) > 0, drop = FALSE]
  if (nrow(weights) == 0 || ncol(weights) == 0) return(NULL)

  top_terms <- top_terms[match(colnames(weights), top_terms$Term), , drop = FALSE]
  list(weights = weights, terms = top_terms)
}

save_enrichment_chord <- function(
  result_df, deg_df, panel, group_name, n_terms, legend_title,
  analysis_id = "4C_vs_2C", analysis_label = "4-cell vs 2-cell"
) {
  chord <- prepare_chord_data(result_df, deg_df, n_terms, chord_top_genes)
  if (is.null(chord)) {
    warning(paste(panel, group_name, "has no eligible gene-term links."), call. = FALSE)
    return(invisible(NULL))
  }

  weights <- chord$weights
  term_ids <- colnames(weights)
  gene_ids <- rownames(weights)
  palette <- c("#E5D2DD", "#53A85F", "#F1BB72", "#F3B1A0", "#D6E7A3")
  term_colors <- setNames(rep(palette, length.out = length(term_ids)), term_ids)
  gene_colors <- setNames(rep("#6CB8D2", length(gene_ids)), gene_ids)
  grid_colors <- c(term_colors, gene_colors)
  link_colors <- matrix(
    rep(unname(term_colors), each = nrow(weights)),
    nrow = nrow(weights), ncol = ncol(weights), dimnames = dimnames(weights)
  )
  sector_order <- c(term_ids, gene_ids)
  gap_degrees <- c(
    rep(2, length(term_ids) - 1), 15,
    rep(2, length(gene_ids) - 1), 15
  )

  prefix <- file.path(output_dir, paste0("Figure_", panel, "_", analysis_id, "_", group_name))
  write.csv(chord$terms, paste0(prefix, "_Terms.csv"), row.names = FALSE)
  write.csv(weights, paste0(prefix, "_Link_Weights.csv"))

  draw_chord <- function() {
    circos.clear()
    on.exit(circos.clear(), add = TRUE)
    par(mar = c(1, 1, 4, 1))
    circos.par(
      start.degree = 90,
      gap.degree = gap_degrees,
      track.margin = c(0.01, 0.01),
      canvas.xlim = c(-1.3, 1.3),
      canvas.ylim = c(-1.3, 1.8)
    )
    chordDiagram(
      weights,
      order = sector_order,
      grid.col = grid_colors,
      col = link_colors,
      transparency = 0.35,
      directional = 1,
      direction.type = "arrows",
      link.arr.type = "big.arrow",
      annotationTrack = "grid",
      annotationTrackHeight = mm_h(4)
    )
    circos.track(track.index = 1, panel.fun = function(x, y) {
      sector <- get.cell.meta.data("sector.index")
      xlim <- get.cell.meta.data("xlim")
      ylim <- get.cell.meta.data("ylim")
      circos.text(
        mean(xlim), mean(ylim), sector,
        facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5), cex = 0.75
      )
    }, bg.border = NA)
    legend(
      "topleft",
      legend = paste0(chord$terms$Term, ": ", chord$terms$Description),
      fill = unname(term_colors),
      title = legend_title,
      bty = "n",
      cex = 0.7
    )
    title(paste0("Figure ", panel, ": ", analysis_label, " (", group_name, ")"), cex.main = 1)
  }

  pdf(paste0(prefix, "_Chord.pdf"), width = 9, height = 9)
  tryCatch(draw_chord(), finally = dev.off())
  png(paste0(prefix, "_Chord.png"), width = 2700, height = 2700, res = 300)
  tryCatch(draw_chord(), finally = dev.off())
  invisible(chord)
}

for (group_name in names(gene_inputs)) {
  current <- gene_inputs[[group_name]]
  if (length(current$entrez) < 5) {
    warning(paste("3E", group_name, "has fewer than five mapped genes."), call. = FALSE)
    next
  }
  kegg_res <- tryCatch(
    enrichKEGG(
      gene = current$entrez,
      organism = "ssc",
      pAdjustMethod = "BH",
      pvalueCutoff = 0.05,
      qvalueCutoff = 0.05
    ),
    error = function(e) {
      warning(paste("KEGG enrichment failed:", conditionMessage(e)), call. = FALSE)
      NULL
    }
  )
  if (is.null(kegg_res) || nrow(as.data.frame(kegg_res)) == 0) {
    warning(paste("3E", group_name, "has no KEGG terms passing the cutoffs."), call. = FALSE)
    next
  }
  kegg_res <- setReadable(kegg_res, OrgDb = org.Ss.eg.db, keyType = "ENTREZID")
  kegg_df <- as.data.frame(kegg_res)
  write.csv(
    kegg_df,
    file.path(output_dir, paste0("Figure_3E_4C_vs_2C_", group_name, "_KEGG_Results.csv")),
    row.names = FALSE
  )
  save_enrichment_chord(
    kegg_df, current$deg_df, "3E", group_name, kegg_top_terms, "KEGG pathways"
  )
}

# Figure 3F ------------------------------------------
if (!file.exists(target_excel_path)) {
  stop(paste("Marker table not found:", target_excel_path))
}
table_markers <- read_excel(target_excel_path)
target_clusters <- c("4_cell", "8_cell")

cleaned_markers <- table_markers %>%
  dplyr::filter(cluster %in% target_clusters) %>%
  dplyr::filter(p_val_adj < 0.05 & avg_log2FC > 0.25) %>%
  dplyr::filter(!grepl("^MT-|^MRP|^RPL|^RPS", gene, ignore.case = TRUE)) %>%
  dplyr::filter(!toupper(gene) %in% c(
    "COX1", "COX2", "COX3", "ATP6", "ATP8", "ND1", "ND2",
    "ND3", "ND4", "ND4L", "ND5", "ND6", "CYTB"
  ))

ZGA_Active_genes <- cleaned_markers %>% pull(gene) %>% toupper() %>% unique()
if (length(ZGA_Active_genes) < 5) {
  stop("Fewer than five genes passed the enrichment filters.")
}

gene_entrez <- bitr(
  ZGA_Active_genes,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Ss.eg.db
)
if (nrow(gene_entrez) == 0) {
  stop("No pig Entrez IDs matched the supplied gene symbols.")
}

ego <- tryCatch({
  enrichGO(
    gene = gene_entrez$ENTREZID,
    OrgDb = org.Ss.eg.db,
    ont = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.1,
    qvalueCutoff = 0.1,
    readable = TRUE
  )
}, error = function(e) {
  warning(paste("GO enrichment failed:", conditionMessage(e)), call. = FALSE)
  NULL
})

if (!is.null(ego) && nrow(as.data.frame(ego)) > 0) {
  go_df <- as.data.frame(ego)
  write.csv(
    go_df,
    file.path(output_dir, "Figure_3F_Combined_4C_8C_GO_Results.csv"),
    row.names = FALSE
  )
  go_chord_genes <- data.frame(
    Gene = toupper(cleaned_markers$gene),
    avg_log2FC = cleaned_markers$avg_log2FC
  )
  go_symbols <- unique(unlist(strsplit(as.character(go_df$geneID), "/", fixed = TRUE)))
  go_chord_genes$Gene <- go_symbols[match(go_chord_genes$Gene, toupper(go_symbols))]
  save_enrichment_chord(
    go_df, go_chord_genes, "3F", "GO", go_top_terms, "GO biological processes",
    analysis_id = "Combined_4C_8C",
    analysis_label = "Combined 4-cell and 8-cell markers"
  )
} else {
  warning("3F has no GO terms passing the cutoffs.", call. = FALSE)
}

# Figure 3G ------------------------------------------
pbmc$Stage_Cleaned <- pbmc$Stage
Idents(pbmc) <- "Stage_Cleaned"

target_genes_ordered <- c(
  "PHF19", "CLOCK", "BSX", "EHF", "OTX2", "GSX2", "NFIB", "NR1H4", "CREB3L3",
  "TET1", "NME2", "TGIF1", "TFAP2C", "SPIC", "AEBP2", "GATA6", "MYCN",
  "XBP1", "BACH2", "ATF3", "MYC", "ZNF114", "KLF4"
)
real_genes <- target_genes_ordered[target_genes_ordered %in% rownames(pbmc)]

avg_exp <- AverageExpression(pbmc, features = real_genes, assays = "RNA", slot = "data")$RNA
mat_raw <- as.matrix(avg_exp)

colnames(mat_raw) <- gsub("MII_Oocyte|MII.Oocyte|MII-Oocyte", "MII-Oocyte", colnames(mat_raw), ignore.case = TRUE)
colnames(mat_raw) <- gsub("g?1.*cell|1C", "1-cell", colnames(mat_raw), ignore.case = TRUE)
colnames(mat_raw) <- gsub("g?2.*cell|2C", "2-cell", colnames(mat_raw), ignore.case = TRUE)
colnames(mat_raw) <- gsub("g?4.*cell|4C", "4-cell", colnames(mat_raw), ignore.case = TRUE)
colnames(mat_raw) <- gsub("g?8.*cell|8C", "8-cell", colnames(mat_raw), ignore.case = TRUE)

desired_stage_order <- c("1-cell", "2-cell", "4-cell", "8-cell")

available_stages <- desired_stage_order[desired_stage_order %in% colnames(mat_raw)]
mat_ordered <- mat_raw[, available_stages, drop = FALSE]

clean_rows <- apply(mat_ordered, 1, function(x) { !any(is.na(x)) && (sd(x) > 0) })
mat_ordered <- mat_ordered[clean_rows, , drop = FALSE]
mat_ordered <- mat_ordered[match(real_genes, rownames(mat_ordered)), , drop = FALSE]
mat_ordered <- na.omit(mat_ordered)

mat_scaled <- t(apply(mat_ordered, 1, function(x) { (x - mean(x)) / sd(x) }))

col_fun = colorRamp2(c(-1.5, 0, 1.5), c("#808080", "#FFFFFF", "#FF0000"))

n_cols <- ncol(mat_scaled)
n_rows <- nrow(mat_scaled)

cell_width_mm <- 4.5
cell_height_mm <- 4.0

ht <- Heatmap(
  mat_scaled,
  name = "Expression",
  col = col_fun,

  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_dend = FALSE,

  width = unit(n_cols * cell_width_mm, "mm"),
  height = unit(n_rows * cell_height_mm, "mm"),

  row_names_gp = gpar(fontface = "plain", size = 8.5),

  column_names_rot = 45,
  column_names_side = "bottom",
  column_names_gp = gpar(fontface = "plain", size = 9, hjust = 1, vjust = 1),

  rect_gp = gpar(col = "white", lwd = 1.0),

  heatmap_legend_param = list(
    direction = "vertical",
    legend_height = unit(3.0, "cm"),
    grid_width = unit(0.35, "cm"),
    title_gp = gpar(fontface = "plain", size = 9),
    labels_gp = gpar(size = 8),
    at = c(-1.5, 0, 1.5),
    labels = c("-1.5", "0.0", "1.5")
  )
)

pdf_width <- 3.9
pdf_height <- 6.5

cairo_pdf(file.path(output_dir, "Figure_3G_Heatmap.pdf"), width = pdf_width, height = pdf_height)
draw(ht, heatmap_legend_side = "right", padding = unit(c(5, 5, 12, 8), "mm"))
dev.off()

png(file.path(output_dir, "Figure_3G_Heatmap.png"), width = pdf_width, height = pdf_height, units = "in", res = 300)
draw(ht, heatmap_legend_side = "right", padding = unit(c(5, 5, 12, 8), "mm"))
dev.off()
