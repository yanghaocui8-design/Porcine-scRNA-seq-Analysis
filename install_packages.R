# R packages required for porcine scRNA-seq analysis

required_packages <- c(
  "Seurat",
  "harmony",
  "monocle3",
  "Mfuzz",
  "ggplot2",
  "patchwork",
  "pheatmap",
  "ComplexHeatmap",
  "circlize",
  "biomaRt",
  "clusterProfiler"
)

install_if_missing <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    message(paste("Installing package:", pkg))
    if (pkg %in% c("Seurat", "harmony", "ggplot2", "patchwork", "pheatmap", "circlize")) {
      install.packages(pkg)
    } else {
      if (!requireNamespace("BiocManager", quietly = TRUE)) {
        install.packages("BiocManager")
      }
      BiocManager::install(pkg)
    }
  }
}

sapply(required_packages, install_if_missing)
