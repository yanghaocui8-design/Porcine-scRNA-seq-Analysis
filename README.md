# Porcine Single-Cell RNA-seq Data Analysis

<p align="left">
  <img src="https://img.shields.io/badge/R-%3E%3D_4.3.0-blue.svg" alt="R Version">
  <img src="https://img.shields.io/badge/License-MIT-yellow.svg" alt="License">
  <img src="https://img.shields.io/badge/SingleCell-Seurat_v5-brightgreen.svg" alt="Seurat">
</p>

This repository contains reproducible R scripts and downstream analysis pipelines for the porcine single-cell RNA-seq (scRNA-seq) dataset, organized sequentially by manuscript figures and panels.

---

## Abstract

Here, we present a comprehensive single-cell transcriptomic analysis of porcine pre-implantation embryonic development. By integrating multi-source scRNA-seq datasets, applying Harmony for batch correction, and utilizing advanced bioinformatics tools (including Monocle3 trajectory inference and Mfuzz temporal clustering), we uncover the key transcriptional dynamics, ZGA (zygotic genome activation) programs, and lineage-specific molecular signatures (ICM, TE, EPI, and PE) during early embryogenesis.

---

## 📁 Repository Structure

```text
.
├── README.md                  # Project overview and documentation
├── install_packages.R         # Script for checking and installing dependencies
├── code/                      # Main analysis scripts organized by modules/figures
│   ├── 01_preprocessing/      # Data integration, QC, and batch correction
│   ├── Fig1/                  # QC evaluation, gene metrics, and transcriptome correlation
│   ├── Fig2/                  # UMAP clustering, marker bubble plots, and feature heatmaps
│   ├── Fig3/                  # DEG statistics, ZGA transition, and transcription factor dynamics
│   ├── Fig4/                  # Monocle3 pseudotime trajectory and lineage dynamics
│   └── Fig5/                  # Mfuzz expression modules (C1-C5) and functional enrichment
└── data/                      # Sample metadata and downsampled objects (git-ignored large files)

## 🛠️ Software & Environment

- **Language**: R (v4.3.0)
- **Core Packages & Versions**:
  - `Seurat` (v5.0.0)
  - `Harmony` (v1.2.0)
  - `Monocle3` (v1.3.1)
  - `Mfuzz` (v2.60.0)
  - `pheatmap` (v1.0.12)
  - `biomaRt` (v2.58.0)
  - `clusterProfiler`, `ComplexHeatmap`, `circlize`, `ggplot2`
