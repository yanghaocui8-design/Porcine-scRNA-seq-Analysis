# Porcine Pre-implantation Embryo scRNA-seq Analysis

<p align="left">
  <img src="https://img.shields.io/badge/Language-R-blue.svg" alt="R">
  <img src="https://img.shields.io/badge/Analysis-scRNA--seq-brightgreen.svg" alt="scRNA-seq">
  <img src="https://img.shields.io/badge/Workflow-Figure--based-orange.svg" alt="Figure-based workflow">
</p>

## Overview

This repository contains R scripts for the downstream analysis and visualization of publicly available single-cell RNA sequencing (scRNA-seq) datasets from porcine pre-implantation embryos.

The analysis integrates seven public GEO datasets and follows the main figures of the associated manuscript. The workflow includes quality control, data integration, Harmony-based batch correction, developmental-stage characterization, differential gene expression analysis, transcriptional dynamics, pseudotime trajectory inference, and temporal expression clustering.

The repository is organized by analysis modules and manuscript figures to facilitate code navigation and reproducibility.

## Abstract

This study investigates transcriptional dynamics during porcine pre-implantation embryonic development using an integrated scRNA-seq dataset spanning oocytes, early cleavage-stage embryos, morulae, and blastocyst-associated cell states.

The analysis characterizes developmental-stage-associated expression patterns, examines transcriptional changes associated with zygotic genome activation (ZGA), and explores lineage-associated molecular features in the inner cell mass (ICM), trophectoderm (TE), epiblast (EPI), and primitive endoderm (PE). Monocle3 is used to infer transcriptional trajectories, while Mfuzz is applied to identify temporal expression patterns across developmental stages and cell states.

These analyses provide a transcriptomic framework for investigating developmental transitions and lineage-associated expression programs in porcine pre-implantation embryos.

## Repository Structure

```text
.
├── README.md
├── install_packages.R
├── code/
│   ├── 01_preprocessing/
│   │   └── # Data preparation, quality control, and integration
│   ├── Fig1/
│   │   └── # Quality control and transcriptome correlation
│   ├── Fig2/
│   │   └── # UMAP and stage-/lineage-associated marker expression
│   ├── Fig3/
│   │   └── # Differential expression and ZGA-associated dynamics
│   ├── Fig4/
│   │   └── # Monocle3 pseudotime and gene expression dynamics
│   └── Fig5/
│       └── # Mfuzz clustering and functional enrichment
└── data/
    └── # Metadata and required input objects, where applicable
```

The directory descriptions above summarize the intended organization. Actual script and data filenames should be added once the repository contents have been finalized.

Large datasets and intermediate objects are not necessarily included in this repository. Please refer to the data availability and input requirements documented below.

## Analysis Workflow

| Module        | Main analyses                                                                                                           |
| ------------- | ----------------------------------------------------------------------------------------------------------------------- |
| Preprocessing | Data preparation, quality control, metadata harmonization, and Harmony-based integration                                |
| Figure 1      | Quality control metrics and transcriptome correlation across developmental stages                                       |
| Figure 2      | UMAP visualization and expression patterns of stage- and lineage-associated genes                                       |
| Figure 3      | Differential gene expression, ZGA-associated transcriptional changes, and candidate regulatory-factor dynamics          |
| Figure 4      | Monocle3 pseudotime trajectory inference and gene expression dynamics                                                   |
| Figure 5      | Mfuzz temporal expression clustering, functional enrichment, and protein–protein interaction analysis, where applicable |

## Data Availability

The analysis uses seven publicly available porcine embryonic scRNA-seq datasets obtained from the Gene Expression Omnibus (GEO).

The corresponding GEO accession numbers, data sources, and input-file requirements should be documented in a dataset manifest or in a table in this section.

Raw sequencing data, processed expression matrices, cell annotations, and integrated analysis objects are distinct data types. The specific inputs required for each analysis module should be stated in the corresponding script or module-level README.

Large input files and intermediate objects may be excluded from version control. Users should obtain the required files from their original data repositories or other documented sources.

## Software and Environment

The analysis is implemented in R. The software environment and package versions should reflect the environment used to generate the reported results.

Key packages used in the workflow may include:

* `Seurat` — single-cell data processing and analysis
* `harmony` — batch correction and integration
* `monocle3` — pseudotime trajectory inference
* `Mfuzz` — temporal expression clustering
* `ggplot2` and `patchwork` — visualization and figure assembly
* `pheatmap` and `ComplexHeatmap` — heatmap visualization
* `circlize` — color mapping and circular visualization
* `biomaRt` — gene annotation, where applicable
* `clusterProfiler` — functional enrichment analysis

Additional packages may be required by individual scripts.

For reproducibility, the exact R version and package versions should be recorded from the analysis environment. A session information file can be generated using `sessionInfo()` after the analysis is completed.

## How to Use

1. Clone or download this repository.
2. Install the required R packages using the provided installation script, after verifying its package list.
3. Obtain the input datasets and metadata specified by the relevant analysis scripts.
4. Run the preprocessing and integration workflow before executing downstream analyses that depend on the integrated object.
5. Execute the figure-specific scripts in the order required by their input dependencies.

The exact execution order, required input files, and output locations should be documented in each module's README or at the beginning of its R scripts.

## Reproducibility Notes

* Analysis results depend on the input datasets, quality-control criteria, metadata annotations, and analysis parameters.
* Differential expression thresholds, clustering parameters, and functional enrichment settings should be explicitly recorded in the corresponding scripts.
* Pseudotime represents a computational ordering of transcriptional states and should not be interpreted as directly measured developmental time.
* Candidate transcription factors and stage-associated expression patterns are descriptive findings unless supported by additional functional evidence.
* Random seeds and software versions should be recorded where applicable.

## License

A license should be added after the intended terms of use have been confirmed. The license for the analysis code does not automatically determine the redistribution rights of third-party datasets.

## Citation

If you use this repository, please cite the associated manuscript once its publication details are available.

A formal citation can be added here after the manuscript has been accepted or published.
