# Cell Annotation by PCA Projection and Cell Label Transfer

This repository provides an R workflow for **annotating cell types in single‑cell RNA‑seq data** by projecting query cells into the PCA space of a reference dataset and transferring labels based on nearest neighbors. It demonstrates an implementation alongside Seurat’s built‑in label transfer functions.

> The workflow uses **PCA projection**, **nearest neighbor search**, and **majority voting or Seurat label transfer** to assign cell type labels from a reference (large) dataset to a smaller query dataset.

---

## 📋 Overview

Single‑cell RNA sequencing reveals cellular heterogeneity across tissues, but assigning biologically meaningful **cell type labels** to clusters of cells requires reference data or marker information. This project illustrates a method where:

1. A larger **reference dataset** (e.g., 10k PBMC cells) is processed and reduced via PCA.
2. A smaller **query dataset** (e.g., 3k PBMC cells) is normalized and projected into the reference PCA space.
3. Nearest neighbors in the PCA space are used to **transfer labels** from the reference to the query.
4. Optional: Seurat’s `FindTransferAnchors` and `TransferData` functions perform label transfer.

---

## 🧠 Key Concepts

- **PCA Projection** — Project query cells into reference PCA space to compare similar expression profiles. :contentReference[oaicite:1]{index=1}  
- **Nearest Neighbors** — Identify closest cells in reference PCA space to assign labels.  
- **Majority Vote** — Labels from nearest neighbors are aggregated to assign a most likely cell type.  
- **Seurat Label Transfer** — Optional built‑in Seurat method that uses anchors for transfer.

---

## 🧰 Dependencies

This analysis uses R and the following packages:

```r
library(Seurat)
library(Matrix)
library(irlba)        # fast PCA
library(RcppAnnoy)    # fast nearest neighbor search
library(dplyr)
library(SeuratData)   # demo datasets
library(ggplot2)
library(ComplexHeatmap)
