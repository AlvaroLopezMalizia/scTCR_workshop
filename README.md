# scTCR_workshop

## Disclaimer 🚨

This script is designed solely for educational purposes. It does not provide an exhaustive analysis of any dataset, and the parameters, thresholds, and procedures outlined should not be directly applied to research tasks without careful consideration. In real-world scenarios, every decision made should be thoroughly evaluated, drawing upon deep knowledge of how the data was collected and the specific questions being addressed. Additionally, while this script reflects the author's preferred methods, it may contain areas for optimization or personalization to suit individual needs and preferences.

## Introduction to the Workshop 📚

At the top-right corner of the script windows, adjacent to the "Source" dropdown menu, you will find the outline button. Clicking it will reveal a clickable index to various parts of the script. It's recommended to organize scripts into sections and subsections using hashtags. For example, a section can be created with a single hashtag, and subsections can be nested with additional hashtags.

## Script Setup 🛠️

Working Directory Setup: The script automatically sets the working directory to the location of the script file using the setwd function. If needed, the directory can be manually set.

Libraries: Several libraries are loaded, including Seurat, dplyr, tidyr, clustree, patchwork, ggplot2, DescTools, divo, ggrepel, scGate, stringdist, and pheatmap, each serving various purposes in data analysis and visualization.

Random Number Generator Seed: The seed for the random number generator is set to ensure reproducibility of results.

02) Load Data 📊
RNAseq Data Loading: RNAseq counts, cell names, and gene names are loaded from cellranger output using the Read10X function from Seurat package.

TCR Data Loading: TCR data is loaded from a CSV file. Only full-length and productive contigs with TRA or TRB chains are selected for further analysis.

Processing TCR Data: TCR data is processed to include V gene information in strings along with sequence information. Clonotypes are defined based on this information. Droplets containing multiple TCR chains are identified as doublets.

Merge TCR Data with Seurat Object: TCR data is merged with the Seurat object, ensuring compatibility and preventing mixing of TCR information between samples with repeated barcodes.
