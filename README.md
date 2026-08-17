# Companion code repository for publication "Helicobacter pylori triggers gastric mucosal remodeling toward a fetal-like transcriptional program via stromal IL-1β signaling" by Beccaceci et al.

This repository contains R-scripts used for analysis of microarray and scRNA data from the manuscript mentioned above. 

R scripts were run using R version R-4.3.2. 
The following packages must be installed: 

```
install.packages(c("cowplot", "data.table", "dplyr", "future", "ggplot2", "grid", "harmony", "knitr",
  "magrittr", "naturalsort", "pheatmap", "plotly", "png", "readxl", "writexl", "BiocManager"))
BiocManager::install(c("BiocParallel","clusterProfiler","EnhancedVolcano","GSEABase","homologene",
  "limma","msigdbr","org.Hs.eg.db","org.Mm.eg.db","ReactomePA","scater","scDblFinder","Seurat",
  "speckle","SummarizedExperiment","UCell"))
```

To reproduce the results, follow the instructions in Data/Raw/README_RAW_MICROARRAY_FILES and Data/scRNA/CellRanger/README to download data from Gene Expression Omnibus first.

Then run the shell scripts Code/Microarray_analysis/run_all.sh and Code/scRNA/run_all.sh.
