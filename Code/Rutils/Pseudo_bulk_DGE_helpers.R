require(Seurat)
require(DESeq2)

# extract pseudo bulk count matrices for each cluster, summing over all cells in each sample
seurat_to_cell_type_matrices <- function(object, group_col, sample_id_col, min_sample_num = 2, min_cells_per_sample = 0, assay="RNA") {
  
  all_orig_groups <- unique(as.character(object@meta.data[[group_col]]))
  names(all_orig_groups) <- paste0("G", 1:length(all_orig_groups))
  all_orig_samples <- unique(as.character(object@meta.data[[sample_id_col]]))
  names(all_orig_samples) <- paste0("G", 1:length(all_orig_samples))
  
  object@meta.data$group_tmp <- names(all_orig_groups)[match(object@meta.data[[group_col]], all_orig_groups)]
  object@meta.data$sample_tmp <- names(all_orig_samples)[match(object@meta.data[[sample_id_col]], all_orig_samples)]
  
  object@meta.data[["cluster_and_sample"]] = paste0(object@meta.data$group_tmp, "...", object@meta.data$sample_tmp)
  cell_counts_per_cluster_and_sample <- table(object@meta.data[["cluster_and_sample"]])
  good_entries = names(cell_counts_per_cluster_and_sample[cell_counts_per_cluster_and_sample >= min_cells_per_sample])
  num_total = length(cell_counts_per_cluster_and_sample)
  num_excluded = num_total - length(good_entries)
  
  if(num_excluded > 0) {
    warning(paste0("Excluding ", num_excluded, "/", num_total, " entries with less than ", min_cells_per_sample, " cells per sample and condition.\n"))  
  }
  
  agg_mat <- AggregateExpression(object, assays=assay, slot="counts", group.by = "cluster_and_sample")[[assay]]
  colnames(agg_mat) <- as.character(colnames(agg_mat))
  agg_mat <- agg_mat[, good_entries]

  tmp <- strsplit(colnames(agg_mat), "...", fixed=T)
  cl <- sapply(tmp, `[`,1) 
  si <- sapply(tmp, `[`,2)
  
  all_groups <- unique(object@meta.data$group_tmp) |> sort()
  
  count_matrices <- list()
  for (cc in all_groups) {
    tmp <- agg_mat[,cl==cc, drop=F]
    if(ncol(tmp)<min_sample_num) {
      warning(paste0("Attempting to select count matrix with less than ", min_sample_num," samples for group ",cc,". Skipping.\n"))
      next
    }
    colnames(tmp) <- all_orig_samples[si[cl==cc]]
    count_matrices[[all_orig_groups[cc]]] <- tmp
  }
  
  return(count_matrices)
}


compute_diff <- function(count_matrices, control_cluster = "control") {
  all_target_clusters <- names(count_matrices)[!names(count_matrices) == control_cluster]
  sel_target_clusters <- all_target_clusters
  
  all_res <- list()
  
  for (tt in sel_target_clusters) {
    ctrl_mat <- count_matrices[[control_cluster]]
    colnames(ctrl_mat) <- paste(colnames(ctrl_mat),"__", control_cluster)
    trgt_mat <- count_matrices[[tt]]
    colnames(trgt_mat) <- paste(colnames(trgt_mat),"__", tt)
    
    cm <- cbind(ctrl_mat, trgt_mat)
    
    col_data <- data.frame(row.names=colnames(cm), group=factor(c(rep(control_cluster,ncol(ctrl_mat)), rep(tt,ncol(trgt_mat)) ), levels=c(control_cluster, tt)) )
    
    dds <- DESeqDataSetFromMatrix(cm, col_data, ~ group)
    dds <- DESeq(dds, fitType = "local")
    #print(resultsNames(dds))
    res <- results(dds, tidy = T)
    res$comparison = paste0(tt, "_vs_", control_cluster)
    res <- res[order(res$log2FoldChange, decreasing = T),]
    
    colnames(res) <- gsub("row","GeneSymbol", colnames(res))

    all_res[[tt]] <- res
  }
  
  return(all_res)
}


