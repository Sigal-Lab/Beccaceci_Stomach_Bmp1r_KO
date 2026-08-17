require(MatrixGenerics)
dotplot_scRNA_multi_feature <- function(data, cluster, features, 
                                        group=NULL, 
                                        assay="RNA", 
                                        scale=F, 
                                        group_order=NULL, 
                                        cluster_order = NULL, 
                                        dot.max=6, 
                                        percentile_scale_upper=0.98, 
                                        min_val_scale = -2.5, 
                                        max_val_scale = 2.5,
                                        log_scale = F) {
  
  valid_features <- features[features %in% Features(data, assay=assay)]
  
  if(length(valid_features)==0) stop("None of the features found in selected assay.")
  
  cluster_orig_names <- unique(unique(data@meta.data[[cluster]]))
  names(cluster_orig_names) <- paste0("C",as.numeric(factor(cluster_orig_names)))
  group_orig_names <- unique(data@meta.data[[group]])
  names(group_orig_names) <- paste0("G",as.numeric(factor(group_orig_names)))
  
  data@meta.data$cluster_fixed_name <- names(cluster_orig_names)[match(data@meta.data[[cluster]], cluster_orig_names)]
  data@meta.data$group_fixed_name <- names(group_orig_names)[match(data@meta.data[[group]], group_orig_names)]

  dd <- LayerData(data, layer="data", assay=assay, features = valid_features)

  a = AverageExpression(data, assays = assay, features = valid_features, group.by = c("cluster_fixed_name", "group_fixed_name"), layer = "data")[[assay]]

  # Seurat::DotPlot() scales averaged values, not raw values. Even if this might not be ideal, let's keep this to make our plot
  # compatible with the Seurat one. Of note, Seurat::AverageExpression applies expm1() to all values before averaging, moving values from the log scale back to linear

  a_df = as.data.frame(t(as.matrix(a)))
  tmp = strsplit(rownames(a_df), "_")
  a_df$g1 = sapply(tmp, `[`,1)
  a_df$g2 = sapply(tmp, `[`,2)
  a_df_ts <- suppressMessages(reshape2::melt(a_df, id_vars = c("g1","g2"), variable.name="gene", value.name = "avg_expression"))
  a_df_ts$id <- paste0(a_df_ts$g1, "__", a_df_ts$g2, "__", a_df_ts$gene)
  
  tmp = apply(dd>0, 1, function(x) tapply(x, list(data@meta.data[["cluster_fixed_name"]], data@meta.data[["group_fixed_name"]]), sum, na.rm=T), simplify = F)
  cnt_expr <- mapply(function(x,n) {d <- reshape2::melt(x, value.name="num_detected"); d$gene <- n; return(d)}, tmp, names(tmp), SIMPLIFY = F) |> do.call(rbind, args=_)
  cnt_expr$num_detected <- ifelse(is.na(cnt_expr$num_detected),  0, cnt_expr$num_detected)
  cnt_expr$id <- paste0(cnt_expr$Var1, "__", cnt_expr$Var2, "__", cnt_expr$gene)
  cnt_expr$id2 <- paste0(cnt_expr$Var1, "__", cnt_expr$Var2)
  
  total_cnts = table(data@meta.data[["cluster_fixed_name"]], data@meta.data[["group_fixed_name"]]) |> as.data.frame()
  total_cnts$Freq <- ifelse(total_cnts$Freq == 0, 1, total_cnts$Freq) # avoid div by 0 problems
  total_cnts$id <- paste0(total_cnts$Var1, "__", total_cnts$Var2)
  colnames(total_cnts) = gsub("Freq","Total",colnames(total_cnts))  
  cnt_expr <- merge(cnt_expr, total_cnts[, c("id","Total")], by.x="id2", by.y="id", all.x=T, sort=F)
  cnt_expr$pct_expr <- cnt_expr$num_detected / cnt_expr$Total
  
  a_df <- merge(a_df_ts, cnt_expr[, c("id", "pct_expr")], by="id", all.x=T, sort=F)
  a_df$cluster_orig <- cluster_orig_names[a_df$g1]
  a_df$group_orig <- group_orig_names[a_df$g2]
  a_df$circle_size=100*a_df$pct_expr
  #y_label_fun = function(x) levels(data$pathway)[x]
  a_df$expr_scaled = split(a_df$avg_expression, a_df$gene) |> lapply(scale) |> unsplit(f=a_df$gene)
  
  expr_col = ifelse(scale, "expr_scaled", "avg_expression")
  #expr_col = "avg_expression"
  # min_scaled = min(a_df$expr_scaled, na.rm=T)
  # max_scaled = max(a_df$expr_scaled, na.rm=T)
  
  expr_label = ifelse(scale, "Expression (scaled)", "Expression")
  
  if(!is.null(group_order)) {
    a_df$group_orig <- factor(a_df$group_orig, levels=group_order)
  }
  
  if(!is.null(cluster_order)) {
    a_df$cluster_orig <- factor(a_df$cluster_orig, levels = cluster_order)
  }
  
  a_df$gene = factor(a_df$gene, levels=valid_features)
  
  p = ggplot(a_df, aes(x=group_orig, y= cluster_orig)) + 
    geom_point(aes(size=circle_size/10, color=.data[[expr_col]])) + 
    facet_grid(~ gene, switch="x") +
    theme_classic() + 
    theme(axis.text.x = element_text(angle=45, hjust=1, size=12 ),
          axis.text.y = element_text(size=14),
          plot.title = element_text(size=20)) + 
    labs(title="") + ylab(cluster) + xlab(group) + 
    scale_size_area(
      #range=c(.01, 6),
      max_size = dot.max,
      limits=c(0,max(2.5, max(a_df$circle_size/10)) ),
      breaks = seq(0,10,2.5),
      labels = paste(c(0,25,50,75,100), "%")) + 
    guides(size = guide_legend(title = "% detected"), color = guide_colorbar(title=expr_label )) #+ 
    #scale_x_discrete(drop = FALSE) + scale_y_discrete(drop = FALSE) # maintain subgroups without expression
  
  max_val = quantile(a_df[[expr_col]], percentile_scale_upper, na.rm=T) |> round(2)
  min_val = quantile(a_df[[expr_col]], 1-percentile_scale_upper, na.rm=T) |> round(2)
  
  if(scale) { 
    min_val <- min_val_scale
    max_val <- max_val_scale
    p = p + scale_colour_gradient2(low="blue", mid="grey95", high="red", midpoint = 0, limits=c(min_val,max_val), breaks=seq(min_val,max_val, by=max_val/4), oob=scales::squish )
  } else {
    p = p + scale_color_distiller(palette="Reds",direction = 1, limits=c(0,max_val), breaks=seq(0, max_val, by=max_val/4), oob=scales::squish )
  }
  
  return(p)
}

