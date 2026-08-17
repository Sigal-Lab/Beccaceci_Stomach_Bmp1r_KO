require(org.Hs.eg.db)
require(org.Mm.eg.db)
require(data.table)
require(dplyr)

enrichGSEA <- function(genes, collection, organism, universe) {
  m_t2g <- msigdbr(species = organism, category = collection)[, c("gs_name", "entrez_gene")]
  m_t2g$gs_name = ifelse(nchar(m_t2g$gs_name)>80, paste0(substr(m_t2g$gs_name,1,80),"~"), m_t2g$gs_name)
  enricher(genes, TERM2GENE=m_t2g, universe = universe) 
}

enrichMAFRES <- function(genes, universe) {
  m_t2g <- maf_ref_top_lists_gs
  m_t2g$gs_name = ifelse(nchar(m_t2g$gs_name)>80, paste0(substr(m_t2g$gs_name,1,80),"~"), m_t2g$gs_name)
  enricher(genes, TERM2GENE=m_t2g, universe=universe) 
}

# enrichGSEA <- function(genes, collection, organism) {
#   m_t2g <- msigdbr(species = organism, category = collection)[, c("gs_name", "entrez_gene")]
#   m_t2g$gs_name = ifelse(nchar(m_t2g$gs_name)>80, paste0(substr(m_t2g$gs_name,1,80),"~"), m_t2g$gs_name)
#   enricher(genes, TERM2GENE=m_t2g) 
# }

enrichCustomSigs <- function(genes, universe, genesets) {
  m_t2g <- genesets
  m_t2g$gs_name = ifelse(nchar(m_t2g$gs_name)>80, paste0(substr(m_t2g$gs_name,1,80),"~"), m_t2g$gs_name)
  enricher(genes, TERM2GENE=m_t2g, universe=universe)
}

# Seurat::FindMarkers and Seurat::FindAllMarkers use a min.pct=0.01 parameter to pre-filter genes that will be tested. 
# This means that the universe of genes is restricted to those detected in at least 1% of samples in at least one group.
# Here we compute the proportion of cells with detection for each gene and each cluster and then define the universe for enrichment accordingly.
get_FindMarkers_universe <- function(scData, cluster_col = "cluster_label", min_pct = 0.01) {
  dcg <- list()
  for (cc in unique(scData@meta.data[[cluster_col]])) {
    tmp <- LayerData(scData, assay="RNA", layer="data")[, scData@meta.data[[cluster_col]] == cc]
    det_counts <- apply(tmp>0, 1, sum)
    det_count_prop <- det_counts / ncol(tmp)
    dcg[[cc]] <- det_count_prop
  }
  dcg_mat <- do.call(cbind, dcg)
  max_dcp <- apply(dcg_mat, 1, max)
  #hist(max_dcp, 100)
  dge_universe = names(max_dcp[max_dcp > min_pct])
  return(dge_universe)
}

fracstr_to_value <- Vectorize(function(x) eval(parse(text=x)), vectorize.args = "x") 

# This function takes a table containing DE calls for one to several comparisons and 
# computes functional enrichment using ClusterProfiler (i.e. overrepresentation analysis)
# Parameters:
# sc_data   The full Seurat object. Used to derive the universe if not provided specifically.
# curr_dge_tab  data frame of prefiltered significant hits containing one row per gene and comparison
# log2FC_col  name of the log2FC column (used to define up/down direction)
# gene_col  name of the gene column
# celltype_col  name of the column defining celltypes or comparisons 
# custom_sets_Hs  list of gene sets with human gene symbols
# custom_sets_Mm  list of gene sets with mouse gene symbols
compute_functional_enrichment <- function(sc_data, 
                                          curr_dge_tab, 
                                          universe=NULL, 
                                          gene_col = "GeneSymbol", 
                                          species = "Hs", 
                                          log2FC_col = "avg_log2FC", 
                                          celltype_col = "celltype",
                                          custom_sets_Hs = NULL,
                                          custom_sets_Mm = NULL) {

  if(! species %in% c("Hs","Mm")) stop("Only Hs and Mm supported for species.")
  
  # assign up/down status to DE results
  curr_dge_tab$direction = ifelse(curr_dge_tab[[log2FC_col]] > 0, "Up", ifelse(curr_dge_tab[[log2FC_col]] < 0, "Down","nc"))
  curr_dge_tab$gene_group = paste0(curr_dge_tab[[celltype_col]],"_",curr_dge_tab$direction)
  
  if(species == "Hs") {
    tmp = AnnotationDbi::select(org.Hs.eg.db, keys=rownames(sc_data@assays$RNA), keytype = "SYMBOL", columns="ENTREZID")  
  } else {
    tmp = AnnotationDbi::select(org.Mm.eg.db, keys=rownames(sc_data@assays$RNA), keytype = "SYMBOL", columns="ENTREZID")  
  }
  
  
  gs2e = tapply(tmp$ENTREZID, tmp$SYMBOL, paste, collapse=",")
  e2gs = tapply(tmp$SYMBOL, tmp$ENTREZID, paste, collapse=",")
  
  genes_by_class = split(gs2e[curr_dge_tab[[gene_col]]], factor(curr_dge_tab$gene_group))
  genes_by_class = lapply(genes_by_class, function(x) unique(x[!is.na(x)]))
  
  if(is.null(universe)) {
    universe = gs2e[Features(sc_data, layer = "data")]
  } else {
    universe = gs2e[universe]
  }
  universe = universe[!is.na(universe)]
  
  all_enrichments = new.env() # futures can only be assigned to an environment, not a list
  
  if(species == "Hs") {
    org_db = "org.Hs.eg.db"
    org_reactome = "human"
    org_name = "Homo sapiens"
  } else {
    org_db = "org.Mm.eg.db"
    org_reactome = "mouse"
    org_name = "Mus musculus"
  }
  
  custom_sets <- list()
  if(!is.null(custom_sets_Hs)) {
    if(species == "Hs") {
      custom_sets <- append(custom_sets, custom_sets_Hs)
    } else {
      custom_sets <- append(custom_sets, lapply(custom_sets_Hs, function(x) human2mouse(x, db=homologeneData2) |> (\(x) {unique(x$mouseGene[!is.na(x$mouseGene)])} )()) )
    }
  }
  if(!is.null(custom_sets_Mm)) {
    if(species == "Mm") {
      custom_sets <- append(custom_sets, custom_sets_Mm)
    } else {
      custom_sets <- append(custom_sets, lapply(custom_sets_Mm, function(x) mouse2human(x, db=homologeneData2) |> (\(x) {unique(x$humanGene[!is.na(x$mouseGene)])} )()) )
    }
  }
  
  
  
  # future has problems identifying the enrich... functions from preexisting packages as global variabales (at least on Windows) - lets define the full path for them
  all_enrichments[["GO_BP"]] %<-% compareCluster(geneCluster = genes_by_class, fun = enrichGO, OrgDb=org_db, ont="BP", universe = universe)
  #all_enrichments[["KEGG"]] %<-% compareCluster(geneCluster = genes_by_class, fun = enrichKEGG, organism = "hsa", use_internal_data=T)
  all_enrichments[["REACTOME"]] %<-% compareCluster(geneCluster = genes_by_class, fun = enrichPathway, organism = org_reactome, universe = universe)
  all_enrichments[["GSEA_H"]] %<-% compareCluster(geneCluster = genes_by_class, fun = enrichGSEA, organism = org_name, collection="H", universe = universe)
  all_enrichments[["GSEA_C2"]] %<-% compareCluster(geneCluster = genes_by_class, fun = enrichGSEA, organism = org_name, collection="C2", universe = universe)
  all_enrichments[["GSEA_C3"]] %<-% compareCluster(geneCluster = genes_by_class, fun = enrichGSEA, organism = org_name, collection="C3", universe = universe)
  all_enrichments[["GSEA_C8"]] %<-% compareCluster(geneCluster = genes_by_class, fun = enrichGSEA, organism = org_name, collection="C8", universe = universe)

  if(length(custom_sets)>0) {
    tmp <- lapply(custom_sets, function(x) gs2e[x]) |> lapply(function(x) x[!is.na(x) & x!=""])
    sig_tab <- data.frame(gs_name = rep(names(tmp), times=sapply(tmp, length)), 
                          entrez_gene = unlist(tmp))
    all_enrichments[["CUSTOM_SETS"]] %<-% compareCluster(geneCluster = genes_by_class, fun = enrichCustomSigs, universe = universe, genesets = sig_tab)
  }
  
  all_enrichments = as.list(all_enrichments)
  
  # write to table
  for (n in names(all_enrichments)) {
    tmp = all_enrichments[[n]]
    tmp@compareClusterResult$database = n
    
    p1 = fracstr_to_value(tmp@compareClusterResult$GeneRatio)
    p2 = fracstr_to_value(tmp@compareClusterResult$BgRatio)
    
    tmp@compareClusterResult$OR = (p1/(1-p1))/(p2/(1-p2))
    tmp@compareClusterResult$log2OR = log2(tmp@compareClusterResult$OR)

    all_enrichments[[n]] = tmp
  }
  
  # build output table
  HM_cluster_enrichment = do.call(rbind, lapply(all_enrichments, function(x) x@compareClusterResult) )
  
  HM_cluster_enrichment$GeneSymbols = Map(function(x) {sapply(strsplit(x,"/"), function(y) sort(paste(e2gs[y], collapse=",")))}, HM_cluster_enrichment$geneID )
  cc = colnames(HM_cluster_enrichment)
  sel_cols = c("Cluster","database")
  sel_cols = append(sel_cols, cc[!cc %in% sel_cols & !cc =="geneID"])

  tmp <- group_split(HM_cluster_enrichment %>% group_by(Cluster))
  tmp2 <- lapply(tmp, function(x) x[order(x$database, x$pvalue), c("Cluster","database","Description","GeneRatio","BgRatio","pvalue","p.adjust","qvalue","OR","log2OR","Count","GeneSymbols")] )
  tmp2 <- lapply(tmp2, function(x) {x$GeneSymbols = unlist(lapply(x$GeneSymbols, paste, collapse=", ")) ; return(x) } )
  names(tmp2) = unlist(lapply(tmp2, function(x) unique(make.names(as.character(x$Cluster)))))
  
  return(list(clustprof_enrichments = all_enrichments, output_tables = tmp2 ) )
}

