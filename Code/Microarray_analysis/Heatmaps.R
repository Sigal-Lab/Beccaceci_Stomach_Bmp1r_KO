setwd("S:/AG/AG-Sigal/Research/Documentation/Projects/Hilmar Berger/Analysen/Giulia_Hpy_antrum_corpus/MicroArray_inf_uninf_antrum_corpus")

d_antrum <- read.table("./Results/Inf_vs_Uninf/FC_average_across_dyeswap.txt", sep="\t", header=T)

m <- as.matrix(d_antrum[, "infected_antrum_vs_uninfected_antrum" ])
rownames(m) <- d_antrum$GeneSymbol
m_avg <- limma::avereps(m, rownames(m))

sel_genes <- strsplit("Cxcl2,Cxcl10,S100A8,Cxcl6,Ccl20,Il1b, S100A9,Il1a, Il6,Ccl2,C4b,C3", split=",") |> unlist()

library(pheatmap)
pheatmap(m_avg[rownames(m_avg) %in% sel_genes,,drop=F], cluster_rows = T, cluster_cols=F)


setwd("S:/AG/AG-Sigal/Research/Documentation/Projects/Hilmar Berger/Analysen/Giulia_Col1a1_Il1r/")

d_il1rko <- read.table("./Results/FC_average_across_dyeswap.txt", sep="\t", header=T)

m <- as.matrix(d_il1rko[, 1:2])
rownames(m) <- d_il1rko$GeneSymbol
m_avg <- limma::avereps(m, rownames(m))

sel_genes <- strsplit("Cxcl2,Cxcl10,S100A8,Cxcl6,Ccl20,Il1b, S100A9,Il1a, Il6,Ccl2,C4b,C3", split=",") |> unlist()

pheatmap(m_avg[rownames(m_avg) %in% sel_genes,,drop=F], cluster_rows = T, cluster_cols=F)



