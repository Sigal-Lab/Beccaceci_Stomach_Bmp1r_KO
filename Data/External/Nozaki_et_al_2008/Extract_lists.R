setwd("S:/AG/AG-Sigal/Research/Documentation/Projects/Hilmar Berger/DataSets/SPEM/Nozaki_et_al_2008")

library(readxl)

up = as.data.frame(read_excel("SPEM_Signature.xlsx", sheet = "Up"))
down = as.data.frame(read_excel("SPEM_Signature.xlsx", sheet = "Down"))

all_probes = unique(c(up$probe, down$probe))

library(mouse4302.db)
library(AnnotationDbi)
symbols = select(mouse4302.db, keys=all_probes, keytype="PROBEID", columns="SYMBOL")
symbols  = subset(symbols, !is.na(SYMBOL))
probe_to_symbol = tapply(symbols$SYMBOL, symbols$PROBEID, paste, collapse=",")

spem_sig = data.frame(probeid = c(up$probe, down$probe), direction = c(rep("Up", nrow(up)), rep("Down",nrow(down))) )
spem_sig$GeneSymbol = probe_to_symbol[spem_sig$probeid]

write.table(spem_sig, file="SPEM_signature.txt", sep="\t", row.names=F, quote=F)
