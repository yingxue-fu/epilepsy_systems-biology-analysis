library(org.Hs.eg.db)
library(dplyr)
library(biomaRt)

# GB compound-target interactions
GBL_compound_target0 = read.csv("08 drug target interactions/GB_C_T_0.65.csv")
# get target entrez ID
load("03 brain cell type enrichment/human_mouse_biomaRt.RData")
GBLtarget_entrezID = getBM(c("hgnc_symbol", "entrezgene"), filters = "hgnc_symbol", values = GBL_compound_target0$Target, mart = human)
# check duplicatin
sum(duplicated(GBLtarget_entrezID$hgnc_symbol))
# GBL targets in interactome
GBLtarget_interome = GBLtarget_entrezID[GBLtarget_entrezID$entrezgene %in% ppi_v$name, ]
# compound target lists
GBL_compound_target = merge(GBL_compound_target0, GBLtarget_interome, by.x = "Target", by.y = "hgnc_symbol")
GBL_compound_target$entrezgene = as.character(GBL_compound_target$entrezgene)
GBL_compound_target_list = split(GBL_compound_target$entrezgene, GBL_compound_target$Compound)

save(GBLtarget_entrezID, GBLtarget_interome, GBL_compound_target, GBL_compound_target_list, file = "08 drug target interactions/GBL_compound_target.RData")
