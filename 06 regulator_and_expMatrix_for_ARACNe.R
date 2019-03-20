library(biomaRt)

load("03 brain cell type enrichment/human_mouse_biomaRt.RData")

# Transcription Factors
TFs <- getBM(c("entrezgene", "mgi_symbol"), filters = "go", values = c("GO:0003700", "GO:0003677", "GO:0140110"), mart = mouse)
TFs = na.omit(TFs) # 1844

# Membrane Receptors
# GO:0045202_synapse_Cellular Component
# GO:0004888_transmembrane signaling receptor activity_Molecular Function
MRs <- getBM(c("entrezgene", "mgi_symbol"), filters = "go", values = c("GO:0045202", "GO:0004888"), mart = mouse)
MRs = na.omit(MRs) # 663

# regulators
TFs$type = "TF"
MRs$type = "MR"
regulators = rbind(TFs, MRs) # 2507

save(regulators, file = "06 master regulator identification/regulators_fromGO.RData")

# relate regulators to probeIDs
load("02 gene coexpression network analysis/probe_modules_anno.RData")
regulators_anno = merge(regulators, probe_modules_anno, by.x = "entrezgene", by.y = "Gene_ID") # 1507

save(regulators_anno, file = "06 master regulator identification/regulators_anno.RData")

# output the regulon probe ID data for ARACNe software
regulon_probe_ID = as.data.frame(regulators_anno$probe_ID)
write.table(regulon_probe_ID, file = "06 master regulator identification/regulon.txt", quote = F, row.names = F, col.names = F)

# output the expression matrix data for ARACNe software
load("01 sample and gene filtering/samples_and_probes_filtering.RData")
write.csv(datExpr, file = "06 master regulator identification/expMatrix.csv", quote = F)
# save it as a tab-seperated .txt file using Excel 
