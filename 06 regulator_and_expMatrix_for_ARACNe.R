library(biomaRt)
mouse_maRt <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")
list_of_filters <- listFilters(mouse_maRt)
list_of_attributes = listAttributes(mouse_maRt)
save(mouse_maRt, file = "06 inferring protein activity/mouse_maRt.RData")

# Transcription Factors GO:0003700∪(GO:0003677∩(GO:0140110∪GO:0006355))

# GO:0003700, ‘DNA binding transcription factor activity’
GO_0003700 <- getBM("mgi_symbol", filters = "go", values = "GO:0003700", mart = mouse_maRt)

# GO:0003677, ‘DNA binding’
GO_0003677 <- getBM("mgi_symbol", filters = "go", values = "GO:0003677", mart = mouse_maRt)

# transcription regulation genes 
# GO:0140110, 'transcription regulator activity'
GO_0140110 <- getBM("mgi_symbol", filters = "go", values = "GO:0140110", mart = mouse_maRt)
# GO:0003712, ‘transcription coregulator activity’
GO_0003712 <- getBM("mgi_symbol", filters = "go", values = "GO:0003712", mart = mouse_maRt)
# GO:0006355, 'regulation of transcription, DNA-templated'
GO_0006355 <- getBM("mgi_symbol", filters = "go", values = "GO:0006355", mart = mouse_maRt)

save(GO_0003677, GO_0003700, GO_0003712, GO_0006355, GO_0140110, file = "06 inferring protein activity/transcription-related GO.RData")

trans_regu <- union(union(GO_0140110$mgi_symbol, GO_0003712$mgi_symbol), GO_0006355$mgi_symbol)
# restrict to DNA-binding
trans_regu_DNA_binding <- intersect(trans_regu, GO_0003677$mgi_symbol)

TFs <- union(GO_0003700$mgi_symbol, trans_regu_DNA_binding)

# Synaptic proteins
# GO:0045202_synapse  # Cellular Component
GO_0045202 <- getBM("mgi_symbol", filters = "go", values = "GO:0045202", mart = mouse_maRt)
# GO:0030424 axon and GO:0030425 dendrite
GO_0030424 <- getBM("mgi_symbol", filters = "go", values = "GO:0030424", mart = mouse_maRt)
GO_0030425 <- getBM("mgi_symbol", filters = "go", values = "GO:0030425", mart = mouse_maRt)

save(GO_0045202,GO_0030424,GO_0030425, file = "06 inferring protein activity/synapse-related GO.RData")

SPs <- union(GO_0045202$mgi_symbol, union(GO_0030424$mgi_symbol, GO_0030425$mgi_symbol))

SPs <- SPs[-match(intersect(SPs, TFs), SPs)]

# Signaling proteins that are not overlaped with above regulators
# GO:0007165 ‘signal transduction’
GO_0007165 <- getBM("mgi_symbol", filters = "go", values = "GO:0007165", mart = mouse_maRt)
# overlapped genes
olap <- union(intersect(GO_0007165$mgi_symbol, SPs), intersect(GO_0007165$mgi_symbol, TFs))
# Final set
SignalP <- GO_0007165$mgi_symbol[-match(olap, GO_0007165$mgi_symbol)]

save(SPs, TFs, SignalP, file = "06 inferring protein activity/regulators_sep.RData")

# regulators dataframe
TFs_df = data.frame(Gene_symbol = TFs)
SPs_df = data.frame(Gene_symbol = SPs)
SignalP_df = data.frame(Gene_symbol = SignalP)

TFs_df$type = "TF"
SPs_df$type = "SP"
SignalP_df$type = "Signal"
regulators = rbind(TFs_df, SPs_df, SignalP_df) # 4564

sum(duplicated(regulators$Gene_symbol))

save(regulators, file = "06 inferring protein activity/regulators_df.RData")

# relate regulators to genes in the expression matrix
load("03 Gene coexpression network/gene_modules.RData")
regulators_module = merge(regulators, gene_modules, by = "Gene_symbol") # 1611 unique regulators

table(regulators_module$type)

save(regulators_module, file = "06 inferring protein activity_a/regulators_module.RData")

# output the regulon probe ID data for ARACNe software
write.table(regulators$Gene_symbol, file = "06 inferring protein activity_a/regulon.txt", quote = T, row.names = F, col.names = F)

# output the expression matrix data for ARACNe software
load("03 Gene coexpression network/datExpr.RData")
write.table(datExpr, file = "06 inferring protein activity_a/expMatrix.txt", quote = T, row.names = T, col.names = T, sep = "\t")
# save it as a tab-seperated .txt file using Excel 




################################################################################
target_GO <- getBM(c("mgi_symbol", "go_id", "name_1006"), filters = "mgi_symbol", values = c("Scgn", "Kirrel3", "Ntng1", "Cdh13", "Gdf10", "Bdnf", "Gfap"), mart = mouse_maRt)

save(target_GO, file = "06 inferring protein activity_a/GO for target genes.RData")
################################################################################
