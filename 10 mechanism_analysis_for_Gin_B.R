library(tidyverse)


# M10 regulators
load("06 master regulator identification/mrs_all_anno.RData")

M10_regulators = mrs_alld_anno %>% filter(Modules == "M10")
write.csv(M10_regulators, "10 mechanism analysis for Gin_B/M10_regulators.csv")

# M10 genes
load("04 differential expression analysis/DEresults.RData")

M10_genes = DE_alld_anno %>% filter(Modules == "M10")

# subgraph of GB targets and M10 regulators and AED targets and epilepsy genes
GB_mech = c(GBL_compound_target_list$Ginkgolide_B, kregulator_interome_list$M10)

a = human_interactome$Gene_A_Entrez_ID %in% GB_mech
b = human_interactome$GeneB_Entrez_ID %in% GB_mech

M10network0 = human_interactome[a|b, ]

x = M10network0$Gene_A_Entrez_ID %in% c(probe_modules_anno_m2h$NCBI.gene.ID.1, 5724)
y = M10network0$GeneB_Entrez_ID %in% c(probe_modules_anno_m2h$NCBI.gene.ID.1, 5724)

M10network = M10network0[x&y, ]

M10_igraph <- graph_from_data_frame(M10network, directed = FALSE)
summary(M10_igraph)

M10_degree = degree(M10_igraph, v = V(M10_igraph), loops = F, normalized = FALSE)
M10_gene_anno = probe_modules_anno_m2h[match(names(M10_degree), as.character(probe_modules_anno_m2h$NCBI.gene.ID.1)), ]
M10_gene_anno$degree = M10_degree

degree1 = read.csv("10 mechanism analysis for Gin_B/degree1.csv", header = F)

s = M10network$Gene_A_Entrez_ID %in% degree1$V1
t = M10network$GeneB_Entrez_ID %in% degree1$V1

M10network_nod1 = M10network[!(s|t), ]
M10network_nod1$geneAsymbol = probe_modules_anno_m2h$HGNC.symbol[match(M10network_nod1$Gene_A_Entrez_ID, probe_modules_anno_m2h$NCBI.gene.ID.1)]
M10network_nod1$geneBsymbol = probe_modules_anno_m2h$HGNC.symbol[match(M10network_nod1$GeneB_Entrez_ID, probe_modules_anno_m2h$NCBI.gene.ID.1)]

write.csv(M10network_nod1, file = "10 mechanism analysis for Gin_B/M10network_nod1.csv")


# DE regulators and their targets
DEregulator_targets = subset.data.frame(acance_network, acance_network$Regulator %in% DEregulator_probe_anno$probe_ID)

save(DEregulator_targets, file = "06 master regulator identification/DEregulator_targets.RData")

# ARR3 ILMN_2717844
arr3_targets = subset.data.frame(DEregulator_targets, Regulator == "ILMN_2717844")
arr3_targets_info = merge(arr3_targets, DEG1, by.x = "Target", by.y = "probe_ID")
arr3_targets_info = merge(arr3_targets_info, probe_modules_anno, by.x = "Target", by.y = "probe_ID")

# GLRA2 ILMN_2729364
glra2_targets = subset.data.frame(DEregulator_targets, Regulator == "ILMN_2729364")
glra2_targets_info = merge(glra2_targets, DEG_info, by.x = "Target", by.y = "probe_ID")


# M13 pathway gene
M13_pathway_gene = read.csv("10 mechanism analysis for Gin_B/M13_pathway_gene.csv")
M13_pathway_gene_anno = merge(M13_pathway_gene, DEregulator_probe_anno, by.x = "Gene", by.y = "Gene_symbol")


# pcr gene primers
pcr_genes = read.table("04 differential expression analysis/PCR_gene_primer.txt", col.names = "Gene_symbol")
load("03 brain cell type enrichment/human_mouse_biomaRt.RData")
pcr_genes_h2m = getLDS(attributes = c("entrezgene", "hgnc_symbol"), 
                       filters = "hgnc_symbol", 
                       values = pcr_genes$Gene_symbol, 
                       mart = human, 
                       attributesL = c("entrezgene", "mgi_symbol"), 
                       martL = mouse) %>% na.omit()
pcr_genes_anno = merge(pcr_genes_h2m, probe_modules_anno, by.x = "NCBI.gene.ID.1", by.y = "Gene_ID")
write.csv(pcr_genes_anno, file = "04 differential expression analysis/pcr_genes_anno.csv")
