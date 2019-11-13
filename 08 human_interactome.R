library(dplyr)
library(biomaRt)
library(igraph)

######################## read PPI data ################################
# protein interactions
human_interactome = read.table(file = "08 human interactome/SD1.txt", header = T)
# whether there are duplicated interactions
sum(duplicated(human_interactome[, c(1,2)]))
# unique proteins
uni_proteins = unique(c(human_interactome$Gene_A_Entrez_ID, human_interactome$GeneB_Entrez_ID)) %>% as.character()
# mapping to gene symbol using biomaRt
load("03 brain cell type enrichment/human_mouse_biomaRt.RData")
uni_gene_symbol = getBM(c("entrezgene", "hgnc_symbol"), filters = "entrezgene", values = uni_proteins, mart = human)

save(human_interactome, uni_proteins, uni_gene_symbol, file = "08 human interactome/human_interactome.RData")


############################# build PPI network ##########################
# creat PPI network using igraph
ppi_network <- graph_from_data_frame(human_interactome, directed = FALSE)
summary(ppi_network)

# Nodes and edges
ppi_v = V(ppi_network) # 15,970
ppi_e = E(ppi_network) # 217,160

save(ppi_network, ppi_v, ppi_e, file = "08 human interactome/interactome_igraph.RData")

# calculate degree (not including self-interaction)
protein_degree = degree(ppi_network, v = ppi_v, loops = F, normalized = FALSE)
degree_breaks = c(0:10, seq(15, 100, by = 5), seq(150, 450, by = 50), 500, 750, 10000)
protein_degree_bins = split(protein_degree, findInterval(protein_degree, degree_breaks))

save(protein_degree, degree_breaks, protein_degree_bins, file = "07 human interactome/protein_degree.RData")

# calculate distances
protein_distance = distances(ppi_network, weights = NA, algorithm = "unweighted")
protein_distance[1:5, 1:5]
hist(protein_distance)

save(protein_distance, file = "07 human interactome/proteins_distance.RData")

########################### mapping desired proteins ###################

# 1. key regulators in interactome
load("06_07 mouse to human/DEregulators_m2h.RData")

kregulator_m2h <- DEregulators_m2h

kregulator_m2h$NCBI.gene.ID.1 = as.character(kregulator_m2h$NCBI.gene.ID.1)
## key regulators in interactome
kregulator_interome = kregulator_m2h[kregulator_m2h$NCBI.gene.ID.1 %in% ppi_v$name, ]

kregulator_interome$Modules = as.character(kregulator_interome$Modules)
kregulator_interome_list = split(kregulator_interome$NCBI.gene.ID.1, kregulator_interome$Modules)

save(kregulator_interome, kregulator_interome_list, file = "08 human interactome/kregulator_interome.RData")

# 2. established epilepsy genes 
epilp_genes0 = read.csv("08 human interactome/epilepsy genes_84.csv")
## get thier entrezID
load("06 inferring protein activity_b/human_maRt.RData")
epilp_genes_entrezID = getBM(attributes = c("hgnc_symbol", "entrezgene_id", "external_gene_name"), 
                             filters = "hgnc_symbol", 
                             values = epilp_genes0$Gene,
                             mart = human_maRt)
epilp_genes = epilp_genes_entrezID[!duplicated(epilp_genes_entrezID$hgnc_symbol), ]
epilp_genes$entrezgene = as.character(epilp_genes$entrezgene)
## epilepsy genes in interactome
epilp_genes_interome = epilp_genes[epilp_genes$entrezgene %in% ppi_v$name, ]

save(epilp_genes, epilp_genes_interome, file = "07 human interactome/epilp_genes.RData")

write.csv(epilp_genes, file = "07 human interactome/epilp_genes.csv")

# 3. targets of AEDs
AED_targets0 = read.csv("08 human interactome/AED-targets.csv")
## convert UniprotID to entrezID
library(org.Hs.eg.db)
AED_targets_entrezID = select(org.Hs.eg.db, as.character(AED_targets0$Target_ID), "ENTREZID", "UNIPROT")
AED_targets_entrezID = AED_targets_entrezID[-84, ] # delete a wrong match
AED_targets = merge(AED_targets0, AED_targets_entrezID, by.x = "Target_ID", by.y = "UNIPROT")
## AED targets in interactome
AED_targets_interome = AED_targets[AED_targets$ENTREZID %in% ppi_v$name, ]
## merge with drug-target interactions
AED_target_net = read.csv("08 human interactome/AED-target-net.csv")
AED_target_net = merge(AED_target_net, AED_targets_interome, by = "Target_ID")
AED_target_net_list = split(AED_target_net$ENTREZID, AED_target_net$DRUG.NAME)

save(AED_targets, AED_targets_interome, AED_target_net, AED_target_net_list, file = "08 human interactome/AED_targets.RData")

# 4. GB compound-target interactions
GBL_compound_target0 = read.csv("08 human interactome/GB_C_T_0.65.csv")
# get target entrez ID
load("06 inferring protein activity_b/human_maRt.RData")
GBLtarget_entrezID = getBM(attributes = c("hgnc_symbol", "entrezgene_id"), filters = "hgnc_symbol", values = GBL_compound_target0$Target, mart = human_maRt)
# check duplicatin
sum(duplicated(GBLtarget_entrezID$hgnc_symbol))
# GBL targets in interactome
GBLtarget_interome = GBLtarget_entrezID[GBLtarget_entrezID$entrezgene %in% ppi_v$name, ]
# compound target lists
GBL_compound_target = merge(GBL_compound_target0, GBLtarget_interome, by.x = "Target", by.y = "hgnc_symbol")
GBL_compound_target$entrezgene = as.character(GBL_compound_target$entrezgene)
GBL_compound_target_list = split(GBL_compound_target$entrezgene, GBL_compound_target$Compound)

save(GBLtarget_entrezID, GBLtarget_interome, GBL_compound_target, GBL_compound_target_list, file = "08 human interactome/GBL_compound_target.RData")


# closest distance
closest_dis <- function(geneset1, geneset2, dismat){
  require(matrixStats)
  disteff <- dismat[geneset1, geneset2]
  if(length(geneset1) == 1){
    result <- min(disteff, na.rm = T)
  } else {
    closd = rowMins(disteff)
    closd[closd==Inf] = NA
    result <- mean(closd, na.rm = T) 
  }
  return(result)
}

closest_dis(epilp_genes_interome$entrezgene, kregulator_interome$NCBI.gene.ID.1, protein_distance) # 1.519481
closest_dis(AED_targets_interome$ENTREZID, kregulator_interome$NCBI.gene.ID.1, protein_distance) # 1.546392
real_closd = data.frame(type = c("AEDTar", "epilGen"), closd = c(1.546392, 1.519481))


# regulator randomness
ran_closd_AEDTar_kregul = vector(mode = "double", length = 1000)
ran_closd_epilGen_kregul = vector(mode = "double", length = 1000)
for (i in 1:1000){
  random_regulators = sample(ppi_v$name, 181)
  ran_closd_AEDTar_kregul[i] = closest_dis(AED_targets_interome$ENTREZID, random_regulators, protein_distance)
  ran_closd_epilGen_kregul[i] = closest_dis(epilp_genes_interome$entrezgene, random_regulators, protein_distance)
}

ran_closd0 = data.frame(AEDTar = ran_closd_AEDTar_kregul, epilGen = ran_closd_epilGen_kregul)
ran_closd = ran_closd0 %>% gather(`AEDTar`, `epilGen`, key = "type", value = "closd")

save(ran_closd0, ran_closd,file = "07 human interactome/ran_closd.RData")

ggplot(ran_closd, aes(x=closd, color=type, fill = type)) +
  geom_density(alpha=0.2)+
  scale_color_manual(values=c("#00BFC4", "#F8766D", "#999999"))+
  scale_fill_manual(values=c("#00BFC4", "#F8766D", "#999999"))+
  theme_classic() + 
  geom_vline(data=real_closd, aes(xintercept=closd, color=type))

sum(ran_closd0$AEDTar < 1.546392)/1000  # 0.003
sum(ran_closd0$epilGen < 1.519481)/1000 # 0.008
