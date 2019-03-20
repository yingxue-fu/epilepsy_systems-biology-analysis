library(matrixStats)
library(tidyverse)


load("06 master regulator identification/DEregulator_probe_anno.RData")
load("07 human interactome/protein_degree.RData")
load("07 human interactome/proteins_distance.RData")
load("07 human interactome/kregulator_interome.RData")
load("07 human interactome/AED_targets.RData")
load("08 drug target interactions/GBL_compound_target.RData")


# GBL compounds vs. epilepsy (all key regulators)
GBLs_epilepsy = matrix(data = NA, nrow = length(GBL_compound_target_list), ncol = 1)
for (i in 1:length(GBL_compound_target_list)){
  GBLs_epilepsy[i] = Zscor(GBL_compound_target_list[[i]], as.character(kregulator_interome$NCBI.gene.ID.1), protein_distance)
  print(i)
}
GBLs_epilepsy_df = data.frame(Compound = names(GBL_compound_target_list), 
                              Z_score = GBLs_epilepsy)
# GBC class
GBC_class = read.csv("09 network proximity calculation/GBC_class.csv")
GBLs_epilepsy_df = merge(GBLs_epilepsy_df, GBC_class, by = "Compound")

# AEDs vs. epilepsy (all key regulators)
AEDs_epilepsy = matrix(data = NA, nrow = length(AED_target_net_list), ncol = 1)
for (i in 1:length(AED_target_net_list)){
  AEDs_epilepsy[i] = Zscor(AED_target_net_list[[i]], as.character(kregulator_interome$NCBI.gene.ID.1), protein_distance)
  print(i)
}
AEDs_epilepsy_df = data.frame(Compound = names(AED_target_net_list), 
                              Z.score = AEDs_epilepsy) %>% na.omit()
AEDs_epilepsy_df$Class = "AEDs"
GBL_AED_epilepsy = rbind(GBLs_epilepsy_df, AEDs_epilepsy_df)

ggplot(GBL_AED_epilepsy, aes(x=Class, y=Z.score, fill=Class)) +
  geom_boxplot(alpha=0.2) + 
  geom_dotplot(binaxis='y', stackdir='center', binwidth = 0.2, dotsize = 0.5, position=position_dodge(1), color="gray", fill = "gray", alpha = 0.7) + 
  scale_x_discrete(limits=c("AEDs", "Carboxylic acids", "Flavonoids", "TTLs")) +
  theme_classic() + 
  theme(legend.position="none") +
  coord_flip()


#  theme(axis.text.x = element_text(face="bold", size=10, angle=30, hjust = 0.6, vjust = 0.8))
save(GBLs_epilepsy_df, AEDs_epilepsy_df, GBL_AED_epilepsy, file = "09 network proximity calculation/compound_epilepsy_Zscor.RData")

# GBL compounds vs. modules
kregulator_interome_list$M14 = NULL
GBLs_modules = matrix(data = NA, nrow = length(GBL_compound_target_list), ncol = length(kregulator_interome_list))
for (i in 1:length(GBL_compound_target_list)){
  for (j in 1:length(kregulator_interome_list)){
    GBLs_modules[i, j] = Zscor(GBL_compound_target_list[[i]], kregulator_interome_list[[j]], protein_distance)
    print(j)
  }
}
rownames(GBLs_modules) = names(GBL_compound_target_list)
colnames(GBLs_modules) = names(kregulator_interome_list)
write.csv(GBLs_modules, file = "09 network proximity calculation/GBLs_modules.csv")

GBLs_modules_sub = GBLs_modules[, c("M2", "M3", "M6", "M7", "M10")]
GBLs_modules_sub[GBLs_modules_sub>0] <- NA

GBLs_modules_df = as.data.frame(GBLs_modules_sub)
GBLs_modules_df$Compound = rownames(GBLs_modules_sub)

GBLs_modules_tbl = GBLs_modules_df %>% gather(`M2`, `M3`, `M6`, `M7`, `M10`, key = "Modules", value = "Z_score") %>% na.omit()

# AEDs vs. modules
AEDs_modules = matrix(data = NA, nrow = length(AED_target_net_list), ncol = length(kregulator_interome_list))
for (i in 1:length(AED_target_net_list)){
  for (j in 1:length(kregulator_interome_list)){
    AEDs_modules[i, j] = Zscor(AED_target_net_list[[i]], kregulator_interome_list[[j]], protein_distance)
    print(j)
  }
}
rownames(AEDs_modules) = names(AED_target_net_list)
colnames(AEDs_modules) = names(kregulator_interome_list)


AEDs_modules_sub = AEDs_modules[, c("M2", "M3", "M6", "M7", "M10")]
AEDs_modules_sub[AEDs_modules_sub>0] <- NA

AEDs_modules_df = as.data.frame(AEDs_modules_sub)
AEDs_modules_df$Compound = rownames(AEDs_modules_sub)

AEDs_modules_tbl = AEDs_modules_df %>% gather(`M2`, `M3`, `M6`, `M7`, `M10`, key = "Modules", value = "Z_score") %>% na.omit()

save(GBLs_modules, GBLs_modules_tbl, AEDs_modules, AEDs_modules_tbl, file = "09 network proximity calculation/compound_module_Zscore.RData")

# volinplot

GBLs_modules_tbl$Compd_type = "GBLs"
AEDs_modules_tbl$Compd_type = "AEDs"
GBL_AED_modules = rbind(GBLs_modules_tbl, AEDs_modules_tbl) 

ggplot(GBL_AED_modules, aes(x=Modules, y=Z_score, color=Compd_type, fill=Compd_type)) +
  geom_violin(alpha=0.7) + 
  geom_dotplot(binaxis='y', stackdir='center', binwidth = 0.16, position=position_dodge(0.9), dotsize = 0.8,color="black") + 
  theme_classic() + 
  scale_x_discrete(limits=c("M2","M3", "M6", "M7", "M10"))



# heatmap
GBLs_modules_sub0 = GBLs_modules[, c("M2", "M3", "M6", "M7", "M10")]
GBLs_modules_sub0[GBLs_modules_sub0>0] <- 0
pheatmap::pheatmap(GBLs_modules_sub0, color = colorRampPalette(c("firebrick1", "white"))(100), cluster_cols = F)





################################# Network proxicimity Function #############################
# calculate drug-disease Z-score
Zscor <- function(drug_target, disease_gene, dismat){
  drugDise_dist <- closest_dis(drug_target, disease_gene, dismat)
  alltemp <- numeric(1000)
  for (i in 1:1000) {
    alltemp[i] = randomdis(drug_target, disease_gene, dismat)
  }
  zscore <- (drugDise_dist - mean(alltemp)) / sd(alltemp)
  return(zscore)
}

# average closest distance between drug targets and disease genes
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

# average closest distance between two random sets of genes
randomdis <- function(drug_target, disease_gene, dismat){
  set1 <- pick_random_genes(drug_target, protein_degree)
  set2 <- pick_random_genes(disease_gene, protein_degree)
  return(closest_dis(set1, set2, dismat))
}

# pick random genes with matched size and degree
pick_random_genes <- function(selected_genes, all_genes_degree){
  degree_breaks = c(0:10, seq(15, 100, by = 5), seq(150, 450, by = 50), 500, 750, 10000)
  selected_degree_hist = as.matrix(table(findInterval(all_genes_degree[selected_genes], degree_breaks)))
  random_nodes_all = numeric()
  for (i in 1:nrow(selected_degree_hist)) {
    n = as.numeric(rownames(selected_degree_hist)[i])
    random_nodes = sample(protein_degree_bins[[n]], selected_degree_hist[i])
    random_nodes_all = c(random_nodes_all, random_nodes)
  }
  return(names(random_nodes_all))
}

# test
all_genes_degree[selected_genes]
all_genes_degree[pick_random_genes(selected_genes, protein_degree)]

######################################################################################




