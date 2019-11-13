library(matrixStats)
library(tidyverse)

load("06_07 mouse to human/DEregulators_m2h.RData")
load("08 human interactome/protein_degree.RData")
load("08 human interactome/proteins_distance.RData")
load("08 human interactome/kregulator_interome.RData")
load("08 human interactome/AED_targets.RData")
load("08 human interactome/epilp_genes.RData")
load("08 human interactome/GBL_compound_target.RData")

# GBL compounds vs. kregulators_all
GBLs_kregulators_all = matrix(data = NA, nrow = length(GBL_compound_target_list), ncol = 1)
for (i in 1:length(GBL_compound_target_list)){
    GBLs_kregulators_all[i, 1] = Zscor(GBL_compound_target_list[[i]], kregulator_interome$NCBI.gene.ID.1, protein_distance)
    print(i)
}
rownames(GBLs_kregulators_all) = names(GBL_compound_target_list)

# GBL compounds vs. modules
kregulator_interome_list$M11 = NULL
kregulator_interome_list$M4 = NULL
GBLs_modules = matrix(data = NA, nrow = length(GBL_compound_target_list), ncol = length(kregulator_interome_list))
for (i in 1:length(GBL_compound_target_list)){
  for (j in 1:length(kregulator_interome_list)){
    GBLs_modules[i, j] = Zscor(GBL_compound_target_list[[i]], kregulator_interome_list[[j]], protein_distance)
    print(j)
  }
}
rownames(GBLs_modules) = names(GBL_compound_target_list)
colnames(GBLs_modules) = names(kregulator_interome_list)

save(GBLs_modules, file = "09 network proximity/GBLs_modules.RData")

write.csv(GBLs_modules, file = "09 network proximity calculation/GBLs_modules.csv")

GBLs_modules_sub = GBLs_modules[, c("M1", "M3", "M6", "M7", "M8")]
# GBLs_modules_sub[GBLs_modules_sub>0] <- NA

GBLs_modules_df = as.data.frame(GBLs_modules_sub)
GBLs_modules_df$Compound = rownames(GBLs_modules_sub)

GBLs_modules_tbl = GBLs_modules_df %>% gather(`M1`, `M3`, `M6`, `M7`, `M8`, key = "Modules", value = "Z_score") 

# GBC class
GBLs_class = read.csv("09 network proximity/GBC_class.csv")
GBLs_modules_tbl = merge(GBLs_modules_tbl, GBLs_class, by = "Compound")

# GBL compounds vs. epilepsy genes
GBLs_epilepsy_genes = matrix(data = NA, nrow = length(GBL_compound_target_list), ncol = 1)
for (i in 1:length(GBL_compound_target_list)){
  GBLs_epilepsy_genes[i, 1] = Zscor(GBL_compound_target_list[[i]], epilp_genes_interome$entrezgene, protein_distance)
  print(i)
}
rownames(GBLs_epilepsy_genes) = names(GBL_compound_target_list)


GBLs_epilepsy_genes <- as.data.frame(GBLs_epilepsy_genes)
GBLs_epilepsy_genes$Modules <- "EpilGene"
GBLs_epilepsy_genes$Compound <- rownames(GBLs_epilepsy_genes)
GBLs_epilepsy_genes = merge(GBLs_epilepsy_genes, GBLs_class, by = "Compound")
GBLs_epilepsy_genes <- GBLs_epilepsy_genes[,c(1,3,2,4)]
colnames(GBLs_epilepsy_genes)[3] <- "Z_score"

GBLs_modules_tbl <- rbind(GBLs_modules_tbl, GBLs_epilepsy_genes)

GBLs_modules_tbl$Class <- factor(GBLs_modules_tbl$Class, levels = c("TTLs", "Flavonoids","flavonoid oligomers","Carboxylic acids"))

save(GBLs_modules_tbl, file = "10 mechanism analysis for Gin_B/GBLs_modules_tbl.RData")

ggplot(GBLs_modules_tbl, aes(x=Modules, y=Z_score, color=Class, fill=Class)) +
  geom_boxplot(alpha=0.3) + 
  geom_dotplot(binaxis='y', stackdir='center', binwidth = 0.2, dotsize = 0.3, position=position_dodge(0.75), color="dimgray", alpha = 0.9) + 
  theme_classic() + 
  scale_x_discrete(limits=c("M1","M3", "M6", "M7", "M8", "EpilGene")) +
 # coord_flip() +
  geom_hline(yintercept = -0.15, linetype="dashed", size = 0.5, color = "black",alpha=0.6)

# GBL compounds vs. AED targets
GBLs_AED_targets = matrix(data = NA, nrow = length(GBL_compound_target_list), ncol = 1)
for (i in 1:length(GBL_compound_target_list)){
  GBLs_AED_targets[i, 1] = Zscor(GBL_compound_target_list[[i]], AED_targets_interome$ENTREZID, protein_distance)
  print(i)
}
rownames(GBLs_AED_targets) = names(GBL_compound_target_list)

GBLs_AED_targets_df <- as.data.frame(GBLs_AED_targets)
GBLs_AED_targets_df$Compound <- rownames(GBLs_AED_targets)
GBLs_AED_targets_df = merge(GBLs_AED_targets_df, GBLs_class, by = "Compound")
GBLs_AED_targets_df$Class <- factor(GBLs_AED_targets_df$Class, levels = c("Carboxylic acids","flavonoid oligomers", "Flavonoids","TTLs"))

ggplot(GBLs_AED_targets_df, aes(x=Class, y=V1, color=Class, fill=Class)) +
  geom_boxplot(alpha=0.3) + 
  geom_dotplot(binaxis='y', stackdir='center', binwidth = 0.4, dotsize = 1, position=position_dodge(0.75), color="dimgray", alpha = 0.9) + 
  theme_classic() + 
  coord_flip() +
  geom_hline(yintercept = -0.15, linetype="dashed", size = 0.5, color = "black",alpha=0.6)



save(GBLs_AED_targets, file = "09 network proximity/GBLs_AED_targets.RData")


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

save(AEDs_modules, file = "09 network proximity/AEDs_modules.RData")

AEDs_modules_df = as.data.frame(AEDs_modules)
AEDs_modules_df$Compound = rownames(AEDs_modules)

AEDs_modules_tbl = AEDs_modules_df %>% gather(`M1`, `M3`, `M6`, `M7`, `M8`, key = "Modules", value = "Z_score") %>% na.omit()

ggplot(AEDs_modules_tbl, aes(x=Modules, y=Z_score)) +
  geom_boxplot(alpha=0.2, color="#64b2cd", fill="#64b2cd") + 
  geom_dotplot(binaxis='y', stackdir='center', binwidth = 0.16, position=position_dodge(0.9), dotsize = 0.8,color="gray", fill="#64b2cd") + 
  theme_classic() + 
  scale_x_discrete(limits=c("M1","M3", "M6", "M7", "M8")) +
  geom_hline(yintercept = -0.15, linetype="dashed", size = 0.5, color = "black",alpha=0.6)

# AEDs vs. epilepsy genes
AEDs_epilepsy_genes = matrix(data = NA, nrow = length(AED_target_net_list), ncol = 1)
for (i in 1:length(AED_target_net_list)){
  AEDs_epilepsy_genes[i, 1] = Zscor(AED_target_net_list[[i]], epilp_genes_interome$entrezgene, protein_distance)
  print(i)
}
rownames(AEDs_epilepsy_genes) = names(AED_target_net_list)



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


# boxplot


ggplot(GBL_AED_epilepsy, aes(x=Class, y=Z.score, color=Class, fill=Class)) +
  geom_boxplot(alpha=0.2) + 
  geom_dotplot(binaxis='y', stackdir='center', binwidth = 0.2, dotsize = 0.5, position=position_dodge(1), color="black", fill = "gray", alpha = 0.7) + 
  scale_x_discrete(limits=c("AEDs", "Carboxylic acids", "flavonoid oligomers", "Flavonoids", "TTLs")) +
  theme_classic() + 
  theme(legend.position="none") +
  coord_flip()


#  theme(axis.text.x = element_text(face="bold", size=10, angle=30, hjust = 0.6, vjust = 0.8))

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




