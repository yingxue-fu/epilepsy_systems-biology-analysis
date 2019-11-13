library(dplyr)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)

# ADME analysis for GBL compounds
GBL_adme <- read.csv(file = "10 mechanism analysis for Gin_B/GBL_ADME.csv")
GBL_adme$Compound_class <- factor(GBL_adme$Compound_class, levels = c("TTLs", "Flavonoids","Flavonoid oligomers","Carboxylic acids"))


gbl_mw <- ggplot(GBL_adme, aes(x=Compound_class, y=MW, color =Compound_class, fill=Compound_class)) +
  geom_violin(alpha=0.7) + 
  geom_dotplot(binaxis='y', stackdir='center', binwidth = 8, dotsize = 1.5, position=position_dodge(1), color="black") + 
  theme_bw() + 
  theme(legend.position="none") +
  scale_fill_manual(values = c("#7189BF", "#df7599", "#ffc785","#72d6c9")) +
  scale_color_manual(values = c("#7189BF", "#df7599", "#ffc785","#72d6c9")) +
  theme(axis.text.x = element_blank()) +
  geom_hline(yintercept=500, linetype="dashed", color = "red",alpha=0.5)

gbl_hbd <- ggplot(GBL_adme, aes(x=Compound_class, y=HBD, color =Compound_class, fill=Compound_class)) +
  geom_violin(alpha=0.7) + 
  geom_dotplot(binaxis='y', stackdir='center', binwidth = 0.18, dotsize = 1.5, position=position_dodge(1), color="black") + 
  theme_bw() + 
  theme(legend.position="none") +
  scale_fill_manual(values = c("#7189BF", "#df7599", "#ffc785","#72d6c9")) +
  scale_color_manual(values = c("#7189BF", "#df7599", "#ffc785","#72d6c9")) +
  theme(axis.text.x = element_blank()) +
  scale_y_continuous(limits = c(2,13)) +
  geom_hline(yintercept=5, linetype="dashed", color = "red",alpha=0.5)

gbl_hba <- ggplot(GBL_adme, aes(x=Compound_class, y=HBA, color =Compound_class, fill=Compound_class)) +
  geom_violin(alpha=0.7) + 
  geom_dotplot(binaxis='y', stackdir='center', binwidth = 0.2, dotsize = 1.5, position=position_dodge(1), color="black") + 
  theme_bw() + 
  theme(legend.position="none") +
  scale_fill_manual(values = c("#7189BF", "#df7599", "#ffc785","#72d6c9")) +
  scale_color_manual(values = c("#7189BF", "#df7599", "#ffc785","#72d6c9")) +
  theme(axis.text.x = element_blank())+
  geom_hline(yintercept=10, linetype="dashed", color = "red",alpha=0.5)

gbl_alogp <- ggplot(GBL_adme, aes(x=Compound_class, y=AlogP, color =Compound_class, fill=Compound_class)) +
  geom_violin(alpha=0.7) + 
  geom_dotplot(binaxis='y', stackdir='center', binwidth = 0.18, dotsize = 1.5, position=position_dodge(1), color="black") + 
  theme_bw() + 
  theme(legend.position="none") +
  scale_fill_manual(values = c("#7189BF", "#df7599", "#ffc785","#72d6c9")) +
  scale_color_manual(values = c("#7189BF", "#df7599", "#ffc785","#72d6c9")) +
  theme(axis.text.x = element_blank()) +
  geom_hline(yintercept=5, linetype="dashed", color = "red",alpha=0.5)

gbl_logs <- ggplot(GBL_adme, aes(x=Compound_class, y=logS, color =Compound_class, fill=Compound_class)) +
  geom_violin(alpha=0.7) + 
  geom_dotplot(binaxis='y', stackdir='center', binwidth = 0.1, dotsize = 1.5, position=position_dodge(1), color="black") + 
  theme_bw() + 
  theme(legend.position="none") +
  scale_fill_manual(values = c("#7189BF", "#df7599", "#ffc785","#72d6c9")) +
  scale_color_manual(values = c("#7189BF", "#df7599", "#ffc785","#72d6c9")) +
  theme(axis.text.x = element_blank()) +
  scale_y_continuous(limits = c(-4,1)) +
  geom_hline(yintercept=-4, linetype="dashed", color = "red",alpha=0.5)

gbl_ppb <- ggplot(GBL_adme, aes(x=Compound_class, y=PPB, color =Compound_class, fill=Compound_class)) +
  geom_violin(alpha=0.7) + 
  geom_dotplot(binaxis='y', stackdir='center', binwidth = 0.021, dotsize = 1.5, position=position_dodge(1), color="black") + 
  theme_bw() + 
  theme(legend.position="none") +
  scale_fill_manual(values = c("#7189BF", "#df7599", "#ffc785","#72d6c9")) +
  scale_color_manual(values = c("#7189BF", "#df7599", "#ffc785","#72d6c9")) +
  theme(axis.text.x = element_blank()) +
  geom_hline(yintercept=0.9, linetype="dashed", color = "red",alpha=0.5)

library(ggpubr)
ggarrange(gbl_mw, gbl_hbd, gbl_hba, gbl_alogp, gbl_logs, gbl_ppb, labels = c("A", "B", "C","D", "E", "F"), ncol = 3, nrow = 2, align = "v")

save.image(file = "10 mechanism analysis for Gin_B/GBL_adme.RData")

# legend
ggplot(GBL_adme, aes(x=Compound_class, y=PPB, color =Compound_class, fill=Compound_class)) +
  geom_violin(alpha=0.7) + 
  geom_dotplot(binaxis='y', stackdir='center', binwidth = 0.021, dotsize = 1.5, position=position_dodge(1), color="black") + 
  theme_bw() + 
  scale_fill_manual(values = c("#7189BF", "#df7599", "#ffc785","#72d6c9")) +
  scale_color_manual(values = c("#7189BF", "#df7599", "#ffc785","#72d6c9")) +
  theme(axis.text.x = element_blank()) +
  geom_hline(yintercept=0.9, linetype="dashed", color = "red",alpha=0.5)

# M1 regulators
load("08 human interactome/kregulator_interome.RData")
load("08 human interactome/GBL_compound_target.RData")
load("08 human interactome/human_interactome.RData")

M1_regulators = kregulator_interome %>% filter(Modules == "M1")
write.csv(M1_regulators, file = "10 mechanism analysis for Gin_B/M1_regulators.csv")

# subgraph of GB targets and M10 regulators
GB_mech = c(GBL_compound_target_list$Ginkgolide_B, kregulator_interome_list$M1)

a = human_interactome$Gene_A_Entrez_ID %in% GB_mech
b = human_interactome$GeneB_Entrez_ID %in% GB_mech

M1_PPI_network = human_interactome[a&b,]

write.csv(M1_PPI_network, file = "10 mechanism analysis for Gin_B/M1_PPI_network.csv")

# ARR3
a = human_interactome$Gene_A_Entrez_ID == "407"
b = human_interactome$GeneB_Entrez_ID == "407"

ARR3_network = human_interactome[a|b,]
ARR3_network_node <- c(ARR3_network$Gene_A_Entrez_ID,ARR3_network$GeneB_Entrez_ID) %>% unique() %>% as.character()

load("06 inferring protein activity_b/human_maRt.RData")
ARR3_network_GENE <- getBM(attributes = c("entrezgene_id", "hgnc_symbol"), 
                           filters = "entrezgene_id",
                           values = ARR3_network_node, 
                           mart = human_maRt)

GB_mech_arr3 = c(GB_mech, ARR3_network_node) %>% unique()

a = human_interactome$Gene_A_Entrez_ID %in% GB_mech_arr3
b = human_interactome$GeneB_Entrez_ID %in% GB_mech_arr3

M1_PPI_network_1 = human_interactome[a&b,]

write.csv(M1_PPI_network_1, file = "10 mechanism analysis for Gin_B/M1_PPI_network_1.csv")


M1_PPI_network_node <- read.csv(file = "10 mechanism analysis for Gin_B/M1_PPI_network_node.csv")

M1R_avtivity <- DEregulators_m2h[match(M1_PPI_network_node$name, DEregulators_m2h$HGNC.symbol),]
write.csv(M1R_avtivity, file = "10 mechanism analysis for Gin_B/M1R_avtivity.csv")

M1R_avtivity_reorder <- read.csv("10 mechanism analysis for Gin_B/M1R_avtivity_reorder.csv")

library(pheatmap)
library(RColorBrewer)
# creat matrix
M1R_matrix = as.matrix(M1R_avtivity_reorder[,2:4])
rownames(M1R_matrix) <- M1R_avtivity_reorder$HGNC.symbol

# generate annotations
annotation_row = data.frame(Type = factor(M1R_avtivity_reorder$type))
rownames(annotation_row) = M1R_avtivity_reorder$HGNC.symbol
# specify colors
ann_colors = list(Type = c(SP = "#C7F0DB", Signal = "#8BBABB", TF = "#6C7B95"))


pheatmap(M1R_matrix, 
         color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100),
         annotation_row = annotation_row,
         annotation_colors = ann_colors, 
         treeheight_col = 10,
         treeheight_row = 15,
         show_rownames = T,
         show_colnames = T,
         fontsize = 8,
         fontsize_col = 6,
         cluster_rows = F,
         cluster_cols = F,
         angle_col = 45) 









other_node1 = human_interactome[a, 2]
other_node2 = human_interactome[b, 1]

M10network_node = unique(c(as.character(intersect(other_node1, other_node2)), GB_mech))

x = human_interactome$Gene_A_Entrez_ID %in% M10network_node
y = human_interactome$GeneB_Entrez_ID %in% M10network_node

M10network = rbind(human_interactome[a&y, ], human_interactome[b&x, ]) %>% unique()

write.csv(M10network, file = "10 mechanism analysis for Gin_B/M10_PPInetwork.csv")

library(biomaRt)
load("05 conservation in human/human_mouse_biomaRt.RData")
M10network_node_anno = getBM(attributes = c("entrezgene", "hgnc_symbol"), 
                             filters = "entrezgene", 
                             values = M10network_node, 
                             mart = human)



