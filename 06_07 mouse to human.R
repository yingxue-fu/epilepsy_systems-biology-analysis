library(biomaRt)
load("03 Gene coexpression network/gene_connect.RData")

human_maRt = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
save(human_maRt, file = "06 inferring protein activity_b/human_maRt.RData")

gene_modules_m2h = getLDS(attributes = c("mgi_symbol", "entrezgene"), 
                            filters = "mgi_symbol", 
                            values = gene_modules_connect$Gene_symbol, 
                            mart = mouse_maRt, 
                            attributesL = c("hgnc_symbol", "entrezgene"), 
                            martL = human_maRt) 

gene_modules_m2h = gene_modules_m2h[!gene_modules_m2h$HGNC.symbol=="",]

gene_modules_m2h <- merge(gene_modules_connect, gene_modules_m2h, by.x="Gene_symbol", by.y="MGI.symbol")

gene_modules_m2h_cq <- gene_modules_m2h %>% group_by(HGNC.symbol) %>% top_n(1,scaled_kWithin)

save(gene_modules_m2h_cq, file = "06_07 mouse to human/gene_modules_m2h_cq.RData")

# regulators mouse to human
load("06 inferring protein activity_b/mrs_all_anno.RData")
load("06 inferring protein activity_b/key_regulators.RData")
regulators_m2h = merge(mrs_alld_anno, gene_modules_m2h_cq[,c(1,7,8,9)], by.x="Regulon", by.y="Gene_symbol")
DEregulators_m2h = merge(DE_regulators, gene_modules_m2h_cq[,c(1,7,8,9)], by.x="Regulon", by.y="Gene_symbol")

save(DEregulators_m2h, file = "06_07 mouse to human/DEregulators_m2h.RData")
write.csv(DEregulators_m2h, file = "06_07 mouse to human/DEregulators_m2h.csv")

# pathway regulators
pathwayR <- read.table(file = "06_07 mouse to human/pathway_regulators.txt")

pathwayR_avtivity <- DEregulators_m2h[match(pathwayR$V1, DEregulators_m2h$HGNC.symbol),]
write.csv(pathwayR_avtivity, file = "06_07 mouse to human/pathwayR_avtivity.csv")

pathwayR_avtivity_reorder <- read.csv("06_07 mouse to human/pathwayR_avtivity_reorderd.csv")
library(pheatmap)
library(RColorBrewer)
# creat matrix
pathwayR_matrix = as.matrix(pathwayR_avtivity_reorder[,2:4])
rownames(pathwayR_matrix) <- pathwayR_avtivity_reorder$HGNC.symbol

# generate annotations
annotation_row = data.frame(Modules = factor(pathwayR_avtivity_reorder$Modules),
                            Type = factor(pathwayR_avtivity_reorder$type))
rownames(annotation_row) = pathwayR_avtivity_reorder$HGNC.symbol
# specify colors
module_color = pathwayR_avtivity_reorder[, 7:8][!duplicated(pathwayR_avtivity_reorder[, 7:8]), ]
m_color = as.character(module_color$Colors); names(m_color) = module_color$Modules
ann_colors = list(Type = c(SP = "#66A61E", Signal = "#E6AB02", TF = "#A6761D"),
                  Modules = m_color)


pheatmap(t(pathwayR_matrix), 
         color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100),
         annotation_col = annotation_row,
         annotation_colors = ann_colors, 
         treeheight_col = 10,
         treeheight_row = 15,
         show_rownames = T,
         show_colnames = T,
         fontsize = 9,
         fontsize_col = 7.5,
         cluster_cols = F,
         angle_col = 45,
         gaps_col = c(17,45)) 
