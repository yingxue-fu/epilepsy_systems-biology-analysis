library(dplyr)

load("03 Gene coexpression network/gene_connect.RData")
# distribution of gene degree
hist(gene_connect$scaled_kWithin, breaks = 100, col = "gray80", border = "gray", freq = F)
quantile(gene_connect$scaled_kWithin, probs = c(0, 0.25,0.5,0.75,1))
abline(v = 0.2, col = "red")
gene_modules_connect_0.2 <- gene_modules_connect %>% filter(scaled_kWithin>0.2)

# distribution of DEGs in modules
load("02 Gene meta-signature/metaRP.DEG.RData")
all.DEG = data.frame()
for (i in 1:length(metaRP.DEG_overlap)){
  x = data.frame(Gene_symbol = metaRP.DEG_overlap[[i]]); x$group = names(metaRP.DEG_overlap[i])
  all.DEG = rbind(all.DEG, x)
}

all.DEG.moduels <- merge(all.DEG, gene_modules_connect, by = "Gene_symbol")

save(all.DEG.moduels, file = "05 Calculate MAS/all.DEG.moduels.RData")
write.csv(all.DEG.moduels, file = "05 Calculate MAS/all.DEG.moduels.csv")

table(all.DEG$group, all.DEG$Modules)

hist(all.DEG$scaled_kWithin, breaks = 100, col = "gray80", border = "gray")

all.DEG %>% group_by(group) %>% summarise(mean_k = mean(scaled_kWithin))




# fgsea
library(fgsea)

# module gene sets
module_genes = gene_modules_connect_0.2[, c(2, 1)]
module_genes_list = split(as.character(module_genes$Gene_symbol), module_genes$Modules)

# cell gene sets
cell_genes <- read.csv("04 brain cell type enrichment/science_singleCell_brainCell_types.csv")
cell_genes_list = split(as.character(cell_genes$Gene), cell_genes$Category)

# rank gene lists
acute_gene_rank <- metaRP_acute_list$logFC; names(acute_gene_rank) <- metaRP_acute_list$ID
acute_gene_rank = sort(acute_gene_rank, decreasing = T)

latent_gene_rank <- metaRP_latent_list$logFC; names(latent_gene_rank) <- metaRP_latent_list$ID
latent_gene_rank = sort(latent_gene_rank, decreasing = T)

chronic_gene_rank <- metaRP_chronic_list$logFC; names(chronic_gene_rank) <- metaRP_chronic_list$ID
chronic_gene_rank = sort(chronic_gene_rank, decreasing = T)

# running module fgsea
fgseaRes_acute <- fgsea(module_genes_list, acute_gene_rank, minSize=15, maxSize=3000, nperm=1000)
fgseaRes_latent <- fgsea(module_genes_list, latent_gene_rank, minSize=15, maxSize=3000, nperm=1000)
fgseaRes_chronic <- fgsea(module_genes_list, chronic_gene_rank, minSize=15, maxSize=3000, nperm=1000)

plotGseaTable(module_genes_list[order(fgseaRes_acute$NES, decreasing = T)+1], acute_gene_rank, fgseaRes_acute, gseaParam = 0.5)
plotGseaTable(module_genes_list[order(fgseaRes_latent$NES, decreasing = T)+1], latent_gene_rank, fgseaRes_latent, gseaParam = 0.5)
plotGseaTable(module_genes_list[order(fgseaRes_chronic$NES, decreasing = T)+1], chronic_gene_rank, fgseaRes_chronic, gseaParam = 0.5)

# running cell fgsea
fgseaRes_acute <- fgsea(cell_genes_list, acute_gene_rank, minSize=15, maxSize=3000, nperm=1000)
fgseaRes_latent <- fgsea(cell_genes_list, latent_gene_rank, minSize=15, maxSize=3000, nperm=1000)
fgseaRes_chronic <- fgsea(cell_genes_list, chronic_gene_rank, minSize=15, maxSize=3000, nperm=1000)

plotGseaTable(cell_genes_list[order(fgseaRes_acute$NES, decreasing = T)], acute_gene_rank, fgseaRes_acute, gseaParam = 0.5)
plotGseaTable(cell_genes_list[order(fgseaRes_latent$NES, decreasing = T)], latent_gene_rank, fgseaRes_latent, gseaParam = 0.5)
plotGseaTable(cell_genes_list[order(fgseaRes_chronic$NES, decreasing = T)], chronic_gene_rank, fgseaRes_chronic, gseaParam = 0.5)

save(fgseaRes_acute, fgseaRes_latent, fgseaRes_chronic, file = "05 Calculate MAS/fgseaRes.RData")




# draw the plot
library(ggrepel)
library(ggpubr)

mas_acute <- 
  ggplot(fgseaRes_acute, aes(x = NES, y = -log10(padj), 
                                    fill = as.character(sign(NES)), size = size)) + 
  geom_point(shape = 21, alpha = 0.5, colour="black") + 
  theme_bw() + 
  geom_text_repel(aes(label = pathway), size = 3, colour= "gray20") + 
  theme(legend.position='none') +
  geom_hline(yintercept = -log10(0.05), linetype="dashed", size = 0.5, alpha=0.5) +
  geom_vline(xintercept = 0, linetype="dashed", size = 0.5, alpha=0.5) +
  scale_x_continuous(name="MAS", limits=c(-4, 4), breaks = c(-4, -2, 0, 2, 4))+
  scale_y_continuous(limits=c(0.25, 2.5)) +
  scale_fill_brewer(palette="Dark2")

mas_latent <- 
  ggplot(fgseaRes_latent, aes(x = NES, y = -log10(padj), 
                                      fill = as.character(sign(NES)),size=size)) + 
  geom_point(shape = 21, alpha = 0.5, colour="black") + 
  theme_bw() + 
  geom_text_repel(aes(label = pathway), size = 3, colour= "gray20") + 
  theme(legend.position='none') +
  geom_hline(yintercept = -log10(0.05), linetype="dashed", size = 0.5, alpha=0.5) +
  geom_vline(xintercept = 0, linetype="dashed", size = 0.5, alpha=0.5) +
  scale_x_continuous(name="MAS", limits=c(-4, 4), breaks = c(-4, -2, 0, 2, 4))+
  scale_y_continuous(limits=c(0.25, 2.5)) +
  scale_fill_brewer(palette="Dark2")

mas_chronic <- 
  ggplot(fgseaRes_chronic, aes(x = NES, y = -log10(padj), 
                                      fill = as.character(sign(NES)), size = size)) + 
  geom_point(shape = 21, alpha = 0.5, colour="black") + 
  theme_bw() + 
  geom_text_repel(aes(label = pathway), size = 3, colour= "gray20") + 
  theme(legend.position='none') +
  geom_hline(yintercept = -log10(0.05), linetype="dashed", size = 0.5, alpha=0.5) +
  geom_vline(xintercept = 0, linetype="dashed", size = 0.5, alpha=0.5) +
  scale_x_continuous(name="MAS", limits=c(-4, 4), breaks = c(-4, -2, 0, 2, 4))+
  scale_y_continuous(limits=c(0.25, 2.5)) +
  scale_fill_brewer(palette="Dark2")

save(mas_acute, mas_latent, mas_chronic, file = "05 Calculate MAS/MAS_stages.RData")

ggarrange(mas_acute, mas_latent, mas_chronic, ncol = 3, nrow = 1)
