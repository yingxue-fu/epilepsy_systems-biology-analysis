library(Biobase)
library(GEOquery)
library(limma)
library(RColorBrewer)
library(WGCNA)
library(tidyverse)
library(Category)
library(pheatmap)
library(biomaRt)
library(clusterProfiler)
library(enrichplot)

load("02 gene coexpression network analysis/probe_modules_anno.RData")

# load series matrix data from GEO
gset <- getGEO(filename = "01 sample and gene filtering/GSE73878_series_matrix.txt.gz", GSEMatrix = TRUE, getGPL = F)
load("01 sample and gene filtering/samples_and_probes_filtering.RData")
dim(datExpr)
gset <- gset[rownames(datExpr), colnames(datExpr)]

# set group names
setgroups <- read.csv("04 differential expression analysis/new_sample_group.csv")
setgroups = setgroups[match(colnames(datExpr), setgroups$GSM_ID), ]
fl <- as.factor(setgroups$new_Group)
gset$description <- fl

design <- model.matrix(~ description + 0, gset)
colnames(design) <- levels(fl)

save(gset, design, file = "04 differential expression analysis/gset_design.RData")

# G0-Sham, G1-KA_7d, G2-KA_28d, G3-KA_60d
fit <- lmFit(gset, design)
cont.matrix <- makeContrasts(G1-G0, G2-G0, G3-G0, levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2, 0.01)

DEresults = decideTests(fit2)
mfrow.old <- par()$mfrow
par(mfrow=c(1,2))
vennDiagram(vennCounts(DEresults), 
            names=c("7 day", "28 day", "60 day"), 
            cex=c(1, 0.8, 0.3), 
            circle.col = c("chocolate1", "green3", "purple"))
vennDiagram(DEresults, 
            names=c("7 day", "28 day", "60 day"), 
            include=c("up", "down"),
            counts.col=c("red", "blue"), 
            cex=c(1, 0.8, 0.3), 
            circle.col = c("chocolate1", "green3", "purple"))

save(fit2, file = "04 differential expression analysis/fit2.RData")

# GSEA analysis of gene logFC on module gene sets
DE_7d <- topTable(fit2, coef = 1, number = 12000)
DE_28d <- topTable(fit2, coef = 2, number = 12000)
DE_60d <- topTable(fit2, coef = 3, number = 12000)
DE_alld <- topTable(fit2, number = 12000)

# match to probe annotation
DE_alld$probe_ID = rownames(DE_alld)
DE_alld_anno = merge(DE_alld, probe_modules_anno, by = "probe_ID")

save(DE_7d, DE_28d, DE_60d, DE_alld, DE_alld_anno, file = "04 differential expression analysis/DEresults.RData")

# t-statistic ranked gene signature
DE_7d_t = DE_7d$t; names(DE_7d_t) = rownames(DE_7d)
DE_7d_t = sort(DE_7d_t, decreasing = T)
DE_28d_t = DE_28d$t; names(DE_28d_t) = rownames(DE_28d)
DE_28d_t = sort(DE_28d_t, decreasing = T)
DE_60d_t = DE_60d$t; names(DE_60d_t) = rownames(DE_60d)
DE_60d_t = sort(DE_60d_t, decreasing = T)

DE_7d_t_GSEA = GSEA(DE_7d_t, TERM2GENE = probe_modules_anno[, c(2, 1)], maxGSSize = 3000, pvalueCutoff = 1, verbose=FALSE)
DE_28d_t_GSEA = GSEA(DE_28d_t, TERM2GENE = probe_modules_anno[, c(2, 1)], maxGSSize = 3000, pvalueCutoff = 1, verbose=FALSE)
DE_60d_t_GSEA = GSEA(DE_60d_t, TERM2GENE = probe_modules_anno[, c(2, 1)], maxGSSize = 3000, pvalueCutoff = 1, verbose=FALSE)

DE_7d_t_GSEA_df = DE_7d_t_GSEA@result %>% na.omit()
DE_28d_t_GSEA_df = DE_28d_t_GSEA@result %>% na.omit()
DE_60d_t_GSEA_df = DE_60d_t_GSEA@result %>% na.omit()

save(DE_7d_t_GSEA, DE_28d_t_GSEA, DE_60d_t_GSEA, DE_7d_t_GSEA_df, DE_28d_t_GSEA_df, DE_60d_t_GSEA_df, file = "04 differential expression analysis/DE_t_GSEA.RData")

# draw the plot
library(ggrepel)
library(ggpubr)

p_7d <- ggplot(DE_7d_t_GSEA_df, aes(x = NES, y = -log10(p.adjust), 
                                  fill = as.character(-sign(NES)))) + 
  geom_point(size = 4, shape = 21, alpha = 0.7) + 
  theme_bw() + 
  geom_text_repel(aes(label = ID), size = 3.5) + 
  scale_fill_discrete() + 
  theme(legend.position='none') +
  geom_hline(yintercept = -log10(0.05), linetype="dashed", size = 0.5) +
  geom_vline(xintercept = 0, linetype="dashed", size = 0.5)

p_28d <- ggplot(DE_28d_t_GSEA_df, aes(x = NES, y = -log10(p.adjust), 
                                    fill = as.character(-sign(NES)))) + 
  geom_point(size = 4, shape = 21, alpha = 0.7) + 
  theme_bw() + 
  geom_text_repel(aes(label = ID), size = 3.5) + 
  scale_fill_discrete() + 
  theme(legend.position='none') +
  geom_hline(yintercept = -log10(0.05), linetype="dashed", size = 0.5) +
  geom_vline(xintercept = 0, linetype="dashed", size = 0.5)

p_60d <- ggplot(DE_60d_t_GSEA_df, aes(x = NES, y = -log10(p.adjust), 
                                    fill = as.character(-sign(NES)))) + 
  geom_point(size = 4, shape = 21, alpha = 0.7) + 
  theme_bw() + 
  geom_text_repel(aes(label = ID), size = 3.5) + 
  scale_fill_discrete() + 
  theme(legend.position='none') +
  geom_hline(yintercept = -log10(0.05), linetype="dashed", size = 0.5) +
  geom_vline(xintercept = 0, linetype="dashed", size = 0.5)

save(p_7d, p_28d, p_60d, file = "04 differential expression analysis/module_evolution.RData")

ggarrange(p_7d, p_28d, p_60d, ncol = 3, nrow = 1)


NES7 = DE_7d_GSEA_df[DE_7d_GSEA_df$Modules %in% c("M2", "M9", "M10", "M11", "M14"), c("NES", "Modules")]; NES7$Stage = "7 day"
NES28 = DE_28d_GSEA_df[DE_28d_GSEA_df$Modules %in% c("M2", "M9", "M10", "M11", "M14"), c("NES", "Modules")]; NES28$Stage = "28 day"
NES60 = DE_60d_GSEA_df[DE_60d_GSEA_df$Modules %in% c("M2", "M9", "M10", "M11", "M14"), c("NES", "Modules")]; NES60$Stage = "60 day"

module_dynamics = rbind(NES7, NES28, NES60)


ggplot(module_dynamics, aes(x = Stage, y = NES, group = Modules, color = Modules, fill = Modules)) +
  geom_line() +
  geom_point(shape = 21, size = 3, alpha = 0.7) +
  scale_x_discrete(limits=c("7 day", "28 day", "60 day")) +
  theme_bw() + 
  geom_hline(yintercept = 0, linetype="dashed")

