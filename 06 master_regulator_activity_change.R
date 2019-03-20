library(viper)
library(pheatmap)
library(limma)
library(dplyr)
library(ggplot2)
library(ggpubr)

library(RColorBrewer)
library(msigdf)
library(clusterProfiler)
library(enrichplot)

load("04 differential expression analysis/gset_design.RData")

# Generating the regulon object form ARACNe output file
regulon_target_net = aracne2regulon("06 master regulator identification/network_3col.txt", gset, format = "3col", verbose = TRUE)
print(regulon_target_net)

save(regulon_target_net, file = "06 master regulator identification/aracne2regulon.RData")

# Master Regulator Analysis performed by msVIPER
## Generating the gene expression signatures
signature_7d <- rowTtest(gset, "description", "G1", "G0")
signature_28d <- rowTtest(gset, "description", "G2", "G0")
signature_60d <- rowTtest(gset, "description", "G3", "G0")

# Z-transformation
signature_7d_z <- (qnorm(signature_7d$p.value/2, lower.tail = FALSE) * sign(signature_7d$statistic))[, 1]
signature_28d_z <- (qnorm(signature_28d$p.value/2, lower.tail = FALSE) * sign(signature_28d$statistic))[, 1]
signature_60d_z <- (qnorm(signature_60d$p.value/2, lower.tail = FALSE) * sign(signature_60d$statistic))[, 1]

## NULL model by sample permutations
nullmodel_7d <- ttestNull(gset, "description", "G1", "G0", per = 1000, 
                       repos = TRUE, verbose = FALSE)
nullmodel_28d <- ttestNull(gset, "description", "G2", "G0", per = 1000, 
                          repos = TRUE, verbose = FALSE)
nullmodel_60d <- ttestNull(gset, "description", "G3", "G0", per = 1000, 
                          repos = TRUE, verbose = FALSE)

## msVIPER
mrs_7d <- msviper(signature_7d_z, regulon_target_net, nullmodel_7d, verbose = FALSE)
mrs_28d <- msviper(signature_28d_z, regulon_target_net, nullmodel_28d, verbose = FALSE)
mrs_60d <- msviper(signature_60d_z, regulon_target_net, nullmodel_60d, verbose = FALSE)

mrs_7d_sum = summary(mrs_7d, mrs = 478)
mrs_28d_sum = summary(mrs_28d, mrs = 478)
mrs_60d_sum = summary(mrs_60d, mrs = 478)

save(mrs_7d_sum, mrs_28d_sum, mrs_60d_sum, file = "06 master regulator identification/mrs_sum.RData")

# relate to other information
load("06 master regulator identification/regulators_anno.RData")

mrs_7d_anno = merge(mrs_7d_sum, regulators_anno, by.x = "Regulon", by.y = "probe_ID")
mrs_28d_anno = merge(mrs_28d_sum, regulators_anno, by.x = "Regulon", by.y = "probe_ID")
mrs_60d_anno = merge(mrs_60d_sum, regulators_anno, by.x = "Regulon", by.y = "probe_ID")

mrs_alld_anno = cbind(mrs_7d_anno$NES, mrs_28d_anno$NES, mrs_60d_anno$NES, mrs_7d_anno[, -c(3:5)])
names(mrs_alld_anno)[1:3] = c("NES_7d", "NES_28d", "NES_60d")

save(mrs_7d_anno, mrs_28d_anno, mrs_60d_anno, mrs_alld_anno, file = "06 master regulator identification/mrs_all_anno.RData")

table(mrs_alld_anno$Modules, mrs_alld_anno$type)

mrs_alld_anno_qc = mrs_alld_anno[!duplicated(mrs_alld_anno$Regulon), ]

# Create a barplot
regulator_distri = mrs_alld_anno_qc %>% group_by(Modules, type) %>% summarise(n = n())

library(plyr)
regulator_distri <- ddply(regulator_distri, "Modules",
                   transform, n_ypos=cumsum(n))


ggplot(data = regulator_distri, aes(x = Modules, y = n, fill = type)) +
  geom_bar(stat = "identity") +
  geom_text(aes(y = n_ypos, label = n), vjust=1.6, 
            color ="white", size = 3.5) +
  scale_fill_brewer(palette = "Paired") +
  xlim("M1", "M2", "M3","M4","M5", "M6", "M7", "M8", "M9", "M10", "M11", "M12", "M13", "M14", "M15") +
  theme_classic()

# draw pheatmap
# creat matrix
nes_matrix = as.matrix(mrs_alld_anno_qc[, 1:3])
rownames(nes_matrix) = mrs_alld_anno_qc$Regulon
nes_matrix = t(nes_matrix)
# generate annotations
annotation_col = data.frame(Modules = factor(mrs_alld_anno_hm$Modules),
                            Type = factor(mrs_alld_anno_hm$type))
rownames(annotation_col) = mrs_alld_anno_hm$Regulon
# specify colors
module_color = mrs_alld_anno[, 9:10][!duplicated(mrs_alld_anno[, 9:10]), ]
m_color = as.character(module_color$Colors); names(m_color) = module_color$Modules
ann_colors = list(Type = c(TF = "#1B9E77", MR = "#D95F02"), 
                  Modules = m_color)
pheatmap(nes_matrix, 
         color = colorRampPalette(rev(brewer.pal(n = 9, name = "RdBu")))(100),
         annotation_col = annotation_col,
         annotation_colors = ann_colors, 
         show_colnames = F,
         treeheight_col = 20,
         cluster_rows = F)


# draw vennplot
top_mrs_7d = mrs_7d_anno[!duplicated(mrs_7d_anno$Regulon), ] %>% filter(abs(NES) >= 2)
top_mrs_28d = mrs_28d_anno[!duplicated(mrs_28d_anno$Regulon), ] %>% filter(abs(NES) >= 2)
top_mrs_60d = mrs_60d_anno[!duplicated(mrs_60d_anno$Regulon), ] %>% filter(abs(NES) >= 2)

save(top_mrs_7d, top_mrs_28d, top_mrs_60d, file = "06 master regulator identification/top_mrs_all.RData")

library(VennDiagram)
overlap <- calculate.overlap(x = list("DE_7id" = as.character(top_mrs_7d$Regulon), "DE_28id" = as.character(top_mrs_28d$Regulon), "DE_60id" = as.character(top_mrs_60d$Regulon)))

venn.plot <- draw.triple.venn(
  area1 = 153,
  area2 = 157,
  area3 = 127,
  n12 = 111,
  n23 = 115,
  n13 = 100,
  n123 = 89,
  lwd = rep(1, 3),
  category = c("7 day", "28 day", "60 day"),
  col = c("chocolate1", "green3", "purple"),
  fill = c("chocolate1", "green3", "purple"),
  alpha = rep(0.05, 3),
  fontfamily = rep("sans", 7),
  cat.fontfamily = rep("sans", 3),
  cex = 1,
  cat.cex = 1.2,
  cat.col = c("chocolate1", "green3", "purple"),
  sep.dist = 0.1)


y = data.frame()
for (i in 1:length(overlap)){
  x = data.frame(probeID = overlap[[i]]); x$group = names(overlap[i])
  y = rbind(y, x)
}

dim(y) # 200   2

# Key regulators (NES>2, <-2) in all three groups
key_regulators = mrs_alld_anno_qc[mrs_alld_anno_qc$Regulon %in% y$probeID, ]
table(key_regulators$Modules, key_regulators$type)



save(key_regulators, file = "06 master regulator identification/key_regulators.RData")



key_regulator_distri = key_regulators %>% group_by(Modules, type) %>% summarise(n = n())

library(plyr)
key_regulator_distri <- ddply(key_regulator_distri, "Modules",
                          transform, n_ypos=cumsum(n))


ggplot(data = key_regulator_distri, aes(x = Modules, y = n, fill = type)) +
  geom_bar(stat = "identity") +
  geom_text(aes(y = n_ypos, label = n), vjust=1.6, 
            color ="white", size = 3.5) +
  scale_fill_brewer(palette = "Paired") +
  xlim("M2", "M3","M6", "M7", "M10", "M12") +
  theme_classic()




# module regulator activity changes during epilepsy progression
regulator_dynamics = mrs_alld_anno_qc %>% group_by(type, Modules) %>% summarize(meanFC_7 = mean(NES_7d), meanFC_28 = mean(NES_28d), meanFC_60 = mean(NES_60d), n = n())
write.csv(regulator_dynamics, file = "06 master regulator identification/regulator_dynamics.csv")


# MR dynamics
mrs_MR = mrs_alld_anno_qc[mrs_alld_anno_qc$type == "MR", ]

MR_meanFC_7 = mrs_MR  %>% group_by(Modules) %>% summarize(NES = mean(NES_7d), sd = sd(NES_7d), n = n())
MR_meanFC_7$Stage = "7day"

MR_meanFC_28 = mrs_MR  %>% group_by(Modules) %>% summarize(NES = mean(NES_28d), sd = sd(NES_28d), n = n())
MR_meanFC_28$Stage = "28day"

MR_meanFC_60 = mrs_MR  %>% group_by(Modules) %>% summarize(NES = mean(NES_60d), sd = sd(NES_60d), n = n())
MR_meanFC_60$Stage = "60day"

MR_FC_dynamics = rbind(MR_meanFC_7, MR_meanFC_28, MR_meanFC_60)

# plot width 600 hight 400
# 
MR_dynamics = ggplot(MR_FC_dynamics[MR_FC_dynamics$Modules %in% c("M2", "M3", "M6", "M7", "M10", "M12"), ], aes(x = Stage, y = NES, group = Modules, color = Modules)) +
  geom_line(size=0.5) +
  geom_point() +
  scale_x_discrete(limits=c("7day", "28day", "60day")) +
  geom_errorbar(aes(ymin=NES-sd, ymax=NES+sd), width=.3,
                position=position_dodge(0.02)) +
  theme_classic() + 
  ggtitle("SPs") + scale_y_continuous(limits=c(-4.1, 4.2), breaks=c(-4, -2, 0, 2, 4))

# TF dynamics
mrs_TF = mrs_alld_anno_qc[mrs_alld_anno_qc$type == "TF", ]

TF_meanFC_7 = mrs_TF  %>% group_by(Modules) %>% summarize(NES = mean(NES_7d), sd = sd(NES_7d), n = n())
TF_meanFC_7$Stage = "7day"

TF_meanFC_28 = mrs_TF  %>% group_by(Modules) %>% summarize(NES = mean(NES_28d), sd = sd(NES_28d), n = n())
TF_meanFC_28$Stage = "28day"

TF_meanFC_60 = mrs_TF  %>% group_by(Modules) %>% summarize(NES = mean(NES_60d), sd = sd(NES_60d), n = n())
TF_meanFC_60$Stage = "60day"

TF_FC_dynamics = rbind(TF_meanFC_7, TF_meanFC_28, TF_meanFC_60)

# plot width 600 hight 400
# 
TF_dynamics = ggplot(TF_FC_dynamics[TF_FC_dynamics$Modules %in% c("M2", "M3", "M6", "M7", "M10", "M12"), ], aes(x = Stage, y = NES, group = Modules, color = Modules)) +
  geom_line(size=0.5) +
  geom_point() +
  scale_x_discrete(limits=c("7day", "28day", "60day")) +
  geom_errorbar(aes(ymin=NES-sd, ymax=NES+sd), width=.3,
                position=position_dodge(0.02)) +
  theme_classic() +
  ggtitle("TFs")

MR_TF_dynamics = ggarrange(MR_dynamics, TF_dynamics, ncol = 2, nrow = 1)
MR_TF_dynamics

# allMRTF_dynamics
all_meanFC_7 = mrs_alld_anno_qc  %>% group_by(Modules) %>% summarize(NES = mean(NES_7d), sd = sd(NES_7d), n = n())
all_meanFC_7$Stage = "7day"

all_meanFC_28 = mrs_alld_anno_qc  %>% group_by(Modules) %>% summarize(NES = mean(NES_28d), sd = sd(NES_28d), n = n())
all_meanFC_28$Stage = "28day"

all_meanFC_60 = mrs_alld_anno_qc  %>% group_by(Modules) %>% summarize(NES = mean(NES_60d), sd = sd(NES_60d), n = n())
all_meanFC_60$Stage = "60day"

allFC_dynamics = rbind(all_meanFC_7, all_meanFC_28, all_meanFC_60)

# plot width 600 hight 400
# 
ggplot(allFC_dynamics[allFC_dynamics$Modules %in% c("M2", "M3", "M6", "M7", "M10", "M12"), ], aes(x = Stage, y = NES, shape = Modules, group = Modules, color = Modules)) +
  geom_line(size=0.5) +
  geom_point() +
  scale_x_discrete(limits=c("7day", "28day", "60day")) +
  geom_errorbar(aes(ymin=NES-sd, ymax=NES+sd), width=.3,
                position=position_dodge(0.02)) +
  theme_classic()
  


save(MR_FC_dynamics, TF_FC_dynamics, file = "06 master regulator identification/regulator_FC_dynamics.RData")





# top 25 regulators in M2,M10 in 60 day
top25_M2M10_60d = filter(mrs_60d_anno, Modules %in% c("M2", "M10"))
# compared to mRNA level changes
load("04 differential expression analysis/DEresults.RData")
top25_mRNA_logFC = DE_alld_anno[match(top25_M2M10_60d$Regulon, DE_alld_anno$probe_ID), c(1, 4, 12)]

top25_M2M10_60d= merge(top25_M2M10_60d, top25_mRNA_logFC, by = "Gene_symbol")

# GSEA analysis for regulators
msig_hallmark = msigdf.mouse %>% filter(category_code=="hallmark")
msig_c5 = msigdf.mouse %>% filter(category_code=="c5")

mrs_60d_anno = mrs_60d_anno[order(mrs_60d_anno$NES, decreasing = T), ]
mrs_60d_NES = mrs_60d_anno$NES; names(mrs_60d_NES) = mrs_60d_anno$mgi_symbol
mrs_60d_gsea = GSEA(mrs_60d_NES, TERM2GENE = msig_c5[, 3:4])
mrs_60d_gsea_re = mrs_60d_gsea@result

save(mrs_60d_gsea, mrs_60d_gsea_re, file = "06 master regulator identification/mrs_60d_gsea.RData")

gsea1 = gseaplot2(mrs_60d_gsea, geneSetID = 23, 
          title = mrs_60d_gsea_re$Description[23], 
          color = "red")

gsea2 = gseaplot2(mrs_60d_gsea, geneSetID = 10, 
          title = mrs_60d_gsea_re$Description[10], 
          color = "red")

gsea3 = gseaplot2(mrs_60d_gsea, geneSetID = 35, 
          title = mrs_60d_gsea_re$Description[35], 
          color = "red")
library(ggpubr)
ggarrange(gsea1, gsea2, gsea3, nrow = 1, ncol = 3)
