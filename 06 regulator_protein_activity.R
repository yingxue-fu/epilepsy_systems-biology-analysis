library(GEOquery)
library(viper)
library(ggplot2)
library(limma)
library(dplyr)
library(ggpubr)


library(clusterProfiler)

load("03 Gene coexpression network/datExpr.RData")
load("06 inferring protein activity_a/regulators_df.RData")

# Regulators and their targets in ARACNe output file
regulator_target_net = read.table(file = "06 inferring protein activity_a/network.txt", header = T)
# 41364 interactions
length(unique(regulator_target_net$Regulator)) # 1493 regulators
length(unique(regulator_target_net$Target)) # 5695 targets
range(regulator_target_net$pvalue)

regulator_target_net$type <- regulators$type[match(regulator_target_net$Regulator, regulators$Gene_symbol)]

table(regulator_target_net$type)

# Signal     SP     TF 
#  11287  21819   8258 

regulator_in_net <- regulator_target_net[!duplicated(regulator_target_net[,c(1,5)]),c(1,5)]

table(regulator_in_net$type)

# Signal     SP     TF 
#   363     696     434 

target_in_net <- regulator_target_net[!duplicated(regulator_target_net[,c(2,5)]),c(2,5)]

table(target_in_net$type)

# Signal     SP     TF 
#  3486     4660   3573 

write.table(regulator_target_net[,-4], file = "06 inferring protein activity_b/network_3col.txt", sep = "\t", row.names = F, col.names = F, quote = F)

# Generating the regulon object form ARACNe output file
regulon = aracne2regulon("06 inferring protein activity_b/network_3col.txt", datExpr, format = "3col")
print(regulon) 
# Object of class regulon with 1307 regulators, 5618 targets and 41178 interactions

regulon_l2d = data.frame()
for (i in 1:1307){
  a = data.frame(regulator = rep(names(regulon[i]), length(regulon[[i]]$likelihood)), target = names(regulon[[i]]$tfmode), tfmode = regulon[[i]]$tfmode, likelihood = regulon[[i]]$likelihood)
  regulon_l2d = rbind(regulon_l2d, a)
}

plot(regulon_l2d$tfmode, regulon_l2d$likelihood)
hist(regulon_l2d$tfmode, breaks = 100)
hist(regulon_l2d$likelihood, breaks = 100)

save(regulon, regulon_l2d, file = "06 inferring protein activity_b/aracne2regulon.RData")

# plot of mode of regulation
ggplot(regulon_l2d, aes(x=tfmode)) + 
  #geom_histogram(aes(y=..density..), colour="gray50", fill="white", binwidth = 0.02)+
  geom_density(alpha=.5, color = "gray50", fill="#FF6666") +
  theme_classic() +
  scale_y_continuous(breaks = c(0, 0.5, 1, 1.5, 2, 2.5))
# plot of interaction confidence
ggplot(regulon_l2d, aes(x=likelihood)) + 
  geom_histogram(aes(y=..density..), colour="gray", fill="white", binwidth = 0.01)+
  geom_density(alpha=.5, color = "gray50", fill="#FF6666") +
  theme_classic()

gset <- ExpressionSet(datExpr)
sample_info_gse73878_IL <- read.csv(file = "00 meta-analysis/sample_info_gse73878_IL.csv")
gset$description <- sample_info_gse73878_IL$Stage

save(gset, file = "06 inferring protein activity_b/gset.RData")

# Master Regulator Analysis performed by msVIPER
## Generating the gene expression signatures
signature_acute <- rowTtest(gset, "description", "acute", "control")
signature_latent <- rowTtest(gset, "description", "latent", "control")
signature_chronic <- rowTtest(gset, "description", "chronic", "control")

# Z-transformation
signature_acute_z <- (qnorm(signature_acute$p.value/2, lower.tail = FALSE) * sign(signature_acute$statistic))[, 1]
signature_latent_z <- (qnorm(signature_latent$p.value/2, lower.tail = FALSE) * sign(signature_latent$statistic))[, 1]
signature_chronic_z <- (qnorm(signature_chronic$p.value/2, lower.tail = FALSE) * sign(signature_chronic$statistic))[, 1]

# compared to RP list
v_signature_acute <- data.frame(Gene_symbol=names(signature_acute_z), zvalue=signature_acute_z)

load("02 Gene meta-signature/metaRP_gene_list.RData")
v_signature_RP_acute <- merge(v_signature_acute, metaRP_acute_list, by.x = "Gene_symbol", by.y = "ID")
plot(v_signature_RP_acute$zvalue, v_signature_RP_acute$logFC)

hist(v_signature_RP_acute$zvalue, breaks = 100)
hist(v_signature_RP_acute$logFC, breaks = 100)
hist(v_signature_RP_acute$NormRank, breaks = 100)

save(v_signature_RP_acute, file = "06 inferring protein activity_b/sample.RData")

## NULL model by sample permutations
nullmodel_acute <- ttestNull(gset, "description", "acute", "control", per = 1000, repos = TRUE, verbose = FALSE)
nullmodel_latent <- ttestNull(gset, "description", "latent", "control", per = 1000, repos = TRUE, verbose = FALSE)
nullmodel_chronic <- ttestNull(gset, "description", "chronic", "control", per = 1000, repos = TRUE, verbose = FALSE)

par(mfrow = c(3,1))
hist(nullmodel_acute[,1], breaks = 100)
hist(nullmodel_latent[,1], breaks = 150)
hist(nullmodel_chronic[,1], breaks = 100)

save(nullmodel_acute, nullmodel_latent, nullmodel_chronic, file = "06 inferring protein activity_b/nullmodel.RData")

## msVIPER
mrs_acute <- msviper(signature_acute_z, regulon, nullmodel_acute, verbose = FALSE)
mrs_latent <- msviper(signature_latent_z, regulon, nullmodel_latent, verbose = FALSE)
mrs_chronic <- msviper(signature_chronic_z, regulon, nullmodel_chronic, verbose = FALSE)

mrs_acute_sum = summary(mrs_acute, mrs = 521)
mrs_latent_sum = summary(mrs_latent, mrs = 521)
mrs_chronic_sum = summary(mrs_chronic, mrs = 521)

a = ggplot(mrs_acute_sum, aes(x=NES, y=p.value)) + 
  geom_point(shape=21,color="black",fill="gray",alpha=0.5) +
  theme_bw() +
  geom_vline(xintercept = c(2,-2),color="red",linetype="dashed") +
  geom_hline(yintercept = 0.05,color="red",linetype="dashed")+
  scale_x_continuous(breaks = c(-4,-2,0,2,4))

l = ggplot(mrs_latent_sum, aes(x=NES, y=p.value)) + 
  geom_point(shape=21,color="black",fill="gray",alpha=0.5) +
  theme_bw() +
  geom_vline(xintercept = c(2,-2),color="red",linetype="dashed") +
  geom_hline(yintercept = 0.05,color="red",linetype="dashed")+
  scale_x_continuous(breaks = c(-4,-2,0,2,4))

c = ggplot(mrs_chronic_sum, aes(x=NES, y=p.value)) + 
  geom_point(shape=21,color="black",fill="gray",alpha=0.5) +
  theme_bw() +
  geom_vline(xintercept = c(2,-2),color="red",linetype="dashed") +
  geom_hline(yintercept = 0.05,color="red",linetype="dashed")+
  scale_x_continuous(breaks = c(-4,-2,0,2,4))

a1 = ggplot(mrs_acute_sum, aes(x=NES)) + 
  geom_histogram(color="black", fill="gray",alpha=0.5) +
  scale_x_continuous(breaks = c(-4,-2,0,2,4)) +
  geom_vline(xintercept = c(2,-2),color="red",linetype="dashed") +
  theme_bw()

l1 = ggplot(mrs_latent_sum, aes(x=NES)) + 
  geom_histogram(color="black", fill="gray",alpha=0.5) +
  scale_x_continuous(breaks = c(-4,-2,0,2,4)) +
  geom_vline(xintercept = c(2,-2),color="red",linetype="dashed") +
  theme_bw()

c1 = ggplot(mrs_chronic_sum, aes(x=NES)) + 
  geom_histogram(color="black", fill="gray",alpha=0.5) +
  scale_x_continuous(breaks = c(-4,-2,0,2,4)) +
  geom_vline(xintercept = c(2,-2),color="red",linetype="dashed") +
  theme_bw()

a2 = ggplot(mrs_acute_sum, aes(x=p.value)) + 
  geom_histogram(color="black", fill="gray",alpha=0.5) +
  geom_vline(xintercept = 0.05,color="red",linetype="dashed") +
  theme_bw()

l2 = ggplot(mrs_latent_sum, aes(x=p.value)) + 
  geom_histogram(color="black", fill="gray",alpha=0.5) +
  geom_vline(xintercept = 0.05,color="red",linetype="dashed") +
  theme_bw()

c2 = ggplot(mrs_chronic_sum, aes(x=p.value)) + 
  geom_histogram(color="black", fill="gray",alpha=0.5) +
  geom_vline(xintercept = 0.05,color="red",linetype="dashed") +
  theme_bw()

ggarrange(a, a1, a2, l, l1,l2, c,c1,c2, ncol = 3, nrow = 3)

plot(mrs_chronic, mrs = c("Arr3", "Scgn"), cex = .9)
plot(mrs_chronic, mrs = 25, cex = .7)

save(mrs_acute, mrs_latent, mrs_chronic, mrs_acute_sum, mrs_latent_sum, mrs_chronic_sum, file = "06 inferring protein activity_b/mrs_sum.RData")

# relate to other information
load("06 inferring protein activity_a/regulators_module.RData")

mrs_acute_anno = merge(mrs_acute_sum, regulators_module, by.x = "Regulon", by.y = "Gene_symbol")
mrs_latent_anno = merge(mrs_latent_sum, regulators_module, by.x = "Regulon", by.y = "Gene_symbol")
mrs_chronic_anno = merge(mrs_chronic_sum, regulators_module, by.x = "Regulon", by.y = "Gene_symbol")

mrs_alld_anno = cbind(mrs_acute_anno$NES, mrs_latent_anno$NES, mrs_chronic_anno$NES, mrs_acute_anno[, -c(3:5)])
names(mrs_alld_anno)[1:3] = c("NES_acute", "NES_latent", "NES_chronic")

save(mrs_acute_anno, mrs_latent_anno, mrs_chronic_anno, mrs_alld_anno, file = "06 inferring protein activity_b/mrs_all_anno.RData")

table(mrs_alld_anno$Modules, mrs_alld_anno$type)

write.table(mrs_alld_anno$Regulon, file = "06 inferring protein activity_b/final_regulator.txt", col.names = F, row.names = F, quote = F)

############################### DE regulators #############################
DE_mrs_acute = mrs_acute_anno %>% filter(abs(NES) >= 2 & p.value < 0.05)
DE_mrs_latent = mrs_latent_anno %>% filter(abs(NES) >= 2 & p.value < 0.05)
DE_mrs_chronic = mrs_chronic_anno %>% filter(abs(NES) >= 2 & p.value < 0.05)

# draw vennplot
library(VennDiagram)
overlap <- calculate.overlap(list("acute" = as.character(DE_mrs_acute$Regulon), "latent" = as.character(DE_mrs_latent$Regulon), "chronic" = as.character(DE_mrs_chronic$Regulon)))

venn.plot <- draw.triple.venn(
  area1 = 214,
  area2 = 198,
  area3 = 156,
  n12 = 148,
  n23 = 149,
  n13 = 119,
  n123 = 113,
  lwd = rep(1, 3),
  category = c("acute", "latent", "chronic"),
  col = c("#c00000", "#0070c0", "#00b050"),
  fill = c("#c00000", "#0070c0", "#00b050"),
  alpha = rep(0.1, 3),
  fontfamily = rep("sans", 7),
  cat.fontfamily = rep("sans", 3),
  cex = 1,
  cat.cex = 1.2,
  cat.col = c("#c00000", "#0070c0", "#00b050"),
  sep.dist = 0.1)

y = data.frame()
for (i in 1:length(overlap)){
  x = data.frame(Gene_symbol = overlap[[i]]); x$group = names(overlap[i])
  y = rbind(y, x)
}

dim(y) # 265   2

# DE regulators in all three groups
DE_regulators = mrs_alld_anno[mrs_alld_anno$Regulon %in% y$Gene_symbol, ]
table(DE_regulators$Modules, DE_regulators$type)

save(DE_mrs_acute, DE_mrs_latent, DE_mrs_chronic, DE_regulators, file = "06 inferring protein activity_b/DE_regulators.RData")

write.csv(DE_regulators, file = "06 inferring protein activity_b/DE_regulators.csv")
write.table(DE_regulators$Regulon, file = "06 inferring protein activity_b/DE_regulators.txt", quote = F, col.names = F, row.names = F)

## DE regulators in each module #############################################
DE_regulators$type <- factor(DE_regulators$type, levels = c("SP","Signal","TF"))
DE_regulator_distri = DE_regulators %>% group_by(Modules, type) %>% summarise(n = n())

ggplot(data = DE_regulator_distri, aes(x = Modules, y = n, fill = type)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = n), position = position_stack(vjust = 0.7), color ="white", size = 2) +
  scale_fill_brewer(palette = "Paired") +
  xlim("M1","M2", "M3","M4", "M6", "M7","M8","M9","M10","M11") +
  theme_classic()


## module regulator activity changes during epilepsy progression
DE_acute_meanFC = DE_regulators %>% group_by(Modules) %>% summarize(NES = mean(NES_acute), sd = sd(NES_acute), n = n())
DE_acute_meanFC$Stage = "acute"

DE_latent_meanFC = DE_regulators  %>% group_by(Modules) %>% summarize(NES = mean(NES_latent), sd = sd(NES_latent), n = n())
DE_latent_meanFC$Stage = "latent"

DE_chronic_meanFC = DE_regulators  %>% group_by(Modules) %>% summarize(NES = mean(NES_chronic), sd = sd(NES_chronic), n = n())
DE_chronic_meanFC$Stage = "chronic"

DE_FC_dynamics = rbind(DE_acute_meanFC, DE_latent_meanFC, DE_chronic_meanFC) %>% na.omit()

ggplot(DE_FC_dynamics[DE_FC_dynamics$Modules %in% c("M1", "M3", "M6", "M7","M8"), ], aes(x = Stage, y = NES, group = Modules, color = Modules)) +
  geom_line(size = 0.75) +
  geom_point(size = 1.2) +
  scale_x_discrete(limits=c("acute", "latent", "chronic")) +
  geom_errorbar(aes(ymin=NES-sd, ymax=NES+sd), width=0.5,
                position=position_dodge(0.01)) +
  theme_classic()+
  scale_color_manual(values = c("#8DD3C7","#BEBADA","#B3DE69","#FCCDE5","#D9D9D9"))
  
save(DE_regulator_distri, DE_FC_dynamics, file = "06 inferring protein activity_b/DEregulator_distri_dynamics.RData")


############################## heatmap for DE regulators #############################
library(pheatmap)
library(RColorBrewer)

# creat matrix
DE_nes_matrix = as.matrix(DE_regulators[,1:3])
rownames(DE_nes_matrix) <- DE_regulators$Regulon

# generate annotations
annotation_row = data.frame(Modules = factor(DE_regulators$Modules),
                            Type = factor(DE_regulators$type))
rownames(annotation_row) = DE_regulators$Regulon
# specify colors
module_color = DE_regulators[, 7:8][!duplicated(DE_regulators[, 7:8]), ]
m_color = as.character(module_color$Colors); names(m_color) = module_color$Modules
ann_colors = list(Type = c(SP = "#66A61E", Signal = "#E6AB02", TF = "#A6761D"),
                  Modules = m_color)


pheatmap(DE_nes_matrix, 
         color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100),
         annotation_row = annotation_row,
         annotation_colors = ann_colors, 
         treeheight_col = 10,
         treeheight_row = 15,
         show_rownames = F) 

###################### KEGG enrichment analysis for DE regulators #########################
library(clusterProfiler)
load("06 inferring protein activity_a/mouse_maRt.RData")

DEregulators_entrezID <- getBM(attributes = c("mgi_symbol","entrezgene"),filters = "mgi_symbol",
                               values = DE_regulators$Regulon, mart = mouse_maRt)

DE_regulator_kegg <- enrichKEGG(gene = DEregulators_entrezID$entrezgene,
                 organism = 'mmu',
                 pvalueCutoff = 0.05)

DE_regulator_kegg_df <- DE_regulator_kegg@result

DE_regulator_kegg_0.05 <- DE_regulator_kegg_df %>% filter(p.adjust<0.05)

###########################################################################
top15kegg <- read.csv("06 inferring protein activity_b/DE_regulator_KEGG - selected.csv")

ggplot(data = top15kegg, aes(x=Term,y=-log10(Benjamini))) +
  geom_bar(stat="identity", fill="gray", width = 0.7) +
  xlim(as.character(top15kegg$Term)) +
  theme_classic() +
  coord_flip() 


###########################################################################
# relations to meta-DEG
mrs_alld_anno$acute_logFC <- metaRP_acute_list$logFC[match(mrs_alld_anno$Regulon, metaRP_acute_list$ID)]

mrs_alld_anno$latent_logFC <- metaRP_latent_list$logFC[match(mrs_alld_anno$Regulon, metaRP_latent_list$ID)]

mrs_alld_anno$chronic_logFC <- metaRP_chronic_list$logFC[match(mrs_alld_anno$Regulon, metaRP_chronic_list$ID)]

acute_lm <- mrs_alld_anno[,c("NES_acute","acute_logFC")] %>% na.omit()
a_lm_eq <- lm_eqn(acute_lm, 'acute_logFC','NES_acute')
a <- ggplot(acute_lm, aes(x=NES_acute, y=acute_logFC))+
  geom_vline(xintercept = 0,linetype="dashed") +
  geom_hline(yintercept = 0,linetype="dashed")+
  geom_point(alpha=0.5,color="gray30") +
  geom_smooth(method=lm, se=FALSE,color="red",formula = y~x,show.legend = TRUE) +
  theme_bw()+
  annotate("text", x=-0.5,y=2,label=a_lm_eq,color='red',parse=T,size=3)

latent_lm <- mrs_alld_anno[,c("NES_latent","latent_logFC")] %>% na.omit()
l_lm_eq <- lm_eqn(latent_lm, 'latent_logFC','NES_latent')
l <- ggplot(latent_lm, aes(x=NES_latent, y=latent_logFC))+
  geom_vline(xintercept = 0,linetype="dashed") +
  geom_hline(yintercept = 0,linetype="dashed")+
  geom_point(alpha=0.5,color="gray30") +
  geom_smooth(method=lm, se=FALSE,color="red") +
  theme_bw()+
  annotate("text", x=0,y=2,label=l_lm_eq,color='red',parse=T,size=3)

chronic_lm <- mrs_alld_anno[,c("NES_chronic","chronic_logFC")] %>% na.omit()
c_lm_eq <- lm_eqn(chronic_lm, 'chronic_logFC','NES_chronic')
c <- ggplot(chronic_lm, aes(x=NES_chronic, y=chronic_logFC))+
  geom_vline(xintercept = 0,linetype="dashed") +
  geom_hline(yintercept = 0,linetype="dashed")+
  geom_point(alpha=0.5,color="gray30") +
  geom_smooth(method=lm, se=FALSE,color="red") +
  theme_bw()+
  annotate("text", x=0,y=1.2,label=c_lm_eq,color='red',parse=T,size=3)

library(ggpubr)
ggarrange(a,l,c, nrow = 1,ncol = 3)

lm_eqn <- function(df, y, x){
  formula = as.formula(sprintf('%s ~ %s', y, x))
  m <- lm(formula, data=df);
  # formating the values into a summary string to print out
  # ~ give some space, but equal size and comma need to be quoted
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2*","~~italic(p)~"="~pvalue, 
                   list(a = format(as.vector(coef(m)[1]), digits = 2), 
                        b = format(as.vector(coef(m)[2]), digits = 2), 
                        r2 = format(summary(m)$r.squared, digits = 2),
                        # getting the pvalue is painful
                        pvalue = format(summary(m)$coefficients[2,'Pr(>|t|)'], digits=2))
                   )
  as.character(as.expression(eq));                 
}



cor.test(x=acute_lm$NES_acute,y=acute_lm$acute_logFC,method = "spearman")
cor.test(x=latent_lm$NES_latent,y=latent_lm$latent_logFC,method = "spearman")
cor.test(x=chronic_lm$NES_chronic,y=chronic_lm$chronic_logFC,method = "spearman")


