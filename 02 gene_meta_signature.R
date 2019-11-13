library(Biobase)
library(GEOquery)
library(RankProd)
library(dplyr)

load("00 meta-analysis/meta_gse.RData")
load("00 meta-analysis/meta_sample_info.RData")
load("01 differential expression-limma/gpl.RData")

################################### GPLs ################################################
library(VennDiagram)
gpl1261_gene <- unique(gpl1261$`Gene symbol`)
gpl2896_gene <- unique(gpl2896$`Gene symbol`)
gpl6247_gene <- unique(gpl6247$`Gene symbol`)
gpl6885_gene <- unique(gpl6885$`Gene symbol`)

gpl_overlap <- calculate.overlap(list("gpl1261" = gpl1261_gene, "gpl2896" = gpl2896_gene, "gpl6247" = gpl6247_gene, "gpl6885" = gpl6885_gene))

write.csv(gpl1261_gene, file = "02 Gene meta-signature/gpl1261_gene.csv")
write.csv(gpl2896_gene, file = "02 Gene meta-signature/gpl2896_gene.csv")
write.csv(gpl6247_gene, file = "02 Gene meta-signature/gpl6247_gene.csv")
write.csv(gpl6885_gene, file = "02 Gene meta-signature/gpl6885_gene.csv")

core_genes <- data.frame(Gene_symbol = gpl_overlap$a6)

save(gpl_overlap, core_genes, gpl1261_gene, gpl2896_gene, gpl6247_gene, gpl6885_gene, file = "02 Gene meta-signature/gpl_overlap.RData")

################## extract common genes in each GSE dataset ###########################
# GSE14763
exprs_gse14763 <- na.omit(log2(exprs(gse14763)+1)) %>% as.data.frame()
exprs_gse14763$ID <- rownames(exprs_gse14763)
exprs_gse14763 <- merge(exprs_gse14763, gpl2896[,c(1,3)], by = "ID")

exprs_gse14763$mean = apply(exprs_gse14763[,-c(1,20)], 1, mean)

exprs_gse14763_qc <- exprs_gse14763 %>% group_by(`Gene symbol`) %>% top_n(1, mean)

exprs_gse14763_core <- exprs_gse14763_qc[match(core_genes$Gene_symbol, exprs_gse14763_qc$`Gene symbol`) %>% na.omit(),-c(1,21,22)] 

# GSE27166
exprs_gse27166 <- na.omit(log2(exprs(gse27166)+1)) %>% as.data.frame()
exprs_gse27166$ID <- rownames(exprs_gse27166)
exprs_gse27166 <- merge(exprs_gse27166, gpl2896[,c(1,3)], by = "ID")

exprs_gse27166$mean = apply(exprs_gse27166[,-c(1,26)], 1, mean)

exprs_gse27166_qc <- exprs_gse27166 %>% group_by(`Gene symbol`) %>% top_n(1, mean)

exprs_gse27166_core <- exprs_gse27166_qc[match(core_genes$Gene_symbol, exprs_gse27166_qc$`Gene symbol`) %>% na.omit(),-c(1,27)] 

# GSE49849
exprs_gse49849 <- na.omit(exprs(gse49849)) %>% as.data.frame()
exprs_gse49849$ID <- rownames(exprs_gse49849)
exprs_gse49849 <- merge(exprs_gse49849, gpl6247[,c(1,3)], by = "ID")

exprs_gse49849$mean = apply(exprs_gse49849[,-c(1,22)], 1, mean)

exprs_gse49849_qc <- exprs_gse49849 %>% group_by(`Gene symbol`) %>% top_n(1, mean)

exprs_gse49849_core <- exprs_gse49849_qc[match(core_genes$Gene_symbol, exprs_gse49849_qc$`Gene symbol`) %>% na.omit(),-c(1,23)] 

# GSE73878
exprs_gse73878 <- na.omit(exprs(gse73878)) %>% as.data.frame()
exprs_gse73878$ID <- rownames(exprs_gse73878)
exprs_gse73878 <- merge(exprs_gse73878, gpl6885[,c(1,3)], by = "ID")

exprs_gse73878$mean = apply(exprs_gse73878[,-c(1,82)], 1, mean)

exprs_gse73878_qc <- exprs_gse73878 %>% group_by(`Gene symbol`) %>% top_n(1, mean)

exprs_gse73878_core <- exprs_gse73878_qc[match(core_genes$Gene_symbol, exprs_gse73878_qc$`Gene symbol`) %>% na.omit(),-c(1,83)] 


# GSE88992
exprs_gse88992 <- na.omit(exprs(gse88992)) %>% as.data.frame()
exprs_gse88992$ID <- rownames(exprs_gse88992)
exprs_gse88992 <- merge(exprs_gse88992, gpl1261[,c(1,3)], by = "ID")

exprs_gse88992$mean = apply(exprs_gse88992[,-c(1,19)], 1, mean)

exprs_gse88992_qc <- exprs_gse88992 %>% group_by(`Gene symbol`) %>% top_n(1, mean)

exprs_gse88992_core <- exprs_gse88992_qc[match(core_genes$Gene_symbol, exprs_gse88992_qc$`Gene symbol`) %>% na.omit(),-c(1,20)] 

save(exprs_gse14763_core, exprs_gse27166_core, exprs_gse49849_core, exprs_gse73878_core, exprs_gse88992_core, file = "02 Gene meta-signature/exprs_gse_core.RData")

###########################################################################################

fl_14763 <- as.factor(sample_info_gse14763$Stage)
fl_27166 <- as.factor(sample_info_gse27166_right$Stage)
fl_49849 <- as.factor(sample_info_gse49849$Stage)
fl_73878 <- as.factor(sample_info_gse73878_IL$Stage)
fl_88992 <- as.factor(sample_info_gse88992$Stage)

####################################### acute meta ##########################################
# acute 14763
exprs_gse14763_acute_core <- exprs_gse14763_core[,c(which(fl_14763==c("acute")), which(fl_14763==c("control")), 19)]
cl_a1 <- rep(c(0,1),c(10,3))

# acute 49849
exprs_gse49849_acute_core <- exprs_gse49849_core[,c(which(fl_49849==c("acute")), which(fl_49849==c("control")),21)]
cl_a2 <- rep(c(0,1),c(5,10))

# acute 73878
exprs_gse73878_core_IL <- exprs_gse73878_core[ ,c(gse73878$`hemisphere:ch1` == "IL",TRUE)]
exprs_gse73878_core_IL <- exprs_gse73878_core_IL[,-c(16,11,29,30,31,32)]

save(exprs_gse73878_core_IL, file = "02 Gene meta-signature/exprs_gse73878_core_IL.RData")

exprs_gse73878_acute_core <- exprs_gse73878_core_IL[,c(which(fl_73878==c("acute")), which(fl_73878==c("control")),35)]
cl_a3 <- rep(c(0,1),c(7,11))

# acute 88992
exprs_gse88992_acute_core <- exprs_gse88992_core
cl_a4 <- rep(c(1,0,1,0,1,0),c(3,3,3,3,3,2))

# merge all datasets

acute_core <- merge(exprs_gse14763_acute_core, exprs_gse49849_acute_core, by = "Gene symbol")
acute_core <- merge(acute_core, exprs_gse73878_acute_core, by = "Gene symbol")
acute_core <- merge(acute_core, exprs_gse88992_acute_core, by = "Gene symbol")

acute_core_matrix <- acute_core[,-1] %>% as.matrix()
rownames(acute_core_matrix) <- acute_core$`Gene symbol`

cl_acute <- c(cl_a1, cl_a2, cl_a3, cl_a4)
origin_acute <- c(rep(1, length(cl_a1)), rep(2, length(cl_a2)), rep(3, length(cl_a3)), rep(4, length(cl_a4)))

metaRP_acute <- RP.advance(acute_core_matrix, cl_acute, origin_acute, logged=TRUE, gene.names = rownames(acute_core_matrix), rand = 123)

plotRP(metaRP_acute, cutoff=0.05)

topgene_acute <- topGene(metaRP_acute, cutoff=0.05, method="pfp", logged=TRUE, logbase=2, gene.names=rownames(acute_core_matrix))

metaRP_acute_df <- RPlist2df(metaRP_acute)

save(acute_core, acute_core_matrix, cl_acute, origin_acute, metaRP_acute, topgene_acute, metaRP_acute_df, file = "02 Gene meta-signature/acute_meta.RData")


#################################### latent meta #############################################

# latent 14763
exprs_gse14763_latent_core <- exprs_gse14763_core[,c(which(fl_14763==c("latent")), which(fl_14763==c("control")), 19)]
cl_l1 <- rep(c(0,1),c(5,3))

# latent 49849
exprs_gse49849_latent_core <- exprs_gse49849_core[,c(which(fl_49849==c("latent")), which(fl_49849==c("control")),21)]
cl_l2 <- rep(c(0,1),c(5,10))

# latent 73878
exprs_gse73878_latent_core <- exprs_gse73878_core_IL[,c(which(fl_73878==c("latent")), which(fl_73878==c("control")),35)]
cl_l3 <- rep(c(0,1),c(8,11))


# merge all datasets

latent_core <- merge(exprs_gse14763_latent_core, exprs_gse49849_latent_core, by = "Gene symbol")
latent_core <- merge(latent_core, exprs_gse73878_latent_core, by = "Gene symbol")


latent_core_matrix <- latent_core[,-1] %>% as.matrix()
rownames(latent_core_matrix) <- latent_core$`Gene symbol`

cl_latent <- c(cl_l1, cl_l2, cl_l3)
origin_latent <- c(rep(1, length(cl_l1)), rep(2, length(cl_l2)), rep(3, length(cl_l3)))

metaRP_latent <- RP.advance(latent_core_matrix, cl_latent, origin_latent, logged=TRUE, gene.names = rownames(latent_core_matrix), rand = 123)

plotRP(metaRP_latent, cutoff=0.05)

topgene_latent <- topGene(metaRP_latent, cutoff=0.05, method="pfp", logged=TRUE, logbase=2, gene.names=rownames(latent_core_matrix))

metaRP_latent_df <- RPlist2df(metaRP_latent)

save(latent_core, latent_core_matrix, cl_latent, origin_latent, metaRP_latent, topgene_latent, metaRP_latent_df, file = "02 Gene meta-signature/latent_meta.RData")

#################################### chronic meta #############################################

# chronic 27166
exprs_gse27166_chronic_core <- exprs_gse27166_core[, c(gse27166$`side:ch1`=="right", TRUE)]
exprs_gse27166_chronic_core <- exprs_gse27166_chronic_core[, -c(7,12)]
cl_c1 <- rep(c(0,1,0,1),c(3,3,2,2))

# chronic 73878
exprs_gse73878_chronic_core <- exprs_gse73878_core_IL[,c(which(fl_73878==c("chronic")), which(fl_73878==c("control")),35)]
cl_c2 <- rep(c(0,1),c(8,11))

# merge all datasets

chronic_core <- merge(exprs_gse27166_chronic_core, exprs_gse73878_chronic_core, by = "Gene symbol")

chronic_core_matrix <- chronic_core[,-1] %>% as.matrix()
rownames(chronic_core_matrix) <- chronic_core$`Gene symbol`

cl_chronic <- c(cl_c1, cl_c2)
origin_chronic <- c(rep(1, length(cl_c1)), rep(2, length(cl_c2)))

metaRP_chronic <- RP.advance(chronic_core_matrix, cl_chronic, origin_chronic, logged=TRUE, gene.names = rownames(chronic_core_matrix), rand = 123)

plotRP(metaRP_chronic, cutoff=0.05)

topgene_chronic <- topGene(metaRP_chronic, cutoff=0.05, method="pfp", logged=TRUE, logbase=2, gene.names=rownames(chronic_core_matrix))

metaRP_chronic_df <- RPlist2df(metaRP_chronic)

save(chronic_core, chronic_core_matrix, cl_chronic, origin_chronic, metaRP_chronic, topgene_chronic, metaRP_chronic_df, file = "02 Gene meta-signature/chronic_meta.RData")

save(metaRP_acute_df, metaRP_chronic_df, metaRP_latent_df, file = "02 Gene meta-signature/metaRP_df.RData")

##################################### gene list ########################################
library(dplyr)
# acute 
metaRP_acute_df_up <- metaRP_acute_df %>% filter(logFC > 0) 
metaRP_acute_df_up$NormRank <- 1-((metaRP_acute_df_up$RPrank_up-1)/max(metaRP_acute_df_up$RPrank_up))

metaRP_acute_df_down <- metaRP_acute_df %>% filter(logFC < 0)
metaRP_acute_df_down$NormRank <- -1+((metaRP_acute_df_down$RPrank_down-1)/max(metaRP_acute_df_down$RPrank_down))

metaRP_acute_list <- rbind(metaRP_acute_df_up[,c(10,9,11)], metaRP_acute_df_down[,c(10,9,11)])
metaRP_acute_list <- metaRP_acute_list[order(metaRP_acute_list$NormRank, decreasing = T),]

# latent 
metaRP_latent_df_up <- metaRP_latent_df %>% filter(logFC > 0) 
metaRP_latent_df_up$NormRank <- 1-((metaRP_latent_df_up$RPrank_up-1)/max(metaRP_latent_df_up$RPrank_up))

metaRP_latent_df_down <- metaRP_latent_df %>% filter(logFC < 0)
metaRP_latent_df_down$NormRank <- -1+((metaRP_latent_df_down$RPrank_down-1)/max(metaRP_latent_df_down$RPrank_down))

metaRP_latent_list <- rbind(metaRP_latent_df_up[,c(10,9,11)], metaRP_latent_df_down[,c(10,9,11)])
metaRP_latent_list <- metaRP_latent_list[order(metaRP_latent_list$NormRank, decreasing = T),]

# chronic 
metaRP_chronic_df_up <- metaRP_chronic_df %>% filter(logFC > 0) 
metaRP_chronic_df_up$NormRank <- 1-((metaRP_chronic_df_up$RPrank_up-1)/max(metaRP_chronic_df_up$RPrank_up))

metaRP_chronic_df_down <- metaRP_chronic_df %>% filter(logFC < 0)
metaRP_chronic_df_down$NormRank <- -1+((metaRP_chronic_df_down$RPrank_down-1)/max(metaRP_chronic_df_down$RPrank_down))

metaRP_chronic_list <- rbind(metaRP_chronic_df_up[,c(10,9,11)], metaRP_chronic_df_down[,c(10,9,11)])
metaRP_chronic_list <- metaRP_chronic_list[order(metaRP_chronic_list$NormRank, decreasing = T),]

save(metaRP_acute_list, metaRP_latent_list, metaRP_chronic_list, file = "02 Gene meta-signature/metaRP_gene_list.RData")


aR <- ggplot(metaRP_acute_list, aes(x=NormRank, y=logFC)) + 
  geom_point(color="gray",size=0.5,alpha=0.8) +
  geom_bin2d(bins=c(200,75)) + 
  theme_bw()

lR <- ggplot(metaRP_latent_list, aes(x=NormRank, y=logFC)) + 
  geom_point(color="gray",size=0.5,alpha=0.8) +
  geom_bin2d(bins=c(200,75)) + 
  theme_bw()

cR <- ggplot(metaRP_chronic_list, aes(x=NormRank, y=logFC)) + 
  geom_point(color="gray",size=0.5,alpha=0.8) +
  geom_bin2d(bins=c(200,75)) + 
  theme_bw()

ggarrange(aR, lR, cR, ncol = 1, nrow = 3)

################################# top 100 genes ######################################
top100_acute = rbind(head(metaRP_acute_list, n = 50L), tail(metaRP_acute_list, n = 50L))
top100_latent = rbind(head(metaRP_latent_list, n = 50L), tail(metaRP_latent_list, n = 50L))
top100_chronic = rbind(head(metaRP_chronic_list, n = 50L), tail(metaRP_chronic_list, n = 50L))

top100.overlap <- calculate.overlap(list("acute" = top100_acute$ID, "latent" = top100_latent$ID, "chronic" = top100_chronic$ID))

top100.overlap.up <- calculate.overlap(list("acute" = top100_acute$ID[1:50], "latent" = top100_latent$ID[1:50], "chronic" = top100_chronic$ID[1:50]))
top100.overlap.down <- calculate.overlap(list("acute" = top100_acute$ID[51:100], "latent" = top100_latent$ID[51:100], "chronic" = top100_chronic$ID[51:100]))

save()

# acute
RP_gse14763_acute_df_cq <- RP_gse14763_acute_df %>% group_by(`Gene symbol`) %>% top_n(1, abs(logFC))
top100_acute_14763_fc <- RP_gse14763_acute_df_cq[match(top100_acute$ID, RP_gse14763_acute_df_cq$`Gene symbol`), c(10,12)]

RP_gse49849_acute_df_cq <- RP_gse49849_acute_df %>% group_by(`Gene symbol`) %>% top_n(1, abs(logFC))
top100_acute_49849_fc <- RP_gse49849_acute_df_cq[match(top100_acute$ID, RP_gse49849_acute_df_cq$`Gene symbol`), c(10,12)]

RP_gse73878_acute_df_cq <- RP_gse73878_acute_df %>% group_by(`Gene symbol`) %>% top_n(1, abs(logFC))
top100_acute_73878_fc <- RP_gse73878_acute_df_cq[match(top100_acute$ID, RP_gse73878_acute_df_cq$`Gene symbol`), c(10,12)]

RP_gse88992_acute_df_cq <- RP_gse88992_acute_df %>% group_by(`Gene symbol`) %>% top_n(1, abs(logFC))
top100_acute_88992_fc <- RP_gse88992_acute_df_cq[match(top100_acute$ID, RP_gse88992_acute_df_cq$`Gene symbol`), c(10,12)]

top100_acute_fcs <- data.frame(Gene_symbol = top100_acute$ID, GSE14763 = top100_acute_14763_fc$logFC, GSE49849 = top100_acute_49849_fc$logFC, GSE73878 = top100_acute_73878_fc$logFC, GSE88992 = top100_acute_88992_fc$logFC)

top100_acute_fcs_matrix <- as.matrix(top100_acute_fcs[,2:5]); rownames(top100_acute_fcs_matrix) <- top100_acute_fcs$Gene_symbol

top100_acute_fcs_matrix_norm <- rbind(top100_acute_fcs_matrix[1:50,]/5.13, top100_acute_fcs_matrix[51:100,]/2.15)

library(pheatmap); library(RColorBrewer)
pheatmap(top100_acute_fcs_matrix_norm, 
         color = colorRampPalette(rev(brewer.pal(n = 9, name = "RdBu")))(100),
         border_color = NA,
         cluster_rows = F,
         cluster_cols = F)


# latent
RP_gse14763_latent_df_cq <- RP_gse14763_latent_df %>% group_by(`Gene symbol`) %>% top_n(1, abs(logFC))
top100_latent_14763_fc <- RP_gse14763_latent_df_cq[match(top100_latent$ID, RP_gse14763_latent_df_cq$`Gene symbol`), c(10,12)]

RP_gse49849_latent_df_cq <- RP_gse49849_latent_df %>% group_by(`Gene symbol`) %>% top_n(1, abs(logFC))
top100_latent_49849_fc <- RP_gse49849_latent_df_cq[match(top100_latent$ID, RP_gse49849_latent_df_cq$`Gene symbol`), c(10,12)]

RP_gse73878_latent_df_cq <- RP_gse73878_latent_df %>% group_by(`Gene symbol`) %>% top_n(1, abs(logFC))
top100_latent_73878_fc <- RP_gse73878_latent_df_cq[match(top100_latent$ID, RP_gse73878_latent_df_cq$`Gene symbol`), c(10,12)]

top100_latent_fcs <- data.frame(Gene_symbol = top100_latent$ID, GSE14763 = top100_latent_14763_fc$logFC, GSE49849 = top100_latent_49849_fc$logFC, GSE73878 = top100_latent_73878_fc$logFC)

top100_latent_fcs_matrix <- as.matrix(top100_latent_fcs[,2:4]); rownames(top100_latent_fcs_matrix) <- top100_latent_fcs$Gene_symbol

top100_latent_fcs_matrix_norm <- rbind(top100_latent_fcs_matrix[1:50,]/5.21, top100_latent_fcs_matrix[51:100,]/1.34)

pheatmap(top100_latent_fcs_matrix_norm, 
         color = colorRampPalette(rev(brewer.pal(n = 9, name = "RdBu")))(100),
         border_color = NA,
         cluster_rows = F,
         cluster_cols = F)

# chronic
RP_gse27166_chronic_df_cq <- RP_gse27166_chronic_df %>% group_by(`Gene symbol`) %>% top_n(1, abs(logFC))
top100_chronic_27166_fc <- RP_gse27166_chronic_df_cq[match(top100_chronic$ID, RP_gse27166_chronic_df_cq$`Gene symbol`), c(10,12)]

RP_gse73878_chronic_df_cq <- RP_gse73878_chronic_df %>% group_by(`Gene symbol`) %>% top_n(1, abs(logFC))
top100_chronic_73878_fc <- RP_gse73878_chronic_df_cq[match(top100_chronic$ID, RP_gse73878_chronic_df_cq$`Gene symbol`), c(10,12)]

top100_chronic_fcs <- data.frame(Gene_symbol = top100_chronic$ID, GSE27166 = top100_chronic_27166_fc$logFC, GSE73878 = top100_chronic_73878_fc$logFC)

top100_chronic_fcs_matrix <- as.matrix(top100_chronic_fcs[,2:3]); rownames(top100_chronic_fcs_matrix) <- top100_chronic_fcs$Gene_symbol

top100_chronic_fcs_matrix_norm <- rbind(top100_chronic_fcs_matrix[1:50,]/3.87, top100_chronic_fcs_matrix[51:100,]/1.05)

pheatmap(top100_chronic_fcs_matrix_norm, 
         color = colorRampPalette(rev(brewer.pal(n = 9, name = "RdBu")))(100),
         border_color = NA,
         cluster_rows = F,
         cluster_cols = F)

save(top100_acute, top100_acute_14763_fc, top100_acute_49849_fc, top100_acute_73878_fc, top100_acute_88992_fc, top100_acute_fcs, top100_chronic_fcs_matrix, top100_acute_fcs_matrix_norm, file = "02 Gene meta-signature/top100_acute.RData")

save(top100_chronic, top100_chronic_27166_fc, top100_chronic_73878_fc, top100_chronic_fcs, top100_chronic_fcs_matrix, top100_chronic_fcs_matrix_norm, top100_chronic, file = "02 Gene meta-signature/top100_chronic.RData")

save(top100_latent, top100_latent_14763_fc, top100_latent_49849_fc, top100_latent_73878_fc, top100_latent_fcs, top100_latent_fcs_matrix, top100_latent_fcs_matrix_norm, file = "02 Gene meta-signature/top100_latent.RData")

################################## meta DEGs ########################################

metaRP_acute_DEG <- metaRP_acute_df %>% filter(pfp_up<0.05 | pfp_down<0.05)
metaRP_acute_DEG_logFC <- metaRP_acute_DEG %>% filter(abs(logFC)>0.5)

metaRP_latent_DEG <- metaRP_latent_df %>% filter(pfp_up<0.05 | pfp_down<0.05)
metaRP_latent_DEG_logFC <- metaRP_latent_DEG %>% filter(abs(logFC)>0.5)

metaRP_chronic_DEG <- metaRP_chronic_df %>% filter(pfp_up<0.05 | pfp_down<0.05)
metaRP_chronic_DEG_logFC <- metaRP_chronic_DEG %>% filter(abs(logFC)>0.5)

write.csv(metaRP_acute_DEG, file = "02 Gene meta-signature/metaRP_acute_DEG.csv")
write.csv(metaRP_latent_DEG, file = "02 Gene meta-signature/metaRP_latent_DEG.csv")
write.csv(metaRP_chronic_DEG, file = "02 Gene meta-signature/metaRP_chronic_DEG.csv")

metaRP.DEG_overlap <- calculate.overlap(x = list("acute" = metaRP_acute_DEG$ID, "latent" = metaRP_latent_DEG$ID, "chronic" = metaRP_chronic_DEG$ID))

venn.plot <- draw.triple.venn(
  area1 = 2404,
  area2 = 1000,
  area3 = 373,
  n12 = 830,
  n23 = 292,
  n13 = 293,
  n123 = 258,
  lwd = rep(1, 3),
  category = c("acute", "latent", "chronic"),
  col = c("chocolate1", "green3", "purple"),
  fill = c("chocolate1", "green3", "purple"),
  alpha = rep(0.05, 3),
  fontfamily = rep("sans", 7),
  cat.fontfamily = rep("sans", 3),
  cex = 1,
  cat.cex = 1.2,
  cat.col = c("chocolate1", "green3", "purple"),
  sep.dist = 0.1)

save(metaRP_acute_DEG, metaRP_chronic_DEG, metaRP_latent_DEG, metaRP.DEG_overlap, file = "02 Gene meta-signature/metaRP.DEG.RData")

######################################################################################
RPlist2df <- function(RPlist){
  a <- RPlist
  RPdf <- cbind(a$RPs, a$RPrank, a$pfp, a$pval, a$AveFC)
  colnames(RPdf) <- c("RPs_down", "RPs_up", "RPrank_down", "RPrank_up", "pfp_down", "pfp_up", "pval_down", "pval_up", "logFC")
  RPdf <- as.data.frame(RPdf)
  RPdf$ID <- rownames(RPdf) 
  return(RPdf)
}

