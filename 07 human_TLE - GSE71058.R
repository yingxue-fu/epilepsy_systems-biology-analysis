library(GEOquery)
library(dplyr)
library(stringr)
library(DESeq2)
library(edgeR)
library(RColorBrewer)

# GSE71085
############################# Data preprocessing ############################
# Sample information
GSE71058_sampleinfo <- getGEO(filename = "07 human TLE/GSE71058_series_matrix.txt.gz", getGPL = F) %>% pData()
GSE71058_sampleinfo <- GSE71058_sampleinfo[,c(1,2,40,41)]
GSE71058_sampleinfo$diagnosis <- rep(c("hs", "no_hs", "hs", "no_hs", "hs", "no_hs", "hs", "no_hs"), c(4,2,2,2,1,2,1,8)) %>% as.factor()
GSE71058_sampleinfo$title <- str_remove(GSE71058_sampleinfo$title, "duke")
rownames(GSE71058_sampleinfo) <- GSE71058_sampleinfo$title

save(GSE71058_sampleinfo, file = "07 human TLE/GSE71058_sampleinfo.RData")

# read count matrix
GSE71058_seqdata <- read.delim("07 human TLE/GSE71058_RNA_rawcounts_12samples_22alignments.txt", stringsAsFactors = FALSE)
GSE71058_countdata <- GSE71058_seqdata[,-1]

# Gene symbols and Sample names
write.csv(GSE71058_seqdata[,1], file = "07 human TLE/GSE71058_Gene_symbol.csv")
GSE71058_genes <- read.csv("07 human TLE/GSE71058_Gene_symbol - noDate.csv")
rownames(GSE71058_countdata) <- GSE71058_genes$x
colnames(GSE71058_countdata) <- str_remove(colnames(GSE71058_countdata), "duke")

save(GSE71058_countdata, file = "07 human TLE/GSE71058_countdata.RData")

############################ Data exploration ############################

# filter genes using CPM 
cpm_matrix <- cpm(GSE71058_countdata)
# or calculate it manually
# cpm_matrix_m <- apply(GSE71058_countdata, 2, function(x)(x/sum(x))*1000000)

mean(apply(GSE71058_countdata, 2, sum))

keep_gene <- rowSums(cpm_matrix > 0.1) >= 8

GSE71058_counts <- GSE71058_countdata[keep_gene,]

save(GSE71058_counts, file = "07 human TLE/GSE71058_counts.RData")

library(RUVSeq)
set <- newSeqExpressionSet(as.matrix(GSE71058_counts), phenoData = GSE71058_sampleinfo)
x <- factor(set$diagnosis, levels = c("no_hs", "hs"))

library(RColorBrewer)
colors <- brewer.pal(3, "Set2")
par(mar=c(7,3,1,1),cex=0.7)
plotRLE(set, outline=FALSE, ylim=c(-8, 6), col=colors[x],las=2)
par(mar=c(7,5,1,1),cex=0.7)
plotPCA(set, col=colors[x])

set <- betweenLaneNormalization(set, which="upper")
par(mar=c(7,3,1,1),cex=0.7)
plotRLE(set, outline=FALSE, ylim=c(-8, 6), col=colors[x],las=2)
par(mar=c(7,5,1,1),cex=0.7)
plotPCA(set, col=colors[x])

design <- model.matrix(~x, data=pData(set))
y <- DGEList(counts=counts(set), group=x)
y <- calcNormFactors(y, method="upperquartile")
y <- estimateGLMCommonDisp(y, design)
y <- estimateGLMTagwiseDisp(y, design)
fit <- glmFit(y, design)
lrt <- glmLRT(fit, coef=2)
top <- topTags(lrt, n=nrow(set))$table
empirical <- rownames(set)[which(!(rownames(set) %in% rownames(top)[1:5000]))]

set2 <- RUVg(set, empirical, k=1)
pData(set2)

par(mar=c(7,3,1,1),cex=0.7)
plotRLE(set2, outline=FALSE, ylim=c(-4, 4), col=colors[x],las=2)
par(mar=c(7,5,1,1),cex=0.7)
plotPCA(set2, col=colors[x])

save(set2, file = "07 human TLE/RUVg_set2.RData")

# DE analysis
library(DESeq2)
dds_71058 <- DESeqDataSetFromMatrix(countData = counts(set2),
                              colData = pData(set2),
                              design = ~ W_1 + diagnosis)
dds_71058$diagnosis <- factor(dds_71058$diagnosis, levels = c("no_hs", "hs"))

rld_71058 <- rlog(dds_71058, blind = FALSE)

save(dds_71058, rld_71058, file = "07 human TLE/dds_71058.RData")

dds_71058 <- DESeq(dds_71058)
res <- results(dds_71058)
res_Sort <- res[order(res$pvalue),] %>% as.data.frame()
res_Sort$Gene_symbol <- rownames(res_Sort)

hs_DEG <- res_Sort %>% filter(padj < 0.05 & abs(log2FoldChange) >0.5)
write.csv(hs_DEG, file = "07 human TLE/hs_DEG.csv")

save(hs_DEG, file = "07 human TLE/hs_DEG.RData")

DER_hs_nohs <- intersect(hs_DEG$Gene_symbol, DEregulators_m2h$HGNC.symbol)

DER_hs_nohs <- hs_DEG[match(DER_hs_nohs, hs_DEG$Gene_symbol),]

DER_hs_nohs_activity <- merge(DER_hs_nohs, DEregulators_m2h, by.x="Gene_symbol", by.y = "HGNC.symbol")

HS_acute <- ggplot(data = DER_hs_nohs_activity, aes(x=NES_acute,y=log2FoldChange))+
  geom_point(shape=21,color="black",alpha=0.5,fill="gray")+
  theme_bw()+
  geom_vline(xintercept = 0, color="black",linetype="dashed") +
  geom_hline(yintercept = 0, color="black",linetype="dashed")

HS_latent <- ggplot(data = DER_hs_nohs_activity, aes(x=NES_latent,y=log2FoldChange))+
  geom_point(shape=21,color="black",alpha=0.5,fill="gray")+
  theme_bw()+
  geom_vline(xintercept = 0, color="black",linetype="dashed") +
  geom_hline(yintercept = 0, color="black",linetype="dashed")

HS_chronic <- ggplot(data = DER_hs_nohs_activity, aes(x=NES_chronic,y=log2FoldChange))+
  geom_point(shape=21,color="black",alpha=0.5,fill="gray")+
  theme_bw()+
  geom_vline(xintercept = 0, color="black",linetype="dashed") +
  geom_hline(yintercept = 0, color="black",linetype="dashed")

library(ggpubr)
ggarrange(SF_acute,SF_latent,SF_chronic,HS_acute,HS_latent,HS_chronic,nrow = 2,ncol = 3,align = "hv")

DER_both <- merge(DER_frequency, DER_hs_nohs, by ="Gene_symbol")

save(DER_hs_nohs, DER_both, file = "07 human TLE/DER_hs_nohs.RData")



for (i in 1:10) {
  geneData <- data.frame(log2counts = assay(rld_71058)[DER_both$Gene_symbol[i],],
                         condition = colData(rld_71058)$diagnosis)
  temp_plot <- ggplot(geneData, aes(x=condition, y=log2counts, fill=condition)) + 
    geom_boxplot(alpha=0.7) + 
    geom_dotplot(binaxis='y', stackdir='center', dotsize=1) +
    theme_bw() + 
    theme(legend.position="none") + 
    scale_fill_brewer(palette = "Paired") +
    labs(title=DER_both$Gene_symbol[i],x="condition", y = "log2(counts)")
  ggsave(temp_plot, file=paste0("plot_", i+11,".png"), width = 1.5, height = 2, units = "in")
}

hs_plot = list()
for (i in 1:10) {
  geneData <- data.frame(log2counts = assay(rld_71058)[DER_both$Gene_symbol[i],],
                         condition = colData(rld_71058)$diagnosis)
  hs_plot[[i]] <- ggplot(geneData, aes(x=condition, y=log2counts, fill=condition)) + 
    geom_boxplot(alpha=0.7) + 
    geom_dotplot(binaxis='y', stackdir='center', dotsize=1) +
    theme_bw() + 
    theme(legend.position="none") + 
    scale_fill_brewer(palette = "Paired") +
    labs(title=DER_both$Gene_symbol[i],x="condition", y = "log2(counts)")
}

library(ggpubr)
ggarrange(hs_plot[[1]], hs_plot[[2]],hs_plot[[3]],hs_plot[[4]],hs_plot[[5]],hs_plot[[6]],hs_plot[[7]],hs_plot[[8]],hs_plot[[9]],hs_plot[[10]],ncol = 5,nrow = 2,align = "hv")

