library(GEOquery)
library(dplyr)
library(stringr)
library(edgeR)
library(DESeq2)
library(ggplot2)

load("06 inferring protein activity_b/human_maRt.RData")

################################ GSE127871 #################################
# Sample information
GSE127871_sampleinfo <- getGEO(filename = "07 human TLE/GSE127871_series_matrix.txt.gz", getGPL = F) %>% pData()
GSE127871_sampleinfo <- GSE127871_sampleinfo[,c(1,2,38)]
colnames(GSE127871_sampleinfo)[3] <- "frequency"

save(GSE127871_sampleinfo, file = "07 human TLE/GSE127871_sampleinfo.RData")

# use readDGE function in edgeR to read count files
GSE127871_DGE <- readDGE(files = str_c(GSE127871_sampleinfo$geo_accession, ".txt"), path = "07 human TLE/GSE127871_RAW")

# get the count matrix
GSE127871_countdata <- GSE127871_DGE$counts
GSE127871_countdata <- GSE127871_countdata[-(63677:63681), ]

save(GSE127871_countdata, file = "07 human TLE/GSE127871_countdata.RData")

# filter genes using CPM 
cpm_matrix <- cpm(GSE127871_countdata)
keep_gene <- rowSums(cpm_matrix > 0.5) >= 6

# get gene symbols
library(biomaRt)
Gene_symbols_127871 <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"), filters = "ensembl_gene_id", values = rownames(GSE127871_countdata)[keep_gene], mart = human_maRt)
Gene_symbols_127871 <- Gene_symbols_127871[!Gene_symbols_127871$hgnc_symbol == "",]

GSE127871_counts <- GSE127871_countdata[Gene_symbols_127871$ensembl_gene_id,]
rownames(GSE127871_counts) <- Gene_symbols_127871$hgnc_symbol

save(Gene_symbols_127871, GSE127871_counts, file = "07 human TLE/GSE127871_counts.RData")

############################################################################
# Create DESeq2Dataset object
dds_127871 <- DESeqDataSetFromMatrix(countData = GSE127871_counts, 
                              colData = GSE127871_sampleinfo, 
                              design = ~ frequency)

rld_127871 <- rlog(dds_127871)

save(dds_127871, rld_127871, file = "07 human TLE/dds_rld_127871.RData")

# Spearman rank correlation 
gene_frequency_cor <- cor(t(assay(rld_127871)), as.numeric(GSE127871_sampleinfo$frequency),method = "spearman") %>% as.data.frame()
gene_frequency_cor$Gene_symbol <- rownames(gene_frequency_cor)

pvalue <- matrix()
for (i in 1:nrow(gene_frequency_cor)){
  pvalue[i] <- cor.test(assay(rld_127871)[i,], as.numeric(GSE127871_sampleinfo$frequency), method = "spearman")$p.value
}

gene_frequency_cor$p.value <- pvalue

gene_frequency_cor$p.adjust <- p.adjust(pvalue, method = "BH", n = length(pvalue))

save(gene_frequency_cor, file = "07 human TLE/gene_frequency_cor.RData")

gene_frequency_cor_0.05 <- gene_frequency_cor %>% filter(p.value<0.05)
write.csv(gene_frequency_cor_0.05, file = "07 human TLE/gene_frequency_cor_0.05.csv")

DER_frequency <- gene_frequency_cor_0.05[match(DEregulators_m2h$HGNC.symbol, gene_frequency_cor_0.05$Gene_symbol),] %>% na.omit()

DER_frequency_activity <- merge(DER_frequency, DEregulators_m2h, by.x="Gene_symbol", by.y = "HGNC.symbol")

SF_acute <- ggplot(data = DER_frequency_activity, aes(x=NES_acute,y=V1))+
  geom_point(shape=21,color="black",alpha=0.5,fill="gray")+
  theme_bw()+
  geom_vline(xintercept = 0, color="black",linetype="dashed") +
  geom_hline(yintercept = 0, color="black",linetype="dashed")

SF_latent <- ggplot(data = DER_frequency_activity, aes(x=NES_latent,y=V1))+
  geom_point(shape=21,color="black",alpha=0.5,fill="gray")+
  theme_bw()+
  geom_vline(xintercept = 0, color="black",linetype="dashed") +
  geom_hline(yintercept = 0, color="black",linetype="dashed")

SF_chronic <- ggplot(data = DER_frequency_activity, aes(x=NES_chronic,y=V1))+
  geom_point(shape=21,color="black",alpha=0.5,fill="gray")+
  theme_bw()+
  geom_vline(xintercept = 0, color="black",linetype="dashed") +
  geom_hline(yintercept = 0, color="black",linetype="dashed")



for (i in 1:10) {
  geneData <- data.frame(log2counts = assay(rld_127871)[DER_both$Gene_symbol[i],],
                         frequency = as.numeric(GSE127871_sampleinfo$frequency))
  temp_plot <- ggplot(geneData, aes(x=frequency, y=log2counts)) + 
    geom_point(alpha=0.5) + 
    geom_smooth(method=lm, se=FALSE) +
    theme_bw() + 
    scale_x_sqrt(breaks=c(0,4,30,60,90,120)) +
    labs(title=DER_both$Gene_symbol[i],x="Seizure frequency", y = "log2(counts)")
  ggsave(temp_plot, file=paste0("plot_", i,".png"), width = 2, height = 2, units = "in")
}

fre_plot = list()
for (i in 1:10) {
  geneData <- data.frame(log2counts = assay(rld_127871)[DER_both$Gene_symbol[i],],
                         frequency = as.numeric(GSE127871_sampleinfo$frequency))
  fre_plot[[i]] <- ggplot(geneData, aes(x=frequency, y=log2counts)) + 
    geom_point(alpha=0.5) + 
    geom_smooth(method=lm, se=FALSE) +
    theme_bw() + 
    scale_x_sqrt(breaks=c(0,4,30,60,90,120)) +
    labs(title=DER_both$Gene_symbol[i],x="Seizure frequency", y = "log2(counts)")
}

library(ggpubr)
ggarrange(fre_plot[[1]], fre_plot[[2]],fre_plot[[3]],fre_plot[[4]],fre_plot[[5]],fre_plot[[6]],fre_plot[[7]],fre_plot[[8]],fre_plot[[9]],fre_plot[[10]],ncol = 5,nrow = 2,align = "hv")

# 7.79,3.89

