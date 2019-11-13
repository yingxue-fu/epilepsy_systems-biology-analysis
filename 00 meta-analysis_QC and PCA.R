library(GEOquery)
library(arrayQualityMetrics)
library(ggfortify)

# control gray # acute #c00000 # latent #0070c0 # chronic #00b050

# GSE14763
gse14763 <- getGEO(filename = "00 meta-analysis/GSE14763_series_matrix.txt.gz", getGPL = F)
sample_info_gse14763 <- pData(gse14763)
# QC
arrayQualityMetrics(expressionset = gse14763,
                    outdir = "00 meta-analysis/Report_for_gse14763",
                    intgroup = "source_name_ch1",
                    force = TRUE)
# PCA
exp_gse14763 <- na.omit(log2(exprs(gse14763)+1))

pca_gse14763 <- prcomp(t(exp_gse14763))
# Define stage
sample_info_gse14763$Stage <- rep(c("acute", "latent", "control"), c(5,10,3))
# draw plot
autoplot(pca_gse14763, data = sample_info_gse14763, colour = 'Stage', alpha = 0.7) +
  scale_x_continuous(limits = c(-0.5, 0.6)) +
  scale_y_continuous(limits = c(-0.6, 0.5)) + 
  theme_bw() +
  scale_color_manual(values = c("#c00000", "gray30", "#0070c0")) +
  geom_hline(yintercept = 0, linetype="dashed", alpha=0.5) +
  geom_vline(xintercept = 0, linetype="dashed", alpha=0.5)

########################################################################################
# GSE27166
gse27166 <- getGEO(filename = "00 meta-analysis/GSE27166_series_matrix.txt.gz", getGPL = F)
dim(gse27166)
sample_info_gse27166 <- pData(gse27166)
# QC
arrayQualityMetrics(expressionset = gse27166,
                    outdir = "00 meta-analysis/Report_for_gse27166_all",
                    intgroup = c("characteristics_ch1.1","characteristics_ch1.2"),
                    force = TRUE)
arrayQualityMetrics(expressionset = 
                      gse27166[, gse27166$characteristics_ch1.2 == "side: right"],
                    outdir = "00 meta-analysis/Report_for_gse27166",
                    intgroup = c("characteristics_ch1.1","characteristics_ch1.2"),
                    force = TRUE)
arrayQualityMetrics(expressionset = 
                      gse27166[, gse27166$characteristics_ch1.2 == "side: left"],
                    outdir = "00 meta-analysis/Report_for_gse27166_left",
                    intgroup = c("characteristics_ch1.1","characteristics_ch1.2"),
                    force = TRUE)
# PCA
exp_gse27166 <- na.omit(log2(exprs(gse27166)+1))
exp_gse27166_right <- exp_gse27166[, gse27166$characteristics_ch1.2 == "side: right"]

pca_gse27166 <- prcomp(t(exp_gse27166_right[,-c(7,12)]))

sample_info_gse27166_right <- sample_info_gse27166[gse27166$characteristics_ch1.2 == "side: right",]
# define stage
sample_info_gse27166_right$Stage <- rep(c("chronic", "control", "chronic", "control"), c(3,3,3,3))
# draw plot
autoplot(pca_gse27166, data = sample_info_gse27166_right[-c(7,12),], colour = 'Stage', alpha = 0.7) +
  scale_x_continuous(limits = c(-0.6, 0.6)) +
  scale_y_continuous(limits = c(-0.5, 0.7)) + theme_bw() +
  scale_color_manual(values = c("#00b050", "gray30")) +
  geom_hline(yintercept = 0, linetype="dashed", alpha=0.5) +
  geom_vline(xintercept = 0, linetype="dashed", alpha=0.5)

######################################################################################
# GSE27015
gse27015 <- getGEO(filename = "00 meta-analysis/GSE27015_series_matrix.txt.gz", getGPL = F)
dim(gse27015)
sample_info_gse27015 <- pData(gse27015)
# QC
arrayQualityMetrics(expressionset = gse27015,
                    outdir = "00 meta-analysis/Report_for_gse27015_all",
                    intgroup = c("characteristics_ch1.2","characteristics_ch1.1"),
                    force = TRUE)
arrayQualityMetrics(expressionset = 
                      gse27015[, gse27015$characteristics_ch1.2 == "side: right"],
                    outdir = "00 meta-analysis/Report_for_gse27015",
                    intgroup = c("characteristics_ch1.1","characteristics_ch1.2"),
                    force = TRUE)

######################################################################################
# GSE49849
gse49849 <- getGEO(filename = "00 meta-analysis/GSE49849_series_matrix.txt.gz", getGPL = F)
dim(gse49849)
sample_info_gse49849 <- pData(gse49849)
# QC
arrayQualityMetrics(expressionset = gse49849,
                    outdir = "00 meta-analysis/Report_for_gse49849",
                    intgroup = c("characteristics_ch1.3","characteristics_ch1.2"),
                    force = TRUE)

# PCA
exp_gse49849 <- na.omit(exprs(gse49849))
pca_gse49849 <- prcomp(t(exp_gse49849))

sample_info_gse49849$Stage <- as.factor(rep(c("control", "acute", "control", "latent"), c(5,5,5,5)))

autoplot(pca_gse49849, data = sample_info_gse49849, colour = 'Stage', alpha = 0.7) +
  scale_x_continuous(limits = c(-0.6, 0.4)) +
  scale_y_continuous(limits = c(-0.5, 0.5)) + theme_bw() +
  scale_color_manual(values = c("#c00000", "gray30", "#0070c0")) +
  geom_hline(yintercept = 0, linetype="dashed", alpha=0.5) +
  geom_vline(xintercept = 0, linetype="dashed", alpha=0.5)

#######################################################################################
# GSE73878
gse73878 <- getGEO(filename = "00 meta-analysis/GSE73878_series_matrix.txt.gz", getGPL = F)
dim(gse73878)
sample_info_gse73878 <- pData(gse73878)
# QC
arrayQualityMetrics(expressionset = gse73878[ ,gse73878$`hemisphere:ch1` == "IL"],
                    outdir = "00 meta-analysis/Report_for_gse73878",
                    intgroup = c("treatment:ch1","treatment time:ch1"),
                    force = TRUE)

# PCA
exp_gse73878 <- na.omit(exprs(gse73878))
exp_gse73878_IL <- exp_gse73878[,gse73878$`hemisphere:ch1` == "IL"]

pca_gse73878 <- prcomp(t(exp_gse73878_IL[,-c(16,11,29,30,31,32)]))

sample_info_gse73878_IL <- sample_info_gse73878[gse73878$`hemisphere:ch1` == "IL",]
sample_info_gse73878_IL <- sample_info_gse73878_IL[-c(16,11,29,30,31,32),]
sample_info_gse73878_IL$group <- str_c(sample_info_gse73878_IL$`treatment:ch1`, sample_info_gse73878_IL$`treatment time:ch1`,sep = "_")
# define stages
write.csv(sample_info_gse73878_IL, file = "00 meta-analysis/sample_info_gse73878_IL.csv")
sample_info_gse73878_IL <- read.csv(file = "00 meta-analysis/sample_info_gse73878_IL.csv")

autoplot(pca_gse73878, data = sample_info_gse73878_IL, colour = 'Stage', alpha = 0.7) +
  scale_x_continuous(limits = c(-0.4, 0.3)) +
  scale_y_continuous(limits = c(-0.4, 0.4)) + theme_bw() +
  scale_color_manual(values = c("#c00000", "#00b050", "gray30", "#0070c0")) +
  geom_hline(yintercept = 0, linetype="dashed", alpha=0.5) +
  geom_vline(xintercept = 0, linetype="dashed", alpha=0.5)

# GSE88992
gse88992 <- getGEO(filename = "00 meta-analysis/GSE88992_series_matrix.txt.gz", getGPL = F)
dim(gse88992)
sample_info_gse88992 <- pData(gse88992)
library(stringr)
library(dplyr)
condition <- str_split(sample_info_gse88992$source_name_ch1, ", ") %>% as.data.frame() %>% t()
gse88992$treatment <- condition[, 2]
gse88992$treatment_time <- condition[, 3]

arrayQualityMetrics(expressionset = gse88992,
                    outdir = "00 meta-analysis/Report_for_gse88992",
                    intgroup = c("treatment","treatment_time"),
                    force = TRUE)

# PCA
exp_gse88992 <- na.omit(exprs(gse88992))

pca_gse88992 <- prcomp(t(exp_gse88992))

sample_info_gse88992$Stage <- rep(c("control","acute","control","acute","control","acute"),c(3,3,3,3,3,2))

autoplot(pca_gse88992, data = sample_info_gse88992, colour = 'Stage', alpha = 0.7) +
  scale_x_continuous(limits = c(-0.3, 0.4)) +
  scale_y_continuous(limits = c(-0.7, 0.4)) + theme_bw() +
  scale_color_manual(values = c("#c00000", "gray30")) +
  geom_hline(yintercept = 0, linetype="dashed", alpha=0.5) +
  geom_vline(xintercept = 0, linetype="dashed", alpha=0.5)




save(gse14763, gse27166, gse49849, gse73878, gse88992, file = "00 meta-analysis/meta_gse.RData")
save(sample_info_gse14763, sample_info_gse27166_right, sample_info_gse49849, sample_info_gse73878_IL, sample_info_gse88992, file = "00 meta-analysis/meta_sample_info.RData")
