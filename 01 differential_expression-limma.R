library(Biobase)
library(GEOquery)
library(limma)

load("00 meta-analysis/meta_gse.RData")
load("00 meta-analysis/meta_sample_info.RData")

####################################### GSE14763 ########################################

# log2 transform
quantile(exprs(gse14763), c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T)
exprs(gse14763) <- log2(exprs(gse14763)+1)

fl <- as.factor(sample_info_gse14763$Stage)
gse14763$description <- fl

design_gse14763 <- model.matrix(~ description + 0, gse14763)
colnames(design_gse14763) <- levels(fl)

fit <- lmFit(gse14763, design_gse14763)

cont.matrix <- makeContrasts(acute-control, latent-control, levels=design_gse14763)

fit2 <- contrasts.fit(fit, cont.matrix)

fit2_gse14763 <- eBayes(fit2, 0.01)

# DEGs
DEG_acute_gse14763 <- topTable(fit2_gse14763, coef = 1, number = 10000, p.value = 0.05, lfc = 0.5)
DEG_latent_gse14763 <- topTable(fit2_gse14763, coef = 2, number = 10000, p.value = 0.05, lfc = 0.5)

# match to probe annotation
gpl2896 <- getGEO(filename = "00 meta-analysis/GPL2896.annot.gz")
gpl2896 <- Table(gpl2896)
gpl2896 <- gpl2896[!gpl2896$`Gene ID` == "", 1:4]

DEG_acute_gse14763$ID = rownames(DEG_acute_gse14763); 
DEG_latent_gse14763$ID = rownames(DEG_latent_gse14763)

DEG_acute_gse14763 = merge(DEG_acute_gse14763, gpl2896, by = "ID")
DEG_latent_gse14763 = merge(DEG_latent_gse14763, gpl2896, by = "ID")



####################################### GSE27166 ########################################

# log2 transform
quantile(exprs(gse27166), c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T)
exprs(gse27166) <- log2(exprs(gse27166)+1)

fl <- as.factor(sample_info_gse27166_right$Stage)
gse27166_right <- gse27166[, gse27166$characteristics_ch1.2 == "side: right"]

# remove outlier 
fl <- fl[-c(7,12)]
gse27166_right <- gse27166_right[, -c(7,12)]

gse27166_right$description <- fl

design_gse27166 <- model.matrix(~ description + 0, gse27166_right)
colnames(design_gse27166) <- levels(fl)

fit <- lmFit(gse27166_right, design_gse27166)

cont.matrix <- makeContrasts(chronic-control, levels=design_gse27166)

fit2 <- contrasts.fit(fit, cont.matrix)

fit2_gse27166 <- eBayes(fit2, 0.01)

# DEGs
DEG_chronic_gse27166 <- topTable(fit2_gse27166, number = 1900, p.value = 0.05)

# no significant DE genes


####################################### GSE49849 #######################################

# log2 transform
quantile(exprs(gse49849), c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T)

fl <- as.factor(sample_info_gse49849$Stage)
gse49849$description <- fl

design_gse49849 <- model.matrix(~ description + 0, gse49849)
colnames(design_gse49849) <- levels(fl)

fit <- lmFit(gse49849, design_gse49849)

cont.matrix <- makeContrasts(acute-control, latent-control, levels=design_gse49849)

fit2 <- contrasts.fit(fit, cont.matrix)

fit2_gse49849 <- eBayes(fit2, 0.01)

# DEGs
DEG_acute_gse49849 <- topTable(fit2_gse49849, coef = 1, number = 10000, p.value = 0.05, lfc = 0.5)
DEG_latent_gse49849 <- topTable(fit2_gse49849, coef = 2, number = 10000, p.value = 0.05, lfc = 0.5)


# match to probe annotation
gpl6247 <- getGEO(filename = "00 meta-analysis/GPL6247.annot.gz")
gpl6247 <- Table(gpl6247)
gpl6247 <- gpl6247[!gpl6247$`Gene ID` == "", 1:4]

DEG_acute_gse49849$ID = rownames(DEG_acute_gse49849)
DEG_latent_gse49849$ID = rownames(DEG_latent_gse49849)

DEG_acute_gse49849 = merge(DEG_acute_gse49849, gpl6247, by = "ID")
DEG_latent_gse49849 = merge(DEG_latent_gse49849, gpl6247, by = "ID")

####################################### GSE73878 ########################################

# log2 transform
quantile(exprs(gse73878), c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T)

fl <- as.factor(sample_info_gse73878_IL$Stage)
gse73878_IL <- gse73878[ ,gse73878$`hemisphere:ch1` == "IL"]
gse73878_IL <- gse73878_IL[,-c(16,11,29,30,31,32)]
gse73878_IL$description <- fl

design_gse73878 <- model.matrix(~ description + 0, gse73878_IL)
colnames(design_gse73878) <- levels(fl)

fit <- lmFit(gse73878_IL, design_gse73878)

cont.matrix <- makeContrasts(acute-control, latent-control, chronic-control, levels=design_gse73878)

fit2 <- contrasts.fit(fit, cont.matrix)

fit2_gse73878 <- eBayes(fit2, 0.01)

# DEGs
DEG_acute_gse73878 <- topTable(fit2_gse73878, coef = 1, number = 10000, p.value = 0.05, lfc = 0.5)
DEG_latent_gse73878 <- topTable(fit2_gse73878, coef = 2, number = 10000, p.value = 0.05, lfc = 0.5)
DEG_chronic_gse73878 <- topTable(fit2_gse73878, coef = 3, number = 10000, p.value = 0.05, lfc = 0.5)

# match to probe annotation
gpl6885 <- getGEO(filename = "00 meta-analysis/GPL6885.annot.gz")
gpl6885 <- Table(gpl6885)
gpl6885 <- na.omit(gpl6885[, 1:4])

DEG_acute_gse73878$ID = rownames(DEG_acute_gse73878)
DEG_latent_gse73878$ID = rownames(DEG_latent_gse73878)
DEG_chronic_gse73878$ID = rownames(DEG_chronic_gse73878)

DEG_acute_gse73878 = merge(DEG_acute_gse73878, gpl6885, by = "ID")
DEG_latent_gse73878 = merge(DEG_latent_gse73878, gpl6885, by = "ID")
DEG_chronic_gse73878 = merge(DEG_chronic_gse73878, gpl6885, by = "ID")

####################################### GSE88992 ########################################

# log2 transform
quantile(exprs(gse88992), c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T)

fl <- as.factor(sample_info_gse88992$Stage)
gse88992$description <- fl

design_gse88992 <- model.matrix(~ description + 0, gse88992)
colnames(design_gse88992) <- levels(fl)

fit <- lmFit(gse88992, design_gse88992)

cont.matrix <- makeContrasts(acute-control, levels=design_gse88992)

fit2 <- contrasts.fit(fit, cont.matrix)

fit2_gse88992 <- eBayes(fit2, 0.01)

# DEGs
DEG_acute_gse88992 <- topTable(fit2_gse88992, number = 20000, p.value = 0.05, lfc = 0.5)

# match to probe annotation
gpl1261 <- getGEO(filename = "00 meta-analysis/GPL1261.annot.gz")
gpl1261 <- Table(gpl1261)
gpl1261 <- gpl1261[!gpl1261$`Gene ID` == "", 1:4]

DEG_acute_gse88992$ID <- rownames(DEG_acute_gse88992)

DEG_acute_gse88992 <- merge(DEG_acute_gse88992, gpl1261, by = "ID")

#############################################################################################
# DEG acute
write.csv(unique(DEG_acute_gse14763$`Gene symbol`), file = "01 differential expression-limma/DEG_acute_gse14763.csv")
write.csv(unique(DEG_acute_gse49849$`Gene symbol`), file = "01 differential expression-limma/DEG_acute_gse49849.csv")
write.csv(unique(DEG_acute_gse73878$`Gene symbol`), file = "01 differential expression-limma/DEG_acute_gse73878.csv")
write.csv(unique(DEG_acute_gse88992_fc$`Gene symbol`), file = "01 differential expression-limma/DEG_acute_gse88992.csv")

# DEG latent
write.csv(unique(DEG_latent_gse14763$`Gene symbol`), file = "01 differential expression-limma/DEG_latent_gse14763.csv")
write.csv(unique(DEG_latent_gse49849$`Gene symbol`), file = "01 differential expression-limma/DEG_latent_gse49849.csv")
write.csv(unique(DEG_latent_gse73878$`Gene symbol`), file = "01 differential expression-limma/DEG_latent_gse73878.csv")

# DEG chronic
write.csv(unique(DEG_chronic_gse27166$`Gene symbol`), file = "01 differential expression-limma/DEG_chronic_gse27166.csv")
write.csv(unique(DEG_chronic_gse73878$`Gene symbol`), file = "01 differential expression-limma/DEG_chronic_gse73878.csv")

############################################################################################

save(DEG_acute_gse14763, DEG_acute_gse49849, DEG_acute_gse73878, DEG_acute_gse88992, file = "01 differential expression-limma/DEG_acute.RData")
save(DEG_latent_gse14763, DEG_latent_gse49849, DEG_latent_gse73878, file = "01 differential expression-limma/DEG_latent.RData")
save(DEG_chronic_gse27166, DEG_chronic_gse73878, file = "01 differential expression-limma/DEG_chronic.RData")
