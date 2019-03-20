library(Biobase)
library(GEOquery)
library(WGCNA)
library(dplyr)

# load series and platform data from GEO
gset <- getGEO(filename = "01 sample and gene filtering/GSE73878_series_matrix.txt.gz", GSEMatrix = TRUE, getGPL = F)
# keep only the "IL" samples
gset_IL = gset[, gset$`hemisphere:ch1` == "IL"]

# extract phenotype information of the samples
ph_IL <- pData(gset_IL)
# have a look on sample numbers of different groups
table(ph_IL[, c("treatment:ch1", "treatment time:ch1")])

# extract gene expresssion matrix
ex_IL = exprs(gset_IL)
dim(ex_IL)
# get a feeling about the expression data
par(las = 2, mar = c(8, 4, 2, 2), cex = 0.6)
boxplot(ex_IL, pch = 20, cex = 0.8, col = "lightgray") # the data have already been log2 transformed and quantile normalized

# 1. cluster the samples to detect outliers
sampleTree = hclust(dist(t(ex_IL)), method = "average")
# Plot the sample tree
par(cex = 0.6, mar = c(1,5,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5, cex.axis = 1.5, cex.main = 1.5) # export the figure
# there seems to be 4 outliers
abline(h = 70, col = "red")
clust = cutreeStatic(sampleTree, cutHeight = 70, minSize = 5)
table(clust)
# clust 0 contains the samples we want to remove.
rmSamples = (clust == 0)
datExpr0 = ex_IL[, !rmSamples]
dim(datExpr0)

# 2. filter low-expression probes
# histogram of expression matrix
ex_hist = hist(datExpr0, breaks = 100, freq = F)
# determine the threshold
th = ex_hist$breaks[ex_hist$counts == max(ex_hist$counts)]
abline(v = th, col = "red")
# filter probes that are lowly expressed in more than 75% of the samples
low_e_probes = rowQuantileC(datExpr0, 0.7) < th
table(low_e_probes)

datExpr0 = datExpr0[!low_e_probes, ]
dim(datExpr0)

# 3. filter low-variation probes based on coefficient of variation (CV)
CV = function(x){
  cv = sd(x)/mean(x)
}
probes_CV = apply(datExpr0, 1, CV)
# histogram of probes CV
hist(probes_CV, breaks = 100, freq = F)
abline(v = quantile(probes_CV, 0.2), col = "red")
# remove 20% of the probes with lowest CV
low_cv_probes = probes_CV < quantile(probes_CV, 0.2)
table(low_cv_probes)

datExpr0 = datExpr0[!low_cv_probes, ]
dim(datExpr0)

# 4. finally, filter the redundant probes of the same gene based on CV
probeCV_r = apply(datExpr0, 1, CV)
head(probeCV_r)

## extract probe annotation from GPL data
gpl = getGEO(filename = "01 sample and gene filtering/GPL6885.annot.gz")
probe_annot = na.omit(Table(gpl)[, 1:4])
# save the probe annotation data for further analysis
save(probe_annot, file = "01 sample and gene filtering/probes_annotation.RData")

# match probe CV to geneID
probeCV_r = as.data.frame(probeCV_r, row.names = names(probeCV_r))
probeCV_r$ID = rownames(probeCV_r)
probe_CV_geneID = merge(probe_annot, probeCV_r, by = "ID")
colnames(probe_CV_geneID) = c("ID", "Gene_title", "Gene_symbol", "Gene_ID", "probeCV")

probes_final <- probe_CV_geneID %>% group_by(Gene_ID) %>% top_n(1, probeCV)

datExpr = datExpr0[probes_final$ID, ]
dim(datExpr)

save(datExpr, probes_final, file = "01 sample and gene filtering/samples_and_probes_filtering.RData")
