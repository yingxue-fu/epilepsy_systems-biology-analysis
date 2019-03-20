library(WGCNA)
library(tidyverse)
library(stringr)
library(ggrepel)

load("01 sample and gene filtering/samples_and_probes_filtering.RData")
datExpr = t(datExpr)

# 1. Calculate biweight midcorrelation between genes
biweight_coref = bicor(datExpr)

# 2. determine the soft-thresholding powers
# convert negative correlation to 0 for "signed hybrid" network type
biweight_coref_signed = biweight_coref 
biweight_coref_signed[biweight_coref_signed < 0] = 0

# Network topology analysis to choose the soft power
powers = c(seq(1, 10, by = 1), seq(12, 20, by = 2))
sft = pickSoftThreshold.fromSimilarity(biweight_coref_signed, blockSize = 5500, verbose = 5)

par(mfrow = c(1,2)); cex1 = 0.9

# Scale-free fit indices as a function of the soft-thresholding power 
plot(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3])*sft$fitIndices[, 2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2", 
     type="n", main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
abline(h = 0.98, col = "red")

# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[, 1], sft$fitIndices[, 5],
     xlab = "Soft Threshold (power)",ylab = "Mean Connectivity", type = "n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[, 1], sft$fitIndices[, 5], labels = powers, cex = cex1,col = "red")
abline(h = 30, col = "red")

# Soft power is chosed to be 6
softPower = 6

# Claculate ajacency matrix
adjacency = adjacency.fromSimilarity(biweight_coref, 
                         type = "signed hybrid", 
                         power = 6)
save(adjacency, file = "02 gene coexpression network analysis/adjacency.RData")

# Calculate TOM matrix
TOM = TOMsimilarity(adjacency)
save(TOM, file = "02 gene coexpression network analysis/TOM.RData")

# Cluster genes based on TOM similarity
dissTOM = 1-TOM
geneTree = hclust(as.dist(dissTOM), method = "average")
plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity", labels = FALSE, hang = 0.04)

# Set the minimum module size
minModuleSize = 30

# Module identification using dynamic tree cut:
original_Modules = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                            deepSplit = 2, pamRespectsDendro = FALSE,
                            minClusterSize = minModuleSize)
head(original_Modules)
table(original_Modules)

# Convert numeric lables into colors
original_Colors = labels2colors(original_Modules, colorSeq = standardColors()[12:35])
head(original_Colors)
table(original_Colors)

# Plot the dendrogram and colors underneath
plotDendroAndColors(geneTree, original_Colors, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors")

# Calculate module eigengenes
MEList = moduleEigengenes(datExpr, colors = original_Colors)
MEs = MEList$eigengenes

colorOrder = data.frame(color = c("grey", standardColors()[12:35]), order1 = 0:24)
MEcolor2module = colorOrder$order1[match(str_remove(names(MEs), "ME"), colorOrder$color)]
names(MEs) = paste("M", MEcolor2module, sep = "")

# Module cluster based on their eigengenes
METree = hclust(as.dist(1-cor(MEs)), method = "average")
plot(METree, main = "Clustering of original module eigengenes", xlab = "", sub = "")
abline(h = 0.2, col = "red")

# Merge similar modules
merge = mergeCloseModules(datExpr, original_Modules, cutHeight = 0.2, colorSeq = c("grey", standardColors()[12:35]), verbose = 3)
mergedColors = original_Colors[match(merge$colors, original_Modules)]

# Comparsion before and after the merge
plotDendroAndColors(geneTree, cbind(original_Colors, mergedColors),
                    c("Dynamic Tree Cut", "Merged dynamic"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)

# New module colors and labels
moduleColors = mergedColors
moduleLabels = match(moduleColors, colorOrder$color)-1
ordered_moduleLabels = as.data.frame(table(moduleLabels))
ordered_moduleLabels$orderedLabels = 0:15

new_moduleLabels = paste("M", ordered_moduleLabels$orderedLabels[match(moduleLabels, ordered_moduleLabels$moduleLabels)], sep = "")
table(new_moduleLabels)

save(moduleColors, moduleLabels, new_moduleLabels, file = "02 gene coexpression network analysis/module_color_labels.RData")

# gene cluster and final modules
plotDendroAndColors(geneTree, moduleColors, groupLabels = "Modules", rowText = new_moduleLabels, rowTextAlignment = "center", addTextGuide = T, dendroLabels = FALSE, hang = 0.03, addGuide = F, guideHang = 0.05)

# save probe_modules data
probe_modules = data.frame(probe_ID = colnames(datExpr), Modules = new_moduleLabels, Colors = moduleColors)
probe_modules_anno = merge(probe_modules, probes_final, by.x = "probe_ID", by.y = "ID")[, -7]
save(probe_modules_anno, file = "02 gene coexpression network analysis/probe_modules_anno.RData")
write.csv(probe_modules_anno, file = "02 gene coexpression network analysis/probe_modules_anno.csv")

# Relating modules to external clinical traits
# Define numbers of genes and samples
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)
# Recalculate MEs with color labels
MEs0 = moduleEigengenes(datExpr, new_moduleLabels)$eigengenes
MEs = orderMEs(MEs0)
names(MEs) = str_remove(names(MEs), "ME")
# recluster
METree = hclust(as.dist(1-cor(MEs)), method = "average");
plot(METree, main = "Clustering of module eigengenes", xlab = "", sub = "")
# divide modules into three groups 
abline(h = 1, col = "red")

# datTraits
datTraits = read.csv("02 gene coexpression network analysis/datTraits.csv")
datTraits = datTraits[match(rownames(datExpr), datTraits$GSM_ID), -2]
datTraits = data.frame(Sham_KA = datTraits$Sham_KA)
rownames(datTraits) = rownames(datExpr)
save(datTraits, file = "02 gene coexpression network analysis/datTraits.RData")

# calculate correlation between module eigen-genes and sample phenotype
MEs_colored = orderMEs(moduleEigengenes(datExpr, moduleColors)$eigengenes)
color_seq = str_remove(names(MEs_colored), "ME")
c2m = paste(color_seq, probe_modules$Modules[match(color_seq, probe_modules$Colors)], sep = ".")

moduleTraitCor = cor(MEs_colored, datTraits, use = "p")
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)

# Display the correlation values within a heatmap plot
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "")
textMatrix = as.data.frame(textMatrix)
par(mar = c(2, 8.5, 1, 3))
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = "S",
               yLabels = names(MEs_colored),
               ySymbols = c2m,
               colorLabels = T,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = NULL)

# Gene-Trait correlation for epilepsy
Gene_Trait_Cor = data.frame(GTC = cor(datExpr, datTraits, use = "p"), 
                            p_GTC = corPvalueStudent(cor(datExpr, datTraits, use = "p"), nSamples))
colnames(Gene_Trait_Cor) = c("GTC", "p_GTC")
save(Gene_Trait_Cor, file = "02 gene coexpression network analysis/Gene_Trait_Cor.RData")

# intramodular connectivity
probe_connect = intramodularConnectivity(adjacency, moduleColors, scaleByMax = FALSE)
connect_scaled = intramodularConnectivity(adjacency, moduleColors, scaleByMax = TRUE)
probe_connect$probe_ID = rownames(probe_connect)
probe_connect$scaled_kWithin = connect_scaled$kWithin

save(probe_connect, file = "02 gene coexpression network analysis/probe_connect.RData")

probe_connect_GTC = cbind(probe_connect, Gene_Trait_Cor)
probe_connect_GTC_anno = merge(probe_connect_GTC, probe_modules_anno, by = "probe_ID")

save(probe_connect_GTC_anno, file = "02 gene coexpression network analysis/probe_connect_GS_anno.RData")

M2_M13_hubs = filter(probe_connect_GTC_anno, Modules %in% c(2, 13), GTC >= 0.5, p_GTC < 0.05, scaled_kWithin >= 0.3)
M2_M13_hubs$Modules = as.factor(M2_M13_hubs$Modules)

write.csv(M2_M13_hubs, file = "02 gene coexpression network analysis/M2_M13_hubs.csv")

ggplot(M2_M13_hubs, aes(x=scaled_kWithin, y=GTC, color = Modules)) + 
  geom_point(shape = 21) +
  geom_smooth(method=lm, se=FALSE, fullrange=TRUE) +
  theme_classic() + 
  annotate(geom="text", x=0.74, y=0.8, label="Bdnf", color="#00BFC4") +
  annotate(geom="text", x=0.42, y=0.82, label="Gfap", color="#F8766D") +
  annotate(geom="text", x=1, y=0.77, label="Scg2", color="#00BFC4") +
  annotate(geom="text", x=0.81, y=0.76, label="Vgf", color="#00BFC4") +
  annotate(geom="text", x=0.98, y=0.65, label="Cxcl16", color="#F8766D") +
  annotate(geom="text", x=0.95, y=0.68, label="Tlr2", color="#F8766D") +
  annotate(geom="text", x=0.47, y=0.88, label="Gpnmb", color="#F8766D") +
  annotate(geom="text", x=0.34, y=0.9, label="C4b", color="#F8766D") +
  annotate(geom="text", x=0.37, y=0.85, label="Gldn", color="#F8766D") +
  annotate(geom="text", x=0.59, y=0.79, label="Nptx2", color="#00BFC4") +
  annotate(geom="text", x=0.87, y=0.74, label="Npy", color="#00BFC4") + 
  annotate(geom="text", x=0.84, y=0.64, label="Kcnk1", color="#00BFC4") +
  annotate(geom="text", x=0.83, y=0.71, label="Ncam2", color="#00BFC4") +
  annotate(geom="text", x=0.95, y=0.68, label="Cd52", color="#F8766D")




