library(WGCNA)
library(dplyr)
library(stringr)
library(ggrepel)

load("02 Gene meta-signature/exprs_gse73878_core_IL.RData")

datExpr0 <- exprs_gse73878_core_IL %>% as.data.frame()
datExpr <- datExpr0[,-35] %>% as.matrix(); rownames(datExpr) <- datExpr0$`Gene symbol`

save(datExpr, file = "03 Gene coexpression network/datExpr.RData")

# 1. Calculate biweight midcorrelation between genes
datExpr <- t(datExpr)
biweight_coref = bicor(datExpr)
person_coref = cor(datExpr, method = "pearson")

hist(biweight_coref[upper.tri(biweight_coref)], breaks = 50, freq = F)
hist(person_coref[upper.tri(person_coref)], breaks = 50, freq = F)

# 2. determine the soft-thresholding powers
# convert negative correlation to 0 for "signed hybrid" network type
biweight_coref_signed = biweight_coref 
biweight_coref_signed[biweight_coref_signed < 0] = 0

# Network topology analysis to choose the soft power
powers = c(seq(1, 10, by = 1), seq(12, 20, by = 2))
sft = pickSoftThreshold.fromSimilarity(biweight_coref_signed, blockSize = 9000, verbose = 5)

par(mfrow = c(1,2)); cex1 = 0.9
# Power   SFT.R.sq    slope    truncated.R.sq   mean.k.   median.k.   max.k.
# 4        0.902     -1.64        0.984         67.700     44.8000    389.0

# Scale-free fit indices as a function of the soft-thresholding power 
plot(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3])*sft$fitIndices[, 2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2", 
     type="n", main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
abline(h = 0.9, col = "red")

# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[, 1], sft$fitIndices[, 5],
     xlab = "Soft Threshold (power)",ylab = "Mean Connectivity", type = "n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[, 1], sft$fitIndices[, 5], labels = powers, cex = cex1,col = "red")
abline(h = 67, col = "red")

# Soft power is chosed to be 4
softPower = 4

# Claculate ajacency matrix
adjacency = adjacency.fromSimilarity(biweight_coref, 
                         type = "signed hybrid", 
                         power = 4)
save(adjacency, file = "03 Gene coexpression network/adjacency.RData")


fNetCon <- fundamentalNetworkConcepts(adj = adjacency)

scaleFreePlot(
  fNetCon$Connectivity, 
  nBreaks = 12, 
  truncated = T,  
  removeFirst = F, 
  pch = 20, col = "gray50")

scaleFreeFitIndex(fNetCon$Connectivity, nBreaks = 10, removeFirst = FALSE)

hist(fNetCon$ScaledConnectivity, breaks = 50)

# Connectivity vs. Cluster Coefficient
fNetCon_df <- data.frame(Connectivity = fNetCon$Connectivity, 
                         ScaledConnectivity = fNetCon$ScaledConnectivity, 
                         ClusterCoef = fNetCon$ClusterCoef)
fNetCon_df$Gene_symbol <- rownames(fNetCon_df)

fNetCon_df <- merge(fNetCon_df, gene_modules, by = "Gene_symbol")
Module_Color <- unique(fNetCon_df[,c(5,6)])
Module_Color <- Module_Color[order(Module_Color$Modules),]

ggplot(fNetCon_df, aes(x=Connectivity, y=ClusterCoef, color=Modules)) +
  geom_point(alpha=0.5) +
  scale_color_manual(values = as.character(Module_Color$Colors))+
  theme_bw()
  

############################## Module Detection #############################
# Calculate TOM matrix
TOM = TOMsimilarity(adjacency)
save(TOM, file = "03 Gene coexpression network/TOM.RData")

# Cluster genes based on TOM similarity
dissTOM = 1-TOM
geneTree = hclust(as.dist(dissTOM), method = "average")
plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity", labels = FALSE, hang = 0.04)

save(geneTree, file = "03 Gene coexpression network/geneTree.RData")

# Set the minimum module size
minModuleSize = 30

# Module identification using dynamic tree cut:
original_Modules = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                            deepSplit = 2, pamRespectsDendro = FALSE,
                            minClusterSize = minModuleSize)
head(original_Modules)
table(original_Modules)

# Convert numeric lables into colors
original_Colors = labels2colors(original_Modules, colorSeq = c(brewer.pal(n = 12, name = "Set3"), brewer.pal(n = 7, name = "Dark2")))
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

colorOrder = data.frame(color = c("grey", brewer.pal(n = 12, name = "Set3"), brewer.pal(n = 7, name = "Dark2")), order1 = 0:19)
MEcolor2module = colorOrder$order1[match(str_remove(names(MEs), "ME"), colorOrder$color)]
names(MEs) = paste("M", MEcolor2module, sep = "")

# Module cluster based on their eigengenes
METree = hclust(as.dist(1-cor(MEs)), method = "average")
plot(METree, main = "Clustering of original module eigengenes", xlab = "", sub = "")
abline(h = 0.2, col = "red")

# Merge similar modules
merge = mergeCloseModules(datExpr, original_Modules, cutHeight = 0.2, colorSeq = c("grey", brewer.pal(n = 12, name = "Set3"), brewer.pal(n = 7, name = "Dark2")), verbose = 3)
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
ordered_moduleLabels$orderedLabels = 0:13

new_moduleLabels = paste("M", ordered_moduleLabels$orderedLabels[match(moduleLabels, ordered_moduleLabels$moduleLabels)], sep = "")
table(new_moduleLabels)

save(moduleColors, moduleLabels, new_moduleLabels, file = "03 Gene coexpression network/module_color_labels.RData")

# gene cluster and final modules
plotDendroAndColors(geneTree, moduleColors, groupLabels = "Modules", rowText = new_moduleLabels, rowTextAlignment = "center", addTextGuide = T, dendroLabels = FALSE, hang = 0.03, addGuide = F, guideHang = 0.1)

# save probe_modules data
gene_modules = data.frame(Gene_symbol = colnames(datExpr), Modules = new_moduleLabels, Colors = moduleColors)
save(gene_modules, file = "03 Gene coexpression network/gene_modules.RData")
write.csv(gene_modules, file = "03 Gene coexpression network/gene_modules.csv")

table(gene_modules$Modules)
#    M0   M1  M10  M11  M12  M13   M2   M3   M4   M5   M6   M7   M8   M9 
#    12 1896  134   70   63   44 1533 1701  866  443  680  383  284  275 

# intramodular connectivity
gene_connect = intramodularConnectivity(adjacency, moduleColors, scaleByMax = FALSE)
connect_scaled = intramodularConnectivity(adjacency, moduleColors, scaleByMax = TRUE)
gene_connect$Gene_symbol = rownames(gene_connect)
gene_connect$scaled_kWithin = connect_scaled$kWithin

gene_modules_connect <- merge(gene_modules, gene_connect[, c(1,2,5,6)])

save(gene_connect, gene_modules_connect, file = "03 Gene coexpression network/gene_connect.RData")


############################### metaDEGs distribution ##################################

load(file = "02 Gene meta-signature/metaRP.DEG.RData")
metaRP_acute_DEG$Modules <- gene_modules$Modules[match(metaRP_acute_DEG$ID, gene_modules$Gene_symbol)]
a_m <- table(metaRP_acute_DEG$Modules) %>% as.data.frame()
a_m$Stage <- "acute"

metaRP_latent_DEG$Modules <- gene_modules$Modules[match(metaRP_latent_DEG$ID, gene_modules$Gene_symbol)]
l_m <- table(metaRP_latent_DEG$Modules) %>% as.data.frame()
l_m$Stage <- "latent"

metaRP_chronic_DEG$Modules <- gene_modules$Modules[match(metaRP_chronic_DEG$ID, gene_modules$Gene_symbol)]
c_m <- table(metaRP_chronic_DEG$Modules) %>% as.data.frame()
c_m$Stage <- "chronic"

stage_module <- rbind(a_m, l_m, c_m)

ggplot(data=stage_module, aes(x=Var1, y=Freq, fill=Stage)) +
  geom_bar(stat="identity", position=position_dodge())

