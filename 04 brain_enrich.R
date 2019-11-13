library(WGCNA)
library(ggplot2)
library(stringr)

load("03 Gene coexpression network/gene_modules.RData")

science_singCell = userListEnrichment(gene_modules$Gene_symbol, labelR = gene_modules$Modules, fnIn = "04 brain cell type enrichment/science_singleCell_brainCell_types.csv", catNmIn = "")

science_singCell_pval = science_singCell$pValues
science_singCell_pval$cell_type = str_remove(science_singCell_pval$UserDefinedCategories, "__")
# significant enrichment with P-value < 0.05
sig_enrich_science = science_singCell_pval[science_singCell_pval$Pvalues < 0.05, ]

sig_enrich_science_0.11 = science_singCell_pval[science_singCell_pval$Pvalues < 0.11, ]
sig_enrich_science_0.11 = sig_enrich_science_0.11[order(sig_enrich_science_0.11$CorrectedPvalues), ]

sig_enrich_science_0.11$Signifcance <- factor(rep(c("yes", "no"), c(11,11)), levels = c("yes", "no"))

save(science_singCell_pval, sig_enrich_science, file = "04 brain cell type enrichment/science_singCell_pval.RData")

# draw the plot
ggplot(sig_enrich_science_0.11, 
       aes(x = InputCategories, y = cell_type, size = NumOverlap, fill = -log10(CorrectedPvalues), color = Signifcance)) + 
  geom_point(shape = 21) + 
  scale_fill_continuous(low = "white", high = "red") + 
  scale_color_manual(values = c("red", "gray10")) +
  scale_size(range = c(1, 10)) + 
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), 
        plot.margin = margin(1, 0, 0.5, 0, "cm")) +  
  scale_y_discrete(limits=c("Interneuron","CA1 Pyramidal Neuron","S1 Pyramidal Neuron", "Astrocyte",  "Microglia", "Ependymal cell", "Oligodendrocyte","Mural cell", "Vascular endothelial cell")) + 
  scale_x_discrete(limits=c("M10", "M11", "M12", "M13", "M1", "M2", "M3", "M4","M5", "M6", "M7", "M8", "M9"))
