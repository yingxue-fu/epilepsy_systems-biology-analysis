library(biomaRt)
library(dplyr)
library(WGCNA)
library(ggplot2)
library(stringr)

load("02 gene coexpression network analysis/probe_modules_anno.RData")

science_singCell = userListEnrichment(probe_modules_anno$Gene_symbol, labelR = probe_modules_anno$Modules, fnIn = "03 brain cell type enrichment/science_singleCell_brainCell_types.csv", catNmIn = "")

science_singCell_pval = science_singCell$pValues
science_singCell_pval$cell_type = str_remove(science_singCell_pval$UserDefinedCategories, "__")
# significant enrichment with P-value < 0.05
sig_enrich_science = science_singCell_pval[science_singCell_pval$Pvalues < 0.05, ]

# draw the plot
ggplot(sig_enrich_science, 
       aes(x = InputCategories, 
           y = cell_type, 
           size = NumOverlap, 
           fill = Pvalues)) + 
  geom_point(shape = 21) + 
  scale_fill_continuous(low = "red", high = "white") + 
  scale_size(range = c(1, 10)) + 
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), 
        plot.margin = margin(1, 0, 0.5, 0, "cm"))

