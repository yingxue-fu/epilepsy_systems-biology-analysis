library(Biobase)
library(GEOquery)
library(WGCNA)
library(dplyr)

load("02 gene coexpression network analysis/probe_modules_anno.RData")

# find the homologous gene in human
human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
save(human, mouse, file = "05 conservation in human/human_mouse_biomaRt.RData")

module_genes_m2h = getLDS(attributes = c("entrezgene", "mgi_symbol"), 
                          filters = "entrezgene", 
                          values = probe_modules_anno$Gene_ID, 
                          mart = mouse, 
                          attributesL = c("entrezgene", "hgnc_symbol"), 
                          martL = human) %>% na.omit()

probe_modules_anno_m2h = merge(probe_modules_anno, module_genes_m2h, by.x = "Gene_ID", by.y = "NCBI.gene.ID")
save(probe_modules_anno_m2h, file = "05 conservation in human/probe_modules_anno_m2h.RData")




# load series and platform data from GEO
gset_h <- getGEO(filename = "05 conservation in human/GSE63808_series_matrix.txt.gz", GSEMatrix = TRUE, getGPL = F)

# extract gene expresssion matrix
ex_h = exprs(gset_h)
dim(ex_h)
ex_h = log2(ex_h)

# extract probe annotation from GPL data
gpl_h = getGEO(filename = "05 conservation in human/GPL6947.annot.gz")
probe_annot_h = Table(gpl_h)[, 1:4]
probe_annot_h[probe_annot_h == ""] = NA
probe_annot_h = probe_annot_h %>% na.omit()

save(probe_annot_h, file = "05 conservation in human/probes_annotation_h.RData")

# expressed probes (p-value < 0.05 in more than 25% of the samples)
load("05 conservation in human/gsm_extract.RData")
ptest <- function(x) {sum(x < 0.05)/129 > 0.25}
expressed_probes <- rownames(pvalMatrix)[apply(pvalMatrix, 1, ptest)]
expressed_probes_anno = probe_annot_h[probe_annot_h$ID %in% expressed_probes, ]

conserve_probe = expressed_probes_anno[expressed_probes_anno$`Gene ID` %in% probe_modules_anno_m2h$NCBI.gene.ID.1, ]

save(conserve_probe, file = "05 conservation in human/conserve_probes.RData")

length(unique(conserve_probe$`Gene ID`))


