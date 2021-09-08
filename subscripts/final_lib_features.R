setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")

lib <- as.data.table(readRDS("Rdata/300bp_enhancer_chromatin_features_and_gene_assignment.rds"))

#---------------------------------#
# Retrieve counts for most informative motifs
#---------------------------------#
mot <- readRDS("Rdata/informatrive_motifs_300bp_enhancers.rds")
lib <- cbind(lib, 
             mot)
motSom <- readRDS("Rdata/informatrive_motifs_300bp_enhancers_SOM.rds")
lib <- cbind(lib, 
             motSom)


# SAVE
setcolorder(lib, c("ID_BA",
                   "ID_twist08",
                   "ID_twist12",
                   "group", 
                   "detail", 
                   "linker_twist08", 
                   "linker_twist12", 
                   "seqnames", 
                   "start", 
                   "end", 
                   "strand",
                   "col",
                   "closest_tss", 
                   "closest_sel_tss",
                   "dev_log2FoldChange", 
                   "hk_log2FoldChange"))
saveRDS(lib, "Rdata/final_300bp_enhancer_features.rds")










