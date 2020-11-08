setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
sapply(list.files("/groups/stark/vloubiere/functions/", ".R", full.names = T), source)
require(data.table)
require(TxDb.Dmelanogaster.UCSC.dm3.ensGene)
require(BSgenome.Dmelanogaster.UCSC.dm3)
require(org.Dm.eg.db)
require(GenomicRanges)

#------------------------------------------#
# Merge BA and my TWIST into a single non-redundant clean object
#------------------------------------------#
vl <- as.data.table(readRDS("Rdata/library/vl_library_112019.rds"))
BA <- fread("Rdata/library/BA_300bp_TWIST_STARRSeq.txt")
lib <- merge(BA[, .(ID_BA= ID, BA_group= Group, BA_enhancer_group= enhancer_group, BA_enh_group_detail= enhancer_group_detail, 
                    dev_log2FoldChange, hk_log2FoldChange, seqnames, start, end, width, strand)], 
             vl[, .(ID_BA= BA_ID, ID_vl, group, detail, linker_ID, seqnames, start, end, width, strand)], all.x= T, all.y= T)
# Uniq IDs
lib[, vl:= ifelse(!is.na(ID_vl), T, F)]
lib[, uniq_id:= ifelse(is.na(ID_vl), as.character(ID_BA), as.character(ID_vl))]
# Uniq groups
lib[BA_enhancer_group=="dev", uniq_group:= "dev"]
lib[BA_enhancer_group=="hk", uniq_group:= "hk"]
lib[BA_enhancer_group=="Inducible", uniq_group:= "inducible"]
lib[BA_enhancer_group=="shared", uniq_group:= "shared"]
lib[BA_enhancer_group=="Controls", uniq_group:= "control"]
lib[group=="OSC", uniq_group:= "OSC"]
lib[group=="control", uniq_group:= "control"]
# Uniq details
lib[uniq_group== "dev", uniq_detail:= cut(dev_log2FoldChange, c(-Inf,4,6,8,10), labels = c("inactive", "weak", "medium", "strong"))]
lib[uniq_group== "shared", uniq_detail:= cut(dev_log2FoldChange, c(-Inf,4,6,8,10), labels = c("inactive", "weak", "medium", "strong"))]
lib[uniq_group== "hk", uniq_detail:= cut(hk_log2FoldChange, quantile(hk_log2FoldChange, seq(0, 1, length.out= 5)), include.lowest = T, labels = c("inactive", "weak", "medium", "strong"))]
lib[uniq_group== "inducible" & BA_enh_group_detail== "Ecdysone", uniq_detail:= "ecdysone"]
lib[uniq_group== "inducible" & BA_enh_group_detail== "HeatShock", uniq_detail:= "heatshock"]
lib[uniq_group== "inducible" & BA_enh_group_detail== "Cadmium", uniq_detail:= "cadmium"]
lib[uniq_group== "inducible" & BA_enh_group_detail== "Paraquat", uniq_detail:= "paraquat"]
lib[uniq_group== "inducible" & BA_enh_group_detail== "PGN", uniq_detail:= "PGN"]
lib[uniq_group== "inducible" & BA_enh_group_detail== "Wnt", uniq_detail:= "wnt"]
lib[uniq_group== "OSC", uniq_detail:= "OSC_specific"]
lib[uniq_group== "control" & detail== "exon", uniq_detail:= "exon"]
lib[uniq_group== "control" & detail== "Ecoli", uniq_detail:= "ecoli"]
lib[uniq_group== "control" & is.na(uniq_detail), uniq_detail:= "flat_genomic_region"]

# colors
clean <- lib[, .(ID= uniq_id, group= uniq_group, detail= uniq_detail, vl, linker_ID, coor= as.character(GRanges(lib[, seqnames:strand])), dev_log2FoldChange, hk_log2FoldChange)]
class_Cc <- data.table(group= c("hk", "shared", "dev", "OSC", "inducible", "control"),
                       col= c("tomato", "royalblue2", "#74C27A", "black", "gold", "lightgrey"))
clean <- clean[class_Cc, , on= "group"]

# SAVE
saveRDS(clean, "Rdata/library/uniq_library_final.rds")
