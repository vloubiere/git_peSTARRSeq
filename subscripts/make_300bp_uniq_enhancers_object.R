# Merge vl BA libraries ####
# Merge BA and my TWIST into a single non-redundant clean object
BA <- fread("Rdata/BA_300bp_TWIST_STARRSeq.txt")
vl8 <- as.data.table(readRDS("Rdata/vl_library_twist008_112019.rds"))
vl12 <- as.data.table(readRDS("Rdata/vl_library_twist12_210610.rds"))
vl12[BA, ID_BA:= i.ID, on= c("seqnames", "start", "end", "strand")]
lib <- merge(BA[, .(ID_BA= ID, 
                    BA_group= Group, 
                    BA_enhancer_group= enhancer_group, 
                    BA_enh_group_detail= enhancer_group_detail, 
                    dev_log2FoldChange, 
                    hk_log2FoldChange, 
                    seqnames, 
                    start, 
                    end, 
                    width, 
                    strand)], 
             vl8[, .(ID_BA= BA_ID, 
                     ID_twist08= ID_vl, 
                     group_twist08= group, 
                     detail_twist08= detail, 
                     linker_twist08= linker_ID, 
                     seqnames, 
                     start, 
                     end, 
                     width, 
                     strand)], all.x= T, all.y= T)
lib <- merge(lib, 
             vl12[, .(ID_BA, 
                      ID_twist12= ID,  
                      group_twist12= group, 
                      detail_twist12= detail, 
                      linker_twist12= linker_ID, 
                      seqnames, 
                      start, 
                      end, 
                      width= 249, 
                      strand)], all.x= T, all.y= T)

# Uniq groups
lib[BA_enhancer_group=="shared", uniq_group:= ifelse(dev_log2FoldChange>hk_log2FoldChange, "dev", "hk")]
lib[BA_enhancer_group=="dev", uniq_group:= "dev"]
lib[BA_enhancer_group=="hk", uniq_group:= "hk"]
lib[BA_enhancer_group=="Inducible", uniq_group:= "inducible"]
lib[BA_enhancer_group=="Controls", uniq_group:= "control"]
lib[group_twist08=="OSC", uniq_group:= "OSC"]
lib[group_twist08=="control", uniq_group:= "control"]
lib[group_twist12=="Repressor", uniq_group:= "Silencer"]
lib[group_twist12=="SUHW_peak", uniq_group:= "SUHW_peak"]
lib[group_twist12=="DHS_peak", uniq_group:= "DHS_peak"]
lib[group_twist12=="CP", uniq_group:= "CP"]

# Uniq details
lib[uniq_group== "dev", 
    uniq_detail:= cut(dev_log2FoldChange, c(-Inf,4,6,8,Inf), labels = c("inactive", "weak", "medium", "strong"))]
lib[uniq_group== "hk", 
    uniq_detail:= cut(hk_log2FoldChange, 
                      quantile(hk_log2FoldChange, seq(0, 1, length.out= 5)), 
                      include.lowest = T, 
                      labels = c("inactive", "weak", "medium", "strong"))]
lib[uniq_group== "inducible" & BA_enh_group_detail== "Ecdysone", uniq_detail:= "ecdysone"]
lib[uniq_group== "inducible" & BA_enh_group_detail== "HeatShock", uniq_detail:= "heatshock"]
lib[uniq_group== "inducible" & BA_enh_group_detail== "Cadmium", uniq_detail:= "cadmium"]
lib[uniq_group== "inducible" & BA_enh_group_detail== "Paraquat", uniq_detail:= "paraquat"]
lib[uniq_group== "inducible" & BA_enh_group_detail== "PGN", uniq_detail:= "PGN"]
lib[uniq_group== "inducible" & BA_enh_group_detail== "Wnt", uniq_detail:= "wnt"]
lib[uniq_group== "OSC", uniq_detail:= "OSC_specific"]
lib[uniq_group== "control" & detail_twist08== "exon", uniq_detail:= "exon"]
lib[uniq_group== "control" & detail_twist08== "Ecoli", uniq_detail:= "ecoli"]
lib[uniq_group== "control" & is.na(uniq_detail), uniq_detail:= "flat_genomic_region"]
lib[uniq_group== "DHS_peak" & is.na(uniq_detail), uniq_detail:= "DHS_peak"]
lib[uniq_group== "Silencer" & is.na(uniq_detail), uniq_detail:= "Silencer"]
lib[uniq_group== "SUHW_peak" & is.na(uniq_detail), uniq_detail:= "SUHW_peak"]
lib[uniq_group== "CP", uniq_detail:= detail_twist12]

# colors
clean <- lib[, .(ID_BA,
                 ID_twist08,
                 ID_twist12,
                 group= uniq_group, 
                 detail= uniq_detail, 
                 linker_twist08,
                 linker_twist12,
                 seqnames, 
                 start, 
                 end, 
                 strand, 
                 dev_log2FoldChange, 
                 hk_log2FoldChange)]
class_Cc <- data.table(group= c("hk", "dev", "OSC", "inducible", "control", "CP", "DHS_peak", "Silencer", "SUHW_peak"),
                       col= c("tomato", "#74C27A", "black", "gold", "lightgrey", "cyan3", "royalblue2", "deeppink2", "darkorchid3"))
clean <- clean[class_Cc, , on= "group"]

# SAVE
saveRDS(clean, "Rdata/uniq_300bp_enhancers.rds")
