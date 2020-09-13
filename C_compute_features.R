setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
source("/groups/stark/vloubiere/scripts/R_functions/R_shell_cmd_wrap_1.0.R")
require(data.table)
require(TxDb.Dmelanogaster.UCSC.dm3.ensGene)
require(BSgenome.Dmelanogaster.UCSC.dm3)
require(org.Dm.eg.db)
require(GenomicRanges)

vl <- as.data.table(readRDS("Rdata/vl_library_112019.rds"))
BA <- fread("Rdata/BA_300bp_TWIST_STARRSeq.txt")
lib <- merge(BA[, .(ID_BA= ID, BA_group= Group, BA_enhancer_group= enhancer_group, BA_enh_group_detail= enhancer_group_detail, 
                    twist_log2FC= dev_log2FoldChange, twist_padj= dev_padj, hk_twist_log2FC= hk_log2FoldChange, hk_twist_padj= hk_padj, 
                    seqnames, start, end, width, strand)], 
             vl[, .(ID_BA= BA_ID, ID_vl, group, detail, seqnames, start, end, width, strand)], all.x= T, all.y= T)
lib <- lib[seqnames %in% c("chr2L", "chr2R", "chr3L", "chr3R", "chr4", "chrX")]
lib <- lib[order(seqnames, start, end)]
lib[, coor:= as.character(GRanges(lib))]
lib[, uniq_ID:= ifelse(is.na(ID_vl), as.character(ID_BA), ID_vl)]
lib <- lib[, .(coor, uniq_ID, ID_BA, BA_group, BA_enhancer_group, BA_enh_group_detail, 
             ID_vl, group, detail, twist_log2FC, twist_padj, hk_twist_log2FC, hk_twist_padj)]

### Gene assignment #################################
tss <- as.data.table(resize(genes(TxDb.Dmelanogaster.UCSC.dm3.ensGene), 1, fix = "start"))
tss <- tss[seqnames %in% c("chr2L", "chr2R", "chr3L", "chr3R", "chr4", "chrX")]
tss[, symbol:= mapIds(org.Dm.eg.db, key= gene_id, column="SYMBOL", keytype="FLYBASE", multiVals="first")]
tss <- GRanges(tss[order(seqnames, start, end)], name= tss$symbol)
tss_mult <- tss[tss$name %in% c("ush", "shn", "ct", "InR", "Eip75B", "Mur2B", "Smr", "brat", "kay", "chinmo")]
# Temp files
tmp1 <- tempfile(fileext = ".bed")
export(GRanges(lib$coor, name= lib$uniq_ID), tmp1)
tmp2 <- tempfile(fileext = ".bed")
export(tss, tmp2)
tmp3 <- tempfile(fileext = ".bed")
export(tss_mult, tmp3)
# All genes
all <- fread(cmd= closestBed(tmp1, tmp2))
all <- all[, .(V10= V10[1]), V1:V4]
mult <- fread(cmd= closestBed(tmp1, tmp3))
mult <- mult[, .(V10= V10[1]), V1:V4]
# Add to lib
lib[all, ctss:= i.V10 , on= "uniq_ID==V4"]
lib[mult, ctss_mult:= i.V10 , on= "uniq_ID==V4"]

### Compute gw STARR-Seq enrichment ###############################
STARR <- data.table(file= c("/groups/stark/vloubiere/data/gw_STARRSeq_bernardo/raw_data/input_DSCP_200bp_cut.bed",
                            "/groups/stark/vloubiere/data/gw_STARRSeq_bernardo/raw_data/DSCP_200bp_gw_Rep1.UMI_cut.bed", 
                            "/groups/stark/vloubiere/data/gw_STARRSeq_bernardo/raw_data/DSCP_200bp_gw_Rep2_A.UMI_cut.bed"), cdition= c("input", "screen", "screen"))
counts <- my_countReads(GRanges(lib$coor, name= lib$uniq_ID), STARR$file)
# Filter out low input counts & ompute enrichment
STARR <- STARR[counts, , on= "file"]
STARR[, filter:= all(.SD[cdition=="input", counts]>50), name]
STARR <- STARR[(filter), !"filter"]
STARR <- STARR[, .(counts= sum(counts)+1, total_reads= sum(total_reads)), .(cdition, name)]
STARR[, DSCP200_log2FC:= log2(counts/total_reads)-log2(.SD[cdition=="input", counts/total_reads]), name]
# Add to lib
lib <- merge(lib, STARR[cdition=="screen", .(uniq_ID= name, DSCP200_log2FC)], all.x= T)

### Add BAC and gw rep-STARR-Seq screens ###############################
rep <- data.table(file= list.files(c("/groups/stark/vloubiere/data/rep_STARRSeq_lorena/BAC32_screens", 
                                     "/groups/stark/vloubiere/data/rep_STARRSeq_lorena/sgl_gw_200bp"), ".bed$", full.names = T))
rep[, c("assay", "cdition", "sample"):= .(ifelse(grepl("gw", file), "gw", "BAC"), 
                                          ifelse(grepl("input", file), "input", "screen"),
                                          unlist(tstrsplit(basename(file), "[.]", keep= 1))), (rep)]
rep <- rep[grepl("sgl|ham", sample)]
counts <- my_countReads(peaks_bed = GRanges(lib$coor, name= lib$uniq_ID), reads_beds_paths = rep$file)
# Filter out low input counts & ompute enrichment
rep <- rep[counts, , on= "file"]
rep[, filter:= all(.SD[cdition=="input", counts]>50), .(assay, sample, name)]
rep <- rep[(filter), !"filter"]
rep <- rep[, .(counts= sum(counts)+1, total_reads= sum(total_reads)), .(assay, cdition, sample, name)]
rep[, log2FC:= log2(counts/total_reads)-log2(.SD[cdition=="input", counts/total_reads]), .(assay, sample, name)]
rep <- rep[cdition=="screen"]
# Compute expected additive score
rep[grepl("ham", sample), stab_enh_act:= lib[uniq_ID=="dev_medium_C_00065", DSCP200_log2FC]]
rep[grepl("sgl", sample), stab_enh_act:= lib[uniq_ID=="dev_medium_B_00587", DSCP200_log2FC]]
rep <- merge(rep, lib[, .(name= uniq_ID, individual_act= DSCP200_log2FC)])
rep[, exp_add:= log2((2^stab_enh_act)+(2^individual_act))-stab_enh_act]
# Clean
rep <- melt(rep, id.vars = c("name", "sample", "assay"), measure.vars = c("log2FC", "exp_add"))
rep[, split:= paste0(c(sample, assay, as.character(variable)), collapse= "_"), (rep)]
rep <- dcast(rep, name~split, value.var = "value")
# Add to lib
lib <- merge(lib, rep, by.x= "uniq_ID", by.y= "name", all.x= T)

### AVAILABLE data ############################
avail <- fread("../../data/available_data/metadata/SampleTable_final.txt")
avail <- rbind(avail, data.table(processed_bed_file= c("/groups/stark/vloubiere/data/PROseq/bed/S2_control_1.uniq.bed",
                                                       "/groups/stark/vloubiere/data/PROseq/bed/S2_control_2.uniq.bed"),
                                                       assay= "PRO", cdition= "parental", sample= "S2", rep= c("rep1", "rep2")), fill= T)
avail[cdition=="GAFRNAi", c("cdition", "sample"):= .("untreated", "input")]
counts <- my_countReads(peaks_bed = resize(GRanges(lib$coor, name= lib$uniq_ID), 1000, "center"), reads_beds_paths = avail$processed_bed_file)
colnames(counts)[5] <- "uniq_ID"
# Filter out low cov
avail <- avail[counts, , on= "processed_bed_file==file"]
avail[, filter:= all(counts>4), .(uniq_ID, center_project_name, assay, cdition, sample)]
avail <- avail[(filter), !"filter"]
# Merge
avail <- avail[, .(counts= sum(counts), cov= sum(total_reads)), .(center_project_name, assay, cdition, sample, uniq_ID)]
# input norm
avail[, input_norm := any(sample=="input"), .(center_project_name, assay, cdition, uniq_ID)]
avail[(input_norm), log2_enr:= log2(counts/cov)-log2(.SD[sample=="input", counts/cov]), .(center_project_name, assay, cdition, uniq_ID)]
# negative regions norm
avail[!(input_norm), log2_enr:= log2(counts)-log2(median(.SD[grepl("NegativeRegions$|^control", uniq_ID), counts])), .(center_project_name, assay, cdition)]
# Clean 
avail <- avail[sample!="input"]
avail[sample=="S2", sample:= assay]
avail[, sample:= paste0(c(sample, cdition, center_project_name), collapse = "_"), .(sample, cdition, center_project_name)]
# Visual inspection
avail <- avail[!grepl("K27me3_dslacz|H3K4me1_WT_GSE41440|H3K27ac_WT|GSE96581", sample)]
# Add to lib
avail <- dcast(avail, uniq_ID~sample, value.var = "log2_enr")
lib <- merge(lib, avail, by= "uniq_ID", all.x= T)

########## AVE
saveRDS(lib, "Rdata/C_features_final_table.rds")










