setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
sapply(list.files("/groups/stark/vloubiere/functions/", ".R$", full.names = T), source)
require(data.table)

# retrieve features and identify sgl short enhancer ID
constructs <- fread("/groups/stark/vloubiere/exp_data/constructs_sequences.txt", key= "name")
feat <- readRDS("Rdata/library/lib_features.rds")
feat[, c("seqnames", "start", "end"):= tstrsplit(coor, ":|-", keep= 1:3)]
cols= c("start", "end");
feat[, (cols):= lapply(.SD, as.numeric), .SDcols= cols]
sgl_ID <- feat[constructs["sgl_LH"], ID, on= c("seqnames", "start<end", "end>start")]

# peSTARR-Seq data
dat<- readRDS("Rdata/processed_peSTARRSeq_data/SCR1_peSTARRSeq_final_table.rds")

# Merged table containing, the activity of sgl/enh in pe-STARR-Seq (sgl in left position!!) and rep-STARR-Seq
res <- dat[enh_L==sgl_ID, .(ID, enh= enh_R, peSTARR_log2FC= log2FoldChange)]
res[feat, repSTARR_log2FC:= i.sgl_repSTARR_log2FoldChange, on= "enh==ID"]


plot(res$peSTARR_log2FC, res$repSTARR_log2FC)
