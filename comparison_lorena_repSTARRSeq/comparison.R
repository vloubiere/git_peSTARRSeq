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
res <- dat[enh_L==sgl_ID, .(ID, enh= enh_R, peSTARR_log2FC= log2FoldChange, peSTARR_exp= log2FC_add)]
res[feat, c("repSTARR_log2FC", "repSTARR_exp"):= .(i.sgl_repSTARR_log2FoldChange, log2(2^i.dev_log2FoldChange+2^feat[ID==sgl_ID, dev_log2FoldChange])-feat[ID==sgl_ID, dev_log2FoldChange]), on= "enh==ID"]
res[, c("peSTARR-Seq_diff", "repSTARR-Seq_diff") := .(peSTARR_log2FC-peSTARR_exp, repSTARR_log2FC-repSTARR_exp)]

# Compare peSTARR-Seq and rep-STARR-Seq
mres <- melt(res, id.vars = "ID", measure.vars = c("peSTARR-Seq_diff", "repSTARR-Seq_diff"))
mres[, variable:= as.character(variable)]
mres[, variable:= switch(variable, "peSTARR-Seq_diff"= "CP, enhA, enhB", "repSTARR-Seq_diff"= "enhA, CP, enhB"), variable]

pdf("pdf/comparison_repSTARRSeq_lorena/boxplot_observed_expected.pdf", 4, 7)
par(mar=c(10,5,5,5))
my_boxplot(value~variable, mres, ylab = "Observed vs Expected additive model (log2)", las= 2, pval_list = list(c(1,2)), ylim= c(-2,3.5))
abline(h=0, lty= 2)
dev.off()
