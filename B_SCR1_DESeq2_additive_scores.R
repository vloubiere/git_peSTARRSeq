setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
require(data.table)
require(DESeq2)

#----------------------------------------------#
# 1- Load counts
#----------------------------------------------#
dat <- data.table(file= list.files("/groups/stark/vloubiere/data/pe_STARRSeq/Rdata", "SCR1.*UMI", full.names= T))
dat[, sample:= gsub("libvl002_SCR1_|_UMI.rds", "", basename(file))]
dat <- dat[, readRDS(file), .(file, sample)]
dat <- dcast(dat, enh_L+enh_R~sample, value.var = "counts", fill = 0)
dat[, input_rep1:= as.integer(rowSums(.SD, na.rm= T)), .SDcols= c("input_rep1", "input_rep2")]
dat$input_rep2 <- NULL
colnames(dat)[11] <- "input_rep2"
dat[, input_merge:= rowSums(.SD, na.rm= T), .SDcols= grep("input", colnames(dat))]
dat[, DSCP_merge:= rowSums(.SD, na.rm= T), .SDcols= grep("DSCP", colnames(dat))]
saveRDS(dat, "/groups/stark/vloubiere/data/pe_STARRSeq/Rdata/libvl002_SCR1_raw_counts.rds")
dat <- dat[input_merge >= 30]
saveRDS(dat, "/groups/stark/vloubiere/data/pe_STARRSeq/Rdata/libvl002_SCR1_clean_counts_table_final.rds")

#----------------------------------------------#
# 2- DESeq 2
#----------------------------------------------#
# Format tables
sampleTable <- data.table(sample= grep("rep", colnames(dat), value= T))
sampleTable[, c("condition", "replicate") := tstrsplit(sample, "_")]
sampleTable <- data.frame(sampleTable[, -1], row.names= sampleTable$sample)
ctls_idx <- which(grepl("^control", dat$enh_L) & grepl("^control", dat$enh_R) & dat$enh_L != dat$enh_R)
cols <- rownames(sampleTable)
DF <- data.frame(dat[, ..cols], row.names = dat[, paste0(enh_L, ".", enh_R), .(enh_L, enh_R)][, V1])
DF <- DF+1
# Differential analysis
dds <- DESeqDataSetFromMatrix(countData= DF, colData= sampleTable, design= ~replicate+condition)
# sizeFactors(dds) <- estimateSizeFactorsForMatrix(as.matrix(DF[ctls_idx,])+0.1)
sizeFactors(dds) <- estimateSizeFactorsForMatrix(as.matrix(DF[ctls_idx,]))
res <- DESeq(dds)
diff <- as.data.table(as.data.frame(lfcShrink(res, contrast= c("condition", "DSCP", "input"))), keep.rownames= T)
diff <- diff[, c("enh_L", "enh_R") := tstrsplit(rn, "[.]")]
diff <- diff[, .(enh_L, enh_R, baseMean, log2FoldChange, padj)]
boxplot(diff[ctls_idx, log2FoldChange])
abline(h= 0, lty= 2)
saveRDS(diff, "Rdata/B_DESeq2_control_pairs.rds")

#----------------------------------------------#
# 3- Robust controls
#----------------------------------------------#
# Identify robust controls producing few outliers
ctls <- diff[grepl("control", enh_L) | grepl("control", enh_R), .(enh_L, enh_R, log2FoldChange)]
ctls[grepl("control", enh_L), outlier_ctl_L:= log2FoldChange<boxplot(log2FoldChange, plot=F)$stats[1,1] | log2FoldChange>boxplot(log2FoldChange, plot=F)$stats[5,1], enh_R]
ctls[grepl("control", enh_L), outlier_freq_L:= length(which(outlier_ctl_L=="TRUE"))/.N, enh_L]
ctls[grepl("control", enh_R), outlier_ctl_R:= log2FoldChange<boxplot(log2FoldChange, plot=F)$stats[1,1] | log2FoldChange>boxplot(log2FoldChange, plot=F)$stats[5,1], enh_L]
ctls[grepl("control", enh_R), outlier_freq_R:= length(which(outlier_ctl_R=="TRUE"))/.N, enh_R]
# PCC LR
ctls[ctls, log2FoldChange_rev:= i.log2FoldChange, on= c("enh_L==enh_R", "enh_R==enh_L")]
ctls[grepl("control", enh_L) & !is.na(log2FoldChange_rev), PCC_L:= if(.N>5){cor.test(log2FoldChange, log2FoldChange_rev)$estimate}else{as.numeric(NA)}, enh_L]
ctls[grepl("control", enh_R) & !is.na(log2FoldChange_rev), PCC_R:= if(.N>5){cor.test(log2FoldChange, log2FoldChange_rev)$estimate}else{as.numeric(NA)}, enh_R]
# sel
L_ctls <- unique(ctls[grepl("control", enh_L) & outlier_freq_L<0.025 & PCC_L>0.8, enh_L])
R_ctls <- unique(ctls[grepl("control", enh_R) & outlier_freq_R<0.025 & PCC_R>0.8, enh_R])

#----------------------------------------------#
# 4- Additive scores
#----------------------------------------------#
# Compute expected scores
exp <- copy(diff)
# Individual activity
exp[!grepl("control", enh_L), c("all_L", "median_L", "N_L"):= .SD[enh_R %in% R_ctls, .(list(log2FoldChange), median(log2FoldChange), .N)], enh_L]
exp[!grepl("control", enh_R), c("all_R", "median_R", "N_R"):= .SD[enh_L %in% L_ctls, .(list(log2FoldChange), median(log2FoldChange), .N)], enh_R]
exp <- na.omit(exp[N_L>4 & N_R>4, !c("N_L", "N_R")])
# Method 1
# exp[, c("log2FCs_add_ls", "log2FCs_add_box"):= 
#       {
#         current <- log2(rowSums(CJ(2^unlist(all_L), (2^unlist(all_R)))))
#         list(.(quantile(current, seq(0, 1, length.out = 101))), .(boxplot(current, plot = F)$stats[, 1]))
#       }, .(enh_L, enh_R)]
# exp[, log2FC_add:= sapply(log2FCs_add_ls, median)]
# exp[, log2FC_add_perc:= mapply(function(x, y){length(which(x>y))/length(y)}, x= log2FoldChange, y= log2FCs_add_ls, SIMPLIFY = T)]
# exp[, diff:= log2FoldChange-log2FC_add]
# Method 2
exp[, log2FC_add:= log2((2^median_L)+(2^median_R)), .(enh_L, enh_R)]
exp[, diff:= log2FoldChange-log2FC_add]
# Add extra stats
exp[exp, log2FoldChange_rev:= i.log2FoldChange, on= c("enh_L==enh_R", "enh_R==enh_L")]
# SAVE short table
clean <- exp[diff, , on=c("enh_L", "enh_R")]
clean <- unique(clean[, .(enh_L, enh_R, baseMean= i.baseMean, log2FoldChange= i.log2FoldChange, padj= i.padj, log2FoldChange_rev,
                          median_L, median_R, log2FC_add, diff)])
saveRDS(clean, "Rdata/B_SCR1_peSTARRSeq_final_table.rds")

