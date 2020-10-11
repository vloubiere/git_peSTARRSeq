setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
require(data.table)
require(DESeq2)

#----------------------------------------------#
# 1- Load counts
#----------------------------------------------#
if(!file.exists("db/Rdata/all_uniq_counts.rds"))
{
  dat <- data.table(file= list.files("Rdata", "SCR1.*uniq.UMI.rds", full.names= T))
  dat[, sample:= gsub("libvl002_SCR1_|.uniq.UMI.rds", "", basename(file))]
  dat <- dat[, readRDS(file), .(file, sample)]
  dat <- dat[, .(counts= .N), c(colnames(dat))]
  saveRDS(dat, "db/Rdata/all_uniq_counts.rds")
}
dat <- readRDS("db/Rdata/all_uniq_counts.rds")
# Merge rep 1 and 2
mer <- data.table(old= paste0("input_rep", 1:6), new= paste0("input_rep", c(1,1,2,3,4,5)))
dat[mer, sample:= i.new, on= "sample==old"]
dat <- dat[, .(counts= sum(counts)), .(sample, enh_L, enh_R)]
# dcast and low counts cutoff
mat <- dcast(dat, enh_L+enh_R~sample, value.var = "counts", fill = 0)
mat[, check:= rowSums(.SD)>20, .SDcols= patterns("input")]
mat <- mat[(check), !"check"]
mat[, rn:= paste0(enh_L, "_vs_", enh_R), .(enh_L, enh_R)]
cmb <- data.table(DSCP_rep1= c("DSCP_rep1", "DSCP_rep3"), input_rep1= c("input_rep2", "input_rep4"),
                  DSCP_rep2= c("DSCP_rep2", "DSCP_rep4"), input_rep2= c("input_rep3", "input_rep5"))

for(i in seq(nrow(cmb)))
{
  #----------------------------------------------#
  # 2- DESeq 2
  #----------------------------------------------#
  sampleTable <- data.table(rn= c(cmb[i, DSCP_rep1:input_rep2]))
  sampleTable[, condition:= tstrsplit(rn, "_", keep= 1)]
  sampleTable[, replicate:= c("rep1", "rep2"), condition]
  sampleTable <- data.frame(sampleTable, row.names = "rn")
  cols <- rownames(sampleTable)
  DF <- data.frame(mat[, ..cols], row.names = mat$rn)
  # DESeq2
  dds <- DESeqDataSetFromMatrix(countData= DF, colData= sampleTable, design= ~replicate+condition)
  # sizeFactors(dds) <- estimateSizeFactorsForMatrix(as.matrix(DF[grep("control.*vs.*control", rownames(DF)),]+0.25)) # Use pseudocounts
  controls <- apply(as.matrix(DF[grep("control.*vs.*control", rownames(DF)),]), 2, sum)
  sizeFactors(dds) <- controls/min(controls)
  res <- DESeq(dds)
  # Differential expression
  diff <- as.data.table(as.data.frame(lfcShrink(res, contrast= c("condition", "DSCP", "input"))), keep.rownames= T)
  diff[, c("enh_L", "enh_R"):= tstrsplit(rn, "_vs_")]
  diff <- diff[, .(enh_L, enh_R, baseMean, log2FoldChange, lfcSE, stat, pvalue, padj)]
  
  #----------------------------------------------#
  # 3- Robust controls from full exp
  #----------------------------------------------#
  pediff <- readRDS("/groups/stark/vloubiere/projects/pe_STARRSeq/Rdata/processed_peSTARRSeq_data/DESeq2_FC_table.rds")
  # Identify robust controls producing few outliers
  ctls <- pediff[grepl("control", enh_L) | grepl("control", enh_R), .(enh_L, enh_R, log2FoldChange)]
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
  clean <- unique(clean[, .(ID= paste0(enh_L, "_vs_", enh_R), enh_L, enh_R, baseMean= i.baseMean, 
                            log2FoldChange= i.log2FoldChange, padj= i.padj, log2FoldChange_rev,
                            median_L, median_R, log2FC_add, log2FC_mult= median_L+median_R, diff)])
  
  file <- paste0("/groups/stark/vloubiere/projects/pe_STARRSeq/Rdata/processed_peSTARRSeq_data/", 
                 gsub("_rep", "r", paste(rownames(sampleTable), collapse= "_")), ".rds")
  saveRDS(clean, file)
}


