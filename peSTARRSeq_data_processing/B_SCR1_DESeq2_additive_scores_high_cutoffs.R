setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
sapply(list.files("/groups/stark/vloubiere/functions/", ".R$", full.names = T), source)
require(data.table)
require(DESeq2)
require(circlize)

#----------------------------------------------#
# 1- Load counts
#----------------------------------------------#
if(!file.exists("db/read_counts/all_uniq_counts.rds"))
{
  dat <- data.table(file= list.files("db/read_counts/", "SCR1.*uniq.UMI.rds", full.names= T))
  dat[, sample:= gsub("libvl002_SCR1_|.uniq.UMI.rds", "", basename(file))]
  dat <- dat[, readRDS(file), .(file, sample)]
  dat <- dat[, .(counts= .N), c(colnames(dat))]
  saveRDS(dat, "db/read_counts/all_uniq_counts.rds")
}
if(!file.exists("Rdata/processed_peSTARRSeq_data/filtered_counts_prior_DESeq2.rds"))
{
  dat <- readRDS("db/read_counts/all_uniq_counts.rds")
  # Merge rep 1 and 2
  mer <- data.table(old= paste0("input_rep", 1:6), new= paste0("input_rep", c(1,1,2,3,4,5)))
  dat[mer, sample:= i.new, on= "sample==old"]
  dat <- dat[, .(counts= sum(counts)), .(sample, enh_L, enh_R)]
  mat <- dcast(dat, enh_L+enh_R~sample, value.var = "counts", fill = 0)
  # low counts cutoff
  check <- mat[, apply(.SD, 1, function(x) all(x>0) & sum(x)>100), .SDcols= patterns("input")]
  mat <- mat[(check)]
  check <- mat[, apply(.SD, 1, function(x) all(x>0)), .SDcols= patterns("DSCP")]
  mat <- mat[(check)]
  mat[, rn:= paste0(enh_L, "_vs_", enh_R), .(enh_L, enh_R)]
  saveRDS(mat, "Rdata/processed_peSTARRSeq_data/filtered_counts_prior_DESeq2_high_cutoff.rds")
}


#----------------------------------------------#
# 2- DESeq 2
#----------------------------------------------#
# Format tables
if(!file.exists("/groups/stark/vloubiere/projects/pe_STARRSeq/Rdata/processed_peSTARRSeq_data/dds_result_object_high_cutoff.rds"))
{
  mat <- readRDS("Rdata/processed_peSTARRSeq_data/filtered_counts_prior_DESeq2_high_cutoff.rds")
  sampleTable <- grep("rep", colnames(mat), value = T)
  sampleTable <- data.frame(condition= sapply(sampleTable, function(x) strsplit(x, "_")[[1]][1]),
                            replicate= sapply(sampleTable, function(x) strsplit(x, "_")[[1]][2]),
                            row.names= sampleTable)
  DF <- data.frame(mat[, DSCP_rep1:input_rep5], row.names = mat$rn)
  # DESeq2
  dds <- DESeqDataSetFromMatrix(countData= DF, colData= sampleTable, design= ~replicate+condition)
  # SizeFactors
  sizeFactors(dds) <- estimateSizeFactorsForMatrix(as.matrix(DF[grep("control.*vs.*control", rownames(DF)),]))
  # Result
  res <- DESeq(dds)
  saveRDS(res, "/groups/stark/vloubiere/projects/pe_STARRSeq/Rdata/processed_peSTARRSeq_data/dds_result_object_high_cutoff.rds")
  
  # Differential expression
  diff <- as.data.table(as.data.frame(results(res, contrast= c("condition", "DSCP", "input"))), keep.rownames= T)
  diff[, c("enh_L", "enh_R"):= tstrsplit(rn, "_vs_")]
  diff <- diff[, .(enh_L, enh_R, baseMean, log2FoldChange, lfcSE, stat, pvalue, padj)]
  
  boxplot(diff[grepl("control", enh_L) & grepl("control", enh_R), log2FoldChange], notch= T)
  abline(h= 0, lty= 2)
  saveRDS(diff, "/groups/stark/vloubiere/projects/pe_STARRSeq/Rdata/processed_peSTARRSeq_data/DESeq2_FC_table_high_cutoff.rds")
}

#----------------------------------------------#
# 3- Robust controls
#----------------------------------------------#
# Identify robust controls with high ctl~X / X~ctl PCCs
diff <- readRDS("/groups/stark/vloubiere/projects/pe_STARRSeq/Rdata/processed_peSTARRSeq_data/DESeq2_FC_table_high_cutoff.rds")
PCC <- data.table(ctl= unique(grep("^control", c(diff$enh_L, diff$enh_R), value = T)))
PCC[, PCC:=
    {
      sub <- diff[enh_L==ctl]
      sub[diff, rev:= i.log2FoldChange, on=c("enh_L==enh_R", "enh_R==enh_L")]
      sub <- na.omit(sub[, .(log2FoldChange, rev)])
      if(nrow(sub)>2)
      {
        cor.test(sub[[1]], sub[[2]])$estimate
      }else
      {
        as.numeric(NA)
      }
    }, ctl]
PCC <- na.omit(PCC)

# Identify robust controls producing few outliers
ctls <- diff[enh_L %in% PCC$ctl | enh_R %in% PCC$ctl, .(enh_L, enh_R, log2FoldChange)]
ctls[enh_L %in% PCC$ctl & !is.na(log2FoldChange), zscore_L:= scale(log2FoldChange), enh_R]
ctls[enh_R %in% PCC$ctl & !is.na(log2FoldChange), zscore_R:= scale(log2FoldChange), enh_L]

lim <- c(-0.75, 0.75)
b1 <- my_boxplot(zscore_L~enh_L, ctls[enh_L %in% PCC$ctl], plot= F)
Cc1 <- as.character(cut(b1$DT_plot$lim3, c(-Inf, lim, Inf), labels = c("cornflowerblue", "white", "red")))
b2 <- my_boxplot(zscore_R~enh_R, ctls[enh_R %in% PCC$ctl], plot= F)
Cc2 <- as.character(cut(b2$DT_plot$lim3, c(-Inf, lim, Inf), labels = c("cornflowerblue", "white", "red")))
Ccgood <- ifelse(between(b1$DT_plot$lim3, lim[1], lim[2]) & between(b2$DT_plot$lim3, lim[1], lim[2]) & PCC$PCC>0.7, "black", "red")

pdf("pdf/peSTARRSeq/diag_plots_selection_negative_controls_high_cutoff.pdf", 10, 17)
layout(matrix(1:3, ncol = 3), widths = c(0.75,0.5,0.5))
par(las= 1, mar= c(5,17,1,1), xaxs= "i", yaxs= "i")

my_boxplot(zscore_L~enh_L, ctls[enh_L %in% PCC$ctl], col= Cc1, points.jit = 0.2, horizontal = T, names= NA)
names <- ctls[enh_L %in% PCC$ctl, enh_L, enh_L]$enh_L
text(par("usr")[1], seq(names), names, pos= 2, xpd= T, offset= 0.7, col= Ccgood)
abline(v= lim)
my_fig_label("A", cex= 2)

par(mar= c(5,5,1,1))
my_boxplot(zscore_R~enh_R, ctls[enh_R %in% PCC$ctl], col= Cc2, names= NA, points.jit = 0.2, 
           horizontal = T, yaxt= "n")
abline(v= lim)
my_fig_label("B", cex= 2)

barplot(PCC$PCC, col= ifelse(PCC$PCC<0.7, "cornflowerblue", "grey"), xlab= "PCC ctl~X / X~ctl", 
        horiz = T, border=NA, xlim= c(0, 1))
abline(v= 0, col= "grey")
abline(v= 0.7)
my_fig_label("C", cex= 2)
dev.off()

# Filtering
ctls <- b1$DT_plot[abs(lim3)<lim[2], .id]
ctls <- ctls[ctls %in% b2$DT_plot[abs(lim3)<lim[2], .id]]
ctls <- ctls[ctls %in% PCC[PCC>0.7, ctl]]

#----------------------------------------------#
# 4- Additive scores
#----------------------------------------------#
# Compute expected scores
exp <- copy(diff)
# Individual activity
exp[, c("all_L", "median_L", "sd_L", "N_L"):= .SD[enh_R %in% ctls, .(list(log2FoldChange), median(log2FoldChange), sd(log2FoldChange), .N)], enh_L]
exp[, c("all_R", "median_R", "sd_R", "N_R"):= .SD[enh_L %in% ctls, .(list(log2FoldChange), median(log2FoldChange), sd(log2FoldChange), .N)], enh_R]
saveRDS(exp, "Rdata/processed_peSTARRSeq_data/all_expected_score_high_cutoff.rds")
exp <- na.omit(exp[!grepl("control", enh_L) & !grepl("control", enh_R) & N_L>4 & N_R>4, !c("N_L", "N_R")])
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
# bg_act <- median(diff[enh_L %in% ctls & enh_R %in% ctls, log2FoldChange])
# exp[, log2FC_add:= log2((2^median_L)+(2^median_R)-2^bg_act), .(enh_L, enh_R)]
exp[, log2FC_add:= log2((2^median_L)+(2^median_R)), .(enh_L, enh_R)]
exp[, diff:= log2FoldChange-log2FC_add]
# Add extra stats
exp[exp, log2FoldChange_rev:= i.log2FoldChange, on= c("enh_L==enh_R", "enh_R==enh_L")]
# SAVE short table
clean <- exp[diff, , on=c("enh_L", "enh_R")]
clean <- unique(clean[, .(ID= paste0(enh_L, "_vs_", enh_R), enh_L, enh_R, baseMean= i.baseMean, log2FoldChange= i.log2FoldChange, padj= i.padj, log2FoldChange_rev,
                          median_L, sd_L, median_R, sd_R, log2FC_add, log2FC_mult= median_L+median_R, diff)])
clean <- clean[!is.na(median_L) & !is.na(median_R) & !is.na(log2FC_add) & !is.na(log2FoldChange)]
saveRDS(clean, "/groups/stark/vloubiere/projects/pe_STARRSeq/Rdata/processed_peSTARRSeq_data/SCR1_peSTARRSeq_final_table_high_cutoff.rds")

