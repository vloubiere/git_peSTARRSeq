setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
sapply(list.files("/groups/stark/vloubiere/functions/", ".R$", full.names = T), source)
require(data.table)
options(datatable.print.topn=1)
require(yarrr)
require(ggplot2)
require(patchwork)

feat <- readRDS("Rdata/library/lib_features.rds")
.c <- readRDS("Rdata/processed_peSTARRSeq_data/filtered_counts_prior_DESeq2.rds")
dat <- readRDS("Rdata/processed_peSTARRSeq_data/expected_score.rds")

#----------------------------------------------------------------#
# 1- format data
#----------------------------------------------------------------#
# Add group
dat[feat, group_L:= i.group, on= "enh_L==ID"]
dat[feat, group_R:= i.group, on= "enh_R==ID"]
# Handle E coli cases
dat[is.na(group_L), group_L:="control"]
dat[is.na(group_R), group_R:="control"]

pld <- dat[group_L==group_R]
pld[, group:= group_L]
pld[, homotypic:= ifelse(enh_L==enh_R, "homotypic", "heterotypic")]
pld <- melt(pld, id.vars = "group", measure.vars = list("homotypic", "log2FoldChange"))

# Add group
.c[feat, group_L:= i.group, on= "enh_L==ID"]
.c[feat, group_R:= i.group, on= "enh_R==ID"]
# Handle E coli cases
.c[is.na(group_L), group_L:="control"]
.c[is.na(group_R), group_R:="control"]
.c[, input:= log2(rowSums(.SD)+1), .SDcols= patterns("input")]
.c[, DSCP:= log2(rowSums(.SD)+1), .SDcols= patterns("DSCP")]

plc <- .c[group_L==group_R]
plc[, group:= group_L]
plc[, homotypic:= ifelse(enh_L==enh_R, "homotypic", "heterotypic")]
plc <- melt(plc, id.vars = "group", measure.vars = list("homotypic", "input", "DSCP"))

#-----------------------------------------------------#
# PLOT
#-----------------------------------------------------#
pdf("pdf/peSTARRSeq/homotypic_pairs_checks.pdf", width = 5, height = 10)
par(mfrow=c(3,1), mar= c(10,5,2,2))
my_boxplot(value2~value1+group, plc, las= 2, at= c(1,2,9,10,7,8,5,6,3,4), ylab= "log2(counts+1)",
           col_box= rep(c("lightgrey", "royalblue2", "gold", "tomato", "#74C27A"), each= 2),
           pval_list = list(c(1,2),c(3,4),c(5,6),c(7,8),c(9,10)), ylim= c(3.5, 10), main= "input counts")

my_boxplot(value3~value1+group, plc, las= 2, at= c(1,2,9,10,7,8,5,6,3,4), ylab= "log2(counts+1)",
           col_box= rep(c("lightgrey", "royalblue2", "gold", "tomato", "#74C27A"), each= 2),
           pval_list = list(c(1,2),c(3,4),c(5,6),c(7,8),c(9,10)), ylim= c(2, 12), main= "DSCP counts")

my_boxplot(value2~value1+group, pld, las= 2, at= c(1,2,9,10,7,8,5,6,3,4), ylab= "activity (log2)",
           col_box= rep(c("lightgrey", "royalblue2", "gold", "tomato", "#74C27A"), each= 2),
           pval_list = list(c(1,2),c(3,4),c(5,6),c(7,8),c(9,10)), ylim= c(-4,9.5), main= "activity")
abline(h= 0, lty= 2)
# legend("topleft", fill= Cc1, legend= c("control", "OSC", "inducible", "hk", "dev"), bty= "n")
dev.off()


