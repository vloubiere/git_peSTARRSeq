setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
require(data.table)

# Import data ----
lib <- readRDS("db/FC_tables/DSCP_large_WT_DESeq2.rds")
luc <- readRDS("Rdata/validations_luciferase_final_table.rds")
dat <- merge(luc,
             lib, 
             by= c("L", "R"),
             suffixes= c("_luc", "_STARR"))
dat[, actClass:= fcase(grepl("^control", L) & grepl("^control", R), "Ctl./Ctl.",
                       grepl("^control", R), "Enh./Ctl.",
                       grepl("^control", L), "Ctl./Enh.",
                       default= "Enh./Enh.")]
dat <- dat[actClass %in% c("Enh./Ctl.", "Ctl./Enh.")]
dat[, actClass:= factor(actClass, c("Enh./Ctl.", "Ctl./Enh."))]
dat <- na.omit(dat)

# Plot
col <- c("royalblue2", "purple")
Cc <- adjustcolor(col, .6)

pdf("pdf/draft/review_luciferase_vs_STARRSeq_differences.pdf", 5, 3)
vl_par(mfrow= c(1, 2),
       cex= .1)
vl_boxplot(log2FoldChange_luc~actClass,
           dat,
           outline= T,
           ylim= c(-1, 7),
           col= Cc,
           tilt.names= T,
           ylab= "Normalized luciferase act. (log2)",
           cex= .5)
vl_boxplot(log2FoldChange_STARR~actClass,
           dat,
           outline= T,
           ylim= c(-1, 7),
           col= Cc,
           tilt.names= T,
           ylab= "STARR-Seq activity (log2)",
           cex= .5)
dev.off()