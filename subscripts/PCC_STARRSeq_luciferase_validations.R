setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")

# Import STARR-Seq data
dat <- readRDS("Rdata/validations_luciferase_final_table.rds")
dat <- dat[!is.na(log2FoldChange_luc) & ! is.na(log2FoldChange_STARR)]
dat[, col:= colorRampPalette(c(col_L, col_R))(3)[2], .(col_L, col_R)]
setkeyv(dat, c("group_L", "group_R"))
dat[.("control", "control"), pch:= 15]
dat[group_L=="control" & is.na(pch), pch:= 16]
dat[group_R=="control" & is.na(pch), pch:= 17]
dat[is.na(pch), pch:= 19]
.lm <- lm(log2FoldChange_luc~log2FoldChange_STARR,
          dat)

#----------------------------------------------#
# PLOT
#----------------------------------------------#
# PCC ---------------#
pdf("pdf/luciferase/PCC_luciferase_validations_peSTARRSeq.pdf", 4.5, 5)
par(las= 1)
plot(dat[, .(log2FoldChange_STARR, log2FoldChange_luc)],
     col= dat$col, 
     pch= dat$pch,
     xlab= "pe-STARR-Seq activity (log2)",
     ylab= "Normalized luciferase activity (log2)")
segments(dat$log2FoldChange_STARR,
         dat$log2FoldChange_luc-dat$sd_luc,
         dat$log2FoldChange_STARR,
         dat$log2FoldChange_luc+dat$sd_luc,
         col= dat$col)
abline(.lm, lty=2)
legend("bottomright", 
       legend = dat[, paste0(group_L, "*", group_R), .(group_L, group_R)]$V1, 
       col= unique(dat[, .(group_L, group_R, col)])$col, 
       bty= "n", 
       pch= unique(dat[, .(group_L, group_R, pch)])$pch,
       cex= 0.8)
legend("topleft", 
       legend = paste("R²=", round(summary(.lm)$r.squared, 2)), 
       lty= 2, 
       bty= "n")
dev.off()