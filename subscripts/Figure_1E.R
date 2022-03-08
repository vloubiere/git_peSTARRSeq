setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")

# Import STARR-Seq data
dat <- readRDS("Rdata/validations_luciferase_final_table.rds")
dat[, inferred:= is.na(log2FoldChange)]
dat[is.na(log2FoldChange), c("log2FoldChange", "act_wilcox_L", "act_wilcox_R"):= .(0, 1, 1)]
dat[, col:= fcase(group_L=="control" & group_R=="control", "lightgrey",
                  act_wilcox_L>=0.001 & group_R=="control", "royalblue2",
                  act_wilcox_L<0.001 & group_R=="control", "royalblue2",
                  group_L=="control" & act_wilcox_R<0.001, "royalblue2",
                  act_wilcox_L<0.001 & act_wilcox_R<0.001, "#74C27A")]
.lm <- lm(log2FoldChange_luc~log2FoldChange,
          dat)

#----------------------------------------------#
# PLOT
#----------------------------------------------#
# PCC ---------------#
pdf("pdf/draft/Figure_1E.pdf", 4.5, 4.5)
par(las= 1,
    mar= c(5.1, 4.1, 2.1, 2.1))
plot(dat[, .(log2FoldChange, log2FoldChange_luc)],
     col= dat$col, 
     pch= 19,
     xlab= "pe-STARR-Seq activity (log2)",
     ylab= "Normalized luciferase activity (log2)")
segments(dat$log2FoldChange,
         dat$log2FoldChange_luc-dat$sd_luc,
         dat$log2FoldChange,
         dat$log2FoldChange_luc+dat$sd_luc,
         col= dat$col)
dat[(inferred), points(log2FoldChange, log2FoldChange_luc)]
abline(.lm, lty=2)
legend("topleft", 
       legend = paste("R²=", round(summary(.lm)$r.squared, 2)), 
       lty= 2, 
       bty= "n")
dev.off()

