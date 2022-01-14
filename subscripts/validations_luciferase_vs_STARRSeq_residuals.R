setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
dat <- readRDS("Rdata/validations_luciferase_final_table.rds")

pdf("pdf/luciferase/PCC_residuals_luc_validations_peSTARRSeq.pdf", 4.5, 5)
x <- dat$log2FoldChange_luc-dat$additive_luc
y <- dat$log2FoldChange_STARR-dat$additive_STARR
plot(x,
     y, 
     xlab= "luciferase residuals (log2 o/e)",
     ylab= "pe-STARR-Seq",
     pch= 19, 
     las= 1)
legend("topleft", 
       bty= "n", 
       legend = paste("PCC=", round(cor.test(x,y)$estimate, 2)))
dev.off()
