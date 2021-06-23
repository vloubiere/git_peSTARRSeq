dat <- readRDS("Rdata/validations_luciferase_final_table.rds")

# Add features
feat <- readRDS("Rdata/master_lib_features.rds")
dat[feat, group_L:= i.group, on= "L==ID"]
dat[feat, group_R:= i.group, on= "R==ID"]

# Collapse per enhancer pair
coll <- dat[, .(mean= mean(luc_norm, na.rm= T), 
                sd= sd(luc_norm, na.rm= T),
                mean_add= mean(luc_mean_L+luc_mean_R, na.rm= T),
                sd_add= sd(luc_mean_L+luc_mean_R, na.rm= T)), .(L, R, group_L, group_R)]
coll <- coll[group_L=="dev" & group_R=="dev"]

# Add STARR-Seq data
STARR <- readRDS("Rdata/")
coll[STARR, diff:= i.diff, on= c("L==enh_L", "R==enh_R")]
# coll[STARR, diff:= i.log2FoldChange-i.add, on= c("L", "R")]

pdf("pdf/PCC_residuals_luc_validations_peSTARRSeq.pdf", 4.5, 5)
x <- log2(coll$mean)-log2(coll$mean_add)
y <- coll$diff
plot(x,
     y, 
     xlab= "luciferase residuals (log2 o/e)",
     ylab= "pe-STARR-Seq",
     pch= 19, 
     las= 1)
legend("bottomleft", 
       bty= "n", 
       legend = paste("PCC=", round(cor.test(x,y)$estimate, 2)))
dev.off()
