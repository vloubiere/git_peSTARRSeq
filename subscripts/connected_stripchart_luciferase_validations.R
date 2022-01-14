setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")

# Import STARR-Seq data
dat <- readRDS("Rdata/validations_luciferase_final_table.rds")
dat[, col:= colorRampPalette(c(col_L, col_R))(3)[2], .(col_L, col_R)]
dat <- dat[(group_L=="control" & group_R=="control") |
           (group_L!="control" & group_R!="control"),]
dat[group_L=="control", x:=1 ]
dat[group_L=="dev" & group_R=="dev", x:= 5]
dat[is.na(x), x:= 3]

#--------------------------------#
# PLOT
#--------------------------------#
pdf("pdf/luciferase/Stripchart_luciferase_validations_peSTARRSeq.pdf", 4, 5)
par(mar= c(7.1, 4.1, 2.1, 2.1),
    las= 2)
plot(dat$x+1, 
     2^dat$log2FoldChange_luc,
     col= dat$col, 
     pch= 19,
     xlim= c(0.5,6.5),
     xaxt= "n",
     ylab= "Normalized luciferase activity",
     xlab= NA)
points(dat$x, 
       2^dat$additive_luc,
       col= dat$col, 
       pch= 19,
       xlim= c(0.5,4.5))
segments(dat$x+1, 
         2^dat$log2FoldChange_luc,
         dat$x, 
         2^dat$additive_luc,
         col= dat$col)
axis(1, 
     at = 1:6, 
     rep(c("Exp. Additive", "Observed"), 3))
legend("topleft", 
       legend = dat[, paste0(group_L, "*", group_R), .(group_L, group_R)]$V1, 
       col= unique(dat[, .(group_L, group_R, col)])$col, 
       bty= "n", 
       pch= 19,
       cex= 0.8)
dev.off()
