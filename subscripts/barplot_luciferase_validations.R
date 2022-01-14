setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")

# Import STARR-Seq data
dat <- readRDS("Rdata/validations_luciferase_final_table.rds")
dat[, col:= colorRampPalette(c(col_L, col_R))(3)[2], .(col_L, col_R)]
setkeyv(dat, c("group_L", "group_R"))
dat <- na.omit(dat[(group_L=="control" & group_R=="control") | (group_L!="control" & group_R!="control")])
dat <- dat[order(factor(col, levels= c("#D3D3D3", "#FF6347", "#B99260", "#74C27A")), additive_luc, log2FoldChange_luc)]

# barplot ---------------#
pdf("pdf/luciferase/Barplot_luciferase_validations.pdf", 9, 4)
par(mar= c(2, 4.1, 2.1, 2.1),
    las= 2)
bar <- barplot(dat[, additive_luc], 
               col= dat$col, 
               space = 3, 
               border = dat$col,
               ylim= c(0,8),
               ylab= "Normalized luciferase activity (log2)")[,1]
barplot(dat[, mean_luc_L], 
        col= "black",
        space= 3,
        add= T,
        density = 20, 
        angle = 45)
barplot(dat[, log2FoldChange_luc], 
        col= dat$col,
        space= c(4.25, rep(3, length(bar)-1)),
        add= T,
        border = "black")
segments(bar+1.25,
         dat$log2FoldChange_luc-dat$sd_luc,
         bar+1.25,
         dat$log2FoldChange_luc+dat$sd_luc)
legend("topleft",
       density = c(0,20,0), 
       angle = c(NA, 45, NA),
       legend= c("Left ind. act.", "Right ind. act.", "Combined"),
       bty= "n")
segments(0.5, 6.5, 0.5, 6.7, xpd= T)

dev.off()