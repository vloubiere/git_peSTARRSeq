setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")

short <- fread("db/final_tables_exp_model/vllib002_short_spacer_peSTARRSeq_max_cutoff_final_oe.txt")
long <- fread("db/final_tables_exp_model/vllib006_long_spacer_peSTARRSeq_high_cutoff_final_oe.txt")

mer_L <- merge(unique(short[, .(L, median_L)]), 
               unique(long[, .(L, median_L)]),
               by= "L",
               suffixes= c("_short", "_long"))
.lm_L <- lm(median_L_long~median_L_short, mer_L)
mer_R <- merge(unique(short[, .(R, median_R)]), 
               unique(long[, .(R, median_R)]),
               by= "R",
               suffixes= c("_short", "_long"))
.lm_R <- lm(median_R_long~median_R_short, mer_R)

dir.create("pdf/comparison_STARR_revSTARR", showWarnings = F)

pdf("pdf/comparison_STARR_revSTARR/comparison_STARR_revSTARR_LR.pdf", width = 7, height = 2.75)
par(las= 1,
    pch= 19)
layout(matrix(1:4, byrow = T, ncol= 4),
       widths = c(1,0.5,1,0.5))
lims <- c(-3, 8)
cex <- 0.5

# LEFT
plot(mer_L$median_L_short, 
     mer_L$median_L_long,
     xlim= lims,
     ylim= lims,
     cex= cex,
     xlab= "A act. short spacer",
     ylab= "A act. long spacer")
abline(.lm_L)
legend("topleft",
       paste0("PCC= ", round(cor.test(mer_L$median_L_short, mer_L$median_L_long)$estimate, 2)), 
       bty="n")
par(mar= c(5.1, 2, 4.1, 2.1))
boxplot(mer_L$median_L_short,
        mer_L$median_L_long, 
        notch= T,
        ylim= c(-2, 7),
        cex= 0.6, 
        names= c("short", "long"),
        las= 2)
segments(1, 6.5, 2, 6.5)
text(1.5, 
     6.5, 
     pos= 3,
     offset= 0.2, 
     cex= 0.6,
     paste0("Wilcox=", 
            formatC(wilcox.test(mer_L$median_L_short, 
                                mer_L$median_L_long)$p.value, format = "e", digits = 1)))


# RIGHT
par(mar= c(5.1, 4.1, 4.1, 2.1))
plot(mer_R$median_R_short, 
     mer_R$median_R_long,
     xlim= lims,
     ylim= lims,
     cex= cex,
     xlab= "B act. short spacer",
     ylab= "B act. long spacer")
abline(.lm_R)
legend("topleft",
       paste0("PCC= ", round(cor.test(mer_R$median_R_short, mer_R$median_R_long)$estimate, 2)), 
       bty="n")
par(mar= c(5.1, 2, 4.1, 2.1))
boxplot(mer_R$median_R_short,
        mer_R$median_R_long,
        notch= T,
        ylim= c(-4, 8),
        cex= 0.6, 
        names= c("short", "long"),
        las= 2)
segments(1, 7.5, 2, 7.5)
text(1.5, 
     7.5, 
     pos= 3,
     offset= 0.2, 
     cex= 0.6,
     paste0("Wilcox=", 
            formatC(wilcox.test(mer_R$median_R_short, 
                                mer_R$median_R_long)$p.value, format = "e", digits = 1)))
dev.off()