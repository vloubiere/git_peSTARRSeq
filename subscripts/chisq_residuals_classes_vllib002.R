setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
require(vlfunctions)

# Import data
dat <- readRDS("db/FC_tables/vllib002_pe-STARR-Seq_DSCP_T8_SCR1_300_counts_norm_final_oe.rds")[!grepl("^control", L) & !grepl("^control", R)]
.m <- merge(unique(dat[, .(ID= L, resGroupL= factor(resGroupL, rev(levels(resGroupL))))]),
            unique(dat[, .(ID= R, resGroupR)]))
chisq <- chisq.test(.m$resGroupL, .m$resGroupR)
chisq_residuals <- matrix(chisq$residuals, nrow= nrow(chisq$residuals))
rownames(chisq_residuals) <- rownames(chisq$residuals)
colnames(chisq_residuals) <- colnames(chisq$residuals)
chisq_obs <- matrix(chisq$observed, nrow= nrow(chisq$observed))

pdf("pdf/draft/density_residuals_vllib002.pdf", 4, 3.6)
par(mar= c(6.5,6.5,2,4),
    las= 2,
    tcl= -0.2)
vl_heatmap(chisq_residuals,
           cluster_rows = F, 
           cluster_cols = F, 
           legend_title = "X2\nresiduals", 
           display_numbers = T,
           display_numbers_matrix = chisq_obs,
           breaks= c(-9,0,9))
title(main= paste0("X2 pval= ", formatC(chisq$p.value, format = "e", digits = 1)),
      line= 1)
dev.off()