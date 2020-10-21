setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
sapply(list.files("/groups/stark/vloubiere/functions/", ".R$", full.names = T), source)
require(data.table)

high <- readRDS("Rdata/processed_peSTARRSeq_data/SCR1_peSTARRSeq_final_table_high_cutoff.rds")
dat <- readRDS("Rdata/processed_peSTARRSeq_data/SCR1_peSTARRSeq_final_table.rds")
dat[high, c("median_Lh", "sd_Lh", "median_Rh", "sd_Rh"):= .(i.median_L, i.sd_L, i.median_R, i.sd_R), on= c("enh_L", "enh_R")]
feat <- readRDS("Rdata/library/lib_features.rds")
luc <- readRDS("Rdata/luciferase_validations/C_luc_validations_final_table.rds")
luc[, luc_L:= mean(log2(.SD[grepl("^control", enh_R), luc_norm]), na.rm= T), enh_L]
luc[, luc_R:= mean(log2(.SD[grepl("^control", enh_L), luc_norm]), na.rm= T), enh_R]

# format
.c <- merge(unique(dat[, .(ID= enh_L, median_L, sd_L, median_Lh, sd_Lh)]),
            unique(dat[, .(ID= enh_R, median_R, sd_R, median_Rh, sd_Rh)]), all.x= T, all.y= T)
.c <- merge(.c, feat[, .(ID, TWIST_L= dev_log2FoldChange,
                             TWIST_Lh= dev_log2FoldChange, 
                             TWIST_R= dev_log2FoldChange,
                             TWIST_Rh= dev_log2FoldChange)], all.x= T)
.c <- merge(.c, unique(luc[, .(ID= enh_L, luc_L, luc_Lh= luc_L)]), all.x= T)
.c <- merge(.c, unique(luc[, .(ID= enh_R, luc_R, luc_Rh= luc_R)]), all.x= T)
pl <- melt(.c, variable.name = "pos",
           measure.vars = patterns(median= "^median", sd= "^sd", TWIST= "TWIST", luc= "luc"))
pl[, pos:= gsub("^median_", "", grep("^median", colnames(.c), value = T))[pos]]
pl <- melt(pl, measure.vars = c("TWIST", "luc"))

setkeyv(pl, c("pos", "variable"))
ord <- data.table(c("L", "L", "Lh", "R", "R", "Rh"), c("TWIST", "luc", "TWIST", "TWIST", "luc", "TWIST")) 
pl <- pl[ord]
pl[ord, lab:= c("A", "C", "E", "B", "D", "F")[.GRP], .(variable, pos)]

# PLOT
pdf("pdf/peSTARRSeq/individual_scores_TWIST_compare.pdf", width= 5, height = 7.3)
par(mfcol= c(3, 2), las= 1, mar= c(5,5,2,2))
col <- adjustcolor("lightgrey", 0.5)
pl[ord, 
   {
     .c <- na.omit(.SD)
     xl <- c("TWIST-STARR-seq activity", "luciferase activity")[match(variable[1], c("TWIST", "luc"))]
     plot(rep(.c$value, 2), c(.c$median-.c$sd, .c$median+.c$sd), pch= NA, ylab= "pe-STARR-seq activity", xlab= xl)
     my_fig_label(lab, cex= 2)
     title <- c("X~control", "X~control (high cutoff)", "control~X", "control~X (high cutoff)")[match(pos[1], c("L", "Lh", "R", "Rh"))]
     mtext(title, line= 0.3, cex= 0.8)
     segments(.c$value, .c$median-.c$sd, .c$value, .c$median+.c$sd, col= col)
     points(.c$value, .c$median, cex= 0.7, pch= 1, lwd= 0.5)
     lines(seq(-5, 15, 0.1), predict(loess(.c$median~.c$value), seq(-5, 15, 0.1)), col= "red")
     abline(0, 1, lty= 2)
     print(nrow(.c))
   }, .(variable, pos, lab)]
dev.off()
