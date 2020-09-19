load("/groups/stark/almeida/data/motifs/motif_collection_v7_no_transfac_SteinAerts/TF_clusters_PWMs.RData")
setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
source("/groups/stark/vloubiere/projects/pe_STARRSeq/scripts/source_data.R")
source("/groups/stark/vloubiere/scripts/R_functions/my_plots.R")
require(ppcor)
require(gridExtra)

#### soi
cur <- feat[BA_enhancer_group=="dev"]
cur[, sgl_DSCP_gw_diff:= sgl_DSCP_gw_log2FC-sgl_DSCP_gw_exp_add]
med <- median(cur$sgl_DSCP_gw_diff, na.rm = T)
cur[between(twist_log2FC, 2, 6), sup_group:= ifelse(sgl_DSCP_gw_diff <= med, "low", "high")]
cur <- cur[between(sgl_DSCP_gw_diff, quantile(sgl_DSCP_gw_diff, 0.0005, na.rm= T), quantile(sgl_DSCP_gw_diff, 0.995, na.rm= T))]

#### PCC motif counts~residuals
c_mot <- mot[cur[!is.na(sup_group), .(uniq_ID, sup_group, sgl_DSCP_gw_diff)], , on= "uniq_ID"]
c_mot[, check:= sum(high_motif_count), motif]
c_mot <- c_mot[check>=50, !"check"]
c_mot[, high_motif_count_bin:= 
        {
          cur <- as.data.table(table(high_motif_count))
          cur[, check:= 1-cumsum(N)/sum(N)]
          if(any(cur$check<0.025))
          {
            sat <- as.numeric(cur[max(which(check>=0.025)), high_motif_count])
            high_motif_count[high_motif_count>=sat] <- sat
          }
          high_motif_count
        }, motif]
c_mot[ , PCC:= if(any(high_motif_count_bin>0)){cor.test(high_motif_count_bin, sgl_DSCP_gw_diff)$estimate}else{as.numeric(NA)}, motif]
c_mot[ , PCC_pval:= if(any(high_motif_count_bin>0)){cor.test(high_motif_count_bin, sgl_DSCP_gw_diff)$p.value}else{as.numeric(NA)}, motif]
#### PCC motif scores~residuals
c_mot[ , PCC_s:= if(any(low_motif_score>0)){cor.test(low_motif_score, sgl_DSCP_gw_diff)$estimate}else{as.numeric(NA)}, motif]
c_mot[ , PCC_s_pval:= if(any(low_motif_score>0)){cor.test(low_motif_score, sgl_DSCP_gw_diff)$p.value}else{as.numeric(NA)}, motif]
#### fisher test high~low SA
c_mot[, fisher_log2OR:= if(any(high_motif_count_bin>0)){log2(fisher.test(table(sup_group=="high", high_motif_count_bin>0))$estimate)}else{as.numeric(NA)}, motif]
c_mot[, fisher_pval:= if(any(high_motif_count_bin>0)){fisher.test(table(sup_group=="high", high_motif_count_bin>0))$p.value}else{as.numeric(NA)}, motif]
setkey(c_mot, motif)

res <- na.omit(unique(c_mot[, .(motif, PCC, PCC_pval, PCC_s, PCC_s_pval, fisher_log2OR, fisher_pval, Dmel_prot)]))
res[, c("PCC_padj", "PCC_s_padj", "fisher_padj") := .(p.adjust(PCC_pval), p.adjust(PCC_s_pval), p.adjust(fisher_pval))]

#### plots
pdf("pdf_wip/motifs_counts_analysis_lm_fisher_repSTARRSeq_1.0.pdf", width = 15, height = 11)
layout(matrix(c(rep(c(1,1,2,2,3,3,4,4), 2), 5:28), ncol= 8, byrow = T))
par(lty= 2, las= 1)

smoothScatter(cur[, .(twist_log2FC, obs_vs_exp= sgl_DSCP_gw_log2FC-sgl_DSCP_gw_exp_add)], ylim= c(-2, 2))
abline(v= c(2, 6))
abline(h= med)
text(c(4, 4), c(-1.8, 2), c("low", "high"), offset= 0)

res <- res[order(PCC)]
plot(res$PCC, pch= NA, xlim= c(-100, nrow(res)+100), ylab= "PCC (motif counts vs residuals)")
text(1:nrow(res), res$PCC, labels = res$Dmel_prot, cex= 0.3, col= ifelse(res$PCC_padj<0.001, "red", "black"))

res <- res[order(fisher_log2OR)]
plot(res$fisher_log2OR, pch= NA, xlim= c(-100, nrow(res)+100), ylab= "fisher log2OR (motif counts vs residuals)")
text(1:nrow(res), res$fisher_log2OR, labels = res$Dmel_prot, cex= 0.3, col= ifelse(res$fisher_padj<0.001, "red", "black"))

res[PCC_padj < 0.001, c("col", "ord"):= .("green", 2)]
res[fisher_padj < 0.001, c("col", "ord"):= .("gold", 2)]
res[PCC_padj < 0.001 & fisher_padj < 0.001, c("col", "ord"):= .("red", 1)]
res[is.na(col), c("col", "ord"):= .("black", 4)]
res <- res[order(ord, decreasing = T)]
plot(res$PCC, res$fisher_log2OR, pch= NA, las= 1, xlab= "PCC", ylab= "fisher log2 OR")
text(res$PCC, res$fisher_log2OR, res$Dmel_prot, offset = 0, col= res$col, cex= 0.5)
legend("topleft", col = c("green", "gold", "red"), pch= 19, legend = c("PCC padj < 0.001", "fisher padj < 0.001" , "both padj < 0.001"), bty= "n", cex= 0.6)

sel <- res[order(PCC)][c(c(1:10),.N-c(9:0)), motif]
c_mot[sel, 
      {
        boxplot(sgl_DSCP_gw_diff~high_motif_count_bin, main= Dmel_prot[1], lty= 1)
        abline(lm(sgl_DSCP_gw_diff~high_motif_count_bin), col= "red")
        legend("topleft", legend = paste0("R2=", formatC(summary(lm(sgl_DSCP_gw_diff~high_motif_count_bin))$r.squared, digits = 2)), bty= "n")
        print(motif)
      }, keyby= .(PCC, motif)]
dev.off()

pdf("pdf_wip/motifs_scores_analysis_lm_fisher_repSTARRSeq_1.0.pdf", width = 15, height = 11)
layout(matrix(c(rep(c(1,1,2,2,3,3,4,4), 2), 5:28), ncol= 8, byrow = T))
par(lty= 2, las= 1)

smoothScatter(cur[, .(twist_log2FC, obs_vs_exp= sgl_DSCP_gw_log2FC-sgl_DSCP_gw_exp_add)], ylim= c(-2, 2))
abline(v= c(2, 6))
abline(h= med)
text(c(4, 4), c(-1.8, 2), c("low", "high"), offset= 0)

res <- res[order(PCC_s)]
plot(res$PCC_s, pch= NA, xlim= c(-100, nrow(res)+100), ylab= "PCC (motif counts vs residuals)")
text(1:nrow(res), res$PCC_s, labels = res$Dmel_prot, cex= 0.3, col= ifelse(res$PCC_s_padj<0.001, "red", "black"))

plot.new()
plot.new()

sel <- res[order(PCC_s)][c(c(1:10),.N-c(9:0)), motif]
c_mot[sel, 
      {
        plot(sgl_DSCP_gw_diff~low_motif_score, main= Dmel_prot[1], lty= 1, cex= 0.5, pch= 19, col= adjustcolor("lightgrey", 0.5))
        abline(lm(sgl_DSCP_gw_diff~low_motif_score), col= "red")
        legend("topleft", legend = paste0("R2=", formatC(summary(lm(sgl_DSCP_gw_diff~low_motif_score))$r.squared, digits = 2)), bty= "n")
        print(motif)
      }, keyby= .(PCC, motif)]
dev.off()

######## Check sgl motif composition
sgl_mot <- mot[uniq_ID==SGL & motif %in% res[col=="red", motif], ]
sgl_mot[res, PCC_counts_residuals:= round(i.PCC, 2), on= "motif"]
setorder(sgl_mot, -PCC_counts_residuals)

pdf("pdf_wip/motifs_counts_sgl_1.0.pdf", width = 12, height = 6)
grid.table(sgl_mot[, .(motif, Dmel_prot, PCC_counts_residuals, low_motif_count, high_motif_count)])
dev.off()

######## Compare with sgl-containing peSTARR-Seq pairs
sub <- dat[act_group=="active~active"]
setkey(c_mot, motif)

pdf("pdf_wip/motifs_effect_compare_rep_peSTARRSeq_1.0.pdf", width = 14, height = 12)
par(las= 1, mfrow= c(4, 5))
sgl_mot[, 
        {
          idx <- c_mot[motif][high_motif_count>0, uniq_ID]
          cur <- copy(sub[grepl("_C_", enh_L) & grepl("_C_", enh_R)])
          cur[!enh_L %in% idx & !enh_R %in% idx, class:= paste0("-/- \n(", .N, ")")]
          cur[enh_L %in% idx & !enh_R %in% idx, class:= paste0("+/- \n(", .N, ")")]
          cur[!enh_L %in% idx & enh_R %in% idx, class:= paste0("-/+ \n(", .N, ")")]
          cur[enh_L %in% idx & enh_R %in% idx, class:= paste0("+/+ \n(", .N, ")")]
          title <- paste0(Dmel_prot, "\n(exp effect= ", PCC_counts_residuals, ")")
          box <- boxplot(diff~class, cur, plot= F)
          boxplot(diff~class, cur, main= title, notch= T, outline= F, ylim= c(min(box$stats), max(box$stats)+1))
          my_pval_plot(diff~class, cur, pairs = list(c(1,2), c(1,3), c(1,4)))
          abline(h= median(cur[grepl("^-/-", class), diff]), lty= 2, lwd= 0.5)
          print(".")
        }, .(motif, Dmel_prot)]
dev.off()
































