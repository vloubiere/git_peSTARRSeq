setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
source("/groups/stark/vloubiere/projects/pe_STARRSeq/scripts/source_data.R")

### PLOTS ############################################################################################################################################
########### Flow chart ###########
# grViz("digraph flowchart {
#       # node definitions with substituted label text
#       node [fontname = Helvetica, shape = rectangle]
#       tab1 [label = '@@1']
#       tab2 [label = '@@2']
#       tab3 [label = '@@3']
#       tab4 [label = '@@4']
#       tab5 [label = '@@5']
#       tab6 [label = '@@6']
#       tab7 [label = '@@7']
#       tab8 [label = '@@8']
#       tab9 [label = '@@9']
#       tab10 [label = '@@10']
#       tab11 [label = '@@11']
#         # edge definitions with the node IDs
#       tab1 -> tab3
#       tab2 -> tab3 -> tab4 -> tab5 -> tab6 -> tab7
#       tab7 -> tab8 -> tab11
#       tab7 -> tab9 -> tab11
#       tab7 -> tab10 -> tab11
#       }
#         [1]: 'trim reads1 - trim_galore v0.6.2 --hardtrim3 34'
#       [2]: 'Build custom bowtie genome index \\n 1000 TWIST sequences + \\n 3 template switching constructs'
#       [3]: 'Align each read separately - bowtie2 v2.3.4.2 -U, default parameters'
#       [4]: 'Filter unique reads with mapq>10 - samtools v1.9'
#       [5]: 'Pair reads back together and collapse UMIs \\n -> counts/combination'
#       [6]: 'Compute expected additive counts \\n L~R exp.= median(Lenh~CTLs) + median(CTLs~Renh)'
#       [7]: 'Compute activity \\n DESeq2 \\n L~R combinations >20 reads \\n center on CTL~CTL combinations'
#       [8]: 'L~R activity \\n STARR-Seq obs./input obs.'
#       [9]: 'L~R expected add. activity \\n STARR-Seq exp./input obs.'
#       [10]: 'obs./exp. score \\n STARR-Seq exp./STARR-Seq obs.'
#       [11]: 'Further analyses \\n ...'
#       ")

########## Reads Stats ###########
# Compute sam stats
stats <- data.table(file= list.files("/groups/stark/vloubiere/data/pe_STARRSeq/Rdata", "all.rds|UMI.rds", full.names = T))
stats[, counts:= 
        {
          current <- readRDS(file)
          sum(current[, ncol(current), with= F])
        }, file]
stats[, c("sample", "stat"):= .(gsub("_UMI.rds|_all.rds", "", basename(file)), ifelse(grepl("all.rds$", file), "raw_counts", "UMI_collapsed"))]
stats <- dcast(stats, sample~stat, value.var= "counts")
stats[, "collapsing %":= round((1- UMI_collapsed/raw_counts)*100, 1)]
cols <- c("UMI_collapsed", "raw_counts")
stats[, (cols):= lapply(.SD, function(x) formatC(x, big.mark = ",")), .SDcols= cols]

# plot
pdf("/groups/stark/vloubiere/projects/pe_STARRSeq/pdf_sup/Read_stats_screen1.pdf", width = 17)
grid.table(stats)
dev.off()

########### Template switching ###########
cur <- copy(counts_raw)
cols <- grep("input_rep", colnames(cur), value = T)
cur[, input:= rowSums(.SD, na.rm = T), .SDcols= cols]
cols <- grep("DSCP_rep", colnames(cur), value = T)
cur[, DSCP:= rowSums(.SD, na.rm = T), .SDcols= cols]
cur <- cur[, .(enh_L, enh_R, input, DSCP)]
# Compute
ts_dat <- ts[, .(ID_vl, name= detail)]
ts_dat <- ts_dat[, .(input_reads= formatC(sum(cur[enh_L==ID_vl | enh_R==ID_vl, input]), format= "d", big.mark = ","),
                     "switch (%)"= round(100-(cur[enh_L==ID_vl & enh_R==ID_vl, input]/sum(cur[enh_L==ID_vl | enh_R==ID_vl, input])*100), 2),
                     DSCP_reads= formatC(sum(cur[enh_L==ID_vl | enh_R==ID_vl, DSCP]), format= "d", big.mark = ","),
                     "switch (%)"= round(100-(cur[enh_L==ID_vl & enh_R==ID_vl, DSCP]/sum(cur[enh_L==ID_vl | enh_R==ID_vl, DSCP])*100), 2)), name]
# plot
pdf("/groups/stark/vloubiere/projects/pe_STARRSeq/pdf_sup/Template_switching_check.pdf", 8, 4)
grid.table(ts_dat)
dev.off()

########### PCC replicates ###########
cols <- grep("rep", colnames(counts_raw))
mat <- as.matrix(counts_raw[!grepl("temp_switch", enh_L) & !grepl("temp_switch", enh_R), ..cols])
mat <- mat[complete.cases(mat),]
PCC <- cor(mat)

pdf("pdf_sup/PCC_replicates.pdf", 7, 7)
par(mar=c(6,6,5,5))
my_pheatmap(PCC, grid_lwd = 1, display_numbers= T, ROUND_FUN = function(x) round(x, 2), legend_title = "PCC")
dev.off()

########### Library complexity ###########
pl <- counts_raw[!grepl("temp_switch", enh_L) & !grepl("temp_switch", enh_R)]

pdf("pdf_sup/complexity_seq_SCR1.pdf", height = 4)
par(mfrow= c(1, 2))

my_empty_plot(xlim=c(-1, 10), ylim= c(0, 0.5), ylab= "Density", xlab= "log2(counts+1)", main= "library complexity")
lines(density(log2(pl$input_merge+1)), lwd= 3, col= adjustcolor("cornflowerblue", 0.5))
lines(density(log2(pl$DSCP_merge+1)), lwd= 3, col= adjustcolor("tomato", 0.5))
legend("topright", lwd=2, legend= c("Input", "STARR-Seq"), col= adjustcolor(c("cornflowerblue", "tomato"), 0.5), bty= "n")

input <- sort(log2(pl$input_merge+1))
input <- c(rep(0, 1e6-length(input)), input)
DSCP <- sort(log2(pl$DSCP_merge+1))
DSCP <- c(rep(0, 1e6-length(DSCP)), DSCP)
my_empty_plot(xlim=c(1, 1e6), ylim= c(0, 15), ylab= "log2(counts+1)", xlab= "combinations", main= "library complexity")
lines(input, lwd= 3, col= adjustcolor("cornflowerblue", 0.5))
lines(DSCP, lwd= 3, col= adjustcolor("tomato", 0.5))
legend("topleft", lwd=2, legend= c("Input", "STARR-Seq"), col= adjustcolor(c("cornflowerblue", "tomato"), 0.5), bty= "n")
dev.off()

########### Library biases sub libraries #############
pl <- counts_clean[!grepl("temp_switch", enh_L) & !grepl("temp_switch", enh_R)]
pl[lib, sub_L:= linker_ID, on= c("enh_L==ID_vl")]
pl[lib, sub_R:= linker_ID, on= c("enh_R==ID_vl")]
pl <- pl[, .(obs= sum(input_merge)), .(sub_L, sub_R)]
pl <- pl[CJ(sub_L= lib$linker_ID, sub_R= lib$linker_ID)[, .(exp= .N), .(sub_L, sub_R)], , on= c("sub_L", "sub_R")]
cols <- c("obs", "exp")
pl[, (cols):= lapply(.SD, function(x) x/sum(x)), .SDcols= cols]
pl[, ratio := log2(obs)-log2(exp)]
mat <- as.matrix(dcast(pl, sub_L~sub_R, value.var = "ratio"), 1)

pdf("pdf/Sub_lib_complex_o_vs_e.pdf", 4.5, 4.5)
par(mar=c(6,6,6,6))
my_pheatmap(mat, cluster_rows= F, cluster_cols= F, lim= c(-1, 1), display_numbers = T, grid_lwd = 1,
            main= "Sub-libraries", legend_title = "log2 obs/exp fraction (input)")
dev.off()

########### Library biases functional subgroups#############
pl <- counts_raw[!grepl("temp_switch", enh_L) & !grepl("temp_switch", enh_R)]
pl[lib, sub_L:= paste(group, detail), on= c("enh_L==ID_vl")]
pl[lib, sub_R:= paste(group, detail), on= c("enh_R==ID_vl")]
pl <- pl[, .(obs= sum(input_merge)), .(sub_L, sub_R)]
pl <- pl[CJ(sub_L= paste(lib$group, lib$detail), sub_R= paste(lib$group, lib$detail))[, .(exp= .N), .(sub_L, sub_R)], , on= c("sub_L", "sub_R")]
cols <- c("obs", "exp")
pl[, (cols):= lapply(.SD, function(x) x/sum(x, na.rm = T)), .SDcols= cols]
pl[, ratio:= log2(obs)-log2(exp)]
mat <- as.matrix(dcast(pl, sub_L~sub_R, value.var = "ratio"), 1)

pdf("pdf_sup/Func_categories_complex_o_vs_e.pdf", 10, 9)
par(mar=c(12,12,6,12))
my_pheatmap(mat, cluster_rows= F, cluster_cols= F, lim= c(-1, 1), display_numbers = T, grid_lwd = 1,
            main= "Functional categories o/e", legend_title = "log2 obs/exp fraction (input)")
dev.off()

########### Library combinations counts #############
pl <- counts_clean[!grepl("temp_switch", enh_L) & !grepl("temp_switch", enh_R)]
pl[lib, sub_L:= group, on= c("enh_L==ID_vl")]
pl[lib, sub_R:= group, on= c("enh_R==ID_vl")]
pl <- pl[, .(obs= .N), .(sub_L, sub_R)]
mat <- as.matrix(dcast(pl, sub_L~sub_R, value.var = "obs"), 1)

pdf("pdf_sup/Func_categories_N_cmbn.pdf", 8, 8)
par(mar= rep(10, 4))
my_pheatmap(mat, display_numbers = T, grid_lwd = 1, ROUND_FUN = function(x) formatC(x, big.mark = ","),
            main= "Functional categories combinations (DESeq2)", legend_title = "log2 obs/exp fraction")
dev.off()

########### Comparison BA screen ################
pl <- merge(unique(dat[!is.na(median_L), .(enh= enh_L, median_L)]), unique(dat[!is.na(median_R), .(enh= enh_R, median_R)]), all.x= T, all.y= T)
pl <- feat[pl, .(uniq_ID, twist_log2FC, median_L= i.median_L, median_R= i.median_R), on= "uniq_ID==enh"]

# Correlations
pdf("pdf_sup/TWIST_screen_BA_correlation.pdf", 4.5, 5)

plot(pl[, .(twist_log2FC, median_L)], xlab= "twist STARR-Seq", ylab= "pe-STARR-Seq Left individual activity",
     pch= 19, cex= 0.8, col= adjustcolor("grey", 0.6), las= 1)
.lm <- lm(median_L~twist_log2FC, pl)
abline(.lm, lty= 2)
rsq <- round(summary(.lm)$r.square, 2)
PCC <- round(cor.test(pl$twist_log2FC, pl$median_L)$estimate, 2)
legend("bottomright", legend = paste0("R2= ", rsq, "\nPCC= ", PCC), bty= "n")

plot(pl[, .(twist_log2FC, median_R)], xlab= "twist STARR-Seq", ylab= "pe-STARR-Seq Right individual activity",
     pch= 19, cex= 0.8, col= adjustcolor("grey", 0.6), las= 1)
.lm <- lm(median_R~twist_log2FC, pl)
abline(.lm, lty= 2)
rsq <- round(summary(.lm)$r.square, 2)
PCC <- round(cor.test(pl$twist_log2FC, pl$median_R)$estimate, 2)
legend("bottomright", legend = paste0("R2= ", rsq, "\nPCC= ", PCC), bty= "n")
dev.off()

# Groups activity
pl <- copy(dat[!grepl("temp_switch", enh_L) & !grepl("temp_switch", enh_R)])
pl[lib, meta_L:= paste(group, detail), on= "enh_L==ID_vl"]
pl[grepl("^hk", enh_L), meta_L:= "hk"]
pl[grepl("^control", enh_L), meta_L:= "control"]
pl[grepl("^OSC|^ecd|^heat", enh_L), meta_L:= "non-S2"]
pl[lib, meta_R:= paste(group, detail), on= "enh_R==ID_vl"]
pl[grepl("^hk", enh_R), meta_R:= "hk"]
pl[grepl("^control", enh_R), meta_R:= "control"]
pl[grepl("^OSC|^ecd|^heat", enh_R), meta_R:= "non-S2"]

subgroups <- as.matrix(dcast(pl, meta_L~meta_R, value.var = "log2FoldChange", fun.aggregate= median, na.rm= T), 1)
subgroups <- subgroups[order(rowMeans(subgroups)),]
subgroups <- subgroups[, order(colMeans(subgroups))]

pdf("pdf_sup/TWIST_screen_BA_categories_act_heatmap.pdf", 9, 9)
par(mar= rep(12, 4), cex= 0.9)
my_pheatmap(subgroups, legend_title = "log2(activity)", main= "Func. groups activity (median)",
            display_numbers = T, cluster_rows= F, cluster_cols= F)
dev.off()

