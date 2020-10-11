setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
sapply(list.files("/groups/stark/vloubiere/functions/", ".R$", full.names = T), source)
require(data.table)
require(readxl)
require(circlize)

#------------------------------------------------------------#
# 1- Import and format data
#------------------------------------------------------------#
# Plate schemes
scheme <- data.table(file= list.files("db/luciferase/peSTARRSeq_validations/", "scheme", recursive = T, full.names = T))
scheme[, c("date", "plate") := tstrsplit(basename(file), "_|[.]", keep= c(1,3)), file]
scheme <- scheme[, melt(fread(file, header = T, colClasses = "character", na.strings = "", fill= T), 
                        id.vars = "rownames", value.name = "Sample_ID"), (scheme)]
colnames(scheme)[4:5] <- c("row", "col")
# Add technical replicates
rep <- fread("Rdata/luciferase_validations/plates_replicates_grid.csv", header = T, colClasses = "character", na.strings = "", fill= T)
rep <- melt(rep, id.vars = "rownames")
colnames(rep) <- c("row", "col", "tech_replicate")
scheme <- merge(scheme, rep)
scheme <- scheme[, .SD[, .(row, col, replicate= .GRP), .(tech_replicate, date, plate)], Sample_ID]
scheme <- scheme[, !"tech_replicate"]
# ID/sample correspondance
ID <- as.data.table(read_excel("Rdata/luciferase_validations/clean_stocks.xlsx"))
scheme[ID, c("enh_L", "enh_R"):= .(i.enh_L, i.enh_R), on= "Sample_ID"]
# luciferase data
luc <- data.table(file= list.files("db/luciferase/peSTARRSeq_validations/", "peSTARRvalid", recursive = T, full.names = T))
luc[, c("date", "plate", "lum") := tstrsplit(basename(file), "_|[.]", keep= c(1,3,4)), file]
luc <- luc[, melt(fread(file, header = T)[, 1:25], id.vars = "V1"), (luc)]
colnames(luc)[5:6] <- c("row", "col")
luc <- dcast(luc, date+plate+row+col~lum)
colnames(luc)[5:6] <- c("luc", "ren")

#------------------------------------------------------------#
# 2- Process luciferase data
#------------------------------------------------------------#
merged <- merge(scheme, luc, by= c("date", "plate", "row", "col"))
# Cutoffs rennilla and N tech replicates
dat_all <- merged[ren>7500]
dat_all[, check := .N>=3 & !is.na(Sample_ID), Sample_ID]
dat_all <- dat_all[(check), !"check"]
# Mean replicates
dat_all <- dat_all[, .(luc_norm= mean(luc/ren)), .(Sample_ID, enh_L, enh_R, replicate)]
# Normalize for negative controls
dat_all[, luc_norm:= luc_norm/mean(dat_all[grepl("^control", enh_L) & grepl("^control", enh_R), luc_norm])]
# Compute additive scores
dat_all[, luc_mean_L:= mean(.SD[grepl("^control", enh_R), luc_norm], na.rm= T), .(enh_L)]
dat_all[, luc_mean_R:= mean(.SD[grepl("^control", enh_L), luc_norm], na.rm= T), .(enh_R)]
dat_all[, luc_add:= luc_mean_L+luc_mean_R]

dat <- dat_all[!grepl("^control", enh_L) & !grepl("^control", enh_R), lapply(.SD, mean, na.rm= T), .(enh_L, enh_R), 
               .SDcols= c("luc_mean_L", "luc_mean_R", "luc_add", "luc_norm")]
dat[, ratio:= luc_norm/luc_add]

# dcast/melt
dmat <- dcast(dat, -luc_mean_L+enh_L~luc_mean_R+enh_R, value.var = "ratio", fill = NA)
mat <- as.matrix(dmat[, -c(1:2)])
pl <- melt(dmat, id.vars = c("luc_mean_L", "enh_L"), value.name = "ratio", variable.name = "enh_R")
pl[, luc_mean_L:= -luc_mean_L]
pl[, "luc_mean_R":= tstrsplit(enh_R, "_", keep = 1)]
pl[, "enh_R":= gsub(paste0(luc_mean_R, "_"), "", enh_R), .(enh_R, luc_mean_R)]
pl[, luc_mean_R:= as.numeric(luc_mean_R)]
pl[grepl("^dev", enh_L), col_L:= "#74C27A"]
pl[grepl("^hk", enh_L), col_L:= "tomato"]
pl[grepl("^dev", enh_R), col_R:= "#74C27A"]
pl[grepl("^hk", enh_R), col_R:= "tomato"]

#----------------------------------------------------------------#
# 3- plot
#----------------------------------------------------------------#
# pdf
pdf("pdf/luciferase_validations/Heatmap_luciferase_validations.pdf", 9.5, 7)
layout(matrix(1:9, ncol=3), widths= c(0.26,0.05,1), heights= c(0.35,0.065,1))

####-------LEFT INDIVIDUAL ACTIVITY, FEATURES and MARGINALS-------####
par(mar=c(0.25,10,0.5,0.25), xaxs="i", yaxs="i", las= 1)
# median
.c <- log2(pl[, luc_mean_L, .(enh_L, luc_mean_L)]$luc_mean_L)
plot.new()
plot.new()
bar <- barplot(-rev(.c), horiz = T, border= NA, col= ifelse(rev(.c)<=0, "blue", "red"), axes=F, space = 0, 
               xlim= c(-6.2,0), names.arg = pl[, enh_L, enh_L]$enh_L)
segments(0, 0, 0, max(bar)+0.5, lty= 1, col= "grey")
segments(-1, 0, -1, max(bar)+0.5, lty= 2)
axis(3, at= c(-6, 0), labels = c(6, 0), line = 1, las= 1)

# group
par(mar=c(0.25,0.5,0.5,0.25), xaxs="i", yaxs="i")
.c <- rev(pl[, col_L, enh_L]$col_L)
plot.new()
plot.new()
barplot(rep(1, length(.c)), horiz = T, border= NA, col= .c, axes=F, space = 0)

####-------RIGHT INDIVIDUAL ACTIVITY, FEATURES and MARGINALS-------####
par(mar= c(0.25,0.5,9,12), xaxs="i", yaxs="i")
# median
.c <- log2(pl[, luc_mean_R, .(enh_R, luc_mean_R)]$luc_mean_R)
bar <- barplot(.c, border= NA, col= ifelse(.c <= 0, "blue", "red"), axes=F, space = 0, ylim= c(0,6.2))
text(bar, 7, labels = pl[, enh_R, enh_R]$enh_R, srt= 45, pos = 4, xpd= T, offset = -0.25)
segments(0, 0, max(bar)+0.5, 0, lty= 1, col= "grey")
segments(0, 1, max(bar)+0.5, 1, lty= 2)
axis(2, at= c(0, 6), labels = c(0, 6), line = 1, las= 1)

# group
par(mar= c(0.25,0.5,0.5,12), xaxs="i", yaxs="i")
.c <- pl[, col_R, enh_R]$col_R
barplot(rep(1, length(.c)), border= NA, col= .c, axes=F, space = 0)

####-------ADDITIVE SCORES-------####
# PAIR
my_pheatmap(log2(mat), cluster_rows = F, cluster_cols = F, breaks= c(-1,0,1), legend_cex = 1.5,
            legend_title= "Additive\nscore (log2)", col= c("cornflowerblue", "white", "tomato"), 
            legend_adj_top = -0.04, legend_adj_left = -0.04, display_numbers = T, row_labels = F, col_labels = F)

####-------LEGEND-------####
points(rep(1.11, 2), seq(0.47, 0.52, length.out = 2), pch=15, xpd= T, cex= 3, col= c("tomato", "#74C27A"))
text(rep(1.11, 2), seq(0.47, 0.52, length.out = 2), pos= 4, c("hk", "dev"), xpd= T, offset = 1.5, cex= 1.5)
text(1.0775, 0.6, "Enhancer\ntype", xpd= T, cex= 1.5, pos= 4)

dev.off()
