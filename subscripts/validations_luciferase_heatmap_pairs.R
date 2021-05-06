lib <- readRDS("Rdata/master_lib_features.rds")
dat <- readRDS("Rdata/validations_luciferase_final_table.rds")
dat <- dat[!grepl("^control", enh_L) & !grepl("^control", enh_R)]
dat[, obs_exp:= luc_norm/luc_add]

cols <- c("luc_norm", "luc_mean_L",  "luc_mean_R",  "luc_add", "obs_exp")
pl <- dat[, lapply(.SD, mean, na.rm= T), .(Sample_ID, enh_L, enh_R), .SDcols= cols]
pl[, luc_mean_L:= mean(luc_mean_L), enh_L]
pl[, luc_mean_R:= mean(luc_mean_R), enh_R]
pl <- na.omit(pl)
pl[lib, col_L:= i.col, on= "enh_L==ID"]
pl[lib, col_R:= i.col, on= "enh_R==ID"]

dmat <- dcast(pl, -luc_mean_L+enh_L~luc_mean_R+enh_R, value.var = "obs_exp", fun.aggregate = function(x) mean(log2(x)), fill = NA)
mat <- as.matrix(dmat[, -1], 1)
mat <- mat[apply(mat, 1, function(x) !all(is.na(x))), ]
mat <- mat[, apply(mat, 2, function(x) !all(is.na(x)))]

#----------------------------------------------------------------#
# 3- plot
#----------------------------------------------------------------#
pdf("pdf/Heatmap_luciferase_validations.pdf", 9.5, 7)
layout(matrix(1:9, ncol=3), widths= c(0.26,0.05,1), heights= c(0.35,0.065,1))

####-------LEFT INDIVIDUAL ACTIVITY, FEATURES and MARGINALS-------####
par(mar=c(0.25,10,0.5,0.25), xaxs="i", yaxs="i", las= 1)
# median
.c <- unique(pl[order(luc_mean_L), .(luc_mean_L= log2(luc_mean_L)), .(enh_L, col_L)])
plot.new()
plot.new()
bar <- barplot(-.c$luc_mean_L, horiz = T, border= NA, col= ifelse(.c$luc_mean_L<=0, "blue", "red"), axes=F, space = 0, 
               xlim= c(-6.2,0), names.arg = .c$enh_L)
segments(0, 0, 0, max(bar)+0.5, lty= 1, col= "grey")
segments(-1, 0, -1, max(bar)+0.5, lty= 2)
axis(3, at= c(-6, 0), labels = c(6, 0), line = 1, las= 1)

# group
par(mar=c(0.25,0.5,0.5,0.25), xaxs="i", yaxs="i")
plot.new()
plot.new()
barplot(rep(1, nrow(.c)), horiz = T, border= NA, col= .c$col_L, axes=F, space = 0)

####-------RIGHT INDIVIDUAL ACTIVITY, FEATURES and MARGINALS-------####
par(mar= c(0.25,0.5,9,12), xaxs="i", yaxs="i")
# median
.c <- unique(pl[order(luc_mean_R), .(luc_mean_R= log2(luc_mean_R)), .(enh_R, col_R)])
bar <- barplot(.c$luc_mean_R, border= NA, col= ifelse(.c$luc_mean_R <= 0, "blue", "red"), axes=F, space = 0, ylim= c(0,6.2))
text(bar, 7, labels = .c$enh_R, srt= 45, pos = 4, xpd= T, offset = -0.25)
segments(0, 0, max(bar)+0.5, 0, lty= 1, col= "grey")
segments(0, 1, max(bar)+0.5, 1, lty= 2)
axis(2, at= c(0, 6), labels = c(0, 6), line = 1, las= 1)

# group
par(mar= c(0.25,0.5,0.5,12), xaxs="i", yaxs="i")
barplot(rep(1, nrow(.c)), border= NA, col= .c$col_R, axes=F, space = 0)

####-------ADDITIVE SCORES-------####
# PAIR
vl_heatmap(mat = mat, 
           cluster_rows = F, 
           cluster_cols = F, 
           breaks= c(-1,0,1), 
           legend_title = "Additive\nscore (log2)", 
           display_numbers = T)

####-------LEGEND-------####
points(rep(1.11, 3), seq(0.445, 0.52, length.out = 3), pch=15, xpd= T, cex= 3, col= unique(c(pl$col_L, pl$col_R)))
text(rep(1.11, 3), seq(0.445, 0.52, length.out = 3), pos= 4, c("hk", "shared", "dev"), xpd= T, offset = 1.5, cex= 1.5)
text(1.0775, 0.6, "Enhancer\ntype", xpd= T, cex= 1.5, pos= 4)

dev.off()
