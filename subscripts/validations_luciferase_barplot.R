# Import
dat <- readRDS("Rdata/validations_luciferase_final_table.rds")
# Add features
feat <- readRDS("Rdata/master_lib_features.rds")
feat <- feat[ID %in% c(dat$L, dat$R)]
feat[, groupID:= paste0(group, formatC(.SD[, .I], width = 3, flag = "00")), group]
dat[, ID:= paste0(feat[match(L, feat$ID), groupID], "~", feat[match(R, feat$ID), groupID])]
dat[feat, c("group_L", "col_L"):= .(i.group, i.col), on= "L==ID"]
dat[feat, c("group_R", "col_R"):= .(i.group, i.col), on= "R==ID"]
dat[, group:= paste0(group_L, "~", group_R)]
dat[group_L==group_R, col:= col_L, .(col_L, col_R)]
dat[group_L!=group_R, col:= colorRampPalette(c(col_L, col_R))(3)[2], .(col_L, col_R)]
# Mean values / replicates
cols <- c("luc_norm", "luc_mean_L", "luc_mean_R", "luc_add")
dat[, paste0(cols, "_mean"):= lapply(.SD, mean, na.rm=T), .(Sample_ID, ID, group, L, R), .SDcols= cols]
dat[, paste0(cols, "_sd"):= lapply(.SD, sd, na.rm=T), .(Sample_ID, ID, group, L, R), .SDcols= cols]
# Compute pval additive and observed
dat[, pval:= t.test(na.omit(.SD[, .(luc_add, luc_norm)])[[1]],
                    na.omit(.SD[, .(luc_add, luc_norm)])[[2]], paired = T)$p.value, .(L, R)]
dat[, pval:= cut(pval, c(-Inf, 0.00001, 0.001, 0.01, 0.05, Inf), c("****", "***", "**", "*", "N.S"))]
# filter and collpase
dat <- dat[!xor(grepl("^control", L), grepl("^control", R))]
sel <- c("Sample_ID", "L", "R", "ID", "group", "col", "pval", grep("mean$|sd$", colnames(dat), value = T))
dat <- dat[, lapply(.SD, function(x) list(na.omit(x))), sel, .SDcols= cols]
# Add class
# STARR <- readRDS("Rdata/master_results_peSTARRSeq.rds")[lib=="libvl002", !"lib"] # Used the old data cause the new calculations are more stringent -> some pairs are lost :(
STARR <- readRDS("../pe_STARRSeq/Rdata/processed_peSTARRSeq_data/SCR1_peSTARRSeq_final_table.rds")
dat[STARR, class:= i.diff, on= c("L==enh_L", "R==enh_R")]
dat[, class:= ifelse(class>1, "D_super-additive", ifelse(class< -1, "B_sub-additive", "C_additive"))]
dat[grepl("^control", L) & grepl("^control", R), class:= "A_control"]
dat <- dat[order(class, luc_norm_mean)]

pdf("pdf/barplot_luc_validations.pdf", 15, 5)
par(las= 2, mar= c(5,4,2,1))

# Barplot
ylim <- c(0, 110)
bar1 <- barplot(t(as.matrix(dat[, .(luc_mean_L_mean, luc_mean_R_mean)])), border= NA, space = 1.6, 
               col= c("grey20", "grey70"), xaxt= "n", ylim= ylim, ylab= "Normalized luciferase activity", xlim= c(6, 190))
bar2 <- barplot(dat$luc_norm_mean, space = c(2.6, rep(1.6, nrow(dat)-1)), border= dat$col, density = 50, angle = 45,
        col= dat$col, xaxt= "n", ylim= ylim, ylab= "Normalized luciferase activity", add= T)
bar2 <- bar2[,1]
abline(h= -0.2, lwd= 2)
text(rowMeans(matrix(c(bar1, bar2), ncol= 2)), -2, srt= 45, labels = dat$ID, pos= 2, xpd= T, cex= 0.7, offset = -0.2)

# sd
y1 <- dat[, luc_add_mean+luc_add_sd]
segments(bar1, y1, bar1, dat[, luc_add_mean-luc_add_sd])
y2 <- dat[, luc_norm_mean+luc_norm_sd]
segments(bar2, y2, bar2, dat[, luc_norm_mean-luc_norm_sd])

# p.value
x <- rowMeans(matrix(c(bar1, bar2), ncol= 2))
y <- apply(matrix(c(y1, y2), ncol= 2), 1, max)
text(x, y, dat$pval, pos= 3, offset= 0, cex= 0.6, 
     col= ifelse(dat[, luc_add_mean<luc_norm_mean], "red", "blue"))

# STARR-Seq prediction
cl <- cumsum(table(dat$class))
x1 <- bar1[c(1, cl[-length(cl)]+1)]
x2 <- bar2[cl]
ycl <- ifelse(unique(dat$class)=="A_control", 10, ylim[2]-10)
segments(x1, ycl, x2, ycl, xpd= T)
text(rowMeans(matrix(c(x1, x2), ncol=2)), ycl, gsub("A_|B_|C_|D_", "", unique(dat$class)), pos=3, xpd= T)
text(20, ylim[2]-10, "STARR-Seq prediction -->", pos=3, xpd= T)

# legend
Cc <- unique(dat[, .(group, col)])
legend(0, ylim[2]-10, legend = c("LEFT~control", "control~RIGHT"), fill= c("grey20", "grey70"), bty= "n")
legend(0, ylim[2]-30, legend = Cc$group, fill= Cc$col, bty= "n", density = 50, angle = 45)

dev.off()
