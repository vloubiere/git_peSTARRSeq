# setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
setwd("/Users/vincent.loubiere/Dropbox (VBC)/untitled folder/")
# sapply(list.files("/groups/stark/vloubiere/functions/", ".R$", full.names = T), source)
sapply(list.files("functions/", ".R$", full.names = T), source)
require(data.table)

# Import
dat <- readRDS("Rdata/luciferase_validations/C_luc_validations_final_table.rds")
feat <- readRDS("Rdata/library/lib_features.rds")
feat <- feat[ID %in% c(dat$enh_L, dat$enh_R)]
feat <- feat[, groupID:= paste0(group, formatC(.SD[, .I], width = 3, flag = "00")), group]
# Add features
dat[, ID:= paste0(feat[match(enh_L, feat$ID), groupID], "~", feat[match(enh_R, feat$ID), groupID])]
dat[, group:= paste0(feat[match(enh_L, feat$ID), group], "~", feat[match(enh_R, feat$ID), group])]
# Mean values / replicates
cols <- c("luc_norm", "luc_mean_L", "luc_mean_R", "luc_add")
dat[, paste0(cols, "_mean"):= lapply(.SD, mean, na.rm=T), .(Sample_ID, ID, group, enh_L, enh_R), .SDcols= cols]
dat[, paste0(cols, "_sd"):= lapply(.SD, sd, na.rm=T), .(Sample_ID, ID, group, enh_L, enh_R), .SDcols= cols]
# filter
dat <- dat[!xor(grepl("^control", enh_L), grepl("^control", enh_R))]
sel <- c("Sample_ID", "enh_L", "enh_R", "ID", "group", grep("mean$|sd$", colnames(dat), value = T))
dat <- dat[, lapply(.SD, function(x) list(na.omit(x))), sel, .SDcols= cols]
# Add colors
dat[data.table(group= c("hk~hk", "control~dev", "control~hk", "dev~control", "control~control", "dev~dev", "hk~dev", "hk~control", "dev~hk"),
                 col= c("red", "cornflowerblue", "chocolate1", "cyan3", "grey", "darkblue", "gold", "coral4", "green")), col:= i.col, on= "group"]
# Add class
STARR <- readRDS("Rdata/processed_peSTARRSeq_data/SCR1_peSTARRSeq_final_table.rds")
dat[STARR, class:= ifelse(i.diff>1, "D_super-additive", ifelse(i.diff<1, "B_sub-additive", "C_additive")), on= c("enh_L", "enh_R")]
dat[grepl("^control", enh_L) & grepl("^control", enh_R), class:= "A_control"]
dat <- dat[order(class, luc_norm_mean)]


pdf("pdf/luciferase_validations/barplot_luc_validations.pdf", 25, 5)
par(las= 2, mar= c(5,10,2,1))
ylim <- c(0, 110)
# barplot
barplot(t(as.matrix(dat[, .(luc_mean_L_mean, luc_mean_R_mean)])), border= NA, space = 1.6, 
               col= c("grey20", "grey70"), xaxt= "n", ylim= ylim, ylab= "Normalized luciferase activity")
barplot(dat$luc_norm_mean, space = c(2.6, rep(1.6, nrow(dat)-1)), border= dat$col, density = 50, angle = 45,
        col= dat$col, xaxt= "n", ylim= ylim, ylab= "Normalized luciferase activity", add= T)

dev.off()


spa <- c(0.1, rep(c(0.1, 0.5), nrow(dat))[-nrow(dat)])
bar <- barplot(mat[1:2,], border= NA, col= c("grey20", "grey70"), space= spa, xaxt= "n", xlim= c(6, 185), ylim= ylim, ylab= "Normalized luciferase activity")
barplot(mat[3,], border= ifelse(seq(ncol(mat)) %% 2 !=0, NA, col), col= col, space= spa, xaxt= "n", yaxt= "n", add= T, density = 50, angle = 45)
# sd
m1 <- dat[unique(dmat$ID), mean(rowSums(.SD), na.rm= T), ID, on= "ID", .SDcols= c("luc_mean_L", "luc_mean_R")]$V1
sd1 <- dat[unique(dmat$ID), sd(rowSums(.SD), na.rm= T), ID, on= "ID", .SDcols= c("luc_mean_L", "luc_mean_R")]$V1
segments(bar[(seq(bar) %% 2)!=0], m1+sd1,
         bar[(seq(bar) %% 2)!=0], m1-sd1)
m2 <- dat[unique(dmat$ID), mean(luc_norm, na.rm= T), ID, on= "ID"]$V1
sd2 <- dat[unique(dmat$ID), sd(luc_norm, na.rm= T), ID, on= "ID"]$V1
segments(bar[(seq(bar) %% 2)==0], m2+sd2,
         bar[(seq(bar) %% 2)==0], m2-sd2)
# pvalue
x <- rowMeans(matrix(bar, ncol=2, byrow = T))
y <- apply(matrix(c(m1+sd1, m2+sd2), ncol= 2),1,max)
pval <- dat[unique(dmat$ID), t.test(.SD[[1]], .SD[[2]], paired = T)$p.value, ID, on= "ID", .SDcols= c("luc_mean_L", "luc_mean_R")]$V1
pval <- cut(pval, c(-Inf, 0.00001, 0.001, 0.01, 0.05, Inf), c("****", "***", "**", "*", "N.S"))
text(x[pval!="N.S"], y[pval!="N.S"], pval[pval!="N.S"], pos= 3, offset= 0, col= ifelse((m1>m2)[pval!="N.S"], "blue", "red"))
text(x[pval=="N.S"], y[pval=="N.S"], pval[pval=="N.S"], cex= 0.6, pos= 3, offset= 0.3)
# peSTARRSeq prediction
text(rowMeans(matrix(bar, ncol= 2, byrow = T)), -2, srt= 45, labels = unique(gsub("_obs|_add", "", colnames(mat))), pos= 2, xpd= T, cex= 0.7, offset = -0.2)
cl <- c(bar[1], bar[cumsum(table(dmat$class))])
segments(cl[-length(cl)]+0.3, ifelse(unique(dmat$class)=="A_control", 10, ylim[2]-10), cl[-1]-0.3, ifelse(unique(dmat$class)=="A_control", 10, ylim[2]-10), xpd= T)
text(cl[-length(cl)]+diff(cl)/2, ifelse(unique(dmat$class)=="A_control", 10, ylim[2]-10), gsub("A_|B_|C_|D_", "", unique(dmat$class)), pos=3, xpd= T)
text(28, ylim[2]-10, "STARR-Seq prediction -->", pos=3, xpd= T)
# legend
legend(bar[1], ylim[2]-20, legend = c("LEFT~control", "control~RIGHT"), fill= c("grey20", "grey70"), bty= "n")
Cc <- Cc[group %in% dmat$group][order(match(group, dmat$group))]
legend(15, ylim[2]-20, legend = Cc$group, fill= Cc$col, bty= "n", density = 50, angle = 45)
