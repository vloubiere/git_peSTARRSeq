setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
sapply(list.files("/groups/stark/vloubiere/functions/", ".R$", full.names = T), source)
require(data.table)

STARR <- readRDS("Rdata/processed_peSTARRSeq_data/SCR1_peSTARRSeq_final_table.rds")
dat <- readRDS("Rdata/luciferase_validations/C_luc_validations_final_table.rds")
feat <- readRDS("Rdata/library/lib_features.rds")
feat <- feat[ID %in% c(dat$enh_L, dat$enh_R)]
feat <- feat[, groupID:= paste0(group, formatC(.SD[, .I], width = 3, flag = "00")), group]
dat[, ID:= paste0(feat[match(enh_L, feat$ID), groupID], "~", feat[match(enh_R, feat$ID), groupID])]
dat[, group:= paste0(feat[match(enh_L, feat$ID), group], "~", feat[match(enh_R, feat$ID), group])]
Cc <- data.table(group= c("hk~hk", "control~dev", "control~hk", "dev~control", "control~control", "dev~dev", "hk~dev", "hk~control", "dev~hk"),
                 col= c("red", "cornflowerblue", "chocolate1", "cyan3", "grey", "darkblue", "gold", "coral4", "green"))

#------------------------------------------------------------#
# 1- Add peSTARR-Seq data, groups and colors
#------------------------------------------------------------#
pl <- unique(dat[(grepl("^control", enh_L) & grepl("^control", enh_R)) | (!grepl("^control", enh_L) & !grepl("^control", enh_R)), .(ID, group, enh_L, enh_R)])
pl[group=="control~control", class:= "A_control"]
pl[STARR, class:= cut(i.diff, c(-Inf, -1, 1, Inf), c("B_sub_add", "C_add", "D_super_add")), on= c("enh_L", "enh_R")]
pl <- dat[pl, .(class, left= mean(luc_mean_L, na.rm= T), right= mean(luc_mean_R, na.rm= T), pair= mean(luc_norm, na.rm= T)), .EACHI, on= c("ID", "group")]
pl<- pl[, .(lab1= rep(c("left", "right", "pair"), 2), 
            lab2= rep(c("add", "obs"), each= 3),
            value= c(left, right, rep(NA, 3), pair)), .(ID, group, class)]
pl[is.na(value) & lab1 %in% c("right", "left"), value:= 0]
dmat <- dcast(pl, ID+group+lab2+class~lab1, value.var = "value")
ord <- dmat[, ord:= max(pair, na.rm= T), ID]$V1
dmat <- dmat[order(class, ord), !"ord"]
dmat[, ID_lab:= paste0(ID, "_", lab2)]
mat <- t(as.matrix(dmat[, .(ID_lab, left, right, pair)], 1))

pdf("pdf/luciferase_validations/barplot_luc_validations.pdf", 20, 7)
par(las= 2, mar= c(5,10,2,1))
Cc <- Cc[dmat[match(gsub("_add|_obs", "", dmat[, ID_lab, ID_lab][, ID_lab]), dmat$ID), group], col, on= "group"]
spa <- c(0.1, rep(c(0.1, 1), ncol(mat)/2)[-ncol(mat)])
bar <- barplot(mat[1:2,], border= NA, col= c("grey20", "grey70"), space= spa, xaxt= "n", ylim= c(0,80))
barplot(mat[3,], border= NA, col= Cc, space= spa, xaxt= "n", yaxt= "n", add= T)
text(rowMeans(matrix(bar, ncol= 2, byrow = T)), -2, srt= 45, labels = unique(gsub("_obs|_add", "", colnames(mat))), pos= 2, xpd= T, cex= 0.7, offset = -0.2)
dev.off()
