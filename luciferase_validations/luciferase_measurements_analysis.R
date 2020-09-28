setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
sapply(list.files("/groups/stark/vloubiere/functions/", ".R$", full.names = T), source)
require(data.table)
require(readxl)

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
dat_all[, luc_mean_L:= mean(.SD[grepl("^control", enh_R), luc_norm], na.rm= T), .(replicate, enh_L)]
dat_all[, luc_mean_R:= mean(.SD[grepl("^control", enh_L), luc_norm], na.rm= T), .(replicate, enh_R)]
dat_all[, luc_add:= luc_mean_L+luc_mean_R]
dat_all[grepl("^control", enh_L)|grepl("^control", enh_R), luc_add:= NA]
# Collapse on enhancer combinations
dat <- dat_all[, .(luc_norm_all= list(luc_norm), luc_norm_mean= mean(luc_norm, na.rm= T),
                   luc_mean_L_all= list(luc_mean_L), luc_mean_L_mean= mean(luc_mean_L, na.rm= T),
                   luc_mean_R_all= list(luc_mean_R), luc_mean_R_mean= mean(luc_mean_R, na.rm= T),
                   luc_add_all= list(unique(.SD[, .(luc_add, replicate)])$luc_add), 
                   luc_add_mean= mean(unique(.SD[, .(luc_add, replicate)])$luc_add, na.rm= T)), .(Sample_ID, enh_L, enh_R)]

#------------------------------------------------------------#
# 3- Add peSTARR-Seq data
#------------------------------------------------------------#
# peSTARR-Seq data
peSTARR <- readRDS("Rdata/processed_peSTARRSeq_data/SCR1_peSTARRSeq_final_table.rds")
cols <- c("peSTARR_FC", "peSTARR_median_L", "peSTARR_median_R", "peSTARR_add", "peSTARR_diff")
dat[peSTARR, (cols):= .(i.log2FoldChange, i.median_L, i.median_R, i.log2FC_add, i.diff), on= c("enh_L", "enh_R")]

#------------------------------------------------------------#
# 4- PCC scatterplot
#------------------------------------------------------------#
pl <- copy(dat)
pl[, luc_norm_mean:= log2(luc_norm_mean)]
pl[, luc_norm_all:= lapply(luc_norm_all, log2)]
pl[, luc_norm_sd:= sapply(luc_norm_all, sd)]
pl[grepl("^control", enh_L) & grepl("^control", enh_R), c("col", "class"):= .("grey", "control~control")]
pl[grepl("^dev", enh_L) & grepl("^control", enh_R), c("col", "class"):= .("cyan3", "dev~control")]
pl[grepl("^control", enh_L) & grepl("^dev", enh_R), c("col", "class"):= .("cornflowerblue", "control~dev")]
pl[grepl("^dev", enh_L) & grepl("^dev", enh_R), c("col", "class"):= .("darkblue", "dev~dev")]
pl[grepl("^hk", enh_L) & grepl("^control", enh_R), c("col", "class"):= .("coral4", "hk~control")]
pl[grepl("^control", enh_L) & grepl("^hk", enh_R), c("col", "class"):= .("chocolate1", "control~hk")]
pl[grepl("^hk", enh_L) & grepl("^hk", enh_R), c("col", "class"):= .("red", "hk~hk")]
pl[grepl("^hk", enh_L) & grepl("^dev", enh_R), c("col", "class"):= .("gold", "hk~dev")]
pl[grepl("^dev", enh_L) & grepl("^hk", enh_R), c("col", "class"):= .("green", "dev~hk")]
leg <- unique(pl[, .(class, col)])[c(4,3,2,5,6,7,8,1)]

pdf("pdf/luciferase_validations/correlation_luc_peSTARRSeq.pdf", 6, 6)
par(pty="s")

plot(pl[, .(peSTARR_FC, luc_norm_mean)], las= 1, col= pl$col, pch= 19,
     xlab= "pe-STARR-Seq activity (log2)", ylab= "Normalized luciferase activity (log2)")

# sd
arrows(pl[, peSTARR_FC], pl[, luc_norm_mean], pl[, peSTARR_FC], pl[, luc_norm_mean+luc_norm_sd], 
       length = 0.025, angle = 90, col= pl$col)
arrows(pl[, peSTARR_FC], pl[, luc_norm_mean], pl[, peSTARR_FC], pl[, luc_norm_mean-luc_norm_sd], 
       length = 0.025, angle = 90, col= pl$col)

# lm
.lm1 <- lm(luc_norm_mean~peSTARR_FC, pl)
abline(.lm1, lty= 1)
rsq <- paste("R²= ", round(summary(.lm1)$r.squared, 2))

legend("topleft", pch= c(rep(19, 8), NA), bty= "n", col= c(leg$col, "black"), legend= c(leg$class, rsq), 
       lty= c(rep(NA, 8), 1), cex= 0.8)
dev.off()

#------------------------------------------------------------#
# 5- Barplot combinations
#------------------------------------------------------------#
pl2 <- copy(dat)
pl2 <- pl2[!xor(grepl("control", enh_L), grepl("control", enh_R))]
pl2[, luc_norm_sd:= sapply(luc_norm_all, sd, na.rm= T)]
pl2[, luc_add_sd:= sapply(luc_add_all, sd, na.rm= T)]
pl2[grepl("^control", enh_L) & grepl("^control", enh_R), c("col", "class"):= .("grey", "control~control")]
pl2[grepl("^dev", enh_L) & grepl("^dev", enh_R), c("col", "class"):= .("darkblue", "dev~dev")]
pl2[grepl("^hk", enh_L) & grepl("^hk", enh_R), c("col", "class"):= .("red", "hk~hk")]
pl2[grepl("^hk", enh_L) & grepl("^dev", enh_R), c("col", "class"):= .("gold", "hk~dev")]
pl2[grepl("^dev", enh_L) & grepl("^hk", enh_R), c("col", "class"):= .("green", "dev~hk")]
pl2[grepl("^dev", enh_L), c("col_L", "class_L"):= .("cyan3", "dev left")]
pl2[grepl("^hk", enh_L), c("col_L", "class_L"):= .("chocolate1", "hk left")]
pl2[grepl("^dev", enh_R), c("col_R", "class_R"):= .("cornflowerblue", "dev right")]
pl2[grepl("^hk", enh_R), c("col_R", "class_R"):= .("brown4", "hk right")]
# add predicted group
pl2[class=="control~control", group:= "1_control"]
pl2[class!="control~control" & peSTARR_diff < 1, group:= "2_additive"]
pl2[class!="control~control" & (peSTARR_diff >= 1 | is.na(peSTARR_diff)), group:= "3_super-additive"]
pl2 <- pl2[order(group, luc_norm_mean)]
pl2[class=="control~control", c("luc_mean_L_mean", "luc_add_mean"):= NA]

leg <- rbind(unique(pl2[, .(col, class)]),
             unique(pl2[, .(col= col_L, class= class_L)]),
             unique(pl2[, .(col= col_R, class= class_R)]))
leg <- na.omit(leg)
leg <- leg[c(1,2,3,4,6,7,5,8)]

#######
pdf("pdf/luciferase_validations/barplot_luc_validations.pdf", 30, 7)
par(mar= c(4,7,2,2))

# Mean luciferase activity
space <- pl2[, .N-1, group][, V1]
spa1 <- c(0, rep(0.5, space[1]), 4, rep(1.5, space[2]), 4, rep(1.5, space[3]))
bar1 <- barplot(pl2$luc_norm_mean, col= pl2$col, las= 1, lwd= 0.25, ylim= c(0, 120), space = spa1, ylab= "Normalized luciferase activity / control pairs")

xpos <- jitter(rep(bar1, lengths(pl2[, luc_norm_all])))
points(xpos, unlist(pl2[, luc_norm_all]), pch= 19, cex= 0.2, col= "grey70")
arrows(bar1, pl2[, luc_norm_mean], bar1, pl2[, luc_norm_mean+luc_norm_sd], length = 0.025, angle = 90)
arrows(bar1, pl2[, luc_norm_mean], bar1, pl2[, luc_norm_mean-luc_norm_sd], length = 0.025, angle = 90)

# Mean Additive score
spa2 <- c(0, rep(0.5, space[1]), 3, rep(1.5, space[2]), 4, rep(1.5, space[3]))
barplot(pl2$luc_add_mean, col= pl2$col_R, las= 1, add= T, lwd= 0.25, axes= F, space = spa2)
# Left enhancer activity
bar2 <- barplot(pl2$luc_mean_L_mean, col= pl2$col_L, las= 1, add= T, axes= F, space = spa2)
barplot(pl2$luc_mean_L_mean, las= 1, border= NA, add= T, density = 30, angle = 45, col= "black", 
        axes= F, space = spa2)

xpos <- rep(bar2, lengths(pl2[, luc_add_all]))
points(jitter(xpos), unlist(pl2[, luc_add_all]), pch= 19, cex= 0.2, col= "grey70")
arrows(bar2, pl2[, luc_add_mean], bar2, pl2[, luc_add_mean+luc_add_sd], length = 0.025, angle = 90)
arrows(bar2, pl2[, luc_add_mean], bar2, pl2[, luc_add_mean-luc_add_sd], length = 0.025, angle = 90)

# xlab legend
text(bar1[1]-8, -5, "pe-STARR-Seq assignment ->", xpd= T, pos= 1)
x1 <- bar2[1]
x2 <- bar1[space[1]+1]
segments(x1, -5, x2, -5, xpd= T)
text(mean(c(x1, x2)), -5, "control pairs", xpd= T, pos= 1)

x3 <- bar2[space[1]+2]
x4 <- bar1[sum(space[1:2])+2]
segments(x3, -5, x4, -5, xpd= T)
text(mean(c(x3, x4)), -5, "additive pairs", xpd= T, pos= 1)

x5 <- bar2[sum(space[1:2])+3]
x6 <- bar1[sum(space[1:3])+3]
segments(x5, -5, x6, -5, xpd= T)
text(mean(c(x5, x6)), -5, "super-additive pairs", xpd= T, pos= 1)

# pval
pval <- pl2[, ifelse(length(na.omit(unlist(luc_norm_all)))>2 & length(na.omit(unlist(luc_add_all)))>2, 
                     t.test(unlist(na.omit(luc_norm_all)), unlist(na.omit(luc_add_all)))$p.value, 1), Sample_ID][, V1]
pval[pl2$class=="control~control"] <- NA
pval_f <- my_pval_format(pval)$p.value
Cc <- ifelse(pl2$luc_add_mean-pl2$luc_norm_mean > 0, "blue", "red")
Cc[pval_f=="N.S"] <- "black"
pval_cex <- my_pval_format(pval)$cex
xpos <- rowMeans(cbind(bar1, bar2))
ypos <- pl2[, max(c(unlist(luc_add_all), unlist(luc_norm_all)), na.rm= T), Sample_ID][, V1]
tab <- data.table(xpos, ypos, pval_f, pval_cex, Cc, offset= ifelse(is.na(pval_f) | pval_f=="N.S", 1.1, 0.75))

tab[, text(xpos, ypos, pval_f, cex= pval_cex, pos= 3, col= Cc, offset= offset), (tab)]

# FC
FC <- round(pl2$luc_norm_mean/pl2$luc_add_mean, 1)
FC <- ifelse(is.na(FC), NA, paste0("x", FC))
text(xpos, ypos, FC, cex= 0.7, pos= 3, offset = 0.5, col= Cc)

# Legend
legend("topleft", bty= "n", fill= leg$col, legend= leg$class)

dev.off()













