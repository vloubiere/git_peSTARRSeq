setwd("/groups/stark/vloubiere/projects/0006_vllib002_SCR1_validations_luc/")
sapply(list.files("/groups/stark/vloubiere/functions/", ".R$", full.names = T), source)
require(data.table)

# Import plate schemes
scheme <- data.table(file= list.files("luciferase", "scheme", recursive = T, full.names = T))
scheme[, c("date", "plate") := tstrsplit(basename(file), "_|[.]", keep= c(1,3)), file]
scheme <- scheme[, melt(fread(file, header = T, colClasses = "character", na.strings = "", fill= T), id.vars = "rownames", value.name = "Sample_ID"), (scheme)]
colnames(scheme)[4:5] <- c("row", "col")

# Import data
luc <- data.table(file= list.files("luciferase", "peSTARRvalid", recursive = T, full.names = T))
luc[, c("date", "plate", "lum") := tstrsplit(basename(file), "_|[.]", keep= c(1,3,4)), file]
luc <- luc[, melt(fread(file, header = T)[, 1:25], id.vars = "V1"), (luc)]
colnames(luc)[5:6] <- c("row", "col")
luc[scheme, Sample_ID:= i.Sample_ID, on= c("date", "plate", "row", "col")]
luc <- dcast(luc, date+plate+row+col+Sample_ID~lum)
# add enh_L and enh_R
ID <- readRDS("Rdata/sum_up.rds")
ID <- ID[, .(Sample_ID= unlist(Sample_ID)), setdiff(colnames(ID), "Sample_ID")]
luc[ID, c("enh_L", "enh_R", "type"):= .(i.enh_L, i.enh_R, i.type), on= "Sample_ID"]
luc[, pair_name:= paste0(enh_L, " x ", enh_R), .(enh_L, enh_R)]
# checks for low renilla / N replicates
luc <- luc[Sample_ID!= "NA"]
luc <- luc[, check := (lum2 >= boxplot(lum2, plot= F)$stats[2,] & lum2 >= 5000), .(date, plate)]
luc <- luc[(check)]
luc <- luc[, check := .N>2, .(date, plate, Sample_ID)]
luc <- luc[(check), !"check"]
# Compute mean activity
luc <- luc[, .(log2FC_luc= mean(log2(lum1)-log2(lum2))), .(date, plate, Sample_ID, enh_L, enh_R, type, pair_name)]
luc[, check:= .N>2, .(enh_L, enh_R)]
luc <- luc[(check), !"check"]
luc <- luc[, .(log2FC_luc= mean(log2FC_luc)), .(Sample_ID, enh_L, enh_R, type, pair_name)]
# Scale on negative controls
luc[, log2FC_luc:= log2FC_luc-median(luc[type=="control", log2FC_luc])]
# Compute expected additive
luc[, exp_luc_L:= mean(.SD[grepl("^control", enh_R), log2FC_luc]), enh_L]
luc[, exp_luc_R:= mean(.SD[grepl("^control", enh_L), log2FC_luc]), enh_R]
luc[, exp_luc:= log2(2^exp_luc_L + 2^exp_luc_R)]
# Add class
dat <- readRDS("/groups/stark/vloubiere/projects/pe_STARRSeq/Rdata/B_SCR1_peSTARRSeq_final_table.rds")
dat[, class:= ifelse(diff>2, "super-additive", "additive")]
luc[dat, class:= i.class, on= c("enh_L", "enh_R")]
luc[type != "enh_pair", class:= NA]
luc[type == "control", class:= "control"]
luc <- luc[(exp_luc_L>2 & exp_luc_R>2) | (grepl("^control", enh_L) & grepl("^control", enh_R))]
# luc <- luc[exp_luc_L>1 & exp_luc_R>1]

# Plot
pl <- melt(luc, id.vars = c("Sample_ID", "enh_L", "enh_R", "type", "pair_name", "class"), measure.vars = c("log2FC_luc", "exp_luc"))
pl <- pl[!grepl("^hk", enh_L) & !grepl("^hk", enh_R)]
pl[type=="control" & variable=="log2FC_luc", x:= 1]
pl[class=="additive" & variable=="exp_luc", x:= 2]
pl[class=="additive" & variable=="log2FC_luc", x:= 3]
pl[class=="super-additive" & variable=="exp_luc", x:= 4]
pl[class=="super-additive" & variable=="log2FC_luc", x:= 5]
pl[, col:="grey"]
pl[grepl("^dev", enh_L) & grepl("^dev", enh_R), col:= "cornflowerblue"]
pl[grepl("^hk", enh_L) | grepl("^hk", enh_R), col:="chocolate1"]
pl[grepl("^hk", enh_L) & grepl("^hk", enh_R), col:="red"]
set.seed(1234)
pl[, xjit:= jitter(x, 0.75)]
# pl[, value:= 10*value]
pl[, value:= 2^value]

pdf("pdf/stripchart_luciferase_validations.pdf", 9, 9)
my_empty_plot(xlim = c(0.5, 5.5), ylim= c(0, 90), xaxt= "n", xlab= "", ylab= "Normalized luciferase (Firefly/Renilla)")
points(pl$xjit, pl$value, pch= 19, cex= 0.8, col= pl$col)
# text(pl$xjit, pl$value, pl$Sample_ID)
pl[, segments(xjit[1], value[1], xjit[2], value[2], col= adjustcolor(col[1], 0.7)), pair_name]
axis(1, at = 1:5, labels = c("control\n pairs", rep(c("expected\n additive score", "observed"), 2)), padj= 0.35)
segments(c(1.75, 3.75), -11, c(3.25, 5.25), -11, xpd= T)
# legend("topleft", pch= 19, col= c("grey", "chocolate1", "cornflowerblue"), c("x2 negative controls pair", "x2 developmental enhancers pair", "hk. x dev. enhancers pair"), bty= "n")
legend("topleft", pch= 19, col= c("grey", "cornflowerblue"), c("x2 negative controls pair", "x2 developmental enhancers pair"), bty= "n")
text(c(2.5, 4.5), -11.2, c("ADDITIVE\n PAIRS", "SUPER-ADDITIVE\n PAIRS"), pos= 1, xpd= T)
# polygon(c(2, 2, 3, 3, 3, 2), c(68, 70, 70, 68, 70, 70))
# polygon(c(4, 4, 5, 5, 5, 4), c(83, 85, 85, 83, 85, 85))
# stat <- na.omit(dcast(pl[class=="additive"], pair_name~variable))
# text(2.5, 70, paste0("x", round(mean(stat$log2FC_luc/stat$exp_luc), 1)), pos= 3)
# text(2.5, 70, my_pval_format(t.test(stat$exp_luc, stat$log2FC_luc, paired= T)$p.value)$p.value, pos= 3, offset = 1.4, cex= 0.7)
# stat <- dcast(pl[class=="super-additive"], pair_name~variable)
# text(4.5, 85, paste0("x", round(mean(stat$log2FC_luc/stat$exp_luc), 1)), pos= 3)
# text(4.5, 85, my_pval_format(t.test(stat$exp_luc, stat$log2FC_luc, paired= T)$p.value)$p.value, pos= 3, offset = 1.1)
dev.off()

dat <- readRDS("/groups/stark/vloubiere/projects/pe_STARRSeq/Rdata/B_SCR1_peSTARRSeq_final_table.rds")
luc <- dat[luc, , on= c("enh_L", "enh_R")]
plot(luc$log2FoldChange, luc$log2FC_luc, xlab= "STARR-Seq activity", ylab= "luciferase activity", las= 1)
legend("topleft", paste0("PCC= ", round(cor.test(luc$log2FoldChange, luc$log2FC_luc)$estimate, 2)), bty= "n")
abline(lm(log2FC_luc~log2FoldChange, luc))


