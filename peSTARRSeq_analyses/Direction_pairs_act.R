setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
sapply(list.files("/groups/stark/vloubiere/functions/", ".R$", full.names = T), source)
require(data.table)
options(datatable.print.topn=1)
require(yarrr)

feat <- as.data.table(readRDS("Rdata/library/vl_library_112019.rds"))
feat[group %in% c("heatshock", "ecdysone"), group:= "inducible"]
dat <- readRDS("Rdata/processed_peSTARRSeq_data/DESeq2_FC_table.rds")
dat <- dat[enh_L!=enh_R]
dat[feat, group_L:= i.group, on= "enh_L==ID_vl"]
dat[feat, group_R:= i.group, on= "enh_R==ID_vl"]

#------------------------------------------------------------#
# 1- All pairs
#------------------------------------------------------------#
pa <- as.matrix(dcast(dat, enh_L~enh_R, value.var = "log2FoldChange"), 1)
pa[lower.tri(pa)] <- NA
pa <- na.omit(as.data.table(reshape2::melt(pa)))
colnames(pa) <- c("enh_L", "enh_R", "log2FoldChange")
pa[dat, rev:= i.log2FoldChange, on= c("enh_L==enh_R", "enh_R==enh_L")]
pa[dat, c("group_L", "group_R"):= .(i.group_L, i.group_R), on= c("enh_L", "enh_R")]
pa <- na.omit(pa)
.lm.pa <- lm(rev~log2FoldChange, pa)

pa.m <- melt(pa, id.vars = c("enh_L", "enh_R", "group_L", "group_R"))
pa.m[, cdition:= paste(group_L, "*", group_R, ":", ifelse(variable=="log2FoldChange", "X~Y", "Y~X"))]
par(mar=c(10,5,2,2))
box <- my_boxplot(formula= value~cdition, data= pa.m, plot= F)$DT_plot
at.all <- box[, rep(median(unlist(var)), 2), gsub(" : Y~X| : X~Y", "", .id)]
at.all <- order(order(at.all$V1, decreasing = T))

#------------------------------------------------------------#
# 2- Control pairs
#------------------------------------------------------------#
ctl <- CJ(control= grep("^control", unique(c(dat$enh_L, dat$enh_R)), value = T), 
          enh= unique(c(dat$enh_L, dat$enh_R), invert = T, value = T))
ctl[dat, right:= i.log2FoldChange, on= c("control==enh_L", "enh==enh_R")]
ctl[dat, left:= i.log2FoldChange, on= c("control==enh_R", "enh==enh_L")]
ctl[feat, group:= i.group, on= "enh==ID_vl"]
ctl <- na.omit(ctl[enh!=control])
.lm.ctl <- lm(ctl$right~ctl$left, ctl)

ctl.m <- melt(ctl, id.vars = c("control", "enh", "group"))
ctl.m[, cdition:= paste(group, variable)]
box <- my_boxplot(formula= value~cdition, data= ctl.m, plot= F)$DT_plot
box[, .id:= gsub(" left| right", "", .id)]
at.ctl <- box[, rep(median(unlist(var)), 2), .id]
at.ctl <- order(order(at.ctl$V1, decreasing = T))

#------------------------------------------------------------#
# 3- PLOT
#------------------------------------------------------------#
pdf("pdf/peSTARRSeq/Direction_impact.pdf", width = 11, height = 9)
layout(matrix(c(1,2,2,3,4,5), ncol = 3, byrow = T))
par(mar= c(10,5,2,2))

smoothScatter(pa[, .(log2FoldChange, rev)], las= 1, xlab= "left~right activity (log2)", ylab= "right~left activity (log2)", bandwidth = 0.1)
abline(.lm.pa, lty= 2)
PCC <- paste("PCC=", round(cor.test(pa$log2FoldChange, pa$rev)$estimate, 2))
rsq <- paste("R²=", round(summary(.lm.pa)$r.squared, 2))
legend("topleft", c(PCC, rsq), bty= "n")
my_fig_label("A", cex= 2)
my_boxplot(formula= value~cdition, data= pa.m, remove_empty = T, las= 2, at= at.all,
           pval_list = lapply(seq(1, length(at.all), 2), function(x) c(x, x+1)), ylim= c(-3, 9), ylab = "activity (log2)")
my_fig_label("B", cex= 2)

smoothScatter(ctl[, .(left, right)], las= 1, xlab= "candidate~control activity (log2)", ylab= "control~candidate activity (log2)", bandwidth = 0.1)
abline(.lm.ctl, lty= 2)
PCC <- paste("PCC=", round(cor.test(ctl$right, ctl$left)$estimate, 2))
rsq <- paste("R²=", round(summary(.lm.ctl)$r.squared, 2))
legend("topleft", c(PCC, rsq), bty= "n")
my_fig_label("C", cex= 2)
Cc <- c("#74C27A", "tomato", "royalblue2", "gold", "lightgrey")
my_boxplot(formula= value~cdition, data= ctl.m, remove_empty = T, las= 2, at= at.ctl, col_box = rep(Cc, each=2),
           pval_list = lapply(seq(1, length(at.ctl), 2), function(x) c(x, x+1)), ylim= c(-3, 8), ylab = "activity (log2)")
my_fig_label("D", cex= 2)
my_boxplot(formula= right-left~group, data= ctl, remove_empty = T, las= 2, col_box = Cc, yalb= "right-left activity (log2)")
abline(h= 0, lty= 2)
my_fig_label("E", cex= 2)
dev.off()


pdf("pdf/peSTARRSeq/linkers_impact_direction_effect.pdf", height = 5)
sub <- ctl.m[group=="dev"]
sub[variable=="right", linker:= paste0(substr(control, nchar(control)-6, nchar(control)-6), substr(enh, nchar(enh)-6, nchar(enh)-6))]
sub[variable=="left", linker:= paste0(substr(enh, nchar(enh)-6, nchar(enh)-6), substr(control, nchar(control)-6, nchar(control)-6))]
my_boxplot(value~variable+linker, sub, las= 2, pval_list = lapply(seq(1, 17, 2), function(x) c(x, x+1)), ylim= c(-3,8), 
           col_box = c(rep(c("lightgrey", "lightgrey", "black", "black"), 4), "lightgrey", "lightgrey"), ylab = "activity (log2)")
dev.off()








