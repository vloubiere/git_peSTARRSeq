setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
sapply(list.files("/groups/stark/vloubiere/functions/", ".R$", full.names = T), source)
require(data.table)

dat <- readRDS("Rdata/SCR1_peSTARRSeq_final_table.rds")
dat <- dat[grepl("dev", enh_L) & grepl("dev", enh_R) & median_L>2 & median_R>2]
dat <- dat[enh_L!=enh_R]
feat <- readRDS("Rdata/lib_features.rds")
feat <- melt(feat, id.vars = "ID", measure.vars = patterns("^motif__"))
feat <- feat[ID %in% dat$enh_L | ID %in% dat$enh_R]
setkeyv(feat, c("variable", "ID"))

#####
combinations <- CJ(mot_L= unique(feat$variable), mot_R= unique(feat$variable))
combinations <- combinations[,
                             {
                               current <- data.table(feat[.(mot_L, dat$enh_L), .(enh_L= ID, counts_L= value), allow.cartesian= T],
                                                     feat[.(mot_R, dat$enh_R), .(enh_R= ID, counts_R= value), allow.cartesian= T])
                             }, .(mot_L, mot_R)]
combinations <- combinations[, c("median_L", "median_R", "diff") := dat[.BY, .(median_L, median_R, diff), on= c("enh_L", "enh_R")], .(enh_L, enh_R)]

res <- combinations[,
                    {
                      res1 <- as.data.table(summary(lm(diff~median_L+median_R+counts_L+counts_R, .SD))$coefficients)
                      res2 <- as.data.table(summary(lm(diff~median_L+median_R+counts_L+counts_R+counts_L*counts_R, .SD))$coefficients)
                      res <- data.table(rn= c("mot_L", "mot_R", "mot_LR"),
                                        tval= c(res1$`t value`[4:5], res2$`t value`[6]),
                                        pval= c(res1$`Pr(>|t|)`[4:5], res2$`Pr(>|t|)`[6]))
                    }, .(mot_L, mot_R)]
res[, padj:= p.adjust(pval), rn]
res[, c("mot_L", "mot_R"):= lapply(.SD, function(x) gsub("^motif__", "", x)), .SDcols= c("mot_L", "mot_R")]

cl_info <- readRDS("Rdata/som_enriched_motifs.rds")$info
cl_info[Dmel=="", Dmel:= best_match]
mat <- dcast(res[rn=="mot_LR"], mot_L~mot_R, value.var = "tval")
# mat <- dcast(res[rn=="mot_R"], mot_L~mot_R, value.var = "tval")
# mat <- dcast(res[rn=="mot_L"], mot_L~mot_R, value.var = "tval")
colnames(mat)[-1] <- cl_info$Dmel[match(colnames(mat)[-1], cl_info$best_match)]
mat$mot_L <- cl_info$Dmel[match(mat$mot_L, cl_info$best_match)]

Cc <- colorRampPalette(c("cornflowerblue", "white", "tomato"))(100)
par(mar= c(20,20,5,7))
mat <- as.matrix(mat, 1)
my_pheatmap(mat, lim = c(-5, 5), col= Cc, col_labels = T, row_labels = T)


