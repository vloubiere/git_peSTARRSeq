setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
require(vlfunctions)

#-----------------------------------------------#
# Import data
#-----------------------------------------------#
lib <- readRDS("Rdata/final_results_table.rds")
feat <- readRDS("Rdata/final_300bp_enhancer_features.rds")

# Import
dat <- lib[vllib=="vllib016"][, diff:= log2FoldChange-additive]
dat <- feat$add_feature(dat, feat$lib)
dat <- feat$add_feature(dat, feat$top_motifs)

# Filtering
dat <- dat[group_L != "repressor" & group_R != "repressor" & median_L<=0]
setorderv(dat, "diff", -1)

# Enriched motifs
vl_boxplot(dat$diff, violin = T)
abline(h= 1, lty= 2)

mot <- readRDS(feat$top_motifs)
cols <- names(mot)[-1]
mot <- as.matrix(mot[match(dat$L, ID), ..cols])+as.matrix(mot[match(dat$R, ID), ..cols])
par(mar= c(3,15,2,5))
enr <- vl_motif_enrich(counts = mot[dat$diff>0,],
                       control_counts = mot[dat$diff<=0,],
                       plot= F)
enr <- enr[padj<0.001 & log2OR>0]
plot(enr)

par(mfrow= c(1,3))
vl_heatmap(as.matrix(dat[, .(median_L, median_R, diff)]),
           cluster_rows= F,
           cluster_cols= F,
           breaks= c(-3,0,3),
           clustering_method = "ward.D2",
           auto_margins = F,
           show_rownames = F)
cols <- enr[c(5,9,11,14), variable]
vl_heatmap(log2(mot[, cols]+1)[dat$diff>1,],
           cluster_rows= F,
           cluster_cols= F, 
           show_rownames = F,
           col= c("blue", "yellow"))
vl_heatmap(log2(mot+1)[dat$diff>1,],
           cluster_rows= F,
           cluster_cols= F, 
           show_rownames = F,
           col= c("blue", "yellow"))

test <- apply(mot, 2, function(x) summary(lm(dat$diff~x)))




# Pairs
P <- as.matrix(dcast(dat, median_L+L~median_R+R, value.var = "diff")[, -1], 1)
P <- P[apply(P, 1, function(x) sum(is.na(x))<70), apply(P, 2, function(x) sum(is.na(x))<70)]
L <- as.matrix(dcast(dat[group=="L"], L~R, value.var = "diff"), 1)
R <- as.matrix(dcast(dat[group=="R"], L~R, value.var = "diff"), 1)

lay <- matrix(9, ncol= 3, nrow=3)
lay[2,1] <- 2
lay[2,2] <- 1
lay[3,2] <- 3
lay[1,2] <- 4
lay[2,3] <- 5
lay[lay==9] <- seq(lay[lay==9])+5

par(mar=c(0.2,0.2,0.2,0.2))
layout(lay, 
       widths = c(1,0.3,0.3),
       heights = c(0.3,0.3,1))
heat <- vl_heatmap(P,
                   cutree_rows = 2,
                   cutree_cols = 2,
                   breaks= c(-3,0,3),
                   clustering_method = "ward.D2",
                   auto_margins = F, 
                   show_row_dendrogram= F, 
                   show_rownames = F, 
                   show_colnames = F)

heatL <- vl_heatmap(L,
                    cluster_cols = F,
                    cluster_rows = F,
                    breaks= c(-3,0,3),
                    show_legend = F, 
                    auto_margins = F, 
                    show_row_dendrogram= F,
                    show_rownames = F, 
                    show_colnames = F)
top_rows <- rownames(L)[rowMeans(L, na.rm= T)>quantile(rowMeans(L, na.rm= T), 0.75)]
abline(h= sum(!rownames(L) %in% top_rows)+0.5)
heatP <- vl_heatmap(P,
                    cluster_cols = F,
                    cluster_rows = F,
                    breaks= c(-3,0,3),
                    auto_margins = F, 
                    show_row_dendrogram= F, 
                    show_rownames = F, 
                    show_colnames = F)
top_cols <- colnames(R)[colMeans(R, na.rm= T)>quantile(colMeans(R, na.rm= T), 0.75)]
abline(h= sum(!rownames(L) %in% top_rows)+0.5)
abline(v= sum(colnames(R) %in% top_cols)+0.5)
heatR <- vl_heatmap(R,
                    cluster_cols = F,
                    cluster_rows = T,
                    breaks= c(-3,0,3),
                    auto_margins = F, 
                    show_row_dendrogram= F,
                    show_rownames = F, 
                    show_colnames = F)
abline(v= sum(colnames(R) %in% top_cols)+0.5)

# Enriched motifs
mot <- readRDS(feat$top_motifs)
cols <- names(mot)[-1]
enrL <- vl_motif_enrich(counts = as.matrix(mot[ID %in% top_rows, ..cols]),
                        control_counts = as.matrix(mot[ID %in% setdiff(rownames(L), top_rows), ..cols]),
                        plot= F)
enrR <- vl_motif_enrich(counts = as.matrix(mot[ID %in% top_cols, ..cols]),
                        control_counts = as.matrix(mot[ID %in% setdiff(colnames(R), top_cols), ..cols]),
                        plot= F)

layout(matrix(c(2,1,4,3),
              ncol= 2, 
              byrow = T), 
       heights = c(1, 0.25),
       widths =  c(0.25, 1))
par(mar= c(0.25,0.25,5,7))
plot(heat, 
     cutree_rows= 2, 
     cutree_cols= 2, 
     show_rownames= F,
     show_colnames= F, 
     breaks= c(-3,0,3), 
     auto_margins= F, 
     legend_title = "Obs./Add. (log2)")
par(mar= c(0.25,7,5,0))
plot(heat, 
     add= log2(mat+1), 
     add_inherit_row_order= T,
     auto_margins= F,
     breaks= c(0,2,4))
vl_heatkey(c(0,4),
           x= -0.5,
           y= 1)
par(mar= c(2,0.25,0,7))
plot(heat, 
     add= log2(matR+1), 
     add_inherit_row_order= F,
     add_inherit_col_order= T,
     auto_margins= F,
     breaks= c(0,2,4))
vl_heatkey(c(0,4), 
           x= -0.5,
           y= 1)

# Motifs pairs
sel <- dat[L %in% rownames(mat) & R %in% colnames(mat) & diff>1]
mot <- melt(sel, 
            id.vars = c("L", "R", "diff"), 
            measure.vars = patterns("top_motif"))
mot[, variable:= gsub("top_motif_|_L$|_R$", "", variable)]
motifs <- data.table(name= unique(mot$variable)) 
motifs[vl_Dmel_motifs_DB_full, motif:= i.motif, on= "name==uniqName_noSpecialChar"]
mot <- mot[, .(value= sum(value)), .(L, R, diff, variable)]
mot <- dcast(mot, 
             L+R+diff~variable, 
             value.var= "value")

# Control pairs
rdm <- vl_control_regions_BSgenome(BSgenome.Dmelanogaster.UCSC.dm3::BSgenome.Dmelanogaster.UCSC.dm3,
                                   5000, 
                                   width = 249*2)
rdm_seq <- getSeq(BSgenome.Dmelanogaster.UCSC.dm3::BSgenome.Dmelanogaster.UCSC.dm3, 
                  rdm$seqnames,
                  rdm$start,
                  rdm$end, 
                  as.character= T)
mot_ctl <- vl_motif_counts(rdm_seq, 
                           sel = motifs$motif)
setnames(mot_ctl, 
         new= motifs[names(mot_ctl), name, on= "motif"])

# Compute enrichment
counts <- rbindlist(list(mot, 
                         mot_ctl),
                    fill= T)

enr <- vl_motif_cl_enrich(as.matrix(counts[, !c("L", "R", "diff")]), 
                          cl_IDs = c(rep("sel", nrow(mot)),
                                     rep("ctl", nrow(mot_ctl))), 
                          control_cl = "ctl")
enr <- na.omit(enr[padj<0.001])
enr[, padj:= -log10(padj)]
enr <- enr[order(log2OR)]
enr[log2OR==Inf, log2OR:= max(enr[is.finite(log2OR), log2OR])]
Cc <- circlize::colorRamp2(range(enr$padj),
                           c("blue", "red"))

#--------------------------#
# PLOT
#--------------------------#
pdf("pdf/analyses/RpS12_inactive_pairs.pdf", 
    width = 8,
    height = 18)
layout(matrix(1:2, nrow= 2), 
       heights = c(1, 1.5))
par(cex= 0.5, 
    mar= c(12,12,3,9),
    cex= 0.4)
vl_heatmap(mat,
           cutree_rows= 3,
           cutree_cols= 2, 
           breaks = c(-2,0,2), 
           auto_margins = F, 
           legend_title = "Obs./Add. (log2)")
par(mar=c(4,10,0,2),
    cex= 0.7)
barplot(enr[, log2OR], 
        horiz = T,
        col= Cc(enr[, padj]),
        border= NA,
        names= enr$motif,
        las= 1,
        xlab= "log2OR")
vl_heatkey(enr$padj, Cc, 0.7, 0.7, main= "padj (-log10)")
dev.off()
