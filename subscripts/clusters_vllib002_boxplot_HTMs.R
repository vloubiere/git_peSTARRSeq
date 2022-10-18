setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
require(vlfunctions)

# Import data
dat <- readRDS("db/FC_tables/vllib002_pe-STARR-Seq_DSCP_T8_SCR1_300_counts_norm_final_oe.rds")[!grepl("^control", L) & !grepl("^control", R)]

pl <- rbind(unique(dat[, .(side= "3'", ID= L, gp= resGroupL, vl_toDTranges(coorL))]),
            unique(dat[, .(side= "5'", ID= R, gp= resGroupR, vl_toDTranges(coorR))]))
files <- c("../available_data_dm3/db/bw/ATAC_merged.bw",
           "/groups/stark/haberle/data/public/dm3/S2_cells/ChIPseq/GSE41440/tracks/S2_H3K27ac_ChIP_Rep1.bw",
           "/groups/stark/haberle/data/public/dm3/S2_cells/ChIPseq/GSE41440/tracks/S2_H3K27ac_ChIP_Rep2.bw",
           "/groups/stark/haberle/data/public/dm3/S2_cells/ChIPseq/GSE41440/tracks/S2_H3K4me1_ChIP_Rep1.bw",
           "/groups/stark/haberle/data/public/dm3/S2_cells/ChIPseq/GSE41440/tracks/S2_H3K4me1_ChIP_Rep2.bw", 
           "/groups/stark/haberle/data/public/dm3/S2_cells/ChIPseq/GSE81795/tracks/S2_H3K4me3_Rep1.bw",
           "/groups/stark/haberle/data/public/dm3/S2_cells/ChIPseq/GSE81795/tracks/S2_H3K4me3_Rep2.bw")
names <- c("ATAC", "H3K27Ac.1", "H3K27Ac.2", "H3K4me1.1", "H3K4me1.2", "H3K4me3.1", "H3K4me3.2")
pl[, (names):= lapply(files, function(x) vl_bw_coverage(pl, x))]
pl[, H3K27Ac:= H3K27Ac.1+H3K27Ac.2]
pl[, H3K4me1:= H3K4me1.1+H3K4me1.2]
pl[, H3K4me3:= H3K4me3.1+H3K4me3.2]
cols <- c("ATAC", "H3K27Ac", "H3K4me1", "H3K4me3")

pdf("pdf/draft/clusters_vllib002_boxplot_HTMs.pdf",
    width= 3.8,
    height= 2.5)
par(mfrow=c(1,4),
    mgp= c(2, 0.5, 0),
    oma= c(0,3,0,0),
    mar= c(4.75,1,1.5,0.5),
    tcl= -0.2,
    las= 1,
    xpd= NA)
pl[, {
  lapply(seq(.SD), function(i) 
    vl_boxplot(.SD[[i]]~gp, 
               tilt.names= T,
               main= names(.SD)[i], 
               compute_pval= list(c(1,2), c(2,3), c(1,3)),
               ylab= ifelse(i==1, "Enrichment", NA),
               col= grey.colors(4, 0.9, 0.5)[i]))
  ""
}, side, .SDcols= cols]
dev.off()