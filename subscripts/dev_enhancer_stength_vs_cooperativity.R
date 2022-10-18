setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
require(vlfunctions)

# Import data
dat <- readRDS("db/FC_tables/vllib002_pe-STARR-Seq_DSCP_T8_SCR1_300_counts_norm_final_oe.rds")
QL <- quantile(dat[actClassL!= "inactive", .(L, indL)]$indL, seq(0,1, length.out= 5))
dat[actClassL!= "inactive", actClassL:= cut(indL, QL, c("Q1", "Q2", "Q3", "Q4"), include.lowest= T)]
dat[actClassL== "inactive", actClassL:= "Inactive"]
dat[, actClassL:= factor(actClassL, rev(c("Inactive", "Q1", "Q2", "Q3", "Q4")))]
QR <- quantile(dat[actClassR!= "inactive", .(R, indR)]$indR, seq(0,1, length.out= 5))
dat[actClassR!= "inactive", actClassR:= cut(indR, QR, c("Q1", "Q2", "Q3", "Q4"), include.lowest= T)]
dat[actClassR== "inactive", actClassR:= "Inactive"]
dat[, actClassR:= factor(actClassR, c("Inactive", "Q1", "Q2", "Q3", "Q4"))]

mat <- dcast(dat, actClassL~actClassR, value.var = "residuals", fun.aggregate = mean)
mat <- as.matrix(mat, 1)

pdf("pdf/draft/individual_strengths_vs_synergy.pdf", 5.5, 4)
par(las= 1,
    mar= c(3,5,3.5,9),
    xpd= T)
pl <- vl_heatmap(mat, 
                 cluster_rows = F, 
                 cluster_cols = F, 
                 breaks = c(-0.75, 0, 0.75), 
                 show_legend = F)
lw <- diff(grconvertX(c(0, 1), "line", "user"))
lh <- diff(grconvertY(c(0, 1), "line", "user"))
polygon(par("usr")[c(1,2,2)], 
        par("usr")[4]+lh*c(0.5,0.5,1.5),
        border= NA,
        col= "grey20")
text(2.5, par("usr")[4]+lh*1.75, "3' Individual activity")
polygon(par("usr")[2]+lw*c(0.5,0.5,1.5),
        par("usr")[c(3,4,4)], 
        border= NA,
        col= "grey20")
text(par("usr")[2]+lw*1.75, 2.5, "5' Individual activity", srt= -90)
vl_heatkey(breaks = pl$breaks, 
           col= pl$col,
           left = par("usr")[4]+lw*2.5,
           top =  par("usr")[2],
           main = "Mean residuals", 
           height = 1.5)
dev.off()

coor <- unique(dat[, .(actClassL, L)])
seq <- as.data.table(readRDS("Rdata/vl_library_twist008_112019.rds"))
coor[seq, sequence:= i.enh_sequence, on= "L==ID_vl"]

sel <- vl_Dmel_motifs_DB_full[!is.na(FBgn), motif_ID]
counts <- vl_motif_counts(coor$sequence, sel= sel)
counts_list <- split(counts, coor$actClassL, drop = T)
enr <- vl_motif_cl_enrich(counts_list, control_cl = "Inactive")
setorderv(enr, "padj")
enr <- enr[variable %in% enr[, variable[1], name]$V1]

par(mar= c(3,15,2,7),
    las= 1)
pl <- plot(enr, padj_cutoff= 0.05, top_enrich= 6)
vl_add_motifs(pl)




