setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
require(vlfunctions)

#-----------------------------------#
# Import data and compute quantiles
#-----------------------------------#
dat <- readRDS("db/FC_tables/vllib002_pe-STARR-Seq_DSCP_T8_SCR1_300_counts_norm_final_oe.rds")
QL <- quantile(dat[actClassL!= "inactive", .(L, indL)]$indL, seq(0,1, length.out= 6))
dat[, actClassL:= cut(indL, 
                      c(-Inf, QL), 
                      c("Inactive", paste0("Q", 1:5)), include.lowest= T)]
dat[, actClassL:= factor(actClassL, rev(levels(actClassL)))]
QR <- quantile(dat[actClassR!= "inactive", .(R, indR)]$indR, seq(0,1, length.out= 6))
dat[, actClassR:= cut(indR, 
                      c(-Inf, QR), 
                      c("Inactive", paste0("Q", 1:5)), include.lowest= T)]

#-----------------------------------#
# Matrix heatmap
#-----------------------------------#
mat <- dcast(dat, actClassL~actClassR, value.var = "residuals", fun.aggregate = mean)
mat <- as.matrix(mat, 1)

#-----------------------------------#
# Motif enrichment
#-----------------------------------#
lib <- as.data.table(readRDS("Rdata/vl_library_twist008_112019.rds"))
seq <- rbind(unique(dat[, .(actClass= actClassL, ID= L, side= "L")]),
             unique(dat[, .(actClass= actClassR, ID= R, side= "R")]))
seq[lib, sequence:= i.enh_sequence, on= "ID==ID_vl"]
sel <- vl_Dmel_motifs_DB_full[!is.na(FBgn), motif_ID]
counts <- cbind(seq, vl_motif_counts(seq$sequence, sel= sel))
# Left
enrL <- vl_motif_cl_enrich(split(counts[side=="L", ..sel], counts[side=="L", actClass], drop = T), 
                           control_cl = "Inactive")
setorderv(enrL, "padj")
enrL <- enrL[variable %in% enrL[, variable[1], name]$V1] # Select top enrichment
# Right
enrR <- vl_motif_cl_enrich(split(counts[side=="R", ..sel], counts[side=="R", actClass], drop = T), 
                           control_cl = "Inactive")
setorderv(enrR, "padj")
enrR <- enrR[variable %in% enrR[, variable[1], name]$V1] # Select top enrichment

#-----------------------------------#
# PLOT
#-----------------------------------#
pdf("pdf/draft/individual_strengths_vs_synergy.pdf", 
    width= 6.5,
    height = 5)
par(las= 1,
    mar= c(6,9,5.5,10),
    xpd= T)
# Quantiles residuals heatmap ---------#
pl <- vl_heatmap(mat, 
                 cluster_rows = F, 
                 cluster_cols = F, 
                 breaks = c(-0.75, 0, 0.75), 
                 show_legend = F, 
                 show_colnames = F)
axis(1, at= seq(ncol(mat)), labels= colnames(mat), lty= 0, gap.axis = 0)
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

# Motif enrichment --------------------#
par(mar= c(3,20,3,6),
    las= 1)
# Left
pl <- plot(enrL, 
           padj_cutoff= 0.05, 
           top_enrich= 6)
title(main= "Left enhancer", line= 1.5)
vl_add_motifs(pl)
# Right
pl <- plot(enrR, 
           padj_cutoff= 0.05, 
           top_enrich= 6)
title(main= "Right enhancer", line= 1.5)
vl_add_motifs(pl)

# Motif counts and variety ------------#
par(mar= c(3,5,3,5),
    mgp= c(2,0.5,0),
    tcl= -0.2,
    mfrow= c(2,2))
enrL[name %in% c("GATA/4", "AP1/1"), {
  .c <- split(counts[side=="L", variable, with= F][[1]], counts[side=="L", actClass], drop = T)
  .c <- .c[c("Q1", "Q5")]
  br <- seq(0, max(unlist(.c)))
  yl <- c(0, max(sapply(.c, function(x) max(hist(x, breaks= br, plot= F)$density))))
  pl <- function(x, col, ...)
    hist(x, 
         freq = F, 
         breaks= br, 
         ylim= yl, 
         col= col,
         main= paste0("5'", name, "counts"),
         xlab= "Motif counts",
         ...)
  Cc <- c("grey80",
          adjustcolor("tomato", 0.3))
  pl(.c$Q1, Cc[1])
  pl(.c$Q5, Cc[2], add= T)
  legend("topright",
         fill= Cc,
         legend= c("Q1", "Q5"),
         xpd= NA,
         bty= "n")
  .SD
}, .(name, variable)]
barplot(t(as.matrix(enrL[padj<0.05, .N, keyby= cl][.N:1], 1)),
        ylab= "Number of motif clusters (padj<0.05)")
dev.off()
