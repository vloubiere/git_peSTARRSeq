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
dat[, actClassR:= factor(actClassR, rev(c("Inactive", "Q1", "Q2", "Q3", "Q4")))]

# Motif enrichment
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

#-------------------------#
# Plot
#-------------------------#
pdf("pdf/draft/dev_quantiles_motifs_enrichment_vllib002.pdf",
    width= 6.5,
    height = 5)
par(mar= c(2,20,3,7),
    las= 1)

# Left enrich
pl <- plot(enrL, 
           padj_cutoff= 0.05, 
           top_enrich= 7)
title(main= "Left enhancer", line= 1.5)
vl_add_motifs(pl)

# Motif counts
par(mar= c(3,18,3,7),
    mfrow=c(2,1),
    tcl= -0.2,
    mgp= c(2,0.5,0))
enrL[name %in% c("GATA/4", "AP1/1"), {
  .c <- split(counts[side=="L", variable, with= F][[1]], counts[side=="L", actClass], drop = T)
  .c <- .c[c("Q1", "Q4")]
  br <- seq(0, max(unlist(.c)))
  yl <- c(0, max(sapply(.c, function(x) max(hist(x, breaks= br, plot= F)$density))))
  pl <- function(x, col, ...)
    hist(x, 
         freq = F, 
         breaks= br, 
         ylim= yl, 
         col= col,
         main= name,
         xlab= "Motif counts",
         ...)
  Cc <- c(adjustcolor("tomato", 0.8),
          adjustcolor("grey50", 0.3))
  pl(.c$Q1, Cc[1])
  pl(.c$Q4, Cc[2], add= T)
  legend(par("usr")[2],
         par("usr")[4],
         fill= rev(Cc),
         legend= c("Q1", "Q4"),
         xpd= NA,
         bty= "n")
  .SD
}, .(name, variable)]

# Right enrichment
par(mar= c(2,20,3,7),
    mfrow= c(1,1),
    mgp= c(3,1,0))
pl <- plot(enrR, 
           padj_cutoff= 0.05, 
           top_enrich= 7)
title(main= "Right enhancer", line= 1.5)
vl_add_motifs(pl)
dev.off()