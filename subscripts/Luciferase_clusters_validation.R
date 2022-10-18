setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
require(data.table)
require(vlfunctions)

#-----------------------------------------------#
# Import data
#-----------------------------------------------#
dat <- readRDS("db/FC_tables/vllib002_pe-STARR-Seq_DSCP_T8_SCR1_300_counts_norm_final_oe.rds")
dat[, resGroupL:= cut(meanResidualsL,
                       quantile(unique(meanResidualsL), c(0, 0.2, 0.9, 1)),
                       include.lowest= T,
                       labels= c("Sub-efficient", "Expected", "Over-efficient"))]
dat[, resGroupR:= cut(meanResidualsR,
                       quantile(unique(meanResidualsR), c(0, 0.3, 0.9, 1)),
                       include.lowest= T,
                       labels= c("Sub-efficient", "Expected", "Over-efficient"))]
dat <- dat[!grepl("^control", L) & !grepl("^control", R) & actClassL=="active" & actClassR=="active"]
luc <- readRDS("Rdata/validations_luciferase_final_table.rds")[, .(L, R, log2FoldChange, additive)]
dat <- dat[luc, on= c("L", "R"), nomatch= NULL]
dat[, group:= fcase(resGroupL=="Sub-efficient" | resGroupR=="Sub-efficient", "Sub-efficient",
                    resGroupL=="Over-efficient" | resGroupR=="Over-efficient", "Over-efficient",
                    default = "Expected")]

pl <- melt(dat, 
           id.vars = c("L", "R", "group"),
           measure.vars = c("log2FoldChange", "additive"))
pl[variable=="log2FoldChange", variable:= "Observed"]
pl[variable=="additive", variable:= "Exp. additive"]
pl[, variable:= factor(variable, c("Exp. additive", "Observed"))]
Cc <- c("tomato", "purple", "limegreen")
Cc <- adjustcolor(Cc, 0.3)

#-----------------------------------------------#
# Plot
#-----------------------------------------------#
pdf("pdf/draft/Luciferase_clusters_validation.pdf", 
    height = 3.5, 
    width = 2.25)
par(mgp= c(1.5, 0.5, 0),
    oma= c(0,0,2.25,0),
    mar= c(4.5,3.5,0,0.25),
    tcl= -0.2,
    las= 1,
    xpd= NA)
box <- vl_boxplot(value~variable+group,
                  pl, 
                  at= rep(seq(1, 4, length.out=3), each= 2)+c(0, 0.6),
                  tilt.names= T,
                  names= function(x) gsub(".Hk$|.Strong$|.Syn$|", "", x),
                  col= rep(Cc, each= 2),
                  compute_pval= list(c(1,2), c(3,4), c(5,6)),
                  ylab= "Luciferase activity (log2)")
box$pval[, {
  .y0 <- unlist(dat0)
  .y1 <- unlist(dat1)
  .x0 <- rep(x0[1], length(.y0))
  .x1 <- rep(x1[1], length(.y1))
  
  segments(.x0, .y0, .x1, .y1, col= Cc[.GRP], lwd= 1.5)
  # points(.x0, .y0, col= Cc[.GRP], pch= 16, cex= 0.6)
  # points(.x1, .y1, col= Cc[.GRP], pch= 16, cex= 0.6)
}, .(x0, x1)]
legend(par("usr")[1]-strwidth("M"),
       par("usr")[4]+strheight("M")*4,
       fill= Cc,
       legend= c(">=1 Hk.", ">=1 strong enh.", "Dev. pairs"),
       bty= "n",
       cex= 0.75,
       border= NA,
       xpd= NA)
dev.off()