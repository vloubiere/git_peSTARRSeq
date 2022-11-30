setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
require(data.table)
require(vlfunctions)

#-----------------------------------------------#
# Import data
#-----------------------------------------------#
dat <- readRDS("db/FC_tables_DESeq2/vllib002_pe-STARR-Seq_DSCP_T8_SCR1_300_counts_norm_final_oe.rds")
QL <- unique(dat[actClassL!= "inactive", .(L, indL)])[indL>quantile(indL, 0.95), L]
QR <- unique(dat[actClassR!= "inactive", .(R, indR)])[indR>quantile(indR, 0.95), R]
luc <- readRDS("Rdata/validations_luciferase_final_table.rds")
setnames(luc, c("mean_L", "mean_R"), c("indL", "indR"))
model <- readRDS("db/linear_models/vllib002_pe-STARR-Seq_DSCP_T8_SCR1_300_lm.rds")
luc[, predicted:= predict(model, .SD)]
# luc[, predicted:= indL+indR]
dat <- dat[luc, on= c("L", "R"), nomatch= NULL]
dat <- dat[!grepl("^control", L) & !grepl("^control", R)]

# Groups
dat[, group:= fcase(grepl("^hk", L) | grepl("^hk", R), "Hk pair",
                    L %in% QL | R %in% QR, "Saturating pair",
                    default= "Enhancer pair")]

# melt
pl <- melt(dat, 
           id.vars = c("L", "R", "group"),
           measure.vars = c("log2FoldChange", "predicted"))
pl[variable=="log2FoldChange", variable:= "Observed"]
pl[variable=="predicted", variable:= "Predicted"]
pl[, variable:= factor(variable, c("Predicted", "Observed"))]
pl[, group:= factor(group, c("Hk pair", "Enhancer pair", "Saturating pair"))]
Cc <- c("tomato", "limegreen", "purple")
Cc <- adjustcolor(Cc, 0.3)

#-----------------------------------------------#
# Plot
#-----------------------------------------------#
pdf("pdf/draft/Luciferase_clusters_validation.pdf", 
    height = 3, 
    width = 2.5)
par(mgp= c(1.5, 0.5, 0),
    oma= c(0,0,2.25,0),
    mar= c(4,3.5,0,0.25),
    tcl= -0.2,
    las= 1,
    xpd= NA)
box <- vl_boxplot(value~variable+group,
                  pl, 
                  at= rep(seq(1, 4, length.out=3), each= 2)+c(0, 0.6),
                  tilt.names= T,
                  names= function(x) gsub(paste0(".", unique(pl$group), collapse= "|"), "", x),
                  col= rep(Cc, each= 2),
                  compute_pval= list(c(1,2), c(3,4), c(5,6)),
                  ylab= "Luciferase activity (log2)")
box$pval[, {
  .y0 <- unlist(dat0)
  .y1 <- unlist(dat1)
  .x0 <- rep(x0[1], length(.y0))
  .x1 <- rep(x1[1], length(.y1))

  segments(.x0, .y0, .x1, .y1, col= Cc[c(1,3,2)][.GRP], lwd= 1.5)
  # points(.x0, .y0, col= Cc[.GRP], pch= 16, cex= 0.6)
  # points(.x1, .y1, col= Cc[.GRP], pch= 16, cex= 0.6)
}, .(x0, x1)]
legend(par("usr")[1]-strwidth("M"),
       par("usr")[4]+strheight("M")*4.25,
       fill= Cc,
       legend= c(">=1 Hk.", "Dev. pairs", ">=1 sat. enh."),
       bty= "n",
       cex= 0.75,
       border= NA,
       xpd= NA)
dev.off()