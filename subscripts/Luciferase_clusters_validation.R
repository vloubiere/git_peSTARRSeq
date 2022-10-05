setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
require(data.table)
require(vlfunctions)

#-----------------------------------------------#
# Import data
#-----------------------------------------------#
lib <- readRDS("Rdata/final_results_table.rds")[vllib=="vllib002" & class_act=="enh./enh.", .(L, R, median_L, median_R)]
dat <- readRDS("Rdata/validations_luciferase_final_table.rds")[, .(L, R, log2FoldChange, additive)]
dat <- lib[dat, on= c("L", "R"), nomatch= NULL]
dat[, group:= fcase(grepl("hk", L) | grepl("hk", R), "Hk",
                    median_L>6 | median_R>6, "Strong",
                    default = "Syn")]

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