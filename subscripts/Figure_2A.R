setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
require(data.table)
require(vlfunctions)

#-----------------------------------------------#
# Import data
#-----------------------------------------------#
if(!exists("lib"))
  lib <- readRDS("Rdata/final_results_table.rds")
dat <- lib[vllib=="vllib002" & class_act== "enh./enh."]
dat <- dat[sample(nrow(dat), nrow(dat))]

# Train CV linear model
set.seed(1)
sel_L <- dat[, 1-.N/dat[,.N], L][, sample(L, round(.N/10), prob = V1)]
set.seed(1)
sel_R <- dat[, 1-.N/dat[,.N], R][, sample(R, round(.N/10), prob = V1)]
dat[, set:= ifelse(L %in% sel_L | R %in% sel_R, "test", "train")]
# Train linear model
model <- lm(formula = log2FoldChange~median_L*median_R, 
            data= dat[set=="train"])
dat[, predicted:= predict(model, newdata = dat)]
rsq <- vl_model_eval(observed = dat[set=="test", log2FoldChange], 
                     predicted = predict(model, new= dat[set=="test"]))
eq <- vl_model_equation(model, digits = 2)
eq <- gsub("median_L", "5'", eq)
eq <- gsub("median_R", "3'", eq)
eq <- gsub("log2FoldChange ", "Act.", eq)
eq <- gsub("\\s\\*\\s", "*", eq)

# Select examples
# sel <- dat[between(predicted, 4, 8), .(L, R, log2FoldChange, predicted, median_L, median_R)]
# sel[, diff:= log2FoldChange-predicted]
# sel <- sel[abs(diff)>=1.5 | between(diff, -0.2, 0.2)]
# sel <- sel[abs(median_L-median_R)<0.5]
# sel[, cut:= cut(predicted, seq(4, 8, length.out= 20), include.lowest= T)]
# sel <- sel[, check:= .N>=2 & any(diff>1.5) & any(diff<(-1.5)), .(cut, L)][(check)]
# sel <- sel[order(L, diff, decreasing = T)]
ex <- dat[list(c("dev_strong_B_00302", "dev_strong_B_00302","dev_strong_B_00302"),
               c("dev_medium_C_00480", "dev_medium_C_00466", "dev_medium_B_00585")), on= c("L", "R")]
ex[, col:= c("tomato", "limegreen", "cornflowerblue")]

# Examples barplot function
ex_barplot <- function(dat, L, R, col)
{
  par(lwd= 0.1)
  .c <- SJ(L= L, R= R)
  for(cdition in c("add", "mult"))
  {
    if(cdition=="add")
    {
      current <- dat[.c, .(Exp= 2^additive, 
                           Activity= 2^log2FoldChange,
                           ind_L= 2^median_L), on= c("L", "R")]
      yl <- "Activity"
      names <- c("Additive", "Observed")
    }else if (cdition=="mult")
    {
      current <- dat[.c, .(Exp= multiplicative, 
                           Activity= log2FoldChange,
                           ind_L= median_L), on= c("L", "R")]
      yl <- "Activity (log2)"
      names <- c("Multiplicative", "Observed")
    }
    bar <- barplot(unlist(current[, .(Exp, Activity)]), 
                   col= c(col, "lightgrey"), 
                   xaxt= "n",
                   ylab= yl,
                   space= 1,
                   xlim= c(0.25,4))
    barplot(current[, ind_L], 
            col= "black",
            density= 40,
            add= T,
            space= 1)
    text(bar, 
         par("usr")[3]-strheight("M", cex= 0.5), 
         names,
         pos= 2, 
         offset= -0.25, 
         xpd= T, 
         srt= 45, 
         cex= 0.8)
  }
  par(lwd= 1)
}

# Plot
Cc <- colorRampPalette(c("white", "grey90", "grey60", "grey20"))

pdf("pdf/draft/Figure_2A.pdf",
    height = 6, 
    width = 6)
layout(matrix(c(1,2,2,2,
                3,4,5,10,
                3,6,7,10,
                3,8,9,10), 
              ncol= 4, 
              byrow = T),
       widths = c(1,0.2,0.2,0.6),
       heights =  c(1,1/3,1/3,1/3))
par(las= 1,
    mar= c(3,4,2,0.5),
    mgp= c(1.75, 0.5, 0),
    tcl= -0.2)
smoothScatter(dat[, .(`Expected additive (log2)`= additive, 
                      `Activity (log2)`= log2FoldChange)],
              main= "Additive",
              colramp = Cc)
abline(0,1,lty=2)
smoothScatter(dat[, .(`Expected multiplicative (log2)`= multiplicative, 
                      `Activity (log2)`=log2FoldChange)],
              main= "Multiplicative",
              colramp = Cc)
abline(0,1,lty=2)
smoothScatter(dat[, .(`Predicted (log2)`= predicted, 
                      `Activity (log2)`= log2FoldChange)],
              main= "Linear model",
              colramp = Cc)
abline(0,1,lty=2)
ex[, points(predicted, 
            log2FoldChange, 
            col= col, 
            pch= 16), .(L, R, col)]
legend("topleft", 
       c(eq, 
         paste0("R²= ", round(rsq$Rsquare, 2))),
       bty= "n")

par(mar= c(3,2.6,0.5,0.1),
    mgp= c(1.4, 0.5, 0),
    cex.axis= 0.6, 
    cex.lab= 0.7)
options(scipen= 1)
ex[, ex_barplot(ex, L, R, col), .(L, R, col)]

par(mar= c(3,4,2,0.5),
    mgp= c(1.75, 0.5, 0),
    cex.axis= 1, 
    cex.lab= 1)
vl_boxplot(dat[, .(log2FoldChange-additive,
                   log2FoldChange-multiplicative,
                   log2FoldChange-predicted)],
           ylab= "Residuals (log2)",
           xaxt= "n")
text(1:3, 
     par("usr")[3]-strheight("M", cex= 0.5), 
     c("Additive", "Multiplicative", "Linear model"),
     pos= 2, 
     offset= -0.25, 
     xpd= T, 
     srt= 45, 
     cex= 0.8)
abline(h= 0, lty= 2)
dev.off()