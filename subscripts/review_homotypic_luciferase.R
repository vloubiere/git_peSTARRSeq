setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
require(vlfunctions)

addFC <- function(value,
                  bar,
                  obsIdx,
                  predIdx,
                  adj= c(1, 2, 1),
                  horiz= F)
{
  width <- if(horiz)
    strwidth("M")*0.25 else
      strheight("M")*0.25
  x <- rep(bar[c(obsIdx, predIdx)], each= 2)
  y <- c(value[obsIdx]+width*adj[1],
         max(value[c(predIdx, obsIdx)])+width*adj[2], 
         max(value[c(predIdx, obsIdx)])+width*adj[2], 
         value[predIdx]+width*adj[3])
  if(horiz)
  {
    y <- x
    x <- y
  }
  lines(x,
        y,
        xpd= T,
        lwd= .75)
  x <- mean(bar[c(obsIdx, predIdx)])
  y <- max(value[c(predIdx, obsIdx)])+width*adj[2]
  if(horiz)
  {
    y <- x
    x <- y
  }
  text(x,
       y,
       pos= ifelse(horiz, 4, 3),
       paste0("x", round(value[predIdx]/value[obsIdx], 1)),
       cex= 5/12,
       offset= .1,
       xpd= T)
}

dat <- readRDS("Rdata/ham_luciferase_final_table.rds")

col1 <- c("#BFFFA4", "#87D873", "#87D873", "limegreen",
          "#FFBFB8", "#D6857F", "#D6857F", "tomato",
          "#D6857F", "tan", "#87D873", "orange4")
col2 <- c(NA, NA, "#BFFFA4", NA,
          NA, NA, "#FFBFB8", NA,
          "#BFFFA4", NA, "#FFBFB8", NA)

var1 <- c(dat[L=="A" & R =="control", log2FoldChange],
          dat[L=="control" & R =="A", log2FoldChange],
          dat[L=="A" & R =="A", log2(2^indL+2^indR-1)],
          dat[L=="A" & R =="A", log2FoldChange],
          dat[L=="B" & R =="control", log2FoldChange],
          dat[L=="control" & R =="B", log2FoldChange],
          dat[L=="B" & R =="B", log2(2^indL+2^indR-1)],
          dat[L=="B" & R =="B", log2FoldChange],
          dat[L=="A" & R =="B", log2(2^indL+2^indR-1)],
          dat[L=="A" & R =="B", log2FoldChange],
          dat[L=="B" & R =="A", log2(2^indL+2^indR-1)],
          dat[L=="B" & R =="A", log2FoldChange])
space <- c(.25,.25,.25,.25,1,
           .25,.25,.25,1,
           .25,1,
           .25)

pdf("pdf/draft/review_homotypic_pairs.pdf", 3.75, 3)
vl_par(lwd= .25)
bar <- barplot(2^var1,
               col= col1,
               ylab= "Normalized luciferase activity",
               space = space)
vl_tilt_xaxis(bar, 
              labels = c("A/Ctl.", "Ctl./A", "A/A exp. additive", "A/A",
                         "B/Ctl.", "Ctl./B", "B/B exp. additive", "B/B",
                         "A/B exp. additive", "A/B", "B/A exp. additive", "B/A"))
barplot(2^c(NA,NA,dat[L=="A" & R =="control", log2FoldChange],
            NA,NA,NA,dat[L=="B" & R =="control", log2FoldChange],
            NA,dat[L=="A" & R =="control", log2FoldChange],
            NA,dat[L=="B" & R =="control", log2FoldChange],
            NA),
        col= col2,
        ylab= "Normalized luciferase activity",
        space= space,
        add= T)
addFC(2^var1, bar, obsIdx = 3, predIdx = 4)
addFC(2^var1, bar, obsIdx = 7, predIdx = 8)
addFC(2^var1, bar, obsIdx = 9, predIdx = 10)
addFC(2^var1, bar, obsIdx = 11, predIdx = 12)
dev.off()