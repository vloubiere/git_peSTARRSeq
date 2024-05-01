setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
require(vlfunctions)

dat <- readRDS("Rdata/ham_luciferase_final_table.rds")
indL <- dat[variable=="indL"]
act <- dat[!(variable %in% c("indL", "indR"))]

sp <- c(.2,.2,1,.2,1,.2,1,.2)
bar <- vl_barplot(height = act$mean,
                  sd = act$sd,
                  individual.var = act$all,
                  col= col1 <- c("#D6857F", "orange4",
                                 "#87D873", "tan",
                                 "#87D873", "limegreen",
                                 "#D6857F", "tomato"),
                  compute.diff= list(c(1,2), c(3,4), c(5,6), c(7,8)),
                  space= sp)
axis(1,
     at= (bar[seq(nrow(bar)) %% 2==1]+bar[seq(nrow(bar)) %% 2==0])/2,
     labels = c("A/B", "B/A", "A/A", "B/B"),
     lty= 0)
vl_barplot(height = c(indL[c(1, NA, 2, NA, 3, NA, 4, NA), mean]),
           sd = c(indL[c(1, NA, 2, NA, 3, NA, 4, NA), sd]),
           individual.var = c(indL[c(1, NA, 2, NA, 3, NA, 4, NA), all]),
           col= col1 <- c("#BFFFA4", NA,
                          "#FFBFB8", NA,
                          "#BFFFA4", NA,
                          "#FFBFB8", NA),
           space= sp,
           add= T)
