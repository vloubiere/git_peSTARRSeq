
luc <- readRDS("Rdata/validations_luciferase_final_table.rds")

setorderv(luc, "additive")
luc[indL>1 & indR>1, {
  plot(additive,
       log2FoldChange,
       # ylim= range(c(additive, log2FoldChange)),
       xlab= "Predicted additive (5'+3')",
       ylab= "Combined activity")
  .lo <- loess(log2FoldChange~additive)
  lines(.lo$x, .lo$fitted, lty= "11")
  abline(0, 1, lty= "11")
  lm <- lm(log2FoldChange ~ 0 + indL * indR)
  print(summary(lm))
  coefs <- formatC(coef(lm), digits = 2)
  eq <- paste0("y= 0+", coefs[1], "*5'+", coefs[2], "*3'", coefs[3], "*5'*3'")
  legend(par("usr")[1],
         par("usr")[4],
         legend= c("y= predicted additive", eq),
         bty= "n",
         cex= .6,
         text.col= c("black", "red"))
  x <- seq(1, 6, length.out= 100)
  lines(sapply(x, function(x) log2(2^x+2^x-1)),
        predict(lm, data.table(indL= x, indR= x)),
        lty= "11",
        col= "red")
  }]
