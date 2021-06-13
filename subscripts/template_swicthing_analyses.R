# Template switching 
ts <- data.table(file= list.files("db/count/", "_umi_counts.txt", full.names = T))
ts[, lib:= tstrsplit(basename(file), "_", keep=1), file]
ts[, cdition:= ifelse(grepl("input", file), "input", "DSCP")]
ts[, rep:= tstrsplit(basename(file), "DSCP_|input_|_umi", keep = 2)]
ts <- ts[, fread(file), ts]

# Barplot percentage ####
res <- ts[, .(perc= round(sum(.SD[spike=="switched", count])/sum(.SD[!is.na(spike), count])*100, 1)), .(lib, cdition, rep)]

pdf("pdf/template_switching_barplot.pdf", height = 6)
par(mar= c(10,4,2,2), las= 2)
bar <- barplot(res$perc, ylim= c(0, 20), ylab= "% template switching", 
               names.arg = res[, paste0(lib, "_", cdition, "_" , rep)])
text(bar, res$perc, round(res$perc, 1), pos= 3)
dev.off()

#---------------------------------------------------------####
# function scatterplots observed vs expected ####
sp <- function(x, y)
{
  lims <- c(0, 20)
  
  layout(matrix(c(2,4,1,3), ncol=2, byrow = T), heights = c(0.25,1), widths = c(1, 0.25))
  par(mar=c(4,4,1,1), pch= 16, las= 1)
  plot(x, y, xlim= lims, ylim= lims, col= adjustcolor("lightgrey", 0.5),
       xlab= "Switched counts", ylab= "Putative correct counts")
  abline(0, 1)
  abline(1, 1, lty= 2)
  abline(-1, 1, lty= 2)
  par(mar=c(1,4,1,1))
  plot(density(x), xlim= lims, xaxt= "n", ylab= "density", main= "")
  .dens <- density(na.omit(y))
  par(mar=c(4,1,1,1))
  plot(.dens$y, .dens$x, type= "l", ylim= lims, yaxt= "n", xlab= "density", main= "")
  plot.new()
}
#---------------------------------------------------------####
# function scatterplots observed vs expected ####
cdition_collapsed <- ts[, .(count= sum(count)), .(lib, cdition, L, R, exist, spike)]
pdf("pdf/template_switching_scatterplots_oe.pdf")
cdition_collapsed[, {
  # Left
  if(lib== "vllib002")
  {
    ok <- .SD[spike== "ok"]
    switched <- .SD[spike== "switched"  & grepl("^ts", L)]
  }else
  {
    ok <- .SD[spike== "ok" & grepl("_C_", L)]
    switched <- .SD[spike== "switched"  & grepl("_C_", L)]
  }
  .c <- merge(switched, ok, by= "L")[!is.na(count.x) & ! is.na(count.y)]
  x <- log2(.c$count.x)
  y <- log2(.c$count.y)
  sp(x, y)
  text(0.5, 0.8, paste0(lib, "_", cdition, "_L"), bty="n", cex= 0.8)
  text(0.5, 0.5, paste0(">2= ", round(length(which(y>(x*2)))/nrow(.c)*100), "%"), bty="n", cex= 0.8)
  text(0.5, 0.2, paste0(">10= ", round(length(which(y>(x*10)))/nrow(.c)*100), "%"), bty="n", cex= 0.8)
  
  # Right
  if(lib== "vllib002")
  {
    ok <- .SD[spike== "ok"]
    switched <- .SD[spike== "switched"  & grepl("^ts", R)]
  }else
  {
    ok <- .SD[spike== "ok" & grepl("_C_", R)]
    switched <- .SD[spike== "switched"  & grepl("_C_", R)]
  }
  .c <- merge(switched, ok, by= "R")[!is.na(count.x) & ! is.na(count.y)]
  x <- log2(.c$count.x)
  y <- log2(.c$count.y)
  sp(x, y)
  text(0.5, 0.8, paste0(lib, "_", cdition, "_R"), bty="n", cex= 0.8)
  text(0.5, 0.5, paste0(">2= ", round(length(which(y>(x*2)))/nrow(.c)*100), "%"), bty="n", cex= 0.8)
  text(0.5, 0.2, paste0(">10= ", round(length(which(y>(x*10)))/nrow(.c)*100), "%"), bty="n", cex= 0.8)
  
  print("")
}, .(lib, cdition)]
dev.off()



