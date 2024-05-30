setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
require(vlfunctions)

# Import ----
dat <- list(hkCP= readRDS("db/linear_models/FC_RpS12_focused_lm_predictions.rds"),
            devCP= readRDS("db/linear_models/FC_DSCP_focused_lm_predictions.rds"))
dat <- rbindlist(dat, idcol = "CP")
dat <- dat[grepl("^hk|^dev", L) & grepl("^hk|^dev", R)]
dat[, CP:= factor(CP, c("hkCP", "devCP"))]
dat[, classL:= tstrsplit(L, "_", keep= 1)]
dat[, classR:= tstrsplit(R, "_", keep= 1)]
dat[, class:= paste0(CP, "/", classL, "/", classR)]
sel <- c("hkCP/hk/hk", "hkCP/hk/dev", "hkCP/dev/hk", "devCP/dev/dev", "devCP/dev/hk", "devCP/hk/dev")
dat <- dat[class %in% sel]
dat[, class:= factor(class, sel)]
dat[, additive:= log2(2^indL+2^indR-1)]
pl <- melt(dat,
           id.vars = "class",
           measure.vars = c("additive", "log2FoldChange"))
setorderv(pl,
          c("class", "variable"))

# Plot ----
ylim <- c(0, 11)
pdf("pdf/draft/boxplot_hkCP_dev_vs_hk_pairs.pdf",
    width = 3, 
    height = 3)
par(mai= c(0.75,.5, 1, 1), 
    mgp= c(0.75, 0.25, 0),
    cex.lab= 8/12,
    cex.axis= 7/12,
    las= 1,
    tcl= -0.1,
    bty= "n",
    lend= 2,
    lwd= .75,
    xpd= NA)
at <- c(1,1.25,2,2.25,3,3.25,5,5.25,6,6.25,7,7.25)
plot(NA,
     ylim= ylim,
     xlim= c(.75, 7),
     frame= F,
     ylab= "Activity (log2)",
     xaxt= "n",
     xlab= NA)
rect(at[1]-.25,
     par("usr")[3],
     at[6]+.25,
     par("usr")[4]+strheight("M")*1.25,
     border= NA,
     col= adjustcolor("tomato", .2))
text(mean(at[1:6]),
     par("usr")[4]+strheight("M")*1.25,
     "Hk. CP",
     pos= 3,
     cex= 8/12)
rect(at[7]-.25,
     par("usr")[3],
     at[12]+.25,
     par("usr")[4]+strheight("M")*1.25,
     border= NA,
     col= adjustcolor("limegreen", .2))
text(mean(at[7:12]),
     par("usr")[4]+strheight("M")*1.25,
     "Dev. CP",
     pos= 3,
     cex= 8/12)
box <- vl_boxplot(value~variable+class,
                  pl,
                  xaxt= "n",
                  col= "white",
                  at= at,
                  xlab= NA,
                  boxwex= .2,
                  compute.pval= list(c(1,2), c(3,4), c(5,6),
                                     c(2,4), c(2,6),
                                     c(7,8), c(9,10), c(11,12)),
                  tilt.names= T,
                  add= T,
                  pval.cex= 6/12,
                  lwd= .5)$stats
legend(par("usr")[2],
       par("usr")[4],
       legend= c("Observed activity", "Expected additive"),
       density= c(0,100),
       bty= "n",
       cex= 6/12)
vl_tilt_xaxis(apply(matrix(at, ncol= 2, byrow = T), 1, mean),
              labels = gsub("hkCP/|devCP/", "", levels(dat$class)),
              srt = 37.5)
rect(at[1]-.1, box[2,1], at[1]+.1, box[4,1], dens= 100, lwd= .25)
rect(at[3]-.1, box[2,3], at[3]+.1, box[4,3], dens= 100, lwd= .25)
rect(at[5]-.1, box[2,5], at[5]+.1, box[4,5], dens= 100, lwd= .25)
rect(at[7]-.1, box[2,7], at[7]+.1, box[4,7], dens= 100, lwd= .25)
rect(at[9]-.1, box[2,9], at[9]+.1, box[4,9], dens= 100, lwd= .25)
rect(at[11]-.1, box[2,11], at[11]+.1, box[4,11], dens= 100, lwd= .25)
dev.off()