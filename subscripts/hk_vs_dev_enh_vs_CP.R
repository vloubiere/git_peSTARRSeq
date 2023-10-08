setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
require(vlfunctions)

# Import ----
dat <- list(hkCP= readRDS("db/FC_tables/vllib016_DESeq2.rds"),
            devCP= readRDS("db/FC_tables/vllib015_DESeq2.rds"))
dat <- rbindlist(dat,
                 idcol = "CP")
dat[, CP:= factor(CP, c("hkCP", "devCP"))]
dat <- dat[(grepl("^hk", L) & grepl("^hk", R)) | (grepl("^dev", L) & grepl("^dev", R))]
dat[, class:= tstrsplit(L, "_", keep= 1)]
dat[, class:= factor(class, c("hk", "dev"))]
dat[, additive:= log2(2^indL+2^indR-1)]
pl <- melt(dat,
           id.vars = c("CP", "class"),
           measure.vars = c("additive", "log2FoldChange"))


pdf("pdf/draft/boxplot_hkCP_dev_vs_hk_pairs.pdf",
    width = 3, 
    height = 3)
par(mai= c(0.75,0.5,1,1.5), 
    mgp= c(0.75, 0.25, 0),
    cex.lab= 8/12,
    cex.axis= 7/12,
    las= 1,
    tcl= -0.1,
    bty= "n",
    lend= 2,
    lwd= .75,
    pty= "s",
    xpd= NA)
at <- c(1,1.25,2,2.25,3.5,3.75,4.5,4.75)
plot(NA,
     ylim= c(-2, 11),
     xlim= c(.75, 5),
     frame= F,
     ylab= "Activity (log2)",
     xaxt= "n",
     xlab= NA)
rect(at[1]-.25,
     -2,
     at[4]+.25,
     14,
     border= NA,
     col= adjustcolor("tomato", .2))
text(mean(at[1:4]),
     14,
     "Hk. CP",
     pos= 3,
     cex= 8/12)
rect(at[5]-.25,
     -2,
     at[8]+.25,
     14,
     border= NA,
     col= adjustcolor("limegreen", .2))
text(mean(at[5:8]),
     14,
     "Dev. CP",
     pos= 3,
     cex= 8/12)
box <- vl_boxplot(value~variable+class+CP,
                  pl,
                  # notch= T,
                  col= "white",
                  at= at,
                  xaxt= "n",
                  xlab= NA,
                  boxwex= .2,
                  compute_pval= list(c(2,4), c(3,4), c(6,8)),
                  tilt.names= T,
                  add= T,
                  pval_cex= 7/12)$stats
legend(par("usr")[2],
       par("usr")[4],
       legend= c("Observed activity (log2)", "Expected additive (log2)"),
       density= c(0,100),
       bty= "n",
       cex= 6/12)
vl_tilt_xaxis(apply(matrix(at, ncol= 2, byrow = T), 1, mean),
              labels = c("Hk. / Hk.", "Dev. / Dev.", "Hk. / Hk.", "Dev. / Dev."),
              srt = 37.5)
rect(at[1]-.1, box[2,1], at[1]+.1, box[4,1], dens= 100, lwd= .25)
rect(at[3]-.1, box[2,3], at[3]+.1, box[4,3], dens= 100, lwd= .25)
rect(at[5]-.1, box[2,5], at[5]+.1, box[4,5], dens= 100, lwd= .25)
rect(at[7]-.1, box[2,7], at[7]+.1, box[4,7], dens= 100, lwd= .25)
dev.off()
