setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
require(vlfunctions)

# motif counts focused library ----
if(!file.exists("db/motif_counts/twist012_motif_counts_Dref_AP1_GATA.rds"))
{
  lib <- readRDS("Rdata/vl_library_twist12_210610.rds")
  sel <- data.table(name= c("GATA", "AP-1", "Dref", "SREBP"),
                    motif_ID= c("cisbp__M4320",
                                "cisbp__M6317",
                                "homer__AVYTATCGATAD_DREF",
                                "cisbp__M2388"))
  counts <- vl_motif_counts(lib$enh_seq,
                            sel$motif_ID,
                            genome= "dm3",
                            p.cutoff = 1e-04)
  names(counts) <- sel$name
  counts[, ID:= lib$ID]
  saveRDS(counts,
          "db/motif_counts/twist012_motif_counts_Dref_AP1_GATA.rds")
}else
  counts <- readRDS("db/motif_counts/twist012_motif_counts_Dref_AP1_GATA.rds")

# Import data ----
dat <- readRDS("db/linear_models/FC_RpS12_focused_lm_predictions.rds")
dat <- dat[grepl("^hk|^shared|^dev", L) & grepl("^hk|^shared|^dev", R)]
hk <- counts[Dref>=1 & `AP-1`+GATA+SREBP==0, ID]
dev1 <- counts[`AP-1`+GATA+SREBP>=1, ID]
dev2 <- counts[`AP-1`+GATA+SREBP>=2, ID]
dev3 <- counts[`AP-1`+GATA+SREBP>=3, ID]

vl_boxplot(dat[L %in% hk & R %in% hk, `Additive model`],
           dat[L %in% hk & R %in% hk, log2FoldChange],
           dat[L %in% hk & R %in% dev1, `Additive model`],
           dat[L %in% hk & R %in% dev1, log2FoldChange],
           dat[L %in% hk & R %in% dev2, `Additive model`],
           dat[L %in% hk & R %in% dev2, log2FoldChange],
           dat[L %in% hk & R %in% dev3, `Additive model`],
           dat[L %in% hk & R %in% dev3, log2FoldChange],
           col= c("white", "grey"),
           compute.pval= list(c(1,2), c(3,4), c(5,6)),
           notch= T)



sel <- data.table(name= c("GATA", "AP-1", "Twist", "Trl",  "Dref", "Dref"),
                  motif_ID= c("cisbp__M4320", "cisbp__M6317",
                              "flyfactorsurvey__CG16778_SANGER_5_FBgn0003715",
                              "homer__CTCTCTCTCY_GAGA-repeat",
                              "flyfactorsurvey__Dref_FlyReg_FBgn0015664",
                              "homer__AVYTATCGATAD_DREF"),
                  col= rep(c("limegreen", "tomato", "cornflowerblue"), each= 2))


dat[, class:= fcase(grepl("^hk", L) & grepl("^hk", R), "hk/hk",
                    grepl("^shared", L) & grepl("^shared", R), "shared/shared",
                    grepl("^hk", L) & grepl("^shared", R), "hk/shared",
                    grepl("^shared", L) & grepl("^hk", R), "shared/hk")]
dat <- dat[!is.na(class)]

vl_boxplot(log2FoldChange~class,
           dat)

dat <- dat[grepl]
dat <- rbindlist(dat, idcol = "CP")
dat <- dat[grepl("^hk|^dev", L) & grepl("^hk|^dev", R)]
dat[, CP:= factor(CP, c("hkCP", "devCP"))]
dat[, classL:= tstrsplit(L, "_", keep= 1)]
dat[, classR:= tstrsplit(R, "_", keep= 1)]
dat[, class:= paste0(classL, "/", classR)]
dat <- dat[!((CP=="hkCP" & class=="dev/dev") | (CP=="devCP" & class=="hk/hk"))]
dat[, class:= factor(class, c("hk/hk", "dev/dev", "dev/hk", "hk/dev"))]
dat[, additive:= log2(2^indL+2^indR-1)]
pl <- melt(dat,
           id.vars = c("CP", "class"),
           measure.vars = c("additive", "log2FoldChange"))
setorderv(pl,
          c("CP", "class", "variable"))
pl <- pl[, .(value= .(value)), .(CP, class, variable)]

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
box <- vl_boxplot(pl$value,
                  # notch= T,
                  col= "white",
                  at= at,
                  xaxt= "n",
                  xlab= NA,
                  boxwex= .2,
                  compute.pval= list(c(1,2), c(3,4), c(5,6),
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
              labels = c("Hk. / Hk.", "Dev. / Dev.", "Hk. / Hk.", "Dev. / Dev."),
              srt = 37.5)
rect(at[1]-.1, box[2,1], at[1]+.1, box[4,1], dens= 100, lwd= .25)
rect(at[3]-.1, box[2,3], at[3]+.1, box[4,3], dens= 100, lwd= .25)
rect(at[5]-.1, box[2,5], at[5]+.1, box[4,5], dens= 100, lwd= .25)
rect(at[7]-.1, box[2,7], at[7]+.1, box[4,7], dens= 100, lwd= .25)
rect(at[9]-.1, box[2,9], at[9]+.1, box[4,9], dens= 100, lwd= .25)
rect(at[11]-.1, box[2,11], at[11]+.1, box[4,11], dens= 100, lwd= .25)
dev.off()