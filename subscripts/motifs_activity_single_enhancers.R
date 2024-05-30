setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")

# Import twist data ----
dat <- fread("/groups/stark/vloubiere/projects/pe_STARRSeq/Rdata/BA_300bp_TWIST_STARRSeq.txt")
dat <- dat[, .(ID, dev= dev_log2FoldChange, hk= hk_log2FoldChange)]
dat[, c("seqnames", "start", "end"):= tstrsplit(ID, "_", keep= 1:3)]
dat[, seq:= vl_getSequence(.SD, genome = "dm3")]

# MOTIF COUNTS ----
sel <- data.table(name= c("SREBP", "GATA", "AP-1", "Twist", "Trl",  "Dref"),
                  motif_ID= c("cisbp__M0187",
                              "cisbp__M4320",
                              "cisbp__M6317",
                              "flyfactorsurvey__CG16778_SANGER_5_FBgn0003715",
                              "homer__CTCTCTCTCY_GAGA-repeat",
                              "homer__AVYTATCGATAD_DREF"))
counts <- vl_motif_counts(dat$seq,
                          sel = sel$motif_ID,
                          genome = "dm3")
colnames(counts) <- sel$name
dat <- cbind(dat, as.data.table(counts))

# Melt data ----
.m <- melt(dat,
           id.vars = c("ID", colnames(counts)),
           measure.vars = c("hk", "dev"),
           variable.name = "CP")
.m <- melt(.m, measure.vars = colnames(counts), value.name = "count")
.m[, cutoff:= min(c(floor(quantile(count, .995)), 8)), .(CP, variable)]
.m[, count:= ifelse(count>cutoff, cutoff, count), .(count, cutoff)]

# Plot ----
ylab <- "Normalized activity (log2)"

pdf("pdf/draft/motif_counts_individual_act.pdf", width = 2.75, height = 5)
vl_par(mfrow= c(3, 2),
       font.main= 1,
       mai= c(0.75, 0.5, 0.05, 0.05),
       mgp= c(1, .35, 0),
       lwd= .75)
# Heterotypic dev enhancers
for(CP in c("hk", "dev"))
{
  sub <- dat[`AP-1`>0 & GATA==0 & SREBP==0 & Trl==0 & Twist==0]
  adj <- median(sub[[CP]])
  box <- vl_boxplot(dat[`AP-1`==0 & GATA==0 & SREBP==0  & Trl==0 & Twist==0][[CP]]-adj,
                    sub[[CP]]-adj,
                    dat[`AP-1`>0 & GATA>0 & SREBP==0 & Trl==0 & Twist==0][[CP]]-adj,
                    dat[`AP-1`>0 & GATA>0 & SREBP>0 & Trl==0 & Twist==0][[CP]]-adj,
                    dat[`AP-1`>0 & GATA>0 & SREBP>0 & Trl>0 & Twist==0][[CP]]-adj,
                    dat[`AP-1`>0 & GATA>0 & SREBP>0 & Trl>0 & Twist>0][[CP]]-adj,
                    tilt.names = T,
                    xaxt= "n",
                    ylab= ylab,
                    col= adjustcolor(switch(as.character(CP), "hk"= "#EC644D", "dev"= "#58B038"), .3),
                    lwd= .5)
  vl_tilt_xaxis(seq(box$n),
                labels= c("None",
                          "AP-1",
                          "AP-1+GATA",
                          "AP-1+GATA+SREBP",
                          "AP-1+GATA+SREBP+Trl",
                          "AP-1+GATA+SREBP+Trl+Twist"))
  abline(h= 0, lty= "13")
  abline(-2, 1, col= "red")
}
.m[variable %in% c("AP-1", "Dref"), {
  center <- median(value[count==1])
  box <- vl_boxplot(value-center~count,
                    .SD,
                    ylab= ylab,
                    # xaxt= "n",
                    col= adjustcolor(switch(as.character(CP), "hk"= "#EC644D", "dev"= "#58B038"), .3),
                    lwd= .5)
  if(CP=="hk")
    title(ylab= ylab, xpd= NA)
  title(xlab= paste0(variable, " motif counts"))
  abline(h= 0, lty= "13")
  abline(-2, 1, col= "red")
  print(".")
}, .(variable, CP)]
dev.off()