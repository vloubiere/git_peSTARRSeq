setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")

# Import data ----
dat <- readRDS("db/FC_tables/DSCP_large_WT_FC_DESeq2.rds")
dat <- dat[(grepl("^dev", L) | ctlL) & (grepl("^dev", R) | ctlR)]
dat[, resL:= mean(log2FoldChange-log2(2^indL+2^indR-1), na.rm= T), L]
dat[, resR:= mean(log2FoldChange-log2(2^indL+2^indR-1), na.rm= T), R]

# Select pairs that are abundant in the input ----
rep1 <- fread("db/umi_counts/DSCP_large_WT_input_rep1.txt")[umi_counts>20]
rep2 <- fread("db/umi_counts/DSCP_large_WT_input_rep2.txt")[umi_counts>20]
sel <- merge(rep1,
             rep2,
             by= c("L", "R"))

# Make sure that their homotypic pairs are also abundant in the input ----
sel <- sel[L %in% reps[L==R, L] & R %in% reps[L==R, R]]

# Only keep pairs that we assessed in STARR-Seq and have limited ind act----
sel <- merge(sel,
             dat[indL<3.5 & indR<3.5, .(L, R, ctlL, ctlR, resL, resR)])

# Select best control, paired with as many enhancers as possible ----
ctlL <- sel[(ctlL) & (!ctlR), .N, .(enh= L)]
ctlR <- sel[(!ctlL) & (ctlR), .N, .(enh= R)]
ctl <- merge(ctlL, ctlR, by= "enh", suffixes= c("_L", "_R"))
plot(ctl$N_L,
     ctl$N_R,
     xlim= c(0, 350))
text(ctl$N_L,
     ctl$N_R,
     labels = ctl$enh,
     cex= .5)
# selCtl <- "control_flat_genomic_region_C_00742"
selCtl <- "control_flat_genomic_region_B_00722"

# Select best enhancers, with high average residuals ----
res <- merge(unique(sel[(!ctlL) & R==selCtl, .(enh= L, res= resL)]),
             unique(sel[L==selCtl & (!ctlR), .(enh= R, res= resR)]),
             by= "enh")
res[, sel:= res.x>1.4 & res.y>1.2 & res.x+res.y>2.9]
res[, Cc:= ifelse(sel, "red", "black")]
res[, {
  plot(res.x,
       res.y,
       col= Cc)
  # text(res.x,
  #      res.y,
  #      labels = enh,
  #      cex= .5,
  #      col= Cc)
}]
table(res$sel)

# Check individual activity
vl_boxplot(unique(dat[, .(L, indL)])$indL,
           unique(dat[, .(R, indR)])$indR)
points(rep(1, sum(res$sel)),
       unique(dat[L %in% res[(sel), enh], .(L, indL)])$indL,
       col= "red")
points(rep(2, sum(res$sel)),
       unique(dat[R %in% res[(sel), enh], .(R, indR)])$indR,
       col= "red")

dat[L %in% res[(sel), enh] & R %in% res[(sel), enh], {
  plot(log2FoldChange-log2(2^indL+2^indR-1),
       log2FoldChange,
       xlim= c(1, 8),
       ylim= c(1, 8))
}]
abline(0, 1)

# Retrieve enhancer sequences for selected ones ----
final <- rbind(data.table(enh= selCtl),
               res[(sel), .(enh)])
lib <- readRDS("Rdata/vl_library_twist008_112019.rds")
lib <- as.data.table(lib)
final[lib, enh_seq:= enh_sequence, on= "enh==ID_vl"]
final[, revSeq:= vl_revComp(enh_seq), enh_seq]
final[, checkL:= FALSE]
final[, checkR:= FALSE]
i <- 0
while(any(!final$checkL))
{
  final[(!checkL), left_primer:= substr(enh_seq, 1, 20+i)]
  final[(!checkL), GC:= vl_oligo_Tm(left_primer)$`GC%`, left_primer]
  final[(!checkL), Tm:= vl_oligo_Tm(left_primer)$Tm, left_primer]
  final[(!checkL), checkL:= Tm>55 & grepl("C$|G$", left_primer)]
  # print(table(final$checkL))
  i <- i+1 
}
i <- 0
while(any(!final$checkR))
{
  final[(!checkR), right_primer:= substr(revSeq, 1, 20+i)]
  final[(!checkR), GC:= vl_oligo_Tm(right_primer)$`GC%`, right_primer]
  final[(!checkR), Tm:= vl_oligo_Tm(right_primer)$Tm, right_primer]
  final[(!checkR), checkR:= Tm>55 & grepl("C$|G$", right_primer)]
  # print(table(final$checkR))
  i <- i+1 
}
primers <- melt(final,
                "enh",
                c("left_primer", "right_primer"))
fwrite(primers,
       "db/PCR_primers/PCR_primers_homotypic_pairs_luc.txt",
       col.names = T,
       row.names = F,
       sep= "\t")








