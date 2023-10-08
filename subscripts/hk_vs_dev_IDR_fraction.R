setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")

# Import IDR fraction database ----
db <- fread("db/public/mobidb_result_2023-08-31T08_01_32.511Z.tsv")
add <- AnnotationDbi::select(org.Dm.eg.db::org.Dm.eg.db,
                             keys = db$acc, 
                             columns = "FLYBASE",
                             keytype = "UNIPROT")
db[as.data.table(add), FBgn:= i.FLYBASE, on= "acc==UNIPROT"]
db <- na.omit(db)

idr <- data.table(FBgn= unique(db$FBgn))
idr <- merge(idr,
             db[feature=="prediction-disorder-mobidb_lite", .(FBgn, content_fraction)],
             # db[feature=="prediction-disorder-alphafold", .(FBgn, content_fraction)],
             by= "FBgn",
             all.x= T)
idr[is.na(content_fraction), content_fraction:= 0]
idr <- idr[, .(content_fraction= max(content_fraction)), FBgn]

# Import data Stampfel 2015 nature (dev/hk preference) ----
LEGO <- as.data.table(readxl::read_xlsx("db/public/LEGO_41586_2015_BFnature15545_MOESM9_ESM.xlsx"))
setnames(LEGO, "FBGN", "FBgn")
dat <- merge(LEGO,
             idr,
             by= "FBgn")
dat[, class:= fcase(`4xUAS upstream hkCP FC-LOG2`>`4xUAS dCP FC-LOG2` & `4xUAS upstream hkCP FC-LOG2`>log2(1.5), "Housekeeping",
                    `4xUAS upstream hkCP FC-LOG2`<`4xUAS dCP FC-LOG2` & `4xUAS dCP FC-LOG2`>log2(1.5), "Developmental",
                    default= "None")]
dat <- dat[class!="None"]
dat[, class:= factor(class, c("Housekeeping","Developmental"))]

# Add FPKMs Drosophila S2 cells ----
RPKM <- rbindlist(list(rep1= fread("db/public/GSM1000412_S2_RNAseq_rep1.rpkm.txt.gz"),
                       rep2= fread("db/public/GSM1000413_S2_RNAseq_rep2.rpkm.txt.gz")),
                  idcol = "rep")
setnames(RPKM,
         old= c("V1", "V4"),
         new= c("FBgn", "rpkm"))
RPKM <- RPKM[, .(rpkm= mean(rpkm, na.rm= T)), FBgn]
dat[RPKM, rpkm:= i.rpkm, on= "FBgn"]
dat[is.na(rpkm), rpkm:= 0]
dat[RPKM, rpkm:= log2(i.rpkm+1), on= "FBgn"]

# Select active genes
dat <- dat[rpkm>0.1]

# Plot ----
pdf("pdf/draft/hk_vs_dev_IDRs.pdf",
    width = 3, 
    height= 3)
par(mai= c(1,1.25,0.75,1.1), 
    mgp= c(0.75, 0.25, 0),
    cex.lab= 8/12,
    cex.axis= 7/12,
    las= 1,
    tcl= -0.1,
    bty= "n",
    lend= 2,
    lwd= .75)
Cc <- c("tomato", "limegreen")
vl_boxplot(content_fraction~class,
           dat,
           col= "white",
           names= dat[, paste0(class, " (n= ", .N, ")"), keyby= class]$V1,
           viocol= adjustcolor(Cc, .3),
           compute_pval= list(c(1,2)),
           tilt.names= T,
           notch= T,
           ylab= "IDR fraction",
           violin= T)
dev.off()

