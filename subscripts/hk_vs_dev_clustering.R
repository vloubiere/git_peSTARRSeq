setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
require(rtracklayer)
require(vlfunctions)

load_peak_function <- function(x, width){
  path<- "/groups/stark/almeida/data/STARRseq/dm3/20190401_200bp_gw_STARRseq/data/"
  # load peaks file and convert to granges
  gr <- makeGRangesFromDataFrame(read.delim(paste0(path, x), header = F),
                                 seqnames.field = "V1", start.field = "V2", end.field = "V2", keep.extra.columns = T)
  # resize peaks
  gr <- resize(gr, width, "center")
  # Treat peak information
  names(mcols(gr))[5:7] <- c("Enrch.", "Corr_enrch", "p_value")
  gr
}

# Import enhancers using Bernardo's function ----
dev <- load_peak_function("DSCP_200bp_gw.UMI_cut_merged.peaks.txt", 1)
dev <- as.data.table(dev)[Corr_enrch>5]
hk <- load_peak_function("RpS12_200bp_gw.UMI_cut_merged.peaks.txt", 1)
hk <- as.data.table(hk)[Corr_enrch>5]

# Import genes ----
genes <- import("/groups/stark/vloubiere/genomes/Drosophila_melanogaster/flybase/dm3/dmel-all-filtered-r5.57_genes_and_transcripts_only.gff")
seqlevelsStyle(genes) <- "UCSC"
genes <- as.data.table(genes)
genes <- genes[type=="gene"]
# Add class
class <- fread("db/public/hk_dev_FB_r5.57_2014_03.txt")
genes[class, class:= i.class, on= "ID==FBgn"]
genes <- genes[, .(class, ID, Name, seqnames, start, end, strand)]
genes <- na.omit(genes)

# Select active genes ----
act <- rbindlist(list(fread("db/public/GSM1000412_S2_RNAseq_rep1.rpkm.txt.gz"),
                      fread("db/public/GSM1000412_S2_RNAseq_rep1.rpkm.txt.gz")),
                 idcol = "rep")
act <- act[, .(FPKM= mean(V4)), .(ID= V1)]
act <- genes[ID %in% act[FPKM>1, ID]]
genes[, print(paste(round(sum(ID %in% act$ID)/.N*100, 1), "%", class, "genes are active")), class]

# Enhancer counts at the loci of active genes ----
body <- vl_resizeBed(act, "start", 5000, act[, end-start]+2000)
body[class=="developmental", count:= vl_covBed(.SD, dev)]
body[class=="housekeeping", count:= vl_covBed(.SD, hk)]

# Enhancer counts around the promoters of active genes ----
ext <- 20e3
prom <- vl_resizeBed(act, "start", ext, ext)
prom[class=="developmental", count:= vl_covBed(.SD, dev)]
prom[class=="housekeeping", count:= vl_covBed(.SD, hk)]

# Plots ----
# Save houseekping genes and promoters as temp files
tmp1 <- tempfile(fileext = ".bed")
fwrite(act[class=="housekeeping", .(seqnames, start, end)],
       tmp1,
       sep= "\t",
       col.names = F)
tmp2 <- tempfile(fileext = ".bed")
fwrite(hk,
       tmp2,
       sep= "\t",
       col.names = F)
# Screenshot region
region <- data.table(seqnames= "chr3R",
                     start= 6.9e6,
                     end= 7.9e6)

pdf("pdf/draft/hk_vs_dev_clustering.pdf",
    width = 7,
    height = 4.25)
par(mar= c(9,10,9,2),
    las= 1,
    tcl= -0.2,
    mgp= c(1.75, 0.5, 0),
    xaxs= "i")

# Screenshot
vl_screenshot(region,
              c(tmp1, tmp2),
              names = c("Houseekping genes",
                        "Housekeeping enhancers"),
              widths =  c(300L, 300L),
              col = adjustcolor("grey20", .5))
left <- par("usr")[1]
w <- par("usr")[2]-left
labs <- seq(7e6, 7.8e6, 2e5)
axis(1,
     at= left+(labs-region$start)/(region[,end-start])*w,
     formatC(labs, format = "d", big.mark = ","))
title(xlab= "chr3R")

# Enhancer counts promoter vicinity
par(mar= c(6,7,6,6),
    mfrow= c(1,2))
vl_boxplot(count~class,
           body,
           notch= F,
           compute_pval= list(c(1,2)),
           col= adjustcolor(c("limegreen", "tomato"), .5),
           names= c("Dev. CP", " Hk. CP"),
           tilt.names = T,
           ylab= paste0("Number of enhancers per gene locus"), outline= T)
vl_boxplot(count~class,
           prom,
           notch= F,
           compute_pval= list(c(1,2)),
           col= adjustcolor(c("limegreen", "tomato"), .5),
           names= c("Dev. CP", " Hk. CP"),
           tilt.names = T,
           ylab= paste0("Number of enhancers +/- ", ext/1000, " kb from TSS"))
dev.off()
