setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
require(vlfunctions)

# Oliver metadata ----
meta <- readxl::read_xlsx("/groups/stark/hendy/handover/sequencing/proseq_metadata.xlsx")
meta <- as.data.table(meta)
setnames(meta, function(x) gsub(" ", "_", x))
meta <- meta[grepl("^S2_control", Sample_name)]
meta <- meta[, .(path= list.files(path_to_bed,
                                  paste0(Sample_name, ".*.uniq.bed$"),
                                  full.names = T)), Sample_name]

# Resize PROseq reads
meta[, file:= paste0("./db/bed/PROSeq/", basename(path))]
meta[, {
    if(!file.exists(file))
    {
        .c <- vl_importBed(path)
        .c <- unique(.c)
        .c <- vl_resizeBed(.c, "start", 0, 0)
        .c[, strand:= ifelse(strand=="-", "+", "-")]
        fwrite(.c,
               file,
               col.names= F,
               row.names= F,
               sep= "\t",
               quote= F,
               na= NA)
    }
    print(file)
}, .(path, file)]

# Import Olivers data ----
Oliver <- fread("db/public/GSE184183_PROseq.tss.counts.tsv.gz")
dat <- Oliver[, .(tss, rep1= Parental_0hrIAA_1, rep2= Parental_0hrIAA_2)]
annot <- fread("db/public/GSE184183_PROseq.tss.annotation.tsv.gz")
annot[, tss:= paste0(seqnames, ":", start, "-", end, ":", strand)]
dat[annot, c("class", "gene_id", "transcript_id"):= .(class, gene_id, transcript_id), on= "tss"]
dat[, c("seqnames", "coor", "strand"):= tstrsplit(tss, ":", type.convert= T)]
dat[, c("start", "end"):= tstrsplit(coor, "-", type.convert= T)]
dat$tss <- NULL
setcolorder(dat,
            c("seqnames", "start", "end", "strand"))

# Add transcript total length ----
gff <- rtracklayer::import("/groups/stark/annotations/dm3/dmel-all-filtered-r5.57_genes_and_transcripts_only.gff")
gff <- as.data.table(gff)
gff <- gff[type=="mRNA"]
dat[gff, c("start.transcript", "end.transcript"):= .(i.start, i.end), on= "transcript_id==ID"]
dat[strand=="+", dist.to.end:= end.transcript-end] # Make sure that extended window does not exceed gene end
dat[strand=="-", dist.to.end:= start-start.transcript]

# Compute coverage ----
dat[, TSS_rep1:= vl_covBed(dat, meta$file[1])]
dat[, TSS_rep2:= vl_covBed(dat, meta$file[2])]
dat <- vl_resizeBed(dat, "end", 0, 249)
# dat[, downstream_rep1:= vl_covBed(dat, meta$file[1])]
# dat[, downstream_rep2:= vl_covBed(dat, meta$file[2])]
dat[, TSS_mean:= rowMeans(.SD), .SDcols= patterns("^TSS_rep")]
dat[, mean:= rowMeans(.SD), .SDcols= patterns("^rep")]

# Select active genes ----
act <- dat[TSS_rep1>0 & TSS_rep2>0 & !is.na(class)]

# Comute enrichment class top 10% ----
count <- table(act[TSS_mean>quantile(TSS_mean, 0.9), class])
enr <- act[, {
  .t <- table(factor(act$TSS_mean>quantile(act$TSS_mean, 0.9), levels = c(T, F)),
              factor(act$class==class, levels= c(T,F)))
  .f <- fisher.test(.t,
                    alternative = "greater")
  c(.f[c("estimate", "p.value")],
    list(label= paste0(class, " (", .t[1,1], "/", sum(.t[1,]), ")")))
}, class]

# Plot ----
pdf("pdf/draft/review_PROSeq_hk_vs_dev_genes.pdf", 6, 3)
vl_par(mai= c(.9, 1, 1.2, .5),
       mfrow= c(1,3))
# vl_boxplot(log2(rep1+rep2)~class,
#            dat[rep1>0 & rep2>0],
#            compute.pval= list(c(1,2)),
#            notch= T,
#            tilt.names= T,
#            ylab= "PRO-Seq counts (log2)",
#            violin= T,
#            col= "red",
#            main= "Oliver")
vl_boxplot(log2(TSS_rep1+TSS_rep2)~class,
           act,
           compute.pval= list(c(1,2)),
           notch= T,
           tilt.names= T,
           ylab= "PRO-Seq counts (log2)",
           violin= T,
           col= "lightgrey",
           main= "Vincent")
vl_barplot(count,
           bar.labels = count,
           ylab= "Number of genes",
           names.arg = paste0(names(count), " (", count, ")"),
           main= "Top 10% act. genes")
vl_barplot(enr$estimate,
           bar.labels = paste0("p=", formatC(enr$p.value, format = "e", digits = 1)),
           bar.labels.cex = .5,
           ylab= "Odd ratio",
           names.arg= enr$label,
           main= "Top 10% act. genes", )
abline(h=1, lty= 3)
# vl_boxplot(log2(downstream_rep1+downstream_rep2)~class,
#            dat[rep1>0 & rep2>0],
#            compute.pval= list(c(1,2)),
#            notch= T,
#            tilt.names= T,
#            ylab= "PRO-Seq counts (log2)",
#            violin= T,
#            col= "lightgrey",
#            main= "Downstream")
dev.off()