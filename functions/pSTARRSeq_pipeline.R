#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

############----------------------------------------------------##############
############  Align pSTARR-Seq data and generate counts table   ##############
############----------------------------------------------------##############
# Takes as input lists of gziped _1.fq and _2.fq files 
# 1/ Align them using the specified library index and saves the corresponding bam file
# 2/ Count reads for each pair
# 3/ Collapse reads based on UMI (<=1 difference) and outputs the corresponding count file

# test if there is at least 2 args: if not, return an error
if (length(args)<6) {
  stop("Please specify:\n
       [required] 1/ The type of assay. Can be one of 'pe-STARR-Seq' or 'rev-pe-STARR-Seq' \n
       [required] 2/ The index to use for the alignment \n
       [required] 3/ A list of coma-separated _1.fq files containing read 1 \n
       [required] 4/ A list of coma-separated _2.fq files containing read 2 \n
       [required] 5/ Output bam file (.bam)\n
       [required] 6/ Output .txt file where UMI-collapsed counts will be stored \n")
}
require(Rsubread)
require(data.table)
require(Rsamtools)
require(stringdist)

# Variables and checks ----
type <- args[1]
index <- args[2]
fq1 <- unlist(tstrsplit(args[3], ","))
fq2 <- unlist(tstrsplit(args[4], ","))
bam <- args[5]
umi_counts <- args[6]
if(!type %in% c("pe-STARR-Seq", "rev-pe-STARR-Seq"))
  stop("First argument 'type' should be one of 'pe-STARR-Seq' or 'rev-pe-STARR-Seq'.")
if(!grepl(".bam$", bam))
  stop("Fifth argument 'bam' should finish with .bam extension.")
if(!grepl(".txt$", umi_counts))
  stop("Sixth argument 'umi_counts' should finish with .txt extension.")

# Alignment ----
tmp1 <- tempfile(tmpdir = "/scratch/stark/vloubiere/fastq/", fileext = "_1.fq.gz")
system(paste(c("cat", fq1, ">", tmp1), collapse= " "))
tmp2 <- tempfile(tmpdir = "/scratch/stark/vloubiere/fastq/", fileext = "_2.fq.gz")
system(paste(c("cat", fq2, ">", tmp2), collapse= " "))
if(!file.exists(bam))
{
  align(index = index,
        readfile1 = tmp1,
        readfile2 = tmp2,
        type = "dna",
        output_format = "BAM",
        output_file = bam,
        maxMismatches = 3,
        unique = T,
        nthreads = getDTthreads())
}else
  print(paste0("bam file ", bam, " existed and alignment was skipped!"))

# Import bam and filter read ----
.c <- Rsamtools::scanBam(bam,
                         param = Rsamtools::ScanBamParam(what=c("qname", "rname", "strand", "pos")))[[1]]
.c <- as.data.table(.c)
setnames(.c, c("read", "seqnames", "strand", "start"))
# Extract first read and potential mates
fw <- na.omit(.c)
fw[, idx:= rowid(read)]
if(type=="pe-STARR-Seq")
{
  fw <- fw[strand=="+" & idx==1]
  .c <- .c[strand=="-", .(read, R= seqnames, start)]
}else if (type=="rev-pe-STARR-Seq")
{
  fw <- fw[strand=="-" & idx==1]
  .c <- .c[strand=="+", .(read, R= seqnames, start)]
}
# Extract paired reads
.c <- merge(fw[, .(read, L= seqnames, start)],
            .c,
            by= "read")
if(type=="pe-STARR-Seq")
{
  .c <- .c[start.x<start.y, .(read, L, R)]
}else if(type=="rev-pe-STARR-Seq")
{
  .c <- .c[start.x>start.y, .(read, L, R)]
}

# Extract UMIs ----
.c[, UMI:= gsub(".*_([A-Z]{10}).*", "\\1", read)]
.c <- .c[, .(umi_N= .N), .(L, R, UMI)]
.c[, total_counts:= sum(umi_N), .(L, R)]
setorderv(.c, "umi_N", order = -1)

# Check whether UMI might be collapsed
# .c[, collapsed:= .N==1, .(L, R)] # Old simple way
.c[, collapsed:= T, .(L, R)]
.c[, idx:= .I]
for(i in 1:10)
{
  .c[, check:= idx[1], .(L, R, gsub(paste0("^(.{", i-1, "})."), "\\1", UMI))]
  potentialDup <- unique(.c[(check<idx), c(check, idx)])
  .c[potentialDup, collapsed:= FALSE]
  print(i)
}
.c$idx <- NULL
paste0(sum(.c$collapsed), " / ", nrow(.c), " pre-collapsed")

# UMI collapsing (>1 diff) ----
while(any(!.c$collapsed))
{
  .c[!(collapsed), c("collapsed", "UMI"):= {
    coll <- stringdist(UMI[1],
                       UMI,
                       method="hamming",
                       nthread= getDTthreads()-1)<=1
    UMI[coll] <- UMI[1]
    .(coll, UMI)
  }, .(L, R)]
}
# Final collapsing
.c <- unique(.c[, .(L, R, total_counts, UMI)])
.c <- .c[, .(umi_counts= .N), .(L, R, total_counts)]

# SAVE ----
fwrite(.c, umi_counts)