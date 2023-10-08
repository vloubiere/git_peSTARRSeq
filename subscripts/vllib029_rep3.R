setwd("/groups/stark/vloubiere/projects/pe_STARRSeq/")
require(data.table)
require(parallel)
require(Biostrings)
require(seqinr)
require(Rsubread)
require(readxl)
require(stringdist)

# Variables ----
vllib <- "vllib029"
type <-  "pe-STARR-Seq"
index <-  "/groups/stark/vloubiere/projects/pe_STARRSeq/db/subread_indexes/twist15_lib/twist15"
meta <- readxl::read_xlsx("/groups/stark/vloubiere/exp_data/vl_sequencing_metadata.xlsx")
meta <- as.data.table(meta)[vllib=="vllib029" & rep=="rep3", .(BAM_path, i5, cdition, rep= DESeq2_pseudo_rep)]

# Checks ----
if(!type %in% c("pe-STARR-Seq", "rev-pe-STARR-Seq"))
  stop("type should be one of 'pe-STARR-Seq' or 'rev-pe-STARR-Seq'")
if(!all(c("BAM_path", "i5", "cdition", "rep") %in% names(meta)))
  stop("meta should be a data.table containing 'BAM_path', 'i5', 'cdition', 'rep' columns")
if(length(args)>4)
  sublibRegexprL <- args[5] else
    sublibRegexprL <- ".*"
if(length(args)>5)
  sublibRegexprR <- args[6] else
    sublibRegexprR <- ".*"

# Generate directories and file names ----
tmpFolder <- paste0("/scratch/stark/vloubiere/", vllib, "/")
dir.create(tmpFolder, showWarnings = F)
fqFolder <- paste0(tmpFolder, "fq/")
dir.create(fqFolder, showWarnings = F)
meta[, fq1:= paste0(fqFolder, gsub(".bam$", "", basename(BAM_path)), "_", i5, "_1.fq.gz"), BAM_path]
meta[, fq2:= paste0(fqFolder, gsub(".bam$", "", basename(BAM_path)), "_", i5, "_2.fq.gz"), BAM_path]
bamFolder <- paste0(tmpFolder, "bam/")
dir.create(bamFolder, showWarnings = F)
meta[, bam:= paste0(bamFolder, vllib, "_", cdition, "_", rep, ".bam"), ]
meta[, umi_counts:= paste0("/groups/stark/vloubiere/projects/pe_STARRSeq/db/umi_counts/", vllib, "_", cdition, "_", rep, ".txt")]
meta[, align:= !file.exists(umi_counts)]

# Extract from VBC bam file ----
# For each sequencing run, extract my reads from the bam containing the full lane
mcmapply(function(b, o, fq1, fq2, i, align){
  if(align & (!file.exists(fq1) | !file.exists(fq2)))
  {
    cmd <- vlfunctions::vl_extract_reads_VBC(bam= b,
                                             output_prefix = o,
                                             rev_comp_i5 = i)
    system(cmd)
  }
  print(".")
},
b= meta$BAM_path,
o= gsub("_1.fq.gz", "", meta$fq1),
fq1= meta$fq1,
fq2= meta$fq2,
align= meta$align,
i= as.character(Biostrings::reverseComplement(Biostrings::DNAStringSet(gsub(".*-(.*)", "\\1", meta$i5)))),
mc.preschedule = F,
mc.cores = getDTthreads())

# Alignment ----
# Align each fastq file and produces SAM ouptut
meta[, {
  if(align & !file.exists(bam))
  {
    tmp1 <- tempfile(tmpdir = fqFolder, fileext = "_1.fq.gz")
    system(paste(c("cat", fq1, ">", tmp1), collapse= " "))
    tmp2 <- tempfile(tmpdir = fqFolder, fileext = "_2.fq.gz")
    system(paste(c("cat", fq2, ">", tmp2), collapse= " "))
    align(index = index,
          readfile1 = tmp1,
          readfile2 = tmp2,
          type = "dna",
          output_format = "BAM",
          output_file = bam,
          maxMismatches = 3,
          unique = T,
          nthreads = getDTthreads())
  }
  print(paste0(bam, " -->> DONE"))
}, .(bam, align)]

# UMI counts ----
# +UMI collapsing with 1nt difference
# +annotation based on library
meta[, {
  if(!file.exists(umi_counts))
  {
    # Import aligned reads
    cmd <- paste("module load  build-env/2020; module load samtools/1.9-foss-2018b; module load bedtools/2.27.1-foss-2018b; samtools view -@")
    cmd <- paste(cmd, getDTthreads()-1, "-b", bam, "| bedtools bamtobed -i stdin -bedpe -mate1")
    .c <- fread(cmd= cmd)
    
    # Extract good reads
    if(type=="pe-STARR-Seq")
      .c <- .c[V9=="+" & V10=="-" & V5>V2] else if(type=="rev-pe-STARR-Seq")
        .c <- .c[V9=="-" & V10=="+" & V2>V5]
    .c <- .c[, .(L= V1, R= V4, UMI= gsub(".*_([A-Z]{10}).*", "\\1", V7))] # Extract UMI
    
    # Filter pairs containing correct patterns
    .c <- .c[grepl(sublibRegexprL, L) & grepl(sublibRegexprR, R)]
    
    ## Advanced UMI collapsing (>1 diff) ----
    .c <- .c[, .(umi_N= .N), .(L, R, UMI)]
    .c[, total_counts:= sum(umi_N), .(L, R)]
    setorderv(.c, "umi_N", order = -1)
    .c[, collapsed:= .N==1, .(L, R)]
    while(any(!.c$collapsed))
    {
      .c[!(collapsed), c("collapsed", "UMI"):= {
        coll <- stringdist(UMI[1],
                           UMI,
                           method="hamming",
                           nthread= getDTthreads()-1)<=1
        UMI[coll] <- UMI[1]
        .(coll, UMI)
      }, .(L, R, total_counts)]
    }
    
    # Final collapsing
    .c <- unique(.c[, .(L, R, total_counts, UMI)])
    .c <- .c[, .(umi_counts= .N), .(L, R, total_counts)]
    
    # SAVE
    fwrite(.c, umi_counts)
  }
  print(paste0(umi_counts, " -->>DONE"))
}, .(bam, umi_counts)]