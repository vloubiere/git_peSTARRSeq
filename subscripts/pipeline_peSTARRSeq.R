#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

# test if there is at least one argument: if not, return an error
if (length(args)<4) {
  stop("Please specify:\n
       [required] 1/ library name\n
       [required] 2/ library type ('pe-STARR-Seq' or 'rev-pe-STARR-Seq')\n
       [required] 3/ Rsubread index prefix\n
       [required] 4/ Sample sheet containint BAM_path, i5, cdition, and rep columns\n
       [optional] 5/ Regular expression used to select correct L oligos (default to '.*')\n
       [optional] 6/ Regular expression used to select correct R oligos (default to '.*')\n")
}

require(data.table)
require(parallel)
require(Biostrings)
require(seqinr)
require(Rsubread)
require(readxl)
require(stringdist)

# Variables
vllib <- args[1]
type <- args[2]
index <- args[3]
meta <- fread(args[4])

# Checks
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

# # Example
# vllib <- lib <- "vllib030"
# type <-  "pe-STARR-Seq"
# index <-  "/groups/stark/vloubiere/projects/pe_STARRSeq/db/subread_indexes/twist15_lib/twist15"
# meta <- read_xlsx("/groups/stark/vloubiere/exp_data/vl_sequencing_metadata.xlsx")
# meta <- as.data.table(meta)[(DESeq2) & vllib==lib]
# meta <- data.table(meta[, .(BAM_path, i5, cdition, rep= DESeq2_pseudo_rep)])
# sublibRegexprL <- "_A_"
# sublibRegexprR <- "_A_"

# Generate directories and file names
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

# Check whether DESeq2 can be applied
check <- meta[, length(unique(rep)), cdition]
diff <- all(c("input", "screen") %in% check$cdition)
DESeq <- all(check$V1>=2)
dds_file <- paste0("/groups/stark/vloubiere/projects/pe_STARRSeq/db/dds/", vllib, ".dds")
FC_file <- paste0("/groups/stark/vloubiere/projects/pe_STARRSeq/db/FC_tables/", 
                  vllib, 
                  ifelse(DESeq, "_DESeq2", "_ratio"),
                  ".rds")

#--------------------------------------------------------------#
# Extract from VBC bam file
# For each sequencing run, extract my reads from the bam containing the full lane
#--------------------------------------------------------------#
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
i= as.character(reverseComplement(DNAStringSet(gsub(".*-(.*)", "\\1", meta$i5)))),
mc.preschedule = F,
mc.cores = getDTthreads())

#--------------------------------------------------------------#
# Alignment
# Align each fastq file and produces SAM ouptut
#--------------------------------------------------------------#
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

#--------------------------------------------------------------#
# UMI counts
# +UMI collapsing with 1nt difference
# +annotation based on library
#--------------------------------------------------------------#
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
    
    # Advanced UMI collapsing (>1 diff)
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

#--------------------------------#
# Compute log2FoldChange
#--------------------------------#
if(!file.exists(FC_file)) 
{
  # Import counts
  dat <- unique(meta[, .(file= umi_counts, cdition, rep)])
  setkeyv(dat, c("cdition", "rep"))
  dat <- dat[, fread(file), (dat)]
  counts <- dcast(dat, 
                  L+R~cdition+rep, 
                  value.var = "umi_counts", 
                  fun.aggregate = sum)
  # Remove homotypic pairs and cutoff low input counts
  check <- apply(counts[, .SD, .SDcols= patterns("input")], 1, function(x) all(x>=5))
  counts <- counts[L!=R & (check)]
  cols <- grep("rep", names(counts), value = T)
  # Add pseudocount
  counts[, (cols):= lapply(.SD, function(x) x+1), .SDcols= cols]
  
  # DESeq2 method ----------#
  if(DESeq) 
  {
    if(!file.exists(dds_file))
    {
      # DF and sampleTable
      DF <- data.frame(counts[, !c("L", "R")], row.names = counts[, paste0(L, "__", R)])
      sampleTable <- SJ(name= colnames(DF))
      sampleTable <- data.frame(sampleTable[, c("cdition", "rep"):= tstrsplit(name, "_rep")], 
                                row.names = "name")
      # DESeq
      dds <- DESeq2::DESeqDataSetFromMatrix(countData= DF,
                                            colData= sampleTable,
                                            design= ~rep+cdition)
      sizeFactors(dds) <- DESeq2::estimateSizeFactorsForMatrix(as.matrix(DF[grepl("^control.*__control.*", rownames(DF)),]))
      dds <- DESeq2::DESeq(dds) # Check if missing replicates -> skip
      saveRDS(dds, dds_file)
    }
    # Differential expression
    dds <- readRDS(dds_file)
    norm <- as.data.frame(DESeq2::results(dds, contrast= c("cdition", "screen", "input")))
    norm <- as.data.table(norm, keep.rownames= T)[, c("L", "R"):= tstrsplit(rn, "__")][, .(L, R, log2FoldChange, padj)]
  }else
  {
    # Ratio method (tolerates one replicate) ----------#
    # Normalize counts
    norm <- copy(counts)
    cols <- grep("^input|^screen", names(norm), value= T)
    norm[, norm_input:= rowSums(.SD), .SDcols= patterns("^input")]
    norm[, norm_input:= norm_input/sum(norm_input)*1e6]
    norm[, norm_screen:= rowSums(.SD), .SDcols= patterns("^screen")]
    norm[, norm_screen:= norm_screen/sum(norm_screen)*1e6]
    # FoldChange
    norm[, log2FoldChange:= log2(norm_screen/norm_input)]
  }

  # Select usable, inactive control pairs
  ctlPairs <- norm[grepl("^control", L) & grepl("^control", R)]
  inactCtlL <- ctlPairs[, mean(log2FoldChange), L][between(scale(V1), -1, 1), L]
  inactCtlR <- ctlPairs[, mean(log2FoldChange), R][between(scale(V1), -1, 1), R]
  norm[, ctlL:= L %in% inactCtlL]
  norm[, ctlR:= R %in% inactCtlR]
  if(!DESeq) # Ratio method -> center on inactive pairs (already done earlier for DESeq2)
    norm[, log2FoldChange:= log2FoldChange-median(norm[ctlL & ctlR, log2FoldChange])]
  # Compute individual act and pval (vs ctlPairs)
  ctlPairs <- norm[ctlL & ctlR, log2FoldChange]
  Left <- norm[(ctlR), 
               .(.N>=10, 
                 wilcox.test(log2FoldChange, ctlPairs, alternative = "greater")$p.value,
                 mean(log2FoldChange)), L][(V1)]
  Left[, padj:= p.adjust(V2, "fdr")]
  norm[Left, c("indL", "padjL"):= .(i.V3, i.padj), on= "L"]
  Right <- norm[(ctlL), 
                .(.N>=10, 
                  wilcox.test(log2FoldChange, ctlPairs, alternative = "greater")$p.value,
                  mean(log2FoldChange)), R][(V1)]
  Right[, padj:= p.adjust(V2, "fdr")]
  norm[Right, c("indR", "padjR"):= .(i.V3, i.padj), on= "R"]

  # Remove pairs for which combined or ind act could not be computed accurately
  norm <- norm[!is.na(indL) & !is.na(indR)]

  # Add activity classes
  norm[, groupL:= fcase(grepl("^control", L), "Random seq.", default = "Candidate seq.")]
  norm[, groupR:= fcase(grepl("^control", R), "Random seq.", default = "Candidate seq.")]
  
  # Split Active enhancers based on their strength
  norm[, actL:= fcase(between(indL, 1, 2) & padjL<0.05, "Low",
                      between(indL, 2, 4) & padjL<0.05, "Medium",
                      between(indL, 4, Inf) & padjL<0.05, "Strong",
                      default = "Inactive")]
  norm[, actR:= fcase(between(indR, 1, 2) & padjR<0.05, "Low",
                      between(indR, 2, 4) & padjR<0.05, "Medium",
                      between(indR, 4, Inf) & padjR<0.05, "Strong",
                      default = "Inactive")]
  cols <- c("actL", "actR")
  norm[, (cols):= lapply(.SD, function(x) factor(x, c("Inactive", "Low", "Medium", "Strong"))), .SDcols= cols]
  
  # Handle missing columns
  cols <- c("L", "R",
            "groupL", "groupR", "actL", "actR",
            "indL", "padjL", "indR", "padjR", 
            "log2FoldChange", "padj",
            "ctlL", "ctlR")
  cols <- cols[cols %in% names(norm)]
  norm <- norm[, cols, with= F]
  # SAVE
  saveRDS(norm, FC_file)
}
