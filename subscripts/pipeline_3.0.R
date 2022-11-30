# test if there is at least one argument: if not, return an error
if (length(args)!=1) {
  stop("meta file group should be provided")
}

require(data.table)
require(parallel)
require(Biostrings)
require(seqinr)
require(Rsubread)
require(readxl)
require(stringdist)

meta <- fread(commandArgs(trailingOnly=TRUE)[1])
# meta <- fread("Rdata/metadata_processed.txt")[vllib=="vllib002" & DESeq2] # Example
# meta <- fread("Rdata/metadata_processed.txt")[vllib=="vllib015" & DESeq2] # Example
print("Sample:")
print(meta)
libs <- as.data.table(read_excel("/groups/stark/vloubiere/exp_data/vl_libraries.xlsx"))

#--------------------------------------------------------------#
# Extract from VBC bam file
# For each sequencing run, extract my reads from the bam containing the full lane
#--------------------------------------------------------------#
mcmapply(function(b, o, i){
  fq_files <- paste0(o, c("_1.fq.gz", "_2.fq.gz"))
  if(any(!file.exists(fq_files)))
  {
    cmd <- vlfunctions::vl_extract_reads_VBC(bam= b,
                                             output_prefix = o,
                                             rev_comp_i5 = i)
    system(cmd)
  }
  print(paste0(fq_files, " -->> DONE!"))
},
b= meta$BAM_path,
o= meta$fq_prefix,
i= as.character(reverseComplement(DNAStringSet(gsub(".*-(.*)", "\\1", meta$i5)))),
mc.preschedule = F,
mc.cores = getDTthreads())

#--------------------------------------------------------------#
# Alignment
# Align each fastq file and produces SAM ouptut
#--------------------------------------------------------------#
meta[, {
  if(!file.exists(bam))
  {
    subread_index <- switch(library,
                            "T8"= "/groups/stark/vloubiere/projects/pe_STARRSeq/db/subread_indexes/twist8_lib/twist8",
                            "T12"= "/groups/stark/vloubiere/projects/pe_STARRSeq/db/subread_indexes/twist12_lib/twist12")
    align(index = subread_index,
          readfile1 = fq1,
          readfile2 = fq2,
          type = "dna",
          output_format = "BAM",
          output_file = bam,
          maxMismatches = 3,
          unique = T,
          nTrim5 = 3,
          nthreads = getDTthreads())
  }
  print(paste0(bam, " -->> DONE"))
}, .(fq1, fq2, bam, library)]

#--------------------------------------------------------------#
# Primary counts
# Takes each sam file and Extract UMI reads
#--------------------------------------------------------------#
meta[, {
  if(!file.exists(umi_counts))
  {
    cmd <- paste("module load  build-env/2020; module load samtools/1.9-foss-2018b; module load bedtools/2.27.1-foss-2018b; samtools view -@")
    cmd <- paste(cmd, getDTthreads()-1, "-b", bam, "| bedtools bamtobed -i stdin -bedpe -mate1")
    .c <- fread(cmd= cmd)
    
    # Extract good reads
    .c <- .c[V8>=10] # mapq cutoff
    if(type=="rev-pe-STARR-Seq") # Orientation
      .c <- .c[V9=="-" & V10=="+" & V2>V6] else
        .c <- .c[V9=="+" & V10=="-" & V5>V3]
    .c <- .c[, .(L= V1, R= V4, UMI= gsub(".*_([A-Z]{10}).*", "\\1", V7))] # Extract UMI
    # SAVE
    fwrite(.c, umi_counts)
  }
  print(paste0(umi_counts, " -->>DONE"))
}, .(bam, umi_counts, type, library)]

#--------------------------------------------------------------#
# Merged counts
# Takes counts from separated runs and merge them per condition
# +UMI collapsing with 1nt difference
# +annotation based on library
#--------------------------------------------------------------#
meta[, {
  if(!file.exists(pairs_counts))
  {
    print(paste0("Start ", pairs_counts))
    # Import all counts file / cdition
    .c <- lapply(umi_counts, fread)
    .c <- rbindlist(.c)
    # Filter pairs containing correct subllib indexes
    .c <- .c[grepl(paste0("_", gsub(",", "_|_", libs[lib_ID==vllib, sub_lib_L]), "_"), L) 
             & grepl(paste0("_", gsub(",", "_|_", libs[lib_ID==vllib, sub_lib_R]), "_"), R)]
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
    fwrite(.c, pairs_counts)
  }
  print(paste0(pairs_counts, " -->> DONE"))
}, .(vllib, pairs_counts)]

#--------------------------------------------------------------#
# Compute activity using counts norm (Tolerates one rep)
#--------------------------------------------------------------#
if(any(meta$DESeq2))
{
  #--------------------------------#
  # DESeq2 based approach
  #--------------------------------#
  meta[, {
    if(!file.exists(dds_file))
    {
      # Import counts
      dat <- SJ(file= pairs_counts, 
                cdition= cdition, 
                rep= DESeq2_pseudo_rep)
      setkeyv(dat, c("cdition", "rep"))
      dat <- dat[, fread(file), (dat)]
      counts <- dcast(dat, 
                      L+R~cdition+rep, 
                      value.var = "umi_counts", 
                      fun.aggregate = sum, 
                      sep= "_rep")
      # Remove homotypic pairs and cutoff low input counts
      check <- apply(counts[, .SD, .SDcols= patterns("input")], 1, function(x) all(x>=5))
      counts <- counts[L!=R & (check)]
      cols <- grep("rep", names(counts), value = T)
      # Add pseudocount
      counts[, (cols):= lapply(.SD, function(x) x+1), .SDcols= cols]
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
      dds <- try(DESeq2::DESeq(dds)) # Check if missing replicates -> skip
      if(class(dds)=="DESeqDataSet")
        saveRDS(dds, dds_file)
    }
    if(!file.exists(FC_file_DESeq2) && file.exists(dds_file))
    {
      dds <- readRDS(dds_file)
      # Differential expression
      norm <- as.data.frame(DESeq2::results(dds, contrast= c("cdition", "screen", "input")))
      norm <- as.data.table(norm, keep.rownames= T)[, c("L", "R"):= tstrsplit(rn, "__")][, .(L, R, log2FoldChange, padj)]
      # Select usable, inactive control pairs
      ctlPairs <- norm[grepl("^control", L) & grepl("^control", R)]
      inactCtlL <- ctlPairs[, mean(log2FoldChange), L][between(scale(V1), -1, 1), L]
      inactCtlR <- ctlPairs[, mean(log2FoldChange), R][between(scale(V1), -1, 1), R]
      norm[, ctlL:= L %in% inactCtlL]
      norm[, ctlR:= R %in% inactCtlR]
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
      
      #-----------------------------------------------------#
      # Define active/inactive individual enhancers and pairs
      #-----------------------------------------------------#
      norm[, actClassL:= fcase(padjL<1e-4 & indL>=log2(1.5), "active", default= "inactive")]
      norm[, actClassR:= fcase(padjR<1e-4 & indR>=log2(1.5), "active", default= "inactive")]
      norm[, actClass:= fcase(grepl("control", L), "ctl.", 
                              actClassL=="active", "enh.",
                              default = "inact.")]
      norm[, actClass:= paste0(actClass, "/")]
      norm[, actClass:= paste0(actClass,
                               fcase(grepl("control", R), "ctl.", 
                                     actClassR=="active", "enh.",
                                     default= "inact."))]
      norm[, actClass:= factor(actClass, 
                               c("ctl./ctl.",
                                 "ctl./inact.",
                                 "inact./ctl.",
                                 "inact./inact.",
                                 "enh./ctl.",
                                 "enh./inact.",
                                 "ctl./enh.",
                                 "inact./enh.",
                                 "enh./enh."))]
      
      # Retrieve genomic coordinates
      .lib <- if(library=="T8")
        readRDS("Rdata/vl_library_twist008_112019.rds")else if(library=="T12")
          readRDS("Rdata/vl_library_twist12_210610.rds")
      .lib <- as.data.table(.lib)
      if(library=="T8")
        setnames(.lib, "ID_vl", "ID")
      norm[.lib, coorL:= paste0(i.seqnames, ":", i.start, "-", i.end, ":", i.strand), on= "L==ID"]
      norm[.lib, coorR:= paste0(i.seqnames, ":", i.start, "-", i.end, ":", i.strand), on= "R==ID"]
      
      # Handle missing rep columns
      cols <- c("L", "R", "indL", "indR", "log2FoldChange", "padj",
                "actClassL", "actClassR", "actClass",
                "ctlL", "ctlR",
                "coorL", "coorR")
      cols <- cols[cols %in% names(norm)]
      # SAVE
      saveRDS(norm[, cols, with= F], FC_file_DESeq2)
    }
    .SD
  }, .(dds_file, FC_file_DESeq2, library)]
  
  #--------------------------------#
  # FC based on ratio (tolerates missing replicates)
  #--------------------------------#
  meta[, {
    if(!file.exists(FC_file_ratio))
    {
      # Import counts
      dat <- SJ(file= pairs_counts, 
                cdition= cdition, 
                rep= DESeq2_pseudo_rep)
      setkeyv(dat, c("cdition", "rep"))
      dat <- dat[, fread(file), (dat)]
      counts <- dcast(dat, 
                      L+R~cdition+rep, 
                      value.var = "umi_counts",
                      fun.aggregate = sum, 
                      sep= "_rep")
      # Remove homotypic pairs and cutoff low input counts
      check <- apply(counts[, .SD, .SDcols= patterns("input")], 1, function(x) all(x>=5))
      counts <- counts[L!=R & (check)]
      # Normalize counts
      norm <- copy(counts)
      cols <- grep("^input|^screen", names(norm), value= T)
      norm[, (cols):= lapply(.SD, function(x) x+0.5), .SDcols= cols]
      norm[, norm_input:= rowSums(.SD), .SDcols= patterns("^input")]
      norm[, norm_input:= norm_input/sum(norm_input)*1e6]
      norm[, norm_screen:= rowSums(.SD), .SDcols= patterns("^screen")]
      norm[, norm_screen:= norm_screen/sum(norm_screen)*1e6]
      # FoldChange
      norm[, log2FoldChange:= log2(norm_screen/norm_input)]
      # Select usable, inactive control pairs
      ctlPairs <- norm[grepl("^control", L) & grepl("^control", R)]
      inactCtlL <- ctlPairs[, mean(log2FoldChange), L][between(scale(V1), -1, 1), L]
      inactCtlR <- ctlPairs[, mean(log2FoldChange), R][between(scale(V1), -1, 1), R]
      norm[, ctlL:= L %in% inactCtlL]
      norm[, ctlR:= R %in% inactCtlR]
      # Center activities on the median of inactive pairs
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
      # Remove pairs for which ind act could not be computed and center log2FC
      norm <- norm[!is.na(indL) & !is.na(indR)]
      
      #-----------------------------------------------------#
      # Define active/inactive individual enhancers and pairs
      #-----------------------------------------------------#
      norm[, actClassL:= fcase(padjL<0.05 & indL>=1, "active", default= "inactive")]
      norm[, actClassR:= fcase(padjR<0.05 & indR>=1, "active", default= "inactive")]
      norm[, actClass:= fcase(grepl("control", L), "ctl.", 
                              actClassL=="active", "enh.",
                              default = "inact.")]
      norm[, actClass:= paste0(actClass, "/")]
      norm[, actClass:= paste0(actClass,
                               fcase(grepl("control", R), "ctl.", 
                                     actClassR=="active", "enh.",
                                     default= "inact."))]
      norm[, actClass:= factor(actClass, 
                               c("ctl./ctl.",
                                 "ctl./inact.",
                                 "inact./ctl.",
                                 "inact./inact.",
                                 "enh./ctl.",
                                 "enh./inact.",
                                 "ctl./enh.",
                                 "inact./enh.",
                                 "enh./enh."))]

      # Retrieve genomic coordinates
      .lib <- if(library=="T8")
        readRDS("Rdata/vl_library_twist008_112019.rds")else if(library=="T12")
          readRDS("Rdata/vl_library_twist12_210610.rds")
      .lib <- as.data.table(.lib)
      if(library=="T8")
        setnames(.lib, "ID_vl", "ID")
      norm[.lib, coorL:= paste0(i.seqnames, ":", i.start, "-", i.end, ":", i.strand), on= "L==ID"]
      norm[.lib, coorR:= paste0(i.seqnames, ":", i.start, "-", i.end, ":", i.strand), on= "R==ID"]
      
      # Handle missing rep columns
      cols <- c("L", "R", "indL", "indR", "log2FoldChange",
                "actClassL", "actClassR", "actClass",
                "ctlL", "ctlR",
                "coorL", "coorR",
                "input_rep1", "input_rep2", "screen_rep1", "screen_rep2", 
                "norm_input", "norm_screen_rep1", "norm_screen_rep2")
      cols <- cols[cols %in% names(norm)]
      # SAVE
      saveRDS(norm[, cols, with= F], FC_file_ratio)
    }
  }, .(FC_file_ratio, library)]
}