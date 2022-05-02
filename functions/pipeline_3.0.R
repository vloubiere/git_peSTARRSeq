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
    print(paste0(fq_files, " -->> DONE!"))
  }else
    print(paste0(fq_files, " -->> ALREADY EXISTS!"))
},
b= meta[, BAM_path],
o= meta$fq_prefix,
i= meta[, as.character(reverseComplement(DNAStringSet(gsub(".*-(.*)", "\\1", i5))))],
mc.preschedule = F,
mc.cores = getDTthreads())

#--------------------------------------------------------------#
# Alignment
# Align each fastq file and produces SAM ouptut
#--------------------------------------------------------------#
meta[, {
  if(file.exists(bam))
    print(paste0(bam, " -->> ALREADY EXISTS")) else
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
      print(paste0(bam, " -->> DONE"))
    }
}, .(fq1, fq2, bam, library)]


#--------------------------------------------------------------#
# Primary counts
# Takes each sam file and Extract UMI reads
#--------------------------------------------------------------#
meta[, {
  if(file.exists(umi_counts))
    print(paste0(umi_counts, " -->> ALREADY EXISTS")) else
    {
      .c <- fread(cmd= paste("module load  build-env/2020; module load samtools/1.9-foss-2018b; samtools view -@", getDTthreads()-1, bam),
                  fill= T,
                  select = 1:5, 
                  col.names = c("ID", "FLAG", "seq", "pos", "mapq"))
      # Extract good reads
      if(type=="rev-pe-STARR-Seq")
      {
        .c <- merge(.c[FLAG %in% c(81, 83) & mapq>=20, .(ID, seq, pos)],
                    .c[FLAG %in% c(161, 163) & mapq>=20, .(ID, seq, pos)], 
                    by= "ID",
                    suffixes= c("_L", "_R"))
      }else
      {
        .c <- merge(.c[FLAG %in% c(97, 99) & mapq>=20, .(ID, seq, pos)],
                    .c[FLAG %in% c(145, 147) & mapq>=20, .(ID, seq, pos)], 
                    by= "ID",
                    suffixes= c("_L", "_R"))
      }
      # Clean
      .c <- .c[, .(L= seq_L,
                   pos_L,
                   R= seq_R,
                   pos_R,
                   UMI= gsub(".*_([A-Z]{10}).*", "\\1", ID))]
      # SAVE
      fwrite(.c, umi_counts)
      print(paste0(umi_counts, " -->>DONE"))
    }
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
    # Filter correct pairs
    .c <- .c[grepl(paste0("_", gsub(",", "_|_", libs[lib_ID==vllib, sub_lib_L]), "_"), L) 
             & grepl(paste0("_", gsub(",", "_|_", libs[lib_ID==vllib, sub_lib_R]), "_"), R)]
    if(vllib=="vllib006")
    {
      .c <- .c[between(pos_L, 260, 264, incbounds = T) 
               & between(pos_R, 4, 7, incbounds = T)]
    }else if(vllib=="vllib002")
    {
      .c <- .c[pos_L==5 & pos_R %in% c(260, 264)]
    }else if(vllib %in% c("vllib013", "vllib014"))
    {
      .c <- .c[between(pos_L, 4, 7, incbounds = T) 
               & between(pos_R, 260, 263, incbounds = T)] 
    }else if(vllib %in% c("vllib015", 
                          "vllib016", 
                          "vllib017", 
                          "vllib018", 
                          "vllib019", 
                          "vllib020", 
                          "vllib021", 
                          "vllib022", 
                          "vllib023", 
                          "vllib024", 
                          "vllib025", 
                          "vllib026", 
                          "vllib027", 
                          "vllib028"))
    {
      .c <- .c[between(pos_L, 4, 7, incbounds = T) 
               & between(pos_R, 261, 265, incbounds = T)] 
    }else
    {
      .c <- .c[0]
    }
    # Count UMIs and order
    .c <- .c[, .(umi_N= .N), .(L, R, UMI)]
    setorderv(.c, "umi_N", order = -1)
    # Advanced UMI collapsing (>1 diff)
    .c[, UMI:= {
      if(.N>1)
      {
        res <- rep(as.character(NA), .N)
        while(anyNA(res))
        {
          check <- is.na(res)
          check[check] <- stringdist(UMI[check][1],
                                     UMI[check],
                                     method="hamming",
                                     nthread= getDTthreads()-1)<=1
          res[check] <- UMI[check][1]
        }
        res
      }else
        UMI
    }, .(L, R)]
    # Final collapsing
    .c <- unique(.c[, .(L, R, UMI)])
    .c <- .c[, .(umi_counts= .N), .(L, R)]
    # SAVE
    fwrite(.c, pairs_counts)
  }
  print(paste0(pairs_counts, " -->> DONE"))
}, .(vllib, pairs_counts)]

#--------------------------------------------------------------#
# Method using counts norm (Tolerates one rep)
#--------------------------------------------------------------#
if(any(meta$DESeq2))
{
  meta[, {
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
    # Remove homotypic pairs and cutoff low counts
    counts <- counts[L!=R & rowSums(counts[, !c("L", "R")])>20]
    
    ####### DESeq2 ########
    if(!file.exists(FC_file_DESeq) & 
       all(c("input_rep1", "input_rep2", "screen_rep1", "screen_rep2") %in% names(counts)))
      {
        # Format sampleTable
        sampleTable <- SJ(name= setdiff(names(counts), c("L","R")))
        sampleTable <- data.frame(sampleTable[, c("cdition", "rep"):= tstrsplit(name, "_rep")], 
                                  row.names = "name")
        # Format DF
        DF <- counts[, name:= paste0(L, "__", R)][, !c("L", "R")]
        DF <- data.frame(DF, 
                         row.names = "name")
        # DESeq
        dds <- DESeq2::DESeqDataSetFromMatrix(countData= DF,
                                              colData= sampleTable,
                                              design= ~rep+cdition)
        sizeFactors(dds) <- DESeq2::estimateSizeFactorsForMatrix(as.matrix(DF[grep("control.*__control.*", rownames(DF)),]))
        dds <- DESeq2::DESeq(dds)

        # Differential expression
        FC <- as.data.frame(DESeq2::results(dds, contrast= c("cdition", "screen", "input")))
        FC <- as.data.table(FC, keep.rownames= T)[, c("L", "R"):= tstrsplit(rn, "__")][, .(L, R, log2FoldChange, padj)]

        # Compute expected
        FC[, median_L:= .SD[grepl("^control", R), ifelse(.N>5, median(log2FoldChange), as.numeric(NA))], L]
        FC[, median_R:= .SD[grepl("^control", L), ifelse(.N>5, median(log2FoldChange), as.numeric(NA))], R]
        FC[, additive:= log2(2^median_L+2^median_R)]
        FC[, multiplicative:= median_L+median_R]
        
        # SAVE
        saveRDS(dds, dds_file)
        fwrite(FC, FC_file_DESeq, sep= "\t", na = NA)
        print(paste0(FC_file_DESeq, "  -->> DONE"))
    }
    
    ####### Ratios ########
    if(!file.exists(FC_file_ratio))
    {
      # Format counts
      norm <- copy(counts)
      inputs <- grep("^input", names(norm), value = T)
      norm[, input:= rowSums(.SD), .SDcols= inputs]
      norm <- norm[, !(inputs), with= F]
      # Pseudocount normalized counts
      cols <- grep("^screen|^input$", names(norm), value= T)
      norm[, (cols):= lapply(.SD, function(x) (x+0.5)/sum(x+0.5)*1e6), .SDcols= cols]
      # FoldChange
      norm[, log2FoldChange:= log2(rowMeans(do.call(cbind, lapply(.SD, function(x) x/input)))), .SDcols= patterns("^screen_rep")]
      # Check if individual enhancer is active
      control_pairs_log2FC <- norm[grepl("control", L) & grepl("control", R), log2FoldChange]
      norm[, act_wilcox_L:= {
        .c <- log2FoldChange[grepl("control", R)]
        if(length(.c)>5)
          wilcox.test(.c, control_pairs_log2FC, alternative = "greater")$p.value else
            as.numeric(NA)
      }, L]
      norm[, act_wilcox_R:= {
        .c <- log2FoldChange[grepl("control", L)]
        if(length(.c)>5)
          wilcox.test(.c, control_pairs_log2FC, alternative = "greater")$p.value else
            as.numeric(NA)
      }, R]
      # Subtract basal activity (center controls on 0)
      norm[, log2FoldChange:= log2FoldChange-median(control_pairs_log2FC)]
      # Compute expected
      norm[, median_L:= .SD[grepl("^control", R), ifelse(.N>5, median(log2FoldChange), as.numeric(NA))], L]
      norm[, median_R:= .SD[grepl("^control", L), ifelse(.N>5, median(log2FoldChange), as.numeric(NA))], R]
      norm[, additive:= log2(2^median_L+2^median_R)]
      norm[, multiplicative:= median_L+median_R]
      # SAVE
      norm <- na.omit(norm[, .(L, R, log2FoldChange, median_L, median_R, act_wilcox_L, act_wilcox_R, additive, multiplicative)])
      fwrite(norm, FC_file_ratio, na = NA, sep= "\t")
      print(paste0(FC_file_ratio, "  -->> DONE"))
    }
  }, .(FC_file_DESeq, dds_file, FC_file_ratio)]
}
