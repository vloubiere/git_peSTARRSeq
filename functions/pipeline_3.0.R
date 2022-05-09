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
    if(!file.exists(FC_file))
    {
      ####### FC based on ratio ########
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
      # Normalize counts
      norm <- copy(counts)
      inputs <- grep("^input", names(norm), value= T)
      screens <- grep("^screen", names(norm), value= T)
      norm[, norm_input:= (rowSums(.SD)+0.5)/sum(rowSums(.SD)+0.5)*1e6, .SDcols= inputs]
      norm[, paste0("norm_", screens):= lapply(.SD, function(x) (x+0.5)/sum(x+0.5)*1e6), .SDcols= screens]
      # FoldChange
      norm[, log2FoldChange:= log2(rowMeans(do.call(cbind, lapply(.SD, function(x) x/norm_input)))), .SDcols= patterns("^norm_screen")]
      # Select suitable control sequences
      ctl_L <- norm[grepl("control", L) & grepl("control", R), median(log2FoldChange, na.rm= T), L]
      ctl_L[, c("min", "max"):= as.list(range(boxplot.stats(V1)$stats))]
      norm[, ctl_L:= L %in% ctl_L[V1>=min && V1<=max, L]]
      ctl_R <- norm[grepl("control", L) & grepl("control", R), median(log2FoldChange, na.rm= T), R]
      ctl_R[, c("min", "max"):= as.list(range(boxplot.stats(V1)$stats))]
      norm[, ctl_R:= R %in% ctl_R[V1>=min && V1<=max, R]]
      # Check if individual enhancer is active
      control_pairs_log2FC <- norm[ctl_L & ctl_R, log2FoldChange]
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
      norm[, median_L:= .SD[(ctl_R), ifelse(.N>=5, median(log2FoldChange), as.numeric(NA))], L]
      norm[, median_R:= .SD[(ctl_L), ifelse(.N>=5, median(log2FoldChange), as.numeric(NA))], R]
      norm[, additive:= log2(2^median_L+2^median_R)]
      norm[, multiplicative:= median_L+median_R]
    }
    
    ####### DESeq2 ########
    if(all(c("input_rep1", "input_rep2", "screen_rep1", "screen_rep2") %in% names(norm)))
    {
      counts_cols <- grep("^input|^screen", names(norm), value = T)
      # Format sampleTable
      sampleTable <- SJ(name= counts_cols)
      sampleTable <- data.frame(sampleTable[, c("cdition", "rep"):= tstrsplit(name, "_rep")], 
                                row.names = "name")
      # Format DF
      DF <- norm[, c("L", "R", counts_cols), with= F]
      DF <- DF[, name:= paste0(L, "__", R)][, !c("L", "R")]
      DF <- data.frame(DF,
                       row.names = "name")
      # DESeq
      dds <- DESeq2::DESeqDataSetFromMatrix(countData= DF,
                                            colData= sampleTable,
                                            design= ~rep+cdition)
      controls <- norm[ctl_L & ctl_R, paste0(L, "__", R)]
      sizeFactors(dds) <- DESeq2::estimateSizeFactorsForMatrix(as.matrix(DF[rownames(DF) %in% controls,]))
      dds <- DESeq2::DESeq(dds)
      
      # Differential expression
      FC <- as.data.frame(DESeq2::results(dds, contrast= c("cdition", "screen", "input")))
      FC <- as.data.table(FC, keep.rownames= T)[, c("L", "R"):= tstrsplit(rn, "__")][, .(L, R, log2FoldChange, padj)]
      setnames(FC, 
               c("log2FoldChange", "padj"), 
               function(x) paste0("DESeq_", x))
      norm <- merge(norm,
                    FC,
                    by= c("L", "R"))
      
      # Compute expected
      norm[, DESeq_median_L:= .SD[(ctl_R), ifelse(.N>=5, median(DESeq_log2FoldChange), as.numeric(NA))], L]
      norm[, DESeq_median_R:= .SD[(ctl_L), ifelse(.N>=5, median(DESeq_log2FoldChange), as.numeric(NA))], R]
      norm[, DESeq_additive:= log2(2^DESeq_median_L+2^DESeq_median_R)]
      norm[, DESeq_multiplicative:= DESeq_median_L+DESeq_median_R]
      
      # SAVE
      saveRDS(dds, dds_file)
    }
    # SAVE
    fwrite(norm, FC_file, na = NA, sep= "\t")
  }, .(FC_file, dds_file)]
}
