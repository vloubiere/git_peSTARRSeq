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

meta <- fread(commandArgs(trailingOnly=TRUE)[1])
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
  if(file.exists(sam))
    print(paste0(sam, " -->> ALREADY EXISTS")) else
    {
      subread_index <- switch(library, 
                              "T8"= "/groups/stark/vloubiere/projects/pe_STARRSeq/db/subread_indexes/twist8_lib/twist8",
                              "T12"= "/groups/stark/vloubiere/projects/pe_STARRSeq/db/subread_indexes/twist12_lib/twist12")
      align(index = subread_index,
            readfile1 = fq1,
            readfile2 = fq2,
            output_format = "SAM",
            output_file = sam,
            maxMismatches = 3,
            unique = T,
            nTrim5 = 3,
            nthreads = getDTthreads())
      print(paste0(sam, " -->> DONE"))
    }
}, .(fq1, fq2, sam, library)]

#--------------------------------------------------------------#
# Primary counts
# Takes each sam file and Extract UMI reads
#--------------------------------------------------------------#
meta[, {
  if(file.exists(umi_counts))
    print(paste0(umi_counts, " -->> ALREADY EXISTS")) else 
    {
      .c <- fread(sam,
                  skip = switch(library, "T8"= 1006, "T12"= 395),
                  header= F,
                  fill= T,
                  select = c(1,3,4,7,8))
      # Keep pairs with two mates are aligned
      .c <- .c[V3!="*" & V7!="*"]
      # Select firt read
      .c <- .c[, .SD[1], V1]
      # Keep pairs with reads on opposite enh ends
      if(type=="pe-STARR-Seq")
        .c <- .c[(V8-V4)>230]
      if(type=="rev-pe-STARR-Seq")
        .c <- .c[(V4-V8)>230]
      # Extract UMI
      .c <- .c[, .(L= V3, R= V7, UMI= gsub(".*_([A-Z]{10}).*", "\\1", V1))]
      # Compute total reads
      stat <- data.table(total_reads= nrow(.c))
      # When red 1 and 2 are frome the same enhancer, read 2 seqnems is "="
      .c[R=="=", R:= L]
      # Compute collapsed reads
      stat[, umi_collapsed_reads:= nrow(unique(.c))]
      # SAVE
      fwrite(.c, umi_counts)
      fwrite(stat, umi_summary)
      print(paste0(umi_counts, " -->>DONE"))
    }
}, .(sam, umi_counts, umi_summary, type, library)]

#--------------------------------------------------------------#
# Merged counts
# Takes counts from separated runs and merge them per condition
# +UMI collapsing with 1nt difference
# +annotation based on library
#--------------------------------------------------------------#
meta[, {
  if(all(file.exists(c(summary_counts, pairs_counts, spike_counts, switched_counts))))
    print(paste0("All merged counts files -->> ALREADY EXISTS")) else
    {
      print(paste0("Start ", pairs_counts))
      # Import all counts file / cdition
      .c <- lapply(umi_counts, fread)
      .c <- rbindlist(.c)
      # Count UMIs and order
      .c <- .c[, .(umi_N= .N), .(L, R, UMI)]
      setorderv(.c, "umi_N", order = -1)
      # Remove problematic UMIs
      .c <- .c[!agrepl("GGGGGGGGGG", UMI)]
      # Identify UMIs that might require collapsing
      .c[, done:= T]
      for(i in 0:9) 
        .c[(done), done:= ifelse(.N>1, F, T), .(L, R, sub(paste0("(.{", i, "})."), "\\1", UMI))]
      # Advanced UMI collapsing (>1 diff)
      while(any(!.c$done))
      {
        .c[!(done), c("UMI", "done") := {
          idx <- agrepl(UMI[1], UMI)
          .(ifelse(idx, UMI[1], UMI), idx)
        }, .(L, R)]
      }
      # Final collapsing
      .c <- unique(.c[, .(L, R, UMI)])
      .c <- .c[, .(umi_counts= .N), .(L, R)]
      # Compute patterns to identify library pairs
      .ex <- libs[lib_ID %in% c(vllib, Spike_in)]
      .ex <- .ex[, .(L_pattern= strsplit(sub_lib_L, ";")[[1]],
                     R_pattern= strsplit(sub_lib_R, ";")[[1]]),(.ex)]
      .ex <- .ex[, CJ(L_pattern= strsplit(L, ",")[[1]],
                      R_pattern= strsplit(R, ",")[[1]]), .(lib_ID, L= L_pattern, R= R_pattern)]
      .ex[, c("L_pattern", "R_pattern"):= .(paste0("_", L_pattern, "_"),
                                            paste0("_", R_pattern, "_"))]
      .ex$L <- NULL
      .ex$R <- NULL
      .ex[, Spike:= ifelse(lib_ID==Spike_in & !is.na(Spike_in), T, F)]
      # Check if pair exists and is spike in
      .ex[, {
        .c[grepl(L_pattern, L) & grepl(R_pattern, R), type:= ifelse(Spike, "spike-in", "pair")]
      }, .(L_pattern, R_pattern, Spike)]
      .c[is.na(type), type:= "switched"]
      # Generate read summary
      sum_files <- data.table(file= sam_summary)
      sum_files[, mapped:= fread(file)[V1=="Mapped_fragments", V2], file]
      .summary <- sum_files[, .(mapped= sum(mapped), collapsed= sum(.c$umi_counts))]
      fwrite(.summary, summary_counts)
      # SAVE
      fwrite(.c[type=="pair", !"type"], pairs_counts)
      fwrite(.c[type=="spike-in", !"type"], spike_counts)
      fwrite(.c[type=="switched", !"type"], switched_counts)
      print(paste0(pairs_counts, " -->> DONE"))
    }
}, .(vllib, Spike_in, summary_counts, pairs_counts, spike_counts, switched_counts)]

#--------------------------------------------------------------#
# Method using counts norm (Tolerates one rep)
#--------------------------------------------------------------#
if(any(meta$DESeq2))
{
  meta[, {
    if(file.exists(FC_file))
      print(paste0(FC_file, "  -->> ALREADY EXISTS")) else
      {
        # Import counts
        dat <- unique(data.table(file= pairs_counts, 
                                 rep= DESeq2_pseudo_rep,
                                 cdition))
        dat[cdition=="screen", cdition:= paste0(cdition, "_rep", rep)]
        dat <- dat[, fread(file), (dat)]
        # Collapse input counts and screen reps
        dat <- dat[, .(umi_counts= sum(umi_counts)), .(cdition, L, R)]
        # Filter
        dat[, check_counts:= sum(umi_counts)>=10, .(L, R)]
        dat <- dat[(check_counts), !"check_counts"]
        # dat <- dat[L!=R] # COULD BE USEFUL (SAFER?)!!!!!!!!!!!!
        # Normalize
        dat[, norm:= (umi_counts+1)/sum(umi_counts)*1e6, cdition]
        # Compute FC
        res <- dcast(dat,
                     L+R~cdition,
                     value.var= "norm")
        res <- na.omit(res)
        res[, log2FoldChange:= log2(rowMeans(do.call(cbind, lapply(.SD, function(x) x/input)))), .SDcols= patterns("^screen_rep")]
        # Subtract basal activity (center controls on 0)
        res[, log2FoldChange:= log2FoldChange-median(res[grepl("control", L) & grepl("control", R), log2FoldChange])]
        # Compute expected
        res[, median_L:= {
          sub <- .SD[grepl("control", R)]
          if(nrow(sub)>=5)
            log2(median(2^sub$log2FoldChange)) else
              as.numeric(NA)
        }, L]
        res[, median_R:= {
          sub <- .SD[grepl("control", L)]
          if(nrow(sub)>=5)
            log2(median(2^sub$log2FoldChange)) else
              as.numeric(NA)
        }, R]
        res[, additive:= log2(2^median_L+2^median_R)]
        # SAVE
        res <- na.omit(res[, .(L, R, log2FoldChange, median_L, median_R, additive)])
        fwrite(res, FC_file)
        print(paste0(FC_file, "  -->> DONE"))
      }
  }, FC_file]
}
