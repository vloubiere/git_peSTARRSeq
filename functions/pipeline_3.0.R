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
      .c <- merge(.c[FLAG %in% c(97, 99) & mapq>=20, .(ID, seq, pos)],
                  .c[FLAG %in% c(145, 147) & mapq>=20, .(ID, seq, pos)], 
                  by= "ID",
                  suffixes= c("_L", "_R"))
      # Compute total reads
      stat <- rbindlist(list(L= .c[(filter_pos_L), .(pos= paste0(sort(unique(pos_L)), collapse= ",")), .(gp= gp_L)], 
                             R= .c[(filter_pos_R), .(pos= paste0(sort(unique(pos_R)), collapse= ",")), .(gp= gp_R)]), 
                        idcol = T)
      .c <- .c[filter_pos_L & filter_pos_R, .(L= seq_L,
                                              R= seq_R,
                                              UMI= gsub(".*_([A-Z]{10}).*", "\\1", ID))]
      stat <- rbindlist(list(stat,
                             data.table(.id= "total_reads", value= nrow(.c))), 
                        fill= T)
      # SAVE
      fwrite(.c, umi_counts)
      fwrite(stat, umi_summary)
      print(paste0(umi_counts, " -->>DONE"))
    }
}, .(bam, umi_counts, umi_summary, type, library)]

#--------------------------------------------------------------#
# Merged counts
# Takes counts from separated runs and merge them per condition
# +UMI collapsing with 1nt difference
# +annotation based on library
#--------------------------------------------------------------#
# meta[, {
#   if(all(file.exists(c(summary_counts, pairs_counts, spike_counts, switched_counts))))
#     print(paste0("All merged counts files -->> ALREADY EXISTS")) else
#     {
#       print(paste0("Start ", pairs_counts))
#       # Import all counts file / cdition
#       .c <- lapply(umi_counts, fread)
#       .c <- rbindlist(.c)
#       # Count UMIs and order
#       .c <- .c[, .(umi_N= .N), .(L, R, UMI)]
#       setorderv(.c, "umi_N", order = -1)
#       # Advanced UMI collapsing (>1 diff)
#       .c[, UMI:= {
#         if(.N>1)
#         {
#           res <- rep(as.character(NA), .N)
#           while(anyNA(res))
#           {
#             check <- is.na(res)
#             check[check] <- stringdist(UMI[check][1], 
#                                        UMI[check], 
#                                        method="hamming", 
#                                        nthread= getDTthreads()-1)<=1
#             res[check] <- UMI[check][1]
#           }
#           res
#         }else
#           UMI
#       }, .(L, R)]
#       # Final collapsing
#       .c <- unique(.c[, .(L, R, UMI)])
#       .c <- .c[, .(umi_counts= .N), .(L, R)]
#       # Compute patterns to identify library pairs
#       .ex <- libs[lib_ID %in% c(vllib, Spike_in)]
#       .ex <- .ex[, .(L_pattern= strsplit(sub_lib_L, ";")[[1]],
#                      R_pattern= strsplit(sub_lib_R, ";")[[1]]),(.ex)]
#       .ex <- .ex[, CJ(L_pattern= strsplit(L, ",")[[1]],
#                       R_pattern= strsplit(R, ",")[[1]]), .(lib_ID, L= L_pattern, R= R_pattern)]
#       .ex[, c("L_pattern", "R_pattern"):= .(paste0("_", L_pattern, "_"),
#                                             paste0("_", R_pattern, "_"))]
#       .ex$L <- NULL
#       .ex$R <- NULL
#       .ex[, Spike:= ifelse(lib_ID==Spike_in & !is.na(Spike_in), T, F)]
#       # Check if pair exists and is spike in
#       .ex[, {
#         .c[grepl(L_pattern, L) & grepl(R_pattern, R), type:= ifelse(Spike, "spike-in", "pair")]
#       }, .(L_pattern, R_pattern, Spike)]
#       .c[is.na(type), type:= "switched"]
#       # Generate read summary
#       sum_files <- data.table(file= sam_summary)
#       sum_files[, mapped:= fread(file)[V1=="Mapped_fragments", V2], file]
#       .summary <- sum_files[, .(mapped= sum(mapped), collapsed= sum(.c$umi_counts))]
#       fwrite(.summary, summary_counts)
#       # SAVE
#       fwrite(.c[type=="pair", !"type"], pairs_counts)
#       fwrite(.c[type=="spike-in", !"type"], spike_counts)
#       fwrite(.c[type=="switched", !"type"], switched_counts)
#       print(paste0(pairs_counts, " -->> DONE"))
#     }
# }, .(vllib, Spike_in, summary_counts, pairs_counts, spike_counts, switched_counts)]

#--------------------------------------------------------------#
# Method using counts norm (Tolerates one rep)
#--------------------------------------------------------------#
# if(any(meta$DESeq2))
# {
#   meta[, {
#     if(file.exists(FC_file))
#       print(paste0(FC_file, "  -->> ALREADY EXISTS")) else
#       {
#         res <- dcast(.SD,
#                      L+R~cdition+rep,
#                      value.var= "value",
#                      fill= 0,
#                      sep = "__")
#         
#         sampleTable <- data.table(rn= setdiff(names(res), c("L", "R")))
#         sampleTable[, c("cdition", "rep"):= tstrsplit(rn, "__")]
#         sampleTable <- data.frame(sampleTable, 
#                                   row.names = "rn")
#         
#         # DF
#         res[, rn:= paste0(L, "__", R)]
#         DF <- data.frame(res[, !c("L", "R")], 
#                          row.names = "rn")
#         DF <- DF[rowSums(DF)>20,]
#         dds <- DESeq2::DESeqDataSetFromMatrix(countData= DF,
#                                               colData= sampleTable,
#                                               design= ~rep+cdition)
#         sizeFactors(dds)= DESeq2::estimateSizeFactorsForMatrix(as.matrix(DF[grep("control.*__control.*", rownames(DF)),]))
#         
#         res <- DESeq2::DESeq(dds)
#         
#         # Differential expression
#         FC <- as.data.frame(DESeq2::results(res, contrast= c("cdition", "screen", "input")), keep.rownames= "rn")
#         FC <- as.data.table(FC, keep.rownames= T)[, c("L", "R"):= tstrsplit(rn, "__")][, .(L, R, log2FoldChange, padj)]
#         
#         # Compute expected
#         median_L <- FC[grepl("control", R), .(check= .N>5, median_L= median(log2FoldChange, na.rm = T)), L][(check)]
#         FC[median_L, median_L:= i.median_L, on= "L"]
#         median_R <- FC[grepl("control", L) , .(check= .N>5, median_R= median(log2FoldChange, na.rm = T)), R][(check)]
#         FC[median_R, median_R:= i.median_R, on= "R"]
#         FC[, additive:= log2(2^median_L+2^median_R)]
#         FC[, multiplicative:= median_L+median_R]
#         
#         # Import counts
#         dat <- unique(data.table(file= pairs_counts, 
#                                  rep= DESeq2_pseudo_rep,
#                                  cdition))
#         dat[cdition=="screen", cdition:= paste0(cdition, "_rep", rep)]
#         dat <- dat[, fread(file), (dat)]
#         # Collapse input counts and screen reps
#         dat <- dat[, .(umi_counts= sum(umi_counts)), .(cdition, L, R)]
#         # Filter
#         dat[, check_counts:= sum(umi_counts)>5, .(L, R)]
#         # Cast counts
#         res <- dcast(dat,
#                      L+R~cdition,
#                      value.var= "umi_counts",
#                      fill= 0)
#         cols <- setdiff(names(res), c("L", "R"))
#         res[, check:= rowSums(.SD), .SDcols= cols]
#         res <- res[check>10, !"check"]
#         # Normalize and compute FC
#         res[, (cols):= lapply(.SD, function(x) (x+1)/sum(x)*1e6), .SDcols= cols]
#         res[, log2FoldChange:= log2(rowMeans(do.call(cbind, lapply(.SD, function(x) x/input)))), .SDcols= patterns("^screen_rep")]
#         # Check if individual enhancer is active
#         control_pairs_log2FC <- res[grepl("control", L) & grepl("control", R), log2FoldChange]
#         res[, act_wilcox_L:= {
#           .c <- log2FoldChange[grepl("control", R)]
#           if(length(.c)>5)
#             wilcox.test(.c, control_pairs_log2FC, alternative = "greater")$p.value else
#               as.numeric(NA)
#         }, L]
#         res[, act_wilcox_R:= {
#           .c <- log2FoldChange[grepl("control", L)]
#           if(length(.c)>5)
#             wilcox.test(.c, control_pairs_log2FC, alternative = "greater")$p.value else
#           as.numeric(NA)
#         }, R]
#         # Subtract basal activity (center controls on 0)
#         res[, log2FoldChange:= log2FoldChange-median(control_pairs_log2FC)]
#         # Compute expected
#         median_L <- res[grepl("control", R) , .(check= .N>5, median_L= median(log2FoldChange)), L][(check)]
#         res[median_L, median_L:= i.median_L, on= "L"]
#         median_R <- res[grepl("control", L) , .(check= .N>5, median_R= median(log2FoldChange)), R][(check)]
#         res[median_R, median_R:= i.median_R, on= "R"]
#         # Expected
#         res[, additive:= log2(2^median_L+2^median_R)]
#         res[, multiplicative:= median_L+median_R]
#         # SAVE
#         res <- na.omit(res[, .(L, R, log2FoldChange, median_L, median_R, act_wilcox_L, act_wilcox_R, additive, multiplicative)])
#         fwrite(res, FC_file)
#         print(paste0(FC_file, "  -->> DONE"))
#       }
#   }, FC_file]
# }
